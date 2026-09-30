#!/usr/bin/env python3
"""Bounded refactor checks; never invokes Pixi tasks or the full pipeline.

Run from the project environment with Python 3.11+ and R dependencies installed:
    python3 tests/run_readability_checks.py
All generated files and subprocess logs go to a new temporary directory.
"""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "3704c99ef3b9018e88219d9b272f5a86bab6a536"
RSCRIPT = ROOT / ".pixi/envs/default/bin/Rscript"
WORK = Path(tempfile.mkdtemp(prefix="rpa-readability-checks-"))
ENV = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
# Test-mode selection must be explicit in each test, independent of the caller.
for key in ("BETA_DIV_TEST_N", "DIFF_AB_TEST_N", "BETA_DIV_N_CORES"):
    ENV.pop(key, None)
RESULTS = []


def baseline(path):
    return subprocess.check_output(["git", "show", f"{BASELINE}:{path}"], cwd=ROOT)


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    RESULTS.append(label)
    print("PASS", label, flush=True)


def sandbox(name, original=False):
    path = WORK / name
    path.mkdir()
    (path / ".here").touch()
    if original:
        (path / "scripts").mkdir()
        for filename in ("8_beta_diversity_16s_migration.R",
                         "9_differential_abundance_16s.R",
                         "9b_differential_abundance_genus_16s.R"):
            (path / "scripts" / filename).write_bytes(baseline("scripts/" + filename))
    else:
        (path / "scripts").symlink_to(ROOT / "scripts", target_is_directory=True)
    return path


def run_r(path, code, label, success=True, env=None):
    runner = path / (label + ".R")
    runner.write_text(code)
    proc = subprocess.run([str(RSCRIPT), "--vanilla", str(runner)], cwd=path,
                          env={**ENV, **(env or {})}, capture_output=True,
                          text=True, timeout=180)
    text = proc.stdout + proc.stderr
    (path / (label + ".log")).write_text(text)
    if success and proc.returncode != 0:
        raise AssertionError(f"{label} failed; see {path / (label + '.log')}\n{text[-1500:]}")
    return proc.returncode, text


def snapshot(path):
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob("*.rds")}


def static_checks():
    config = tomllib.loads((ROOT / "pixi.toml").read_text())
    old = tomllib.loads(baseline("pixi.toml").decode())
    tasks = config["tasks"]
    check("public tasks retained", old["tasks"].keys() <= tasks.keys())
    check("dependencies unchanged", config["dependencies"] == old["dependencies"])
    check("lockfile unchanged", (ROOT / "pixi.lock").read_bytes() == baseline("pixi.lock"))
    visited = set()

    def visit(name, active):
        assert name not in active, ("task cycle", name)
        if name in visited:
            return
        for dep in tasks[name].get("depends-on", []):
            visit(dep, active | {name})
        visited.add(name)
    for name in tasks:
        visit(name, set())
        task = tasks[name]
        for inp in task.get("inputs", []):
            if inp.startswith("scripts/"):
                assert (ROOT / inp).is_file(), (name, inp)
        cmd = task.get("cmd", "")
        if cmd.startswith("Rscript --vanilla scripts/"):
            script = cmd.removeprefix("Rscript --vanilla ")
            assert script in task["inputs"]
            for helper in re.findall(r'source\(here::here\("scripts", "lib", "([^\"]+)"\)\)',
                                     (ROOT / script).read_text()):
                assert f"scripts/lib/{helper}" in task["inputs"], (name, helper)
    check("task graph, script paths and shared-helper inputs valid", True)
    d = sandbox("parse")
    run_r(d, 'for (f in list.files("scripts", "[.]R$", recursive=TRUE, full.names=TRUE)) parse(f)\n', "parse")
    check("all R files parse", True)
    (d / "baseline").mkdir()
    for name in ("1a_datacleaning_helius.R", "1b_datacleaning_biome.R",
                 "6_alpha_diversity_16s_ethnicity.R", "7a_beta_diversity_16s_compute.R",
                 "8_beta_diversity_16s_migration.R", "9_differential_abundance_16s.R",
                 "9b_differential_abundance_genus_16s.R", "11_alpha_diversity_shotgun.R"):
        (d / "baseline" / name).write_bytes(baseline("scripts/" + name))
    run_r(d, (ROOT / "tests/readability_structure.R").read_text(), "structure")
    check("cleaning and statistical expressions match baseline", True)
    # Reports/figures must not reintroduce the statistical entry points moved out.
    reports = [p for p in (ROOT / "scripts").glob("*_report.R") if not p.name.startswith("07_")]
    reports += list((ROOT / "scripts").glob("*_figure*.R"))
    forbidden = r"\b(?:adonis2|permutest|Maaslin2|kruskal\.test|pairwise\.wilcox\.test|estimate_richness|lm)\s*\("
    for p in reports:
        code = "\n".join(line for line in p.read_text().splitlines() if not line.lstrip().startswith("#"))
        assert not re.search(forbidden, code), p
    check("reports and figures contain no model/test entry points", True)


def missing_cache_checks():
    cases = [
        ("04_air_pollution_participants_report.R", "airpollution-participants-compute"),
        ("06_air_pollution_amsterdam_report.R", "airpollution-amsterdam-prepare"),
        ("08_alpha_diversity_16s_ethnicity_report.R", "alpha-16s-compute"),
        ("09_alpha_diversity_shotgun_ethnicity_report.R", "alpha-shotgun-compute"),
        ("10_beta_diversity_16s_ethnicity_report.R", "10_beta_diversity_16s_ethnicity_compute.R"),
        ("11_beta_diversity_16s_migration_report.R", "beta-16s-migration-compute"),
        ("12_differential_abundance_16s_asv_ethnicity_report.R", "diffabund-16s-compute"),
        ("13_differential_abundance_16s_genus_ethnicity_report.R", "diffabund-16s-genus-compute"),
    ]
    def test(case):
        script, producer = case
        d = sandbox("missing-" + script[:2])
        rc, log = run_r(d, f'source("scripts/{script}")\n', "missing", success=False)
        assert rc != 0 and producer in log and "cache" in log.lower(), log[-1500:]
        return script
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for script in pool.map(test, cases):
            check("missing cache identifies producer: " + script, True)


def migration_check():
    before = sandbox("migration-before", original=True)
    after = sandbox("migration-after")
    run_r(before, (ROOT / "tests/readability_fixture.R").read_text(), "fixture")
    shutil.copytree(before / "data", after / "data")
    def go(item):
        d, files = item
        code = 'set.seed(20260928)\n' + "\n".join(f'source("scripts/{f}", print.eval=TRUE)' for f in files)
        run_r(d, code, "migration")
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(go, [(before, ["8_beta_diversity_16s_migration.R"]),
                           (after, ["11_beta_diversity_16s_migration_compute.R",
                                    "11_beta_diversity_16s_migration_report.R"])]))
    a = before / "results/beta_diversity_migration"
    b = after / "results/beta_diversity_migration"
    names = {p.name for p in a.iterdir() if p.is_file()}
    assert names == {p.name for p in b.iterdir() if p.is_file()}
    for name in names:
        norm = lambda p: re.sub(rb"/[A-Za-z]*Date \([^)]*\)", b"", p.read_bytes()) if p.suffix == ".pdf" else p.read_bytes()
        assert norm(a / name) == norm(b / name), name
    check(f"synthetic migration: {len(names)} baseline CSV/PDF outputs match", True)
    hashes = snapshot(b / "cache")
    # A report run in a fresh process must work without the processed inputs.
    (after / "data").rename(after / "hidden-inputs")
    run_r(after, 'source("scripts/11_beta_diversity_16s_migration_report.R", print.eval=TRUE)\n', "report-only")
    check("migration report reuses unmodified caches without raw/processed inputs", hashes == snapshot(b / "cache"))


def da_reporting_checks():
    fixture = WORK / "migration-before/data"
    for level in ("asv", "genus"):
        for scenario in ("none", "single", "multiple"):
            dirs = [sandbox(f"da-{level}-{scenario}-{v}", original=v == "before")
                    for v in ("before", "after")]
            for d in dirs:
                shutil.copytree(fixture, d / "data")
                shutil.copy2(ROOT / "tests/readability_da_case.R", d / "case.R")
            def go(item):
                d, version = item
                proc = subprocess.run([str(RSCRIPT), "--vanilla", "case.R", version, level, scenario],
                                      cwd=d, env=ENV, capture_output=True, text=True, timeout=120)
                (d / "case.log").write_text(proc.stdout + proc.stderr)
                assert proc.returncode == 0, str(d / "case.log")
            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
                list(pool.map(go, zip(dirs, ("before", "after"))))
            before, after = dirs
            compare_code = f"""
stopifnot(identical(readRDS({json.dumps(str(before / 'compute_rng.rds'))}), readRDS("compute_rng.rds")))
stopifnot(identical(readRDS({json.dumps(str(before / 'summary.rds'))}), readRDS("summary.rds")))
"""
            run_r(after, compare_code, "rng-summary")
            originals = [p for p in (before / "results").rglob("*") if p.suffix in (".pdf", ".tsv", ".csv", ".rds")]
            for p in originals:
                q = after / p.relative_to(before)
                assert q.exists(), q
                if p.suffix == ".rds":
                    run_r(after, f"stopifnot(identical(readRDS({json.dumps(str(p))}),readRDS({json.dumps(str(q))})))", "model-inputs-" + str(len(RESULTS)))
                else:
                    norm = lambda path: re.sub(rb"/[A-Za-z]*Date \([^)]*\)", b"", path.read_bytes()) if path.suffix == ".pdf" else path.read_bytes()
                    assert norm(p) == norm(q), (p, q)
            # The new report may add summary/UpSet tables, but no extra plots.
            assert {p.relative_to(before / "results") for p in (before / "results").rglob("*.pdf")} == {p.relative_to(after / "results") for p in (after / "results").rglob("*.pdf")}
            check(f"{level} {scenario}: model inputs, summary, RNG and reporting match", True)


def skipped_group_checks():
    fixture = WORK / "migration-before/data"
    for level, prefix in (("asv", "12"), ("genus", "13")):
        d = sandbox("skip-groups-" + level)
        shutil.copytree(fixture, d / "data")
        code = f"""
library(phyloseq)
for (site in c("throat", "nose")) {{
    path <- paste0("data/processed/ps_", site, "_rarefied.RDS")
    ps <- readRDS(path)
    ps <- prune_samples(sample_data(ps)$EthnicityTotal == "Turkish", ps)
    saveRDS(ps, path)
}}
# Fail immediately if the group-count guard permits a model fit.
Maaslin2 <- function(...) stop("Unexpected model fit")
source("scripts/{prefix}_differential_abundance_16s_{level}_ethnicity_compute.R")
source("scripts/{prefix}_differential_abundance_16s_{level}_ethnicity_report.R")
stopifnot(length(cache$pairs) == 0L, nrow(cache$summary) == 0L)
stopifnot(length(list.files("results", "[.]pdf$", recursive=TRUE)) == 0L)
"""
        run_r(d, code, "skip-groups")
        check(f"{level}: one qualifying group skips fitting and plots", True)


if __name__ == "__main__":
    print("Artifacts:", WORK, flush=True)
    static_checks()
    helper_dir = sandbox("shared-plot-helpers")
    run_r(helper_dir, (ROOT / "tests/shared_plot_helpers.R").read_text(), "helpers")
    check("shared plot helper edge cases", True)
    missing_cache_checks()
    migration_check()
    da_reporting_checks()
    skipped_group_checks()
    (WORK / "summary.json").write_text(json.dumps(RESULTS, indent=2))
    print("All bounded checks passed. No full pipeline was run.")
