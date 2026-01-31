#!/usr/bin/env bash
pixi install && pixi run install-bioc-data && pixi run install-decontam
