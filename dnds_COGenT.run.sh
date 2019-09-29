#!/bin/bash

source activate dnds
snakemake --snakefile dnds_COGenT.Snakefile
