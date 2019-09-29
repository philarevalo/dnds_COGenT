# dnds_COGenT

# Purpose
Simple pipeline for calculating dn/ds between groups of microbial genomes identified with PopCOGenT.

# Dependencies
* A linux-based system
* [Miniconda with python 3.7](https://docs.conda.io/en/latest/miniconda.html)

The required packages can be installed by creating a conda environment with the included `dnds_Environment.yml` file as follows:

`conda env create -f dnds_Environment.yml`


This also requires the `seqinr` and `reshape2` packages in R. To sintall To install, please follow the instructions under above to create the conda environment. Then, activate the environment (`source activate dnds`). Finally, run the `Rscript instal_R_dependencies.R`.

# Usage

To run, set the parameters in the `dnds_COGenT_config.yml` file. Then, run the runscript as follows: `bash dnds_COGenT.run.sh`

Your input sequence files should be in `fasta` format and there should only be one entry per strain in each file. The names of the sequences should just be the strain name as in the test example and these strain names should match the strain names in the population file.

# References

This package uses the [`seqinr`](https://rdrr.io/rforge/seqinr/man/kaks.html) R library to calculate dn/ds which in turn uses the method of Li (1993) and Pamilo & Bianchi (1993). To generate codon-based alignments of proteins, we use [`pal2nal`](http://www.bork.embl.de/pal2nal/) developed by Suyama, Torrents, and Bork (2006). Population calls in the test example were generated with [`PopCOGenT`](https://github.com/philarevalo/PopCOGenT).

Arevalo P., VanInsberghe D., Elsherbini J., Gore J., Polz M.F. (2019). A reverse ecology approachbased on a biological definition of microbial populations. Cell, 178(4).doi:10.1016/j.cell.2019.06.03

Li, W.-H. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J. Mol. Evol., 36:96-99.

Pamilo, P., Bianchi, N.O. (1993) Evolution of the Zfx and Zfy genes: Rates and interdependence between genes. Mol. Biol. Evol, 10:271-281

Suyama, M., Torrents, D., Bork, P. (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Res. 34, W609-W612.
