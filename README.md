# TE-clock
This repository contains code for the paper "Transposable element 5mC methylation state of blood cells predicts age and disease" </br>
<<<<<<< HEAD
Due to the large file sizes, this repository contains only code. <br>
Run the script download_data.sh to collect most files required (If download fails, try links in source_data.tsv manually)</br>
Feel free to contact the authors is help is required.

### Script description (processing order)
- download_data.sh:  downloads required infinium array data from GEO
- scripts/utils.R:  includes various utility functions
- scripts/prepare_re_ids.R:  small scripts to associate a convenient unique name to each repeat family found by RepeatMasker, outputs to data/annotation.
- scripts/prepare_array.R:  preprocesses and standardizes array data before the main analysis.
- scripts/prepare_whi.R:  preprocesses and standardizes WHI array data, inclued for transparency however the data is under controlled access and is not downloaded in this code.
- scripts/prepare_atac.R:  preprocesses and standardizes ATAC-seq data before the main analysis.
- scripts/prepare_rna.R:  preprocesses and standardizes RNA-seq data before the main analysis.
- scripts/prepare_rrbs.R:  preprocesses and standardizes RRBS-seq data before the main analysis. Data not downloaded with this code, however the script is included for transparency.
- scripts/paper_part1v2.R:  analysis for part of the paper related to associations between age, disease and methylation changes at TEs in humans, excluing predictors.
- scripts/paper_part2v2.R:  analysis related to age/disease predictors in humans.
- scripts/paper_part3.R:  analysis related to the mouse RRBS clock.
=======
Due to the large file sizes, this repository contains only code. </br>
Feel free to contact the authors is help is required.
>>>>>>> 3411ae85f9f1b07b8cc8ba38894f9b12407add01
