# _In vitro_ model of the human small intestine

This repository contains code for the analysis of scRNA-seq of two tissue models combined by hashtag oligos (HTO) in one dataset (control).

## Directory structure

* `bin/` - Scripts used for the analysis
* `envs/` - **conda** environment YAML files
* `data/` - Data output (not included in **git**)
* `analysis/` - Analysis output (not included in **git**)

## Analysis workflow

* Control
   - Clone the **git** repository
   ```
   git clone https://github.com/saliba-lab/tissue-model-human-intestine.git
   ```
   
   - Install [conda](https://docs.conda.io/en/latest/miniconda.html)
   
   - Create **conda** environment
   ```
   conda env create -f envs/default.yml
   ```
   
   - Run the scripts
   ```
   Rscript bin/dataset-control.R
   Rscript bin/control-analysis-qc.R
   Rscript bin/control-analysis-core.R
   Rscript bin/control-analysis-celltypes.R
   ```
