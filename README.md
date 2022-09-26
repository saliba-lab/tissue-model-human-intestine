# tissue model human intestine
An _in vitro_ model of the human small intestine

## Directory structure

* `bin/` - Scripts used for the analysis
* `envs/` - **conda** environment YAML files
* `data/` - Data output (not included in **git**)
* `analysis/` - Analysis output (not included in **git**)

## Analysis workflow

* Control
   - Clone the **git** repository
   
   - Install [conda](https://docs.conda.io/en/latest/miniconda.html)
   
   - Create **conda** environment
   > conda env create -f envs/default.yml
   
   - Run the scripts
   ```
   Rscript bin/dataset-control.R
   Rscript bin/control-analysis-qc.R
   Rscript bin/control-analysis-core.R
   Rscript bin/control-analysis-celltypes.R
   ```
