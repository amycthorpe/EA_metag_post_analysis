# Metagenome Analysis

Repository containing Snakemake workflow for:  

* Dereplicating assembled bins
* Assessing bin quality with [CheckM2](https://github.com/chklovski/CheckM2)
* Assigning taxonomy with [GTDBtk](https://github.com/Ecogenomics/GTDBTk)
* Identification of functional traits with [microTrait](https://github.com/ukaraoz/microtrait)
* Identification of metabolic traits with [METABOLIC](https://github.com/AnantharamanLab/METABOLIC) and pathways with [metabolisHMM](https://github.com/elizabethmcd/metabolisHMM)

## Initial processing
Raw metagenomic reads are first processed according to [this workflow](https://github.com/amycthorpe/metag_analysis_EA) to assemble MAGs. 

## Conda

[Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions
```

Getting the repository including sub-modules
```bash
git clone --recurse-submodules git@github.com:amycthorpe/EA_metag_post_analysis
```

Create the main `snakemake` environment

```bash
# create venv
conda env create -f envs/snakemake.yaml -n "snakemake"
conda activate snakemake
```

## How to run

The workflow can be launched as follows

```bash
snake_cores=48    # adjust as needed
snake_jobs=4   # adjust as needed
conda_prefix="/hdd0/susbus/tools/conda_env"

# dry-run
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${conda_prefix} --cores ${snake_cores} --jobs ${snake_jobs} --conda-frontend conda --rerun-trigger mtime -rpn  

# full run
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${conda_prefix} --cores ${snake_cores} --jobs ${snake_jobs} --conda-frontend conda --rerun-trigger mtime -rp
```


## Configs

All config files are stored in the folder `config/`

**Important Note(s)**: 
- Edit the paths to 
    - `data_dir`: `/prj/DECODE/EA_metag_post_analysis/data`
    - `results_dir`: `"/prj/DECODE/ea_analysis/results_20240915`
    - `env_dir`: `/prj/DECODE/EA_metag_post_analysis/envs`
    - ***`databases`***: `/prj/DECODE/ea_analysis/database`

- Provide a `sample_list.tsv` in the *config* folder

## Downstream analysis and visualisation
Outputs from the present workflow can be visualised using the R scripts provided in the [downstream analysis](https://github.com/amycthorpe/biofilm_MAG_analysis) workflow.
