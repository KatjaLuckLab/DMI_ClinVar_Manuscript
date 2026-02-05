
This repository contains scripts used in our manuscript to process mutation data from ClinVar, analyze domain-motif interface (DMI) prediction data, process and analyze BRET data, set up AlphaFold (AF) runs, and analyze the resulting structural models. 

## ClinVar and AlphaMissense data
### Data preprocessing and assembly
`ClinVar_preprocessing_1.ipynb`, `ClinVar_preprocessing_2.ipynb`

### Overlap of ClinVar mutations to DMI fragments
`ClinVar_DMI_overlap.ipynb`

### Mapping of AlphaMissense data to DMI ClinVar data
`AlphaMissense_mapping.Rmd`

## Structural models
### Setting up AF runs
- Selection of DMIs for AF modelling `DMI_AFmodellingSelection.Rmd`
- Definition of domain boundaries `domain_pLDDT_extension.py`
- Fasta file generation for AF modelling `DMI_fastaFiles.py`
- AF sample batch run parameters `HuRI_DMI_GPU_script_GPU1.0_commands.sh`

### Compilation of data from AF structural models
- Retrieval of values from .cif, .pdf and dssp files `Model_data_retrieval.txt`
- Combination of all model metrics `AF_DMI_combine_model_metrics.Rmd`

## BRET data
### Fitting BRET data
`titration_analysis_upload_script_v2.py`
### Computing statistics on BRET data
`DMI_var_BRET_quantification.ipynb`

## Random Forest retests
### Functions used for retesting RF, modified from the original DMI_predictor repository
`DMI_RF_functions.py`
### Re-tests of the RF with variations of features and testing sets
`RandomForest_retests.py`

## Manuscript
### Analysis and graphic generation
`Figures.Rmd`
- DMI predictions and ClinVar data analysis
- BRET results analysis
- AF model metrics analysis
- Variant structural classification

### Preparation of supplementary tables
`Supp_Tables.Rmd`
