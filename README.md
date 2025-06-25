
This repository contains scripts used in our manuscript to harvest mutation data, handle DMI prediction data, assemble AF modelling input data, compile it and analyse it. 

## ClinVar and AlphaMissense data
### Data preprocessing
ClinVar_preprocessing_1.ipynb
ClinVar_preprocessing_2.ipynb

### Overlap of ClinVar mutations to DMI fragments
ClinVar_DMI_overlap.ipynb

### Mapping of AlphaMissense data to DMI ClinVar data
AlphaMissense_mapping.Rmd

## Structural models
### AF model data assembly
- Selection of DMI for AF modelling `DMI_AFmodellingSelection.Rmd`
- Define extended domain boundaries by pLDDT `domain_pLDDT_extension.py`
- Produce fastas for AF modelling `DMI_fastaFiles.py`

### AF model info data compilation
- Retrieval of values of .cif, .pdf and dssp files `Model_data_retrieval.txt`
- Combine all model metrics `AF_DMI_combine_model_metrics.Rmd`

## Manuscript
### Analysis and graphic generation
`Figures.Rmd`

### Preparation of supplementary tables
`Supp_Tables.Rmd`
