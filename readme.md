# Pipeline for drug ranking based on computed pathway scores of disease and healthy samples

Python3 pipeline inspired in the [simdrugs](https://github.com/sepehrgolriz/simdrugs/tree/main) repository, to structure and automatize the data processing, model training, drug-based calibrated pathway scores computation and drug ranking.

## Summary

This pipeline contains the following functions: 
(1) Data processing to handle the tansformations needed to obtain the original pathway scores of the samples according to single sample analysis GSEA
(2) Model training based on the disease and healthy sample pathway scores, to classify them; and computation of the calibrated disease samples pathwa scores according to the interaction among drug and targets found in these pathways
(3) Drug ranking based on the disease samples whose calibrated matrix were responsible to change the trained model decision from disease to healthy state.
            
## Input configuration file:
* The pipeline only needs a configuration file as a parameter:
- Configuration file parameters:
    - folder: working directory
    - expression_file: compressed gene expression file for the desired icgc project, it must be separated by tabulation. The following columns are mandatory: submitted_file_id (sample names), raw_read_count (the read counts without normalization) and gene_id (genes in ensembl or hgnc symbol)
    - type_norm: normalization type (tpm, fpkm or fpkm_uq (upper quartile))
    - identifier: project identifier to be used in the result files
    - labels_file (optional for part 1): file with two columns, one named 'sample' corresponding to the unique values of submitted_sample_id; the second named 'label' corresponding to a disease (or confirmed tumour) (1) or a healthy (0) case
    - trained_model (optional for part 1): file with the trained model to separate healthy and disease cases
    
    - * The "labels_file" parameter is mandatory for the scoring matrix calculation, model traning and drug ranking 
    - * The "trained_model" parameter will skip the model training step in the drug ranking
    
## Usage Instructions
### Preparation:
1. ````git clone https://github.com/YasCoMa/caliscoma_pipeline.git````
2. ````cd caliscoma_pipeline````
3. Create conda environment to handle dependencies: ````conda env create -f drugresponse_env.yml````
4. ````conda activate drugresponse_env````
5. Setup an environment variable named "path_workflow" with the full path to this workflow folder

### Getting data for the running example in the LICA-FR and LIRI-JP projects from ICGC
1. Download the [expression file for LICA-FR](https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LICA-FR/exp_seq.LICA-FR.tsv.gz) and put it in data_icgc folder
2. Download the [expression file for LIRI-JP](https://dcc.icgc.org/api/v1/download?fn=/current/Projects/LIRI-JP/exp_seq.LIRI-JP.tsv.gz) and put it in data_icgc folder
3. For the liri-jp project, the labels file is already processed, to given an example of a project that run all steps proposed by this workflow

### Run analysis
- Run all steps: ````python3 main.py 0 config.json````
- Run only data processing: ````python3 main.py 1 config.json````
- Run only model training & modified pathway score matrix: ````python3 main.py 2 config.json````
- Run only drug ranking: ````python3 main.py 3 config.json````

## Reference

## Bug Report
Please, use the [Issues](https://github.com/YasCoMa/caliscoma_pipeline/issues) tab to report any bug.
