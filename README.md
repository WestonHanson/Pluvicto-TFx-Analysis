<img src="misc/genome.png" width="140" align="left">

# Copy Number Analysis via Chromosome Shannon Entropy

A script designed to process ichorCNA Copy Number files to find the Shannon Entropy of each chromosome of each patient systematically. 
<br/><br/>

## Description
This pipeline processes *[ichorCNA](https://github.com/GavinHaLab/ichorCNA)* copy number data to quantify genomic instability at the chromosome level and correlate it with clinical outcomes (survival, tumor fraction, treatment response). It computes Shannon entropy scores for each chromosome per patient and enables longitudinal analysis across treatment cycles.

*[Shannon entropy](https://en.wikipedia.org/wiki/Entropy_(information_theory))* quantifies disorder in a data distribution. Applied to copy number states, higher entropy indicates greater chromosomal instability — a marker of genomic complexity.

**Workflow**: The pipeline follows two stages: **(1)** Create a copy number lookup table (Process 1), then compute entropy scores using multiple entropy functions (Process 2); **(2)** Perform modular analyses — visualization, statistical modeling, survival analysis, and clinical correlations (Processes 3-16)—using the entropy scores and input data. Processes are independent and can be run individually for flexible workflows. Different entropy score iterations are available in `entropy_equation_functions.py` and custom implementations can be added as needed.

## Inputs
*Inputs can be specified in the input file*
| Input | Definition |
| ----- | ---------- |
| genomic_instability_ulp_files | Absolute path to a csv file containing sample ids (column name `Identifier`; depending on ULP or deep) and absoulte paths of copy number files (column name `curated_solution_path`) |
| genomic_instability_deep_ichor_files | Absolute path to a csv file containing sample ids (column name `Identifier`; depending on ULP or deep) and absoulte paths of copy number files (column name `curated_solution_path`) |
| all_samples_tfx_and_ctdnaq_sheet | Absolute path to a csv file containing columns `Identifier`, `Sample_ID`, `Cycle`, `TFx`, `ploidy`, `yield`, `cfdna-q`, `TFx_category`, `ctdnaq_category` |
| pluvicto_master_sheet | Absolute path to a csv file containing all patient metadata like sample id, treatement cycles, PSA, TFx, etc. (This file is the biggest out of the clinical data files and would have to be updated the most throughout the code depending on what columns you want to use and or name them) |
| complex_sv_calls_and_chr8_entropy | Absolute path to a csv file containing structural varient calls from JaBbA (Jabba) and AmpliconArchitect (AA) and entropy values for chr8 (columns: `sample_id`, `caller`, `SV_type`, `location`, `chr8_entropy`) |
| gene_annotation_list | Absolute path to a txt file containing hg38 gene annotations (columns: `Ensg_ID`, `Enst_ID`, `Gene`, `HGNC_symbol`, `Chr`, `cdsStart`, `cdsEnd`, `txStart`, `txEnd`, `TSS`, `Karyotype_band`, `strand`) |
| gene_matrix_from_titan | Absolute path to a csv file containing the deep whole genome copy number of each gene in each patient (columns: `Sample`, *gene_names) |
| purity_ploidy_from_titan | Absolute path to a txt file containing the deep whole genome purity and ploidy of each sample (columns: `Sample`, `Purity`, `Ploidy`) |
| prostate_specific_cancer_gene_list | Absolute path to a txt file containing a list of prostate cancer specific oncogenes (no column headers, but in order from left to right would be gene name, chromosome, gene start, gene end) |

## Outputs
(Does not include all veriations of outputs) <p>
*Output can be specifed in the input file and processes can be found and edited in main.py*
| Process    | Output Description |
|----------- | ------------------ |
| Process 1 | Outputs the look up tables used in process 2 to create entropy scores. 
| Process 2  | Outputs a CSV file per entropy function added to `entropy_functions_dict`. Saved to `/outputs/data-tables/entropy-tables-<sequencing_type>-<metric_to_use>`.
| Process 3 | Outputs entropy distribution and time sequence plots like entropy fold change or waterfall plots from C1-C2.  
| Process 4 | Outputs a distribution plot for one or more entropy equations.
| Process 5 | Outputs Cox Hazard Ratio Plots and Kaplin-Meier Curves for specified chromosome and any covariates specified. 
| Process 6 | Visualizes the relationship between complex structural varients (sv) and chromosome 8 entropy through different qualitative plots. 
| Process 7 | Logistic regression correlating gene level Titan copy number (CN) gain/loss vs CN neutral across all deep samples using entropy high/low accounting for TFx as predictor values - creates a CSV (sorted by entropy adj p-value). For significant genes, create a bar plot of top 20 genes by entropy coefficients and forrest plot of top 20 genes by entropy odds ratio. Does GSEA on gene list from logistic regression using odds ratio as weighting. 
| Process 8 | Logistic regression correlating entropy and SVs. 
| Process 9 | Compares different datasets between each other to determine if entropy is prognostic or potentially predictive of Pluvicto.
| Process 10 | Linear regression correlating entropy and SVs. 
| Process 11 | Linear regression correlating gene level Titan copy number (CN) gain/loss vs CN neutral across all deep samples using entropy high/low accounting for TFx as predictor values - creates a CSV (sorted by entropy adj p-value). For significant genes, create a bar plot of top 20 genes by entropy coefficients and forrest plot of top 20 genes by entropy odds ratio. Does GSEA on gene list from logistic regression using odds ratio as weighting. 
| Process 12 | Linear mixed effects model for entropy. 
| Process 13 | Creating Manuscript supplemental tables.
| Process 14 | More specific visualizations - an extention of process 3. 
| Process 15 | Linear mixed effects model for specific genes CN increase/decrease over time.
| Process 16 | Forest plot with all chromosomes on it (seperate models on one plot and all in one model).

## Mutable Variables in Processes
- "Global" Mutable variabels:
<p>&emsp;&emsp;&emsp; 1. `sequencing_type`: Determines what files will be used to build entropy.
<p>&emsp;&emsp;&emsp; All subdirectories and function parameters are for the most part mutable to your liking. 

- Mutable variables almost every process:
<p>&emsp;&emsp;&emsp; 1. `metric_to_use`: Copy number column from files used to build entropy.
<p>&emsp;&emsp;&emsp; 2. `entropy_functions_dict`: Dictionary of entropy functions from `python_scripts/entropy_on_other_datasets.py`.
<p>&emsp;&emsp;&emsp; 3. `entropy_csv_list`: List of entropy csvs saved from process 2.

- Process 3:
<p>&emsp;&emsp;&emsp; 1. `split_method`: Creates a TFx cutoff (3% or 10% only).

- Process 4:
<p>&emsp;&emsp;&emsp; 1. `chromosome`: Chromosome to analyize.
<p>&emsp;&emsp;&emsp; 2. `tfx_thresholding_tag`: Tag to add onto files and directory names.
<p>&emsp;&emsp;&emsp; 3. `tfx_threshold`: Threshold to remove samples at.

- Process 5:
<p>&emsp;&emsp;&emsp; 1. `chromosome`: Chromosome to analyize.
<p>&emsp;&emsp;&emsp; 2. `tfx_thresholding_tag`: Tag to add onto files and directory names.
<p>&emsp;&emsp;&emsp; 3. `tfx_threshold`: Threshold to remove samples at.

- Process 7:
<p>&emsp;&emsp;&emsp; 1. Any of the function parameters that has chr8 can be changed to any chromosome.
<p>&emsp;&emsp;&emsp; 2. `normalize_by_ploidy`: If true it normalized by the ploidy of the sample, else it normalizes by the average of the chromosome.
<p>&emsp;&emsp;&emsp; 3. `cn_threshold`: The threshold a states has to be greater than or less than to be concidered a gain/loss respectively (This is the normalized number so 1 is relative to the ploidy or median of the sample meaning that bin does not deviate from the ploidy/median).
<p>&emsp;&emsp;&emsp; 4. `responder_grouping`: Either groups the samples into what is specified or treats them as separate groups.
<p>&emsp;&emsp;&emsp; 1. `direction`: Determines if the function calculates bins that increases or decrease from `cn_threshold`.
<p>&emsp;&emsp;&emsp; 1. `gene_matrix_id_column`: Name of ID column from the Titan gene matrix.
<p>&emsp;&emsp;&emsp; 1. `pluvicto_id_column`: Name of ID column from the pluvicto master sheet.
<p>&emsp;&emsp;&emsp; 1. `tfx_thresholding`: Determines if a 10% TFx threshold gets applied or not.

- Process 8:
<p>&emsp;&emsp;&emsp; 1. `sv_to_filter_out`: A list of structural variants to filter out.
<p>&emsp;&emsp;&emsp; 2. `all_sv_tag`: Tag to add onto files and directory names.
<p>&emsp;&emsp;&emsp; 3. `group_svs`: Determines if SVs are grouped into biologically similar groups.
<p>&emsp;&emsp;&emsp; 4. `full_cohort`: Determines if full cohort is used or only patients at or above the median entropy value.
<p>&emsp;&emsp;&emsp; 5. `tfx_thresholding`: Determines if a 10% TFx threshold gets applied or not.

- Process 9:
<p>&emsp;&emsp;&emsp; 1. `tfx_threshold`: Determines if a 10% TFx threshold gets applied or not.
<p>&emsp;&emsp;&emsp; 2. List of boolean variables to determine which areas of the process will be run.

- Process 10:
<p>&emsp;&emsp;&emsp; 1. `sv_to_filter_out`: A list of structural variants to filter out.
<p>&emsp;&emsp;&emsp; 2. `all_sv_tag`: Tag to add onto files and directory names.
<p>&emsp;&emsp;&emsp; 3. `group_svs`: Determines if SVs are grouped into biologically similar groups.
<p>&emsp;&emsp;&emsp; 4. `full_cohort`: Determines if full cohort is used or only patients at or above the median entropy value.
<p>&emsp;&emsp;&emsp; 5. `tfx_thresholding`: Determines if a 10% TFx threshold gets applied or not.

- Process 11:
<p>&emsp;&emsp;&emsp; 1. Any of the function parameters that has chr8 can be changed to any chromosome.
<p>&emsp;&emsp;&emsp; 2. `normalize_by_ploidy`: If true it normalized by the ploidy of the sample, else it normalizes by the average of the chromosome.
<p>&emsp;&emsp;&emsp; 3. `cn_threshold`: The threshold a states has to be greater than or less than to be concidered a gain/loss respectively (This is the normalized number so 1 is relative to the ploidy or median of the sample meaning that bin does not deviate from the ploidy/median).
<p>&emsp;&emsp;&emsp; 4. `responder_grouping`: Either groups the samples into what is specified or treats them as separate groups.
<p>&emsp;&emsp;&emsp; 1. `direction`: Determines if the function calculates bins that increases or decrease from `cn_threshold`.
<p>&emsp;&emsp;&emsp; 1. `gene_matrix_id_column`: Name of ID column from the Titan gene matrix.
<p>&emsp;&emsp;&emsp; 1. `pluvicto_id_column`: Name of ID column from the pluvicto master sheet.
<p>&emsp;&emsp;&emsp; 1. `tfx_thresholding`: Determines if a 10% TFx threshold gets applied or not.

- Process 12:
<p>&emsp;&emsp;&emsp; 1. `extreme_responders_categories_csv`: Absolute path to sample classification csv (columns: `Sample_ID`, `Group`).
<p>&emsp;&emsp;&emsp; 2. `cycles_to_use`: A list of cycles to use from `pluvicto_master_sheet`.
<p>&emsp;&emsp;&emsp; 3. `tfx_cutoff`: Determines what TFx threshold gets applied.
<p>&emsp;&emsp;&emsp; 4. `response_category`: A grouping of samples based on `extreme_responders_categories_csv` or overall survival from `pluvicto_master_sheet`.
<p>&emsp;&emsp;&emsp; 5. `cycles`: A list of cycles to use from the entropy table, procuded from process 2, used in the analysis.

- Process 13:
<p>&emsp;&emsp;&emsp; 1. 'baseline_sample_entropy_characteristic': A boolean variable to determine which area of the process will be run.
<p>&emsp;&emsp;&emsp; 2. 'complex_sv_csv': Absolute path to a complex SV csv that has each sv call by tool type for each patient (columns: `patient_id`, `AA`, `jabba`).
<p>&emsp;&emsp;&emsp; 3. 'patrick_annotations': Absolute path to a pre-annotated csv file that will be added to the final dataframe  (columns: `Sample_ID`, *any_other_columns).
<p>&emsp;&emsp;&emsp; 4. 'ulp_curated_weighted_by_bin_cohort_normalized_column':  A function where the absolute path to an entropy table, procuded from process 2, can be changed, as well as all other parameters.
<p>&emsp;&emsp;&emsp; 5. 'deep_ichor_v1_weighted_by_bin_cohort_normalized_column': A function where the absolute path to an entropy table, procuded from process 2, can be changed, as well as all other parameters.
<p>&emsp;&emsp;&emsp; 6. 'entropy_dfs': A dictionary of tuples containing entropy dataframes and its unique name which will become that entropies column name.

- Process 14:
<p>&emsp;&emsp;&emsp; 1. `tfx_entropy_scatter_plot`: A boolean variable to determine which area of the process will be run.
<p>&emsp;&emsp;&emsp; 2. `tfx_cutoff`: Determines what TFx threshold gets applied.
<p>&emsp;&emsp;&emsp; 3. 'ulp_curated_weighted_by_bin_cohort_normalized_column':  A function where the absolute path to an entropy table, procuded from process 2, can be changed, as well as all other parameters.
<p>&emsp;&emsp;&emsp; 4. 'deep_ichor_v1_weighted_by_bin_cohort_normalized_column': A function where the absolute path to an entropy table, procuded from process 2, can be changed, as well as all other parameters.
<p>&emsp;&emsp;&emsp; 5. 'entropy_dfs': A dictionary of tuples containing entropy dataframes and its unique name which will become that entropies column name.

- Process 15:
<p>&emsp;&emsp;&emsp; 1. `gene_list`: A dictionary of tuples with three variables (chromosome, start, end) for each gene.
<p>&emsp;&emsp;&emsp; 2. `extreme_responders_categories_csv`: Absolute path to sample classification csv (columns: `Sample_ID`, `Group`).
<p>&emsp;&emsp;&emsp; 3. `cycles_to_use`: A list of cycles to use from `pluvicto_master_sheet`.
<p>&emsp;&emsp;&emsp; 4. `tfx_cutoff`: Determines what TFx threshold gets applied.
<p>&emsp;&emsp;&emsp; 5. `response_category`: A grouping of samples based on `extreme_responders_categories_csv` or overall survival from `pluvicto_master_sheet`.

- Process 16:
<p>&emsp;&emsp;&emsp; 1. `tfx_cutoff`: Determines what TFx threshold gets applied.
<p>&emsp;&emsp;&emsp; 2. `z_score`: Determines if entropy data will be normalized by z-score or not.

## Scripts
This pipeline has R scripts that are not build into the pipeline but was used to understand / explore areas of the dataset.
| File Name  | File Description |
|----------- | ---------------- |
| `main.py` | The pipelines main file that is called when using the pipeline. |
| `python_scripts/analysis_of_pluvicto_vs_other_dataset.py` | Functions to analyize different datasets together focusing on entropy. |
| `python_scripts/complex_sv_analysis.py` | Functions to analyize correlations between complex structural variants (SVs) and entropy or SVs and survival. |
| `python_scripts/create_copy_number_dataframe.py` | A function to build copy number profile look up tables for each patient. |
| `python_scripts/entropy_by_chromosome.py` | Functions to calculate entropy files and analyize them. |
| `python_scripts/entropy_equation_functions.py` | Entropy equations as unique callable functions. |
| `python_scripts/entropy_on_other_datasets.py` | Functions to calculate entropy on other datasets and to analyize them. |
| `python_scripts/imports.py` | List of all imports throughout the pipeline. |
| `python_scripts/mock_chr8_for_entropy.py` | A small pipeline where you can list out a chromosome level CN profile and see step by step the entropy calculation. |
*Files that are not listed are not important*

## Instalations
This script can be run on any machine that has the python v.3.12.3 or higher (on rhino use Python/3.12.3-GCCcore-13.3.0).

## Usage
```
pythyon main.py --process[#] [True/False]

# Example (Processes are set to False by default so this case would run process 1 and 3 only )
pythyon main.py --process1 True --process3 True
```

## Requirements
This repo has a docker image attached that can be used if needed (it might be depreciated by now).

## Acknowledgments
This script was developed by Weston Hanson in the Gavin Ha Lab, Fred Hutchinson Cancer Center, under the supervision of Robbert D. Patton and Patrick McDeed.

## License
Copyright 2025 Fred Hutchinson Cancer Center

Permission is hereby granted, free of charge, to any government or not-for-profit entity, or to any person employed at one of the foregoing (each, an "Academic Licensee") who obtains a copy of this software and associated documentation files (the ÒSoftwareÓ), to deal in the Software purely for non-commercial research and educational purposes, including the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or share copies of the Software, and to permit other Academic Licensees to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.Ê
No Academic Licensee shall be permitted to sell or use the Software or derivatives thereof in any service for commercial benefit. For the avoidance of doubt, any use by or transfer to a commercial entity shall be considered a commercial use and will require a separate license with Fred Hutchinson Cancer Center.

THE SOFTWARE IS PROVIDED ÒAS ISÓ, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.