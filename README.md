<img src="misc/genome.png" width="140" align="left">

# Genomic Instibility vs Tumor Fraction or Central Depth Over Time Pipeline for cfDNA ULP-WGS Data 

A script developed to process and give insight into genomic instibilty (fraction genome altered : FGA) vs tumor fraction (TFx) or central depth data over time points for cell-free DNA (cfDNA) ultra-low pass whole genome sequencing (ULP-WGS) data.
<br/><br/>

## Desprition
This script integrates features from *[ichorCNA](https://github.com/GavinHaLab/ichorCNA)*, *[Triton](https://github.com/GavinHaLab/TritonNP)*, and clinical progression data to derive logitudinal changes in genomic instabily over time, specifically understanding its relationship to tumor fraction and transciption factor binding site expression. 

## Inputs
*Inputs can be specified in the input file*
| Input | Definition |
| ----- | ---------- |
| genomic instability file paths: | Paths to ichorCNA batch run *directories* |
| tumor fraction file paths: | Path to tumor fraction meta data |
| patient treatment progression: | Path to progression cycle meta data |
| central depth file path: | Path to Triton central depth file |
| complex sv calls and chr8 entropy: | Path to csv file containing structural varient calls and entropy values for chr8 |
| time point: | Sets what time points will be used in FGA vs TFX scatter plots (seperate by commas) |


## Outputs 
(Does not include all veriations of outputs) <p>
*Output can be specifed in the input file and processes can be found and edited in main.py*
| Process    | Output Description |
|----------- | ------------------ |
| Process 1  | Outputs a CSV file per entropy function added to `entropy_functions_dict`. Saved to `/outputs/data-tables/entropy-tables-<sequencing_type>-<metric_to_use>`. |

| Process 2 | Outputs entropy distribution and time sequence plots like entropy fold change or waterfall plots from C1-C2.  | 

| Process 3 | Outputs a distribution plot for one or more entropy equations. |

| Process 5 | Outputs Cox Hazard Ratio Plots and Kaplin-Meier Curves for specified chromosome and any covariates specified. |

| Process 4 | Visualizes the relationship between complex structural varients (sv) and chromosome 8 entropy through box and swarm plots (`significant_svs` subdir for individual histograms of significant relationships between svs). |

|&#x21B3;| Visualizes the relationship between complex structural varients (sv) and chromosome 8 entropy through histograms (stacked and individual histograms). |

|&#x21B3;| Visualizes the relationship between complex structural varients (sv) and chromosome 8 entropy through stacked boxplot per patient. |

| Process 5 | Logistic regression predicting gene level Titan copy number (CN) gain/loss vs CN neutral across all deep samples using entropy high/low accounting for TFx as predictor values - creates a CSV (sorted by entropy adj p-value). For significant genes, create a bar plot of top 20 genes by entropy coefficients and forrest plot of top 20 genes by entropy odds ratio. Does GSEA on gene list from logistic regression using odds ratio as weighting. 

| Process 6 | Through a binary dataframe of presence of complex structural variances (SVs) by patient (used throughout Process 20), box plots of presence of any SVs vs no presense, and indivudual SVs, in relation to entropy scores are created.

|&#x21B3;| Plots forest plot of odds ratios of SVs that are significantly correlated with entropy values from logistic regression.

|&#x21B3;| Cox regression model forest plot accounting for TFx, quantifying the added prognostic information provided by SV presence on survival.

|&#x21B3;| Kaplan-Meier curves seperated on different SVs.

|&#x21B3;| Above analyses includes sub-directories corresponding to different data subsets, which can be selected in main.py by modifying the full_cohort and tfx_thresholding parameters.

!! Most scatter plots or heatmaps are "responder" vs "non-responder" and that is a user defined variable based on cycles to progression that you can change in each process in main.py !!

## Usage
This script can be run on any machine that has the python v.3.12.3 or higher (on rhino use Python/3.12.3-GCCcore-13.3.0).

## Instalations
None at the moment.

### Inputs to script:
```
pythyon main.py
```

## Requirements
This repo has a docker image attached that will need to be used.

## Acknowledgments
This script was developed by Weston Hanson in the Gavin Ha Lab, Fred Hutchinson Cancer Center, under the supervision of Robbert D. Patton and Patrick McDeed.

## License
Copyright 2025 Fred Hutchinson Cancer Center

Permission is hereby granted, free of charge, to any government or not-for-profit entity, or to any person employed at one of the foregoing (each, an "Academic Licensee") who obtains a copy of this software and associated documentation files (the ÒSoftwareÓ), to deal in the Software purely for non-commercial research and educational purposes, including the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or share copies of the Software, and to permit other Academic Licensees to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.Ê
No Academic Licensee shall be permitted to sell or use the Software or derivatives thereof in any service for commercial benefit. For the avoidance of doubt, any use by or transfer to a commercial entity shall be considered a commercial use and will require a separate license with Fred Hutchinson Cancer Center.

THE SOFTWARE IS PROVIDED ÒAS ISÓ, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
