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
| time point: | Sets what time points will be used in FGA vs TFX scatter plots (seperate by commas) |
| genomic instability and tumor fraction data processing: | Calculates genomic instability and tumor fraction of each sample
| central depth data processing: | Grabs central depth value for each sample's TFBS site
| combine data tables: | Combines FGA/TFx and central depth data tables *(necessary for downstream analysis)*
| genomic instability and tumor fraction scatter plots: | Creates a scatter plot of FGA vs TFx for all patients at C1 (above 3% TFX)
| genomic instability and tumor fraction line plots: | Creates line plots for FGA and/or TFx split into 2 plots (no samples below 3% TFx and if any sample is below 3% TFx) - grouped by patients last time point
| central depth pca plot: | Creates two PCA plots for central depth values - one for sample and feature normalizations, respectively. 
| central depth heatmap: | Creates a heatmap of central depth values from parameters selected in main.py. Dendrogram shows clustering and colored bars on top show TFx grouping and responder vs non-responder groups.
| save heatmap data frames: | Saves the heatmap dataframe into a TSV in ```/data-tables```
| central depth differential activity: | Creates a heatmap similar to *central depth heatmap* above, but only shows significantly differencially expressed TFBS's.
| central depth and tumor fraction density plot: | Density plot of TFx vs central depth values across C1 patients.
| central depth and tumor fraction regression: | Creates 3 plots: 1) a heatmap of signficicantly differentially expressed TFBS's with TFx regressed out. 2) A scatter plot of the residuals vs TFx. And 3) Central depth vs TFx for a specific gene that you can select in main.py.
| central depth and tumor fraction spearman rho plot: | Creates two directories: low and high spearman rho statistics. It calculates spearman rho on each TFBS vs TFx and puts it in its respective directory. Then it creates a distribution plot of spearman rho values.
| cox forest plots of sig tfbs against clinical data: | In progress...


## Outputs 
(Does not include all veriations of outputs) <p>
*Output can be specifed in the input file and processes can be found and edited in main.py*
| Process    | Output File Name                        | Output Description |
|----------- | --------------------------------------- | ------------------ |
| Process 1 | `outputs/data-tables/FGA_data_table.tsv` | Data table showing TFx and FGA for each patient at each time point (as well as the patients progression cycle group) |
| Process 2 | `outputs/data-tables/central_depth_data_table.tsv` | Per-sample central depth table created from the raw "central-depth" feature. Columns include patient_id, time_point, site, and central_depth (also saved as mean_depth for compatibility).
| Process 3 | `outputs/data-tables/combined_data_table.tsv` | Combined data table which merges FGA and central depth on patient_id and time_point. Includes metadata such as progression_cycle, tumor_fraction, and genomic_instability.
| Process 4 | `outputs/tumor%_vs_FGA_grt_0.03tfx.png` | A scatter plot of TFx by FGA across all patients only above 3% TFx.
| Process 5 | `outputs/mean_<FGA/tfx>_over_time_by_progression_cycle_SEM_filtered.png` | A line plot of the mean of patients at each cycle for FGA or TFx.
|&#x21B3;|`outputs/treatment_cycle_<cycle>_split_by_<fga/tfx>_over_time.png` | A line plot of TFx or FGA for each cycle split by if any samples went below 3% TFx.
| Process 6 | `outputs/central_depth_PCA_normalized_by_<normalization>.png` | PCA scatter plot - Plot shows PC1 vs PC2 for C1 samples and is color-coded by TFx group (low/medium/high). 
| Process 7 | `outputs/normalized_<normalization>_<number/all>_tfbs_by_variance_heatmap.png` | Clustered central depth heatmap image produced. Heatmap includes top TFBS by variance and annotation bars for progression group and tumor group. 
|&#x21B3;| `outputs/data-tables/heatmap_data_table.tsv` | Heatmap-ready table saved; contains normalized central depth values (per chosen normalization level), patient metadata, and progression/tumor group annotations. 
| Process 8 | `outputs/sig_diff_exp_sites_normalized_<normalization>_from_<number/all>_by_variance_heatmap.png` | Clustered central depth heatmap image produced. Heatmap includes only significantly differentially expressed TFBS's. Has annotation bars for progression group and tumor group. 
| Process 9 | `outputs/tfx_vs_CD_scatter_plot.png` | A density plot of TFx vs CD to see how correlated TFx is with central depth.
| Process 10 | `outputs/diff_sig_exp_residuals_normalized_<normalization>_<number/all>_tfbs_CI_<pvalue_cutoff>_heatmap.png` | Clustered central depth heatmap image produced. Heatmap includes only significantly differentially expressed TFBS residuals - meaning TFx has been regressed out of central depth. Has annotation bars for progression group and tumor group.
|&#x21B3;| `outputs/tfx_central_depth_residual_plot.png` | A scatter plot of TFx vs central depth with a pearson coefficient value to interpret if TFx is correlated with central depth.
|&#x21B3;|  `outputs/tfx_vs_CD_scatter_plot_specific_genes.png` | A scatter plot for a specific gene of TFx vs central depth with a pearson coefficient value to interpret if TFx is correlated with central depth.
| Process 11 | `outputs/rho-stat-high` | A scatter plot for every TFBS that had a rho stat greater than or equal to 0.5 - TFx vs central depth - only significant values.
|&#x21B3;| `outputs/rho-stat-low` | A scatter plot for every TFBS that had a rho stat less than 0.5 - TFx vs central depth - only significant values.
|&#x21B3;| `outputs/temp_histogram.png` | A histogram of the distribution of spearman rho values acorss significant values. 
| Process 12 | *In process...* |

!! Most scatter plots or heatmaps are "responder" vs "non-responder" and that is a user defined variable based on cycles to progression that you can change in each process in main.py !!

## Usage
This script can be run on any machine that has the python v. 3.12.3 or higher.

### Inputs to script:
```
pythyon main.py
```

## Requirements
This repo has a docker image attached that will need to be used.

## Acknowledgments
This script was developed by Weston Hanson in the Gavin Ha Lab, Fred Hutchinson Cancer Center, under the supervision of Robbert D. Patton and Patrick McDeed.

## License
The MIT License (MIT)

Copyright (c) 2025 Fred Hutchinson Cancer Center

Permission is hereby granted, free of charge, to any government or not-for-profit entity, or to any person employed at one of the foregoing (each, an "Academic Licensee") who obtains a copy of this software and associated documentation files (the “Software”), to deal in the Software purely for non-commercial research and educational purposes, including the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or share copies of the Software, and to permit other Academic Licensees to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

No Academic Licensee shall be permitted to sell or use the Software or derivatives thereof in any service for commercial benefit. For the avoidance of doubt, any use by or transfer to a commercial entity shall be considered a commercial use and will require a separate license with Fred Hutchinson Cancer Center.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
