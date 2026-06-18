# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 04/14/26
# Purpose: Calculates a penalized Shannon entropy score per chromosome and outputs a tsv listing patients in decending order by FGA per chromosome.
# Note: Functions adopted from / uses from entropy_by_chromosome.py.

from python_scripts.imports import *

def create_list_of_files_from_directory(directory, sequencing_type, tool):
    list_of_files = []
    if sequencing_type == "NA":
        return
    else:
        # Finds each patient in the directory
        for patient in os.listdir(directory):
            # Creates a full path with directory and patient identifier
            # If patient is not in the patient dictionary add the patient, else print that patient is already in the list
            if "titan" in tool:
                if os.path.isdir(os.path.join(directory, patient)):
                    full_path = os.path.join(directory, patient)
                    for dir in os.listdir(full_path):
                        if "solution" in dir:
                            full_path = os.path.join(full_path, dir)
                    if os.path.isdir(full_path):
                       list_of_files.append(f"{full_path}/{patient}.cna.seg")
            
            else:
                full_path = os.path.join(directory, patient)
                # If full path is true directory append patient identifier and genome altered to different lists
                if os.path.isdir(full_path):
                    list_of_files.append(f"{full_path}/{patient}.cna.seg")

    return list_of_files

def process_entropy_per_chromosome_modified(directory, patient_dict, metric_to_use, sequencing_type, tool, entropy_equation_function=None):
    '''
    Parameters:
    -----------
        directory: String
            Genomic instability directory path.

        patient_dict: Dictionary
            Dictionary to store patient data with patient identifiers as keys.

        metric_to_use: String
            Parameter to determine which metric_to_use from cna.seg file to use for entropy.

        sequencing_type: String
            Determins what sequencying type to select for (e.g. 'deep', 'ulp', or 'ulp_curated').

        tool: String
            Determins what tool was used to process CNs (to seperate how to process file structure).

        entropy_equation_function: function
            Function to calculate entropy.

    Function:
    ---------
        Note: Modifed from process_entropy_per_chromosome function in entropy_by_chomosome.py
        - Loops through each directory in the genomic instability directory array.
        - For each patient found, extracts entropy per chromosome data using find_entropy_per_chromosome.
        - Reads tumor fraction data from the tumor fraction file.
        - Updates the patient dictionary with tumor fraction and entropy values.

    Returns:
    --------
        None (modifies patient_dict in place)

    '''
    print("Processing FGA_TFX...")

    if sequencing_type == "NA":
        curated_files_df = pd.read_csv(directory[0])
        for index in curated_files_df.index:
            full_path = curated_files_df.iloc[index]["curated_solution_path"]
            patient = curated_files_df.iloc[index]["Identifier"]
            patient_dict[patient] = (None, find_entropy_per_chromosome(full_path, patient, metric_to_use, sequencing_type, None, entropy_equation_function, None), None)
    else:
        # Finds each patient in the directory
        for patient in os.listdir(directory):
            # Creates a full path with directory and patient identifier
            # If patient is not in the patient dictionary add the patient, else print that patient is already in the list
            if patient not in patient_dict:
                if "titan" in tool:
                    if os.path.isdir(os.path.join(directory, patient)):
                        full_path = os.path.join(directory, patient)
                        for dir in os.listdir(full_path):
                            if "solution" in dir:
                                full_path = os.path.join(full_path, dir)
                        if os.path.isdir(full_path):
                            patient_dict[patient] = (None, find_entropy_per_chromosome(full_path, patient, metric_to_use, sequencing_type, None, entropy_equation_function, ".cna.seg"), None)
                
                else:
                    full_path = os.path.join(directory, patient)
                    # If full path is true directory append patient identifier and genome altered to different lists
                    if os.path.isdir(full_path):
                        if "ulp_curated_seg" in sequencing_type:
                            patient_dict[patient] = (None, find_entropy_per_chromosome_seg_file(full_path, patient, metric_to_use, sequencing_type, None, entropy_equation_function, ".seg.txt"), None)
                        else:
                            patient_dict[patient] = (None, find_entropy_per_chromosome(full_path, patient, metric_to_use, sequencing_type, None, entropy_equation_function, ".cna.seg"), None)

            else:
                print(f"{patient} already in list")

def create_entropy_data_table_per_chromosome_modified(patient_dict, output_sub_dir, csv_file_name, output_path, metric_to_use, sequencing_type, cycle_bool):
    '''
    Parameters:
    -----------
        patient_dict: Dictionary
            Dictionary containing patient data with tumor fraction, entropy per chromosome, and cycle information.

        csv_file_name: String
            Name of the output CSV/TSV file.

        output_path: String
            Path to the output directory where the data table will be saved.

        metric_to_use: String
            Parameter to show which metric was used in entropy equation.

        cycle_bool: boolean
            Determins if you split the patient id by an underscore.
    
    Function:
    ---------
        Note: Modifed from create_entropy_data_table_per_chromosome function in entropy_by_chomosome.py
        - Extracts entropy values from the patient dictionary for each chromosome.
        - Creates a header row with chromosome columns (chr1-chr22, chrX).
        - Writes a CSV file with patient id, cycle, and entropy values per chromosome.
        - Saves the file to the output path under data-tables directory.

    Returns:
    --------
        None (creates and writes CSV file)

    '''
    # Create directory
    file_dir = os.path.join(output_path, f'data-tables/{output_sub_dir}')
    os.makedirs(file_dir, exist_ok=True)

    # Sets output path and name
    csvFileName = os.path.join(file_dir, f'{csv_file_name}')

    # patient_pregression_cycle_dict = find_T_cycles(patient_dict)

    # Creates headers for data table
    tsv_header_name = [f"chr{i}" for i in range(1,23)]
    tsv_header_name.append("chrX")
    tsv_header_name.insert(0, "patient_id")
    tsv_header_name.insert(1, "cycle")

    dict_of_patients = {}
    for patient_id, (tumor_fraction, genomic_instability, _) in patient_dict.items():
        for chrom, value in genomic_instability.items():
            if patient_id not in dict_of_patients.keys():
                dict_of_patients[patient_id] = [(chrom, value)]
            else:
                dict_of_patients[patient_id].append((chrom, value))

    with open(csvFileName, 'w', ) as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(tsv_header_name)

        for patient_id, values in dict_of_patients.items():
            if cycle_bool:
                if "_" in patient_id:
                    patient_id_specific = patient_id[:patient_id.find("_")]
                else:
                    patient_id_specific = patient_id
            else:
                patient_id_specific = patient_id

            # Remove trailing cluster value
            if sequencing_type == "deep":
                patient_id = patient_id[:patient_id.rfind("_")] 
            
            patient_time_point = patient_id[patient_id.rfind("_")+1:]
            list_of_entropies = [value[1] for value in values]
            csv_writer.writerow([patient_id_specific, patient_time_point, *list_of_entropies])

    print(f"Saved to: {csvFileName}")

def process_entropy_for_csv(csv_path, patient_dict, id_col, chr_col, metric_to_use, sequencing_type, tool, entropy_equation_function=None, sub_set_col=None, sub_set_filter=None):
    '''
    Parameters:
    -----------
        csv_path: String
            Path to csv with CN data.

        patient_dict: Dictionary
            Dictionary to store patient data with patient identifiers as keys.

        id_col: String
            Name of id column. 
        
        chr_col: String
            Name of chromosome column.

        metric_to_use: String
            Parameter to determine which metric_to_use from cna.seg file to use for entropy (and name of copy number column).

        sequencing_type: String
            Determins what sequencying type to select for (e.g. 'deep', 'ulp', or 'ulp_curated').

        tool: String
            Determins what tool was used to process CNs (to seperate how to process file structure).

        entropy_equation_function: function
            Function to calculate entropy (default: None).

        sub_set_col: String
            Name of the column to subset down to (default: None).
        
        sub_set_filter: String
            String to filter sub_set_col down to (default: None).

    Function:
    ---------
        - Creates an entropy score for every patient, chromosome pair.

    Returns:
    --------
        None (modifies patient_dict in place)
    '''
    print("Processing FGA_TFX...")

    # Import csv file
    df = pd.read_csv(csv_path)

    # Subset df
    if sub_set_col is not None and sub_set_filter is not None:
        df = df[df[sub_set_col] == sub_set_filter].copy()

    # Group df by id-chr pair
    df_grouped = df.groupby([id_col, chr_col], as_index=False)[metric_to_use].agg(list)

    # Update patient_dict 
    for idx, row in df_grouped.iterrows():
        id = row[id_col]
        chrom = row[chr_col]
        cn_list = row[metric_to_use]        
        
        # Calculate entropy
        clean_values = [v for v in cn_list if v != "NA"]
        len_of_chrom = len(clean_values)
        counts = Counter(clean_values)
        counts_probs = {int(key): count / len_of_chrom for key, count in counts.items()}

        patient_dict.setdefault(id, {})[chrom] = entropy_equation_function(counts_probs)

    # Sorts the dictionary by keys so it is in natural order
    def chr_key(chrom):
        """Extract number from chromosome name for natural sorting"""
        chrom_str = str(chrom).replace('chr', '').replace('Chr', '')
        try:
            return int(chrom_str)
        except ValueError:
            # For X, Y, M - put them at the end
            return float('inf')

    # Sorts patient_dict chr's into natural order and adds Nones in a triplet to match formatting needed to process into csv
    for id, dictionary in patient_dict.items():
        patient_dict[id] = (None, dict(sorted(dictionary.items(), key=lambda x: chr_key(x[0]))), None)

def add_column_to_dataframe(df_1_path, df_1_id_col, df_2_path, df_2_id_col, df_2_value_col):
    '''
    Parameters:
    -----------
        df_1_path: String
            Path to csv.

        df_1_id_col: String
            Name of column to merge on.

        df_2_path: String
            Path to csv. 
        
        df_2_id_col: String
            Name of column to merge on.

        df_2_value_col: String
            Name of column to add onto df_1.

    Function:
    ---------
        - Merges df_1 and df_2 on id columns and keeps only df_1 columns + df_2_value_col
        - Returns merged df.

    Returns:
    --------
        df_merged_filtered: pandas DataFrame
            df_1 with df_2_value_column merged on id columns.
    '''
    df_1 = pd.read_csv(df_1_path)
    df_2 = pd.read_csv(df_2_path)

    df_merged = pd.merge(df_1, df_2, left_on=df_1_id_col, right_on=df_2_id_col, how='inner')

    df_merged_filtered = df_merged.loc[:, df_merged.columns.isin(list(df_1.columns) + [df_2_value_col])]

    # Format column names to have no spaces
    df_merged_filtered.columns = [col.replace(' ', '_') for col in df_merged_filtered.columns]
    
    return df_merged_filtered

def main():

    ##########################
    # DEFINE INPUTS
    ##########################

    # Define output path
    output_path = "./outputs"

    # List of entropy functions
    def base_entropy(counts_probs):
        return sum(val * math.log2((1/val)) for i, val in counts_probs.items())
        
    def base_entropy_hn_normalized(counts_probs):
        entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
        return pow(2, entropy) / math.log2(len(counts_probs) + 1)
    
    def base_entropy_hn_normalized_logged(counts_probs):
            entropy = sum(val * math.log2((1/val)) for i, val in counts_probs.items())
            entropy_score = pow(2, entropy) / math.log2(len(counts_probs) + 1)
            return math.log2(entropy_score)

    # Build dictionary for other projects (a list of path to directory, tool used (for file structure - None means its a csv), and outcome data)
    projects = {
        # 'radium223': ['/fh/fast/ha_g/projects/Collaborations/Choudhury_Lab/Radium223_2025/GenomicsAnalysis/ULP_ichorCNA_curated/Radium-223/', 'titan', '/fh/fast/ha_g/projects/Collaborations/Choudhury_Lab/Radium223_2025/sample_manifest.csv'],
        'docetaxel': ['/fh/fast/ha_g/user/dchen4/ichorCNA_analysis/docetaxel_cabazitaxel_results/results/ichorCNA/', 'ichor', '/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/2022_09_13_Docetaxel_cabazitaxel_database_with_TFx_for_Gavin.csv'],
        # 'WCDT': ['/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/jci_insight_161370_sdt1_S1I.csv', None, '/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/012326_WCDT_mCRPC_outcomes_w_assay_info.csv']
    }

    sequencing_type = "ulp_curated_seg_entropy_mod"
    metric_to_use = "Corrected_Copy_Number"

    # ----------------------------
    # Process dataframes if needed 
    # ----------------------------

    # --- WCDT ---
    if 'WCDT' in projects:
        # Add purity column to WCDT outcome data
        wcdt_df = add_column_to_dataframe(df_1_path=projects['WCDT'][2], df_1_id_col='patient_id', df_2_path='/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/jci_insight_161370_sdt1_S1B.csv', df_2_id_col='Sample', df_2_value_col='Tumor purity')

        save_path = '/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/012326_WCDT_mCRPC_outcomes_w_assay_info_with_overall_survival.csv'
        if not os.path.isfile(save_path):
            wcdt_df.to_csv(save_path, index=False)
        
        # Update projects dictionary 
        projects['WCDT'][2] = save_path

    # --- radium223 ---
    if 'radium223' in projects:
    # Add purity column to WCDT outcome data
        radium_df = add_column_to_dataframe(df_1_path=projects['radium223'][2], df_1_id_col='Pre_treatment_specimen', df_2_path='/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/radium223_metadata.csv', df_2_id_col='Pre_ID', df_2_value_col='Date_Death')
        
        # Add an event column
        # radium_df["Completed_Tx_binary"] = (radium_df['Completed_Tx'] == 'Y').astype(int)
        radium_df["Completed_Tx_binary"] = np.where(
            radium_df["Completed_Tx"] == "Censored",
            np.nan,
            (radium_df["Completed_Tx"] == "Y").astype(int)
        )

        save_path = '/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/radium223_manifest_with_overall_survival.csv'
        if not os.path.isfile(save_path):
            radium_df.to_csv(save_path, index=False)
        
        # Update projects dictionary 
        projects['radium223'][2] = save_path

    #########################################
    # PROCESS PROJECT DIRECTORIES FOR ENTROPY 
    #########################################

    if args.process1 == "True":

        # Build entropy function dictionary 
        entropy_functions_dict = {
            "base_entropy_hn_normalized_per_chr_table": base_entropy_hn_normalized,
        }

        for projects_idx, (project, arr) in enumerate(projects.items()):
            directory = arr[0]
            tool = arr[1]
            for entropy_idx, (function_label, function) in enumerate(entropy_functions_dict.items()):
                patient_dict = {}

                entropy_equation_function = function

                if tool is None:
                    process_entropy_for_csv(csv_path=directory, patient_dict=patient_dict, id_col="Sample", chr_col="Chromosome", metric_to_use=metric_to_use, sequencing_type=sequencing_type, tool=tool, entropy_equation_function=entropy_equation_function, sub_set_col='Dataset', sub_set_filter='WCDT101')

                else:
                    process_entropy_per_chromosome_modified(directory=directory, patient_dict=patient_dict, metric_to_use=metric_to_use, sequencing_type=sequencing_type, tool=tool, entropy_equation_function=entropy_equation_function)

                if project == 'docetaxel':
                    create_entropy_data_table_per_chromosome_modified(patient_dict=patient_dict, output_sub_dir=f'other-datasets/{project}/entropy-tables-{sequencing_type}-{metric_to_use}', csv_file_name=f"{function_label}.csv", output_path=output_path, metric_to_use="Corrected_Copy_Number", sequencing_type=sequencing_type, cycle_bool=False)
                    
                else:
                    create_entropy_data_table_per_chromosome(patient_dict=patient_dict, output_sub_dir=f'other-datasets/{project}/entropy-tables-{sequencing_type}-{metric_to_use}', csv_file_name=f"{function_label}.csv", output_path=output_path, metric_to_use="Corrected_Copy_Number", sequencing_type=sequencing_type)

    ##########################
    # VISUALIZE ENTROPY SCORES
    ##########################

    if args.process2 == "True":

        # Define entropy functions you want to process
        entropy_functions_dict = {
            "base_entropy_hn_normalized": "base_entropy_hn_normalized_per_chr_table",
        }

        for projects_idx, (project, arr) in enumerate(projects.items()):

            # Some datasets dont actually have a cycle but they do have a cycle column
            if project in ['radium223']:
                cycle_filter = "pre"
            else:
                cycle_filter = None

            for idx, (label, csv_file) in enumerate(entropy_functions_dict.items()):

                sub_dir = f"other-datasets/{project}/process2_outputs/distribution-of-{label.replace("_", "-")}"

                print("Creating plots...")

                create_distribuition_plot_by_chromosomes(
                    output_path, 
                    f'other-datasets/{project}/entropy-tables-{sequencing_type}-{metric_to_use}',
                    sub_dir, 
                    "Corrected_Copy_Number",
                    sequencing_type,
                    f"Entropy Distribution ({label.replace("_", " ")}) per Patient", 
                    f"{csv_file}.csv", 
                    cycle_filter
                )

    ##########################
    # SURVIVAL ANALYSIS
    ##########################
    if args.process3 == "True":
        # Define entropy functions you want to process
        entropy_functions_list = [
            "base_entropy_hn_normalized_per_chr_table",
        ]

        for projects_idx, (project, arr) in enumerate(projects.items()):
            outcome_data_path = arr[2]

            if outcome_data_path is not None and project != "docetaxel":

                 # Some datasets dont actually have a cycle but they do have a cycle column
                if project == 'radium223':
                    patient_id_col = 'Pre_treatment_specimen'
                    cycle_filter = "pre"
                    event_col = 'Completed_Tx_binary'
                    time_to_event_col = 'Date_Death'
                    covariate_list = ["preTFx"]

                elif project == 'WCDT':
                    patient_id_col = 'patient_id'
                    cycle_filter = None
                    event_col = 'OS_status'
                    time_to_event_col = 'OS_months_from_Bx'
                    covariate_list = []

                for csv_path in entropy_functions_list:
                    cox_forest_output_path = f'entropy-analysis-{sequencing_type}-Corrected_Copy_Number/other-datasets/{project}/cox_hardard_ratios/entropy_plots/{csv_path}'
                    km_curve_output_path = f'entropy-analysis-{sequencing_type}-{metric_to_use}/other-datasets/{project}/kaplan-meier-curves/entropy_km_curves/{csv_path}'

                    _, _, chr_column, cox_model_summary = cox_proportional_hazard_model_from_csvs(
                        output_path,
                        cox_forest_output_path,
                        f"data-tables/other-datasets/{project}/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", 
                        outcome_data_path, 
                        "chr8", 
                        cycle_filter, 
                        event_col, 
                        time_to_event_col, 
                        covariate_list,
                        "Corrected_Copy_Number",
                        sequencing_type,
                        None,
                        None,
                        'TFx_C1',
                        0.10
                    )

                    # Format HR and CI from cox model to add to KM curve
                    row = cox_model_summary.index[cox_model_summary['covariate'] == 'chr8'][0]
                    row = cox_model_summary.iloc[row]
                    cox_model_hr_text = f'HR: {row["HR"]:.2f} (CI 95%: {row["HR_lower"]:.2f}–{row["HR_upper"]:.2f})'

                    df = pd.read_csv(outcome_data_path)
                    # Set the index to patient_id
                    df.set_index(patient_id_col, inplace=True)

                    kaplan_meier_plot(
                        output_path, 
                        km_curve_output_path,
                        df,
                        event_col, 
                        time_to_event_col,
                        "chr8", 
                        chr_column,
                        cycle_filter, 
                        "median",
                        "Corrected_Copy_Number",
                        sequencing_type,
                        cox_model_hr_text
                    )

            # If project does not have OS a logistic regression needs to be computed instead of a survival analysis
            elif project == "docetaxel":
                # Read in outcome_data and entropy data and mege
                outcome_data = pd.read_csv(outcome_data_path, index_col=0, encoding="latin1")

                for csv_path in entropy_functions_list:
                    # output path
                    log_reg_output_path = f'entropy-analysis-{sequencing_type}-{metric_to_use}/other-datasets/{project}/process3_outputs/logistic_regression/{csv_path}'

                    # Get entropy column and merge it with outcome_data
                    entropy_column = get_entropy_column(output_path, f"data-tables/other-datasets/{project}/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "chr8", None)
                    merged_df = pd.merge(outcome_data, entropy_column, left_on='Sample label', right_on='patient_id', how='inner')

                    # Get the TFx column from outcome_data
                    tfx_column = merged_df[['patient_id', 'tumor_fraction']]

                    # Subset merged_df to only columns you want and format
                    merged_df.rename(columns={'Resistance; Y = YES; N = NO': 'resistance', 'PSA decline ³ 50% within 16 weeks after docetaxel start; Y = YES; N = NO': 'psa_50_response'}, inplace=True)
                    filtered_df = merged_df[['patient_id', 'resistance', 'psa_50_response', 'chr8']]

                    filtered_df['resistance'] = (filtered_df['resistance'] == 'Y').astype(int)
                    filtered_df['psa_50_response'] = (filtered_df['psa_50_response'] == 'Y').astype(int)

                    # To correct for direction of the logistic regression you need to drop one of the columns and flip 
                    model_results = logistic_regression_with_covariates(filtered_df, 'patient_id', tfx_column, 'chr8', 'tumor_fraction', 'patient_id', None)

                    odds_ratio_forest_plot(
                        output_path, 
                        log_reg_output_path, 
                        model_results['coefficients'], 
                        "entropy_response_odds_ratio", 
                        "entropy_response_adj_pvalue", 
                        None, 
                        None, 
                        None, 
                        f"Odds Ratios for {project} Survival Associated with Chr8 Entorpy (controlled for TFx)", 
                        f"survival_associated_with_chr8_entropy_odds_ratio_forest_plot_FDR_corrected.png"
                    )

    #################################
    # GENE CN AND ENTROPY CORRELATION
    #################################
    if args.process4 == "True":

        for projects_idx, (project, arr) in enumerate(projects.items()):
            directory = arr[0]
            tool = arr[1]
            patient_dict = {}

            if tool is not None:
                list_of_files = create_list_of_files_from_directory(directory=directory, sequencing_type=sequencing_type, tool=tool)

                # Create list of files dataframe
                df = pd.DataFrame({
                    'curated_solution_path': list_of_files
                })

                save_path = f'/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/outputs/data-tables/gene_matrix_from_ichor/{project}_ichor_file_list.csv'
                if not os.path.isfile(save_path):
                    df.to_csv(save_path, index=False)


        # List of csv files to process
        entropy_csv_list = [
            'base_entropy_hn_normalized_per_chr_table',
        ]

        # Dictionary of ichor gene CN matrix path and project meta data
        gene_matrix_path_list = {
            'radium223': ['/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/outputs/data-tables/gene_matrix_from_ichor/radium223_ichor_geneCN.txt', '/fh/fast/ha_g/projects/Collaborations/Choudhury_Lab/Radium223_2025/sample_manifest.csv', '/fh/fast/ha_g/projects/Collaborations/Choudhury_Lab/Radium223_2025/GenomicsAnalysis/ULP_ichorCNA_curated/Radium-223/ichor_curation_summary.csv'],
            'docetaxel': ['/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/scripts/outputs/data-tables/gene_matrix_from_ichor/docetaxel_ichor_geneCN.txt', '/fh/fast/ha_g/user/whanson/PSMA_Lutetium_whanson/genome_instability/data-files/2022_09_13_Docetaxel_cabazitaxel_database_with_TFx_for_Gavin.csv'],
        }

        gene_annotation_list = '/fh/fast/ha_g/grp/reference/GRCh38/gene_lists/GRCh38.p12.ensembl.gene.annotations.sorted.include.AREnhancer.txt'
        prostate_specific_cancer_gene_list = '/fh/fast/ha_g/projects/ProstateTAN/analysis_CN-SV/JASMINE/Plots/data/20230823_Composite_PC_GeneList.txt'

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:

            for projects_idx, (project, arr) in enumerate(gene_matrix_path_list.items()):

                # Process Gene matrix
                filtered_gene_list = filter_gene_annotation_list(gene_annotation_list, "Chr", [8], "Gene")
                if project == 'radium223':
                    cycle = "pre"
                    full_id = True
                    
                    gene_matrix = final_gene_matrix_modified(
                        annotation_df=filtered_gene_list, 
                        annotation_gene_column="Gene", 
                        CN_matrix_path=arr[0], 
                        CN_matrix_id_column="Sample", 
                        cycle=cycle, 
                        meta_data_path=None, 
                        meta_data_id_column=None,
                        ploidy_column=None, 
                        filter_ids=True
                    )

                    # Read ichor ploidy
                    ploidy_df = pd.read_csv(arr[2])
                    ploidy_df.rename(columns={'Sample _ID': 'Sample_ID'}, inplace=True)

                    # Add ploidy column
                    gene_matrix = add_column_to_df(
                        df_1=gene_matrix, 
                        id_column_1="Sample", 
                        df_2=ploidy_df, 
                        id_column_2='Sample_ID', 
                        column_to_add='Ploidy'
                    )

                    # Rename ploidy column
                    gene_matrix.rename(columns={'Ploidy': 'ploidy'}, inplace=True)

                elif project == 'docetaxel':
                    cycle = None
                    full_id = False
                    gene_matrix = final_gene_matrix_modified(
                        annotation_df=filtered_gene_list, 
                        annotation_gene_column="Gene", 
                        CN_matrix_path=arr[0], 
                        CN_matrix_id_column="Sample", 
                        cycle=cycle, 
                        meta_data_path=None, 
                        meta_data_id_column=None,
                        ploidy_column=None, 
                        filter_ids=False
                    )

                    # Read ploidy
                    ploidy_df = pd.read_csv(arr[1], encoding='latin1')

                    # Add ploidy column
                    gene_matrix = add_column_to_df(
                        df_1=gene_matrix, 
                        id_column_1="Sample", 
                        df_2=ploidy_df, 
                        id_column_2='Sample label', 
                        column_to_add='ploidy'
                    )

                # Process Entropy table 
                entropy_column = get_entropy_column(
                    output_path=output_path, 
                    entropy_file_path=f"data-tables/other-datasets/{project}/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", 
                    chr_column="chr8", 
                    cycle=cycle, 
                    full_id=full_id
                )
  
                # Set Parameters for next analyses
                gene_matrix_id_column = "Sample"
                cn_threshold = 1
                responder_grouping = "no_grouping" # 'median' or 'quartiles_extreme'
                direction = "gain" # 'gain' or 'loss'
                tfx_thresholding = False
                encoding=None
                if project == 'docetaxel':
                    meta_data_id_column = "Sample label"
                    tfx_col = 'tumor_fraction'
                    encoding="latin1"
                    normalize_by_ploidy = True
                elif project == 'radium223':
                    meta_data_id_column = "Pre_treatment_specimen"
                    tfx_col = 'preTFx'
                    normalize_by_ploidy = True

                # ----------------------------------
                # Prep DataFrame for Analyses
                # ----------------------------------
                # Create gene matrix/dataframe
                gene_matrix_binary = threshold_gene_matrix(
                    gene_matrix=gene_matrix, 
                    gene_matrix_id_column=gene_matrix_id_column, 
                    gene_matrix_ploidy_column="ploidy", 
                    entropy_df=entropy_column, 
                    entropy_df_id_column="patient_id", 
                    chr_column="chr8", 
                    cn_threshold=cn_threshold, 
                    direction=direction, 
                    normalize_by_ploidy=normalize_by_ploidy
                )

                tfx_df = get_tfx_column(
                    output_path=output_path, 
                    pluvicto_master_sheet_path=arr[1], 
                    id_column=meta_data_id_column, 
                    tfx_column=tfx_col,
                    encoding=encoding
                )

                # Set 10% TFx threshold if tfx_thresholding is True
                tfx_thresholding_tag = ''
                tfx_suffix = ''
                if tfx_thresholding:
                    tfx_df = tfx_df[tfx_df[tfx_col] >= 0.10].reset_index(drop=True)
                    tfx_thresholding_tag = "_tfx_10_percent_threshold"
                    tfx_suffix = f" - {tfx_thresholding_tag.replace('_', ' ').title()}" if tfx_thresholding_tag else ""

                normalization_tag = ''
                normalization_suffix = ''
                if normalize_by_ploidy:
                    normalization_tag = '_normalized_by_ploidy'
                    normalization_suffix = f'(Normalized by Ploidy)'

                # ----------------------------------
                # Logistic Regression
                # ----------------------------------
                model_results = logistic_regression_with_covariates(
                    gene_matrix_binary=gene_matrix_binary, 
                    gene_matrix_id_column=gene_matrix_id_column, 
                    tfx_df=tfx_df, 
                    chr_column='chr8', 
                    tfx_column=tfx_col, 
                    tfx_id_column=meta_data_id_column, 
                    responder_grouping=responder_grouping
                )

                # --------------------------------------
                # Plot metrics from logistic regression
                # --------------------------------------
                # Plot coefficients on bar plots
                coefficient_barplot(
                    output_path=output_path, 
                    sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/other-datasets/{project}/process4_outputs/gene_matrix_and_entropy/{csv_path}", 
                    model_results=model_results['coefficients'], 
                    coefficent_column="entropy_response_coefficient", 
                    adj_pvalue_column="entropy_response_adj_pvalue", 
                    responder_grouping=responder_grouping, 
                    cn_threshold=cn_threshold, 
                    direction=direction, 
                    normalize_by_ploidy=normalize_by_ploidy
                )

                # Plot odds ratios on forest plots
                odds_ratio_forest_plot(
                    output_path=output_path, 
                    sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/other-datasets/{project}/process4_outputs/gene_matrix_and_entropy/{csv_path}", 
                    model_results=model_results['coefficients'], 
                    coefficent_column="entropy_response_odds_ratio", 
                    adj_pvalue_column="entropy_response_adj_pvalue", 
                    responder_grouping=responder_grouping, 
                    cn_threshold=cn_threshold, 
                    direction=direction, 
                    plot_title=f'Odds Ratios for Treatment Response (Median Chr8 Entropy) - CN {direction.title()} {tfx_suffix} {normalization_suffix}', 
                    file_name=f"gene_matrix_and_entropy_response_odds_ratio_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}{tfx_thresholding_tag}{normalization_tag}_{direction}_vs_other_forest_plot.png"
                )

                # Create a gene list from a prostate specific oncogene list
                gene_list = get_prostate_specific_gene_list(gene_list_path=prostate_specific_cancer_gene_list, chr_column="chr8")
                
                #Plot odds ratios on forest plots subsetted to gene list
                odds_ratio_forest_plot(
                    output_path=output_path, 
                    sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/other-datasets/{project}/process4_outputs/gene_matrix_and_entropy/{csv_path}", 
                    model_results=model_results['coefficients'], 
                    coefficent_column="entropy_response_odds_ratio", 
                    adj_pvalue_column="entropy_response_adj_pvalue", 
                    responder_grouping=responder_grouping, 
                    cn_threshold=cn_threshold, 
                    direction=direction, 
                    plot_title=f"Odds Ratios for Treatment Response (Median Chr8 Entropy) - CN {direction.title()} {tfx_suffix} {normalization_suffix}", 
                    file_name=f"gene_matrix_and_entropy_response_odds_ratio_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}{tfx_thresholding_tag}{normalization_tag}_{direction}_vs_other_prostate_gene_list_forest_plot.png", 
                    gene_list=gene_list
                )

                # Create dot plots for GSEA pathways
                try:
                    GSEA(
                        output_path=output_path, 
                        sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/other-datasets/{project}/process4_outputs/gene_matrix_and_entropy/{csv_path}", 
                        model_results=model_results['coefficients'], 
                        gsea_ranking_column="entropy_response_odds_ratio", 
                        adj_pvalue_column="entropy_response_adj_pvalue", 
                        gene_set="MSigDB_Hallmark_2020", 
                        responder_grouping=responder_grouping, 
                        cn_threshold=cn_threshold, 
                        direction=direction
                    )
                except Exception as e:
                    print(f"GSEA could not run: {e}")

                

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plots tumor fraction vs fraction genome altered.")

    parser.add_argument("--process1", default="False", help="Process entropy from CN files for other datasets.")
    parser.add_argument("--process2", default="False", help="Show entropy distributions.")
    parser.add_argument("--process3", default="False", help="Survival Analysis.")
    parser.add_argument("--process4", default="False", help="Gene CN and entorpy correlation.")

    args = parser.parse_args()

    main()