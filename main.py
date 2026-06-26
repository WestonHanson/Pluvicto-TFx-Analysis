# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 12/01/25
# Purpose: Main script to analyize Pluvicto samples.

from python_scripts.imports import *

# get_input_file function
# Creates empty arrays to hold directories and time points of patients
# Returns genomic instability directory array, tumor fraction directory array, patient progression directory array, and time point array
def get_input_file(input_path):
    # Creates empty dictionary to hold directories and time points of patients
    input_map = {}

    input_file = os.path.join(input_path, "input.csv")

    # Opens input file
    with open(input_file, "r") as file:
        for line in file:
            delimited_line = line.strip().split(",")

            # Saves tags (eg. "file paths:")
            tag = delimited_line[0].strip().lower()

            # Adds file paths to directory array and time points to time point array
            if tag == "genomic_instability_ulp_files:":
                input_map["genomic_instability_ulp_files"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == 'genomic_instability_deep_ichor_files:':
                input_map["genomic_instability_deep_ichor_files"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == "all_samples_tfx_and_ctdnaq_sheet:":
                input_map["all_samples_tfx_and_ctdnaq_sheet"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == "pluvicto_master_sheet:":
                input_map["pluvicto_master_sheet"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == "complex_sv_calls_and_chr8_entropy:":
                input_map["complex_sv_calls_and_chr8_entropy"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == "gene_annotation_list:":
                input_map["gene_annotation_list"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == "gene_matrix_from_titan:":
                input_map["gene_matrix_from_titan"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == "purity_ploidy_from_titan:":
                input_map["purity_ploidy_from_titan"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == "prostate_specific_cancer_gene_list:":
                input_map["prostate_specific_cancer_gene_list"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

    return input_map

# main function
# Initializes the patient dictionary and calls the get_input_file function to get the directories and time points
# If the process argument is True, it processes all samples and creates a data table
# Always creates a scatter plot and a line plot over time points
def main(args):
    # Dictionary with key as patient name and values as (genome instability, tumor fraction)
    patient_dict = {}

    # Inpput directory
    input_path = "./inputs"
    
    # Output directory and plot
    output_path = "./outputs"

    # Gets all directories for searching ulp, deep, tumor fraction, patient progression, and time point data from input file
    # genomic_instability_directory_array, tumor_fraction_directory_array, patient_progression_directory_array, central_depth_directory_array,time_point_array = get_input_file(input_path)
    input_map = get_input_file(input_path)

    sequencing_type = "deep_ichor_seg_transitions"

    # Processes samples if user selects True
    if sequencing_type not in ["ulp", "deep", "ulp_curated", "ulp_curated_seg", "ulp_curated_seg_entropy_mod", "ulp_curated_seg_subclone", "ulp_curated_seg_transitions", "ulp_curated_seg_transitions_mod", "deep_seg", 'deep_seg_entropy_mod', 'deep_seg_500kb', 'deep_seg_10snp', 'deep_seg_20snp', 'deep_ichor_seg_v1', 'deep_ichor_seg_v2', 'deep_ichor_seg_v1_entropy_mod', 'deep_ichor_seg_v2_entropy_mod', 'deep_ichor_v1', 'deep_ichor_seg', 'deep_ichor_seg_transitions']:
        print("Wrong sequencing type!")
        return    

    # ==============
    # Process 1:
    # ==============
    if args.process1 == "True":

        metric_to_use = "Corrected_Copy_Number"

        input_dirs = input_map[f"genomic_instability_{sequencing_type.split('_seg')[0]}_files"]

        patient_df = create_copy_number_profile_dataframe(
            dir_list=input_dirs, 
            metric_to_use=metric_to_use, 
            sequencing_type=sequencing_type
        )

        save_path = f'{output_path}/data-tables/cohort-copy-number-profiles-dataframes'
        file_name = f'{sequencing_type}-{metric_to_use}-cohort-copy_number_profiles'

        os.makedirs(save_path, exist_ok=True)

        patient_df.to_pickle(
            f'{save_path}/{file_name}.pkl',
        )

        print(f'Saved cohort copy number profiles to {save_path}/{file_name}.pkl')

    # ===============
    # Process 2:
    # ===============
    if args.process2 == "True":

        # Dictionary of functions (list of entropy functions in entropy_equation_functions.py)
        entropy_functions_dict = {
            "base_entropy_per_chr_table": base_entropy,
            "base_entropy_exponentiate_per_chr_table": base_entropy_exponentiate,
            "base_entropy_hn_cohort_normalized_per_chr_table": base_entropy_hn_cohort_normalized,
        }

        metric_to_use = "Corrected_Copy_Number"

        # Create variable to find pickle 
        if 'seg' in sequencing_type and 'deep' not in sequencing_type:
            seq_type = 'ulp_curated_seg'
        elif 'seg' not in sequencing_type:
            seq_type = 'ulp_curated'
        else:
            seq_type = 'deep_ichor_seg'

        # Read in cohort copy number profile pickle as a dataframe
        patient_cn_df = pd.read_pickle(
            f'{output_path}/data-tables/cohort-copy-number-profiles-dataframes/{seq_type}-{metric_to_use}-cohort-copy_number_profiles.pkl'
        )

        for index, (function_label, function) in enumerate(entropy_functions_dict.items()):

            # Reset patient_dict
            patient_dict = {}

            print(f"Processing {function_label}")

            max_count = 0
            
            for (patient, chromosome), group in patient_cn_df.groupby(['sample_id', 'chr']):
                # Instantiate matrix with 0s
                count_matrix = np.zeros((1000, 1000), dtype=int)
                
                # Reset group index
                group.reset_index(drop=True, inplace=True)

                prev = None
                for i in range(0, len(group)):
                    # Get current copy number
                    curr = group.iloc[i]['Corrected_Copy_Number']

                    # Add number of self transitioning bins
                    count_matrix[curr, curr] += group.iloc[i]['num.mark'] 
                    
                    if 'transitions' in sequencing_type:
                        count_matrix[curr, curr] -= 1
                        if prev is not None and prev != curr:
                            count_matrix[prev, curr] += 1

                    prev = curr

                # Find the probabilities of the matrix
                total = count_matrix.sum()
                if total != 0:
                    probability_matrix = count_matrix / total
                else: 
                    probability_matrix = count_matrix

                # Find the amount non-zero cells in the matrix
                count = np.count_nonzero(count_matrix)
                max_count = count if count > max_count else max_count

                # Set the entropy of the patient/chromosome pair
                if patient not in patient_dict:
                    patient_dict[patient] = {}
                patient_dict[patient][chromosome] = float(function(probability_matrix))

            # Sort dict
            for patient, chrom_dict in patient_dict.items():
                def chr_key(chrom):
                    """Extract number from chromosome name for natural sorting"""
                    chrom_str = str(chrom).replace('chr', '').replace('Chr', '')
                    try:
                        return int(chrom_str)
                    except ValueError:
                        # For X, Y, M - put them at the end
                        return float('inf')
                chrom_dict = dict(sorted(chrom_dict.items(), key=lambda x: chr_key(x[0])))
                patient_dict[patient] = chrom_dict

            # Update patient_dict if base_entropy_hn_cohort_normalized_per_chr_table
            if function_label == "base_entropy_hn_cohort_normalized_per_chr_table":
                for patient_data in patient_dict.values():
                    for chrom in patient_data:
                        patient_data[chrom] /= max_count

            create_entropy_data_table_per_chromosome_v2(patient_dict, f'entropy-tables-{sequencing_type}-{metric_to_use}', f"{function_label}.csv", output_path, "Corrected_Copy_Number", sequencing_type)

    # ===============
    # Process 3:
    # ===============
    if args.process3 == "True":

        entropy_functions_dict = {
            "base_entropy_hn_cohort_normalized": "base_entropy_hn_cohort_normalized_per_chr_table",
        }

        for idx, (label, csv_file) in enumerate(entropy_functions_dict.items()):
            
            metric_to_use = "Corrected_Copy_Number"
            # split_method = "favorable_unfavorable_ctdnaq_grt_than_3pct_TFx"
            split_method = "favorable_unfavorable_ctdnaq_grt_than_10pct_TFx"

            print(f"Visualizing {label}")
            
            # *************************************************************
            # Create violin plot & histograms for entropy equation
            # *************************************************************

            print(f"\t↳ Creating chr level histogram and violin plots")

            sub_dir = f'process3_outputs/distribution-of-{label.replace("_", "-")}'

            create_distribuition_plot_by_chromosomes(
                output_path, 
                f'entropy-tables-{sequencing_type}-{metric_to_use}',
                sub_dir, 
                "Corrected_Copy_Number",
                sequencing_type,
                f'Entropy Distribution ({label.replace("_", " ")}) per Patient', 
                f"{csv_file}.csv", 
                "C1"
            )

            # *************************************************************
            # Create line plot across cycles for entropy equation
            # *************************************************************

            print(f"\t↳ Creating line plot across cycles")

            sub_dir = f"process3_outputs/entropy_across_cycles/{label}"
            file_name = f'{label}_across_cycles_line_plot'

            # Combine pluvicto_master_sheet, all_samples_tfx_and_ctdnaq_sheet, and entropy_df to pass into line plot function
            master_sheet_df = pd.read_csv(input_map["pluvicto_master_sheet"][0], index_col=0)
            ctdna_sheet_df = pd.read_csv(input_map["all_samples_tfx_and_ctdnaq_sheet"][0], index_col=0)
            entropy_df = pd.read_csv(f'{output_path}/data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_file}.csv')

            # Subset ctdna_sheet_df to C1 and remove '-' from column names
            ctdna_sheet_df_filtered = ctdna_sheet_df[ctdna_sheet_df['Cycle'] == 'C1']
            ctdna_sheet_df_filtered = ctdna_sheet_df_filtered.rename(columns={'cfdna-q': 'cfdna_q'})

            df_combined = pd.merge(master_sheet_df, ctdna_sheet_df_filtered, left_on='Sample_ID', right_on='Sample_ID')
            df_combined = pd.merge(df_combined, entropy_df, left_on='Sample_ID', right_on='patient_id')

            df_combined = df_combined.drop(columns={'Sample_ID'})

            if '_grt_than_10pct_TFx' in split_method:
                file_name = file_name + "_grt_than_10pct_TFx"
                df_combined = df_combined[df_combined['TFx_C1'] > 0.10]

            # Filter down to desired columns
            df_filtered = df_combined[['patient_id', 'chr8', 'cycle', 'ctdnaq_category']]

            if '_grt_than_3pct_TFx' in split_method:
                file_name = file_name + "_grt_than_3pct_TFx"
            
            line_plot_across_cycles(
                output_path,
                sub_dir,
                "Corrected_Copy_Number",
                sequencing_type,
                "chr8",
                'ctdnaq_category',
                split_method,
                f'Entropy Distribution ({label.replace("_", " ")}) per Patient', 
                file_name,
                df_filtered,
                'cycle',
                ["C1", "C2", "C3", "C4", "C5", "C6"]
            )

            # *************************************************************
            # Create bar plot of C1-C2 entropy change
            # *************************************************************

            print(f"\t↳ Creating C1-C2 stacked bar plot")

            # ! Sub dir defined in outer for loop !

            # Read in csvs and combine
            ctdna_sheet_df = pd.read_csv(input_map["all_samples_tfx_and_ctdnaq_sheet"][0], index_col=0)
            entropy_df = pd.read_csv(f'{output_path}/data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_file}.csv')
            combined_df = pd.merge(ctdna_sheet_df, entropy_df, left_on=['Sample_ID', 'Cycle'], right_on=['patient_id', 'cycle'], how='inner')

            # Drop unnessasary columns
            combined_df = combined_df.drop(columns={'Sample_ID', 'Cycle'})

            if '_grt_than_10pct_TFx' in split_method:
                file_name = file_name + "_grt_than_10pct_TFx"
                df_combined = df_combined[df_combined['TFx'] > 0.10]

            # Rename columns
            combined_df = combined_df.rename(columns={'cfdna-q': 'cfdna_q'})

            for values in ['cfdna_q', 'chr8', 'TFx']:

                if values == 'chr8':
                    sub_dir = f"process3_outputs/entropy_across_cycles/{label}"
                    file_name = f'{label}_change_barplot'
                else:
                    sub_dir = f"process3_outputs/entropy_across_cycles/TFx_and_ctDNAQ_reference"
                    file_name = f'{values}_change_barplot'
                    if os.path.isfile(f'{output_path}/entropy-analysis-{sequencing_type}-{metric_to_use}/{sub_dir}/{file_name}'):
                        continue

                if '_grt_than_3pct_TFx' in split_method:
                    file_name = file_name + "_grt_than_3pct_TFx"

                if '_grt_than_10pct_TFx' in split_method:
                    file_name = file_name + "_grt_than_10pct_TFx"

                for direction in ["both", "only_increasing", "only_decreasing", "only_the_same"]:

                    bar_plot_comparing_cycles(
                        output_path,
                        sub_dir,
                        "Corrected_Copy_Number",
                        sequencing_type,
                        values,
                        'patient_id',
                        'ctdnaq_category',
                        split_method,
                        f'Comparing Entropy Over Cycles ({label.replace("_", " ")}) per Patient', 
                        file_name,
                        combined_df,
                        'cycle',
                        ["C1", "C2"],
                        direction
                    )

            # *************************************************************
            # Create violin plot looking at FC from C1-C2
            # *************************************************************

            print(f"\t↳ Creating C1-C2 fold change violin plot (favorable vs unfavorable)")

            sub_dir = f"process3_outputs/entropy_across_cycles/{label}"

            file_name = f'{label}_C2_response_log2_FC_violin_plot'

            if '_grt_than_3pct_TFx' in split_method:
                file_name = file_name + "_grt_than_3pct_TFx"

            # Read in csvs and combine
            ctdna_sheet_df = pd.read_csv(input_map["all_samples_tfx_and_ctdnaq_sheet"][0], index_col=0)
            entropy_df = pd.read_csv(f'{output_path}/data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_file}.csv')
            combined_df = pd.merge(ctdna_sheet_df, entropy_df, left_on=['Sample_ID', 'Cycle'], right_on=['patient_id', 'cycle'], how='inner')

            if '_grt_than_10pct_TFx' in split_method:
                file_name = file_name + "_grt_than_10pct_TFx"
                combined_df = combined_df[combined_df['TFx'] > 0.10]

            # Drop unnessasary columns
            combined_df = combined_df.drop(columns={'Sample_ID', 'Cycle'})

            groups = ["C1", "C2"]

            entropy_fc_df = calculate_fold_change_across_groups_per_patient(
                "chr8",
                'patient_id',
                'cycle',
                groups,
                combined_df
            )

            # Filter to the last cycle
            ctdna_sheet_df_filtered = ctdna_sheet_df[ctdna_sheet_df['Cycle'] == groups[1]]

            # Add response grouping to combined_df and drop redundent columns
            entropy_fc_df = pd.merge(entropy_fc_df, ctdna_sheet_df_filtered[['Sample_ID', 'ctdnaq_category']], left_on='patient_id', right_on='Sample_ID', how="left")
            entropy_fc_df = entropy_fc_df.drop(columns={'Sample_ID'})

            # Pivot so column names are ctdna_categorys
            entropy_fc_df_pivot = entropy_fc_df.pivot(
                index='patient_id',
                columns='ctdnaq_category',
                values='log2_fc'
            )

            fc_groups_combined = pd.DataFrame()
            if 'grt_than_3pct_TFx' in split_method:
                fc_groups_combined['Favorable'] = entropy_fc_df_pivot['Low']
                fc_groups_combined['Unfavorable'] = entropy_fc_df_pivot['Moderate'].combine_first(entropy_fc_df_pivot['High'])
            else:
                fc_groups_combined['Favorable'] = entropy_fc_df_pivot['Undetectable'].combine_first(entropy_fc_df_pivot['Low'])
                fc_groups_combined['Unfavorable'] = entropy_fc_df_pivot['Moderate'].combine_first(entropy_fc_df_pivot['High'])

            # Plot on a violin plot
            publication_multi_violin_plot(
                fc_groups_combined, 
                ['Favorable', 'Unfavorable'], 
                title=f'Log2-FC C2/C1 ({label.replace("_", " ")}) ', 
                xlabel=None, 
                ylabel=None, 
                cycle=None, 
                out_dir=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/{sub_dir}', 
                file_name=file_name, 
                add_significance=True
            )

            # ********************************************************************************************************
            # Create violin plot looking at percent change from C1-C2 (different responder vs non-responder groupings)
            # ********************************************************************************************************

            print(f"\t↳ Creating C1-C2 percent change violin plot (different responder vs non-responder groupings)")

            exteme_response_categories_csv = '/path/to/proteus_deg_sample_classification.csv'

            # Read in csvs and combine
            ctdna_sheet_df = pd.read_csv(input_map["all_samples_tfx_and_ctdnaq_sheet"][0], index_col=0)
            entropy_df = pd.read_csv(f'{output_path}/data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_file}.csv')
            combined_df = pd.merge(ctdna_sheet_df, entropy_df, left_on=['Sample_ID', 'Cycle'], right_on=['patient_id', 'cycle'], how='inner')

            if '_grt_than_10pct_TFx' in split_method:
                file_name = file_name + "_grt_than_10pct_TFx"
                combined_df = combined_df[combined_df['TFx'] > 0.10]

            # Drop unnessasary columns
            combined_df = combined_df.drop(columns={'Sample_ID', 'Cycle'})

            # Rename columns to not include dashes
            combined_df = combined_df.rename(columns={'cfdna-q': 'cfdna_q'})

            cycles = ["C1", "C2"]

            sub_dir = f"process3_outputs/entropy_across_cycles/{label}"
            file_name_cycles_tag = f"cycles_{cycles[0]}_through_{cycles[len(cycles)-1]}"
            file_name = f'{label}_entropy_percent_change_{file_name_cycles_tag}'

            if '_grt_than_3pct_TFx' in split_method:
                file_name = file_name + "_grt_than_3pct_TFx"

            df_pivot = combined_df.pivot(
                index='patient_id',
                columns='cycle',
                values='chr8'
            ).reset_index()

            df_pivot = add_max_change_columns(df_pivot, cycles)

            # Drop NANs in df_pivot
            df_pivot = df_pivot.dropna(subset=['pct_change'])

            # Add new response groupings
            exteme_response_categories_df = pd.read_csv(exteme_response_categories_csv)
            exteme_response_categories_df = exteme_response_categories_df[['Sample_ID', 'Group']]
            exteme_response_categories_df = exteme_response_categories_df.rename(columns={'Group':'Response'})
            df_pivot_groupings = pd.merge(df_pivot, exteme_response_categories_df, left_on='patient_id', right_on='Sample_ID', how='inner').drop(columns={'Sample_ID'})

            df_wide_groupings = df_pivot_groupings.pivot(
                index='patient_id',
                columns='Response',
                values='pct_change'
            ).reset_index()

            df_wide_groupings.set_index('patient_id', inplace=True)

            # Plot on a violin plot
            publication_multi_violin_plot(
                df_wide_groupings, 
                ['Responder', 'Non-Responder'], 
                title=f'Percent Change from C1 to C2 ({label.replace("_", " ")}) ', 
                xlabel=None, 
                ylabel=None, 
                cycle=None, 
                out_dir=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/{sub_dir}', 
                file_name=file_name, 
                add_significance=True
            )

            # *************************************************************
            # Create water fall plot analyizing relative change from C1-C6 
            # *************************************************************

            print(f"\t↳ Creating percent change waterfall plot")

            # ! Sub dir is defined in the for loop !

            # Read in csvs and combine
            ctdna_sheet_df = pd.read_csv(input_map["all_samples_tfx_and_ctdnaq_sheet"][0], index_col=0)
            entropy_df = pd.read_csv(f'{output_path}/data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_file}.csv')
            combined_df = pd.merge(ctdna_sheet_df, entropy_df, left_on=['Sample_ID', 'Cycle'], right_on=['patient_id', 'cycle'], how='inner')

            if '_grt_than_10pct_TFx' in split_method:
                file_name = file_name + "_grt_than_10pct_TFx"
                combined_df = combined_df[combined_df['TFx'] > 0.10]

            # Drop unnessasary columns
            combined_df = combined_df.drop(columns={'Sample_ID', 'Cycle'})

            # Rename columns to not include dashes
            combined_df = combined_df.rename(columns={'cfdna-q': 'cfdna_q'})

            cycles = ["C1", "C2"]

            for values in ['cfdna_q', 'chr8', 'TFx']:

                if values == 'chr8':
                    sub_dir = f"process3_outputs/entropy_across_cycles/{label}"
                    file_name_cycles_tag = f"cycles_{cycles[0]}_through_{cycles[len(cycles)-1]}"
                    file_name = f'{label}_absolute_entropy_change_waterfall_plot_{file_name_cycles_tag}'
                else:
                    sub_dir = f"process3_outputs/entropy_across_cycles/TFx_and_ctDNAQ_reference"
                    file_name_cycles_tag = f"cycles_{cycles[0]}_through_{cycles[len(cycles)-1]}"
                    file_name = f'absolute_{values}_change_waterfall_plot_{file_name_cycles_tag}'
                    if os.path.isfile(f'{output_path}/entropy-analysis-{sequencing_type}-{metric_to_use}/{sub_dir}/{file_name}'):
                        continue

                if '_grt_than_3pct_TFx' in split_method:
                    file_name = file_name + "_grt_than_3pct_TFx"

                df_pivot = combined_df.pivot(
                    index='patient_id',
                    columns='cycle',
                    values=values
                ).reset_index()

                df_pivot = add_max_change_columns(df_pivot, cycles)

                # Drop NANs in df_pivot
                df_pivot = df_pivot.dropna(subset=['pct_change'])

                # Add response groupings
                combined_df_filtered = combined_df[combined_df['cycle'] == "C1"] # Filter down to only C1 responder groupings
                df_pivot_groupings = pd.merge(df_pivot, combined_df_filtered[['patient_id', 'ctdnaq_category']], left_on='patient_id', right_on='patient_id')

                # Add grouping label (and filter out undetectable group if so)
                df_pivot_groupings, group_labels, group_n  = stratify_data(df_pivot_groupings, 'ctdnaq_category', split_method, None)

                # Combine groupings into favorable and unfavorable 
                df_pivot_groupings['group'] = np.where(
                    (df_pivot_groupings['ctdnaq_category'] == 'Undetectable') |
                    (df_pivot_groupings['ctdnaq_category'] == 'Low'),
                    'Favorable',
                    'Unfavorable'
                )

                create_waterfall_plot(
                    df=df_pivot_groupings,
                    value_col='pct_change', 
                    group_col='group', 
                    title=f"Biggest Percent Change of {'Entropy' if values == 'chr8' else values} ({label.replace('_', ' ')}) per Patient", 
                    xlabel="Patients (sorted)", 
                    ylabel="Percent Change", 
                    out_dir=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/{sub_dir}', 
                    file_name=file_name
                )


            print("")
        
    # ===============
    # Process 4:
    # ===============
    if args.process4 == "True":

        csv_file_comapision_dict = {
            "base_entropy_hn_cohort_normalized": ["base_entropy_hn_cohort_normalized_per_chr_table.csv"],
        }

        sub_dir = f'process4_outputs/comparison_of_entropy_equations'

        metric_to_use = 'Corrected_Copy_Number'

        chromosome = 20

        tfx_thresholding_tag = '_10prct_tfx'
        tfx_threshold = 0.1

        pluvicto_df = pd.read_csv(input_map["pluvicto_master_sheet"][0])

        list_of_patients = list(pluvicto_df[pluvicto_df['TFx_C1'] >= tfx_threshold]['Sample_ID'])
        
        for idx, (plot_title, csv_files_list) in enumerate(csv_file_comapision_dict.items()):
            all_data = []
            for file in csv_files_list:
                csv_file_path = os.path.join(output_path, f"data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{file}")
                df = pd.read_csv(csv_file_path)
                df_filtered = df[df['patient_id'].isin(list_of_patients)]
                data = df_filtered[f'chr{chromosome}'].dropna().values
                all_data.append(data)

            plot_title = plot_title + tfx_thresholding_tag

            distribution_plot_comparing_entorpy_equations(
                output_path=f"{output_path}/entropy-analysis-{sequencing_type}-{metric_to_use}/{sub_dir}/chr{chromosome}/tfx_thresholding_{tfx_thresholding_tag}", 
                all_data=all_data,
                data_labels=csv_files_list,
                column=f'chr{chromosome}', 
                sequencing_type=sequencing_type, 
                metric_to_use="Corrected_Copy_Number", 
                file_name=plot_title
            )

        # Can use code, but isnt very informative - saving it here so it doesnt get lost
        # qualitative_top_and_bottom_40_patients(
        #     output_path,
        #     "entropy_per_chromosome_table.csv",
        #     "standard",
        #     "entropy_per_chromosome_table_ploidy_2_corrected.csv",
        #     "corrected",
        #     "Corrected_Copy_Number",
        #     sequencing_type,
        #     "chr8",
        #     "C1",
        #     40,
        #     True
        # )

    # ===============
    # Process 5:
    # ===============
    if args.process5 == "True":

        # List of csv files to process
        entropy_csv_list = [
            "base_entropy_exponentiate_per_chr_table",
            "base_entropy_hn_cohort_normalized_per_chr_table",
        ]

        tfx_thresholding_tag = '_10_prct'
        tfx_threshold = 0.10

        chromosome = 8

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:

            cox_forest_output_path = f'entropy-analysis-{sequencing_type}-Corrected_Copy_Number/cox_hardard_ratios/entropy_plots/{csv_path}/chr{chromosome}/tfx_threshold_{tfx_thresholding_tag}'
            
            _, _, chr_column, cox_model_summary = cox_proportional_hazard_model_from_csvs(
                output_path, 
                cox_forest_output_path,
                f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", 
                input_map["pluvicto_master_sheet"][0], 
                f'chr{chromosome}', 
                "C1", 
                "Death", 
                "survival_days", 
                ["TFx_C1"],
                "Corrected_Copy_Number",
                sequencing_type,
                None,
                None,
                'TFx_C1',
                tfx_threshold,
            )
            
            # Grabs and format HR and CI from cox model 
            cox_model_summary.reset_index(drop=True, inplace=True)
            row = cox_model_summary.index[cox_model_summary['covariate'] == f'chr{chromosome}'][0]
            row = cox_model_summary.iloc[row]
            cox_model_hr_text = f'HR: {row["HR"]:.2f} (CI 95%: {row["HR_lower"]:.2f}–{row["HR_upper"]:.2f})'

            metric_to_use = "Corrected_Copy_Number"
            km_curve_output_path = f'entropy-analysis-{sequencing_type}-{metric_to_use}/kaplan-meier-curves/entropy_km_curves/{csv_path}/chr{chromosome}/tfx_threshold_{tfx_thresholding_tag}'

            master_sheet = input_map["pluvicto_master_sheet"][0]
            ctdna_sheet = input_map["all_samples_tfx_and_ctdnaq_sheet"][0]
            master_sheet_df = pd.read_csv(master_sheet, index_col=0)
            ctdna_sheet_df = pd.read_csv(ctdna_sheet, index_col=0)

            # Subset ctdna_sheet_df to C1 and remove '-' from column names
            ctdna_sheet_df = ctdna_sheet_df[ctdna_sheet_df['Cycle'] == 'C1']
            ctdna_sheet_df = ctdna_sheet_df.rename(columns={'cfdna-q': 'cfdna_q'})

            df_combined = pd.merge(master_sheet_df, ctdna_sheet_df, left_on='Sample_ID', right_on='Sample_ID')

            # Filter by TFx if it is deep
            df_combined = df_combined[df_combined['TFx_C1'] > tfx_threshold]

            # KM curve based on (high TFx, high entropy) vs (high TFx, low entropy) vs (low TFx, high entropy) vs (low TFx, low entropy)
            kaplan_meier_plot(
                output_path, 
                km_curve_output_path,
                df_combined,
                "Death", 
                "survival_days", 
                f'chr{chromosome}',
                chr_column, 
                "TFx_C1",
                "tfx_and_entropy_median",
                metric_to_use,
                sequencing_type,
                None,
                [(0, 1), (2, 3)]
            )

            # Creates a km curve for each of the TFx groups seperated on TFx + low entropy (0), TFx + high entropy (1), and base line TFx ('overall')
            tfx_tags = ['low', 'medium', 'high']
            for tag in tfx_tags:
                try:
                    kaplan_meier_plot(
                        output_path, 
                        km_curve_output_path,
                        df_combined,
                        "Death", 
                        "survival_days", 
                        f'chr{chromosome}',
                        chr_column, 
                        "TFx_C1",
                        f"tfx_and_entropy_custom_{tag}",
                        metric_to_use,
                        sequencing_type,
                        None,
                        [(0, 'overall'), ('overall', 1)],
                        plot_title = None,
                        split_params = None,
                        group_labels  = None,
                        complex_sv = None,
                        high_entropy = None,
                        include_overall = True
                    )
                except:
                    print('Error splitting data')

            
            kaplan_meier_plot(
                output_path, 
                km_curve_output_path,
                df_combined,
                "Death", 
                "survival_days", 
                f'chr{chromosome}',
                chr_column, 
                "TFx_C1",
                "quartile_extremes",
                "Corrected_Copy_Number",
                sequencing_type
            )

            kaplan_meier_plot(
                output_path, 
                km_curve_output_path,
                df_combined,
                "Death", 
                "survival_days", 
                f'chr{chromosome}', 
                chr_column, 
                "TFx_C1",
                "quartiles",
                "Corrected_Copy_Number",
                sequencing_type
            )

            kaplan_meier_plot(
                output_path, 
                km_curve_output_path,
                df_combined,
                "Death", 
                "survival_days", 
                f'chr{chromosome}', 
                chr_column,
                "TFx_C1", 
                "median",
                "Corrected_Copy_Number",
                sequencing_type,
                cox_model_hr_text
            )

            kaplan_meier_plot(
                output_path, 
                km_curve_output_path,
                df_combined,
                "Death", 
                "survival_days", 
                f'chr{chromosome}',
                chr_column, 
                "TFx_C1",
                "tumor_fraction",
                "Corrected_Copy_Number",
                sequencing_type
            )

            # KM Curve based on favorable and unfavorable ctDNA-Q
            kaplan_meier_plot(
                output_path, 
                km_curve_output_path,
                df_combined,
                "Death", 
                "survival_days", 
                f'chr{chromosome}',
                chr_column, 
                "cfdna_q",
                "favorable_ctdnaq_median_entropy",
                "Corrected_Copy_Number",
                sequencing_type
            )
            kaplan_meier_plot(
                output_path, 
                km_curve_output_path,
                df_combined,
                "Death", 
                "survival_days", 
                f'chr{chromosome}',
                chr_column, 
                "cfdna_q",
                "unfavorable_ctdnaq_median_entropy",
                "Corrected_Copy_Number",
                sequencing_type
            )

    # ==============
    # Process 6:
    # ==============
    if args.process6 == "True":
        # List of csv files to process
        entropy_csv_list = [
            'base_entropy_hn_normalized_per_chr_table'
        ]

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:
            complex_svs_and_entropy_df = clean_and_fix_complex_sv_csv(
                input_map["complex_sv_calls_and_chr8_entropy"][0], 
                "sample_id", 
                "caller", 
                "location", 
                "SV_type", 
                "chr8_entropy", 
                ['dup', 'del', 'inv', 'tra', 'invdup']
            )

            complex_svs_and_new_entropy_df = replace_chr8_entropy_column(
                output_path, 
                complex_svs_and_entropy_df, 
                f"{csv_path}.csv", 
                "Corrected_Copy_Number", 
                sequencing_type, 
                "chr8", 
                "C1",
            )

            dir_path = f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process6_outputs/{csv_path}'

            complex_sv_and_entropy_visualization(
                complex_svs_and_new_entropy_df, 
                dir_path, 
                "SV_type", 
                "chr8_entropy", 
                "Corrected_Copy_Number", 
                sequencing_type
            )
            complex_sv_and_entropy_visualization(
                complex_svs_and_new_entropy_df, 
                dir_path, 
                "SV_type", 
                "chr8_entropy", 
                "Corrected_Copy_Number", 
                sequencing_type, 
                "Jabba"
            )
            complex_sv_and_entropy_visualization(
                complex_svs_and_new_entropy_df, 
                dir_path, 
                "SV_type", 
                "chr8_entropy", 
                "Corrected_Copy_Number", 
                sequencing_type, 
                "AA"
            )
            tukey_hsd(
                complex_svs_and_new_entropy_df, 
                dir_path, 
                "SV_type", 
                "chr8_entropy", 
                "Corrected_Copy_Number", 
                sequencing_type
            )
            complex_sv_and_entropy_stacked_histogram(
                complex_svs_and_new_entropy_df, 
                dir_path, 
                "SV_type", 
                "chr8_entropy", 
                "Corrected_Copy_Number", 
                sequencing_type
            )
            complex_sv_and_entropy_stacked_histogram_patient_level(
                complex_svs_and_new_entropy_df, 
                dir_path, 
                "SV_type", 
                "patient_id", 
                "chr8_entropy", 
                "Corrected_Copy_Number",
                sequencing_type
            )

    # ==============
    # Process 7:
    # ==============
    if args.process7 == "True":

        # List of csv files to process
        entropy_csv_list = [
            'base_entropy_hn_normalized_per_chr_table',
        ]

        normalize_by_ploidy = True

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:
                
            # Process Gene matrix
            filtered_gene_list = filter_gene_annotation_list(input_map["gene_annotation_list"][0], "Chr", [8], "Gene")
            gene_matrix = final_gene_matrix(filtered_gene_list, "Gene", input_map["gene_matrix_from_titan"], "Sample", "C1", input_map["purity_ploidy_from_titan"], "Ploidy")

            # Process Entropy table 
            entropy_column = get_entropy_column(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "chr8", "C1")

            # Set Parameters for next analyses
            cn_threshold = 1
            responder_grouping = "no_grouping" # 'median' or 'quartiles_extreme'
            direction = "gain" # 'gain' or 'loss'
            gene_matrix_id_column = "Sample"
            pluvicto_id_column = "Sample_ID"
            tfx_thresholding = True

            # ----------------------------------
            # Prep DataFrame for Analyses
            # ----------------------------------
            # Create gene matrix/dataframe
            gene_matrix_binary = threshold_gene_matrix(gene_matrix, gene_matrix_id_column, "Ploidy", entropy_column, "patient_id", "chr8",  cn_threshold, direction, normalize_by_ploidy)
            tfx_df = get_tfx_column(output_path, input_map["pluvicto_master_sheet"][0], pluvicto_id_column, "TFx_C1")

            # Set 10% TFx threshold if tfx_thresholding is True
            tfx_thresholding_tag = ''
            tfx_suffix = ''
            if tfx_thresholding:
                tfx_df = tfx_df[tfx_df['TFx_C1'] >= 0.10].reset_index(drop=True)
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
                tfx_column='TFx_C1',
                tfx_id_column=pluvicto_id_column,
                responder_grouping=responder_grouping
            )

            # --------------------------------------
            # Plot metrics from logistic regression
            # --------------------------------------
            # Plot coefficients on bar plots
            coefficient_barplot(
                output_path=output_path,
                sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process7_outputs/gene_matrix_and_entropy/{csv_path}/logistic_regression",
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
                sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process7_outputs/gene_matrix_and_entropy/{csv_path}/logistic_regression",
                model_results=model_results['coefficients'],
                coefficent_column="entropy_response_odds_ratio",
                adj_pvalue_column="entropy_response_adj_pvalue",
                responder_grouping=responder_grouping,
                cn_threshold=cn_threshold,
                direction=direction,
                plot_title=f'Odds Ratios for Treatment Response (Median Chr8 Entropy) - CN {direction.title()} {tfx_suffix} {normalization_suffix}',
                file_name=f"gene_matrix_and_entropy_response_odds_ratio_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}{tfx_thresholding_tag}{normalization_tag}_{direction}_vs_other_forest_plot"
            )

            # Create a gene list from a prostate specific oncogene list
            gene_list = get_prostate_specific_gene_list(
                gene_list_path=input_map["prostate_specific_cancer_gene_list"][0],
                chr_column="chr8"
            )
            
            #Plot odds ratios on forest plots subsetted to gene list
            odds_ratio_forest_plot(
                output_path=output_path,
                sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process7_outputs/gene_matrix_and_entropy/{csv_path}/logistic_regression",
                model_results=model_results['coefficients'],
                coefficent_column="entropy_response_odds_ratio",
                adj_pvalue_column="entropy_response_adj_pvalue",
                responder_grouping=responder_grouping,
                cn_threshold=cn_threshold,
                direction=direction,
                plot_title=f"Odds Ratios for Treatment Response (Median Chr8 Entropy) - CN {direction.title()} {tfx_suffix} {normalization_suffix}",
                file_name=f"gene_matrix_and_entropy_response_odds_ratio_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}{tfx_thresholding_tag}{normalization_tag}_{direction}_vs_other_prostate_gene_list_forest_plot",
                gene_list=gene_list
            )

            # Create dot plots for GSEA pathways
            try:
                GSEA(
                    output_path=output_path,
                    sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process7_outputs/gene_matrix_and_entropy/{csv_path}/logistic_regression",
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

    # ==============
    # Process 8:
    # ==============
    if args.process8 == "True":

        # List of csv files to process
        entropy_csv_list = [
            "base_entropy_hn_cohort_normalized_per_chr_table",
        ]

        sv_to_filter_out = ['dup', 'del', 'inv', 'tra', 'invdup']
        all_sv_tag = ''
        # sv_to_filter_out = []
        # all_sv_tag = '_all_svs'

        chromosome_list = list(range(1, 23)) + ['X']

        group_svs = False

        for chromosome in chromosome_list:
            # Set CSV file variable to use in process
            for csv_path in entropy_csv_list:

                complex_svs_and_entropy_df = clean_and_fix_complex_sv_csv(input_map["complex_sv_calls_and_chr8_entropy"][0], "sample_id", "caller", "location", "SV_type", "chr8_entropy", sv_to_filter_out)
                # complex_svs_and_new_entropy_df = replace_chr8_entropy_column(output_path, complex_svs_and_entropy_df, f"{csv_path}.csv", "Corrected_Copy_Number", sequencing_type, "chr8", "C1")

                # Filteres dataframe down to only structural variants found in chr8
                complex_svs_and_new_entropy_df_filtered = filter_to_location_in_complex_sv_df(
                    df=complex_svs_and_entropy_df, 
                    location_column='location', 
                    chr_filter=f'{chromosome}'
                )


                if complex_svs_and_new_entropy_df_filtered.empty:
                    print(f'chr{chromosome} does not have any SVs.')
                    continue

                # --- Create a table that shows a binary call for each patient with respect to each tool that calls the SV --- 

                # Create binary dataframes for each caller
                AA_binary_df = create_binary_dataframe(complex_svs_and_new_entropy_df_filtered[complex_svs_and_new_entropy_df_filtered['caller'] == 'AA'].copy())
                AA_calls = AA_binary_df.columns.drop('patient_id')
                AA_binary_df['AA'] = AA_binary_df[AA_calls].apply(
                    lambda row: ', '.join(row.index[row == 1]),
                    axis=1
                )
                AA_col = AA_binary_df[['patient_id', 'AA']]

                jabba_binary_df = create_binary_dataframe(complex_svs_and_new_entropy_df_filtered[complex_svs_and_new_entropy_df_filtered['caller'] == 'Jabba'].copy())
                print(jabba_binary_df)
                jabba_calls = jabba_binary_df.columns.drop('patient_id')
                jabba_binary_df['jabba'] = jabba_binary_df[jabba_calls].apply(
                    lambda row: ', '.join(row.index[row == 1]),
                    axis=1
                )
                jabba_col = jabba_binary_df[['patient_id', 'jabba']]

                # Merge dataframes
                AA_jabba_combined_df = pd.merge(AA_col, jabba_col, on='patient_id', how='outer')

                # Save combined df
                csv_dir_path = f'{output_path}/data-tables/complex_sv_binary_tables/{csv_path}'
                os.makedirs(csv_dir_path, exist_ok=True)
                csv_file_path = os.path.join(csv_dir_path, f'chr{chromosome}_complex_svs_type_by_tool_table.csv')
                AA_jabba_combined_df.to_csv(csv_file_path, index=False)

                # location_array = [(127735434, 127742951)]
                location_array = [(127000000, 128000000)]
                complex_sv_filtered = filter_by_gene_location(complex_svs_and_new_entropy_df_filtered.copy(), 'location', location_array)

                print(f'Samples with a complex structural varient in: {location_array}')
                print(complex_sv_filtered)

                binary_df = create_binary_dataframe(complex_svs_and_new_entropy_df_filtered)
                binary_df = prepare_sv_df_for_survival_analysis(binary_df, input_map["pluvicto_master_sheet"][0])

                # Create new groupings based off biology and drop all other column (except patient_id)
                if group_svs:
                    binary_df_columns = binary_df.columns
                    binary_df['ecDNA_dm'] = ((binary_df.get('ecDNA', 0) == 1.0) | (binary_df.get('dm', 0) == 1.0)).astype(float)
                    binary_df['bfb_tyfonas_pyrgo'] = ((binary_df.get('bfb', 0) == 1.0) | (binary_df.get('tyfonas', 0) == 1.0) | (binary_df.get('pyrgo', 0) == 1.0)).astype(float)
                    binary_df['chromoplexy_Complex_non_cyclic'] = ((binary_df.get('chromoplexy', 0) == 1.0) | (binary_df.get('Complex_non_cyclic', 0) == 1.0)).astype(float)
                    binary_df['Linear_tic'] = ((binary_df.get('Linear', 0) == 1.0) | (binary_df.get('tic', 0) == 1.0)).astype(float)
                    binary_df = binary_df.drop(columns=[col for col in binary_df_columns if col != "patient_id"])
                    sv_grouped_tag = '_grouped'
                else:
                    if 'ecDNA' in binary_df.columns or 'dm' in binary_df.columns:
                        binary_df['ecDNA_dm'] = ((binary_df.get('ecDNA') == 1.0) | (binary_df.get('dm', 0) == 1.0)).astype(float)
                        binary_df.drop(columns=['ecDNA', 'dm'], errors='ignore', inplace=True)
                    sv_grouped_tag = ''
                
                # Variable to use the whole cohort or only pateints with top 50% entropy 
                full_cohort = True
                if full_cohort != True:
                    high_entropy = 'at_or_above_median_entropy_patients'
                else:
                    high_entropy = None

                # Get TFx column from master sheet
                tfx_df = get_tfx_column(output_path, input_map["pluvicto_master_sheet"][0], 'Sample_ID', "TFx_C1")

                # Threshold TFx tp ≥ 10% TFx if tfx_thresholding is True
                tfx_thresholding = False
                tfx_thresholding_tag = ''
                if tfx_thresholding:
                    tfx_df = tfx_df[tfx_df['TFx_C1'] >= 0.10].reset_index(drop=True)
                    tfx_thresholding_tag = "_tfx_10_perct_threshold"
                else:
                    tfx_df = tfx_df.dropna(axis=0)
                    
                tfx_suffix = f" - {tfx_thresholding_tag.replace('_', ' ').title()}" if tfx_thresholding_tag else ""

                # ----------------------------------
                # Prep DataFrame for Analyses
                # ----------------------------------
                # Filter binary_df
                binary_df = filter_binary_df(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "C1", f'chr{chromosome}', binary_df, tfx_df, high_entropy)
                # Get entropy column and merge it with binary_df
                entropy_column = get_entropy_column(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", f'chr{chromosome}', "C1")
                binary_df = pd.merge(binary_df, entropy_column, left_on='patient_id', right_on='patient_id', how='left')

                temp_df = pd.merge(binary_df.copy(), tfx_df, left_on='patient_id', right_on='Sample_ID', how='inner').drop(columns={'Sample_ID'})
                print(temp_df.sort_values(f'chr{chromosome}', ascending=False).to_string())

                # Save binary_df only on whole cohort
                if full_cohort and not tfx_thresholding:
                    binary_df_save = binary_df.sort_values(f'chr{chromosome}', ascending=False)
                    # binary_df_save = binary_df.drop(columns={'chr8'})
                    csv_dir_path = f'{output_path}/data-tables/complex_sv_binary_tables/{csv_path}'
                    os.makedirs(csv_dir_path, exist_ok=True)
                    csv_file_path = os.path.join(csv_dir_path, f'chr{chromosome}_complex_svs{sv_grouped_tag}_binary_table.csv')
                    binary_df_save.to_csv(csv_file_path, index=False)

                # Remove rows with NA
                binary_df.dropna(axis=0, how='any', inplace=True)

                # ----------------------------------
                # Visualize Stack Barplot of each SV with respect to entropy
                # ----------------------------------
                
                complex_sv_and_entropy_stacked_histogram_per_sv(
                    output_path=output_path, 
                    sub_dir=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process8_outputs/{csv_path.replace('_per_chr_table', '')}/complex_svs_vs_entropy_barplots", 
                    df=binary_df, 
                    entropy_column=f'chr{chromosome}',
                    title="Stacked Barplot of Complex SVs",
                    file_name=f'complex_svs{sv_grouped_tag}_stacked_barplot{all_sv_tag}',
                    ylabel='Number of Patients',
                    xlabel='Structural Variant'
                )

                binary_df_stripped = binary_df[[col for col in binary_df.columns if col not in ['entropy_group', f'chr{chromosome}']]]
                horizontal_bar_plot(
                    output_path=output_path, 
                    sub_dir=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process8_outputs/{csv_path.replace('_per_chr_table', '')}/complex_svs_vs_entropy_barplots", 
                    df=binary_df_stripped, 
                    title="Barplot of Complex SVs Counts", 
                    file_name=f'counts_of_complex_svs{sv_grouped_tag}_horizontal_barplot{all_sv_tag}', 
                    ylabel='Structural Variant', 
                    xlabel='Number of Events'
                )

                proportional_stacked_histogram_per_sv(
                    output_path=f"{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process8_outputs/{csv_path.replace('_per_chr_table', '')}/complex_svs_vs_entropy_barplots", 
                    df=binary_df, 
                    id_column='patient_id',
                    entropy_column=f'chr{chromosome}',
                    title="Structural variant frequency stratified by entropy",
                    file_name=f'complex_svs{sv_grouped_tag}_proportional_stacked_barplot{all_sv_tag}',
                    ylabel='Proportion of patients',
                    xlabel='Structural Variant'
                )

                binary_df = binary_df.drop(columns=['entropy_group'])

                # ----------------------------------
                # SV and Entropy Visualization
                # ----------------------------------

                # Add a column to test the correlation between presence of any SVs and entropy
                cols_to_check = binary_df.columns.difference(['patient_id', f'chr{chromosome}'])
                binary_df['any_sv'] = (binary_df[cols_to_check].any(axis=1)).astype(int)

                # Tag for plots title and filename
                if full_cohort == True:
                    high_entropy = 'whole_cohort'

                filter_step = None # either 'high_vs_low_entropy' or None - selects all patients or patients that have high entropy, no sv and low entropy, with sv

                # SV and high entropy correlation
                sv_entorpy_correlation(
                    output_path, 
                    f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process8_outputs/{csv_path.replace('_per_chr_table', '')}/complex_svs{sv_grouped_tag}_vs_entropy_boxplots", 
                    binary_df.copy(), 
                    f'chr{chromosome}', 
                    high_entropy, 
                    tfx_thresholding_tag, 
                    filter_step
                )

                # Drop any sv presence column
                binary_df.drop(columns=['any_sv'], inplace=True)
                binary_df.reset_index(drop=True)

                # --- Simple scatter plot of entropy vs number of complex svs per patient ---
                cols_to_sum = binary_df.columns.difference(["patient_id", f'chr{chromosome}'])
                binary_df['sv_count'] = binary_df[cols_to_sum].sum(axis=1)

                # Whole cohort
                scatter_plot_two_columns_one_dataframe(
                    df=binary_df, 
                    col_1_name='sv_count', 
                    col_2_name=f'chr{chromosome}', 
                    title='Chr8 Entropy vs Number of Complex SVs Per Patient', 
                    xlabel='Number of Complex SVs', 
                    ylabel='Chr8 Entropy', 
                    out_dir=f"{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process8_outputs/{csv_path.replace('_per_chr_table', '')}/scatter_plots", 
                    file_name=f"chr{chromosome}_entropy_vs_sv_count_per_patient_scatter_plot{all_sv_tag}",
                    jitter=True
                )

                # Patients with at least one complex sv
                binary_df_filtered = binary_df[binary_df["sv_count"] >= 1]
                scatter_plot_two_columns_one_dataframe(
                    df=binary_df_filtered, 
                    col_1_name='sv_count', 
                    col_2_name=f'chr{chromosome}', 
                    title='Chr8 Entropy vs Number of Complex SVs Per Patient (at least one SV)', 
                    xlabel='Number of Complex SVs', 
                    ylabel='Chr8 Entropy', 
                    out_dir=f"{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process8_outputs/{csv_path.replace('_per_chr_table', '')}/scatter_plots", 
                    file_name=f"chr{chromosome}_entropy_vs_sv_count_filtered_per_patient_scatter_plot{all_sv_tag}",
                    jitter=True
                )

                # --- bar plot of entropy vs number of complex svs per patient --- 
                bar_plot_groups=None

                box_plot_two_columns_one_dataframe(
                    df=binary_df, 
                    col_1_name='sv_count', 
                    col_2_name=f'chr{chromosome}', 
                    title='Chr8 Entropy vs Number of Complex SVs Per Patient', 
                    xlabel='Number of Complex SVs', 
                    ylabel='Chr8 Entropy', 
                    out_dir=f"{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process8_outputs/{csv_path.replace('_per_chr_table', '')}/box_plots", 
                    file_name=f"chr{chromosome}_entropy_vs_sv_count_per_patient_box_plot{all_sv_tag}",
                    jitter=True,
                    compare_groups=bar_plot_groups
                )

                # Drop any sv_count column if group_svs (when grouped group_svs columm becomes like any_sv column)
                binary_df_temp = binary_df.drop(columns=['sv_count'])
                binary_df_temp.reset_index(drop=True)

                # -----------------------------------------
                # Logistic regression to correlate with svs
                # -----------------------------------------
                model_results = logistic_regression_with_covariates(binary_df_temp.copy(), 'patient_id', tfx_df, 'chr8', 'TFx_C1', 'Sample_ID', None)

                odds_ratio_forest_plot(
                    output_path, 
                    f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process8_outputs/{csv_path.replace('_per_chr_table', '')}/complex_svs_vs_entropy_forest_plots/logistic_regression", 
                    model_results['coefficients'], 
                    "entropy_response_odds_ratio", 
                    "entropy_response_adj_pvalue", 
                    None, 
                    None, 
                    None, 
                    f"Odds Ratios for SVs per Unit Increase in Chr8 Entorpy (controlled for TFx) - ({high_entropy.replace('_', '').title()}) {tfx_suffix}", 
                    f"complex_sv{sv_grouped_tag}_and_entropy_response_odds_ratio_{high_entropy}{tfx_thresholding_tag}{all_sv_tag}_forest_plot_FDR_corrected"
                )

                # ----------------------------------
                # SV Survival Analysis
                # ----------------------------------
                # Drop entropy column
                binary_df = binary_df.drop(columns='chr8')
                # Drop any sv_count column if group_svs (when grouped group_svs columm becomes like any_sv column)
                if group_svs:
                    binary_df.drop(columns=['sv_count'], inplace=True)
                    binary_df.reset_index(drop=True)

                create_sv_kaplan_meier_plots(output_path, input_map["pluvicto_master_sheet"][0], csv_path, binary_df, sequencing_type, high_entropy, tfx_thresholding_tag, sv_grouped_tag)

    # ==============
    # Process 9:
    # ==============
    if args.process9 == "True":

        # Dictionary of other projects (a list of: outcome data)
        projects = {
            'projectA': ['path/to/projectA/clinical/data'],
            'projectB': ['path/to/projectB/clinical/data'],
        }

        # List of csv files to process
        entropy_csv_dict = {
            'base_entropy_hn_normalized': 'base_entropy_hn_normalized_per_chr_table',
        }

        tfx_threshold = True

        # Booleans to run different sections of the code
        psa50_response_logistic_regression = False
        treatment_comopletion_logistic_regression = False
        met_progression_logistic_regression = False
        cohort_entropy_distribution_plots = False
        docetaxel_and_pluvicto_hazard_ratio = True

        # Set CSV file variable to use in process
        for idx, (label, csv_path) in enumerate(entropy_csv_dict.items()):

            # =========================================
            # Logistic regression for psa50 and entropy
            # =========================================

            if psa50_response_logistic_regression:
                # --- Pluvicto model ---
                # Get entropy column
                entropy_column = get_entropy_column(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "chr8", "C1")

                # Read in pluvicto survival data
                master_sheet_df = pd.read_csv(input_map["pluvicto_master_sheet"][0], index_col=0)
                
                # Collect TFx column
                tfx_column = master_sheet_df[['Sample_ID','TFx_C1']]

                # Merge and filter pluvicto master sheet and entropy
                merged_df = pd.merge(entropy_column, master_sheet_df, left_on='patient_id', right_on='Sample_ID', how='inner').drop(columns='TFx_C1')
                filtered_df = merged_df[['patient_id', 'PSA50', 'chr8']]
                filtered_df = filtered_df.rename(columns={'PSA50': 'psa_50_response_pluvicto'})

                # Model logistic regression
                model_results_pluvicto = logistic_regression_with_covariates(filtered_df, 'patient_id', tfx_column, 'chr8', 'TFx_C1', 'Sample_ID', None)

                # --- docetaxel model ---
                project = 'docetaxel'

                # Read in outcome data
                outcome_data = pd.read_csv(projects[project][0], index_col=0, encoding="latin1")

                psa_50_col = 'psa_50_response'

                # Get entropy column and merge it with outcome_data
                entropy_column_temp = get_entropy_column(output_path, f'data-tables/other-datasets/{project}/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv', "chr8", None)
                merged_df_temp = pd.merge(outcome_data, entropy_column_temp, left_on='Sample label', right_on='patient_id', how='inner')

                # Get the TFx column from outcome_data
                tfx_column_temp = merged_df_temp[['patient_id', 'tumor_fraction']]

                # Subset merged_df to only columns you want and format
                merged_df_temp.rename(columns={'PSA decline ³ 50% within 16 weeks after docetaxel start; Y = YES; N = NO': psa_50_col}, inplace=True)
                filtered_df_temp = merged_df_temp[['patient_id', psa_50_col, 'chr8']]

                # To supress warnings
                filtered_df_temp = filtered_df_temp.copy()

                filtered_df_temp[psa_50_col] = (filtered_df_temp[psa_50_col] == 'Y').astype(int)

                # Rename column
                filtered_df_temp = filtered_df_temp.rename(columns={psa_50_col: f'{psa_50_col}_{project}'})

                # To correct for direction of the logistic regression you need to drop one of the columns and flip 
                model_results_project = logistic_regression_with_covariates(filtered_df_temp, 'patient_id', tfx_column_temp, 'chr8', 'tumor_fraction', 'patient_id', None)

                # Merge models
                merged_model_coefficients = pd.concat(
                    [model_results_pluvicto['coefficients'],
                    model_results_project['coefficients']],
                    axis=0
                )

                # Plot merged models
                odds_ratio_forest_plot(
                    output_path, 
                    f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process9_outputs/pluvicto_vs_docetaxel/{label}/logistic-regression", 
                    merged_model_coefficients, 
                    "entropy_response_odds_ratio", 
                    "entropy_response_adj_pvalue", 
                    None, 
                    None, 
                    None, 
                    f"Odds Ratios for PSA 50 Response Associated with Chr8 Entorpy (controlled for TFx) in Pluvicto and Docetaxel", 
                    f"pluvicto_vs_docetaxel_survival_associated_with_chr8_entropy_odds_ratio_forest_plot_FDR_corrected"
                )

            # ========================================================
            # Logistic regression for treatment completion and entropy
            # ========================================================
            if treatment_comopletion_logistic_regression:
                # --- Pluvicto model ---
                # Get entropy column
                entropy_column = get_entropy_column(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "chr8", "C1")

                # Read in pluvicto survival data
                master_sheet_df = pd.read_csv(input_map["pluvicto_master_sheet"][0], index_col=0)
                
                # Collect TFx column
                tfx_column = master_sheet_df[['Sample_ID','TFx_C1']]

                # Merge and filter pluvicto master sheet and entropy
                merged_df = pd.merge(entropy_column, master_sheet_df, left_on='patient_id', right_on='Sample_ID', how='inner').drop(columns='TFx_C1')

                # Create binary cycles column
                merged_df['Completed_Tx_binary'] = (merged_df['T_cycles'] == 6).astype(int)
                filtered_df = merged_df[['patient_id', 'Completed_Tx_binary', 'chr8']]
                filtered_df = filtered_df.rename(columns={'Completed_Tx_binary': 'Completed_Tx_binary_pluvicto'})

                # Model logistic regression
                model_results_pluvicto = logistic_regression_with_covariates(filtered_df, 'patient_id', tfx_column, 'chr8', 'TFx_C1', 'Sample_ID', None)

                # --- Radium223 model ---
                project = 'radium223'

                # Read in outcome data
                outcome_data = pd.read_csv(projects[project][0])

                completed_Tx_col = 'Completed_Tx_binary'
                tfx_col = 'preTFx'

                # Get entropy column and merge it with outcome_data
                entropy_column_temp = get_entropy_column(output_path, f'data-tables/other-datasets/{project}/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv', "chr8", 'pre')
                
                # extract patient_id from the left column
                merged_df_temp = outcome_data.copy()

                # Filter to pre treatment
                merged_df_temp = merged_df_temp[merged_df_temp['Pre_treatment_specimen'].str.contains('pre')]

                # Get the first part of the patient id
                merged_df_temp["patient_id"] = merged_df_temp['Pre_treatment_specimen'].str.split("_").str[0]

                # Drop duplicates 
                merged_df_temp = merged_df_temp.drop_duplicates(subset="patient_id")
                entropy_column_temp = entropy_column_temp.drop_duplicates(subset="patient_id")

                # now do a normal merge
                merged_df_temp = pd.merge(
                    merged_df_temp,
                    entropy_column_temp,
                    on="patient_id",
                    how="inner"
                )

                # Get the TFx column from outcome_data
                tfx_column_temp = merged_df_temp[['patient_id', tfx_col]]

                # Subset merged_df to only columns you want
                filtered_df_temp = merged_df_temp[['patient_id', completed_Tx_col, 'chr8']]

                # To supress warnings
                filtered_df_temp = filtered_df_temp.copy()

                # Rename column
                filtered_df_temp = filtered_df_temp.rename(columns={completed_Tx_col: f'{completed_Tx_col}_{project}'})

                # To correct for direction of the logistic regression you need to drop one of the columns and flip 
                model_results_project = logistic_regression_with_covariates(filtered_df_temp, 'patient_id', tfx_column_temp, 'chr8', tfx_col, 'patient_id', None)


                # Merge models
                merged_model_coefficients = pd.concat(
                    [model_results_pluvicto['coefficients'],
                    model_results_project['coefficients']],
                    axis=0
                )

                # Plot merged models
                odds_ratio_forest_plot(
                    output_path, 
                    f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process9_outputs/pluvicto_vs_radium223/{label}/logistic-regression", 
                    merged_model_coefficients, 
                    "entropy_response_odds_ratio", 
                    "entropy_response_adj_pvalue", 
                    None, 
                    None, 
                    None, 
                    f"Odds Ratios for Treatment Completion Associated with Chr8 Entorpy (controlled for TFx) in Pluvicto and Docetaxel", 
                    f"pluvicto_vs_docetaxel_Tx_completion_associated_with_chr8_entropy_odds_ratio_forest_plot_FDR_corrected"
                )

            # ==========================================================
            # Logistic regression for metastasis progression and entropy (Radium vs Docetaxel)
            # ==========================================================
            if met_progression_logistic_regression:

                # --- Other Projects model ---
                model_results_project = None
                for project in projects:
                    # Read in outcome data
                    if project == 'docetaxel':
                        outcome_data = pd.read_csv(projects[project][0], encoding="latin1")
                        met_cols = ['LN involvement', 'Visceral involvement', 'Liver involvement']
                        tfx_col = 'tumor_fraction'
                        id_col = 'Sample label'
                        cycle_filter = None
                    else:
                        outcome_data = pd.read_csv(projects[project][0])
                        met_cols = ['Bone_Prog', 'LN_Prog', 'Visc_Prog', 'Liver_Prog']
                        tfx_col = 'preTFx'
                        id_col = 'Pre_treatment_specimen'
                        cycle_filter = 'pre'

                    if not all(col in outcome_data.columns for col in met_cols):
                        continue

                    # Create new metastasis column
                    met_bool = outcome_data[met_cols].eq("Y")
                    outcome_data["met_prog"] = met_bool.any(axis=1)
                    outcome_data["met_prog"] = outcome_data["met_prog"].where(
                        ~outcome_data[met_cols].isna().all(axis=1),
                        pd.NA
                    ).astype("Int64")

                    # Get entropy column and merge it with outcome_data
                    entropy_column_temp = get_entropy_column(output_path, f'data-tables/other-datasets/{project}/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv', "chr8", cycle_filter)
                    
                    # extract patient_id from the left column
                    merged_df_temp = outcome_data.copy()

                    if project == 'radium223':
                        # Filter to pre treatment
                        merged_df_temp = merged_df_temp[merged_df_temp[id_col].str.contains('pre')]

                        # Get the first part of the patient id
                        merged_df_temp["patient_id"] = merged_df_temp[id_col].str.split("_").str[0]

                    else:
                        merged_df_temp["patient_id"] = merged_df_temp[id_col]


                    # Drop duplicates 
                    merged_df_temp = merged_df_temp.drop_duplicates(subset='patient_id')
                    entropy_column_temp = entropy_column_temp.drop_duplicates(subset="patient_id")

                    # now do a normal merge
                    merged_df_temp = pd.merge(
                        merged_df_temp,
                        entropy_column_temp,
                        on="patient_id",
                        how="inner"
                    )

                    # Get the TFx column from outcome_data
                    tfx_column_temp = merged_df_temp[['patient_id', tfx_col]]

                    # Subset merged_df to only columns you want
                    filtered_df_temp = merged_df_temp[['patient_id', 'met_prog', 'chr8']]

                    # To supress warnings
                    filtered_df_temp = filtered_df_temp.copy()
    
                    # Rename column
                    filtered_df_temp = filtered_df_temp.rename(columns={'met_prog': f'met_prog_{project}'})

                    # To correct for direction of the logistic regression you need to drop one of the columns and flip 
                    model_results_project_temp = logistic_regression_with_covariates(filtered_df_temp, 'patient_id', tfx_column_temp, 'chr8', tfx_col, 'patient_id', None)


                    # Merge with former model results
                    if model_results_project is not None:
                        model_results_project = pd.concat(
                            [model_results_project['coefficients'],
                            model_results_project_temp['coefficients']],
                            axis=0
                        )
                    else:
                        model_results_project = model_results_project_temp

                # Plot merged models
                odds_ratio_forest_plot(
                    output_path, 
                    f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process9_outputs/radium223_vs_docetaxel/{label}/logistic-regression", 
                    model_results_project, 
                    "entropy_response_odds_ratio", 
                    "entropy_response_adj_pvalue", 
                    None, 
                    None, 
                    None, 
                    f"Odds Ratios for Metatasis Progression Associated with Chr8 Entorpy (controlled for TFx) in Radium23 and Docetaxel", 
                    f"radium223_vs_docetaxel_metatasis_progression_associated_with_chr8_entropy_odds_ratio_forest_plot_FDR_corrected"
                )


            # ==========================================================================
            # Creates cohort wide entropy distribution plots (pluvicto vs other dataset)
            # ==========================================================================

            if cohort_entropy_distribution_plots:
                metric_to_use = "Corrected_Copy_Number"
                
                # ***************************************************
                # Create violin plots for entropy(dataset vs dataset)
                # ***************************************************
                df_1, df_1_tag = prepare_distribution_data(
                    csv_path=f'{output_path}/data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_path}.csv', 
                    id_col='patient_id',
                    list_of_cols=['chr8'],
                    cycle_col='cycle',
                    cycle_filter="C1",
                    col_tag='pluvicto'
                )

                df_2, df_2_tag = prepare_distribution_data(
                    csv_path=f'{output_path}/data-tables/other-datasets/radium223/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv',
                    id_col='patient_id',
                    list_of_cols=['chr8'],
                    cycle_col=None,
                    cycle_filter=None,
                    col_tag='radium223'
                )

                df_3, df_3_tag = prepare_distribution_data(
                    csv_path=f'{output_path}/data-tables/other-datasets/docetaxel/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv',
                    id_col='patient_id',
                    list_of_cols=['chr8'],
                    cycle_col=None,
                    cycle_filter=None,
                    col_tag='docetaxel'
                )

                # Create violin plot
                compare_column_violin_between_dfs(
                    df_1=df_1, 
                    col_name_1=f'chr8_{df_1_tag}', 
                    df_2=df_2, 
                    col_name_2=f'chr8_{df_2_tag}', 
                    df_3=df_3, 
                    col_name_3=f'chr8_{df_3_tag}', 
                    title=f'{df_1_tag.title()} vs {df_2_tag.title()}: Chr8 Entropy', 
                    xlabel=f'Datasets', 
                    ylabel=f'Chr8 Entropy', 
                    out_dir=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process9_outputs/{df_1_tag}_vs_{df_2_tag}_vs_{df_3_tag}/{label}/distribution_plots', 
                    file_name=f'{df_1_tag}_vs_{df_2_tag}_vs_{df_3_tag}_chr8_entropy_volin_plot'
                )

                # Create histograms for entropy
                compare_column_density_distribution_between_dfs(
                    df_1=df_1, 
                    col_name_1=f'chr8_{df_1_tag}', 
                    df_2=df_2, 
                    col_name_2=f'chr8_{df_2_tag}', 
                    df_3=df_3, 
                    col_name_3=f'chr8_{df_3_tag}', 
                    title=f'{df_1_tag.title()} vs {df_2_tag.title()}: Chr8 Entropy', 
                    xlabel=f'Chr8 Entropy', 
                    ylabel=f'Density', 
                    out_dir=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process9_outputs/{df_1_tag}_vs_{df_2_tag}_vs_{df_3_tag}/{label}/distribution_plots', 
                    file_name=f'{df_1_tag}_vs_{df_2_tag}_vs_{df_3_tag}_chr8_entropy_density_distribution_plot'
                )

            # ==========================================================================
            # Cox Proportional Hazard Model (pluvicto vs docetaxel)
            # ==========================================================================
            if docetaxel_and_pluvicto_hazard_ratio:
                # --- Pluvicto model ---
                # Get entropy column
                entropy_column = get_entropy_column(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "chr8", "C1")

                # Read in pluvicto survival data
                master_sheet_df = pd.read_csv(input_map["pluvicto_master_sheet"][0], index_col=0)

                # If deep filter out samples below 10% TFx
                if 'deep' in sequencing_type:
                    master_sheet_df = master_sheet_df[master_sheet_df['TFx_C1'] > 0.20]
                else: 
                    master_sheet_df = master_sheet_df[master_sheet_df['TFx_C1'] > 0.10]
                
                # Collect TFx column
                tfx_column = master_sheet_df[['Sample_ID','TFx_C1', 'Death', 'survival_days']]

                # Merge and filter pluvicto master sheet and entropy
                merged_df = pd.merge(entropy_column, master_sheet_df, left_on='patient_id', right_on='Sample_ID', how='inner')
                filtered_df = merged_df[['patient_id', 'TFx_C1', 'Death', 'survival_days', 'chr8']]

                pluvicto_cox_model = cox_proportional_hazards_model(
                    cox_data=filtered_df, 
                    time_to_event_column='survival_days', 
                    event_col='Death', 
                    formula='chr8 + TFx_C1', 
                    chr_to_pick='chr8'
                )

                # --- docetaxel model ---
                project = 'docetaxel'
                
                # Read in outcome data
                outcome_data = pd.read_csv(projects[project][0], encoding="latin1")

                # If deep filter out samples below 10% TFx
                if 'deep' in sequencing_type:
                    outcome_data = outcome_data[outcome_data['tumor_fraction'] > 0.20]
                else: 
                    outcome_data = outcome_data[outcome_data['tumor_fraction'] > 0.10]

                survival_col_name = 'survival_days'

                # Get entropy column and merge it with outcome_data
                entropy_column_temp = get_entropy_column(output_path, f'data-tables/other-datasets/{project}/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv', "chr8", None)
                merged_df_temp = pd.merge(outcome_data, entropy_column_temp, left_on='ï»¿Docetaxel Sample label', right_on='patient_id', how='inner')

                # Create death column
                merged_df_temp['Death'] = merged_df_temp[
                    'Still alive; Y = YES; N = NO; N/A = Not available'
                ].map({
                    'N': 1,
                    'Y': 0,
                    'N/A': np.nan
                })

                # Get the TFx column from outcome_data
                filtered_df_temp = merged_df_temp[['patient_id', 'Days from sample draw to death', 'Death', 'chr8', 'tumor_fraction']]

                # Rename columns
                filtered_df_temp = filtered_df_temp.rename(columns={'Days from sample draw to death': 'survival_days'})
                
                # Remove NAs
                filtered_df_temp.dropna(inplace=True)
                filtered_df_temp = filtered_df_temp.reset_index(drop=True)

                # To supress warnings
                filtered_df_temp = filtered_df_temp.copy()

                docetaxel_cox_model = cox_proportional_hazards_model(
                    cox_data=filtered_df_temp, 
                    time_to_event_column='survival_days', 
                    event_col='Death', 
                    formula='chr8 + tumor_fraction', 
                    chr_to_pick='chr8'
                )

                # --- Combine HRs onto one plot ---

                # Extract only chr8 rows
                # pluvicto_chr8 = pluvicto_cox_model.summary.loc[['chr8']].copy()
                # docetaxel_chr8 = docetaxel_cox_model.summary.loc[['chr8']].copy()
                pluvicto_chr8 = pluvicto_cox_model.summary.copy()
                docetaxel_chr8 = docetaxel_cox_model.summary.copy()
                

                # print(pluvicto_chr8)
                # print(docetaxel_chr8)
                # sys.exit(0)

                # # Rename index to model names
                # pluvicto_chr8.index = ['Pluvicto_chr8_entropy']
                # docetaxel_chr8.index = ['Docetaxel_chr8_entropy']

                pluvicto_chr8.index = ['Pluvicto_chr8_entropy', 'Pluvicto_TFx']
                docetaxel_chr8.index = ['Docetaxel_chr8_entropy', 'Docetaxel_TFx']

                # Combine
                forest_df = pd.concat([pluvicto_chr8, docetaxel_chr8])

                # FDR correction
                pvalues = forest_df['p'].values
                _, pvalues_corrected, _, _ = multipletests(pvalues, method='fdr_bh')

                # Add corrected pvalues back to dataframe
                forest_df['pvalues_corrected'] = pvalues_corrected

                # Filter froest_df
                forest_df = forest_df.loc[forest_df.index.isin(['Pluvicto_chr8_entropy', 'Docetaxel_chr8_entropy'])]

                file_path = f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process9_outputs/pluvicto_vs_docetaxel/cox_hardard_ratios/entropy_plots/{csv_path}'
                if tfx_threshold:
                    file_path = file_path + '/tfx_threshold'
                file_name=f'forest_plot_combined_entropy'

                out_path, summary = cox_forest_plot_pub(
                    forest_df, 
                    file_path, 
                    'pvalues_corrected',
                    filename=file_name,
                    title='Cox Proportional Hazards Model (Ctrl for TFx - FDR Corrected)'
                )

    # ==============
    # Process 10:
    # ==============
    if args.process10 == "True":
        # List of csv files to process
        entropy_csv_list = [
            "base_entropy_hn_normalized_per_chr_table",
        ]

        sv_to_filter_out = ['dup', 'del', 'inv', 'tra', 'invdup']
        all_sv_tag = ''
        # sv_to_filter_out = []
        # all_sv_tag = '_all_svs'

        group_svs = True

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:

            complex_svs_and_entropy_df = clean_and_fix_complex_sv_csv(input_map["complex_sv_calls_and_chr8_entropy"][0], "sample_id", "caller", "location", "SV_type", "chr8_entropy", sv_to_filter_out)
            # complex_svs_and_new_entropy_df = replace_chr8_entropy_column(output_path, complex_svs_and_entropy_df, f"{csv_path}.csv", "Corrected_Copy_Number", sequencing_type, "chr8", "C1")

            # Filteres dataframe down to only structural variants found in chr8
            complex_svs_and_new_entropy_df_filtered = filter_to_location_in_complex_sv_df(
                df=complex_svs_and_entropy_df, 
                location_column='location', 
                chr_filter='8'
            )

            binary_df = create_binary_dataframe(complex_svs_and_new_entropy_df_filtered)
            binary_df = prepare_sv_df_for_survival_analysis(binary_df, input_map["pluvicto_master_sheet"][0])

            # Create new groupings based off biology and drop all other column (except patient_id)
            if group_svs:
                binary_df_columns = binary_df.columns
                binary_df['ecDNA_dm'] = ((binary_df['ecDNA'] == 1.0) | (binary_df['dm'] == 1.0)).astype(float)
                binary_df['bfb_tyfonas_pyrgo'] = ((binary_df['bfb'] == 1.0) | (binary_df['tyfonas'] == 1.0) | (binary_df['pyrgo'] == 1.0)).astype(float)
                binary_df['chromoplexy_Complex_non_cyclic'] = ((binary_df['chromoplexy'] == 1.0) | (binary_df['Complex_non_cyclic'] == 1.0)).astype(float)
                binary_df['Linear_tic'] = ((binary_df['Linear'] == 1.0) | (binary_df['tic'] == 1.0)).astype(float)
                binary_df = binary_df.drop(columns=[col for col in binary_df_columns if col != "patient_id"])
                sv_grouped_tag = '_grouped'
            else:
                binary_df['ecDNA_dm'] = ((binary_df['ecDNA'] == 1.0) | (binary_df['dm'] == 1.0)).astype(float)
                binary_df.drop(columns=['ecDNA', 'dm'], inplace=True)
                sv_grouped_tag = ''
            
            # Variable to use the whole cohort or only pateints with top 50% entropy 
            full_cohort = True
            if full_cohort != True:
                high_entropy = '_at_or_above_median_entropy_patients'
            else:
                high_entropy = ''

            # Get TFx column from master sheet
            tfx_df = get_tfx_column(output_path, input_map["pluvicto_master_sheet"][0], 'Sample_ID', "TFx_C1")

            # Threshold TFx tp ≥ 10% TFx if tfx_thresholding is True
            tfx_thresholding = False
            tfx_thresholding_tag = ''
            if tfx_thresholding:
                tfx_df = tfx_df[tfx_df['TFx_C1'] >= 0.10].reset_index(drop=True)
                tfx_thresholding_tag = "_tfx_10_perct_threshold"
            else:
                tfx_df = tfx_df.dropna(axis=0)
                
            tfx_suffix = f" - {tfx_thresholding_tag.replace('_', ' ').title()}" if tfx_thresholding_tag else ""

            # ----------------------------------
            # Prep DataFrame for Analyses
            # ----------------------------------
            # Filter binary_df
            binary_df = filter_binary_df(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "C1", "chr8", binary_df, tfx_df, high_entropy)

            # Get entropy column and merge it with binary_df
            entropy_column = get_entropy_column(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "chr8", "C1")
            binary_df = pd.merge(binary_df, entropy_column, left_on='patient_id', right_on='patient_id', how='left')

            # ----------------------------------
            # Add some extra columns
            # ----------------------------------

            # Add a column to test the correlation between presence of any SVs and entropy
            cols_to_check = binary_df.columns.difference(['patient_id', 'chr8'])
            binary_df['any_sv'] = (binary_df[cols_to_check].any(axis=1)).astype(int)

            # # --- Simple scatter plot of entropy vs number of complex svs per patient ---
            cols_to_sum = binary_df.columns.difference(["patient_id", "chr8", 'any_sv'])
            binary_df['sv_count'] = binary_df[cols_to_sum].sum(axis=1)

            # -----------------------------------------
            # Logistic regression to correlate with svs
            # -----------------------------------------
            model_results = linear_regression_with_covariates(binary_df.copy(), 'patient_id', tfx_df, 'chr8', 'TFx_C1', 'Sample_ID', None)
            
            # Take out TFx rows 
            model_df = model_results['coefficients'][~model_results['coefficients']['feature'].str.contains('TFx')].copy().reset_index(drop=True)

            # Create a forest plot with each structural variant 
            coefficients_forest_plot(
                output_path=f"{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process10_outputs/complex_svs_vs_entropy_forest_plots/{csv_path}/linear_regression",
                df=model_df,
                feature_column='feature',
                coefficent_column="coefficient",
                ci_lower_column='ci_lower',
                ci_upper_column='ci_upper',
                adj_pvalue_column="adj_pvalue",
                x_axis_title='Coefficient value (β)',
                y_axis_title='SV',
                plot_title=f"Change in Chr8 Entropy Per Unit Increase in SV (controlled for TFx) -({high_entropy.replace('_', '').title()}) {tfx_suffix}",
                file_name=f"complex_sv{sv_grouped_tag}_and_entropy_response_coefficients{high_entropy}{tfx_thresholding_tag}{all_sv_tag}_forest_plot_FDR_corrected",
            )

    # ==============
    # Process 11:
    # ==============
    if args.process11 == "True":

        # List of csv files to process
        entropy_csv_list = [
            'base_entropy_hn_normalized_per_chr_table',
        ]

        normalize_by_ploidy = True

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:
                
            # Process Gene matrix
            filtered_gene_list = filter_gene_annotation_list(input_map["gene_annotation_list"][0], "Chr", [8], "Gene")
            gene_matrix = final_gene_matrix(filtered_gene_list, "Gene", input_map["gene_matrix_from_titan"], "Sample", "C1", input_map["purity_ploidy_from_titan"], "Ploidy")

            # Process Entropy table 
            entropy_column = get_entropy_column(output_path, f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", "chr8", "C1")

            # Set Parameters for next analyses
            cn_threshold = 1
            responder_grouping = "no_grouping" # 'median' or 'quartiles_extreme'
            direction = "gain" # 'gain' or 'loss'
            gene_matrix_id_column = "Sample"
            pluvicto_id_column = "Sample_ID"
            tfx_thresholding = True

            # ----------------------------------
            # Prep DataFrame for Analyses
            # ----------------------------------
            # Create gene matrix/dataframe
            gene_matrix_binary = threshold_gene_matrix(gene_matrix, gene_matrix_id_column, "Ploidy", entropy_column, "patient_id", "chr8",  cn_threshold, direction, normalize_by_ploidy)
            tfx_df = get_tfx_column(output_path, input_map["pluvicto_master_sheet"][0], pluvicto_id_column, "TFx_C1")

            # Set 10% TFx threshold if tfx_thresholding is True
            tfx_thresholding_tag = ''
            tfx_suffix = ''
            if tfx_thresholding:
                tfx_df = tfx_df[tfx_df['TFx_C1'] >= 0.10].reset_index(drop=True)
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
            model_results = linear_regression_with_covariates(
                gene_matrix_binary=gene_matrix_binary,
                gene_matrix_id_column=gene_matrix_id_column,
                tfx_df=tfx_df,
                chr_column='chr8',
                tfx_column='TFx_C1',
                tfx_id_column=pluvicto_id_column,
                responder_grouping=responder_grouping
            )

            # --------------------------------------
            # Plot metrics from logistic regression
            # --------------------------------------
            # Plot coefficients on bar plots
            coefficient_barplot(
                output_path=output_path,
                sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process11_outputs/gene_matrix_and_entropy/{csv_path}/linear_regression",
                model_results=model_results['coefficients'],
                coefficent_column="coefficient",
                adj_pvalue_column="adj_pvalue",
                responder_grouping=responder_grouping,
                cn_threshold=cn_threshold,
                direction=direction,
                normalize_by_ploidy=normalize_by_ploidy
            )

            # Take out TFx rows 
            model_df = model_results['coefficients'][~model_results['coefficients']['feature'].str.contains('TFx')].copy().reset_index(drop=True)

            # Create a forest plot with each gene 
            coefficients_forest_plot(
                output_path=f"{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process11_outputs/gene_matrix_and_entropy/{csv_path}/linear_regression",
                df=model_df,
                feature_column='feature',
                coefficent_column="coefficient",
                ci_lower_column='ci_lower',
                ci_upper_column='ci_upper',
                adj_pvalue_column="adj_pvalue",
                x_axis_title='Coefficients',
                y_axis_title='Gene',
                plot_title=f'Change in Chr8 Entropy Per Unit Increase in Gene Copy Number (controlled for TFx) - CN {direction.title()} {tfx_suffix} {normalization_suffix}',
                file_name=f'gene_matrix_and_entropy_response_odds_ratio_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}{tfx_thresholding_tag}{normalization_tag}_{direction}_vs_other_forest_plot',
            )

            # Create a gene list from a prostate specific oncogene list
            gene_list = get_prostate_specific_gene_list(
                gene_list_path=input_map["prostate_specific_cancer_gene_list"][0],
                chr_column="chr8"
            )
            
            #Plot odds ratios on forest plots subsetted to gene list
            coefficients_forest_plot(
                output_path=f"{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process11_outputs/gene_matrix_and_entropy/{csv_path}/linear_regression",
                df=model_df,
                feature_column='feature',
                coefficent_column="coefficient",
                ci_lower_column='ci_lower',
                ci_upper_column='ci_upper',
                adj_pvalue_column="adj_pvalue",
                x_axis_title='Coefficients',
                y_axis_title='Gene',
                plot_title=f'Change in Chr8 Entropy Per Unit Increase in Gene Copy Number (controlled for TFx) - CN {direction.title()} {tfx_suffix} {normalization_suffix}',
                file_name=f"gene_matrix_and_entropy_response_odds_ratio_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}{tfx_thresholding_tag}{normalization_tag}_{direction}_vs_other_prostate_gene_list_forest_plot",
                gene_list=gene_list
            )

            # Create dot plots for GSEA pathways
            try:
                GSEA(
                    output_path=output_path,
                    sub_dir_path=f"entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process11_outputs/gene_matrix_and_entropy/{csv_path}/linear_regression",
                    model_results=model_results['coefficients'],
                    gsea_ranking_column="coefficient",
                    adj_pvalue_column="adj_pvalue",
                    gene_set="MSigDB_Hallmark_2020",
                    responder_grouping=responder_grouping,
                    cn_threshold=cn_threshold,
                    direction=direction
                )
            except Exception as e:
                print(f"GSEA could not run: {e}")
    
    # ==============
    # Process 12:
    # ==============
    if args.process12 == "True":
        # List of csv files to process
        entropy_csv_list = [
            "base_entropy_hn_cohort_normalized_per_chr_table",
        ]

        extreme_responders_categories_csv = '/path/to/proteus_deg_sample_classification.csv'

        cycles_to_use = ['TFx_C1', 'TFx_C2']#, 'TFx_C3', 'TFx_C4', 'TFx_C5', 'TFx_C6']
        tfx_cutoff = 0.10
        response_category = 'extreme_responders'

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:
            # Read entropy_df and drop columns
            entropy_df = pd.read_csv(f"{output_path}/data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv")
            entropy_df = entropy_df[['patient_id', 'cycle', 'chr8']]

            # Read outcome_data
            outcome_data = pd.read_csv(input_map["pluvicto_master_sheet"][0])

            # Remove patients less than TFx cutoff
            outcome_data_filtered = outcome_data[outcome_data[cycles_to_use[0]] >= tfx_cutoff]

            # Filter columns
            outcome_data_slice = outcome_data_filtered[['Sample_ID', *cycles_to_use]]

            # Convert outcome_data from wide to long
            outcome_data_T = pd.melt(
                frame=outcome_data_slice,
                id_vars=['Sample_ID'],
                value_vars=cycles_to_use,
                var_name='TFx_cycle',
                value_name='TFx'
            )
            # Remove NAs
            outcome_data_T = outcome_data_T.dropna(how='any')

            # Convert cycle columns into an int
            entropy_df["cycle"] = entropy_df["cycle"].str.extract(r"C(\d+)").astype(int)
            outcome_data_T["TFx_cycle"] = outcome_data_T["TFx_cycle"].str.extract(r"TFx_C(\d+)").astype(int)

            # merge data
            merged_df = pd.merge(entropy_df, outcome_data_T, left_on=['patient_id', 'cycle'], right_on=['Sample_ID', 'TFx_cycle'], how='inner').drop(columns={'Sample_ID', 'TFx_cycle'})
            merged_df.columns = merged_df.columns.str.strip()

            # Convert TFx from object to numeric
            merged_df["TFx"] = pd.to_numeric(merged_df["TFx"], errors="coerce")

            # Drop NAs
            merged_df = merged_df.dropna(how='any')

            # Add response category and code them as 0/1 (responder = 0)
            if response_category == 'extreme_responders':
                response_categories_df = pd.read_csv(extreme_responders_categories_csv)
                response_categories_df = response_categories_df[['Sample_ID', 'Group']]
                response_categories_df['Response'] = (response_categories_df['Group'] == 'Non-Responder').astype(int)
            elif response_category == 'overall_survival':
                response_categories_df = outcome_data_filtered[['Sample_ID', 'survival_days']]
                response_categories_df['Response'] = (response_categories_df['survival_days'] < 252).astype(int)
                response_categories_df.rename(columns={'survival_days': 'Group'}, inplace=True)
            else:
                print('Error: no response category given (e.g. \'extreme_responders\', \'overall_survival\')')
                return

            merged_df = pd.merge(merged_df, response_categories_df, left_on='patient_id', right_on='Sample_ID', how='inner').drop(columns={'Sample_ID', 'Group'})

            # Run linear mixd effects model
            model_result = linear_mixed_effects_model(
                df=merged_df,
                model_equation=f'chr8 ~ cycle*Response + TFx',
                id_col='patient_id',
            )

            linear_mixed_effects_forest_plot(
                model_result=model_result,
                output_path=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process12_outputs/{csv_path}/linear_mixed_effects_model/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct',
                file_name=f'linear_mixed_effects_forest_plot',
            )

            linear_mixed_effects_spaghetti_plot(
                df=merged_df,
                id_col='patient_id',
                x_col='cycle',
                y_col='chr8',
                response_col='Response',
                x_label='cycles',
                y_label='chr8 entropy',
                model_result=model_result,
                model_pval_col='cycle:Response',
                output_path=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process12_outputs/{csv_path}/linear_mixed_effects_model/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct',
                file_name=f'linear_mixed_effects_spaghetti_plot',
            )

            linear_mixed_effects_split_spaghetti_plot(
                df=merged_df,
                id_col='patient_id',
                x_col='cycle',
                y_col='chr8',
                response_col='Response',
                x_label='cycles',
                y_label='chr8 entropy',
                model_result=model_result,
                model_pval_col='cycle:Response',
                output_path=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process12_outputs/{csv_path}/linear_mixed_effects_model/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct',
                file_name=f'linear_mixed_effects_split_spaghetti_plot',
            )

            longitudinal_violin_split_plot(
                df=merged_df,
                id_col='patient_id',
                x_col='cycle',
                y_col='chr8',
                response_col='Response',
                x_label='cycles',
                y_label='chr8 entropy',
                model_result=model_result,
                model_pval_col='cycle:Response',
                output_path=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process12_outputs/{csv_path}/linear_mixed_effects_model/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct',
                file_name=f'linear_mixed_effects_split_violin_plot',
            )

            # Plot on KM curve
            metric_to_use = "Corrected_Copy_Number"

            # Read in csvs and combine
            ctdna_sheet_df = pd.read_csv(input_map["all_samples_tfx_and_ctdnaq_sheet"][0], index_col=0)
            entropy_df = pd.read_csv(f'{output_path}/data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_path}.csv')
            combined_df = pd.merge(ctdna_sheet_df, entropy_df, left_on=['Sample_ID', 'Cycle'], right_on=['patient_id', 'cycle'], how='inner')

            # Drop unnessasary columns
            combined_df = combined_df.drop(columns={'Sample_ID', 'Cycle'})

            # Rename columns to not include dashes
            combined_df = combined_df.rename(columns={'cfdna-q': 'cfdna_q'})

            # Remove patients less than TFx cutoff
            combined_df = combined_df[
                (combined_df['cycle'] != 'C1') |
                (combined_df['TFx'] >= tfx_cutoff)
            ]

            cycles = ["C1", "C2"]

            df_pivot = combined_df.pivot(
                index='patient_id',
                columns='cycle',
                values='chr8'
            ).reset_index()

            # Remove patients with NAs for C1
            df_pivot = df_pivot[df_pivot[cycles[0]].notna()]

            df_pivot = add_max_change_columns(df_pivot, cycles)

            # Drop NANs in df_pivot
            df_pivot = df_pivot.dropna(subset=['pct_change'])

            prct_change_df = df_pivot[['patient_id', 'pct_change']]

            # Grabs and format HR and CI from cox model 
            cox_model_hr_text = f'N/a'

            km_curve_output_path = f'entropy-analysis-{sequencing_type}-{metric_to_use}/process12_outputs/kaplan-meier-curves/entropy_km_curves/{csv_path}/tfx_threshold_{tfx_cutoff*100}_prct'

            master_sheet = input_map["pluvicto_master_sheet"][0]
            ctdna_sheet = input_map["all_samples_tfx_and_ctdnaq_sheet"][0]
            master_sheet_df = pd.read_csv(master_sheet, index_col=0)
            ctdna_sheet_df = pd.read_csv(ctdna_sheet, index_col=0)

            # Subset ctdna_sheet_df to C1 and remove '-' from column names
            ctdna_sheet_df = ctdna_sheet_df[ctdna_sheet_df['Cycle'] == 'C1']
            ctdna_sheet_df = ctdna_sheet_df.rename(columns={'cfdna-q': 'cfdna_q'})

            df_combined = pd.merge(master_sheet_df, ctdna_sheet_df, left_on='Sample_ID', right_on='Sample_ID')

            # KM curve based on (high TFx, high entropy) vs (high TFx, low entropy) vs (low TFx, high entropy) vs (low TFx, low entropy)
            kaplan_meier_plot(
                output_path, 
                km_curve_output_path,
                df_combined,
                "Death", 
                "survival_days", 
                "pct_change",
                prct_change_df, 
                "TFx_C1",
                "median",
                metric_to_use,
                sequencing_type,
                None,
                None,
            )

    # ==============
    # Process 13:
    # ==============
    if args.process13 == "True":
        # List of csv files to process
        entropy_csv_list = [
            "base_entropy_hn_cohort_normalized_per_chr_table",
        ]

        baseline_sample_entropy_characteristic = True

        # =======================================================================
        # Creating a supplemental table with entropy, TFx, ctDNA-Q, and SV status
        # =======================================================================
        if baseline_sample_entropy_characteristic:
            # Load data

            complex_sv_csv = 'path/to/complex_svs_type_by_tool_table.csv'
            patrick_annotations = 'path/to/entropy_curations.csv'

            complex_svs_df = pd.read_csv(complex_sv_csv)
            ctdna_sheet_df = pd.read_csv(input_map["all_samples_tfx_and_ctdnaq_sheet"][0], index_col=0)
            patrick_annotations_df = pd.read_csv(patrick_annotations)

            # ULP Entropy base_entropy_hn_cohort_normalized
            ulp_curated_weighted_by_bin_cohort_normalized_column = get_entropy_column(output_path, f"path/to/scripts/outputs/data-tables/entropy-tables-ulp_curated_seg-Corrected_Copy_Number/base_entropy_hn_cohort_normalized_per_chr_table.csv", "chr8", "C1")

            # Deep ichor base_entropy_hn_cohort_normalized
            deep_ichor_v1_weighted_by_bin_cohort_normalized_column = get_entropy_column(output_path, f"path/to/scripts/outputs/data-tables/entropy-tables-deep_ichor_seg_v1-Corrected_Copy_Number/base_entropy_hn_cohort_normalized_per_chr_table.csv", "chr8", "C1")

            entropy_dfs = [
                (ulp_curated_weighted_by_bin_cohort_normalized_column, 'chr8_ulp_curated_weighted_by_bin_cohort_normalized'),
                (deep_ichor_v1_weighted_by_bin_cohort_normalized_column, 'chr8_deep_ichor_v1_weighted_by_bin_cohort_normalized'),
            ]
            
            # Filter before merging
            # complex_svs_df = complex_svs_df.drop(columns=['chr8'])
            ctdna_sheet_df = ctdna_sheet_df[ctdna_sheet_df['Cycle'] == 'C1']
            ctdna_sheet_df = ctdna_sheet_df[['Sample_ID', 'TFx', 'cfdna-q']].rename(columns={'cfdna-q': 'cfdna_q'}).reset_index(drop=True)

            # Remove tag on id if there are no repeats
            if len(list(complex_svs_df['patient_id'])) == len(set(list(complex_svs_df['patient_id']))):
                complex_svs_df['patient_id'] = complex_svs_df['patient_id'].str.split('_').str[0]
            else:
                print('error')
                sys.exit(0)

            # Merge
            merge_df = pd.merge(complex_svs_df, ctdna_sheet_df, left_on='patient_id', right_on='Sample_ID', how='outer').drop(columns=['patient_id']).rename(columns={'Sample_ID': 'patient_id'})
            merge_df = pd.merge(merge_df, patrick_annotations_df, left_on='patient_id', right_on='Sample_ID', how='outer').drop(columns=['Sample_ID'])

            # Merge all entropy dataframes
            for entropy_df, col_name in entropy_dfs:
                merge_df = pd.merge(merge_df, entropy_df[['patient_id', 'chr8']], on='patient_id', how='outer')
                merge_df = merge_df.rename(columns={'chr8': col_name})

            # Add binary column to show if patient has deep data or not
            merge_df['has_deep_data'] = merge_df['chr8_deep_ichor_v1_weighted_by_bin_cohort_normalized'].apply(lambda x: 'Y' if pd.notna(x) else 'N')
            
            # Format columns
            final_df = merge_df[['patient_id', 'has_deep_data', 'Entropy?', *[col for _, col in entropy_dfs], 'TFx', 'cfdna_q', 'AA', 'jabba', 'Notes']]

            final_df.to_csv(f'{output_path}/data-tables/entropy_ULP_deepTitan_deepIchor_and_sv_types.csv', index=False)

            print(f'csv saved to: {output_path}/data-tables/entropy_deep_ichor_and_sv_types.csv')

    # ==============
    # Process 14:
    # ==============
    if args.process14 == "True":
        # List of csv files to process
        entropy_csv_list = [
            "base_entropy_hn_cohort_normalized_per_chr_table",
        ]
    
        tfx_entropy_scatter_plot = True

        tfx_cutoff = 0.0

        for csv_path in entropy_csv_list:

            if tfx_entropy_scatter_plot:
                
                # Load data
                ctdna_sheet_df = pd.read_csv(input_map["all_samples_tfx_and_ctdnaq_sheet"][0], index_col=0)
                
                # ULP Entropy base_entropy_hn_cohort_normalized
                ulp_curated_weighted_by_bin_cohort_normalized_column = get_entropy_column(output_path, f"path/to/scripts/outputs/data-tables/entropy-tables-ulp_curated_seg-Corrected_Copy_Number/base_entropy_hn_cohort_normalized_per_chr_table.csv", "chr8", "C1")

                # Deep ichor base_entropy_hn_cohort_normalized
                deep_ichor_v1_weighted_by_bin_cohort_normalized_column = get_entropy_column(output_path, f"/path/to/scripts/outputs/data-tables/entropy-tables-deep_ichor_seg_v1-Corrected_Copy_Number/base_entropy_hn_cohort_normalized_per_chr_table.csv", "chr8", "C1")
                
                # List of entropys used
                entropy_dfs = [
                    (ulp_curated_weighted_by_bin_cohort_normalized_column, 'chr8_ulp_curated_weighted_by_bin_cohort_normalized'),
                    (deep_ichor_v1_weighted_by_bin_cohort_normalized_column, 'chr8_deep_ichor_v1_weighted_by_bin_cohort_normalized'),
                ]
                
                # Filter before merging
                # complex_svs_df = complex_svs_df.drop(columns=['chr8'])
                ctdna_sheet_df = ctdna_sheet_df[ctdna_sheet_df['Cycle'] == 'C1']
                ctdna_sheet_df = ctdna_sheet_df[['Sample_ID', 'TFx', 'cfdna-q']].rename(columns={'cfdna-q': 'cfdna_q'}).reset_index(drop=True).rename(columns={'Sample_ID': 'patient_id'})

                # Merge
                merge_df = ctdna_sheet_df

                # Merge all entropy dataframes
                for entropy_df, col_name in entropy_dfs:
                    merge_df = pd.merge(merge_df, entropy_df[['patient_id', 'chr8']], on='patient_id', how='outer')
                    merge_df = merge_df.rename(columns={'chr8': col_name})

                # TFx cutoff
                merge_df = merge_df[merge_df['TFx'] >= tfx_cutoff]
                
                # Scatter plot of ULP entropy vs TFx
                scatter_plot_two_columns_one_dataframe(
                    df=merge_df, 
                    col_1_name='TFx', 
                    col_2_name='chr8_ulp_curated_weighted_by_bin_cohort_normalized', 
                    title='TFx vs ULP Entropy', 
                    xlabel='TFx', 
                    ylabel='Entropy', 
                    out_dir=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process14_outputs/{csv_path}/TFx_vs_entropy_scatter_plots/tfx_threshold_{tfx_cutoff*100}_prct', 
                    file_name="TFx_vs_ULP_entropy_scatter_plot", 
                    jitter=False
                )

                # Scatter plot of deep ichor entropy vs TFx
                scatter_plot_two_columns_one_dataframe(
                    df=merge_df, 
                    col_1_name='TFx', 
                    col_2_name='chr8_deep_ichor_v1_weighted_by_bin_cohort_normalized', 
                    title='TFx vs Deep Ichor Entropy', 
                    xlabel='TFx', 
                    ylabel='Entropy', 
                    out_dir=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process14_outputs/{csv_path}/TFx_vs_entropy_scatter_plots/tfx_threshold_{tfx_cutoff*100}_prct', 
                    file_name="TFx_vs_deep_ichor_entropy_scatter_plot", 
                    jitter=False
                )
        
    # ==============
    # Process 15:
    # ==============
    if args.process15 == "True":

        if 'seg' in sequencing_type:
            print('This process uses bins not segments. Please swithc sequencing_type.')
            sys.exit(0)

        # List of csv files to process
        gene_list = {
            "MYC":('chr8', 127735434, 127742951),
            "RAD21":('chr8', 116845934, 116874776)
        }

        entropy_csv_list = [
            "base_entropy_hn_cohort_normalized_per_chr_table",
        ]

        extreme_responders_categories_csv = 'path/to/proteus_deg_sample_classification.csv'

        metric_to_use = "Corrected_Copy_Number"

        # Read cohort wide patient copy number dataframe
        patient_cn_df = pd.read_pickle(
            f'{output_path}/data-tables/cohort-copy-number-profiles-dataframes/{sequencing_type}-{metric_to_use}-cohort-copy_number_profiles.pkl'
        )

        cycles_to_use = ['TFx_C1', 'TFx_C2']#, 'TFx_C3', 'TFx_C4', 'TFx_C5', 'TFx_C6']
        tfx_cutoff = 0.10
        response_category = 'extreme_responders'

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:
            for gene, (chrom, start, end) in gene_list.items():

                # Filter df
                gene_df = patient_cn_df[
                    (patient_cn_df["chr"] == chrom) &
                    (patient_cn_df["start"] <= start) &
                    (patient_cn_df["end"] >= end)
                ].copy()

                # Read outcome_data
                outcome_data = pd.read_csv(input_map["pluvicto_master_sheet"][0])

                # Remove patients less than TFx cutoff
                outcome_data_filtered = outcome_data[outcome_data[cycles_to_use[0]] >= tfx_cutoff]

                # Filter columns
                outcome_data_slice = outcome_data_filtered[['Sample_ID', *cycles_to_use]]

                # Convert outcome_data from wide to long
                outcome_data_T = pd.melt(
                    frame=outcome_data_slice,
                    id_vars=['Sample_ID'],
                    value_vars=cycles_to_use,
                    var_name='TFx_cycle',
                    value_name='TFx'
                )
                # Remove NAs
                outcome_data_T = outcome_data_T.dropna(how='any')

                # Convert cycle columns into an int
                gene_df["cycle"] = gene_df["cycle"].str.extract(r"C(\d+)").astype(int)
                outcome_data_T["TFx_cycle"] = outcome_data_T["TFx_cycle"].str.extract(r"TFx_C(\d+)").astype(int)

                # merge data
                merged_df = pd.merge(gene_df, outcome_data_T, left_on=['patient', 'cycle'], right_on=['Sample_ID', 'TFx_cycle'], how='inner').drop(columns={'Sample_ID', 'TFx_cycle'})
                merged_df.columns = merged_df.columns.str.strip()

                # Convert TFx from object to numeric
                merged_df["TFx"] = pd.to_numeric(merged_df["TFx"], errors="coerce")

                # Drop NAs
                merged_df = merged_df.dropna(how='any')

                # Add response category and code them as 0/1 (responder = 0)
                if response_category == 'extreme_responders':
                    response_categories_df = pd.read_csv(extreme_responders_categories_csv)
                    response_categories_df = response_categories_df[['Sample_ID', 'Group']].copy()
                    response_categories_df['Response'] = (response_categories_df['Group'] == 'Non-Responder').astype(int)
                elif response_category == 'overall_survival':
                    response_categories_df = outcome_data_filtered[['Sample_ID', 'survival_days']].copy()
                    response_categories_df['Response'] = (response_categories_df['survival_days'] < 252).astype(int)
                    response_categories_df.rename(columns={'survival_days': 'Group'}, inplace=True)
                else:
                    print('Error: no response category given (e.g. \'extreme_responders\', \'overall_survival\')')
                    return

                merged_df = pd.merge(merged_df, response_categories_df, left_on='patient', right_on='Sample_ID', how='inner').drop(columns={'Sample_ID', 'Group'})

                # Run linear mixd effects model
                model_result = linear_mixed_effects_model(
                    df=merged_df,
                    model_equation=f'Corrected_Copy_Number ~ cycle*Response + TFx',
                    id_col='patient',
                )

                # Save model_results
                model_result_df = pd.DataFrame({
                    'coef': model_result.params,
                    'std_err': model_result.bse,
                    'z': model_result.tvalues,
                    'p_value': model_result.pvalues
                })
                model_result_df[['CI_lower', 'CI_upper']] = model_result.conf_int()
                
                model_result_df.to_csv(f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process15_outputs/{csv_path}/linear_mixed_effects_model/{gene}/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct/linear_mixed_effects_model_results.csv')

                linear_mixed_effects_forest_plot(
                    model_result=model_result,
                    output_path=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process15_outputs/{csv_path}/linear_mixed_effects_model/{gene}/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct',
                    file_name=f'linear_mixed_effects_forest_plot',
                )

                linear_mixed_effects_spaghetti_plot(
                    df=merged_df,
                    id_col='patient',
                    x_col='cycle',
                    y_col='Corrected_Copy_Number',
                    response_col='Response',
                    x_label='cycles',
                    y_label=f'{gene} Copy Number',
                    model_result=model_result,
                    model_pval_col='cycle:Response',
                    output_path=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process15_outputs/{csv_path}/linear_mixed_effects_model/{gene}/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct',
                    file_name=f'linear_mixed_effects_spaghetti_plot',
                )

                linear_mixed_effects_split_spaghetti_plot(
                    df=merged_df,
                    id_col='patient',
                    x_col='cycle',
                    y_col='Corrected_Copy_Number',
                    response_col='Response',
                    x_label='cycles',
                    y_label=f'{gene} Copy Number',
                    model_result=model_result,
                    model_pval_col='cycle:Response',
                    output_path=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process15_outputs/{csv_path}/linear_mixed_effects_model/{gene}/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct',
                    file_name=f'linear_mixed_effects_split_spaghetti_plot',
                )

                longitudinal_violin_split_plot_v2(
                    df=merged_df,
                    id_col='patient',
                    x_col='cycle',
                    y_col='Corrected_Copy_Number',
                    response_col='Response',
                    x_label='cycles',
                    y_label=f'{gene} Copy Number',
                    model_result=model_result,
                    model_pval_col='cycle:Response',
                    output_path=f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process15_outputs/{csv_path}/linear_mixed_effects_model/{gene}/response_category_{response_category}/cycles_{cycles_to_use[0][-2:]}_to_{cycles_to_use[len(cycles_to_use)-1][-2:]}/tfx_threshold_{tfx_cutoff*100}_prct',
                    file_name=f'linear_mixed_effects_split_violin_plot',
                )

    # ==============
    # Process 16:
    # ==============
    if args.process16 == "True":
        entropy_csv_list = [
            "base_entropy_hn_cohort_normalized_per_chr_table",
        ]

        metric_to_use = "Corrected_Copy_Number"

        tfx_cutoff = 0.10

        z_score = True

        # Set CSV file variable to use in process
        for csv_path in entropy_csv_list:
            # Initialize combined model dataframe
            forest_df = pd.DataFrame({})

            for x in list(range(1, 23)) + ['X']: 
                # Get entropy column
                entropy_column = get_entropy_column(
                    output_path=output_path, 
                    entropy_file_path=f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv", 
                    chr_column=f'chr{x}', 
                    cycle="C1"
                )

                # Read in pluvicto survival data
                master_sheet_df = pd.read_csv(input_map["pluvicto_master_sheet"][0], index_col=0)

                # Filter out samples below 10% TFx
                master_sheet_df = master_sheet_df[master_sheet_df['TFx_C1'] > tfx_cutoff]
                
                # Collect TFx column
                tfx_column = master_sheet_df[['Sample_ID','TFx_C1', 'Death', 'survival_days']]

                # Merge and filter pluvicto master sheet and entropy
                merged_df = pd.merge(entropy_column, master_sheet_df, left_on='patient_id', right_on='Sample_ID', how='inner')
                filtered_df = merged_df[['patient_id', 'TFx_C1', 'Death', 'survival_days', f'chr{x}']].copy()

                # Z-score data so it plots by a one unit increase in standard deviation
                if z_score:
                    filtered_df[f'zscore_chr{x}'] = (filtered_df[f'chr{x}'] - filtered_df[f'chr{x}'].mean()) / filtered_df[f'chr{x}'].std()
                    filtered_df.drop(columns=[f'chr{x}'], inplace=True)
                    filtered_df.rename(columns={f'zscore_chr{x}': f'chr{x}'}, inplace=True)

                cox_model = cox_proportional_hazards_model(
                    cox_data=filtered_df, 
                    time_to_event_column='survival_days', 
                    event_col='Death', 
                    formula=f'chr{x} + TFx_C1', 
                    chr_to_pick=f'chr{x}'
                )

                # Prepare df for combine with forest_df
                df = cox_model.summary.copy()
                if df.index[0] != f'chr{x}':
                    print('index does not match.')
                    sys.exit(0)
                df.index = [f'chr{x}_entropy', f'chr{x}_TFx']

                # Combine
                forest_df = pd.concat([forest_df, df])

            # Filter froest_df
            forest_df = forest_df.loc[forest_df.index.isin([f'chr{x}_entropy' for x in list(range(1, 23)) + ['X']])]

            # If z_score is true save in different directory 
            if z_score:
                file_path = f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process16_outputs/{csv_path}/cox_hardard_ratios/entropy_plots/tfx_threshold_{tfx_cutoff*100}_prct/zscored_data'
            else:
                file_path = f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process16_outputs/{csv_path}/cox_hardard_ratios/entropy_plots/tfx_threshold_{tfx_cutoff*100}_prct'

            file_name=f'all_chromosomes_entropy_cox_model_forest_plot_no_FDR'

            out_path, summary = cox_forest_plot_pub(
                forest_df, 
                file_path, 
                'p',
                filename=file_name,
                title='Cox Proportional Hazards Model (Ctrl for TFx - No FDR)'
            )

            print(f'\nFile saved at: {file_path}/{file_name}\n')

            # ---------------------------
            # Model with all chrs entropy
            # ---------------------------
            entropy_df = pd.read_csv(f'{output_path}/data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{csv_path}.csv')
            entropy_df = entropy_df[entropy_df['cycle'] == 'C1'].copy()
            entropy_df.drop(columns=['cycle'], inplace=True)
            
            # Read in pluvicto survival data
            master_sheet_df = pd.read_csv(input_map["pluvicto_master_sheet"][0], index_col=0)

            # Filter out samples below 10% TFx
            master_sheet_df = master_sheet_df[master_sheet_df['TFx_C1'] > tfx_cutoff]
            
            # Collect TFx column
            tfx_column = master_sheet_df[['Sample_ID','TFx_C1', 'Death', 'survival_days']]

            # Merge and filter pluvicto master sheet and entropy
            merged_df = pd.merge(entropy_df, master_sheet_df, left_on='patient_id', right_on='Sample_ID', how='inner')
            filtered_df = merged_df[['patient_id', 'TFx_C1', 'Death', 'survival_days', *[f'chr{x}' for x in list(range(1, 23)) + ['X']]]].copy()

            # Z-score data so it plots by a one unit increase in standard deviation
            if z_score:
                for x in list(range(1, 23)) + ['X']:
                    filtered_df[f'zscore_chr{x}'] = (filtered_df[f'chr{x}'] - filtered_df[f'chr{x}'].mean()) / filtered_df[f'chr{x}'].std()
                    filtered_df.drop(columns=[f'chr{x}'], inplace=True)
                    filtered_df.rename(columns={f'zscore_chr{x}': f'chr{x}'}, inplace=True)

            chr_formual_part = " + ".join([f"chr{x}" for x in list(range(1, 23)) + ["X"]])

            cox_model = cox_proportional_hazards_model(
                cox_data=filtered_df, 
                time_to_event_column='survival_days', 
                event_col='Death', 
                formula=f'{chr_formual_part} + TFx_C1', 
                chr_to_pick=f'chr{x}'
            )

            # If z_score is true save in different directory 
            if z_score:
                file_path = f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process16_outputs/{csv_path}/cox_hardard_ratios/entropy_plots/tfx_threshold_{tfx_cutoff*100}_prct/zscored_data'
            else:
                file_path = f'{output_path}/entropy-analysis-{sequencing_type}-Corrected_Copy_Number/process16_outputs/{csv_path}/cox_hardard_ratios/entropy_plots/tfx_threshold_{tfx_cutoff*100}_prct'

            file_name=f'grouped_entropy_cox_model_forest_plot_no_FDR'

            out_path, summary = cox_forest_plot_pub(
                cox_model.summary.copy(), 
                file_path, 
                'p',
                filename=file_name,
                title='Cox Proportional Hazards Model (Ctrl for TFx - No FDR)'
            )

            print(f'\nFile saved at: {file_path}/{file_name}\n')

# If this script is run directly, it will execute the main function
# This allows the script to be run from the command line with arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plots tumor fraction vs fraction genome altered.")

    parser.add_argument("--process1", default="False", help="Creates and saves look up table for copy number profiles.")
    parser.add_argument("--process2", default="False", help="Entropy per chromosome data processing.")
    parser.add_argument("--process3", default="False", help="Entropy per chromosome data plotting.")
    parser.add_argument("--process4", default="False", help="Entropy per chromosome data qualitative analysis.")
    parser.add_argument("--process5", default="False", help="Entropy per chromosome data cox model analysis.")
    parser.add_argument("--process6", default="False", help="Chromosome 8 complex sv events and entropy visualization.")
    parser.add_argument("--process7", default="False", help="Chromosome 8 Gene Matrix and entropy correlation (logistic regression).")
    parser.add_argument("--process8", default="False", help="Chromosome 8 complex sv analysis (logistic regression).")
    parser.add_argument("--process9", default="False", help="Comparing chr8 survival analysis between different datasets.")
    parser.add_argument("--process10", default="False", help="Chromosome 8 complex sv analysis (linear regression).")
    parser.add_argument("--process11", default="False", help="Chromosome 8 Gene Matrix and entropy correlation (linear regression).")
    parser.add_argument("--process12", default="False", help="Linear mixed effects model for entropy.")
    parser.add_argument("--process13", default="False", help="Creating Manuscript supplemental tables.")
    parser.add_argument("--process14", default="False", help="Entropy Visualizations - extention of process3.")
    parser.add_argument("--process15", default="False", help="Linear mixed effects model for specific genes CN increase/decrease over time.")
    parser.add_argument("--process16", default="False", help="Forest plot with all chromosomes on it.")

    args = parser.parse_args()
    main(args)