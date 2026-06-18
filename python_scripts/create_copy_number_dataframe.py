# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 06/04/26
# Purpose: Creates a csv of all patients copy number profile and related meta data.
# Note: Functions adopted from entropy_by_chromosome.py

from python_scripts.imports import *

def create_copy_number_profile_dataframe(dir_list, metric_to_use, sequencing_type):
    '''
    Parameters:
    -----------
        dir_list: List
            Array holding directory paths.

        metric_to_use: String
            Parameter to determine which metric_to_use from cna.seg file to use for entropy.

        sequencing_type: String
            Determins what sequencying type to select for (e.g. 'deep', 'ulp', or 'ulp_curated').

    Function:
    ---------
        - Loops through each directory in the genomic instability directory array.
        - For each patient found, extracts the copy number profile for each chromosome. 

    Returns:
    --------
        None (modifies patient_dict in place)

    '''
    print("Processing copy number profiles")
    
    patient_df = []

    if 'ulp_curated_seg' in sequencing_type:
        curated_files_df = pd.read_csv(dir_list[0])
        for index in curated_files_df.index:
            full_path = curated_files_df.iloc[index]["curated_solution_path"].replace('.cna','') + '.txt'
            
            # read file and append it to patient_df
            df = pd.read_csv(full_path, sep='\t')

            # Rename ID column to sample_id
            df.rename(columns={'ID': 'sample_id'}, inplace=True)

            # Find parts of the the sample_id
            sample_id = df.loc[0]['sample_id']
            parts = sample_id.split("_")
            patient = parts[0]
            timepoint = parts[-1]

            # Add new columns into the second and third columns
            df.insert(1, "patient", patient)
            df.insert(2, "cycle", timepoint)

            # Rename columns
            df.rename(columns={'chrom': 'chr'}, inplace=True)

            patient_df.append(df)

    elif 'ulp_curated' in sequencing_type:
        curated_files_df = pd.read_csv(dir_list[0])
        for index in curated_files_df.index:
            full_path = curated_files_df.iloc[index]["curated_solution_path"]
            
            # read file and append it to patient_df
            df = pd.read_csv(full_path, sep='\t')

            # Find parts of the the sample_id
            sample_id = curated_files_df.iloc[index]["Identifier"]
            parts = sample_id.split("_")
            patient = parts[0]
            timepoint = parts[-1]

            # Add new columns
            df.insert(0, "sample_id", sample_id)
            df.insert(1, "patient", patient)
            df.insert(2, "cycle", timepoint)

            # Rename columns to remove sample id
            df.rename(
                columns=lambda col: col.replace(f'{sample_id}.', ''),
                inplace=True
            )

            patient_df.append(df)

    elif 'deep_ichor_seg' in sequencing_type:
        curated_files_df = pd.read_csv(dir_list[0])

        for index in curated_files_df.index:
            full_path = curated_files_df.iloc[index]["curated_solution_path"]
            
            # read file and append it to patient_df
            df = pd.read_csv(full_path, sep='\t')

            # Find parts of the the sample_id
            sample_id = curated_files_df.iloc[index]["Identifier"]
            print(f'{sample_id}')
            parts = sample_id.split("_")
            patient = parts[0]
            timepoint = parts[-1]

            # Add new columns
            df.insert(0, "sample_id", sample_id)
            df.insert(1, "patient", patient)
            df.insert(2, "cycle", timepoint)

            # Rename columns to remove sample id
            df.rename(
                columns=lambda col: col.replace(f'{sample_id}.', ''),
                inplace=True
            )

            patient_df.append(df)

    master_df = pd.concat(patient_df, ignore_index=True)

    return master_df