# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 12/4/2025
# Purpose: Calculates a penalized Shannon entropy score per chromosome and outputs a tsv listing patients in decending order by FGA per chromosome.
# Note: Functions adopted from genome_instability_v1.py

from python_scripts.imports import *

warnings.filterwarnings('ignore', category=pd.errors.DtypeWarning)

# Convert both warnings to exceptions
warnings.simplefilter('error', PerfectSeparationWarning)
warnings.simplefilter('error', ConvergenceWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings('error', category=RuntimeWarning, message='overflow encountered in exp')
warnings.filterwarnings('error', category=RuntimeWarning, message='divide by zero encountered in log')

def find_entropy_per_chromosome(full_path, patient, metric_to_use, sequencing_type, bin_range = None, entropy_equation_function=None, max_K=None, file_ext=".cna.seg"):
    '''
    Parameters:
    -----------
        full_path: String
            A path to the patient directory.

        patient: String
            Id of the pateint.

        metric_to_use: String
            Parameter to determine which metric_to_use from cna.seg file to use for entropy.

        sequencing_type: String
            Determins what sequencying type to select for (e.g. 'deep', 'ulp', or 'ulp_curated').

        bin_range: int (defult = None)
            Used to set a range on how many bins to select in deep data to manage entropy score (not used for ulp).

        file_ext: String (default = ".cna.seg")
            File extention if parameter is not None. If file_ext is None, assume full_path parameter has an extention. 
    
    Function:
    ---------
        - Collects the cna.seg file from the patient directory.
        - Collects the corrected copy number for each bin across the patient.
        - Calculates entropy per chromosome (with/without error term and/or normalization), adds it to a dictionary, and returns it.

    Returns:
    --------
        copy_number_corrected_dict: Dictionary
            A dictionary with chromosome keys and entropy values.

    '''
    if file_ext is not None:
        # Finds seg file and opens it
        full_path = os.path.join(full_path, patient + file_ext)

    if 'deep_ichor' in sequencing_type:
        full_path = list(Path(full_path).rglob(f"{patient}/*_optimal/*.cna.seg"))[0]

    # Initialize an empty list to store the values from copy number corrected column
    copy_number_corrected_dict = {}

    # Sets the index for the metric and the chromosome id column (Titan has patient id as index 1)
    if metric_to_use == "Corrected_Copy_Number":
        if sequencing_type == "deep":
            n = 23
            chr_coloumn = 1
        else:
            n = 7
            chr_coloumn = 0
    elif metric_to_use == "logR":
        if sequencing_type == "deep":
            n = 12
            chr_coloumn = 1
        else:
            n = 5
            chr_coloumn = 0

    # Process sequencing types differently because Titan has dynamic column lengths, ichor is static
    if 'deep' in sequencing_type and 'ichor' not in sequencing_type:
        df = pd.read_csv(full_path, sep=r"\s+")

        # for idx, row in df.iterrows():
        #     patient_copy_number_corrected = row[metric_to_use]
        #     # Skips over header
        #     if metric_to_use == patient_copy_number_corrected:
        #         continue
        #     # Append the value from copy number corrected column to the list
        #     copy_number_corrected_dict.setdefault(row['Chromosome'], []).append(patient_copy_number_corrected)

        # start, end = 0, 0
        # temp_copy_number_list = []
        # last_seen_chr = ""
        # for _, row in df_filtered.iterrows():
        #     curr_place = row["Start"]
        #     val = row[metric_to_use]
        #     temp_row = row
        #     chr_val = row["Chr"]

        #     if chr_val != last_seen_chr:
        #         last_seen_chr = chr_val
        #         start = curr_place
        #         if start < 500000:
        #             end = bin_range
        #         elif start < bin_range:
        #             end = bin_range + bin_range
        #         else:
        #             end = start - (start % bin_range) + bin_range

        #         print(f"Chr: {chr_val}")    
        #         print(f"start: {start}")    
        #         print(f"end: {end}")    

        #     # If bin_range is not None it will take the median across the bin size (so it doesnt have so much noise) - else it will take every non NA bin into the calculation
        #     if bin_range is not None:
        #         if curr_place < end:
        #             if pd.isna(val):
        #                 continue
        #             temp_copy_number_list.append(float(val))
        #         else:
        #             if not temp_copy_number_list:
        #                 end = curr_place + bin_range
        #                 if pd.isna(val):
        #                     temp_copy_number_list = []
        #                 else:
        #                     temp_copy_number_list = [float(val)]
        #                 continue

        #             # Find median value from temp list and reset counter and list
        #             median_val = statistics.median(temp_copy_number_list)
        #             end = curr_place + bin_range
        #             copy_number_corrected_dict.setdefault(chr_val, []).append(median_val)

        #             print(f"median_val: {median_val}")
        #             print(f"end: {end}")

        #             if pd.isna(val):
        #                 temp_copy_number_list = []
        #             else:
        #                 temp_copy_number_list = [float(val)]   
        #     else:
        #         copy_number_corrected_dict.setdefault(chr_val, []).append(float(val))

        # if temp_copy_number_list:
        #     print("here")
        #     median_val = statistics.median(temp_copy_number_list)
        #     copy_number_corrected_dict.setdefault(temp_row["Chr"], []).append(median_val)

    else:
        file = open(full_path, "r")
        # Read the file line by line
        for line in file:
            patient_data = line.strip().split("\t")
            patient_copy_number_corrected = patient_data[n]
            # Skips over header
            if metric_to_use in patient_copy_number_corrected:
                continue
            # Append the value from copy number corrected column to the list
            copy_number_corrected_dict.setdefault(patient_data[chr_coloumn], []).append(patient_copy_number_corrected)

        # Close the file
        file.close()

    # Sorts the dictionary by keys so it is in natural order
    def chr_key(chrom):
        """Extract number from chromosome name for natural sorting"""
        chrom_str = str(chrom).replace('chr', '').replace('Chr', '')
        try:
            return int(chrom_str)
        except ValueError:
            # For X, Y, M - put them at the end
            return float('inf')

    copy_number_corrected_dict = dict(sorted(copy_number_corrected_dict.items(), key=lambda x: chr_key(x[0])))

    # Check if every copy number value is 2 or not.
    # If value != 2 count it as a copy number alteration
    for chrom, values in copy_number_corrected_dict.items():

        # Calculates entropy per chromosome (with an error term for expected ploidy 2)
        len_of_chrom = sum(values)
        counts = Counter(values)
        counts_probs = {int(key): count / len_of_chrom for key, count in counts.items()}
        
        copy_number_corrected_dict[chrom] = entropy_equation_function(counts_probs)

        # Update max_K 
        max_K[0] = len(counts_probs) if len(counts_probs) > max_K[0] else max_K[0]

        print(copy_number_corrected_dict)
        
    # Return genomic instability per chrom
    return copy_number_corrected_dict

def find_entropy_per_chromosome_seg_file(full_path, patient, metric_to_use, sequencing_type, bin_range=None, entropy_equation_function=None, max_K=None, file_ext=".cna.seg"):
    '''
    Parameters:
    -----------
        full_path: String
            A path to the patient directory.

        patient: String
            Id of the pateint.

        metric_to_use: String
            Parameter to determine which metric_to_use from cna.seg file to use for entropy.

        sequencing_type: String
            Determins what sequencying type to select for (e.g. 'deep', 'ulp', or 'ulp_curated').

        bin_range: int (defult = None)
            Used to set a range on how many bins to select in deep data to manage entropy score (not used for ulp).

        file_ext: String (default = ".cna.seg")
            File extention if parameter is not None. If file_ext is None, assume full_path parameter has an extention. 
    
    Function:
    ---------
        - Collects the .seg.txt file from the patient directory.
        - Collects the corrected copy number for each bin across the patient.
        - Calculates entropy per chromosome (with/without error term and/or normalization), adds it to a dictionary, and returns it.

    Returns:
    --------
        copy_number_corrected_dict: Dictionary
            A dictionary with chromosome keys and entropy values.

    '''
    if file_ext is not None:
        # Finds seg file and opens it
        full_path = os.path.join(full_path, patient + file_ext)

    if 'deep_ichor' in sequencing_type:
        full_path = list(Path(full_path).rglob(f"{patient}/*_optimal/*.seg.txt"))[0]

    # Initialize an empty list to store the values from copy number corrected column
    copy_number_corrected_dict = {}

    # Sets the index for the metric and the chromosome id column (Titan has patient id as index 1)
    if metric_to_use == "Corrected_Copy_Number":
        if '_seg' in sequencing_type:
            n = 10
            subclone_status = 8
            bin_count_column = 4
            chr_coloumn = 1
        else:
            print("error")
            return
    elif metric_to_use == "logR":
        if sequencing_type == "deep":
            n = 12
            chr_coloumn = 1
        else:
            n = 5
            chr_coloumn = 0

    # Process sequencing types differently because Titan has dynamic column lengths, ichor is static
    if 'deep' in sequencing_type and 'ichor' not in sequencing_type:
        df = pd.read_csv(full_path, sep=r"\s+")

        for idx, row in df.iterrows():            
            patient_copy_number_corrected = row[metric_to_use]
            # Skips over header
            if metric_to_use == patient_copy_number_corrected:
                continue

            # # Filter out any conal clusters that have Cellular_Prevalence less then 0.8 
            # if row['Cellular_Prevalence'] is not None and row['Cellular_Prevalence'] < 0.8:
            #     continue

            chrom = row['Chromosome']
            copy_num = patient_copy_number_corrected
            bin_count = int(row['End']) - int(row['Start'])

            if '500kb' in sequencing_type:
                if bin_count < 500000:
                    continue

            if '10snp' in sequencing_type:
                bin_count = int(row['Length.snp.'])
                if bin_count < 10:
                    continue

            if 'entropy_mod' in sequencing_type:
                bin_count = 1

            # Create inner dictionary if chromosome does not exist
            copy_number_corrected_dict.setdefault(chrom, {})

            # If key already exists, add to existing value
            if copy_num in copy_number_corrected_dict[chrom]:
                copy_number_corrected_dict[chrom][copy_num] += int(bin_count)
            else:
                copy_number_corrected_dict[chrom][copy_num] = int(bin_count)

    else:
        file = open(full_path, "r")
        # Read the file line by line
        for line in file:
            patient_data = line.strip().split("\t")
            patient_copy_number_corrected = patient_data[n]
            # Skips over header
            if metric_to_use in patient_copy_number_corrected:
                continue
            # Append the value from copy number corrected column to the list
            # copy_number_corrected_dict.setdefault(patient_data[chr_coloumn], []).append(patient_copy_number_corrected)
            chrom = patient_data[chr_coloumn]
            copy_num = patient_copy_number_corrected
            bin_count = patient_data[bin_count_column]

            if 'entropy_mod' in sequencing_type:
                bin_count = 1

            if 'subclone' in sequencing_type:

                if patient_data[subclone_status] == "TRUE":
                    if int(copy_num) == 1:
                        copy_num = int(copy_num) + 0.5
                    elif int(copy_num) == 3:
                        copy_num = int(copy_num) - 0.5

            # Create inner dictionary if chromosome does not exist
            copy_number_corrected_dict.setdefault(chrom, {})

            # If key already exists, add to existing value
            if copy_num in copy_number_corrected_dict[chrom]:
                copy_number_corrected_dict[chrom][copy_num] += int(bin_count)
            else:
                copy_number_corrected_dict[chrom][copy_num] = int(bin_count)

        # Close the file
        file.close()

    # Sorts the dictionary by keys so it is in natural order
    def chr_key(chrom):
        """Extract number from chromosome name for natural sorting"""
        chrom_str = str(chrom).replace('chr', '').replace('Chr', '')
        try:
            return int(chrom_str)
        except ValueError:
            # For X, Y, M - put them at the end
            return float('inf')

    copy_number_corrected_dict = dict(sorted(copy_number_corrected_dict.items(), key=lambda x: chr_key(x[0])))
    
    # Check if every copy number value is 2 or not.
    # If value != 2 count it as a copy number alteration
    for chrom, values in copy_number_corrected_dict.items():

        # Calculates entropy per chromosome (with an error term for expected ploidy 2)
        len_of_chrom = sum(values.values())
        counts = Counter(values)
        counts_probs = {float(key): count / len_of_chrom for key, count in counts.items()}
        
        copy_number_corrected_dict[chrom] = entropy_equation_function(counts_probs)

        # Update max_K 
        max_K[0] = len(counts_probs) if len(counts_probs) > max_K[0] else max_K[0]
        
    # Return genomic instability per chrom
    return copy_number_corrected_dict

def process_entropy_per_chromosome(genomic_instability_directory_array, tumor_fraction_directory_array, patient_dict, metric_to_use, sequencing_type, bin_range=None, entropy_equation_function=None):
    '''
    Parameters:
    -----------
        genomic_instability_directory_array: List
            Array holding genomic instability directory paths.

        tumor_fraction_directory_array: List
            Array holding tumor fraction directory paths.

        patient_dict: Dictionary
            Dictionary to store patient data with patient identifiers as keys.

        metric_to_use: String
            Parameter to determine which metric_to_use from cna.seg file to use for entropy.

        sequencing_type: String
            Determins what sequencying type to select for (e.g. 'deep', 'ulp', or 'ulp_curated').

        bin_range: int (default = None)
            Used to set a range on how many bins to select in deep data to manage entropy score (not used for ulp).

    Function:
    ---------
        - Loops through each directory in the genomic instability directory array.
        - For each patient found, extracts entropy per chromosome data using find_entropy_per_chromosome.
        - Reads tumor fraction data from the tumor fraction file.
        - Updates the patient dictionary with tumor fraction and entropy values.

    Returns:
    --------
        None (modifies patient_dict in place)

    '''
    print("Processing FGA_TFX...")

    # Max number of copy number states found in the cohort
    max_K = [0]

    # Tumor fraction data
    tumor_fraction_data = open(tumor_fraction_directory_array[0], "r")

    if sequencing_type == "ulp_curated":
        curated_files_df = pd.read_csv(genomic_instability_directory_array[0])
        for index in curated_files_df.index:
            full_path = curated_files_df.iloc[index]["curated_solution_path"]
            patient = curated_files_df.iloc[index]["Identifier"]
            patient_dict[patient] = (None, find_entropy_per_chromosome(full_path, patient, metric_to_use, sequencing_type, bin_range, entropy_equation_function, None), None)
    elif "ulp_curated_seg" in sequencing_type:
        curated_files_df = pd.read_csv(genomic_instability_directory_array[0])
        for index in curated_files_df.index:
            full_path = curated_files_df.iloc[index]["curated_solution_path"].replace('.cna','') + '.txt'
            patient = curated_files_df.iloc[index]["Identifier"]
            patient_dict[patient] = (None, find_entropy_per_chromosome_seg_file(full_path, patient, metric_to_use, sequencing_type, bin_range, entropy_equation_function, max_K, None), None)
    else:
        # For every directory in the array find the patients ulp data
        for directory in genomic_instability_directory_array:
            # For loop to find each patient in each batch
            for patient in os.listdir(directory):
                # Creates a full path with directory and patient identifier
                # If patient is not in the patient dictionary add the patient, else print that patient is already in the list
                if patient not in patient_dict:
                    if "deep" in sequencing_type and 'ichor' not in sequencing_type:
                        if os.path.isdir(os.path.join(directory, patient)):
                            full_path = directory
                            if os.path.isdir(full_path):
                                if "_seg" in sequencing_type:
                                    patient_dict[patient] = (None, find_entropy_per_chromosome_seg_file(full_path, patient, metric_to_use, sequencing_type, bin_range, entropy_equation_function, max_K, file_ext=".titan.ichor.seg.txt"), None)
                                else:
                                    patient_dict[patient] = (None, find_entropy_per_chromosome(full_path, patient, metric_to_use, sequencing_type, bin_range, max_K, file_ext=".titan.ichor.seg.txt"), None)
                    
                    else:
                        full_path = os.path.join(directory, patient)
                        # If full path is true directory append patient identifier and genome altered to different lists
                        if os.path.isdir(full_path):
                            if "Batch5" in full_path:
                                for dir in os.listdir(full_path):
                                    if "solution" in dir:
                                        full_path = os.path.join(full_path, dir)
                            # Creates a dictionary of patient names with values genome altered and None
                            number = patient.split('C')[-1]
                            if "_seg" in sequencing_type:
                                patient_dict[patient] = (None, find_entropy_per_chromosome_seg_file(full_path, patient, metric_to_use, sequencing_type, bin_range, entropy_equation_function, max_K, None), None)
                            else:
                                patient_dict[patient] = (None, find_entropy_per_chromosome(full_path, patient, metric_to_use, sequencing_type, None, entropy_equation_function, max_K, None), number)

                else:
                    print(f"{patient} already in list")


    # Creates an array for each tumor fraction for each patient 
    for line in tumor_fraction_data:
        patient_data = line.strip().split("\t")
        if patient_data[1] == "tfx": # Skips the first row
            continue
        # Update None value for tumor fraction for each patient id
        patient_id = patient_data[0]
        tumor_fraction = patient_data[1]
        if patient_id in patient_dict:
            patient_dict[patient_id] = (float(tumor_fraction), patient_dict[patient_id][1], patient_dict[patient_id][2])

    return max_K[0]

def create_entropy_data_table_per_chromosome(patient_dict, output_sub_dir, csv_file_name, output_path, metric_to_use, sequencing_type):
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
    
    Function:
    ---------
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
            if "_" in patient_id:
                patient_id_specific = patient_id[:patient_id.find("_")]
            else:
                patient_id_specific = patient_id

            # Remove trailing cluster value
            if 'deep' in sequencing_type and 'ichor' not in sequencing_type:
                patient_id = patient_id[:patient_id.rfind("_")] 
            
            patient_time_point = patient_id[patient_id.rfind("_")+1:]
            list_of_entropies = [value[1] for value in values]
            csv_writer.writerow([patient_id_specific, patient_time_point, *list_of_entropies])

def create_entropy_data_table_per_chromosome_v2(patient_dict, output_sub_dir, csv_file_name, output_path, metric_to_use, sequencing_type):
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
    
    Function:
    ---------
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
    for patient_id, (genomic_instability) in patient_dict.items():
        for chrom, value in genomic_instability.items():
            if patient_id not in dict_of_patients.keys():
                dict_of_patients[patient_id] = [(chrom, value)]
            else:
                dict_of_patients[patient_id].append((chrom, value))

    with open(csvFileName, 'w', ) as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        csv_writer.writerow(tsv_header_name)

        for patient_id, values in dict_of_patients.items():
            if "_" in patient_id:
                patient_id_specific = patient_id[:patient_id.find("_")]
            else:
                patient_id_specific = patient_id

            # Remove trailing cluster value
            if 'deep' in sequencing_type and 'ichor' not in sequencing_type:
                patient_id = patient_id[:patient_id.rfind("_")] 
            
            patient_time_point = patient_id[patient_id.rfind("_")+1:]
            list_of_entropies = [value[1] for value in values]
            csv_writer.writerow([patient_id_specific, patient_time_point, *list_of_entropies])

def find_T_cycles(patient_dict):
    '''
    Parameters:
    -----------
        patient_dict: Dictionary
            Dictionary with patient identifiers as keys and tuple values (tumor_fraction, entropy_dict, cycle).
    
    Function:
    ---------
        - Sorts the patient dictionary by patient identifier keys.
        - Iterates through the sorted dictionary to find the maximum treatment cycle for each patient.
        - Creates a new dictionary mapping patient IDs to their maximum cycle number.

    Returns:
    --------
        patient_progression_cycle_dict: Dictionary
            Dictionary with patient IDs as keys and their maximum cycle numbers as values.

    '''
    # Sort dict by key
    patient_dict = dict(sorted(patient_dict.items()))

    # Sort dict by index 2 of tuple value
    {k: v for k, v in sorted(patient_dict.items(), key = lambda item: item[1][2])}
    
    patient_pregression_cycle_dict = {}
    last_patient = None
    for item in patient_dict.keys():
        if last_patient != item:
            patient_id = item.split("_")[0]
            patient_pregression_cycle_dict[patient_id] = patient_dict[item][2]
        last_patient = item

    return patient_pregression_cycle_dict
        
def create_distribuition_plot_by_chromosomes(output_path, data_table_dir, output_sub_dir_name, metric_to_use, sequencing_type, plot_titles, csv_file_name, cycle):
    '''
    Parameters:
    -----------
        output_path: String
            Path to the output directory.

        data_table_dir: String
            Path to directory with data table.

        output_sub_dir_name: String
            Name of the subdirectory to store plots.

        plot_titles: String
            Title for the distribution plots.

        csv_file_name: String
            Name of the CSV file containing entropy data.

        cycle: String
            Treatment cycle to filter data (e.g., 'C1', 'C2').
    
    Function:
    ---------
        - Reads entropy data from the specified CSV file.
        - Filters the data for the specified treatment cycle.
        - Calculates appropriate bin ranges based on maximum values.
        - Creates and saves distribution plots for each chromosome.

    Returns:
    --------
        None (creates and saves histogram plots)

    '''
    # Create directory
    file_dir = os.path.join(output_path, f'data-tables/{data_table_dir}')
    os.makedirs(file_dir, exist_ok=True)

    # Sets output path and name    
    csvFileName = os.path.join(file_dir, f'{csv_file_name}')

    df = pd.read_csv(csvFileName)

    # Find top value of sorted column for chromosome 8 for bin range
    sorted_df = sort_df_by_column(df, "chr8", False)
    top_value = sorted_df.loc[0, "chr8"]

    # Remove all other cycles
    # mask = df["cycle"] == cycle
    # df_filtered = df[mask]
    if cycle is not None:
        df_filtered = df[df["cycle"].str.contains(cycle, na=False)]
    else:
        df_filtered = df

    out_dir = os.path.join(output_path, f"entropy-analysis-{sequencing_type}-{metric_to_use}/{output_sub_dir_name}")
    os.makedirs(out_dir, exist_ok=True)

    # Create bin range (for top_value < 10 it will increment by .5 instead of 1)
    if top_value < 5:
        bin_range = [x * .25 for x in range(0, math.ceil(top_value)*4)]
    elif top_value < 10:
        bin_range = [x * .5 for x in range(0, math.ceil(top_value)*2)]
    else:
        if top_value <= 100:
            bin_range = list(range(0, math.ceil(top_value)))

    # Create histogram for each chromosome
    for column in df_filtered.columns:
        if column == "patient_id" or column == "cycle":
            continue

        if "bin_range" in locals():
            distribution_plot_from_df_column(
                df_filtered, 
                column, 
                plot_titles,
                "Entropy Value",
                "Number of Patients",
                cycle, 
                bin_range, 
                out_dir,
                column
            )

        publication_violin_plot(
            df_filtered, 
            column, 
            plot_titles,
            "Chromosome",
            "Entropy Value",
            cycle, 
            out_dir,
            column
        )
    
    # Create a multi-violin plot 
    chrom_columns = [col for col in df_filtered.columns if col not in ["patient_id", "cycle"]]

    publication_multi_violin_plot(
        df_filtered,
        chrom_columns,
        plot_titles,
        "Chromosome",
        "Entropy Value",
        cycle,
        out_dir,
        file_name=f"all_chromosomes_violin_{cycle}"
    )

def distribution_plot_from_df_column(df, column, title, xlabel, ylabel, cycle, bin_list, out_dir, file_name):
    '''
    Parameters:
    -----------
        df: DataFrame
            Pandas DataFrame containing the data to plot.

        column: String
            Column name from the DataFrame to create histogram from.

        title: String
            Title for the histogram plot.

        xlabel: String
            Label for the x-axis.

        ylabel: String
            Label for the y-axis.

        cycle: String
            Treatment cycle identifier for the plot.

        bin_list: List
            List of bin edges for the histogram.

        out_dir: String
            Output directory path for saving the plot.

        file_name: String
            Name of the output plot file.
    
    Function:
    ---------
        - Creates a histogram from the specified DataFrame column.
        - Applies formatting with black edges and specified color.
        - Adds grid lines and labels for clarity.
        - Saves the plot as a PNG file to the specified output directory.

    Returns:
    --------
        None (creates and saves histogram PNG file)

    '''

    plt.figure(figsize=(8, 6)) # Optional: Adjust plot size

    data = df[column]

    counts, bin, patches = plt.hist(data, bins=bin_list, edgecolor="black", color="#4C72B0", alpha=0.8)

    plt.xticks(bin)

    # Add labels and title for clarity
    plt.xlabel(f"{xlabel}")
    plt.ylabel(f"{ylabel}")
    plt.title(f"{title} ({cycle}) for {column.title()} ")

    plt.grid(axis="y", alpha=0.3)

    plt.tight_layout()

    plot_file_name = os.path.join(out_dir, f'{file_name}_histogram.png')

    plt.savefig(plot_file_name, dpi=300)
    plt.close()

def publication_violin_plot(df, column, title=None, xlabel=None, ylabel=None, cycle=None, out_dir=".", file_name="violin_plot",):
    """
    Create a publication-ready violin plot with:
      - Median
      - Quartiles (25th, 75th percentiles)
      - Whiskers (1.5 * IQR)
    
    Parameters
    ----------
    df : pandas.DataFrame
    column : str
        Column name to plot
    title : str, optional
    xlabel : str, optional
    ylabel : str, optional
    cycle : str/int, optional
        Added to title if provided
    out_dir : str
    file_name : str
    """

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Drop NA values
    data = df[column].dropna().values

    if len(data) == 0:
        raise ValueError(f"No valid data found in column '{column}'")

    # Compute statistics
    q1 = np.percentile(data, 25)
    median = np.percentile(data, 50)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1

    lower_whisker = np.min(data[data >= q1 - 1.5 * iqr])
    upper_whisker = np.max(data[data <= q3 + 1.5 * iqr])

    # Figure setup (journal style)
    plt.rcParams.update({
        "font.size": 12,
        "font.family": "sans-serif",
        "axes.linewidth": 1.2,
        "pdf.fonttype": 42,   # editable text in Illustrator
        "ps.fonttype": 42
    })

    fig, ax = plt.subplots(figsize=(4, 6))

    # Violin
    parts = ax.violinplot(
        data,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )

    # Style violin
    for pc in parts['bodies']:
        pc.set_facecolor("#82CAFF")
        pc.set_edgecolor("black")
        pc.set_alpha(0.8)
        pc.set_linewidth(1.2)

    # Add quartile box
    ax.vlines(1, q1, q3, color='black', linewidth=6)

    # Add median
    ax.scatter(1, median, color='white', edgecolor='black', zorder=3, s=60)

    # Add whiskers
    ax.vlines(1, lower_whisker, upper_whisker, color='black', linewidth=1.5)

    # Whisker caps
    ax.hlines(lower_whisker, 0.95, 1.05, color='black', linewidth=1.5)
    ax.hlines(upper_whisker, 0.95, 1.05, color='black', linewidth=1.5)

    # Formatting
    ax.set_xticks([1])
    ax.set_xticklabels([column])

    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)

    if title:
        if cycle is not None:
            ax.set_title(f"{title} (Cycle {cycle})")
        else:
            ax.set_title(title)

    # Remove top/right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    # Save
    png_path = os.path.join(out_dir, f"{file_name}_violin_plot.png")

    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.close()

def publication_multi_violin_plot(df, columns, title=None, xlabel=None, ylabel=None, cycle=None, out_dir=".", file_name="multi_violin_plot", add_significance=True):

    os.makedirs(out_dir, exist_ok=True)

    # Collect data
    data = [df[col].dropna().values for col in columns]

    # Remove empty columns
    valid = [(col, d) for col, d in zip(columns, data) if len(d) > 0]
    columns = [v[0] for v in valid]
    data = [v[1] for v in valid]

    if len(data) == 0:
        raise ValueError("No valid columns to plot.")

    plt.rcParams.update({
        "font.size": 11,
        "axes.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42
    })

    fig, ax = plt.subplots(figsize=(max(8, len(columns)*0.6), 6))

    parts = ax.violinplot(
        data,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )

    # Style violins
    for pc in parts['bodies']:
        pc.set_facecolor('#82CAFF')
        pc.set_edgecolor("black")
        pc.set_alpha(0.80)
        pc.set_linewidth(1.2)

    medians = {}

    # Add quartiles, medians, whiskers
    for i, d in enumerate(data, start=1):

        q1 = np.percentile(d, 25)
        median = np.percentile(d, 50)
        q3 = np.percentile(d, 75)
        iqr = q3 - q1

        medians[columns[i-1]] = median

        lower = np.min(d[d >= q1 - 1.5 * iqr])
        upper = np.max(d[d <= q3 + 1.5 * iqr])

        ax.vlines(i, q1, q3, color='black', linewidth=4)
        ax.scatter(i, median, color='white', edgecolor='black', s=40, zorder=3)
        ax.vlines(i, lower, upper, color='black', linewidth=1.2)
        ax.hlines(lower, i-0.15, i+0.15, color='black', linewidth=1.2)
        ax.hlines(upper, i-0.15, i+0.15, color='black', linewidth=1.2)

    # ---------------------------
    # Add significance vs highest
    # ---------------------------
    if add_significance and len(columns) > 1:

        # Identify highest median chromosome
        top_chr = max(medians, key=medians.get)
        top_index = columns.index(top_chr)
        top_data = data[top_index]

        p_values = []
        comparison_indices = []

        for i, col in enumerate(columns):
            if col == top_chr:
                continue

            stat, p = mannwhitneyu(
                top_data,
                data[i],
                alternative='two-sided'  # one-sided (top > others)
            )
            p_values.append(p)
            comparison_indices.append(i)

        # FDR correction
        reject, pvals_corrected, _, _ = multipletests(
            p_values,
            method='fdr_bh'
        )

        # Add stars directly above each significant violin
        for idx, sig in zip(comparison_indices, reject):

            if not sig:
                continue

            d = data[idx]
            local_max = max(d)

            offset = (max(d) - min(d)) * 0.08  # small vertical offset
            if offset == 0:
                offset = 0.05  # fallback for flat distributions

            ax.text(
                idx + 1,
                local_max + offset,
                "*",
                ha='center',
                va='bottom',
                fontsize=14,
                fontweight='bold'
            )

    ax.set_xticks(range(1, len(columns) + 1))
    labels = ax.set_xticklabels(columns, rotation=45, ha="right")

    # Bold reference chromosome
    labels[top_index].set_fontweight("bold")

    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)

    if title:
        if cycle is not None:
            ax.set_title(f"{title} (Cycle {cycle})", fontweight="bold")
        else:
            ax.set_title(title, fontweight="bold")

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    png_path = os.path.join(out_dir, file_name)

    plt.savefig(f"{png_path}.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"{png_path}.pdf", dpi=600, bbox_inches="tight")
    plt.close()

def line_plot_across_cycles(output_path, sub_dir, metric_to_use, sequencing_type, entropy_col, stratify_col, split_method, plot_title, file_name, df_combined, cycle_col, cycles):
    '''
    Parameters:
    -----------
        output_path: String
            Path to the output directory.

        sub_dir: String
            Name of the subdirectory to store plots.

        metric_to_use: String
            Parameter to show what metric from cna.seg file was used for entropy.

        sequencing_type: String
            Type of sequencing data (e.g., 'deep', 'ulp_curated').

        entropy_col: String
            Name of chromosome column in entorpy file.

        stratify_col: String
            Name of column to split data on.

        split_method: String
            Method to split stratify_col by.

        plot_title: String
            Title for the distribution plots.

        file_name: String
            File name (no file extention).

        df_combined: DataFrame
            DataFrame that has entropy scores and data to split on.

        cycle_col: String
            Name of cycle column in csv_file. 

        cycles: List
            Treatment cycles to show on plot (e.g., 'C1', 'C2').
    
    Function:
    ---------
        - Create a line plot of the median entropy value across each cycle for the all patients, split between split_method.

    Returns:
    --------
        None (creates and saves plot)

    '''

    df_stratified, group_labels, group_n  = stratify_data(df_combined, stratify_col, split_method, None)

    df_stratified[cycle_col] = pd.Categorical(df_stratified[cycle_col], categories=cycles, ordered=True)
    
    # group + average
    df_stats = (
        df_stratified.groupby(['group', cycle_col], observed=True)[entropy_col]
        .agg(['median', 'sem'])
        .reset_index()
    )

    # Map group labels (optional but recommended)
    # assumes stratify_data already returned this
    df_stats['group_label'] = df_stats['group'].map(group_labels)

    # Pivot for cleaner plotting
    median_df = df_stats.pivot(index=cycle_col, columns='group_label', values='median')
    sem_df = df_stats.pivot(index=cycle_col, columns='group_label', values='sem')

    # --- Plot ---
    plt.figure(figsize=(5.5, 3.8), dpi=300)

    colors = {
        'Favorable': '#1f77b4',
        'Unfavorable': '#d62728'
    }

    for group, color in colors.items():
        y = median_df[group]
        ci = 1.96 * sem_df[group]

        # Clip lower bound at 0
        lower = np.maximum(y - ci, 0)
        upper = y + ci

        # Convert to asymmetric error bars
        yerr = np.vstack([y - lower, upper - y])

        plt.errorbar(
            median_df.index,
            y,
            yerr=yerr,
            fmt='o-',            # line + marker
            linewidth=1.5,
            markersize=4,
            capsize=3,           # little horizontal caps
            capthick=1,
            elinewidth=1,
            color=color,
            label=group
        )

    # Labels (slightly larger)
    plt.xlabel("Treatment Cycle", fontsize=12)
    plt.ylabel(f"Median {entropy_col}", fontsize=12)
    plt.title(plot_title, fontsize=12)

    # Y-axis formatting (0.25 increments)
    ax = plt.gca()
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(MultipleLocator(0.125))

    # Optional: cleaner tick labels
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # Grid (subtle)
    plt.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.grid(True, which='minor', linestyle=':', linewidth=0.3, alpha=0.3)

    # Remove top/right spines (publication style)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Slightly thicker left/bottom spines
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)

    # Tick params
    ax.tick_params(axis='both', which='major', labelsize=10, length=4)
    ax.tick_params(axis='both', which='minor', length=2)

    # Legend (clean + smaller)
    plt.legend(
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        frameon=False,
        fontsize=9,
        handlelength=1.5
    )

    plt.tight_layout(rect=[0, 0, 0.82, 1])

    # Save
    out_dir = os.path.join(output_path, f"entropy-analysis-{sequencing_type}-{metric_to_use}/{sub_dir}")
    os.makedirs(out_dir, exist_ok=True)

    save_path = os.path.join(out_dir, f"{file_name}.png")

    plt.savefig(save_path, dpi=300)
    plt.close()

def bar_plot_comparing_cycles(output_path, sub_dir, metric_to_use, sequencing_type, entropy_col, id_col, stratify_col, split_method, plot_title, file_name, df, cycle_col, cycles, direction_cuttoff):
    '''
    Parameters:
    -----------
        output_path: String
            Path to the output directory.

        sub_dir: String
            Name of the subdirectory to store plots.

        metric_to_use: String
            Parameter to show what metric from cna.seg file was used for entropy.

        sequencing_type: String
            Type of sequencing data (e.g., 'deep', 'ulp_curated').

        entropy_col: String
            Name of chromosome column in dataframe.

        id_col: String
            Name of id column in dataframe.

        stratify_col: String
            Name of column to split data on.

        split_method: String
            Method to split stratify_col by.

        plot_title: String
            Title for the distribution plots.

        file_name: String
            File name (no file extention).

        df: DataFrame
            DataFrame that has entropy scores and data to split on.

        cycle_col: String
            Name of cycle column in csv_file. 

        cycles: List
            Treatment cycles to show on plot (e.g., 'C1', 'C2').
    
    Function:
    ---------
        - Create a line plot of the median entropy value across each cycle for the all patients, split between split_method.

    Returns:
    --------
        None (creates and saves plot)

    '''

    # Create grouping column
    df_stratified, group_labels, group_n  = stratify_data(df, stratify_col, split_method, None)

    # Add groupings
    df_groups = df_stratified[df_stratified['cycle'] == cycles[0]][[id_col, 'group']]
    
    # Transform df to work for plotting
    df_cycles = (
        df[df[cycle_col].isin(cycles)]
        .pivot(
            index=id_col,
            columns=cycle_col,
            values=entropy_col
        ).reset_index()
    )

    # Add groupings
    df_cycles = df_cycles.merge(df_groups, on=id_col)

    # Remove rows with NaNs
    df_cycles = df_cycles.dropna()

    if direction_cuttoff == "both":
        df_cycles = df_cycles
    elif direction_cuttoff == "only_increasing":
        df_cycles = df_cycles[df_cycles[cycles[1]] > df_cycles[cycles[0]]]
    elif direction_cuttoff == "only_decreasing":
        df_cycles = df_cycles[df_cycles[cycles[0]] > df_cycles[cycles[1]]]
    elif direction_cuttoff == "only_the_same":
        df_cycles = df_cycles[df_cycles[cycles[0]] == df_cycles[cycles[1]]]

    file_name = file_name + "_" + direction_cuttoff

    # ----------------
    # Plot
    # ----------------
    
    # --- Sort by top cycle (C1) ---
    df_cycles = df_cycles.sort_values(by=cycles[0], ascending=False).reset_index(drop=True)

    x = np.arange(len(df_cycles))

    fig, ax = plt.subplots(figsize=(7, 4), dpi=300)

    # Colors for cycles (more muted / publication friendly)
    c1_color = '#357A8E'
    c2_color = '#4D2A4D'

    group_colors = {
        0: '#3375B6',   # or "Favorable"
        1: '#BC4024'    # or "Unfavorable"
    }

    # --- seperate columns --
    c1_vals = df_cycles[cycles[0]].abs()
    c2_vals = df_cycles[cycles[1]].abs()


    # --- Bars (touching, with outline) ---
    ax.bar(
        x,
        c1_vals,
        width=1.0,
        color=c1_color,
        edgecolor='black',
        linewidth=0.3
    )

    ax.bar(
        x,
        -c2_vals,
        width=1.0,
        color=c2_color,
        edgecolor='black',
        linewidth=0.3
    )

    # --- Center line ---
    ax.axhline(0, color='black', linewidth=1)

    # --- Remove x-axis completely ---
    ax.set_xticks([])
    ax.set_xlabel('')
    ax.set_xlim(-0.5, len(df_cycles) - 0.5) # Remove padding between axis and first bar

    # --- Y-axis formatting ---
    ax.set_ylabel(entropy_col, fontsize=11)

    # Make y ticks symmetric and clean
    ymax = max(df_cycles[cycles[0]].max(), df_cycles[cycles[1]].max())
    if not np.isfinite(ymax) or ymax == 0:
        ymax = 1e-6
    ax.set_ylim(-ymax * 1.1, ymax * 1.1)

    # Show absolute values on y-axis
    yticks = ax.get_yticks()
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{abs(x):.2f}"))

    # --- Clean spines ---
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Slightly thicker left spine
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_visible(False)

    # --- Keep y grid lines, remove x grid lines ---
    ax.yaxis.grid(True, linestyle='--', linewidth=0.5, color='grey', alpha=0.5)
    ax.xaxis.grid(False)  # keep x-grid off

    # --- Direct labels instead of legend ---
    ax.text(
        0.5, 0.95,
        f'Cycle {cycles[0][1]}',
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment='top',
        color='grey'
    )

    ax.text(
        0.50, 0.05,
        f'Cycle {cycles[1][1]}',
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment='bottom',
        color='grey'
    )

    # --- Add sample number ---
    ax.text(
        0.50, 0.01,
        f'(n = {len(df_cycles)})',
        transform=ax.transAxes,
        fontsize=6,
        verticalalignment='bottom',
        color='grey'
    )

    # --- Group annotation strip ---
    strip_height = ymax * 0.035
    strip_y = ymax * 1.01

    for i, g in enumerate(df_cycles['group']):
        ax.add_patch(
            plt.Rectangle(
                (i - 0.5, strip_y),
                1.0,
                strip_height,
                color=group_colors[g],
                linewidth=0
            )
        )

    # --- Legend ---
    legend_elements = [
        Patch(facecolor=group_colors[0], label='Favorable'),
        Patch(facecolor=group_colors[1], label='Unfavorable')
    ]

    ax.legend(
        handles=legend_elements,
        loc='upper left',
        bbox_to_anchor=(1, 1),
        frameon=False,
        fontsize=6,
        title='Group',
        title_fontsize=6
    )

    # --- Title ---
    ax.set_title(plot_title, fontsize=12)

    plt.tight_layout()

    # Save
    out_dir = os.path.join(output_path, f"entropy-analysis-{sequencing_type}-{metric_to_use}/{sub_dir}")
    os.makedirs(out_dir, exist_ok=True)

    save_path = os.path.join(out_dir, f"{file_name}.png")

    plt.savefig(save_path, dpi=300)
    plt.close()

def calculate_fold_change_across_groups_per_patient(value_col, id_col, group_col, groups, df):
    '''
    Parameters:
    -----------
        output_path: String
            Path to the output directory.
        sub_dir: String
            Name of the subdirectory to store plots.
        metric_to_use: String
            Parameter to show what metric from cna.seg file was used for entropy.
        sequencing_type: String
            Type of sequencing data (e.g., 'deep', 'ulp_curated').
        value_col: String
            Name of column in dataframe that contains value to calculate fold change from.
        id_col: String
            Name of id column in dataframe.
        group_col: String
            Column indicating the condition, timepoint, or group within each patient used to compute fold change (e.g., baseline vs treatment, cycle number).
        groups: List
            Group names in group_col to build ratio (Needs 2 values; fist value is numerator, second value is denominator).
        plot_title: String
            Title for the plot.
        file_name: String
            File name (no file extention).
        df: pandas DataFrame
            DataFrame that has parameters as column names.

    Function:
    ---------
        - Creates a dataframe that has the log2 fold change (with pseudocount of 1) of each patient based on groups.

    Returns:
    --------
        fc_df: pandas DataFrame
            A DataFrame with patient id and fold change column.

    '''
    # Pivot to create columns for each group
    df_groups = df.pivot(
        index=id_col,
        columns=group_col,
        values=value_col
    ).reset_index()

    # Filter to groups parameter
    df_filtered = df_groups[[id_col] + groups]

    # Drop rows with NA
    fc_df = df_filtered.dropna().copy()

    fc_df['log2_fc'] = np.log2((fc_df[groups[0]] + 1) / (fc_df[groups[1]] + 1))

    return fc_df

def add_max_change_columns(df, cols):
    """
    Adds two columns to the DataFrame:
    - max_abs_change: largest absolute difference between any two columns
    - max_pct_change: percent change corresponding to that pair

    Parameters
    ----------
    df : pd.DataFrame
    cols : list of str
        Columns to compare

    Returns
    -------
    pd.DataFrame
    """

    df = df.copy()

    max_diff = np.full(len(df), 0.0)
    max_ratio = np.full(len(df), np.nan)

    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):            
            # Grabs correct rows 
            c1 = cols[i]
            c2 = cols[j]
            v1 = df[c1].values
            v2 = df[c2].values

            # valid (non-NaN) mask - only keps rows where both values exist
            mask = ~np.isnan(v1) & ~np.isnan(v2)

            # Compute difference
            diff = v2 - v1

            # percent change only where v1 != 0
            pct = np.full(len(df), np.nan)
            valid_div = (v1 != 0) & mask
            pct[valid_div] = diff[valid_div] / v1[valid_div]

            # update where this pair's difference is larger than the saved difference
            update_mask = (np.abs(diff) > np.abs(max_diff)) & mask

            # Update stored results
            max_diff[update_mask] = diff[update_mask]
            max_ratio[update_mask] = pct[update_mask]

    df['pct_change'] = max_ratio
    df['pct_change'] = df['pct_change'] * 100

    return df

def create_waterfall_plot(df, value_col, group_col, title=None, xlabel=None, ylabel=None, out_dir=".", file_name="waterfall_plot"):
    groups = sorted(df[group_col].dropna().unique())

    # Assign colors
    group_colors_map = {
        0: '#3375B6',   # Favorable
        1: '#BC4024'    # Unfavorable
    }
    group_colors = {group: group_colors_map[i] for i, group in enumerate(groups)}

    fig, axes = plt.subplots(1, len(groups), figsize=(6 * len(groups), 5), sharey=False)

    if len(groups) == 1:
        axes = [axes]

    for i, (ax, group) in enumerate(zip(axes, groups)):
        df_sub = df[df[group_col] == group][[value_col]].dropna().copy()

        # Sort within group
        df_sub = df_sub.sort_values(by=value_col, ascending=False).reset_index(drop=True)

        values = df_sub[value_col].values
        x = np.arange(len(values))

        # Draw bars individually to control transparency
        for xi, v in zip(x, values):
            color = group_colors[group]
            alpha = 0.6 if v >= 0 else 1.0  # translucent above 0

            ax.bar(xi, v, color=color, alpha=alpha)

        # Zero line (black)
        ax.axhline(0, color='black', linewidth=1.2)

        # Make bars touch edges
        ax.set_xlim(-0.5, len(values) - 0.5)

        # Find how many patients in group to add to x label
        n = len(values)

        # Titles and labels
        ax.set_title(f"{group}")
        ax.set_xlabel(f"{xlabel} (n = {n})")

        if i == 0:
            ax.set_ylabel(ylabel)

        ax.set_xticks([])

        # Set y-axis from 100 (top) to -100 (bottom)
        ax.set_ylim(-100, 100)

    fig.suptitle(title, fontsize=16)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    os.makedirs(out_dir, exist_ok=True)

    save_path = os.path.join(out_dir, f"{file_name}.png")

    plt.savefig(save_path, dpi=300)
    plt.close()

def distribution_plot_comparing_entorpy_equations(output_path, all_data, data_labels, column, sequencing_type, metric_to_use, file_name, labels=None, figsize=(6,4), linewidth=2, dpi=300):
    """
    Plot kernel density distributions for the same column across multiple CSV files.

    Parameters
    ----------
    all_data: list
        List of column data.
    data_labels : list
        List of CSV file paths
    column : str
        Column name to plot
    labels : list, optional
        Labels for each dataset (defaults to file names)
    figsize : tuple
        Figure size
    linewidth : int
        Width of distribution curves
    dpi : int
        Figure resolution

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """

    if labels is None:
        labels = [f.split("/")[-1].replace("_per_chr_table.csv","") for f in data_labels]

    plt.rcParams.update({
        "font.size": 12,
        "axes.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42
    })

    # ---------------------------------------------------------
    # Figure
    # ---------------------------------------------------------

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # ---------------------------------------------------------
    # Shared x range
    # ---------------------------------------------------------

    xmin = min(np.min(d) for d in all_data)
    xmax = max(np.max(d) for d in all_data)

    padding = (xmax - xmin) * 0.05

    x = np.linspace(xmin - padding, xmax + padding, 1000)

    # ---------------------------------------------------------
    # Nature-style color palette
    # ---------------------------------------------------------

    colors = [
        "#000000",
        "#D55E00",
        "#0072B2",
        "#009E73",
        "#CC79A7",
        "#E69F00",
        "#56B4E9"
    ]

    # ---------------------------------------------------------
    # Plot KDEs
    # ---------------------------------------------------------

    for i, (data, label) in enumerate(zip(all_data, labels)):

        data = np.asarray(data)
        data = data[~np.isnan(data)]

        color = colors[i % len(colors)]

        # KDE
        kde = gaussian_kde(data)
        y = kde(x)

        # Main density curve
        ax.plot(
            x,
            y,
            color=color,
            linewidth=linewidth,
            label=f"{label}  (n={len(data)})",
            zorder=3
        )

        # Light fill under curve
        ax.fill_between(
            x,
            y,
            color=color,
            alpha=0.12,
            zorder=2
        )

        # -------------------------------------------------
        # Statistics
        # -------------------------------------------------

        median = np.median(data)

        # Median
        ax.axvline(
            median,
            color=color,
            linestyle="--",
            linewidth=1.5,
            alpha=0.9,
            zorder=1
        )

    # ---------------------------------------------------------
    # Labels
    # ---------------------------------------------------------

    pretty_column = column.replace("_", " ").title()

    ax.set_xlabel(pretty_column)
    ax.set_ylabel("Density")

    # ---------------------------------------------------------
    # Clean Nature-style aesthetics
    # ---------------------------------------------------------

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    ax.tick_params(
        axis="both",
        which="major",
        direction="out"
    )

    # Subtle grid
    ax.grid(
        axis="y",
        linestyle="-",
        linewidth=0.5,
        alpha=0.15
    )

    # ---------------------------------------------------------
    # Legend
    # ---------------------------------------------------------

    stat_legend = [
        Line2D(
            [0],
            [0],
            color="black",
            linestyle="--",
            linewidth=1.5,
            label="Median"
        )
    ]

    handles, existing_labels = ax.get_legend_handles_labels()

    ax.legend(
        handles + stat_legend,
        [f"Median: {median:.3f}"],
        frameon=False,
        loc="upper left",
        # bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
        handlelength=2.5
    )

    # ---------------------------------------------------------
    # Layout
    # ---------------------------------------------------------

    plt.tight_layout()


    # out_dir = os.path.join(output_path, f"entropy-analysis-{sequencing_type}-{metric_to_use}/{sub_dir_name}")
    # os.makedirs(out_dir, exist_ok=True)
    os.makedirs(output_path, exist_ok=True)

    png_path = os.path.join(output_path, f"{file_name}.png")

    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.close()

    print("distribution plot comparing entorpy equations saved")

def sort_df_by_column(df, column, ascending):
    '''
    Parameters:
    -----------
        df: DataFrame
            Pandas DataFrame to sort.

        column: String
            Column name to sort by.

        ascending: Boolean
            True to sort in ascending order, False for descending order.
    
    Function:
    ---------
        - Sorts the DataFrame by the specified column.
        - Resets the DataFrame index after sorting.
        - Extracts the first value from the sorted column.

    Returns:
    --------
        value: Float or Placeholder
            The first (top) value from the sorted column.

    '''
    sorted_df = df.sort_values(by=column, ascending=ascending)
    sorted_df = sorted_df.reset_index()
    return sorted_df


def qualitative_top_and_bottom_40_patients(output_path, csv_path_1, csv_title_1, csv_path_2, csv_title_2, metric_to_use, sequencing_type, column_to_sort_by, cycle, n, mask=False):
    '''
    Parameters:
    -----------
        output_path: String
            Output directory path for saving the plot.

        csv_path_1: String
            Path to csv.

        csv_path_2: String
            Path to csv.

        metric_to_use: String
            Parameter to show what metric from cna.seg file was used for entropy.

        column_to_sort_by: String
            Column passed to sort_df_by_column.

        cycle: String
            Treatment cycle to filter data (e.g., 'C1', 'C2').

        n: int
            Number for top and bottom range.

        mask: Boolean
            Parameter to determine if a mask should be applied to dataframes.
    
    Function:
    ---------
        - Converts csv files to dataframes and sorts them based on chromosome passed.
        - Applies mask if mask else skips.
        - Calculates the jaccard index of each list's top and bottom 40 patients and prints it out.

    Returns:
    --------
        None (prints out jaccard index values)

    '''
    csv_file_name_1 = os.path.join(output_path, f'data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_path_1}')
    csv_file_name_2 = os.path.join(output_path, f'data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{csv_path_2}')

    df_1 = pd.read_csv(csv_file_name_1)
    df_2 = pd.read_csv(csv_file_name_2)

    df_1 = sort_df_by_column(df_1, column_to_sort_by, False)
    df_2 = sort_df_by_column(df_2, column_to_sort_by, False)

    if mask:
        df_1 = df_1[df_1["cycle"] == cycle]
        df_2 = df_2[df_2["cycle"] == cycle]

    df_1_top = list(df_1["patient_id"])[:n]
    df_2_top = list(df_2["patient_id"])[:n]

    df_1_bottom = list(df_1["patient_id"])[-n:]
    df_2_bottom = list(df_2["patient_id"])[-n:]

    jaccared_index_1_vs_2_top = jaccared_index(df_1_top, df_2_top)

    jaccared_index_1_vs_2_bottom = jaccared_index(df_1_bottom, df_2_bottom)

    venn_diagram(
        df_1_top, df_2_top, 
        csv_title_1, csv_title_2, 
        f"{csv_title_1} vs {csv_title_2} ({column_to_sort_by} - {cycle} - top {n})",  
        f"{csv_title_1}_vs_{csv_title_2}_{column_to_sort_by}_{cycle}_top_{n}", 
        metric_to_use,
        sequencing_type,
        output_path,
        jaccared_index_1_vs_2_top
    )

    venn_diagram(
        df_1_bottom, df_2_bottom, 
        csv_title_1, csv_title_2, 
        f"{csv_title_1} vs {csv_title_2} ({column_to_sort_by} - {cycle} - bottom {n})",  
        f"{csv_title_1}_vs_{csv_title_2}_{column_to_sort_by}_{cycle}_bottom_{n}", 
        metric_to_use,
        sequencing_type,
        output_path,
        jaccared_index_1_vs_2_bottom
    )

    print()
    print(f"{csv_title_1} top: {df_1_top}")
    print(f"{csv_title_2} top: {df_2_top}")
    print(f"top jaccard index: {jaccared_index_1_vs_2_top}")
    print()
    print(f"{csv_title_1} bottom: {df_1_bottom}")
    print(f"{csv_title_2} bottom: {df_2_bottom}")
    print(f"bottom jaccard index: {jaccared_index_1_vs_2_bottom}")
    print()

def jaccared_index(c1, c2):
    '''
    Parameters:
    -----------
        c1: List
            List of strings.

        c2: List
            List of strings.
    
    Function:
    ---------
        - Turns lists into sets (use for union and intersection not for duplicates).
        - Calculates union and intersection and returns jaccard index (intersection/union).

    Returns:
    --------
        None (prints out jaccard index values)

    '''
    s1, s2 = set(c1), set(c2)
    if not s1 and not s2:
        return None
    union = s1.union(s2)
    intersetction = s1.intersection(s2)
    return len(intersetction)/len(union)

def venn_diagram(c1, c2, c1_id, c2_id, title, file_name, metric_to_use, sequencing_type, output_path, jaccard_index=None):
    s1, s2 = set(c1), set(c2)
    v = venn2([s1, s2], (c1_id, c2_id))

    v.get_patch_by_id('10').set_color('#1f77b4')   # left-only
    v.get_patch_by_id('01').set_color('#ff7f0e')   # right-only
    v.get_patch_by_id('11').set_color('#2ca02c')   # intersection
    for subset in ('10', '01', '11'):
        v.get_patch_by_id(subset).set_alpha(0.6)

    # Style subset labels (numbers)
    for subset in ('10', '01', '11'):
        if v.get_label_by_id(subset):
            v.get_label_by_id(subset).set_fontsize(14)

    # Style set labels
    left_label, right_label = v.set_labels

    if left_label is not None:
        x, y = left_label.get_position()
        left_label.set_position((x - 0.2, y))
        left_label.set_fontsize(14)

    if right_label is not None:
        x, y = right_label.get_position()
        right_label.set_position((x + 0.1, y))
        right_label.set_fontsize(14)
    
    if jaccard_index:
        plt.text(
            0.5,
            0.52,
            f"Jaccard index: {jaccard_index:.3f}",
            ha="center",
            va="center"
        )
    
    plt.title(title, fontsize=12, fontweight='bold')

    plt.tight_layout()

    file_dir = os.path.join(output_path, f"entropy-analysis-{sequencing_type}-{metric_to_use}/csv_comparision_venn_diagrams")
    os.makedirs(file_dir, exist_ok=True)
    file_path = os.path.join(file_dir, f"{file_name}.png")
    plt.savefig(file_path)

    plt.close()

def cox_proportional_hazard_model_from_csvs(
        output_path, sub_dir_path, entropy_csv_path, pluvicto_csv_path, chr_to_pick, 
        cycle, event_col, time_to_event_column, other_pluvicto_columns,
        metric_to_use, sequencing_type, complex_sv_binary=None, high_entropy=None, tfx_col=None, tfx_cutoff=None):
    """
    Parameters:
    -----------
        output_path: String
            Output directory path for saving plots.
            
        sub_dir_path: String
            Subdirectory path within output_path for saving forest plots.

        entropy_csv_path: String
            Path to entropy per chromosome data CSV file.
        
        pluvicto_csv_path: String
            Path to pluvicto master sheet CSV file containing clinical metadata.

        chr_to_pick: String
            Column name of chromosome to filter entropy dataframe.
        
        cycle: String
            Treatment cycle to filter entropy data (e.g., 'C1', 'C2').

        event_col: String
            Column name for event indicator in pluvicto data.

        time_to_event_column: String
            Column name for time-to-event or censoring data in pluvicto data.

        other_pluvicto_columns: List
            List of additional column names from pluvicto to include as covariates in Cox model.

        metric_to_use: String
            Parameter to show which metric was used for entropy calculation.

        sequencing_type: String
            Type of sequencing data (e.g., 'deep', 'ulp_curated').

        complex_sv_binary: pandas DataFrame, optional
            Binary matrix of complex SVs (alternative to entropy for Cox model). Default is None.

        high_entropy: String, optional
            Tag to filter dataframe down to high entropy group or to keep the whole cohort. Default is None.

    Function:
    ---------
        - Reads entropy and pluvicto CSV files.
        - Filters entropy data to specified chromosome and cycle.
        - Adds responder groupings to pluvicto data.
        - Merges entropy/SV data with clinical data.
        - Fits Cox proportional hazards model with specified covariates.
        - Generates and saves forest plot (PDF and PNG).
        - Checks proportional hazards assumption.
        
    Returns:
    --------
        cox_data: pandas DataFrame
            Merged dataframe used for Cox model fitting.
        
        cox_model: lifelines.CoxPHFitter
            Fitted Cox proportional hazards model.
        
        chr_column: pandas DataFrame
            Dataframe containing patient_id and entropy/SV values used in model.
        
        summary: pandas DataFrame
            Summary table with hazard ratios, confidence intervals, and p-values for all covariates.
    """
    entropy_csv_file = os.path.join(output_path, entropy_csv_path)

    entropy_df = pd.read_csv(entropy_csv_file)

    pluvicto_df = pd.read_csv(pluvicto_csv_path, index_col=0)

    # Add tfx cutoff if not None
    pluvicto_df = pluvicto_df[pluvicto_df[tfx_col] >= tfx_cutoff]

    # Adds all responder groupings to pluvicto sheet
    try:
        pluvicto_groupings_df = add_responder_groupings(pluvicto_df)
    except:
        print("Groupings were not added\n")
        pluvicto_groupings_df = pluvicto_df

    # Filteres entropy df down to cycle then to specific chromosome column
    if cycle is not None:
        entropy_df = entropy_df[entropy_df["cycle"].str.contains(cycle, na=False)]
    else:
        entropy_df = entropy_df

    if cycle == "pre":
        entropy_df['patient_id'] = (
            entropy_df['patient_id'].fillna('NA').astype(str) + '_' +
            entropy_df['cycle'].fillna('NA').astype(str)
        )

    chr_column = entropy_df[["patient_id", chr_to_pick]]

    # if complex_sv_binary is not None use complex sv instead of entropy
    if complex_sv_binary is not None:
        chr_column = complex_sv_binary
        chr_to_pick = complex_sv_binary.columns[1]

    # Filteres pluvicto down to responder group and other columns
    pluvicto_filtered_df = pluvicto_groupings_df[[event_col, time_to_event_column, *other_pluvicto_columns]]

    # Merges dataframes
    if high_entropy:
        cox_data = pd.merge(pluvicto_filtered_df, chr_column, left_index=True, right_on="patient_id", how='right')
    else:
        cox_data = pd.merge(pluvicto_filtered_df, chr_column, left_index=True, right_on="patient_id")

    cox_data.dropna(inplace=True)
    cox_data = cox_data.reset_index(drop=True)

    # Build covariate list
    covariates = [chr_to_pick] + other_pluvicto_columns
    formula = ' + '.join(covariates)

    # z-score data
    cox_data[chr_to_pick] = (cox_data[chr_to_pick] - cox_data[chr_to_pick].mean()) / cox_data[chr_to_pick].std()

    cox_model = cox_proportional_hazards_model(cox_data, time_to_event_column, event_col, formula, chr_to_pick)

    file_dir = os.path.join(output_path, sub_dir_path)
    os.makedirs(file_dir, exist_ok=True)

    if complex_sv_binary is not None:
        title = f'Cox PH Model: Complex SV {chr_column.columns[1]} {f'- {high_entropy}' if high_entropy is not None else ''})'
        filename=f'forest_plot_complex_sv_{chr_column.columns[1]}{f'_{high_entropy}' if high_entropy is not None else ''}'
    else:
        # Create custom title based on parameters
        title = f'Cox PH Model: {chr_to_pick} (Cycle {cycle})'
        filename=f'forest_plot_{cycle}_{chr_to_pick}_entropy'

    out_path, summary = cox_forest_plot_pub(
        cox_model, 
        file_dir, 
        'p',
        filename=filename,
        title=title
    )
    
    # Optional: Check proportional hazards assumption
    print("\nChecking proportional hazards assumption...")
    try:
        cox_model.check_assumptions(cox_data, p_value_threshold=0.05, show_plots=False)
        print("Proportional hazards assumption check complete.")
    except Exception as e:
        print(f"Could not check assumptions: {e}")

    return cox_data, cox_model, chr_column, summary

def cox_proportional_hazard_model_from_df(
        output_path, sub_dir_path, df, chr_to_pick, 
        cycle, event_col, time_to_event_column, other_pluvicto_columns,
        metric_to_use, sequencing_type, complex_sv_binary=None, high_entropy=None):
    """
    Parameters:
    -----------
        output_path: String
            Output directory path for saving plots.
            
        sub_dir_path: String
            Subdirectory path within output_path for saving forest plots.

        df: pandas DataFrame
            

        chr_to_pick: String
            Column name of chromosome to filter entropy dataframe.
        
        cycle: String
            Treatment cycle to filter entropy data (e.g., 'C1', 'C2').

        event_col: String
            Column name for event indicator in pluvicto data.

        time_to_event_column: String
            Column name for time-to-event or censoring data in pluvicto data.

        other_pluvicto_columns: List
            List of additional column names from pluvicto to include as covariates in Cox model.

        metric_to_use: String
            Parameter to show which metric was used for entropy calculation.

        sequencing_type: String
            Type of sequencing data (e.g., 'deep', 'ulp_curated').

        complex_sv_binary: pandas DataFrame, optional
            Binary matrix of complex SVs (alternative to entropy for Cox model). Default is None.

        high_entropy: String, optional
            Tag to filter dataframe down to high entropy group or to keep the whole cohort. Default is None.

    Function:
    ---------
        - Reads entropy and pluvicto CSV files.
        - Filters entropy data to specified chromosome and cycle.
        - Adds responder groupings to pluvicto data.
        - Merges entropy/SV data with clinical data.
        - Fits Cox proportional hazards model with specified covariates.
        - Generates and saves forest plot (PDF and PNG).
        - Checks proportional hazards assumption.
        
    Returns:
    --------
        cox_data: pandas DataFrame
            Merged dataframe used for Cox model fitting.
        
        cox_model: lifelines.CoxPHFitter
            Fitted Cox proportional hazards model.
        
        chr_column: pandas DataFrame
            Dataframe containing patient_id and entropy/SV values used in model.
        
        summary: pandas DataFrame
            Summary table with hazard ratios, confidence intervals, and p-values for all covariates.
    """
    cox_data = df

    cox_data.dropna(inplace=True)
    cox_data = cox_data.reset_index(drop=True)

    # Build covariate list
    covariates = [chr_to_pick] + other_pluvicto_columns
    formula = ' + '.join(covariates)

    cox_model = cox_proportional_hazards_model(cox_data, time_to_event_column, event_col, formula, chr_to_pick)

    file_dir = os.path.join(output_path, sub_dir_path)
    os.makedirs(file_dir, exist_ok=True)

    if complex_sv_binary is not None:
        # title = f'Cox PH Model: Complex SV {chr_column.columns[1]} {f'- {high_entropy}' if high_entropy is not None else ''})'
        # filename=f'forest_plot_complex_sv_{chr_column.columns[1]}{f'_{high_entropy}' if high_entropy is not None else ''}'
        return
    else:
        # Create custom title based on parameters
        title = f'Cox PH Model: {chr_to_pick} (Cycle {cycle})'
        filename=f'forest_plot_{cycle}_{chr_to_pick}_entropy'

    out_path, summary = cox_forest_plot_pub(
        cox_model, 
        file_dir,
        'p',
        filename=filename,
        title=title
    )
    
    # Optional: Check proportional hazards assumption
    print("\nChecking proportional hazards assumption...")
    try:
        cox_model.check_assumptions(cox_data, p_value_threshold=0.05, show_plots=False)
        print("Proportional hazards assumption check complete.")
    except Exception as e:
        print(f"Could not check assumptions: {e}")

    return cox_data, cox_model, None, summary

def cox_proportional_hazards_model(cox_data, time_to_event_column, event_col, formula, chr_to_pick=None):
    if cox_data[chr_to_pick].sum() == 1:
        cox_model = CoxPHFitter(penalizer=0.1)
    else:
        cox_model = CoxPHFitter()

    cox_model.fit(cox_data, duration_col=time_to_event_column, event_col=event_col, formula=formula)
    return cox_model

def cox_forest_plot_pub(
    cox_model,
    output_dir,
    pval_col,
    filename="cox_forest_plot",
    title="Cox Proportional Hazards Model",
    x_lower=0.1,
    x_upper=None,
    fontsize=11
):

    # ------------------------------------------------------------------
    # Prepare summary
    # ------------------------------------------------------------------
    try:
        summary = cox_model.summary.reset_index()
        summary = summary.rename(columns={"index": "covariate"})
    except:
        summary = cox_model

    summary["HR"] = summary["exp(coef)"]
    summary["HR_lower"] = summary["exp(coef) lower 95%"]
    summary["HR_upper"] = summary["exp(coef) upper 95%"]
    summary[pval_col] = summary[pval_col]

    # Sort for nicer plotting (optional but common)
    summary = summary.sort_values(by='HR', ascending=False)

    n_vars = summary.shape[0]
    y_pos = np.arange(n_vars)

    # Axis limits
    if x_upper is None:
        x_upper = summary["HR_upper"].max() * 1.3

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 0.55 * n_vars + 1.5))

    # Confidence intervals
    ax.hlines(
        y=y_pos,
        xmin=summary["HR_lower"],
        xmax=summary["HR_upper"],
        linewidth=2,
        color="black"
    )

    # Point estimates
    ax.scatter(
        summary["HR"],
        y_pos,
        marker="s",
        s=15,
        color="black",
        zorder=3
    )

    # Reference line at HR = 1
    ax.axvline(1.0, linestyle="--", color="gray", linewidth=1)

    # Log scale on x-axis
    ax.set_xscale("log")

    # Major ticks at powers of 10 (10^-1, 10^0, 10^1, ...)
    ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=10))
    ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10.0))

    # Minor ticks without labels (clean)
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs="auto"))
    ax.xaxis.set_minor_formatter(NullFormatter())

    # Add padding to the left of 10^-1
    ax.set_xlim(min(x_lower * 0.8, 0.8), x_upper)

    # Y-axis labels
    ax.set_yticks(y_pos)
    try:
        ax.set_yticklabels(summary["covariate"], fontsize=fontsize)
    except:
        ax.set_yticklabels(summary.index, fontsize=fontsize)

    ax.invert_yaxis()

    # Labels & title
    ax.set_xlabel("Hazard Ratio (log scale)", fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize + 2, pad=10)

    # ------------------------------------------------------------------
    # Add HR (95% CI) + p-value text column
    # ------------------------------------------------------------------
    for i, row in summary.iterrows():
        hr_text = f'HR: {row["HR"]:.2f} [{row["HR_lower"]:.2f}, {row["HR_upper"]:.2f}]'
        p_text = f'p = {row[pval_col]:.2e}' if row[pval_col] < 0.001 else f'p = {row[pval_col]:.3f}'

        ax.text(
            1.02,   # x-position in AXES coordinates (just right of plot)
            y_pos[list(summary.index).index(i)],
            f"{hr_text}\n{p_text}",
            va="center",
            ha="left",
            fontsize=fontsize - 1,
            transform=ax.get_yaxis_transform(),  # <-- key line
            clip_on=False
        )

    # Make room on the right for text
    # ax.set_xlim(x_lower, x_upper)

    # Clean up spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{filename}.pdf")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    out_path = os.path.join(output_dir, f"{filename}.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    return out_path, summary

def add_responder_groupings(pluvicto_master_sheet):
    """
    Note: Adopted from '/prognostic-biomarker-model/scripts/data_processing_functions.py' add_responder_groupings function
    Parameters:
    -----------
        pluvicto_master_sheet: pandas DataFrame
            Data containing meta data.

    Function:
    ---------
        - Adds responder/non-responder label based on certain conditions listed in columns list.
        
    Returns:
    --------
        pluvicto_master_sheet: pandas DataFrame
            Data containing meta data and responder groupings.
    """
    columns = ["survival_days", "PSA_prog_days", "PSA_Progression", "tfx_prog_days", "T_cycles"]

    # Set Sample_ID as index
    try:
        pluvicto_master_sheet.set_index('Sample_ID', inplace=True)
    except Exception as e:
        print("\nSample_ID is already the index\nContinuing with adding responder groupings")
    
    for column in columns:
        # For binary indicators
        if column == "PSA_Progression":
            pluvicto_master_sheet[f"progression_group_{column}"] = pluvicto_master_sheet[column].apply(
                lambda x: "non-responder" if x == 1 else "responder"
            )
            continue
        
        if column == "T_cycles":
            column_temp = f"{column}_1_vs_6"
            pluvicto_master_sheet[f"progression_group_{column_temp}"] = pluvicto_master_sheet['T_cycles'].apply(
                lambda x: "responder" if x == 6 else ("non-responder" if x == 1 else np.nan)
            )
            
            column_temp = f"{column}_2_vs_6"
            pluvicto_master_sheet[f"progression_group_{column_temp}"] = pluvicto_master_sheet['T_cycles'].apply(
                lambda x: "responder" if x == 6 else ("non-responder" if x == 2 else np.nan)
            )
            
            column_temp = f"{column}_1_2_vs_5_6"
            pluvicto_master_sheet[f"progression_group_{column_temp}"] = pluvicto_master_sheet['T_cycles'].apply(
                lambda x: "responder" if x in [5, 6] else ("non-responder" if x in [1, 2] else np.nan)
            )
            
            column_temp = f"{column}_1-5_vs_6"
            pluvicto_master_sheet[f"progression_group_{column_temp}"] = pluvicto_master_sheet['T_cycles'].apply(
                lambda x: "responder" if x == 6 else ("non-responder" if x in [1, 2, 3, 4, 5] else np.nan)
            )
            continue
        
        if column == "survival_days":
            pluvicto_master_sheet[f"progression_group_{column}_252_cutoff"] = pluvicto_master_sheet[column].apply(
                lambda x: "non-responder" if x < 252 else "responder"
            )
        
        # Calculate quantiles on non-NA values
        quantile_df = pluvicto_master_sheet[pluvicto_master_sheet[column].notna()]
        quantile_vector = quantile_df[column].quantile([0.25, 0.5, 0.75])
        print(f"{column} {quantile_vector.values}")
        
        pluvicto_master_sheet[f"progression_group_{column}_median"] = pluvicto_master_sheet[column].apply(
            lambda x: "non-responder" if x < quantile_vector[0.5] else "responder"
        )
        
        pluvicto_master_sheet[f"progression_group_{column}_quartile"] = pluvicto_master_sheet[column].apply(
            lambda x: "non-responder" if x <= quantile_vector[0.25] else ("responder" if x >= quantile_vector[0.75] else np.nan)
        )
    
    return pluvicto_master_sheet

def kaplan_meier_plot(output_path, sub_dir_path, pluvicto_csv_path, event_col, time_to_event_column, stratify_col, chr_column, tfx_cycle, split_method, metric_to_use, sequencing_type, cox_model_hr_text=None, pairwise_comparisons=None, plot_title=None, split_params=None, group_labels=None, complex_sv=None, high_entropy=None, include_overall=False):
    """
    Parameters:
    -----------
        output_path : String
            Base directory for saving output files

        sub_dir_path : String
            Subdirectory path within output_path for saving Kaplan-Meier plots.

        pluvicto_csv_path : String
            Path to CSV file containing patient data and clinical metadata.

        event_col : String
            Column name for event indicator (1=event, 0=censored).

        time_to_event_column : String
            Column name for time-to-event or censoring duration.

        stratify_col : String
            Column name to use for stratification (e.g., entropy column).

        chr_column : pandas DataFrame
            DataFrame with patient_id and stratification feature columns.

        tfx_cycle : String
            Column name for tumor fraction cycle data to include as covariate.

        split_method : String
            Stratification method ('quartile_extremes', 'median', 'tertiles', 'quartiles', 'custom_quantiles', 'threshold', 'top_bottom_pct', 'tumor_fraction', 'tfx_and_entropy_median', 'tfx_and_entropy_custom', etc.).

        metric_to_use : String
            Parameter to show which metric was used for entropy calculation.

        sequencing_type : String
            Type of sequencing data (e.g., 'deep', 'ulp_curated').

        cox_model_hr_text : String, optional
            Pre-calculated Cox model hazard ratio text to display on plot. Default is None.

        pairwise_comparisons : List of tuples, optional
            List of group pairs for pairwise log-rank tests. Default is None.

        plot_title : String, optional
            Custom plot title. Default is None (auto-generated).

        split_params : dict, optional
            Additional parameters for splitting strategy (e.g., 'quantiles', 'value', 'percentage'). Default is None.

        group_labels : dict, optional
            Custom labels for groups mapping group numbers to descriptive labels. Overrides defaults. Default is None.

        complex_sv : String, optional
            Column name for complex structural variant binary matrix (alternative to entropy for stratification). Default is None.

        high_entropy : String, optional
            Tag to filter dataframe to high entropy group or keep whole cohort. Default is None.

    Function:
    ---------
        - Reads pluvicto CSV file and adds responder groupings.
        - Filters data for relevant event and time columns, handles missing values.
        - Merges patient data with entropy/SV data and tumor fraction data.
        - Applies stratification method to create patient groups.
        - Fits Kaplan-Meier curves for each group with confidence intervals.
        - Calculates survival statistics (log-rank test, hazard ratios).
        - Generates Kaplan-Meier survival plot with group labels, statistics, and legend.
        - Saves plot as PNG and PDF files to specified output directory.
        
    Returns:
    --------
        None (saves Kaplan-Meier curves to files).
    """

    # Read data
    pluvicto_df = pluvicto_csv_path
    try:
        pluvicto_groupings_df = add_responder_groupings(pluvicto_df)
    except:
        pluvicto_groupings_df = pluvicto_df

    # Get relevant columns and drop missing values
    if "ctdnaq" in split_method:
        pluvicto_filtered_df = pluvicto_groupings_df[[event_col, time_to_event_column, "ctdnaq_category"]].dropna()
    else:
        pluvicto_filtered_df = pluvicto_groupings_df[[event_col, time_to_event_column]].dropna()

    try:
        pluvicto_tfx_col = pluvicto_groupings_df[[tfx_cycle]].dropna()
    except:
        pluvicto_tfx_col = None

    # Merge dataframes
    km_data_temp = pd.merge(pluvicto_filtered_df, chr_column, left_index=True, right_on="patient_id")

    if pluvicto_tfx_col is not None:
        km_data = pd.merge(km_data_temp, pluvicto_tfx_col, left_on="patient_id", right_index=True)
    else:
        km_data = km_data_temp

    km_data_stratified, default_labels, n_groups = stratify_data(
        km_data, 
        stratify_col, 
        split_method, 
        split_params
    )

    top = km_data[stratify_col].max()
    median = km_data[stratify_col].median()
    bottom = km_data[stratify_col].min()

    # Update km_data_stratified so the sv becomes the grouping, default_labels reflects sv presence, and stratify_col is the sv column
    if complex_sv:
        km_data_stratified = km_data.rename(columns={complex_sv: 'group'})
        default_labels = {0: 'Not present', 1: 'Present'}

    if group_labels is not None:
        final_labels = group_labels
    else:
        final_labels = default_labels

    # Setup plot
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Fit KM curves for each responder group
    kmf = KaplanMeierFitter()
    groups = sorted(km_data_stratified['group'].unique())  # Sort to ensure consistent order
    color_palette = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33']

    if include_overall:
        groups = ['overall'] + groups
    
    for i, group in enumerate(groups):

        if group == 'overall':
            group_data = km_data_stratified
            label_name = "All Patients"
        else:
            group_data = km_data_stratified[km_data_stratified['group'] == group]
            label_name = final_labels[group]

        n_events = group_data[event_col].sum()
        n_censored = len(group_data) - n_events

        label = f'{label_name}: n={len(group_data)}, events={int(n_events)}, censored={int(n_censored)}'

        kmf.fit(group_data[time_to_event_column], group_data[event_col], label=label)
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color_palette[i % len(color_palette)], linewidth=2.5)

        ax.scatter(
            kmf.event_table.index[kmf.event_table['censored'] > 0],
            kmf.survival_function_.loc[kmf.event_table.index[kmf.event_table['censored'] > 0]].values,
            marker='|',
            s=100,
            color=color_palette[i % len(color_palette)],
            zorder=5,
            linewidths=2
        )
    
    if complex_sv:
        covariates = ['group'] + [tfx_cycle]
        hr_formula = ' + '.join(covariates)
    else:
        if tfx_cycle is not None:
            covariates = [stratify_col] + [tfx_cycle]
        else:
            covariates = [stratify_col]
        hr_formula = ' + '.join(covariates)

    # Calculate statistics 
    stats_text = calculate_survival_statistics(
        km_data_stratified,
        time_to_event_column,
        event_col,
        stratify_col,
        hr_formula,
        n_groups,
        cox_model_hr_text,
        pairwise_comparisons
    )

    if cox_model_hr_text:
        stats_text = stats_text + cox_model_hr_text

    y_pos = 0.90 - (n_groups * 0.05)
     
    if stats_text:
        ax.text(0.98, y_pos, stats_text, transform=ax.transAxes, fontsize=11,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                verticalalignment='bottom', horizontalalignment='right')
    
    # Format plot
    ax.set_xlabel('Time (days)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Overall Survival Probability', fontsize=12, fontweight='bold')

    if plot_title is None:
        plot_title = f'Kaplan-Meier Survival Curves based on {stratify_col} Entropy ({split_method.replace("_", " ").title()}'

    if split_method == "median":
        plot_title = f'{plot_title}: {median:3f} [{bottom:3f}, {top:3f}])'
    else:
        plot_title = f'{plot_title})'

    ax.set_title(plot_title, fontsize=13, fontweight='bold')

    ax.legend(loc='best', frameon=True, fontsize=11)
    ax.grid(alpha=0.3, linestyle=':')
    ax.set_ylim([0, 1.05])
    plt.tight_layout()
    
    # Save
    file_dir = os.path.join(output_path, sub_dir_path)
    os.makedirs(file_dir, exist_ok=True)

    plt.savefig(os.path.join(file_dir, f'km_curve_{stratify_col}_{split_method}{f'_{high_entropy}' if high_entropy is not None else ''}.png'), dpi=600, bbox_inches='tight')
    plt.savefig(os.path.join(file_dir, f'km_curve_{stratify_col}_{split_method}{f'_{high_entropy}' if high_entropy is not None else ''}.pdf'), bbox_inches='tight')
    plt.close()
    
    print(f"KM plot saved to: {file_dir}")

def stratify_data(km_data, stratify_col, split_method, split_params):
    """
    Parameters:
    -----------
        km_data: pandas DataFrame
            DataFrame with at least one column in it (stratify_col).

        stratify_col: String
            Name of column to split data on.

        split_method: String
            Method to split stratify_col by.

        split_params: Not exactly sure because I never use it.

    Function:
    ---------
        - Creates a column called 'group' that binarizes the stratify_col depending on the split_method.
        - Returns the df with the 'group' column, along with labels, and the number of groups. s
    
    Returns:
    --------
        km_data_stratified: pandas DataFrame
            DataFrame with 'group' column added

        group_labels: dictionary
            dict mapping group numbers to descriptive labels

        n_groups: int
            int, number of groups created
    """
    km_data_copy = km_data.copy()
    
    if split_method == 'quartile_extremes':
        # Original method: Q1 vs Q4
        q1 = km_data_copy[stratify_col].quantile(0.25)
        q3 = km_data_copy[stratify_col].quantile(0.75)
        km_data_copy = km_data_copy[(km_data_copy[stratify_col] <= q1) | (km_data_copy[stratify_col] >= q3)]
        km_data_copy['group'] = (km_data_copy[stratify_col] >= q3).astype(int)
        group_labels = {0: 'Low (Q1)', 1: 'High (Q4)'}
        n_groups = 2
        
    elif split_method == 'median':
        # Split at median
        median = km_data_copy[stratify_col].median()
        km_data_copy['group'] = (km_data_copy[stratify_col] > median).astype(int)
        group_labels = {0: 'Low (< median)', 1: 'High (≥ median)'}
        n_groups = 2
        
    elif split_method == 'tertiles':
        # Split into three equal groups
        km_data_copy['group'] = pd.qcut(km_data_copy[stratify_col], q=3, labels=[0, 1, 2], duplicates='drop')
        km_data_copy['group'] = km_data_copy['group'].astype(int)
        group_labels = {0: 'Low (T1)', 1: 'Medium (T2)', 2: 'High (T3)'}
        n_groups = 3
        
    elif split_method == 'quartiles':
        # Split into four equal groups
        categories = pd.qcut(km_data_copy[stratify_col], q=4, duplicates='drop')
        km_data_copy['group'] = pd.factorize(categories)[0]
        n_groups = km_data_copy['group'].nunique()
        group_labels = {i: f'Q{i+1}' for i in range(n_groups)}
        
    elif split_method == 'custom_quantiles':
        # Use custom quantile cutoffs
        quantiles = split_params.get('quantiles', [0.33, 0.67])
        km_data_copy['group'] = pd.cut(
            km_data_copy[stratify_col], 
            bins=[-np.inf] + [km_data_copy[stratify_col].quantile(q) for q in quantiles] + [np.inf],
            labels=range(len(quantiles) + 1)
        )
        km_data_copy['group'] = km_data_copy['group'].astype(int)
        n_groups = len(quantiles) + 1
        group_labels = {i: f'Group {i+1}' for i in range(n_groups)}
        
    elif split_method == 'threshold':
        # Split at custom threshold
        threshold = split_params.get('value', km_data_copy[stratify_col].median())
        km_data_copy['group'] = (km_data_copy[stratify_col] >= threshold).astype(int)
        group_labels = {0: f'Low (< {threshold:.3f})', 1: f'High (≥ {threshold:.3f})'}
        n_groups = 2
        
    elif split_method == 'top_bottom_pct':
        # Compare top X% vs bottom X%
        pct = split_params.get('percentage', 0.25)
        lower_cutoff = km_data_copy[stratify_col].quantile(pct)
        upper_cutoff = km_data_copy[stratify_col].quantile(1 - pct)
        km_data_copy = km_data_copy[
            (km_data_copy[stratify_col] <= lower_cutoff) | 
            (km_data_copy[stratify_col] >= upper_cutoff)
        ]
        km_data_copy['group'] = (km_data_copy[stratify_col] >= upper_cutoff).astype(int)
        pct_label = int(pct * 100)
        group_labels = {0: f'Bottom {pct_label}%', 1: f'Top {pct_label}%'}
        n_groups = 2

    elif split_method == "tumor_fraction":
        median = km_data_copy["TFx_C1"].median()
        km_data_copy['group'] = (km_data_copy["TFx_C1"] >= median).astype(int)
        group_labels = {0: 'Low (< median)', 1: 'High (≥ median)'}
        n_groups = 2

    elif split_method == "tfx_and_entropy_median":
        tfx_median = km_data_copy["TFx_C1"].median()
        entropy_median = km_data_copy[stratify_col].median()
        conditions = [
            (km_data_copy["TFx_C1"] >= tfx_median) & (km_data_copy[stratify_col] >= entropy_median), 
            (km_data_copy["TFx_C1"] >= tfx_median) & (km_data_copy[stratify_col] < entropy_median), 
            (km_data_copy["TFx_C1"] <= tfx_median) & (km_data_copy[stratify_col] >= entropy_median), 
            (km_data_copy["TFx_C1"] <= tfx_median) & (km_data_copy[stratify_col] < entropy_median), 
        ]
        condition_labels = [0, 1, 2, 3]
        km_data_copy['group'] = np.select(conditions, condition_labels, default=-1)
        group_labels = {0: 'high TFx, high entropy', 1: 'high TFx, low entropy', 2: 'low TFx, high entropy', 3: 'low TFx, low entropy'}
        n_groups = 4

    elif "tfx_and_entropy_custom" in split_method :
        # TFx tag on the end of split_method
        tag = split_method[split_method.rfind('_')+1:]

        if tag == 'low':
            km_data_copy = km_data_copy[(km_data_copy['TFx_C1'] > 0.03) & (km_data_copy['TFx_C1'] <= 0.15)]
            entropy_median = km_data_copy[stratify_col].median()
            conditions = [
                (km_data_copy[stratify_col] < entropy_median),
                (km_data_copy[stratify_col] >= entropy_median)
            ]

        elif tag == 'medium':
            km_data_copy = km_data_copy[(km_data_copy['TFx_C1'] > 0.15) & (km_data_copy['TFx_C1'] < 0.45)]
            entropy_median = km_data_copy[stratify_col].median()
            conditions = [
                (km_data_copy[stratify_col] < entropy_median),
                (km_data_copy[stratify_col] >= entropy_median)
            ]

        elif tag == 'high':
            km_data_copy = km_data_copy[(km_data_copy['TFx_C1'] >= 0.45)]
            entropy_median = km_data_copy[stratify_col].median()
            conditions = [
                (km_data_copy[stratify_col] < entropy_median),
                (km_data_copy[stratify_col] >= entropy_median)
            ]
        
        condition_labels = [0, 1]
        km_data_copy['group'] = np.select(conditions, condition_labels, default=-1)
        group_labels = {0: f'{tag} TFx, low entropy', 1: f'{tag} TFx, high entropy'}
        n_groups = 3

    elif "ctdnaq" in split_method:
        if "entropy" in split_method:
            if split_method == "favorable_ctdnaq_median_entropy":
                km_data_copy = km_data_copy[(km_data_copy["ctdnaq_category"] == "Undetectable") | (km_data_copy["ctdnaq_category"] == "Low")]
            elif split_method == "unfavorable_ctdnaq_median_entropy":
                km_data_copy = km_data_copy[(km_data_copy["ctdnaq_category"] == "Moderate") | (km_data_copy["ctdnaq_category"] == "High")]
        
            entropy_median = km_data_copy[stratify_col].median()
            conditions = [
                km_data_copy[stratify_col] < entropy_median, 
                km_data_copy[stratify_col] >= entropy_median, 
            ]
            # Drop ctdnaq_category as it will mess up downstream analysis
            km_data_copy = km_data_copy.drop("ctdnaq_category", axis=1)
            condition_labels = [0, 1]
            km_data_copy['group'] = np.select(conditions, condition_labels, default=-1)
            group_labels = {0: 'Low entropy', 1: 'High entropy'}
            n_groups = 2

        elif '_grt_than_3pct_TFx' in split_method:
            conditions = [
                (km_data_copy["ctdnaq_category"] == "Low"), 
                (km_data_copy["ctdnaq_category"] == "Moderate") | (km_data_copy["ctdnaq_category"] == "High"), 
            ]

            condition_labels = [0, 1]
            km_data_copy['group'] = np.select(conditions, condition_labels, default=-1)
            group_labels = {0: 'Favorable', 1: 'Unfavorable'}
            n_groups = 2

            # Drop undetectable group
            drop_cols = km_data_copy['group'] == -1
            km_data_copy = km_data_copy[~drop_cols]
            
        else:
            conditions = [
                (km_data_copy["ctdnaq_category"] == "Undetectable") | 
                (km_data_copy["ctdnaq_category"] == "Low"), 
                (km_data_copy["ctdnaq_category"] == "Moderate") | (km_data_copy["ctdnaq_category"] == "High"), 
            ]

            condition_labels = [0, 1]
            km_data_copy['group'] = np.select(conditions, condition_labels, default=-1)
            group_labels = {0: 'Favorable', 1: 'Unfavorable'}
            n_groups = 2
        
    else:
        raise ValueError(f"Unknown split_method: {split_method}")
    
    return km_data_copy, group_labels, n_groups

def calculate_survival_statistics(km_data, time_col, event_col, stratify_col, hr_formula, n_groups, cox_model_hr_text=None, pairwise_comparisons=None):
    '''
    Parameters:
    -----------
        km_data: pandas DataFrame
            Dataframe with the kaplin Meier data.
        
        time_col: String
            Name of column with time data.

        event_col: String
            Name of column with event data.

        stratify_col: String
            Name of column that data is being stratified on.
        
        hr_formula: String
            List of covariates (Seperated by '+' if multiple).

        n_groups: int
            Number of groups to split data into.

        cox_model_hr_text: String (Default: None)
            The HR and p-value of the cox model run previously. 

        pairwise_comparisons: List of tuples (Default: None)
            Groups that will be compared together from 'group' column. 

    Function:
    ---------
        - If data is split into 2 groups then do a log-rank test on two groups of data, finds HR value, formats string, and returns it.
        - If data is split into more than 2 groups it does a multivariate log-rank test, formats string, and returns it.
        - If data is split into more than 2 groups and pairwise_comparisons is not None, it does log-rank tests on each tuple group, formats string, and returns it.
        
    Returns:
    --------
        stats_text: String
            String for km curve plot (contains Log-rank and HR value).
    '''
    stats_text = ""
    
    # For 2 groups: pairwise log-rank and HR
    if n_groups == 2:
        g0 = km_data[km_data['group'] == 0]
        g1 = km_data[km_data['group'] == 1]
        
        lr_result = logrank_test(g0[time_col], g1[time_col], g0[event_col], g1[event_col])
        
        if lr_result.p_value < 0.00001:
            p_value = "< 0.00001"
        else:
            p_value = f"{lr_result.p_value:.5f}"

        stats_text = f'Log-rank p = {p_value}\n'

        if not cox_model_hr_text:
            cox_model = cox_proportional_hazards_model(km_data, time_col, event_col, hr_formula, 'group')
            row = cox_model.summary.loc[stratify_col]
            hr_ci_text = f'HR: {row["exp(coef)"]:.2f} (CI 95%: {row["exp(coef) lower 95%"]:.2f}–{row["exp(coef) upper 95%"]:.2f})'
                        
            stats_text = stats_text + hr_ci_text

    elif n_groups > 2 and pairwise_comparisons is not None:
        comparison_results = []

        for group1, group2 in pairwise_comparisons:
            if group1 == "overall":
                g1 = km_data
            else:
                g1 = km_data[km_data['group'] == group1]
            
            if group2 == "overall":
                g2 = km_data
            else:
                g2 = km_data[km_data['group'] == group2]

            lr_result = logrank_test(g1[time_col], g2[time_col], g1[event_col], g2[event_col])

            if lr_result.p_value < 0.001:
                p_value = "< 0.001"
            else:
                p_value = f"{lr_result.p_value:.4f}"

            comparison_results.append(f'Groups {group1} vs {group2}: p={p_value}')

            stats_text = stats_text + f'{group1} vs {group2}: p = {p_value}\n'

        stats_text = stats_text.rstrip('\n')
    
    # For >2 groups: multivariate log-rank test
    elif n_groups > 2:
        lr_result = multivariate_logrank_test(
            km_data[time_col], 
            km_data['group'], 
            km_data[event_col]
        )
        
        if lr_result.p_value < 0.001:
            p_value = "< 0.001"
        else:
            p_value = f"{lr_result.p_value:.4f}"
        
        stats_text = f'Log-rank p = {p_value}'

    return stats_text


def complex_sv_and_entropy_visualization(df, output_path, sv_column, entropy_column, metric_to_use, sequencing_type, caller=None):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            Dataframe with a column of structural variant (sv) calls and a column of chr8 entropy values from ulp data.
        
        output_path: String
            Path output directory.

        sv_column: String
            Name of sv column. 
        
        entropy_column: String
            Name of entropy column.

        metric_to_use: String
            Parameter to show which metric was used from cna.seg file to use for entropy.

        sequencing_type: String
            Parameter to show which type of sequencing was used.

        caller: String
            Default: None
            Name of sv caller to filter sv column down to. 

    Function:
    ---------
        - If caller is passed - df is filtered to that caller.
        - Calls boxplot_or_swarmplot twice, once for both types of plots. 
        
    Returns:
    --------
        None
    '''
    if caller:
        df = df[df["caller"] == caller]

    boxplot_or_swarmplot(df, output_path, sv_column, entropy_column, "boxplot", metric_to_use, sequencing_type, caller)
    boxplot_or_swarmplot(df, output_path, sv_column, entropy_column, "swarmplot", metric_to_use, sequencing_type, caller)

def boxplot_or_swarmplot(df, output_path, sv_column, entropy_column, plot_type, metric_to_use, sequencing_type, caller):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            Dataframe with a column of structural variant (sv) calls and a column of chr8 entropy values from ulp data.
        
        output_path: String
            Path output directory.

        sv_column: String
            Name of sv column. 
        
        entropy_column: String
            Name of entropy column.

        plot_type: String
            Type of plot to create: swarm or box plot. 

        caller: String
            Default: None
            Name of sv caller to filter sv column down to. 

        metric_to_use: String
            Parameter to show which metric was used from cna.seg file to use for entropy.

        sequencing_type: String
            Parameter to show which type of sequencing was used.

    Function:
    ---------
        - If caller is passed - create tag on title and file name signifying using 'all' callers.
        - Creates a plot based on plot_type and saves plot to correct output path. 
        
    Returns:
    --------
        None
    '''
    # If caller is None then use 'all' in title and file path
    if caller is None:
        caller = "all"
    
    plt.figure(figsize=(15, 10))
    if plot_type == "boxplot":
        sns.boxplot(data=df, x=sv_column, y=entropy_column)

        # # Add significance starts and pvalues
        # ax = sns.boxplot(data=df, x=sv_column, y=entropy_column)
        # pairs = list(combinations(df[sv_column].unique(), 2))
        # annotator = Annotator(ax, pairs, data=df, x=sv_column, y=entropy_column)
        # annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
        # annotator.apply_and_annotate()
    else:
        sns.swarmplot(data=df, x=sv_column, y=entropy_column, size=3)

    plt.title(f"{caller} Complex SV Events vs Entropy Values (chr8)")
    plt.xlabel("Complex SV Events")
    plt.xticks(rotation=45)
    plt.ylabel("Entropy Values")

    plot_sub_dir = os.path.join(output_path, "complex_svs_vs_entropy_boxplots")
    os.makedirs(plot_sub_dir, exist_ok=True)
    boxplot_file = os.path.join(plot_sub_dir, f"chr8_{caller}_complex_sv_events_vs_entropy_{plot_type}.png")
    plt.tight_layout()
    plt.savefig(boxplot_file)
    plt.close()

def tukey_hsd(df, output_path, sv_column, entropy_column, metric_to_use, sequencing_type):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            Dataframe with a column of structural variant (sv) calls and a column of chr8 entropy values from ulp data.
        
        output_path: String
            Path output directory.

        sv_column: String
            Name of sv column. 
        
        entropy_column: String
            Name of entropy column.

        metric_to_use: String
            Parameter to show which metric was used from cna.seg file to use for entropy.

        sequencing_type: String
            Parameter to show which type of sequencing was used.

    Function:
    ---------
        - Turns entropy column to numeric.
        - Does a pairwise tukey HSD test on entropy and sv columns (does family error correction).
        - Then creates a boxplot for all significant tests.
        
    Returns:
    --------
        None
    '''
    df[entropy_column] = pd.to_numeric(df[entropy_column], errors='coerce')

    tukey = pairwise_tukeyhsd(endog=df[entropy_column], 
                          groups=df[sv_column], 
                          alpha=0.05)
    
    tukey_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
    tukey_df_filtered = tukey_df[tukey_df["reject"] == True]

    for i, row in tukey_df_filtered.iterrows():
        df_masked = df[ df[sv_column].isin([row["group1"], row["group2"]])]
        sns.boxplot(data=df_masked, x=sv_column, y=entropy_column)

        plt.title("Complex SV Events vs Entropy Values (chr8)")
        plt.xlabel("Complex SV Events")
        plt.xticks(rotation=45)
        plt.ylabel("Entropy Values")

        if row['p-adj'] < 0.001:
            p_value = "< 0.001"
        else:
            p_value = row['p-adj']

        plt.text(
            0.95, 0.95, f"p-adj: {p_value}", 
            transform=plt.gca().transAxes,
            ha='right', va='top',
            fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )

        boxplot_dir = os.path.join(output_path, f"entropy-analysis-{sequencing_type}-{metric_to_use}")
        boxplot_sub_dir = os.path.join(boxplot_dir, "complex_svs_vs_entropy_boxplots")
        significant_dir = os.path.join(boxplot_sub_dir, "significant_svs")
        os.makedirs(significant_dir, exist_ok=True)
        boxplot_file = os.path.join(significant_dir, f"{row["group1"]}_vs_{row["group2"]}_boxplot.png")
        plt.tight_layout()
        plt.savefig(boxplot_file)
        plt.close()

def complex_sv_and_entropy_stacked_histogram(df, output_path, sv_column, entropy_column, metric_to_use, sequencing_type):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            Dataframe with a column of structural variant (sv) calls and a column of chr8 entropy values from ulp data.
        
        output_path: String
            Path output directory.

        sv_column: String
            Name of sv column. 
        
        entropy_column: String
            Name of entropy column.

        metric_to_use: String
            Parameter to show which metric was used from cna.seg file to use for entropy.

        sequencing_type: String
            Parameter to show which type of sequencing was used.

    Function:
    ---------
        - Builds a dictionary that makes a list of entropy values for each sv. 
        - Find the top entropy value and builds a bin range based on it.
        - Builds lists of data and labels as well as a colormap.
        - Creates a temp plot to standardize x- and y-axies across all histograms.
        - Creates a histogram of entropy values for each svs.
        - Creates a stacked histogram of all svs. 
        
    Returns:
    --------
        None
    '''
    # Build sv dictionary
    sv_dict = {}
    for i, row in df.iterrows():
        sv_dict.setdefault(row[sv_column], []).append(row[entropy_column])

    top_value = df[entropy_column].max()

    # Create bin range (for top_value < 10 it will increment by .5 instead of 1)
    if top_value < 5:
        bin_range = [x * .25 for x in range(0, math.ceil(top_value)*4)]
    elif top_value < 10:
        bin_range = [x * .5 for x in range(0, math.ceil(top_value)*2)]
    else:
        bin_range = list(range(0, math.ceil(top_value)))

    # Build list of data and labels
    list_of_data = []
    labels = []
    for i, item in sv_dict.items():
        list_of_data.append(item)
        labels.append(i)

    colorblind_colors = [
        '#332288',  # Indigo
        '#117733',  # Green
        '#44AA99',  # Teal
        '#88CCEE',  # Cyan
        '#DDCC77',  # Sand
        '#CC6677',  # Rose
        '#AA4499',  # Purple
        '#882255',  # Wine
        '#661100',  # Brown
        '#999933'   # Olive
    ]
    
    # Create color map for each SV type
    color_map = {label: colorblind_colors[i % len(colorblind_colors)] 
                 for i, label in enumerate(labels)}

    # Get histogram counts to determine y-axis max from Temporary histogram
    counts, _, _ = plt.hist(list_of_data, bins=bin_range, stacked=True, color=[color_map[label] for label in labels])
    max_y = counts.max()
    plt.close()  # Close the temporary figure

    for i, item in sv_dict.items():
        create_histogram(
            output_path, 
            item, 
            i,
            max_y,
            color_map[i],
            'Chr8 Entropy', 
            'Frequency', 
            f'Stacked Histogram of Chr8 Complex SVs {i} vs Entropy',
            metric_to_use, 
            sequencing_type,
            f"{i}_complex_svs_vs_entropy_stacked_histogram", 
            bin_range
        )

    create_histogram(
        output_path, 
        list_of_data, 
        labels,
        max_y,
        [color_map[label] for label in labels],
        'Chr8 Entropy', 
        'Frequency', 
        'Stacked Histogram of Chr8 Complex SVs vs Entropy', 
        metric_to_use, 
        sequencing_type,
        "all_complex_svs_vs_entropy_stacked_histogram", 
        bin_range
    )

def complex_sv_and_entropy_stacked_histogram_per_sv(output_path, sub_dir, df, entropy_column, title, file_name, ylabel='Number of Patients', xlabel='Structural Variant'):
    """
    Publication-ready stacked bar plot of SV frequency by high/low entropy.

    Parameters:
    -----------
    output_path: str
        Path to save the figure.
    sub_dir: str
        Subdirectory for saving the figure.
    df: pd.DataFrame
        DataFrame with SV binary columns and entropy column.
    entropy_column: str
        Column name for entropy values.
    title: str
        Plot title.
    file_name: str
        Name of the output file (without extension).
    ylabel: str
        Y-axis label.
    xlabel: str
        X-axis label.
    """

    # ---------------------------------
    # Define high vs low entropy
    # ---------------------------------
    median_entropy = df[entropy_column].median()
    df['entropy_group'] = ['High' if x > median_entropy else 'Low' for x in df[entropy_column]]

    # ---------------------------------
    # Identify SV columns (binary presence/absence)
    # ---------------------------------

    sv_cols = [col for col in df.columns if col not in ['patient_id', entropy_column, "entropy_group"]]

    # ---------------------------------
    # Count number of patients with each SV in high vs low entropy
    # ---------------------------------

    stacked_counts = pd.DataFrame()
    for sv in sv_cols:
        counts = df.groupby('entropy_group')[sv].sum()
        stacked_counts[sv] = counts

    # Make sure rows are ordered as Low / High
    stacked_counts = stacked_counts.reindex(['Low', 'High'])

    # ---------------------------------
    # Order bars by total SV counts
    # ---------------------------------

    total_counts = stacked_counts.sum(axis=0)
    stacked_counts = stacked_counts[total_counts.sort_values(ascending=False).index]

    # ---------------------------------
    # Plot stacked bar
    # ---------------------------------

    colors = ['#1f77b4', '#d62728']  # Low = blue, High = red
    ax = stacked_counts.T.plot(kind='bar', stacked=True, color=colors, figsize=(10,6), edgecolor='black')

    # ---------------------------------
    # Style
    # ---------------------------------

    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(title='Entropy Group', fontsize=10, title_fontsize=11, loc='upper right')

    plt.tight_layout()

    # ---------------------------------
    # Save figure
    # ---------------------------------

    plot_dir = os.path.join(output_path, sub_dir)
    os.makedirs(plot_dir, exist_ok=True)
    plot_file = os.path.join(plot_dir, file_name)
    plt.savefig(f"{plot_file}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{plot_file}.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def horizontal_bar_plot(output_path, sub_dir, df, title, file_name, ylabel='Structural Variant', xlabel='Number of Patients'):
    """
    Horizontal bar plot of SV frequencies.

    Parameters:
    -----------
    output_path : str
        Path to save the figure.
    sub_dir : str
        Subdirectory for saving the figure.
    df : pd.DataFrame
        DataFrame with SV binary columns.
    title : str
        Plot title.
    file_name : str
        Name of output file (without extension).
    ylabel : str
        Y-axis label.
    xlabel : str
        X-axis label.
    """

    # ---------------------------------
    # Identify SV columns
    # ---------------------------------

    sv_cols = [col for col in df.columns if col != 'patient_id']

    # ---------------------------------
    # Count SV frequency
    # ---------------------------------

    sv_counts = df[sv_cols].sum(axis=0)

    # Sort ascending for clean horizontal display
    sv_counts = sv_counts.sort_values(ascending=True)

    # ---------------------------------
    # Create figure
    # ---------------------------------

    fig_height = max(4, len(sv_counts) * 0.45)

    fig, ax = plt.subplots(
        figsize=(8, fig_height),
        facecolor='white'
    )

    # ---------------------------------
    # Bar plot
    # ---------------------------------

    bar_color = '#3E5871'  # dark blue-grey

    bars = ax.barh(
        y=np.arange(len(sv_counts)),
        width=sv_counts.values,
        height=0.72,
        color=bar_color,
        edgecolor='none',
        zorder=3
    )

    # Add y tick labels
    ax.set_yticks(np.arange(len(sv_counts)))
    ax.set_yticklabels(sv_counts.index)

    # ---------------------------------
    # Add count labels
    # ---------------------------------

    x_offset = max(sv_counts.values) * 0.015

    for i, value in enumerate(sv_counts.values):
        ax.text(
            value * 0.5,
            i,
            str(int(value)),
            va='center',
            ha='left',
            fontsize=6,
            color='white'
        )

    # ---------------------------------
    # Axes labels and title
    # ---------------------------------

    ax.set_xlabel(xlabel, fontsize=12, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=12, labelpad=10)

    ax.set_title(
        title,
        fontsize=14,
        fontweight='bold',
        pad=15
    )

    # ---------------------------------
    # Clean Nature-style formatting
    # ---------------------------------

    # Remove unnecessary spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Light left/bottom spines
    ax.spines['left'].set_color('#B0B0B0')
    ax.spines['bottom'].set_color('#B0B0B0')

    # Make bars not touch y-axis
    ax.set_xlim(left=-max(sv_counts.values) * 0.03)

    # Vertical dashed gridlines
    ax.xaxis.grid(
        True,
        linestyle='--',
        linewidth=0.8,
        color='#D3D3D3',
        alpha=0.8,
        zorder=0
    )

    ax.yaxis.grid(False)

    # Tick styling
    ax.tick_params(
        axis='x',
        labelsize=10,
        colors='#333333'
    )

    ax.tick_params(
        axis='y',
        labelsize=10,
        colors='#333333',
        length=0
    )

    # ---------------------------------
    # Layout
    # ---------------------------------

    plt.tight_layout()

    # ---------------------------------
    # Save figure
    # ---------------------------------

    plot_dir = os.path.join(output_path, sub_dir)
    os.makedirs(plot_dir, exist_ok=True)
    plot_file = os.path.join(plot_dir, file_name)
    plt.savefig(f"{plot_file}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{plot_file}.pdf", dpi=300, bbox_inches='tight')
    plt.close()

def proportional_stacked_histogram_per_sv(output_path, df, id_column, entropy_column, title, file_name, ylabel='Patients (%)', xlabel='Structural Variant', panel_label=None, correction_method="fdr_bh", min_expected=5):
    """
    Publication-ready stacked bar plot of SV frequency by high/low entropy.

    Parameters:
    -----------
    output_path: str
        Path to save the figure.
    sub_dir: str
        Subdirectory for saving the figure.
    df: pd.DataFrame
        DataFrame with SV binary columns and entropy column.
    id_column:
        Column name for ids.
    entropy_column: str
        Column name for entropy values.
    title: str
        Plot title.
    file_name: str
        Name of the output file (without extension).
    ylabel: str
        Y-axis label.
    xlabel: str
        X-axis label.
    panel_label : str, optional
        Single letter/number placed at upper-left (e.g. "a", "b").
    correction_method : str
        Multiple-testing correction passed to
        statsmodels.stats.multitest.multipletests.
        Common options: "fdr_bh" (Benjamini-Hochberg), "bonferroni", "holm".
    min_expected : float
        Minimum expected cell count below which Fisher's exact test is used
        instead of chi-squared.
    """
    # ------------------------------------------------------------------
    # 1. Prepare data
    # ------------------------------------------------------------------
    plot_df = df.copy()
 
    median_val = plot_df[entropy_column].median()
    plot_df["_group"] = np.where(
        plot_df[entropy_column] >= median_val, "high", "low"
    )
 
    excluded = {id_column, entropy_column, "_group"}
    binary_cols = [
        col
        for col in plot_df.columns
        if col not in excluded
        and set(plot_df[col].dropna().unique()).issubset({0, 1})
    ]
 
    total = len(plot_df)
    proportions = []
    for col in binary_cols:
        counts = {
            "0_low":  ((plot_df[col] == 0) & (plot_df["_group"] == "low")).sum()  / total * 100,
            "0_high": ((plot_df[col] == 0) & (plot_df["_group"] == "high")).sum() / total * 100,
            "1_low":  ((plot_df[col] == 1) & (plot_df["_group"] == "low")).sum()  / total * 100,
            "1_high": ((plot_df[col] == 1) & (plot_df["_group"] == "high")).sum() / total * 100,
        }
        proportions.append(counts)
 
    prop_df = pd.DataFrame(proportions, index=binary_cols)

    # ------------------------------------------------------------------
    # 2. Statistical testing — 2x2 contingency table per SV
    #
    #    Table layout:
    #                  Low entropy   High entropy
    #    SV present        a              b
    #    SV absent         c              d
    #
    #    Test choice:
    #      - Any expected cell < min_expected  →  Fisher's exact (two-sided)
    #      - All expected cells >= min_expected →  Chi-squared + Yates correction
    # ------------------------------------------------------------------
    raw_pvals = []
    test_used = []
 
    for col in binary_cols:
        a = int(((plot_df[col] == 1) & (plot_df["_group"] == "low")).sum())
        b = int(((plot_df[col] == 1) & (plot_df["_group"] == "high")).sum())
        c = int(((plot_df[col] == 0) & (plot_df["_group"] == "low")).sum())
        d = int(((plot_df[col] == 0) & (plot_df["_group"] == "high")).sum())
 
        table = np.array([[a, b], [c, d]])
 
        # Expected counts under independence
        row_sums = table.sum(axis=1, keepdims=True)
        col_sums = table.sum(axis=0, keepdims=True)
        n = table.sum()
        expected = (row_sums * col_sums) / n if n > 0 else table.astype(float)

        if (expected < min_expected).any():
            _, p = fisher_exact(table, alternative="two-sided")
            test_used.append("fisher")
        else:
            _, p, _, _ = chi2_contingency(table, correction=True)
            test_used.append("chi2")

        raw_pvals.append(p)
 
    # ------------------------------------------------------------------
    # 3. Multiple-testing correction
    # ------------------------------------------------------------------
    if len(raw_pvals) > 1:
        _, corrected_pvals, _, _ = multipletests(
            raw_pvals, method=correction_method
        )
    else:
        corrected_pvals = np.array(raw_pvals)

    def _significance_label(p: float) -> str:
        """Convert p-value to standard annotation string."""
        if p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"
 
    sig_labels = [_significance_label(p) for p in corrected_pvals]
 
    # ------------------------------------------------------------------
    # 4. Nature figure style
    #    Single-column width: 88 mm  (~3.46 in)
    #    Double-column width: 180 mm (~7.09 in)
    #    Use double-column for > 8 SVs so x-labels don't collide.
    # ------------------------------------------------------------------
    n_bars = len(prop_df)
    fig_width_in = 3.46 if n_bars <= 6 else 7.09
 
    matplotlib.rcParams.update({
        # --- font (Nature requires Helvetica or Arial) ---
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size":         7,
        "axes.titlesize":    8,
        "axes.labelsize":    7,
        "xtick.labelsize":   6,
        "ytick.labelsize":   6,
        "legend.fontsize":   6,
        # --- lines ---
        "axes.linewidth":    0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size":  2.5,
        "ytick.major.size":  2.5,
        "xtick.direction":   "out",
        "ytick.direction":   "out",
        # --- spines ---
        "axes.spines.top":   False,
        "axes.spines.right": False,
        # --- misc ---
        "figure.dpi":        300,
        "savefig.dpi":       300,
        "pdf.fonttype":      42,   # embed fonts as Type 1 (editable in Illustrator)
        "ps.fonttype":       42,
    })
 
    # ------------------------------------------------------------------
    # 5. CMYK-safe, colorblind-friendly palette
    #
    #    Absent  · Low entropy  → light grey
    #    Absent  · High entropy → mid grey
    #    Present · Low entropy  → muted blue
    #    Present · High entropy → deep blue
    #
    #    Blues: #1A4F8A (Nature dark) / #7BAFD4 (light)
    #    Greys: #C8C8C8 / #666666
    # ------------------------------------------------------------------
    colors = {
        "0_low":  "#C8C8C8",
        "0_high": "#666666",
        "1_low":  "#7BAFD4",
        "1_high": "#1A4F8A",
    }
 
    stack_order = ["0_low", "0_high", "1_low", "1_high"]
 
    # ------------------------------------------------------------------
    # 6. Draw figure
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(fig_width_in, fig_width_in * 0.55))
 
    x = np.arange(n_bars)
    bar_width = 0.65
    bottom = np.zeros(n_bars)

    # Raw counts for labels
    count_df = []

    for col in binary_cols:
        counts = {
            "0_low":  ((plot_df[col] == 0) & (plot_df["_group"] == "low")).sum(),
            "0_high": ((plot_df[col] == 0) & (plot_df["_group"] == "high")).sum(),
            "1_low":  ((plot_df[col] == 1) & (plot_df["_group"] == "low")).sum(),
            "1_high": ((plot_df[col] == 1) & (plot_df["_group"] == "high")).sum(),
        }
        count_df.append(counts)

    count_df = pd.DataFrame(count_df, index=binary_cols)

    for key in stack_order:

        bars = ax.bar(
            x,
            prop_df[key].values,
            bottom=bottom,
            color=colors[key],
            edgecolor="white",
            linewidth=0.3,
            width=bar_width,
        )

        # Add n labels inside each segment
        labels = [
            str(int(v)) if v > 3 else ""
            for v in count_df[key].values
        ]

        ax.bar_label(
            bars,
            labels=labels,
            label_type="center",
            fontsize=5,
            color="white" if "high" in key else "black",
        )

        bottom += prop_df[key].values

    # ------------------------------------------------------------------
    # 7. Significance annotations above each bar
    #    Asterisks bold in black; "ns" lighter grey to de-emphasise
    # ------------------------------------------------------------------
    for i, label in enumerate(sig_labels):
        bar_top = prop_df.loc[binary_cols[i], stack_order].sum()
        y_annot = bar_top + 1.5     # 1.5 pp clearance above bar top
 
        ax.text(
            i, y_annot,
            label,
            ha="center", va="bottom",
            fontsize=6,
            fontweight="bold" if label != "ns" else "normal",
            color="black"    if label != "ns" else "#999999",
        )
 
    # ------------------------------------------------------------------
    # 8. Axes formatting
    # ------------------------------------------------------------------
    ax.set_xticks(x)
    ax.set_xticklabels(prop_df.index, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.6, n_bars - 0.4)
    ax.set_ylim(0, 105)
 
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.set_ylabel(ylabel, labelpad=4)
    ax.set_xlabel(xlabel, labelpad=4)
 
    # Remove bottom spine gap at zero
    ax.spines["left"].set_bounds(0, 100)
    ax.spines["bottom"].set_position(("outward", 2))
    ax.spines["left"].set_position(("outward", 2))
 
    # ------------------------------------------------------------------
    # 9. Title
    # ------------------------------------------------------------------
    ax.set_title(title, fontsize=8, fontweight="normal", pad=5, loc="left")
 
    # ------------------------------------------------------------------
    # 10. Panel label (a, b, …) — bold, slightly larger, upper-left
    # ------------------------------------------------------------------
    if panel_label is not None:
        ax.text(
            -0.12, 1.08,
            panel_label,
            transform=ax.transAxes,
            fontsize=9,
            fontweight="bold",
            va="top",
            ha="left",
        )
 
    # ------------------------------------------------------------------
    # 11. Legend — Nature style: no box, tight, upper-right
    # ------------------------------------------------------------------
    legend_labels = {
        "1_high": "Present, high entropy",
        "1_low":  "Present, low entropy",
        "0_high": "Absent, high entropy",
        "0_low":  "Absent, low entropy",
    }
 
    handles = [
        mpatches.Patch(facecolor=colors[k], edgecolor="none", label=legend_labels[k])
        for k in reversed(stack_order)   # present on top in legend
    ]
 
    ax.legend(
        handles=handles,
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        handlelength=1.2,
        handleheight=1.0,
        borderpad=0,
        labelspacing=0.3,
    )
 
    # ------------------------------------------------------------------
    # 12. Save — PNG (300 dpi) and TIFF (300 dpi, Nature submission)
    # ------------------------------------------------------------------
    os.makedirs(output_path, exist_ok=True)
 
    png_path  = os.path.join(output_path, f"{file_name}.png")
    pdf_path = os.path.join(output_path, f"{file_name}.pdf")
 
    plt.tight_layout(pad=0.4)

    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path, dpi=300)
    plt.close()


def create_histogram(output_path, list_of_data, labels, max_y, colors, xlabel, ylabel, title, metric_to_use, sequencing_type, file_name, bin_range=None):
    '''
    Parameters:
    -----------
        output_path: String
            Path output directory.

        list_of_data: List
            List of lists of entropy values.

        labels: List
            List of structural varient (sv) labels.

        bin_range: List
            List of evenly distributed integers.

        max_y: int
            Sets the top of the y-axis.

        colors: dict
            Dictionary of keys (color names) and values (hex codes).

        xlabel: String
            x-axis label.

        ylabel: String
            y-axis label.

        title: String
            Plot title.

        metric_to_use: String
            Parameter to show which metric was used from cna.seg file to use for entropy.

        sequencing_type: String
            Parameter to show which type of sequencing was used.

        file_name: String
            Plot file name.
         
    Function:
    ---------
        - Creates a histogram from list_of_data and labels and saves it to output_path.
        
    Returns:
    --------
        None
    '''
    plt.figure(figsize=(10,12))
    plt.hist(list_of_data, bins=bin_range, stacked=True, label=labels, color=colors)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()

    # Adds a little bit of padding to the y-axis
    plt.ylim(0, max_y * 1.05)
    
    boxplot_sub_dir = os.path.join(output_path, "complex_svs_vs_entropy_histograms")
    os.makedirs(boxplot_sub_dir, exist_ok=True)
    boxplot_file = os.path.join(boxplot_sub_dir, f"{file_name}.png")
    plt.tight_layout()
    plt.savefig(boxplot_file)
    plt.close()

def complex_sv_and_entropy_stacked_histogram_patient_level(df, output_path, sv_column, id_column, entropy_column, metric_to_use, sequencing_type):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            Dataframe with a column of structural variant (sv) calls and a column of chr8 entropy values from ulp data.
        
        output_path: String
            Path output directory.

        sv_column: String
            Name of sv column. 
        
        entropy_column: String
            Name of entropy column.

        metric_to_use: String
            Parameter to show which metric was used from cna.seg file to use for entropy.

        sequencing_type: String
            Parameter to show which type of sequencing was used.

    Function:
    ---------
        - Builds lists of data and labels as well as a colormap.
        - Creates a temp plot to standardize x- and y-axies across all histograms.
        - Creates a stacked histogram of all svs. 
        
    Returns:
    --------
        None
    '''
    # Build sv data set 
    sv_counts = df.groupby([id_column, sv_column,]).size().unstack(fill_value=0)
    entropy_values = df[[id_column, entropy_column]]
    entropy_values = entropy_values.drop_duplicates().set_index(id_column)
    entropy_values = entropy_values.sort_values(by=entropy_column)

    # Sort sv_counts based on entropy_values
    key = pd.Series({k: v for v, k in enumerate(entropy_values.index.unique())})
    sv_counts_sorted = sv_counts.sort_index(key=lambda x: key.reindex(x).values, na_position='last')

    print(f"sv_counts_sorted: \n{sv_counts_sorted}")

    sv_counts_sorted.index = entropy_values[entropy_column].round(2)

    print(f"sv_counts_sorted: \n{sv_counts_sorted}")


    colorblind_colors = [
        '#332288',  # Indigo
        '#117733',  # Green
        '#44AA99',  # Teal
        '#88CCEE',  # Cyan
        '#DDCC77',  # Sand
        '#CC6677',  # Rose
        '#AA4499',  # Purple
        '#882255',  # Wine
        '#661100',  # Brown
        '#999933'   # Olive
    ]
    
    # Create color map for each SV type
    sv_types = sv_counts.columns.tolist()
    color_map = {colorblind_colors[i % len(colorblind_colors)] for i in range(len(sv_types))}

    stacked_barplot(
        output_path, 
        sv_counts_sorted,
        entropy_values,
        color_map,
        'Patient ID', 
        'Frequency', 
        'Stacked Barplot of Chr8 Complex SVs per Patient', 
        metric_to_use, 
        sequencing_type,
        "patient_level_complex_svs_stacked_barplot"
    )

def stacked_barplot(output_path, data, txt_values, colors, xlabel, ylabel, title, metric_to_use, sequencing_type, file_name):
    '''
    Parameters:
    -----------
        output_path: String
            Path output directory.

        list_of_data: List
            List of lists of entropy values.

        labels: List
            List of structural varient (sv) labels.

        bin_range: List
            List of evenly distributed integers.

        max_y: int
            Sets the top of the y-axis.

        colors: dict
            Dictionary of keys (color names) and values (hex codes).

        xlabel: String
            x-axis label.

        ylabel: String
            y-axis label.

        title: String
            Plot title.

        metric_to_use: String
            Parameter to show which metric was used from cna.seg file to use for entropy.

        sequencing_type: String
            Parameter to show which type of sequencing was used.

        file_name: String
            Plot file name.
         
    Function:
    ---------
        - Creates a histogram from list_of_data and labels and saves it to output_path.
        
    Returns:
    --------
        None
    '''
    fig, ax = plt.subplots(figsize=(12,12))
    data.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.8)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()

    # Add value above each bar
    bar_totals = data.sum(axis=1)
    max_bar = bar_totals.max()

    # Set y-axis limit to 10% above highest bar
    ax.set_ylim(0, max_bar * 1.1)
    
    # Add text above each bar
    for i, (idx, value) in enumerate(txt_values.iterrows()):
        # Get the value to display (assuming single column in txt_values)
        text_val = idx.iloc[0] if hasattr(idx, 'iloc') else idx
        
        # Position text above the bar
        ax.text(i, bar_totals.iloc[i] + 0.02, f'{text_val}', 
                ha='center', va='bottom', fontsize=9, rotation=90)

    plt.xticks(rotation=45, ha='center')
    
    barplot_sub_dir = os.path.join(output_path, "complex_svs_vs_entropy_barplot")
    os.makedirs(barplot_sub_dir, exist_ok=True)
    boxplot_file = os.path.join(barplot_sub_dir, f"{file_name}.png")
    plt.tight_layout()
    plt.savefig(boxplot_file, dpi=300)
    plt.close()

def clean_and_fix_complex_sv_csv(csv_file, sample_id, caller, location_column, sv_column, entropy_column, mask):
    '''
    Parameters:
    -----------
        csv_file: String
            Dataframe with a column of structural variant (sv) calls and a column of chr8 entropy values from ulp data.
        
        sample_id: String
            Name of sample id column.

        caller: String
            Name of caller column. 

        location_column: String
            Name of location column. 

        sv_column: String
            Name of sv column. 
        
        entropy_column: String
            Name of entropy column.

        mask: List
            List of Strings, simple structural variants (sv), to filter out. 

    Function:
    ---------
        - Reads csv as a dataframe.
        - Filteres df based on mask so only complex svs are left (anything in mask will get filtered out).
        - Creates a new df with column names from parameters.
        - Because the location column in the old df has multiple locations of the call per patient, the for loop and creates a row in the new df for each location that starts with 8 (chromosome 8).
        
    Returns:
    --------
        df_new: pandas Dataframe
            Dataframe with only complex svs and a row for each sv event.
    '''
    # Read in csv file
    df = pd.read_csv(csv_file)

    # Filters out mask from df
    df_filtered = df[~df[sv_column].isin(mask)]
    df_filtered = df_filtered.reset_index(drop=True)

    return df_filtered        

def replace_chr8_entropy_column(output_path, complex_sv_df, entropy_file_name, metric_to_use, sequencing_type, entropy_column, cycle):
    '''
    Parameters:
    -----------
        output_path: String
            Path to the directory containing entropy tables.
        
        complex_sv_df: pandas DataFrame
            DataFrame with structural variant calls and a sample_id column.
        
        entropy_file_name: String
            Name of the entropy CSV file to read.
        
        metric_to_use: String
            Metric type used in the entropy table directory path.
        
        sequencing_type: String
            Sequencing type used in the entropy table directory path.
        
        entropy_column: String
            Name of the entropy column to extract from the entropy table.
        
        cycle: int or String
            Specific cycle number to filter entropy data.
    
    Function:
    ---------
        - Reads entropy table from the specified path.
        - Filters entropy data for the specified cycle and extracts patient_id, cycle, and entropy columns.
        - Splits sample_id column in complex_sv_df into patient_id and cycle components.
        - Merges complex_sv_df with filtered entropy data based on patient_id and cycle.
        - Replaces the chr8_entropy column with matched entropy values.
    
    Returns:
    --------
        complex_sv_df: pandas DataFrame
            Original DataFrame with updated chr8_entropy column containing merged entropy values.
    '''
    # Read in entropy df
    entropy_df_path = os.path.join(output_path, f'data-tables/entropy-tables-{sequencing_type}-{metric_to_use}/{entropy_file_name}')
    entropy_df = pd.read_csv(entropy_df_path)

    # Filter entropy df to patient id, specific cycle, and specific column
    entropy_df_filtered = entropy_df[["patient_id", "cycle", entropy_column]]
    entropy_df_filtered = entropy_df_filtered.drop(entropy_df_filtered[entropy_df_filtered.cycle != cycle].index)

    # Split patient_id column into sample id and cycle
    # Merge dfs together based on patient id and take chr column
    split_cols = complex_sv_df["sample_id"].str.split('_', expand=True)
    complex_sv_df = (
        complex_sv_df.assign(patient_id=split_cols[0], cycle=split_cols[2])
        .merge(entropy_df_filtered, on=["patient_id", "cycle"], how="left")
        .drop("sample_id", axis=1)
        )
    
    # Set column to new entropy values
    complex_sv_df["chr8_entropy"] = complex_sv_df[entropy_column]
    
    return complex_sv_df

def filter_gene_annotation_list(annotation_path, chr_column, chr_filter, gene_column_filter):
    '''
    Parameters:
    -----------
        annotation_path: String
            File path to annotation gene list.
        
        chr_column: String
            Name of column for chromosome in annotation gene list.

        chr_filter: List
            List of Specific chromosomes (int) to filter down to.

        gene_column_filter: String
            Name of column with gene names in it.

    Function:
    ---------
        - Read in annotation gene list.
        - Filter dataframe down to specific chromosome list and return.

    Return:
    -------
        filtered_df: pandas DataFrame
            filtered gene annotation list. 
    '''
    annotation_df = pd.read_table(annotation_path, sep='\t')
    mask = annotation_df[chr_column].isin(chr_filter)
    filtered_df = annotation_df[mask]
    filtered_unique_df = filtered_df.drop_duplicates(subset=[gene_column_filter])
    filtered_unique_df.reset_index(inplace=True, drop=True)
    return filtered_unique_df

def final_gene_matrix(annotation_df, annotation_gene_column, titan_matrix_path_list, id_column, cycle, titan_ploidy_path_list, ploidy_column):
    '''
    Parameters:
    -----------
        annotation_df: pandas DataFrame
            Filtered annotation gene list.

        annotation_gene_column: String
            Name of gene id column. 
        
        titan_matrix_path_list: List
            List of file paths to Titan gene matrix.

        id_column: String
            Name of patient id column for both files.

        cycle: String
            Cycle to filter patient ids down to.

        titan_ploidy_path_list: List
            List of file paths to Titan purity ploidy list.

        ploidy_column: String
            Name of column to get ploidy values in ploidy file.

        * Note: Make sure titan_matrix_path_list and titan_ploidy_path_list are respective to each other (eg. idexes of paths should align) *

    Function:
    ---------
        - Creates gene list to filter to.
        - For every matrix path in titan_matrix_path_list filter down to gene_list, append ploidy as a column, and add to joined_matrix.
        - Return joined_matrix

    Return:
    -------
        joined_matrix: pandas DataFrame
            A dataframe with all matrixes and ploidy values.
    '''
    # Create list of genes in annotated gene list
    gene_list = list(annotation_df[annotation_gene_column])

    # Instantiate empty dataframe
    joined_matrix = pd.DataFrame()

    # For each Titan matrix, filter down to gene_list, add ploidy as column, and add to joined_matrix
    for i, path in enumerate(titan_matrix_path_list):
        # Filter matrix to gene list
        matrix = pd.read_table(path, sep='\t')
        mask = matrix.columns.isin(gene_list)
        mask[0] = True # Make first index True to keep sample id
        matrix_filtered = matrix.loc[:, mask]

        # Remove samples not containing cycle and clean trailing tags on sample ids
        matrix_filtered_cleaned = clean_df_ids(matrix_filtered, id_column, cycle)

        # Add ploidy as a column to matrix
        ploidy_df = pd.read_table(titan_ploidy_path_list[i], sep='\t')
        ploidy_df_cleaned = clean_df_ids(ploidy_df, id_column, cycle)
        ploidy_df_filtered = ploidy_df_cleaned[[id_column, ploidy_column]]

        matrix_and_ploidy = pd.merge(matrix_filtered_cleaned, ploidy_df_filtered, on=id_column, how='left')

        joined_matrix = pd.concat([joined_matrix, matrix_and_ploidy], ignore_index=True)

    return joined_matrix

def normalize_rows_by_column(df, id_col_name, ref_col_name):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            DataFrame with columns matching id_col_name and ref_col_name.

        id_col_name: String
            Name of id column.

        ref_col_name: String
            Name of column to divide others by.

    Function:
    ---------
        - Sets id column to index.
        - Normalizes all rows by ref_col_name
        - Restore id column.
        - Return normalized df.

    Return:
    -------
        df: pandas DataFrame
            Dataframe with normalized rows.
    '''
    # Find all columns to normalize
    cols_to_normalize = df.columns.drop([id_col_name, ref_col_name])
    # Normalize columns
    df[cols_to_normalize] = df[cols_to_normalize].div(df[ref_col_name], axis=0).copy()

    return df
        
def clean_df_ids(df, id_column, cycle):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            Dataframe with and id column that has a tag with an underscore (eg. patientid_tag).
        
        id_column: String
            Name of id column.

        cycle: String
            Substring in tag of cycle to filter down to. 

    Function:
    ---------
        - Removes rows that doesnt contain cycle in id_column.
        - Removes tailing tag on id column and returns dataframe.

    Return:
    -------
        df: pandas DataFrame
            Dataframe with cleaned id column.
    '''
    # Remove rows that dont contain cycle
    df = df[df[id_column].str.contains(cycle, na=False)]

    # Removes trailing tags on id column
    df.loc[:, id_column] = df[id_column].str.split("_").str[0]

    return df

def get_entropy_column(output_path, entropy_file_path, chr_column, cycle, full_id=False):
    '''
    Parameters:
    -----------
        output_path: String
            Output path.

        entropy_file_path: String
            Path to entropy file.

        chr_column: String
            Name of entropy column to select.

        cycle: String
            Cycle to select for.

        full_id: boolean
            Combines id and cycle into id column.

    Function: 
    ---------
        - Load the entropy table.
        - Filter out any cycle that is not passed in through parameter.
        - Select for chr_column and return.

    Return:
    -------
        entropy_df_chr_column: pandas DataFrame
            A df of patient ids and respective entropy values.

    '''
    entropy_full_path = os.path.join(output_path, entropy_file_path)
    entropy_df = pd.read_csv(entropy_full_path)

    if cycle is not None:
        mask = entropy_df['cycle'].str.contains(cycle)
        entropy_df_filtered = entropy_df[mask]
    else:
        entropy_df_filtered = entropy_df

    if full_id:
        entropy_df_filtered['patient_id'] = entropy_df_filtered['patient_id'] + "_" + entropy_df_filtered['cycle']

    entropy_df_chr_column = entropy_df_filtered[['patient_id', chr_column]]

    entropy_df_chr_column.reset_index(inplace=True, drop=True)

    return entropy_df_chr_column

def find_significantly_different_genes_in_matrix(gene_matrix, gene_matrix_id_column, gene_matrix_ploidy_column, entropy_df, entropy_df_id_column, chr_column, responder_grouping):
    '''
    Parameters:
    -----------
        gene_matrix: pandas DataFrame
            Filtered Titan gene matrix as a df (patient_id x genes).

        gene_matrix_id_column: String
            Name of patient id column in gene_matrix.

        gene_matrix_ploidy_column: String
            Name of ploidy column in gene_matrix.

        entropy_df: pandas DataFrame
            DataFrame containing patient IDs and their corresponding entropy values as separate columns (default integar index).

        entropy_df_id_column: String
            Name of patient id column in entropy_df.

        chr_column: String
            Name of entropy column in gene_matrix.

        responder_grouping: String
            Specifies the grouping method to use for classifying responders and nonresponders.

    Function: 
    ---------
        - For gene_matrix divide all genes for each patient by respective ploidy and drop ploidy column.
        - Merge gene_matrix with entropy_df and remove genes with >20% missing values.
        - Then split the data into seperate dataframes depending on responder groupings.
        - Do a Mann Whitney-U test on responder vs nonresponder for each gene (only numeric columns).
        - Do FDR correction, keep only significant genes, and return.

    Return:
    -------
        results_df: pandas DataFrame
            Df with FDR corrected, significant genes (patient_id x genes).
        
        gene_matrix_with_groupings: pandas DataFrame
            Df used in Mann Whitney-U test. 

    '''
    # Divide all columns by ploidy column
    gene_matrix_normalized = pd.concat([
        gene_matrix[[gene_matrix_id_column]], 
        gene_matrix.drop(columns=[gene_matrix_id_column, gene_matrix_ploidy_column]).div(gene_matrix[gene_matrix_ploidy_column], axis=0)
    ], axis=1)

    # Merge dataframes
    gene_matrix_with_groupings = pd.merge(gene_matrix_normalized, entropy_df, left_on=gene_matrix_id_column, right_on=entropy_df_id_column, how='left').drop(columns=entropy_df_id_column)

    # Remove columns with > 20% NAs
    gene_matrix_with_groupings = gene_matrix_with_groupings.dropna(axis=1, thresh=int(0.8 * gene_matrix_with_groupings.shape[0]))

    # Build respnder and non responder groups
    if responder_grouping == "median":
        median_val = gene_matrix_with_groupings[chr_column].median()
        responder_group = gene_matrix_with_groupings[gene_matrix_with_groupings[chr_column] < median_val]
        non_responder_group = gene_matrix_with_groupings[gene_matrix_with_groupings[chr_column] >= median_val]

    # Identify numeric columns
    gene_columns = gene_matrix_with_groupings.select_dtypes(include='number').columns
    gene_columns = gene_columns.drop(chr_column, errors='ignore')

    # Do a Mann Whitney-U test on responders vs non responders for each gene
    results = {}
    for gene in gene_columns:
        responder = responder_group[gene].dropna()
        non_responder = non_responder_group[gene].dropna()

        if len(responder) > 1 and len(non_responder) > 1:
            mann_statistic, mann_pvalue = stats.mannwhitneyu(responder, non_responder, alternative='two-sided', nan_policy='omit')
            results[gene] = {"pvalue": mann_pvalue}

    # Convert results to dataframe
    results_df = pd.DataFrame(results).T

    # FDR correction
    pvalues = results_df['pvalue'].values
    _, pvalues_corrected, _, _ = multipletests(pvalues, method='fdr_bh')

    # Add corrected pvalues back to dataframe
    results_df['pvalues_corrected'] = pvalues_corrected

    genes_to_drop = []
    for gene, row in results_df.iterrows():
        if row["pvalues_corrected"] > 0.05:
            genes_to_drop.append(gene)
            
    results_df = results_df.drop(genes_to_drop)

    return results_df, gene_matrix_with_groupings

def get_tfx_column(output_path, pluvicto_master_sheet_path, id_column, tfx_column, encoding=None):
    '''
    Parameters:
    -----------
        output_path: String
            File path to output directory.

        pluvicto_master_sheet_path: String
            File path to pluvicto master sheet.

        id_column: String
            Name of patinet id column in pluvicto master sheet.

        tfx_column: String
            Name of TFx column in pluvicto master sheet.

        encoding: String
            To read csvs that need encoding.

    Funciton:
    ---------
        - Filter pluvicto master sheet down to id column and tfx column.

    Return:
    -------
        pluvicto_df[[id_column, tfx_column]]: pandas Dataframe
            Df with all patients and TFx as serperate columns (default integer index).

    '''
    pluvicto_full_path = os.path.join(output_path, pluvicto_master_sheet_path)

    pluvicto_df = pd.read_csv(pluvicto_full_path, encoding=encoding)

    return pluvicto_df[[id_column, tfx_column]]

def threshold_gene_matrix(gene_matrix, gene_matrix_id_column, gene_matrix_ploidy_column, entropy_df, entropy_df_id_column, chr_column, cn_threshold, direction, normalize_by_ploidy): 
    '''
    Parameters:
    -----------
        gene_matrix: pandas DataFrame
            Dataframe with gene CN levels and ploidy (patient_id x genes).
        
        gene_matrix_id_column: String
            Name of id column.

        gene_matrix_ploidy_column: String
            Name of ploidy column.
            
        entropy_df: pandas DataFrame
            Dataframe of patient id and respective entropy of specific chromosome as serperate columns (default integer index).

        entropy_df_id_column: String
            Name of id column.
 
        chr_column: String
            Name of entropy column.
 
        cn_threshold: int
            Threshold limit to be above or under the responder_grouping selected.

        direction: String
            Identifies gain or loss. 

        normalize_by_ploidy: boolean
            Choses if rows will be normalized by median of row or ploidy.

    Function:
    ---------
        - Drop ploidy column (its not needed and is an artifact from an older version of the function.)
        - Divide gene matrix by median (or ploidy) across each sample.
        - Merge entropy column onto gene matrix.
        - Remove all columns that has more than 20% NAs.
        - Binarize matrix on cn_threshold value (gain/loss vs other). 
        - Correct matrix and return. 

    Return:
    -------
        gene_matrix_binary: pandas DataFrame
            Binary matrix of genes with entropy as a column. 
    '''
    # Normalize rows
    if normalize_by_ploidy:
        # Divide rows by ploidy
        gene_matrix_normalized = normalize_rows_by_column(gene_matrix, gene_matrix_id_column, gene_matrix_ploidy_column)

        # Drop ploidy column
        gene_matrix_normalized.drop(columns=[gene_matrix_ploidy_column], inplace=True)
    else:
        # Drop ploidy column
        if gene_matrix_ploidy_column in gene_matrix.columns:
            gene_matrix.drop(columns=[gene_matrix_ploidy_column], inplace=True)

        # Divide all columns by median of row
        gene_matrix_normalized = pd.concat([
            gene_matrix[[gene_matrix_id_column]], 
            gene_matrix.drop(columns=[gene_matrix_id_column]).div(gene_matrix.drop(columns=[gene_matrix_id_column]).median(axis=1), axis=0)
        ], axis=1)

    # Merge dataframes
    gene_matrix_with_groupings = pd.merge(gene_matrix_normalized, entropy_df, left_on=gene_matrix_id_column, right_on=entropy_df_id_column, how='left').drop(columns=entropy_df_id_column)

    # Remove columns with > 20% NAs
    gene_matrix_with_groupings = gene_matrix_with_groupings.dropna(axis=1, thresh=int(0.8 * gene_matrix_with_groupings.shape[0]))
    
    # Make all columns numeric and then set all above/below threshold to 1 and copy neutral and all others to 0
    if direction == "gain":
        gene_matrix_binary = (gene_matrix_with_groupings.drop(columns=[gene_matrix_id_column, chr_column]).apply(pd.to_numeric, errors='coerce') > cn_threshold).astype(int)
    else:
        gene_matrix_binary = (gene_matrix_with_groupings.drop(columns=[gene_matrix_id_column, chr_column]).apply(pd.to_numeric, errors='coerce') < cn_threshold).astype(int)

    # ! Commented out code is for gain/loss vs copy neutral only
    # numeric_data = gene_matrix_with_groupings.drop(columns=[gene_matrix_id_column, chr_column]).apply(pd.to_numeric, errors='coerce')
    # gene_matrix_binary = pd.DataFrame(-1, index=numeric_data.index, columns=numeric_data.columns)
    # gene_matrix_binary[numeric_data == 1.0] = 0

    # if direction == "gain":
    #     gene_matrix_binary[numeric_data > cn_threshold] = 1
    # else:
        # gene_matrix_binary[numeric_data < cn_threshold] = 1
    
    # Merge to get correct order, then drop id column
    # gene_matrix_binary = pd.concat([gene_matrix_with_groupings[[gene_matrix_id_column]], gene_matrix_binary], axis=1)


    # Add back id and entorpy columns
    gene_matrix_binary[gene_matrix_id_column] = gene_matrix_with_groupings[gene_matrix_id_column].values
    gene_matrix_binary[chr_column] = gene_matrix_with_groupings[chr_column].values
    
    # reorder columns to put id column first
    cols = [gene_matrix_id_column] + [col for col in gene_matrix_binary.columns if col != gene_matrix_id_column]
    gene_matrix_binary = gene_matrix_binary[cols]

    return gene_matrix_binary

def logistic_regression_with_covariates(gene_matrix_binary, gene_matrix_id_column, tfx_df, chr_column, tfx_column, tfx_id_column, responder_grouping):
    '''
    Parameters:
    -----------
        gene_matrix_binary: pandas DataFrame
            Dataframe with gene expression and clinical data.
        
        gene_matrix_id_column: pandas DataFrame
            Dataframe with genes as index and pvalues as the column.

        tfx_df: pandas DataFrame
            Dataframe with sample ids as index and TFx as the column.

        chr_column: String
            Name of column for covariate.

        tfx_column: String
            Name of column for covariate.

        tfx_id_column: String
            Name of id column.

        responder_grouping: String
            Specifies the grouping method to use for classifying responders and nonresponders.

    Function:
    ---------
        - Binarizes gene_matrix_binary's 'response' column based on responder grouping method of entropy column, then drops entropy column.
        - Merges gene_matrix_binary and tfx_df.
        - Fits a logistic regression model for every gene. The model uses entropy value + TFx to predict CN gain/loss or not. 
        - Model removes genes that have one class for predictor values or response values.
        - Do FDR for both variable's pvalues. 
        - Return map with results and model specifics.

    Return:
    -------
        {
            'coefficients': results_df,
            'feature_columns': feature_columns,
            'X': X,
            'y': y,
            'sample_ids': modeling_df[gene_matrix_id_column].values
        }

    '''
    if responder_grouping == "median":
        # Add response label high entropy is 0 and low entropy is 1
        gene_matrix_binary['response'] = (gene_matrix_binary[chr_column] < gene_matrix_binary[chr_column].median()).astype(int)
    elif responder_grouping == "quartiles_extreme":
        q1 = gene_matrix_binary[chr_column].quantile(0.25)
        q3 = gene_matrix_binary[chr_column].quantile(0.75)
        gene_matrix_binary = gene_matrix_binary[(gene_matrix_binary[chr_column] <= q1) | (gene_matrix_binary[chr_column] >= q3)].copy()
        gene_matrix_binary['response'] = (gene_matrix_binary[chr_column] >= q3).astype(int)
    else:
        gene_matrix_binary['response'] = gene_matrix_binary[chr_column]

    # Create new df without entropy column
    modeling_df = gene_matrix_binary.drop(columns=chr_column)

    # Merge TFx column 
    how = 'left' if modeling_df.shape[0] < tfx_df.shape[0] else 'right'


    modeling_df = pd.merge(modeling_df, tfx_df, left_on=gene_matrix_id_column, right_on=tfx_id_column, how=how)

    if gene_matrix_id_column != tfx_id_column:
        modeling_df = modeling_df.drop(columns=tfx_id_column)

    # Drop NaNs if there are any 
    modeling_df.dropna(axis=0, inplace=True)
    modeling_df.reset_index(drop=True, inplace=True)

    # Create feature columns but remove entorpy column
    feature_columns = [col for col in modeling_df.columns if col not in [gene_matrix_id_column, 'response', tfx_column]]

    # Fit the regression with response predicting features
    results = []
    for feature in feature_columns:
        # Filter out all -1 rows 
        mask = modeling_df[feature] != -1
        modeling_df_filtered = modeling_df[mask]

        X = modeling_df_filtered[['response', tfx_column]]
        X = sm.add_constant(X)
        y = modeling_df_filtered[feature]

        # Try to fit the model, but if there is an error with perfect seperation (like y being all 1s or 0s) catch and log that result
        try:
            model = sm.Logit(y, X).fit(disp=0)

            # Get confidence intervals for coefficients
            conf_int = model.conf_int()

            # calculate accuracy for this feature
            y_pred = (model.predict(X) > 0.5).astype(int)
            accuracy = (y_pred == y).mean()

            entropy_response_coefficient = model.params['response']
            entropy_response_odds_ratio = np.exp(model.params['response'])
            entropy_response_or_ci_lower = np.exp(conf_int.loc['response', 0])
            entropy_response_or_ci_upper = np.exp(conf_int.loc['response', 1])

            entropy_response_pvalue = model.pvalues['response']
            tfx_coefficient = model.params[tfx_column]
            tfx_odds_ratio = np.exp(model.params[tfx_column])
            tfx_or_ci_lower = np.exp(conf_int.loc[tfx_column, 0])
            tfx_or_ci_upper = np.exp(conf_int.loc[tfx_column, 1])
            tfx_response_pvalue = model.pvalues[tfx_column]
            note = ""

        except (PerfectSeparationWarning, ConvergenceWarning, LinAlgError, RuntimeWarning)  as e:
            entropy_response_coefficient = np.nan
            entropy_response_odds_ratio = np.nan
            entropy_response_or_ci_lower = np.nan
            entropy_response_or_ci_upper = np.nan
            entropy_response_pvalue = np.nan
            tfx_coefficient = np.nan
            tfx_odds_ratio = np.nan
            tfx_or_ci_lower = np.nan
            tfx_or_ci_upper = np.nan
            tfx_response_pvalue = np.nan
            accuracy = np.nan
            note = f"{type(e).__name__}: feature has only one class"

        results.append({
            'feature': feature,
            'entropy_response_coefficient': entropy_response_coefficient,
            'entropy_response_odds_ratio': entropy_response_odds_ratio,
            'entropy_response_or_ci_lower': entropy_response_or_ci_lower,
            'entropy_response_or_ci_upper': entropy_response_or_ci_upper,
            'entropy_response_pvalue': entropy_response_pvalue,
            'entropy_response_adj_pvalue': np.nan,
            'tfx_coefficient': tfx_coefficient,
            'tfx_odds_ratio': tfx_odds_ratio,
            'tfx_or_ci_lower': tfx_or_ci_lower,
            'tfx_or_ci_upper': tfx_or_ci_upper,
            'tfx_response_pvalue': tfx_response_pvalue,
            'tfx_response_adj_pvalue': np.nan,
            'accuracy': accuracy,
            'note': note
        })

    results_df = pd.DataFrame(results)

    # Create an adj pvalue for entropy and tfx pvalues
    for col in ['entropy_response_pvalue', 'tfx_response_pvalue']:
        # Create mask for non NA values
        mask = results_df[col].notna()

        # Apply FDR correction to valid pvalues
        valid_pvalues = results_df.loc[mask, col]
        adj_pvalues = multipletests(valid_pvalues, alpha=0.05, method='fdr_bh')[1]

        if col == 'entropy_response_pvalue':
            adj_col = 'entropy_response_adj_pvalue'
        else:
            adj_col = 'tfx_response_adj_pvalue'

        # Map adj pvalues back to origional df
        results_df.loc[mask, adj_col] = adj_pvalues

    # Sort by entropy_response_adj_pvalue
    results_df = results_df.sort_values('entropy_response_adj_pvalue', ascending=True)

    print(f"Average accuracy: {results_df['accuracy'].mean():.3f}")

    return {
        'coefficients': results_df,
        'feature_columns': feature_columns,
        'X': X,
        'y': y,
        'sample_ids': modeling_df[gene_matrix_id_column].values
    }
    
def linear_regression_with_covariates(gene_matrix_binary, gene_matrix_id_column, tfx_df, chr_column, tfx_column, tfx_id_column, responder_grouping):
    '''
    Parameters:
    -----------
        gene_matrix_binary: pandas DataFrame
            Dataframe with gene expression and clinical data.
        
        gene_matrix_id_column: pandas DataFrame
            Dataframe with genes as index and pvalues as the column.

        tfx_df: pandas DataFrame
            Dataframe with sample ids as index and TFx as the column.

        chr_column: String
            Name of column for covariate.

        tfx_column: String
            Name of column for covariate.

        tfx_id_column: String
            Name of id column.

        responder_grouping: String
            Specifies the grouping method to use for classifying responders and nonresponders.

    Function:
    ---------
        - Binarizes gene_matrix_binary's 'response' column based on responder grouping method of entropy column, then drops entropy column.
        - Merges gene_matrix_binary and tfx_df.
        - Fits a linear regression model for every gene. The model uses entropy value + TFx to predict CN gain/loss or not. 
        - Model removes genes that have one class for predictor values or response values.
        - Do FDR for both variable's pvalues. 
        - Return map with results and model specifics.

    Return:
    -------
        {
            'coefficients': results_df,
            'feature_columns': feature_columns,
            'X': X,
            'y': y,
            'sample_ids': modeling_df[gene_matrix_id_column].values
        }

    '''
    if responder_grouping == "median":
        # Add response label high entropy is 0 and low entropy is 1
        gene_matrix_binary['response'] = (gene_matrix_binary[chr_column] < gene_matrix_binary[chr_column].median()).astype(int)
    elif responder_grouping == "quartiles_extreme":
        q1 = gene_matrix_binary[chr_column].quantile(0.25)
        q3 = gene_matrix_binary[chr_column].quantile(0.75)
        gene_matrix_binary = gene_matrix_binary[(gene_matrix_binary[chr_column] <= q1) | (gene_matrix_binary[chr_column] >= q3)].copy()
        gene_matrix_binary['response'] = (gene_matrix_binary[chr_column] >= q3).astype(int)
    else:
        gene_matrix_binary['response'] = gene_matrix_binary[chr_column]

    # Create new df without entropy column
    modeling_df = gene_matrix_binary.drop(columns=chr_column)

    # Merge TFx column 
    how = 'left' if modeling_df.shape[0] < tfx_df.shape[0] else 'right'

    modeling_df = pd.merge(modeling_df, tfx_df, left_on=gene_matrix_id_column, right_on=tfx_id_column, how=how)

    if gene_matrix_id_column != tfx_id_column:
        modeling_df = modeling_df.drop(columns=tfx_id_column)

    # Drop NaNs if there are any 
    modeling_df.dropna(axis=0, inplace=True)
    modeling_df.reset_index(drop=True, inplace=True)

    # Create feature columns but remove entorpy column
    feature_columns = [col for col in modeling_df.columns if col not in [gene_matrix_id_column, 'response', tfx_column]]

    # Fit the regression with response predicting features
    results = []
    for feature in feature_columns:
        # Filter out all -1 rows 
        mask = modeling_df[feature] != -1
        modeling_df_filtered = modeling_df[mask]

        X = modeling_df_filtered[[feature, tfx_column]]
        X = sm.add_constant(X)
        y = modeling_df_filtered['response']

        # Try to fit the model, but if there is an error with perfect seperation (like y being all 1s or 0s) catch and log that result
        try:
            model = sm.OLS(y, X).fit(disp=0)

            # Get confidence intervals for coefficients
            conf_int = model.conf_int()

            feature_response_coefficient = model.params[feature]
            feature_response_or_ci_lower = conf_int.loc[feature, 0]
            feature_response_or_ci_upper = conf_int.loc[feature, 1]
            feature_response_pvalue = model.pvalues[feature]

            tfx_coefficient = model.params[tfx_column]
            tfx_or_ci_lower = conf_int.loc[tfx_column, 0]
            tfx_or_ci_upper = conf_int.loc[tfx_column, 1]
            tfx_response_pvalue = model.pvalues[tfx_column]

        except (PerfectSeparationWarning, ConvergenceWarning, LinAlgError, RuntimeWarning)  as e:
            feature_response_coefficient = np.nan
            feature_response_odds_ratio = np.nan
            feature_response_or_ci_lower = np.nan
            feature_response_or_ci_upper = np.nan
            feature_response_pvalue = np.nan
            tfx_coefficient = np.nan
            tfx_odds_ratio = np.nan
            tfx_or_ci_lower = np.nan
            tfx_or_ci_upper = np.nan
            tfx_response_pvalue = np.nan
            accuracy = np.nan
            note = f"{type(e).__name__}: feature has only one class"

        results.append({
            'feature': feature,
            f'coefficient': feature_response_coefficient,
            f'ci_lower': feature_response_or_ci_lower,
            f'ci_upper': feature_response_or_ci_upper,
            f'pvalue': feature_response_pvalue,
            f'adj_pvalue': np.nan,
        })

        results.append({
            'feature': f"{feature}_TFx",
            'coefficient': tfx_coefficient,
            'ci_lower': tfx_or_ci_lower,
            'ci_upper': tfx_or_ci_upper,
            'pvalue': tfx_response_pvalue,
            'adj_pvalue': np.nan,
        })

    results_df = pd.DataFrame(results)

    # Create an adj pvalue for feature and tfx pvalues
    for col in ['pvalue']:
        # Create mask for non NA values
        mask = results_df[col].notna()

        # Apply FDR correction to valid pvalues
        valid_pvalues = results_df.loc[mask, col]
        adj_pvalues = multipletests(valid_pvalues, alpha=0.05, method='fdr_bh')[1]

        adj_col = 'adj_pvalue'
        
        # Map adj pvalues back to origional df
        results_df.loc[mask, adj_col] = adj_pvalues

    # Sort by entropy_response_adj_pvalue
    results_df = results_df.sort_values('adj_pvalue', ascending=True)

    return {
        'coefficients': results_df,
        'feature_columns': feature_columns,
        'X': X,
        'y': y,
        'sample_ids': modeling_df[gene_matrix_id_column].values
    }

def coefficient_barplot(output_path, sub_dir_path, model_results, coefficent_column, adj_pvalue_column, responder_grouping, cn_threshold, direction, normalize_by_ploidy):
    '''
    Parameters:
    -----------
        output_path: String
            Path to output dir.

        sub_dir_path: String
            Path to specific sub dir.

        model_results: pandas DataFrame
            Dataframe with coefficient information. 

        coefficent_column: String
            Name of coefficent column in model_results.

        adj_pvalue_column: String 
            Name of adjusted pvalue column in model_results.

        responder_grouping: String
            Specifies the grouping method to use for classifying responders and nonresponders.

        cn_threshold: int
            Threshold limit to be above or under the responder_grouping selected.

        direction: String
            Identifies gain or loss. 

        normalize_by_ploidy: boolean
            Shows if rows were normalized by median of row or ploidy.
        
    Function:
    ---------
        - Saves dataframe to csv.
        - Builds horizontal barplot for coefficients. 

    Return:
    -------
        None
    '''
    # Save model_results
    plot_output_path = os.path.join(output_path, sub_dir_path)
    os.makedirs(plot_output_path, exist_ok=True)
    model_results_path = os.path.join(plot_output_path, f"gene_matrix_and_entropy_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}_{direction}_vs_other.csv")
    model_results.to_csv(model_results_path)

    # Remove rows that arent significant
    model_results = model_results[model_results[adj_pvalue_column] < 0.05]
    
    top_features = model_results.head(20).copy()

    print(top_features)

    if top_features.empty:
        print("No significant features.")
        return

    # Create horizontal bar plot
    plt.figure(figsize=(10, 8))
    colors = ['red' if x < 0 else 'blue' for x in top_features[coefficent_column]]
    plt.barh(range(len(top_features)), top_features[coefficent_column], color=colors)
    plt.yticks(range(len(top_features)), top_features['feature'])
    plt.xlabel('Coefficient (Standardized) - high entropy (left) vs low entropy (right)')
    plt.title(f'Top 20 Features Associated with Treatment Response (Median Chr8 Entropy) - Copy Number {direction.title()}')
    plt.axvline(x=0, color='black', linestyle='--', linewidth=0.8)
    plt.gca().invert_yaxis()

    ax = plt.gca()

    # Get data range and ensure symmetric limits with at least 3 ticks per side of 0
    coef_min = top_features[coefficent_column].min()
    coef_max = top_features[coefficent_column].max()

    # Make sure limits are at least -1.0 to 1.0
    limit_min = min(coef_min - 0.2, -1.0)
    limit_max = max(coef_max + 0.2, 1.0)

    # Round to nearest 0.5
    limit_min = np.floor(limit_min / 0.5) * 0.5
    limit_max = np.ceil(limit_max / 0.5) * 0.5
    ax.set_xlim(limit_min, limit_max)

    # Set tick marks every 0.5
    ticks = np.arange(limit_min, limit_max + 0.1, 0.5)
    ax.set_xticks(ticks)

    # Position text boxes at bottom and within axis limits
    y_pos = len(top_features) - 1 
    x_min, x_max = ax.get_xlim()

    ax.text(x_min * 0.85, y_pos, 'Associated with\nHigh Entropy',
            ha='center', va='center', fontsize=11, color='darkred', 
            fontstyle='italic', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7)
        )
    
    ax.text(x_max * 0.6, y_pos, 'Associated with\nLow Entropy', 
            ha='center', va='center', fontsize=11, color='darkblue', 
            fontstyle='italic', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7)
            )

    plt.tight_layout()

    normalize_tag = ''
    if normalize_by_ploidy:
        normalize_tag = '_normalized_by_ploidy'

    plot_file_name = os.path.join(plot_output_path, f"gene_matrix_and_{coefficent_column}_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}{normalize_tag}_{direction}_vs_other_barplot.png")
    plt.savefig(plot_file_name)
    plt.close()

def odds_ratio_forest_plot(output_path, sub_dir_path, model_results, coefficent_column, adj_pvalue_column, responder_grouping, cn_threshold, direction, plot_title=None, file_name=None, gene_list=None):
    '''
    Parameters:
    -----------
        output_path: String
            Path to output dir.

        sub_dir_path: String
            Path to specific sub dir.

        model_results: pandas DataFrame
            Dataframe with coefficient information, sorted in ascending order of entropy adj pvalue.
            High entropy is coded as 0 and low entropy is coded as 1, meaning < 1 OR is high entorpy and > 1 OR is low entropy.

        coefficent_column: String
            Name of coefficent column in model_results.

        adj_pvalue_column: String 
            Name of adjusted pvalue column in model_results.

        responder_grouping: String
            Specifies the grouping method to use for classifying responders and nonresponders.

        cn_threshold: int
            Threshold limit to be above or under the responder_grouping selected.

        direction: String
            Identifies gain or loss. 

        plot_title: String (Default = None)
            Title of the plot.

        file_name: String (Default = None)
            Name of the file.

        gene_list: String (Default = None)
            List of genes to filter down to.

    Function:
    ---------
        - Builds a forest plot for odds ratios and confidence intervals. 

    Return:
    -------
        None
    '''

    # ----------------------------
    # Filter results
    # ----------------------------
    if gene_list:
        gene_list_title = "_prostate_gene_list"
        top_features = model_results[model_results['feature'].isin(gene_list)].copy()
    else:
        gene_list_title = ""
        top_features = model_results.nsmallest(20, adj_pvalue_column).copy()

    if top_features.empty:
        print("No significant features.")
        return

    top_features = top_features.sort_values(by=coefficent_column, ascending=True)

    top_features = top_features[
        (top_features['entropy_response_odds_ratio'] > 0) &
        (top_features['entropy_response_or_ci_lower'] > 0) &
        (top_features['entropy_response_or_ci_upper'] > 0)
    ].copy()

    if top_features.empty:
        print("No valid features after filtering invalid OR/CI values.")
        return

    # ----------------------------
    # Extract + transform values
    # ----------------------------
    or_vals  = top_features['entropy_response_odds_ratio'].values
    ci_lower = top_features['entropy_response_or_ci_lower'].values
    ci_upper = top_features['entropy_response_or_ci_upper'].values

    # lower_error = or_vals - ci_lower
    # upper_error = ci_upper - or_vals

    log_or    = np.log(or_vals)
    log_lower = np.log(ci_lower)
    log_upper = np.log(ci_upper)

    lower_error = log_or - log_lower
    upper_error = log_upper - log_or

    y_pos = np.arange(len(top_features))
    pvals = top_features[adj_pvalue_column].values

    # ----------------------------
    # Two-panel figure
    # ----------------------------
    fig, (ax_plot, ax_table) = plt.subplots(
        1, 2,
        figsize=(12, max(5, len(top_features) * 0.45 + 1.5)),
        gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.05}
    )

    # ── Left panel: forest plot ──
    for i, (lor, le, ue, p) in enumerate(zip(log_or, lower_error, upper_error, pvals)):
        facecolor = 'black' if p < 0.05 else 'white'
        ax_plot.errorbar(
            lor, i,
            xerr=[[le], [ue]],
            fmt='o',
            markersize=6,
            markerfacecolor=facecolor,
            markeredgecolor='black',
            ecolor='black',
            elinewidth=1.5,
            capsize=3
        )

    ax_plot.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax_plot.set_yticks(y_pos)
    ax_plot.set_yticklabels(top_features['feature'], fontsize=10)
    ax_plot.set_ylim(-0.8, len(top_features) - 0.2)
    ax_plot.set_xlabel('log(OR)', fontsize=11)
    ax_plot.set_ylabel('Gene', fontsize=11)
    ax_plot.spines['top'].set_visible(False)
    ax_plot.spines['right'].set_visible(False)
    ax_plot.grid(axis='x', linestyle=':', alpha=0.4)

    # Legend inside the plot panel (bottom right — now unobstructed)
    legend_elements = [
        Line2D([0], [0], marker='o', color='black', linestyle='None',
               markersize=6, markerfacecolor='black', label='adj p < 0.05 (FDR)'),
        Line2D([0], [0], marker='o', color='black', linestyle='None',
               markersize=6, markerfacecolor='white', label='not significant'),
    ]
    ax_plot.legend(handles=legend_elements, frameon=False, loc='best', fontsize=9)

    if plot_title:
        ax_plot.set_title(plot_title, fontsize=12, fontweight='bold')

    # ── Right panel: annotation table ──
    ax_table.set_xlim(0, 1)
    ax_table.set_ylim(-0.8, len(top_features))
    ax_table.axis('off')

    # Column headers as axis title — always sits just above the panel
    ax_table.text(0.05, len(top_features) - 0.2, 'log(OR)', fontsize=9, fontweight='bold', va='bottom')
    ax_table.text(0.20, len(top_features) - 0.2, '(95% CI)', fontsize=9, fontweight='bold', va='bottom')
    # ax_table.text(0.35, len(top_features) - 0.2, 'adj p', fontsize=9, fontweight='bold', va='bottom')

    # Plot values in OR space (not log(OR)) for each row
    for i, (p, or_val, ci_l, ci_u) in enumerate(zip(pvals, log_or, log_lower, log_upper)):
        weight = 'bold' if p < 0.05 else 'normal'
        ax_table.text(0.05, i, f'{or_val:.2f}', va='center', fontsize=9, fontweight=weight)
        ax_table.text(0.20, i, f'[{ci_l:.2f}, {ci_u:.2f}]', va='center', fontsize=9, fontweight=weight)
        # ax_table.text(0.35, i, f'{p:.2e}', va='center', fontsize=9, fontweight=weight)

    # ----------------------------
    # Save
    # ----------------------------
    plot_output_path = os.path.join(output_path, sub_dir_path)
    os.makedirs(plot_output_path, exist_ok=True)
    plt.savefig(os.path.join(plot_output_path, f"{file_name}.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(plot_output_path, f"{file_name}.pdf"), dpi=300, bbox_inches='tight')
    plt.close()

def coefficients_forest_plot(output_path, df, feature_column, coefficent_column, ci_lower_column, ci_upper_column, adj_pvalue_column, x_axis_title, y_axis_title, plot_title, file_name, gene_list=None):
    '''
    Parameters:
    -----------
        output_path: String
            Path to full output dir.

        df: pandas DataFrame
            Dataframe with coefficient information, sorted in ascending order of entropy adj pvalue.

        feature_column: String
            Name of feature column in model_results.

        coefficent_column: String
            Name of coefficent column in model_results.

        ci_lower_column: String
            Name of lower confidence interval column in model_results.

        ci_upper_column: String
            Name of upper confidence interval column in model_results.

        adj_pvalue_column: String 
            Name of adjusted pvalue column in model_results.

        plot_title: String
            Title of the plot.

        file_name: String
            Name of the file.

    Function:
    ---------
        - Builds a forest plot for odds ratios and confidence intervals. 

    Return:
    -------
        None
    '''

    # ----------------------------
    # Filter results
    # ----------------------------
    if gene_list:
        top_features = df[df['feature'].isin(gene_list)].copy()
    else:
        top_features = df.nsmallest(20, adj_pvalue_column).copy()

    if top_features.empty:
        print("No significant features.")
        return

    top_features = top_features.sort_values(by=coefficent_column, ascending=True)

    # ----------------------------
    # Extract + transform values
    # ----------------------------
    coef_vals  = top_features[coefficent_column].values
    ci_lower = top_features[ci_lower_column].values
    ci_upper = top_features[ci_upper_column].values

    lower_error = coef_vals - ci_lower
    upper_error = ci_upper - coef_vals

    y_pos = np.arange(len(top_features))
    pvals = top_features[adj_pvalue_column].values

    # ----------------------------
    # Two-panel figure
    # ----------------------------
    fig, (ax_plot, ax_table) = plt.subplots(
        1, 2,
        figsize=(12, max(5, len(top_features) * 0.45 + 1.5)),
        gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.05}
    )

    # ── Left panel: forest plot ──
    for i, (lor, le, ue, p) in enumerate(zip(coef_vals, lower_error, upper_error, pvals)):
        facecolor = 'black' if p < 0.05 else 'white'
        ax_plot.errorbar(
            lor, i,
            xerr=[[le], [ue]],
            fmt='o',
            markersize=6,
            markerfacecolor=facecolor,
            markeredgecolor='black',
            ecolor='black',
            elinewidth=1.5,
            capsize=3
        )

    ax_plot.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax_plot.set_yticks(y_pos)
    ax_plot.set_yticklabels(top_features[feature_column], fontsize=10)
    ax_plot.set_ylim(-0.8, len(top_features) - 0.2)
    ax_plot.set_xlabel(x_axis_title, fontsize=11)
    ax_plot.set_ylabel(y_axis_title, fontsize=11)
    ax_plot.spines['top'].set_visible(False)
    ax_plot.spines['right'].set_visible(False)
    ax_plot.grid(axis='x', linestyle=':', alpha=0.4)

    # Legend inside the plot panel (bottom right — now unobstructed)
    legend_elements = [
        Line2D([0], [0], marker='o', color='black', linestyle='None',
               markersize=6, markerfacecolor='black', label='adj p < 0.05 (FDR)'),
        Line2D([0], [0], marker='o', color='black', linestyle='None',
               markersize=6, markerfacecolor='white', label='not significant'),
    ]
    ax_plot.legend(handles=legend_elements, frameon=False, loc='best', fontsize=9)

    if plot_title:
        ax_plot.set_title(plot_title, fontsize=12, fontweight='bold')

    # ── Right panel: annotation table ──
    ax_table.set_xlim(0, 1)
    ax_table.set_ylim(-0.8, len(top_features))
    ax_table.axis('off')

    # Column headers as axis title — always sits just above the panel
    ax_table.text(0.05, len(top_features) - 0.2, 'β', fontsize=9, fontweight='bold', va='bottom')
    ax_table.text(0.20, len(top_features) - 0.2, '(95% CI)', fontsize=9, fontweight='bold', va='bottom')

    # Plot values in OR space (not log(OR)) for each row
    for i, (p, coef_val, ci_l, ci_u) in enumerate(zip(pvals, coef_vals, ci_lower, ci_upper)):
        weight = 'bold' if p < 0.05 else 'normal'
        ax_table.text(0.05, i, f'{coef_val:.2f}', va='center', fontsize=9, fontweight=weight)
        ax_table.text(0.20, i, f'[{ci_l:.2f}, {ci_u:.2f}]', va='center', fontsize=9, fontweight=weight)

    # ----------------------------
    # Save
    # ----------------------------
    os.makedirs(output_path, exist_ok=True)
    plt.savefig(os.path.join(output_path, f"{file_name}.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_path, f"{file_name}.pdf"), dpi=300, bbox_inches='tight')
    plt.close()

def get_prostate_specific_gene_list(gene_list_path, chr_column):
    '''
    Parameters:
    -----------
        gene_list_path: String
            Path to prostate specific gene list.

        chr_column: String
            Chromosome to filter down to.

    Function:
    ---------
        - Filters gene list to specified chromosome and returns it

    Return:
    -------
        gene_list: pandas DataFrame
            Dataframe with only specified chromosome.

    '''
    gene_list_pd = pd.read_table(gene_list_path, sep="\t", header=None)

    gene_list = gene_list_pd[gene_list_pd[1] == chr_column][0].to_list()

    return gene_list

def GSEA(output_path, sub_dir_path, model_results, gsea_ranking_column, adj_pvalue_column, gene_set, responder_grouping, cn_threshold, direction):
    '''
    Parameters:
    -----------
        output_path: String
            Path to output directory.

        sub_dir_path: String
            Path to specific sub dir to save plot to. 

        model_results: pandas DataFrame
            Dataframe with coefficient information, sorted in ascending order of entropy adj pvalue.

        coefficent_column: String
            Name of coefficient column to use.

        adj_pvalue_column: String
            Name of adj pvalue column to use.

        gene_set: String
            Gene set library to use.

        responder_grouping: String
            Responder grouping method.

        cn_threshold: String
            CN threshold used.

        direction: String
            CN gain or loss identifier. 

    Function: 
    ---------
        - Filter model to significant genes.
        - Create a gsea score of the 
    '''
    # Filter out NAs
    model_results = model_results[model_results[adj_pvalue_column].notna()]
    
    # Sort by gsea ranking column 
    ranked_genes = model_results[['feature', gsea_ranking_column]].sort_values(gsea_ranking_column, ascending=False)

    # names = gp.get_library_name(organism='Human')
    # names_filtered = [x for x in names if "hallmark" in x.lower()]
    # print(names_filtered)
    
    # Using gseapy prerank module
    pre_res = gp.prerank(
        rnk=ranked_genes,
        gene_sets=gene_set,
        organism="Human",
        threads=4,
        min_size=5,
        max_size=1000,
        permutation_num=1000,
        outdir=None,
        seed=6,
        verbose=True,
    )

    print(pre_res.res2d.head(5))

    # Save results to a csv
    plot_output_path = os.path.join(output_path, sub_dir_path)
    os.makedirs(plot_output_path, exist_ok=True)
    pre_res.res2d.sort_values("FDR q-val").to_csv(os.path.join(plot_output_path, f"gsea_on_gene_matrix_{gsea_ranking_column}_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}_{direction}_vs_other_{gene_set}_results.csv"))

    plot_file_name = os.path.join(plot_output_path, f"test_dotplot.png")


    # Plot gsea results
    ax = dotplot(
        pre_res.res2d,
        column="FDR q-val",
        title='Top Pathways in gene_set',
        cmap=plt.cm.viridis,
        size=6, # adjust dot size
        figsize=(10,10), 
        cutoff=0.5, 
        show_ring=False
    )
    
    fig = ax.get_figure()
    fig.tight_layout()
    plot_file_name = os.path.join(plot_output_path, f"gsea_on_gene_matrix_{gsea_ranking_column}_threshold_{cn_threshold}_abv_or_bel_{responder_grouping}_{direction}_vs_other_{gene_set}_dotplot.png")
    fig.savefig(plot_file_name)
    plt.close(fig)

def linear_mixed_effects_model(df, model_equation, id_col):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            Dataframe with coefficient information, sorted in ascending order of entropy adj pvalue.

        model_equation: String
            Linear mixed effects model equation (eg. 'chr8 ~ cycle + TFx')

        id_col: String
            Name of coefficient column to use.

        cycles_col: String
            Name of adj pvalue column to use.

        outcome_col: String
            Gene set library to use.

        fixed_effect_col: String
            Responder grouping method.

    Function: 
    ---------
        - 

    Return:
    --------
        model: statsmodel
            Linear mixed effects model.
    '''
    df = df.reset_index(drop=True)
    
    model = MixedLM.from_formula(
        model_equation,
        data=df.copy(),
        groups=df[id_col],
    )

    result = model.fit()

    return result

def linear_mixed_effects_forest_plot(model_result, output_path, file_name):
    '''
    Parameters:
    -----------
        model_result: stats model object

        output_path: String
            Path to output folder

        file_name: String
            Name of file.

    Function: 
    ---------
        - 

    Return:
    --------
        model: statsmodel
            Linear mixed effects model.
    '''

    # ----------------------------
    # Extract + transform values
    # ----------------------------
    print(model_result.summary())
    coef_vals = model_result.params
    cis       = model_result.conf_int()
    pvals     = model_result.pvalues

    # exclude intercept and Group Var
    mask   = ~coef_vals.index.str.contains("Intercept|Group")
    names  = coef_vals.index[mask].tolist()
    vals   = coef_vals[mask].values
    lo     = cis[0][mask].values
    hi     = cis[1][mask].values
    pvals  = pvals[mask].values

    lower_error = vals - lo
    upper_error = hi - vals
    y_pos       = np.arange(len(names))

    # ----------------------------
    # Two-panel figure
    # ----------------------------
    fig, (ax_plot, ax_table) = plt.subplots(
        1, 2,
        figsize=(12, max(5, len(names) * 0.45 + 1.5)),
        gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.05}
    )

    # ── Left panel: forest plot ──
    for i, (val, le, ue, p) in enumerate(zip(vals, lower_error, upper_error, pvals)):
        facecolor = 'black' if p < 0.05 else 'white'
        ax_plot.errorbar(
            val, i,
            xerr=[[le], [ue]],
            fmt='o',
            markersize=6,
            markerfacecolor=facecolor,
            markeredgecolor='black',
            ecolor='black',
            elinewidth=1.5,
            capsize=3
        )

    ax_plot.axvline(x=0, color='gray', linestyle='--', linewidth=1)
    ax_plot.set_yticks(y_pos)
    ax_plot.set_yticklabels(names, fontsize=10)
    ax_plot.set_ylim(-0.8, len(names) - 0.2)
    ax_plot.set_xlabel("Coefficient (β)", fontsize=11)
    ax_plot.set_ylabel("Predictor", fontsize=11)
    ax_plot.spines['top'].set_visible(False)
    ax_plot.spines['right'].set_visible(False)
    ax_plot.grid(axis='x', linestyle=':', alpha=0.4)

    legend_elements = [
        Line2D([0], [0], marker='o', color='black', linestyle='None',
               markersize=6, markerfacecolor='black', label='p < 0.05'),
        Line2D([0], [0], marker='o', color='black', linestyle='None',
               markersize=6, markerfacecolor='white', label='not significant'),
    ]
    ax_plot.legend(handles=legend_elements, frameon=False, loc='best', fontsize=9)
    ax_plot.set_title("Fixed effects forest plot", fontsize=12, fontweight='bold')

    # ── Right panel: annotation table ──
    ax_table.set_xlim(0, 1)
    ax_table.set_ylim(-0.8, len(names))
    ax_table.axis('off')

    ax_table.text(0.05, len(names) - 0.2, 'β',         fontsize=9, fontweight='bold', va='bottom')
    ax_table.text(0.25, len(names) - 0.2, '(95% CI)',  fontsize=9, fontweight='bold', va='bottom')
    ax_table.text(0.65, len(names) - 0.2, 'p-value',   fontsize=9, fontweight='bold', va='bottom')

    for i, (p, val, ci_l, ci_u) in enumerate(zip(pvals, vals, lo, hi)):
        weight = 'bold' if p < 0.05 else 'normal'
        ax_table.text(0.05, i, f'{val:.3f}',              va='center', fontsize=9, fontweight=weight)
        ax_table.text(0.25, i, f'[{ci_l:.3f}, {ci_u:.3f}]', va='center', fontsize=9, fontweight=weight)
        ax_table.text(0.65, i, f'{p:.3f}',                va='center', fontsize=9, fontweight=weight)

    # ----------------------------
    # Save
    # ----------------------------
    os.makedirs(output_path, exist_ok=True)
    final_path = os.path.join(output_path, file_name)
    plt.savefig(f'{final_path}.png', dpi=300)
    plt.savefig(f'{final_path}.pdf', dpi=300)
    plt.close()

def linear_mixed_effects_spaghetti_plot(df, id_col, x_col, y_col, response_col, x_label, y_label, model_result, model_pval_col, output_path, file_name):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            DataFrame model_results come from.

        id_col: String
            Name of id column.

        x_col: String
            Name of column to use on x-axis.

        y_col: String
            Name of column to use on y-axis.

        response_col: String
            Name of response column.

        x_label: String
            X-axis label.

        y_label: String
            Y-axis label.
        
        model_result: stats model object
            Object with model_pval_col in it.

        model_pval_col: String
            P-value for non-responder group compared to responder group.

        output_path: String
            Path to output folder

        file_name: String
            Name of file.

    Function: 
    ---------
        - 

    Return:
    --------
        model: statsmodel
            Linear mixed effects model.
    '''

    # -------------------------------------------------------------------------
    # Nature typography & style constants
    # -------------------------------------------------------------------------
    AXIS_LW       = 0.8              # spine / tick linewidth
    TICK_LEN      = 3.5
    TICK_LEN_MIN  = 2.0
 
    # Palette — Nature-approved, colorblind-safe (Wong 2011)
    COLOR_R   = "#0072B2"   # blue       → Responders
    COLOR_NR  = "#D55E00"   # red        → Non-responders
    TRAJ_ALPHA      = 0.12
    TRAJ_LW         = 0.55
    SCATTER_S       = 10
    SCATTER_ALPHA   = 0.22
    MEAN_LW         = 1.8
    CI_ALPHA        = 0.20
 
    # -------------------------------------------------------------------------
    # Data prep
    # -------------------------------------------------------------------------
    df = df.copy()
 
    if df[x_col].dtype == object:
        extracted = df[x_col].astype(str).str.extract(r"(\d+)").astype(float)
        df["_x_numeric"] = extracted
        plot_x = "_x_numeric"
        x_tick_labels = (
            df[[plot_x, x_col]].drop_duplicates().sort_values(plot_x)
        )
    else:
        plot_x = x_col
        x_tick_labels = None
 
    df = df.sort_values([id_col, plot_x])
 
    # -------------------------------------------------------------------------
    # Global rcParams — Nature style
    # -------------------------------------------------------------------------
    plt.rcParams.update({
        "font.size":          7,
        "axes.labelsize":     8,
        "axes.titlesize":     8,
        "axes.linewidth":     AXIS_LW,
        "xtick.labelsize":    7,
        "ytick.labelsize":    7,
        "xtick.major.width":  AXIS_LW,
        "ytick.major.width":  AXIS_LW,
        "xtick.minor.width":  AXIS_LW * 0.7,
        "ytick.minor.width":  AXIS_LW * 0.7,
        "xtick.major.size":   TICK_LEN,
        "ytick.major.size":   TICK_LEN,
        "xtick.minor.size":   TICK_LEN_MIN,
        "ytick.minor.size":   TICK_LEN_MIN,
        "xtick.direction":    "out",
        "ytick.direction":    "out",
        "legend.fontsize":    7,
        "legend.frameon":     False,
        "pdf.fonttype":       42,   # embeds fonts as TrueType in PDF
        "ps.fonttype":        42,
        "savefig.dpi":        600,
        "figure.dpi":         150,
    })
 
    # Nature single-column width = 89 mm ≈ 3.5 in; double = 183 mm ≈ 7.2 in
    fig, ax = plt.subplots(figsize=(3.5, 2.9))
 
    # -------------------------------------------------------------------------
    # Individual patient trajectories (coloured by response)
    # -------------------------------------------------------------------------
    response_colors = {0: COLOR_R, 1: COLOR_NR}
 
    for response_value, response_df in df.groupby(response_col):
        col = response_colors.get(response_value, "#555555")
 
        for _, patient_df in response_df.groupby(id_col):
            ax.plot(
                patient_df[plot_x],
                patient_df[y_col],
                lw=TRAJ_LW,
                alpha=TRAJ_ALPHA,
                color=col,
                solid_capstyle="round",
                zorder=2,
            )
            ax.scatter(
                patient_df[plot_x],
                patient_df[y_col],
                s=SCATTER_S,
                alpha=SCATTER_ALPHA,
                color=col,
                edgecolors="none",
                zorder=3,
            )
 
    # -------------------------------------------------------------------------
    # Group mean ± SEM trajectories
    # -------------------------------------------------------------------------
    response_labels = {0: "Responders", 1: "Non-responders"}
 
    for response_value, response_df in df.groupby(response_col):
        col   = response_colors.get(response_value, "#555555")
        label = response_labels.get(response_value, str(response_value))
 
        mean_df = (
            response_df.groupby(plot_x)[y_col]
            .agg(["mean", "sem"])
            .reset_index()
        )
 
        ax.plot(
            mean_df[plot_x],
            mean_df["mean"],
            lw=MEAN_LW,
            color=col,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=10,
        )
        ax.fill_between(
            mean_df[plot_x],
            mean_df["mean"] - mean_df["sem"],
            mean_df["mean"] + mean_df["sem"],
            alpha=CI_ALPHA,
            color=col,
            lw=0,
            zorder=9,
        )
 
    # -------------------------------------------------------------------------
    # Axis formatting — Nature style (no top / right spines)
    # -------------------------------------------------------------------------
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(AXIS_LW)
    ax.spines["bottom"].set_linewidth(AXIS_LW)
 
    ax.set_xlabel(x_label, labelpad=4)
    ax.set_ylabel(y_label, labelpad=4)
 
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=False))
    ax.tick_params(axis="both", which="both", pad=2)
 
    if x_tick_labels is not None:
        ax.set_xticks(x_tick_labels[plot_x])
        ax.set_xticklabels(x_tick_labels[x_col])
    else:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
 
    # Extend axes slightly beyond data range (no clip)
    ax.autoscale_view()
 
    # -------------------------------------------------------------------------
    # Statistics annotation — clean inset box
    # -------------------------------------------------------------------------
    try:
        pval = model_result.pvalues[model_pval_col]
 
        if pval < 0.0001:
            p_str = "p < 0.0001"
        elif pval < 0.001:
            p_str = "p < 0.001"
        elif pval < 0.01:
            p_str = "p < 0.01"
        else:
            p_str = f"p = {pval:.3f}"
 
        # Coefficient (β) with 95 % CI from model
        coef   = model_result.params[model_pval_col]
        ci_low, ci_high = model_result.conf_int().loc[model_pval_col]
 
        stat_lines = [
            f"Linear Mixed-Effects Model:",
            f"β = {coef:.3f} [{ci_low:.3f}, {ci_high:.3f}]",
            p_str,
        ]
        stat_text = "\n".join(stat_lines)
 
        ax.text(
            0.97, 0.97,
            stat_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=6,
            linespacing=1.5,
            color="#222222",
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="white",
                edgecolor="#CCCCCC",
                linewidth=0.5,
                alpha=0.92,
            ),
            zorder=20,
        )
 
    except Exception:
        pass
 
    # -------------------------------------------------------------------------
    # Legend — compact, Nature-style
    # -------------------------------------------------------------------------
    legend_handles = [
        Line2D(
            [0], [0],
            color=COLOR_R,
            lw=MEAN_LW,
            label="Responders",
            solid_capstyle="round",
        ),
        Line2D(
            [0], [0],
            color=COLOR_NR,
            lw=MEAN_LW,
            label="Non-responders",
            solid_capstyle="round",
        ),
    ]
 
    leg = ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.01, 0.86),
        handlelength=1.4,
        handleheight=0.6,
        handletextpad=0.4,
        labelspacing=0.3,
        borderpad=0,
        frameon=False,
        fontsize=7,
    )

    # -----------------------------
    # Tight layout
    # -----------------------------
    plt.tight_layout()

    # ----------------------------
    # Save
    # ----------------------------
    os.makedirs(output_path, exist_ok=True)
    final_path = os.path.join(output_path, file_name)
    plt.savefig(f'{final_path}.png', dpi=300)
    plt.savefig(f'{final_path}.pdf', dpi=300)
    plt.close()

def linear_mixed_effects_split_spaghetti_plot(df, id_col, x_col, y_col, response_col, x_label, y_label, model_result, model_pval_col, output_path, file_name):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            DataFrame model_results come from.

        id_col: String
            Name of id column.

        x_col: String
            Name of column to use on x-axis.

        y_col: String
            Name of column to use on y-axis.

        response_col: String
            Name of response column.

        x_label: String
            X-axis label.

        y_label: String
            Y-axis label.
        
        model_result: stats model object
            Object with model_pval_col in it.

        model_pval_col: String
            P-value for non-responder group compared to responder group.

        output_path: String
            Path to output folder

        file_name: String
            Name of file.

    Function: 
    ---------
        - 

    Return:
    --------
        model: statsmodel
            Linear mixed effects model.
    '''

    # -------------------------------------------------------------------------
    # Nature typography & style constants
    # -------------------------------------------------------------------------
    AXIS_LW       = 0.8              # spine / tick linewidth
    TICK_LEN      = 3.5
    TICK_LEN_MIN  = 2.0
 
    # Palette — Nature-approved, colorblind-safe (Wong 2011)
    COLOR_R   = "#0072B2"   # blue       → Responders
    COLOR_NR  = "#D55E00"   # red        → Non-responders
    TRAJ_ALPHA      = 0.12
    TRAJ_LW         = 0.55
    SCATTER_S       = 10
    SCATTER_ALPHA   = 0.22
    MEAN_LW         = 1.8
    CI_ALPHA        = 0.20
 
    # -------------------------------------------------------------------------
    # Data prep
    # -------------------------------------------------------------------------
    df = df.copy()
 
    if df[x_col].dtype == object:
        extracted = df[x_col].astype(str).str.extract(r"(\d+)").astype(float)
        df["_x_numeric"] = extracted
        plot_x = "_x_numeric"
        x_tick_labels = (
            df[[plot_x, x_col]].drop_duplicates().sort_values(plot_x)
        )
    else:
        plot_x = x_col
        x_tick_labels = None
 
    df = df.sort_values([id_col, plot_x])
 
    # -------------------------------------------------------------------------
    # Global rcParams — Nature style
    # -------------------------------------------------------------------------
    plt.rcParams.update({
        "font.size":          7,
        "axes.labelsize":     8,
        "axes.titlesize":     8,
        "axes.linewidth":     AXIS_LW,
        "xtick.labelsize":    7,
        "ytick.labelsize":    7,
        "xtick.major.width":  AXIS_LW,
        "ytick.major.width":  AXIS_LW,
        "xtick.minor.width":  AXIS_LW * 0.7,
        "ytick.minor.width":  AXIS_LW * 0.7,
        "xtick.major.size":   TICK_LEN,
        "ytick.major.size":   TICK_LEN,
        "xtick.minor.size":   TICK_LEN_MIN,
        "ytick.minor.size":   TICK_LEN_MIN,
        "xtick.direction":    "out",
        "ytick.direction":    "out",
        "legend.fontsize":    7,
        "legend.frameon":     False,
        "pdf.fonttype":       42,   # embeds fonts as TrueType in PDF
        "ps.fonttype":        42,
        "savefig.dpi":        600,
        "figure.dpi":         150,
    })
 
    # -------------------------------------------------------------------------
    # Split figure into responder and non-responder plots
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(5.6, 2.7),   # was 5.6
        sharey=True,
        constrained_layout=False
    )

    ax_map = {
        0: axes[0],   # Responders
        1: axes[1],   # Non-responders
    }

    # -------------------------------------------------------------------------
    # Individual patient trajectories + mean trajectories
    # -------------------------------------------------------------------------
    response_colors = {
        0: COLOR_R,
        1: COLOR_NR,
    }

    response_labels = {
        0: "Responders",
        1: "Non-responders",
    }

    for response_value, response_df in df.groupby(response_col):

        ax = ax_map[response_value]

        col   = response_colors.get(response_value, "#555555")
        label = response_labels.get(response_value, str(response_value))

        # -----------------------------------------------------------------
        # Individual patient trajectories
        # -----------------------------------------------------------------
        for _, patient_df in response_df.groupby(id_col):

            ax.plot(
                patient_df[plot_x],
                patient_df[y_col],
                lw=0.45,
                alpha=0.10,
                color=col,
                solid_capstyle="round",
                zorder=1,
            )

            ax.scatter(
                patient_df[plot_x],
                patient_df[y_col],
                s=8,
                alpha=0.16,
                color=col,
                edgecolors="none",
                zorder=2,
            )

        # -----------------------------------------------------------------
        # Group mean ± SEM
        # -----------------------------------------------------------------
        mean_df = (
            response_df.groupby(plot_x)[y_col]
            .agg(["mean", "sem"])
            .reset_index()
        )

        ax.plot(
            mean_df[plot_x],
            mean_df["mean"],
            lw=2.0,
            color=col,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=10,
        )

        # -----------------------------------------------------------------
        # Median value labels at line endpoints
        # -----------------------------------------------------------------
        start_x = mean_df[plot_x].iloc[0]
        end_x   = mean_df[plot_x].iloc[-1]

        start_y = mean_df["mean"].iloc[0]
        end_y   = mean_df["mean"].iloc[-1]

        ax.annotate(
            f"{start_y:.2f}",
            (start_x, start_y),
            xytext=(-4, 0),
            textcoords="offset points",
            fontsize=5.8,
            color=col,
            ha="right",
            va="center",
            fontweight="medium",
            zorder=15,
        )

        ax.annotate(
            f"{end_y:.2f}",
            (end_x, end_y),
            xytext=(4, 0),
            textcoords="offset points",
            fontsize=5.8,
            color=col,
            ha="left",
            va="center",
            fontweight="medium",
            zorder=15,
        )

        ax.fill_between(
            mean_df[plot_x],
            mean_df["mean"] - mean_df["sem"],
            mean_df["mean"] + mean_df["sem"],
            alpha=0.18,
            color=col,
            lw=0,
            zorder=9,
        )

        # -----------------------------------------------------------------
        # Panel titles
        # -----------------------------------------------------------------
        ax.set_title(
            label,
            fontsize=8,
            fontweight="bold",
            pad=6,
        )

        # -----------------------------------------------------------------
        # Axis styling
        # -----------------------------------------------------------------
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.spines["left"].set_linewidth(AXIS_LW)
        ax.spines["bottom"].set_linewidth(AXIS_LW)

        ax.tick_params(
            axis="both",
            which="both",
            pad=2,
        )

        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))

        # -----------------------------------------------------------------
        # Discrete x-axis values
        # -----------------------------------------------------------------
        if x_tick_labels is not None:

            ax.set_xticks(x_tick_labels[plot_x])
            ax.set_xticklabels(x_tick_labels[x_col])

        else:

            unique_x = sorted(df[plot_x].dropna().unique())
            ax.set_xticks(unique_x)

        # -----------------------------------------------------------------
        # Extra horizontal padding for endpoint labels
        # -----------------------------------------------------------------
        x_min = mean_df[plot_x].min()
        x_max = mean_df[plot_x].max()

        x_pad = (x_max - x_min) * 0.2

        ax.set_xlim(x_min - x_pad, x_max + x_pad)

        # -----------------------------------------------------------------
        # Remove duplicate y-axis clutter and removes right y-axis
        # -----------------------------------------------------------------
        if response_value == 1:

            # remove y-axis spine entirely
            ax.spines["left"].set_visible(False)

            # remove y ticks + labels
            ax.tick_params(
                axis="y",
                which="both",
                left=False,
                labelleft=False,
            )

        # -----------------------------------------------------------------
        # Slight padding around data
        # -----------------------------------------------------------------
        ax.margins(x=0.02)

    # -------------------------------------------------------------------------
    # Shared labels
    # -------------------------------------------------------------------------
    fig.supxlabel(
        x_label,
        fontsize=8,
        y=0.02,
    )

    fig.supylabel(
        y_label,
        fontsize=8,
        x=0.03,
    )

    # -------------------------------------------------------------------------
    # Shared title
    # -------------------------------------------------------------------------
    fig.suptitle(
        file_name.replace("_", " "),
        fontsize=9,
        fontweight="bold",
        y=0.98,
    )

    # -------------------------------------------------------------------------
    # Statistics annotation
    # -------------------------------------------------------------------------
    try:

        pval = model_result.pvalues[model_pval_col]

        if pval < 0.0001:
            p_str = "p < 0.0001"
        elif pval < 0.001:
            p_str = "p < 0.001"
        elif pval < 0.01:
            p_str = "p < 0.01"
        else:
            p_str = f"p = {pval:.3f}"

        coef = model_result.params[model_pval_col]
        ci_low, ci_high = model_result.conf_int().loc[model_pval_col]

        stat_text = (
            "Linear mixed-effects model\n"
            f"β = {coef:.3f} "
            f"[{ci_low:.3f}, {ci_high:.3f}]\n"
            f"{p_str}"
        )

        axes[1].text(
            0.98,
            0.98,
            stat_text,
            transform=axes[1].transAxes,
            ha="right",
            va="top",
            fontsize=5.8,
            linespacing=1.35,
            color="#222222",
            bbox=dict(
                boxstyle="round,pad=0.28",
                facecolor="white",
                edgecolor="#D0D0D0",
                linewidth=0.45,
                alpha=0.95,
            ),
            zorder=20,
        )

    except Exception:
        pass

    # -------------------------------------------------------------------------
    # Final layout tuning
    # -------------------------------------------------------------------------
    plt.subplots_adjust(
        # left=0.12,
        # right=0.98,
        bottom=0.22,
        top=0.84,
        wspace=0.01,
    )

    sns.despine(
        trim=True,
        left=False,
    )

    # ----------------------------
    # Save
    # ----------------------------
    os.makedirs(output_path, exist_ok=True)
    final_path = os.path.join(output_path, file_name)
    plt.savefig(f'{final_path}.png', dpi=300)
    plt.savefig(f'{final_path}.pdf', dpi=300)
    plt.close()

def longitudinal_violin_split_plot( df, id_col, x_col, y_col, response_col, x_label, y_label, model_result, model_pval_col, output_path, file_name, violin_width=0.35, bw_method="scott"):
    """
    Nature Medicine-style split violin plot for longitudinal data.
 
    Parameters
    ----------
    df : pandas DataFrame
    id_col : str        Subject identifier column.
    x_col : str         Timepoint column (string or numeric).
    y_col : str         Outcome variable column.
    response_col : str  Binary column: 0 = Responders, 1 = Non-responders.
    x_label : str       Shared x-axis label.
    y_label : str       Shared y-axis label.
    model_result        statsmodels result object (.pvalues, .params, .conf_int()).
    model_pval_col : str  Key for group-contrast coefficient in model_result.
    output_path : str   Output directory.
    file_name : str     File stem; underscores become spaces in the suptitle.
    violin_width : float  Half-width of each violin in data units (default 0.28).
    bw_method           KDE bandwidth passed to scipy gaussian_kde.
 
    Returns
    -------
    fig : matplotlib Figure
    """
 
    # ── Nature Medicine typographic constants ────────────────────────────────
    AXIS_LW      = 1
    TICK_LEN     = 4
    TICK_LEN_MIN = 1
 
    COLOR_R  = "#0072B2"   # blue  — Responders  
    COLOR_NR = "#D55E00"   # red   — Non-responders
 
    BOX_ALPHA    = 0.32
    MEDIAN_LW    = 1
    TREND_LW     = 1.4
    TREND_DOT_S  = 18
 
    # ── rcParams ─────────────────────────────────────────────────────────────
    plt.rcParams.update({
        "font.size":          7,
        "axes.labelsize":     7,
        "axes.titlesize":     7,
        "axes.linewidth":     AXIS_LW,
        "xtick.labelsize":    7,
        "ytick.labelsize":    7, 
        "xtick.major.width":  AXIS_LW,
        "ytick.major.width":  AXIS_LW,
        "xtick.minor.width":  AXIS_LW * 0.7,
        "ytick.minor.width":  AXIS_LW * 0.7,
        "xtick.major.size":   TICK_LEN,
        "ytick.major.size":   TICK_LEN,
        "xtick.minor.size":   TICK_LEN_MIN,
        "ytick.minor.size":   TICK_LEN_MIN,
        "xtick.direction":    "out",
        "ytick.direction":    "out",
        "legend.fontsize":    7,
        "legend.frameon":     False,
        "pdf.fonttype":       42,
        "ps.fonttype":        42,
        "savefig.dpi":        600,
        "figure.dpi":         150,
    })
 
    # ── Data prep ─────────────────────────────────────────────────────────────
    df = df.copy()
 
    if df[x_col].dtype == object:
        extracted = df[x_col].astype(str).str.extract(r"(\d+)").astype(float)
        df["_x_num"] = extracted[0]
        plot_x = "_x_num"
        x_tick_map = (
            df[[plot_x, x_col]]
            .drop_duplicates()
            .sort_values(plot_x)
            .set_index(plot_x)[x_col]
            .to_dict()
        )
    else:
        plot_x = x_col
        x_tick_map = None
 
    df = df.sort_values([id_col, plot_x])
    timepoints = sorted(df[plot_x].dropna().unique())
    
    # Compute dynamic y-limits based on actual data range
    all_y_values = df[y_col].dropna()
    y_min = all_y_values.min()
    y_max = all_y_values.max()
    y_lim = (y_min, y_max)
    
    def _draw_box(ax, x_center, values, color):
        values = np.asarray(values, dtype=float)
        values = values[~np.isnan(values)]
        if len(values) == 0:
            return np.nan

        q1, med, q3 = np.percentile(values, [25, 50, 75])
        iqr = q3 - q1

        box_hw = violin_width * 0.28

        # Box
        ax.add_patch(mpatches.Rectangle(
            (x_center - box_hw, q1),
            2 * box_hw,
            iqr,
            linewidth=0.8,
            edgecolor=color,
            facecolor=color,
            alpha=BOX_ALPHA,
            zorder=5,
        ))

        # Whiskers
        w_lo = max(values.min(), q1 - 1.5 * iqr)
        w_hi = min(values.max(), q3 + 1.5 * iqr)

        ax.plot([x_center, x_center], [w_lo, q1], color=color, lw=0.8, zorder=5)
        ax.plot([x_center, x_center], [q3, w_hi], color=color, lw=0.8, zorder=5)

        cap_hw = box_hw * 0.8
        ax.plot([x_center - cap_hw, x_center + cap_hw], [w_lo, w_lo], color=color, lw=0.8)
        ax.plot([x_center - cap_hw, x_center + cap_hw], [w_hi, w_hi], color=color, lw=0.8)

        # Median line
        ax.plot(
            [x_center - box_hw, x_center + box_hw],
            [med, med],
            color=color,
            lw=MEDIAN_LW,
            zorder=6,
        )

        # Jittered points
        rng = np.random.default_rng(42)  # reproducible jitter
        jitter_strength = box_hw * 0.35

        x_jitter = rng.uniform(
            x_center - jitter_strength,
            x_center + jitter_strength,
            size=len(values)
        )

        ax.scatter(
            x_jitter,
            values,
            s=7,                 # small, clean dots
            color=color,
            alpha=0.35,          # subtle
            edgecolors="none",
            zorder=4             # behind median line, above box fill
        )

        return med
 
    # ── Figure ────────────────────────────────────────────────────────────────
    # Double-column width (7.0 in), modest height
    fig, axes = plt.subplots(
        1, 2,
        figsize=(5, 2),
        sharey=True,
        constrained_layout=False,
    )
 
    ax_map         = {0: axes[0], 1: axes[1]}
    response_colors = {0: COLOR_R, 1: COLOR_NR}
    response_labels = {0: "Responders", 1: "Non-responders"}
 
    # ── Per-panel drawing ─────────────────────────────────────────────────────
    for resp_val, resp_df in df.groupby(response_col):
 
        ax    = ax_map[resp_val]
        col   = response_colors.get(resp_val, "#555555")
        label = response_labels.get(resp_val, str(resp_val))
 
        med_x, med_y = [], []
 
        for tp in timepoints:
            vals = resp_df.loc[resp_df[plot_x] == tp, y_col]
            med  = _draw_box(ax, tp, vals, col)
            if not np.isnan(med):
                med_x.append(tp)
                med_y.append(med)
 
        # Median trend line
        if len(med_x) >= 2:
            ax.plot(med_x, med_y,
                    color=col, lw=TREND_LW,
                    solid_capstyle="round", solid_joinstyle="round",
                    zorder=10)
            ax.scatter(med_x, med_y,
                       s=TREND_DOT_S, color=col,
                       zorder=11, edgecolors="none")
 
        # Endpoint labels (first and last median only)
        if med_x:
            ax.annotate(
                f"{med_y[0]:.2f}",
                (med_x[0], med_y[0]),
                xytext=(-12, 2), textcoords="offset points",
                fontsize=5.5, color=col,
                ha="right", va="center", zorder=15,
                fontweight="bold",
            )
            ax.annotate(
                f"{med_y[-1]:.2f}",
                (med_x[-1], med_y[-1]),
                xytext=(12, 2), textcoords="offset points",
                fontsize=5.5, color=col,
                ha="left", va="center", zorder=15,
                fontweight="bold",
            )
 
        # n label beneath panel title
        n = resp_df[id_col].nunique()
        ax.set_title(
            f"{label}\n$n$ = {n}",
            fontsize=7,
            fontweight="bold",
            color=col,
            pad=4,
            linespacing=1.5,
        )
 
        # Spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(AXIS_LW)
        ax.spines["bottom"].set_linewidth(AXIS_LW)
 
        # Ticks
        ax.tick_params(axis="both", which="both", pad=2)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=False))
        
        # Apply dynamic y-limits
        ax.set_ylim(y_lim)
 
        # X ticks
        ax.set_xticks(timepoints)
        if x_tick_map is not None:
            ax.set_xticklabels([x_tick_map.get(t, str(t)) for t in timepoints])
 
        # X limits — breathing room for endpoint labels
        x_range = (timepoints[-1] - timepoints[0]) if len(timepoints) > 1 else 1
        # x_pad  = x_range * 0.10
        ax.set_xlim(
            timepoints[0]  - violin_width,# - x_pad,
            timepoints[-1] + violin_width #+ x_pad,
        )
 
        # Remove right-panel y-axis clutter and middle y-axis line
        if resp_val == 1:
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", which="both", left=False, labelleft=False)
        
        # Ensure axes connect properly at origin
        ax.set_axisbelow(False)
        
        ax.margins(x=0.02)
 
    # ── Shared axis labels ────────────────────────────────────────────────────
    fig.supxlabel(x_label, fontsize=7, y=0.01, fontweight='bold')
    fig.supylabel(y_label, fontsize=7, x=0.02, fontweight='bold')
 
    # ── Layout and despine ────────────────────────────────────────────────────
    plt.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.18,
        top=0.85,
        wspace=0,
    )
    sns.despine(trim=True, left=False)

    # ----------------------------
    # Save
    # ----------------------------
    os.makedirs(output_path, exist_ok=True)
    final_path = os.path.join(output_path, file_name)
    plt.savefig(f'{final_path}.png', dpi=300)
    plt.savefig(f'{final_path}.pdf', dpi=300)
    plt.close()

def longitudinal_violin_split_plot_v2(
    df, id_col, x_col, y_col, response_col,
    x_label, y_label, model_result, model_pval_col,
    output_path, file_name, violin_width=0.35, bw_method="scott"
):
    """
    Nature Medicine-style split box plot for longitudinal data.
    response_col: 0 = Responders, 1 = Non-responders
    """
    COLOR_R, COLOR_NR, LW = "#0072B2", "#D55E00", 1

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7, "axes.linewidth": LW,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "xtick.major.size": 4, "ytick.major.size": 4,
        "xtick.major.width": LW, "ytick.major.width": LW, "xtick.direction": "out",
        "ytick.direction": "out", "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    # data prep
    df = df.copy()
    if df[x_col].dtype == object:
        df["_xn"] = df[x_col].str.extract(r"(\d+)").astype(float)
        x_tick_map = df[["_xn", x_col]].drop_duplicates().sort_values("_xn").set_index("_xn")[x_col].to_dict()
        plot_x = "_xn"
    else:
        plot_x, x_tick_map = x_col, None

    tps  = sorted(df[plot_x].dropna().unique())
    pos  = [i * 0.7 for i in range(len(tps))]
    df_r, df_nr = df[df[response_col] == 0], df[df[response_col] == 1]

    # y-axis
    v = df[y_col].dropna().values
    yr = v.max() - v.min()
    y_lo, y_hi = v.min() - 0.05 * yr, v.max() + 0.04 * yr
    step = min([0.05,0.1,0.15,0.2,0.25,0.5,1.0], key=lambda c: abs(c - yr/4))
    y_ticks = np.arange(np.ceil(v.min()/step)*step, y_hi, step)
    x_lo, x_hi = pos[0] - 0.42, pos[-1] + 0.42

    # p-value
    pval_str = ""
    if model_result is not None:
        p = model_result.pvalues[model_pval_col]
        pval_str = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"

    # figure
    fig = plt.figure(figsize=(3.5, 2.0))
    gs  = GridSpec(1, 2, figure=fig, left=0.14, right=0.97, top=0.84, bottom=0.16, wspace=0.28)
    ax_r  = fig.add_subplot(gs[0, 0])
    ax_nr = fig.add_subplot(gs[0, 1], sharey=ax_r)

    for ax, data, color, seed in [(ax_r, df_r, COLOR_R, 42), (ax_nr, df_nr, COLOR_NR, 43)]:
        rng, bw, medians = np.random.default_rng(seed), violin_width * 0.45, []

        for tp, px in zip(tps, pos):
            vals = data.loc[data[plot_x] == tp, y_col].dropna().values
            if not len(vals): medians.append(np.nan); continue

            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            iqr = q3 - q1
            w_lo = vals[vals >= q1 - 1.5*iqr].min()
            w_hi = vals[vals <= q3 + 1.5*iqr].max()
            medians.append(med)

            ax.add_patch(mpatches.FancyBboxPatch((px-bw, q1), 2*bw, iqr,
                boxstyle="square,pad=0", linewidth=0, facecolor=color, alpha=0.32, zorder=2))
            for y0, y1 in [(w_lo, q1), (q3, w_hi)]:
                ax.plot([px, px], [y0, y1], color=color, lw=LW, zorder=3)
            for yc in [w_lo, w_hi]:
                ax.plot([px-bw*0.7, px+bw*0.7], [yc, yc], color=color, lw=LW, zorder=3)
            ax.plot([px-bw, px+bw], [med, med], color=color, lw=LW, zorder=5)
            ax.scatter(px + rng.uniform(-bw*.55, bw*.55, len(vals)), vals,
                       s=5, color=color, alpha=0.45, linewidths=0, zorder=4)

        xs, ys = zip(*[(p,m) for p,m in zip(pos, medians) if not np.isnan(m)])
        ax.plot(xs, ys, color=color, lw=1.4, solid_capstyle="round", zorder=6)
        ax.scatter(xs, ys, s=18, color=color, linewidths=0, zorder=7)
        for i, (px, med) in enumerate(zip(pos, medians)):
            if not np.isnan(med):
                ax.text(px + (-0.06 if i==0 else 0.06), med, f"{med:.2f}",
                        ha=("right" if i==0 else "left"), va="center",
                        fontsize=7, fontweight="bold", color=color, zorder=8)

        ax.set_xticks(pos)
        ax.set_xticklabels([x_tick_map[t] for t in tps] if x_tick_map else [str(t) for t in tps])
        ax.set_yticks(y_ticks)
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
        ax.set_xlim(x_lo, x_hi); ax.set_ylim(y_lo, y_hi)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.spines["left"].set_position(("data", x_lo)); ax.spines["bottom"].set_position(("data", y_lo))
        ax.spines["left"].set_bounds(y_ticks[0], y_ticks[-1])
        ax.spines["bottom"].set_bounds(pos[0], pos[-1])
        ax.tick_params(axis="both", length=4, width=LW, direction="out", pad=2)

    ax_r.set_ylabel(y_label, fontsize=7, labelpad=4)
    ax_nr.tick_params(labelleft=False)
    if x_label: fig.text(0.555, 0.01, x_label, ha="center", va="bottom", fontsize=7)
    ax_r.set_title(f"Responders\n$n$ = {df_r[id_col].nunique()}",
                   fontsize=8, fontweight="bold", color=COLOR_R, pad=6)
    ax_nr.set_title(f"Non-responders\n$n$ = {df_nr[id_col].nunique()}",
                    fontsize=8, fontweight="bold", color=COLOR_NR, pad=6)
    fig.suptitle(file_name.replace("_"," ") + (f"  ({pval_str})" if pval_str else ""), fontsize=8, y=0.99)

    # ----------------------------
    # Save
    # ----------------------------
    os.makedirs(output_path, exist_ok=True)
    final_path = os.path.join(output_path, file_name)
    plt.savefig(f'{final_path}.png', dpi=300)
    plt.savefig(f'{final_path}.pdf', dpi=300)
    plt.close()