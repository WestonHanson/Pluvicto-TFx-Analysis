# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 12/12/2025
# Purpose: Prints the directory of a list of patients

from python_scripts.imports import *

# Function adopted from 'python_scripts/genome_instability_v1.py' process_FGA_TFx function
def find_patients_directory(genomic_instability_directory_array, patient_dict, target):
    '''
    Parameters:
    -----------
        genomic_instability_directory_array: List
            Array holding genomic instability directory paths.

        patient_dict: Dictionary
            Dictionary to store patient data with patient identifiers as keys.

        target: List
            List of patient ids (strings).
    
    Function:
    ---------
        - Loops through each directory in the genomic instability directory array.
        - If patient is in target it prints the directory and returns.

    Returns:
    --------
        None (prints directories)

    '''
    # For every directory in the array find the patients directory
    for directory in genomic_instability_directory_array:
        # For loop to find each patient in each batch
        for patient in os.listdir(directory):
            # print("Directory: " + patient)
            # Creates a full path with directory and patient identifier
            # If patient is not in the patient dictionary add the patient, else print that patient is already in the list
            if patient not in patient_dict:
                full_path = os.path.join(directory, patient)
                # If full path is true directory append patient identifier and genome altered to different lists
                if os.path.isdir(full_path):
                    if "Batch5" in full_path:
                        for dir in os.listdir(full_path):
                            if "solution" in dir:
                                full_path = os.path.join(full_path, dir)
                            
                
                    # Creates a dictionary of patient names with values genome altered and None
                    if patient in target:
                        print(f"\n{patient} dir: {full_path}\n")
                        target.remove(patient)
                        if len(target) == 0:
                            return
            else:
                print(f"{patient} already in list")

    print(f"{target} does not exist.")