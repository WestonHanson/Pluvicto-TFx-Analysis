from python_scripts.genome_instability_v1 import *
from python_scripts.central_depth_v2 import *

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
            if tag == "genomic instability file paths:":
                input_map["genomic_instability_files"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "tumor fraction file paths:":
                input_map["tumor_fraction_files"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "patient treatment progression:":
                input_map["patient_treatment_progression_files"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "central depth file path:":
                input_map["central_depth_files"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "time point:":
                input_map["time_point"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]

            elif tag == "genomic instability and tumor fraction data processing:":
                input_map["fga_ftx_processing"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "central depth data processing:":
                input_map["central_depth_processing"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "combine data tables:":
                input_map["combine_data_tables"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "genomic instability and tumor fraction scatter plots:":
                input_map["fga_ftx_scatter_plots"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "genomic instability and tumor fraction line plots:":
                input_map["fga_ftx_line_plots"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "central depth pca plot:":
                input_map["central_depth_pca_plots"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "central depth heatmap:":
                input_map["central_depth_heatmap"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "save heatmap data frames:":
                input_map["save_cd_heatmap_data_frames"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "central depth differential activity:":
                input_map["central_depth_differential_activity"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "central depth and tumor fraction density plot:":
                input_map["cd_tfx_density_plot"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "central depth and tumor fraction regression:":
                input_map["cd_tfx_regression"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "central depth and tumor fraction spearman rho plot:":
                input_map["cd_tfx_spearman_rho"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]
            elif tag == "cox forest plots of sig tfbs against clinical data:":
                input_map["cox_forest_tfbs_vs_clinical"] = [x.strip() for x in delimited_line[1:] if x.strip() != ""]


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

    # Gets all directories for searching ULP , tumor fraction, patient progression, and time point data from input file
    # genomic_instability_directory_array, tumor_fraction_directory_array, patient_progression_directory_array, central_depth_directory_array,time_point_array = get_input_file(input_path)
    input_map = get_input_file(input_path)

    # Processes samples if user selects True
    
    # ===============
    # Process 1:
    # ===============
    if input_map["fga_ftx_processing"][0] == "True":
        process_FGA_TFx(input_map["genomic_instability_files"], input_map["tumor_fraction_files"], patient_dict)
        create_FGA_data_table(input_map["patient_treatment_progression_files"], output_path, patient_dict)  

    # ===============
    # Process 2:
    # ===============
    if input_map["central_depth_processing"][0] == "True":
        process_central_depth(input_map["central_depth_files"], output_path)
    
    # ===============
    # Process 3:
    # ===============
    # Combines data tables if user selects True
    if input_map["combine_data_tables"][0] == "True":
        combine_data_tables(output_path)

    # ===============
    # Process 4:
    # ===============
    # Curreates plots based on users preference
    if input_map["fga_ftx_scatter_plots"][0] == "True":
        create_scatter_plot_from_table(output_path, input_map["time_point"], True, True, False, True)

    # ===============
    # Process 5:
    # ===============
    if input_map["fga_ftx_line_plots"][0] == "True":
        create_line_plot_over_time_points(output_path, input_map["time_point"])
    
    # ===============
    # Process 6:
    # ===============
    if input_map["central_depth_pca_plots"][0] == "True":
        create_central_depth_PCA_plot(output_path, 0, "central_depth_PCA_normalized_by_features")
        create_central_depth_PCA_plot(output_path, 1, "central_depth_PCA_normalized_by_samples")

    # ===============
    # Process 7:
    # ===============
    if input_map["central_depth_heatmap"][0] == "True":
        create_central_depth_heatmap(output_path, [], "None", -1, {1.0: "non-responder", 2.0: "non-responder", 3.0: "non-responder", 4.0: "responder", 5.0: "responder", 6.0: "responder"})
    
    # ===============
    # Process 8:
    # ===============
    if input_map["central_depth_differential_activity"][0] == "True":
        central_depth_differential_activity(output_path, 0.05, "sample level", -1, {1.0: "non-responder", 2.0: "non-responder", 3.0: "None", 4.0: "None", 5.0: "responder", 6.0: "responder"}, False)
    
    # ===============
    # Process 9:
    # ===============
    if input_map["cd_tfx_density_plot"][0] == "True":
        central_depth_vs_tumor_fraction_density_plot(output_path)
    
    # ===============
    # Process 10:
    # ===============
    if input_map["cd_tfx_regression"][0] == "True":
        central_depth_linear_regression(output_path, 0.2, "sample level", -1, {1.0: "non-responder", 2.0: "non-responder", 3.0: "None", 4.0: "None", 5.0: "responder", 6.0: "responder"})
    
    # ===============
    # Process 11:
    # ===============
    if input_map["cd_tfx_spearman_rho"][0] == "True":
        central_depth_and_tumor_fraction_spearman_rho(output_path, 0.05, "sample level", -1, {1.0: "non-responder", 2.0: "non-responder", 3.0: "None", 4.0: "None", 5.0: "responder", 6.0: "responder"})
    
    # ===============
    # Process 12:
    # ===============
    if input_map["cox_forest_tfbs_vs_clinical"][0] == "True":
        tfbs_list_vs_clinical_data_cox_forest_plots()


# If this script is run directly, it will execute the main function
# This allows the script to be run from the command line with arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plots tumor fraction vs fraction genome altered.")
    # parser.add_argument("--process1", default="True", help="Creates data table for tfx and FGA if True, else uses existing data table.")
    args = parser.parse_args()
    main(args)