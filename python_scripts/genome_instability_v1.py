# from central_depth_v1 import *

import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import csv
import pandas as pd
import seaborn as sns
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import sys
import re
from scipy import stats
from scipy.stats import iqr
from scipy.stats import gaussian_kde
from statsmodels.sandbox.stats.multicomp import multipletests
from statsmodels.stats.nonparametric import rank_compare_2indep
import statistics
import mpl_scatter_density
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LinearSegmentedColormap

# find_genome_alt function
# Takes in full path to patient directory and patient identifier
# Finds genomic instability per patient by counting bins that are not equal to 2 and divides that by total bins
# Returns genomic instability percent
def find_genome_alt(full_path, patient):

    # Finds seg file and opens it
    full_path = os.path.join(full_path, patient + ".cna.seg")
    file = open(full_path, "r")

    # Initialize an empty list to store the values from copy number corrected column
    copy_number_corrected_array = []

    # Read the file line by line
    for line in file:
        patient_data = line.strip().split("\t")
        patient_copy_number_corrected = patient_data[7]
        if "Corrected_Copy_Number" in patient_copy_number_corrected:
            continue
        # Append the value from copy number corrected column to the list
        copy_number_corrected_array.append(patient_copy_number_corrected)

    # Close the file
    file.close()

    # Counters for all bins and copy number alterations
    bins = 0
    CN_alts = 0

    # Check if every copy number value is 2 or not.
    # If value != 2 count it as a copy number alteration
    for value in copy_number_corrected_array:
        if value != "2":
            CN_alts += 1
        bins += 1

    # Return genomic instability percent 
    return CN_alts/bins


# process_FGA_TFx function
# Takes in genomic instability directory array (array holding genomic instability directory), 
# tumor fraction directory array (array holding tumor fraction directory), and patient dictionary
# Loops through each directory in the genomic instability directory array and adds patient identifier, genome altered, and tumor fraction values to the patient dictionary
def process_FGA_TFx(genomic_instability_directory_array, tumor_fraction_directory_array, patient_dict):
    print("Processing FGA_TFX...")
    # Tumor fraction data
    tumor_fraction_data = open(tumor_fraction_directory_array[0], "r")

    # For every directory in the array find the patients ULP data
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
                    # Creates a dictionary of patient names with values genome altered and None
                    number = patient.split('C')[-1]
                    patient_dict[patient] = (None, find_genome_alt(full_path, patient), number)
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


# create_FGA_data_table function
# Takes in patient progression cycle directory array (array holding patient progression cycle directory), output path, and patient dictionary
# Creates a data table with patient id (before first '_'), time point (C1, C2, etc.), tumor fraction, and genomic instability
# Writes the data table to a tsv file in the output path
def create_FGA_data_table(patient_progression_cycle_directory_array, output_path, patient_dict):
    # Sets output path and name
    tsvFileName = os.path.join(output_path, 'data-tables/FGA_data_table.tsv')

    # Open patient progression data and add it to a dictionary
    # patient_pregression_cycle_dict = {}
    # with open(patient_progression_cycle_directory_array[0], "r") as patient_progression:
    #     for line in patient_progression:
    #         patient = line.strip().split(",")
    #         if patient[3]:
    #             patient_id = patient[2].strip()
    #             patient_progression_cycle = patient[3]
    #             patient_pregression_cycle_dict[patient_id] = patient_progression_cycle

    # print(f"patient_pregression_cycle_dict: \n{patient_pregression_cycle_dict}")

    patient_pregression_cycle_dict = find_T_cycles(patient_dict)

    # Creates headers for data table
    tsv_header_name = ['patient_id', 'time_point', 'progression_cycle', 'tumor_fraction', 'genomic_instability']
    
    with open(tsvFileName, 'w') as tsvfile:
        tsv_writer = csv.writer(tsvfile, delimiter='\t')
        tsv_writer.writerow(tsv_header_name)

        # Breaks up dictionary into three parts and also splits patient id into the patient name and time point
        for patient_id, (tumor_fraction, genomic_instability, _) in patient_dict.items():
            patient_id_specific = patient_id[:patient_id.find("_")]
            patient_time_point = patient_id[patient_id.rfind("_")+1:]
            progression_cycle = patient_pregression_cycle_dict.get(patient_id_specific, "")
            tsv_writer.writerow([patient_id_specific, patient_time_point, progression_cycle, tumor_fraction, genomic_instability])
    
# find_T_cycles function
# Sorts patient_dict by key then by cycle value and finds the largest cycle each patient went to
def find_T_cycles(patient_dict):
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

# create_scatter_plot_from_table function
# Takes in output path and time point array (array of time points from input file to select from data table)
#  creates a scatter plot with tumor fraction on the x-axis and genome instability on the y-axis
# Filters out patients with tumor fraction less than 0.03
# Saves the plot as a png file in the output path
def create_scatter_plot_from_table(output_path, time_point_array, TFx, FGA, CD, above_3_percent):
    # Selects data table
    FGA_data_table = os.path.join(output_path, "data-tables/FGA_data_table.tsv")
    CD_data_table = os.path.join(output_path, "data-tables/central_depth_data_table.tsv")

    # Sets output plot name
    if TFx and FGA:
        output_plot = os.path.join(output_path, "tfx_vs_FGA_grt_0.03tfx.png")
    if TFx and CD:
       output_plot = os.path.join(output_path, "tfx_vs_CD.png")
    if FGA and CD:
        output_plot = os.path.join(output_path, "FGA_vs_CD.png")
    

    # Initialize empty lists to store patient id, tumor%, genome instability, and central depth
    patient_ids = []
    tumor_fraction_axis = []
    genome_altered_axis = []
    cd_axis = []

    # A switch to select all time points if time point array is empty
    switch = False
    if not time_point_array:
        switch = True
    
    with open(FGA_data_table, "r") as input_file:
        # Read the file line by line
        for line in input_file:
            patient_data = line.strip().split("\t")

            # Skips header
            if  patient_data[0] == "patient_id":
                continue


            # If the selected row is a time point measurement in the time point array or there are no time points selected
            if patient_data[1] in time_point_array or switch:
                patient_ids.append(patient_data[0])
                if float(patient_data[3]) <= 0.03 and above_3_percent: # Filters out all patients with tumor fraction less than specified amout
                    continue
                tumor_fraction_axis.append(float(patient_data[3]))
                genome_altered_axis.append(float(patient_data[4]))

    tfbs_counter = 0
    first_tfbs = ""
    if CD:
        with open(CD_data_table, "r") as input_file:
            for line in input_file:
                patient_data = line.strip().split("\t")
            
                # Skips header
                if  patient_data[0] == "patient_id":
                    continue

                if tfbs_counter == 0:
                    first_tfbs = patient_data[2]
                else:
                    if first_tfbs == patient_data[2]:
                        tfbs_counter = 0

                tfbs_counter += 1

                if patient_data[1] in time_point_array or switch:
                    cd_axis.append(float(patient_data[3]))

    print(f"first tfbs: {first_tfbs}")
    print(f"counter: {tfbs_counter}")


    x_axis = []
    y_axis = []
    if TFx and FGA:
        x_axis = tumor_fraction_axis
        y_axis = genome_altered_axis
    if TFx and CD:
        x_axis = tumor_fraction_axis
        y_axis = cd_axis
        x_axis = [x_axis[i // tfbs_counter] for i in range(len(y_axis))]
    if FGA and CD:
        x_axis = genome_altered_axis
        y_axis = cd_axis

    # Creates the scatter plot with tumor% as x-axis and genome instability as y-axis
    fig, ax = plt.subplots()
    ax.scatter(x_axis, y_axis)
    ax.set_xlabel("Tumor Fraction (%)") 
    ax.set_ylabel("Fraction Genome Altered (%)")
    ax.set_title("Fraction Genome Altered vs Tumor Fraction")
    ax.axvline(x=0.03, color="red", linestyle="--")
    ax.annotate(
        "0.03",
        xy=(0.03, 0), # x in data coords, y at bottom of y-axis
        xycoords='data',
        xytext=(0, -25), # shift downward in pixels
        textcoords='offset points',
        ha='center',
        va='top',
        fontsize=8,
        color='red'
    )
    ax.grid(True)

    # Save plot
    plt.savefig(output_plot, dpi=200) 
    plt.close() 


# create_line_plot_over_time_points function
# Takes in output path
# Reads the data table and filters out N/A values for tumor fraction and genomic instability
# Drops all rows with tumor fraction less than or equal to 0.03 and time point greater than progression cycle
# Groups table by treatment cycle then by time point
# Find mean of all tumor fraction and genomic instability values per time point
# Plot the 6 treatment cycle lines over time points with standard error of the mean confidence intervals
# Calls line_plot_creation function to create the line plot
def create_line_plot_over_time_points(output_path, time_point_array):
    data_table = os.path.join(output_path, "data-tables/FGA_data_table.tsv")

    df = pd.read_table(data_table)

    # Takes out N/A values
    df = df.dropna(subset=['tumor_fraction', 'genomic_instability', 'progression_cycle'])

    # Converts time_point to a numeric value
    df["time_point"] = df["time_point"].str.extract(r'C(\d+)').astype(int)

    # Creates individual line plots for each patient progression cycle
    create_individual_line_plot(df, time_point_array, output_path)

    # Remove all rows with tumor fraction less than or equal to 0.03
    df = df[df['tumor_fraction'] > 0.03]

    # Drops all rows where time point is greater than progression cycle
    df = df[df['time_point'].astype(float) <= df['progression_cycle']]

    # Creates a table for tumor fraction to plot
    tumor_fraction_summary = (
        df.groupby(['progression_cycle', 'time_point'])
        .agg(mean_tumor_fraction=('tumor_fraction', 'mean'),
            sem_tumor_fraction=('tumor_fraction', lambda x: np.std(x, ddof=1) / np.sqrt(len(x))),
            ci_upper=('tumor_fraction', lambda x: np.mean(x) + 1.96 * (np.std(x, ddof=1) / np.sqrt(len(x)))),
            ci_lower=('tumor_fraction', lambda x: np.mean(x) - 1.96 * (np.std(x, ddof=1) / np.sqrt(len(x))))
            )
        .reset_index()
    )

    print(f"tumor_fraction_summary: \n{tumor_fraction_summary}")

    # tumor_fraction_summary = (
    #     df.groupby(["progression_cycle", "time_point"])["tumor_fraction"]
    #     .mean()
    #     .reset_index()
    #     .sort_values(by=["progression_cycle", "time_point"])
    # )
    
    for index, row in df.iterrows():
        if row["patient_id"] == "patient_id":
            print("Found you")
        if row["progression_cycle"] == 1 and row["time_point"] == 1:
            print(row["patient_id"])
    print(df[(df["progression_cycle"] == 1) & (df["time_point"] == 1)]["tumor_fraction"].mean())

    print(tumor_fraction_summary)

    genomic_instability_summary = (
        df.groupby(['progression_cycle', 'time_point'])
        .agg(mean_genomic_instability=('genomic_instability', 'mean'),
            sem_genomic_instability=('genomic_instability', lambda x: np.std(x, ddof=1) / np.sqrt(len(x))),
            ci_upper=('genomic_instability', lambda x: np.mean(x) + 1.96 * (np.std(x, ddof=1) / np.sqrt(len(x)))),
            ci_lower=('genomic_instability', lambda x: np.mean(x) - 1.96 * (np.std(x, ddof=1) / np.sqrt(len(x))))
            )
        .reset_index()
    )
    
    print(f"genomic_instability_summary: \n{genomic_instability_summary}")


    line_plot_creation(tumor_fraction_summary, output_path, 'mean_tumor_fraction')
    line_plot_creation(genomic_instability_summary, output_path, 'mean_genomic_instability')


# line_plot_creation function
# Takes in data frame, summary data frame, output path, and mean column name to plot (either 'mean_tumor_fraction' or 'mean_genomic_instability')
# Creates a the line plot with the mean tumor fraction or genomic instability over time by treatment cycle
# Uses the summary data frame to plot the mean and confidence intervals
# Saves the plot as a png file in the output path with the appropriate name based on the mean column
def line_plot_creation(summary_df, output_path, mean_column):

    # Sets figure size and color palette
    plt.figure(figsize=(10, 6))
    palette = sns.color_palette("tab10", n_colors=summary_df['progression_cycle'].nunique())

    # Creates markers for each cycle
    marker_styles = ['o', 's', 'D', '^', 'v', 'x'] 

    legend_handles = []

    plt.rc("axes", labelsize=20)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rc('legend', fontsize=15)

    # Sorts the summary dataframe by progression cycle and creates an index i and a tuple of cycle and group
    for i, (cycle, group) in enumerate(summary_df.groupby('progression_cycle')):

        # Sorts the group by time point
        sorted_group = group.sort_values('time_point')

        # Sets the color by progression cycle and marker style
        color = palette[i]
        marker = marker_styles[i % len(marker_styles)]

        # Plots the mean line over time points with the specified color
        plt.plot(sorted_group["time_point"], 
                 sorted_group[mean_column],
                 label=f"Cycle {cycle}", 
                 color=color,
                 linewidth=2)
        
        # Plots the mean points at each time point with the marker style
        plt.scatter(sorted_group["time_point"], 
                    sorted_group[mean_column],
                    color=palette[i], 
                    s=50, 
                    zorder=5,
                    marker=marker)  # Use first character of time point as marker

        # Plots the confidence intervals as a filled area between the upper and lower bounds
        plt.fill_between(
            sorted_group["time_point"],
            sorted_group["ci_lower"],
            sorted_group["ci_upper"],
            alpha=0.2, color=palette[i]
        )

        # Creates a legend handle for the cycle
        handle = Line2D([0], [0],
                        color=color,
                        marker=marker,
                        linestyle='-',
                        linewidth=2,
                        markersize=8,
                        label=f"Cycle {cycle}")
        
        legend_handles.append(handle)

    # Sets the x and y labels, title, and legend
    plt.xlabel("Treatment Cycle")

    if mean_column == 'mean_tumor_fraction':
        plt.ylabel("Mean Tumor Fraction (%)", labelpad=5.0)
    else:
        plt.ylabel("Mean Fraction Genome Altered (%)",  labelpad=5.0)

    

    if mean_column == 'mean_tumor_fraction':
        plt.title("Tumor Fraction over Time by Number of Treatment Cycles Completed", fontsize=20)
    else:
        plt.title("Fraction Genome Altered over Time by Number of Treatment Cycles Completed", fontsize=23, pad=10)

    plt.legend(handles=legend_handles, title="Treatment Cycles", loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xticks(sorted(summary_df["time_point"].unique()))
    plt.grid(True)

    if mean_column == 'mean_tumor_fraction':
        plot_path = os.path.join(output_path, "mean_tfx_over_time_by_progression_cycle_SEM_filtered.png")
    else:
        plot_path = os.path.join(output_path, "mean_FGA_over_time_by_progression_cycle_SEM_filtered.png")

    # Put a legend to the right of the current axis
    # plt.legend()
    plt.tight_layout(rect=[0.05, 0, 1, 1.05])
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()


# create_individual_plot function
# Take in data frame (raw data table), data fram of low tumor fraction samples, time point array from input, and output path
# Creates individual line plots for each patient progression cycle
# Splits the data into two subplots: one for patients above 3% tumor fraction
# and one for patients at or below 3% tumor fraction
# Saves the plots in the output path with the appropriate name based on the progression
def create_individual_line_plot(df, time_point_array, output_path):
    if not time_point_array:
        time_point_array = ["C1", "C2", "C3", "C4", "C5", "C6"]

    # Convert time points to integers
    time_point_array = [int(tp[1:]) for tp in time_point_array]  

    # Get unique progression cycles
    unique_cycles = df['progression_cycle'].unique()

    print(f"unique_cycles: \n{unique_cycles}")

    column_of_interest_array = ['tumor_fraction', 'genomic_instability']

    # Classify patients into low and high tumor fraction based on the threshold of 0.03
    low_tfx_patients = set(df[df['tumor_fraction'] <= 0.03]['patient_id'].unique())
    high_tfx_patients = set(df['patient_id'].unique()) - low_tfx_patients
    print(f"low_tfx_patients: \n{low_tfx_patients}")
    print(f"high_tfx_patients: \n{high_tfx_patients}")


    plt.rc("axes", labelsize=20)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rc('legend', fontsize=12)

    # For each column of interest, create a separate figure
    # Splits the data into two subplots: one for patients above 3% tumor fraction, and one for patients at or below 3% tumor fraction
    for column in column_of_interest_array:

        # Creates a figure for each progression cycle
        for cycle in sorted(unique_cycles):

            # Creates two subplots for high and low tumor fraction patients
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

            # Filter the DataFrame for the current cycle
            cycle_df = df[df['progression_cycle'] == cycle]
            palette = sns.color_palette("tab10", len(cycle_df['patient_id'].unique()))

            # For each patient in the cycle, plot the data
            for i, (patient_id, patient_df) in enumerate(cycle_df.groupby('patient_id')):
                # Sorts the patient data frame by time point
                patient_df_sorted = patient_df.sort_values('time_point')
                x_vals = patient_df_sorted['time_point']
                y_vals = patient_df_sorted[column]

                color = palette[i % len(palette)]  # to avoid IndexError

                # Sets the subplot based on whether the patient is in high or low tumor fraction
                ax = ax1 if patient_id in high_tfx_patients else ax1

                ax.plot(x_vals, y_vals, label=patient_id, color=color, marker='o', linewidth=2)

                ax.annotate(
                    patient_id,
                    xy=(x_vals.iloc[-1], y_vals.iloc[-1]),
                    textcoords="offset points",
                    xytext=(0, 10),
                    ha='center',
                    fontsize=8,
                    color=color
                )

            # Sets the y-axis label and legends for each subplot (if figure is tfx then it adds a horizontal line at 0.03)
            for ax in [ax1]:
                ax.axhline(0.03, color='red', linestyle='--', linewidth=1, label='3% threshold')
                ax.set_ylabel(column.replace('_', ' ').title() + " (%)")

                # Only adds legend if there are handles (i.e., plotted lines) in the axis - this is to prevent warning messages
                handles, labels = ax.get_legend_handles_labels()
                if handles:
                    ax.legend(title="Patient ID", bbox_to_anchor=(1.05, 1), loc='upper left')
                
            fig.suptitle(f"{column.replace('_', ' ').title()} Over Time by Treatment Cycle {int(cycle)}", fontsize=16)
            ax1.set_title(f"Patients With Samples Above 3% Tumor Fraction", fontsize=18)
            # ax2.set_title(f"Patients With Samples At or Below 3% Tumor Fraction", fontsize=18)
            # ax2.set_xlabel("Treatment Cycle")

            plt.tight_layout()
            output_file = os.path.join(output_path, f"treatment_cycle_{cycle}_split_by_tfx_{column}_over_time.png")
            plt.savefig(output_file)
            plt.close()

