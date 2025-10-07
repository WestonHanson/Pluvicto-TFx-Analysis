from scripts.python_scripts.genome_instability_v1 import *

import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
import sys
import re
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
import statistics
from lifelines import KaplanMeierFitter, CoxPHFitter

# process_central_depth function
# Create a data frame from output path and save it to a tsv file
def process_central_depth(central_depth_directory_array, output_path):

    # Selects central depth directory
    central_depth_file = central_depth_directory_array[0]

    # Set output path
    central_depth_output = os.path.join(output_path, "data-tables/central_depth_data_table.tsv")

    # Reads table in 
    df = pd.read_table(central_depth_file)

    # Filter for site and feature
    df_flattened = df[df["feature"] == "central-depth"].copy()
    
    # Create df with columns 
    df_flattened['patient_id'] = df_flattened['sample'].apply(lambda x: re.match(r"(FH\d+|FHL\d+)_", x).group(1))
    df_flattened['time_point'] = df_flattened['sample'].apply(lambda x: re.search(r"_(C\d+)", x).group(1))
    df_flattened['central_depth'] = df_flattened['value']
    
    df_flattened['mean_depth'] = df_flattened['value']

    central_depth_df = df_flattened[["patient_id", "time_point", "site", "central_depth"]]

    central_depth_df.to_csv(central_depth_output, sep='\t', index=False)


# combine_data_tables function
# Takes in fraction genome altered table and central depth table
# Combines tables on patinet id and drops unnessarry columns
# Returns the combined data table
def combine_data_tables(output_path):

    # Reads in data table into data frame
    FGA_data_table_file_name = os.path.join(output_path, 'data-tables/FGA_data_table.tsv')
    CD_data_table_file_name = os.path.join(output_path, 'data-tables/central_depth_data_table.tsv')

    FGA_df = pd.read_table(FGA_data_table_file_name)
    CD_df = pd.read_table(CD_data_table_file_name)

    print(f"FGA_df: \n{FGA_df}")
    print(f"CD_df: \n{CD_df}")

    # Merges tables on patient_id the drops unnessary columns
    df = pd.merge(FGA_df, CD_df, on=['patient_id', 'time_point'], how="left")

    print(f"df: \n{df}")

    # Takes out N/A values
    df = df.dropna(subset=['progression_cycle']).reset_index(drop=True)
    df = df.sort_values(by=['patient_id', 'time_point']).reset_index(drop=True)

    # Converts time_point to a numeric value
    df["time_point"] = df["time_point"].str.extract(r'C(\d+)').astype(int)

    # Saves combined data frame
    df_output = os.path.join(output_path, "data-tables/combined_data_table.tsv")
    df.to_csv(df_output, sep='\t', index=False)


def create_central_depth_PCA_plot(output_path, normalize_axis, output_title):
    """
    Parameters:
    -----------
        output_path: str
            Path to output directory
        normalize_axis: str
            Specifies what axis to normalize at: "sample level", "feature level", or "None"
        output_title: str
            File name for PCA plot

    Function:
    ---------
        - Creates a PCA plot from combined_data_table.tsv.
        - Only uses samples at time point 1. 
        - Pivots data table, normalizes it, then does PCA on it
        - Creates a new data frame with meta data
        - Conditions tumor fraction into low, medium, and high
        - Calls scatter plot on new data frame and saves plot
        
    Returns:
    --------
        Void
    """
    
    # Creates name of combined data table
    combined_data_table_path = os.path.join(output_path, 'data-tables/combined_data_table.tsv')

    # Reads in data table into data frame
    df = pd.read_table(combined_data_table_path)

    # Selects for just C1 patients
    df = df[df["time_point"] == 1].reset_index(drop=True)

    # Pivots table around sites so that the patient_id are the rows and sites are the columns
    df_pivot = _pivot_table_(df, "central_depth", "patient_id", "site", 'mean', 0)

    # turns pivot table into data frame by dropping the index and column name
    pca_input = df_pivot.reset_index()

    # Standard normailization across features -> returns a data frame
    df_normalized, patient_id_list = _normalize_(pca_input, normalize_axis)

    pca_df = _PCA_(df_normalized)

    # Appiles principle components to the normalized data, and creates a new data frame from it
    df_normalized_fit = pd.DataFrame(data=pca_df, columns=['PC1', 'PC2'])

    df_normalized_fit.insert(0, "patient_id", patient_id_list)

    # Creates new patient id column in the data frame
    # Then merges the origional data frame on patient_id to add the time point and progression cycle columns
    meta_df = df[['patient_id', 'time_point', 'tumor_fraction','progression_cycle']].drop_duplicates()
    df_normalized_fit = df_normalized_fit.merge(meta_df, on="patient_id", how="left")

    # Conditions and values for conditions
    conditions = [
        (df_normalized_fit["tumor_fraction"] <= 0.03),
        (df_normalized_fit["tumor_fraction"] > 0.03) & (df_normalized_fit["tumor_fraction"] < 0.3),
        (df_normalized_fit["tumor_fraction"] >= 0.3)
    ]
    values = ['low', 'medium', 'high']

    df_normalized_fit["tumor_group"] = np.select(conditions, values, default='unknown')

    # Creatse scatter plot for normalized data 
    sns.set_theme()
    fig, ax = plt.subplots(figsize=(8, 6))
    
    title = ''
    if normalize_axis == 0:
        title = "Cental Depth PCA Plot (per patient) - C1 (Normalized Across Features)"
    else:
        title = "Cental Depth PCA Plot (per patient) - C1 (Normalized Across Samples)"

    _scatter_plot_(
        fig,
        ax,
        x=df_normalized_fit['PC1'],
        y=df_normalized_fit['PC2'],
        color=df_normalized_fit['tumor_group'],
        title=title,
        xlabel='PC1',
        ylabel='PC2',
        legend_title='Tumor Group'
    )
        
    # Saves plot to output
    output_file = os.path.join(output_path, output_title)
    fig.savefig(output_file)


def _save_df_as_tsv_(df, output_path, file_name):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Data containing information to save
        output_path: str
            Path to output directory
        file_path: str
            File name for tsv

    Function:
    ---------
        - Saves data frame as a tsv in output path
        
    Returns: 
    --------
        Void
    """
    # Saves data frame
    df_output_path = os.path.join(output_path, file_name)
    df.to_csv(df_output_path, sep='\t', index=False)


def _pivot_table_(df, values, index, columns, aggfunc, fill_value):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Data containing information to save
        values: str
            Selects column that will be the values of the pivot table
        index: str
            Selects column that will be the indexes/rows of the pivot table
        columns: str
            Selects column that will be the columns of the pivot table
        aggfunc: str
            Specifies how to aggregate the values
        fill_value: int
            Specifies how to fill if there is a NaN value

    Function:
    ---------
        - Pivots a data frame around a certain index and column
        
    Returns:
    --------
        pivoted data frame
    """
    df_pivot = pd.pivot_table(
        df,
        values=values,
        index=index,
        columns=columns,
        aggfunc=aggfunc,
        fill_value=fill_value)
    
    return df_pivot


def _normalize_(df, axis=0):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Data containing information to normalize
        axis: int
            default: 0
            Specifies what axis to normalize at: 0 or 1

    Function:
    ---------
        - If there is a patient_id column it strips and saves it
        - Then zscores the data frame
        
    Returns:
    --------
        normalized data frame, list of patient ids
    """
    # If patient_id strip and save column
    if "patient_id" in df.columns:
        patient_ids = df['patient_id']
        df_temp = df.drop(columns=['patient_id']).copy()
    else:
        patient_ids = None
        df_temp = df.copy()

    normalized_df = df_temp.copy()
    
    # Normalize across columns/tfbs
    if axis == 0:
        for row in df_temp.index:
            std = statistics.pstdev(df_temp.loc[row])
            mean = df_temp.loc[row].mean()
            normalized_df.loc[row] = (df_temp.loc[row] - mean) / std
    # Normalize across rows/patient_ids
    elif axis == 1:
        for column in df_temp.columns:
            std = statistics.pstdev(df_temp[column])
            mean = df_temp[column].mean()
            normalized_df[column] = (df_temp[column] - mean) / std

    return normalized_df, patient_ids


def _variance_(df, axis=0):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Data containing information to save
        axis: int
            Default: 0
            Specifies which axis to find the variance for: 0 is column, 1 is row

    Function:
    ---------
        - If there is a patient_id column it strips and saves it
        - Then finds the variance the data frame at the specified axis
        
    Returns:
    --------
        decending list of vairances, list of patient ids
    """
    if "patient_id" in df.columns:
        patient_ids = df['patient_id']
        df_temp = df.drop(columns=['patient_id'])
    else:
        patient_ids = None
        df_temp = df.copy()

    variance = df_temp.var(ddof=0, axis=axis).sort_values(ascending=False)

    return variance, patient_ids


def _PCA_(df):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Data containing information to save

    Function:
    ---------
        - Creates PCA model
        - Fits data to PCA
        - Transforms data
        
    Returns:
    --------
        2D array of principle components
    """
    # Selects amount of principle componenets
    pca = PCA(n_components=2)

    # Fit normalized data to PCA then applies the PCA transformation to the data
    pca.fit(df)
    pca_df = pca.transform(df)

    return pca_df


def _scatter_plot_(fig, ax, x, y, color=None, title=None, xlabel=None, ylabel=None, palette='tab10', legend_title=None, legend_loc='upper left', legend_outside=True):
    """
    Parameters:
    -----------
        listed parameters

    Function:
    ---------
        - If color is provided it assigns a color pallete for each unique color group and plots each subset with respective color
        - Else it does a normal scatter plot
        - Creates title, axis titles, etc.
        
    Returns: 
    --------
        scatter plot as fig
    """

    if color is not None:
        color = np.array(color)
        groups = np.unique(color)
        colors = sns.color_palette(palette, len(groups))
        color_map = dict(zip(groups, colors))

        for group in groups:
            mask = color == group
            ax.scatter(
                np.array(x)[mask],
                np.array(y)[mask],
                label=group,
                color=color_map[group]
            )
    else:
        ax.scatter(x, y)

    if title:
        ax.set_title(title, fontsize=14)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    if color is not None:
        if legend_outside:
            ax.legend(title=legend_title, bbox_to_anchor=(1.05, 1), loc=legend_loc)
        else:
            ax.legend(title=legend_title, loc=legend_loc)

    fig.tight_layout()
    return fig


def create_central_depth_heatmap(output_path, gene_list, normalization_level, tfbs_amount, progression_grouping_dict):
    """
    Parameters:
    -----------
        output_path: str
            Path to directory where output files will be saved
        gene_list: list
            List of genes to annotate or highlight in the heatmap
        normalization_level: str
            Normalization method to apply ("feature level", "sample level", or "None")
        tfbs_amount: int
            Number of top TFBS sites (by variance) to include; -1 uses all sites
        progression_grouping_dict: dict
            Dictionary mapping patient IDs to responder/non-responder groups

    Function:
    ---------
        - Reads combined data table from output_path
        - Selects only C1 patients
        - Creates heatmap-ready DataFrame with normalization and variance filtering
        - Saves intermediate heatmap data table
        - Generates clustered heatmap of central depth values

    Returns:
    --------
        None
    """

    # Creates name of combined data table
    combined_data_table_path = os.path.join(output_path, 'data-tables/combined_data_table.tsv')

    # Selects the data frame
    # Reads in data table into data frame
    df = pd.read_table(combined_data_table_path)

    # Selects for just C1 patients
    df = df[df["time_point"] == 1].reset_index(drop=True)

    print(f"df: \n{df}")

    df_filtered_patients_with_meta = _create_heatmap_data_frame_(df, normalization_level, tfbs_amount, progression_grouping_dict)

    print(f"df_filtered_patients_with_meta: \n{df_filtered_patients_with_meta}")

    _save_df_as_tsv_(df_filtered_patients_with_meta, output_path, "data-tables/heatmap_data_table.tsv")

    df_filtered_patients_no_meta = heatmap_data_frame_strip_meta_data(df_filtered_patients_with_meta)

    print(f"df_filtered_patients_no_meta: {df_filtered_patients_no_meta}")

    # Shape title of plot and name of file based on tfbs_amount
    if tfbs_amount > -1:
        tfbs_amount_title = tfbs_amount
        number_of_sites_title = f"Top {tfbs_amount} Sites by Variance"

    else:
        tfbs_amount_title = "all"
        number_of_sites_title = f"All Sites"

    # Shape title of plot and name of file based on  normalization_level
    normalization_level_title = f"Normalized at {normalization_level}"
        
    plot_title = f"Central Depth Clustered Heatmap - {number_of_sites_title} - {normalization_level_title}"
    file_name = f"normalized_{normalization_level.replace(" ", "_")}_{tfbs_amount_title}_tfbs_by_variance_heatmap.png"

    _heatmap_(
        df_filtered_patients_no_meta, 
        output_path, 
        plot_title, 
        file_name, 
        gene_list
    )


def _create_heatmap_data_frame_(df, normalization_level, tfbs_amount, progression_grouping_dict):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Input patient/site central depth data
        normalization_level: str
            Normalization method to apply ("feature level", "sample level", or "None")
        tfbs_amount: int
            Number of top TFBS sites (by variance) to include; -1 uses all sites
        progression_grouping_dict: dict
            Dictionary mapping patient IDs to responder/non-responder groups

    Function:
    ---------
        - Pivots the input DataFrame so rows = patients and columns = sites
        - Resets index after pivoting
        - Applies normalization and selects high-variance TFBS sites
        - Merges in patient metadata (time point, progression cycle, tumor fraction, genomic instability)
        - Filters patients by progression group (responder vs nonresponder)
        - Adds progression group and tumor group classification columns
        - Reorders columns so that metadata appears directly after patient_id

    Returns:
    --------
        df_filtered_patients: pandas DataFrame
            Heatmap-ready DataFrame with normalized central depth values, metadata, and progression/tumor group annotations
    """

    # Pivots table around sites so that the patient_id are the rows and sites are the columns
    df_pivot = _pivot_table_(
        df, 
        values="central_depth",
        index="patient_id",
        columns="site",
        aggfunc='first',
        fill_value=0
    )

    # Resetting the index counters on the side of the pivot table
    df_pivot.reset_index(inplace=True)

    # Normalizes based on normalization_level, and filteres out data frame based on tfbs_amount (if "None" and -1 then nothing happens to dataframe)
    df_normalized, high_variance, _ = normalization_and_variance_selection(df_pivot, normalization_level, tfbs_amount)
    df_normalized_filtered = df_normalized.loc[:, df_normalized.columns.intersection(high_variance.index)]

    # Adding in metadata for each patient
    df_meta = df[[
        "patient_id",
        "time_point",
        "progression_cycle",
        "tumor_fraction",
        "genomic_instability"
    ]].drop_duplicates(subset="patient_id")
    # Merge metadata into pivot table
    df_normalized_filtered = df_normalized_filtered.merge(df_meta, on="patient_id", how="left")

    # Prep df by filtering by responder and nonresponder and adding 
    df_filtered_patients = filter_df_by_progression_group(df_normalized_filtered, progression_grouping_dict)
    df_filtered_patients = add_progression_group_column(df_filtered_patients, progression_grouping_dict)
    df_filtered_patients = add_tumor_group_column(df_filtered_patients)

    # Reorder so metadata is right after patient_id
    columns_list = df_filtered_patients.columns.tolist()
    columns_reordered = columns_list[:1] + columns_list[-6:] + columns_list[1:-6]
    df_filtered_patients = df_filtered_patients[columns_reordered]

    return df_filtered_patients


def heatmap_data_frame_strip_meta_data(df):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Input patient/site central depth data

    Function:
    ---------
        - Strips df of meta data

    Returns:
    --------
        df:
            pandas DataFrame with no meta data 
    """

    df = df.drop(columns=["time_point", "progression_cycle", "tumor_fraction", "genomic_instability"])
    return df


def _heatmap_(df, output_path, plot_title, file_name, gene_list):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Input DataFrame containing patient gene/site expression values 
            along with progression_group and tumor_group annotations
        output_path: str
            Directory path where the heatmap image will be saved
        plot_title: str
            Title text to display above the heatmap
        file_name: str
            File name for saving the heatmap (e.g. "heatmap.png")
        gene_list: list or None
            List of gene/site names to include in the heatmap;
            if None, all available columns (except metadata) are used

    Function:
    ---------
        - Extracts progression group and tumor group metadata per patient
        - Creates unique color palettes and lookup tables for each group
        - Builds a color annotation bar (`col_colors`) for the heatmap
        - Optionally filters the DataFrame to include only a subset of genes/sites
        - Generates a hierarchical clustermap with seaborn:
            * Rows = genes/sites
            * Columns = patients
            * Colors represent z-scored expression values
        - Adds custom legends for progression groups and tumor groups
        - Improves aesthetics (thicker dendrogram lines, tick labels, colorbar title)
        - Adjusts y-axis tick marks if fewer than 60 genes/sites
        - Exports the final heatmap figure to the specified output path

    Returns:
    --------
        None
            Saves the heatmap plot to file; does not return a DataFrame or figure
    """

    df_progression_cycle = df[["patient_id", "progression_group"]]
    df_progression_cycle = df_progression_cycle.set_index("patient_id")
    df_tumor_fraction = df[["patient_id", "tumor_group"]]
    df_tumor_fraction = df_tumor_fraction.set_index("patient_id")

    # Create color mapping for progression cycles
    # Grabs only unique progression cycle values
    # Creates a list of unique colors per cycle
    # Creates a look up table that assigns each cycle with a color
    unique_cycles = sorted(df_progression_cycle["progression_group"].dropna().unique())
    cycle_palette = [(210/255, 16/255, 72/255), (244/255, 179/255, 1/255)]#sns.color_palette("colorblind", len(unique_cycles))
    cycles_lookup_table = dict(zip(unique_cycles, cycle_palette))
    # row_colors = df_progression_cycle["progression_group"].map(cycles_lookup_table)

    unique_tumor_group = sorted(df_tumor_fraction["tumor_group"].dropna().unique())
    tumor_palette = sns.color_palette("colorblind", len(unique_tumor_group))
    tumor_group_lookup_table = dict(zip(unique_tumor_group, tumor_palette))
    # row_colors = df_tumor_fraction["tumor_group"].map(tumor_group_lookup_table)

    # Assign the colors to the patient id of progression cycle and tumor group data frames
    col_colors = pd.DataFrame({
        "Responder Group": df_progression_cycle["progression_group"].map(cycles_lookup_table),
        "Tumor Group": df_tumor_fraction["tumor_group"].map(tumor_group_lookup_table)
    })

    if gene_list:
        df = df.set_index("patient_id")
        df = df[df.columns.intersection(gene_list)]
        print(f"df_after_gene_list_filter: \n{df}")
    else:
        print("nothing")
        df = df.set_index("patient_id")
        df = df.drop(columns=["progression_group", "tumor_group"])

    # Creates clustermap with single bar above to show progression cycle
    heatmap = sns.clustermap(
        df.T, 
        col_colors=col_colors, 
        cmap="coolwarm", 
        linecolor='gray', 
        figsize=(30,20),
    )
    
    # Colors in the specific progression cycle bar
    for label in unique_cycles:
        heatmap.ax_col_dendrogram.bar(0, 0, color=cycles_lookup_table[label],
                                label=f'{label}', linewidth=0)
    
    # Creates and places progression cycle legend
    heatmap.ax_col_dendrogram.legend(loc="center", ncol=6)

    # Colors in the specific tumor group bar
    for label in unique_tumor_group:
        heatmap.ax_col_dendrogram.bar(0, 0, color=tumor_group_lookup_table[label],
                                label=f'{label.title()} TFx', linewidth=0)
    
    # Creates and places tumor group legend
    heatmap.ax_col_dendrogram.legend(loc="center", ncol=6, bbox_to_anchor=(.5,1.1), fontsize=15)

    # Add a title with padding
    fig = plt.gcf()

    fig.suptitle(
            plot_title,
            fontsize=16,
            y=0.99
        )
   
    # Trying to update the look of the heatmap

    # Make dendrogram lines thicker
    for ax in [heatmap.ax_row_dendrogram, heatmap.ax_col_dendrogram]:
        for line in ax.collections:
            line.set_linewidth(2.5)

    heatmap.ax_heatmap.tick_params(axis='y', labelsize=15)

    heatmap.cax.set_title('Zscore - Expression (Low to high)', fontsize=16)

    plt.rc('legend', fontsize=20)

    # [left, bottom, width, height]
    heatmap.cax.set_position([.205, 1, -.6, 0])

    # Only reset tick marks if there are less than 100 tfbs
    if len(df.columns) < 60:
        heatmap.ax_heatmap.set_yticks(
            [x + 0.5 for x in range(len(df.columns))]  # adjust offset
        )

        heatmap.ax_heatmap.set_yticklabels(
            df.columns, rotation=0, fontsize=15, ha="left", va="center"
        )

    # End Trying to update the look of the heatmap

    # (left, bottom, right, top)
    fig.tight_layout(rect=[0, 0, .95, 0.95])

    # Creates output file for plot
    output_file = os.path.join(output_path, file_name)
    fig.savefig(output_file)
    plt.close()


def normalization_and_variance_selection(df, normalization_level, tfbs_amount):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Input DataFrame containing patient features (TFBS values)
        normalization_level: str
            Normalization strategy; one of:
                - "feature level" → normalize across features
                - "sample level" → normalize across samples
                - "None" → no normalization
        tfbs_amount: int
            Number of top TFBS to select based on variance.
            If -1, use all TFBS.

    Function:
    ---------
        - Computes variance of features (TFBS) across patients
        - Selects either all TFBS or the top `tfbs_amount` TFBS by variance
        - Applies the requested normalization strategy
        - Reattaches patient IDs to the normalized DataFrame

    Returns:
    --------
        df_normalized: pandas DataFrame
            Normalized feature matrix with patient IDs as index
        variance_series: pandas Series
            Variance values of selected TFBS, sorted in descending order
        patient_id_list: list
            List of patient IDs from the input DataFrame
    """

    # Takes all TFBS
    if tfbs_amount == -1:
        variance_series, patient_id_list = _variance_(df, 0)
    # Takes top tfbs
    else:
        variance_series, patient_id_list = _variance_(df, 0)
        variance_series = variance_series.head(tfbs_amount)

    if normalization_level == "feature level":
        # Standard normailization across features -> returns a data frame
        df_normalized, patient_id_list = _normalize_(df, axis=0)
    elif normalization_level == "sample level":
        # Standard normalization across samples -> returns a data frame
        df_normalized, patient_id_list = _normalize_(df, 1)
    elif normalization_level == "None":
        df_normalized = df
    else:
        print("Provide a normalization level: \"feature level\", \"sample level\", \"None\" ")
        sys.exit(1)

    # Renaming and adding back the patient_id column
    df_normalized.index.name = 'patient_id'
    df_normalized.index = patient_id_list

    return df_normalized, variance_series, patient_id_list


def filter_df_by_progression_group(df, progression_grouping_dict):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Input DataFrame containing patient metadata including "progression_cycle"
        progression_grouping_dict: dict
            Mapping of progression cycle values to group names.
            Values with "None" are excluded.

    Function:
    ---------
        - Identifies valid progression cycles from the grouping dictionary
        - Filters the DataFrame to include only patients belonging to those cycles

    Returns:
    --------
        pandas DataFrame
            Filtered DataFrame containing only patients in selected progression cycles
    """

    key_list = [k for k, v in progression_grouping_dict.items() if v != "None"]
    return df[df["progression_cycle"].isin(key_list)].copy()


def add_progression_group_column(df, progression_grouping_dict):
    """
    Parameters:
    -----------
        df: pandas DataFrame
            Input DataFrame containing a "progression_cycle" column
        progression_grouping_dict: dict
            Mapping from cycle index (1-6) to group names (e.g. responder/non-responder)

    Function:
    ---------
        - Builds a mapping from numerical cycle values to group names
        - Creates a new column "progression_group" in the DataFrame

    Returns:
    --------
        pandas DataFrame
            DataFrame with an additional "progression_group" column
    """

    value_list = list(progression_grouping_dict.values())
    # Map to responders/non-responders
    progression_map = {
        1.0: value_list[0],
        2.0: value_list[1],
        3.0: value_list[2],
        4.0: value_list[3],
        5.0: value_list[4],
        6.0: value_list[5],
    }

    df["progression_group"] = df["progression_cycle"].map(progression_map)

    return df


def add_tumor_group_column(df):
    """
    add_tumor_group_column

    Parameters:
    -----------
        df: pandas DataFrame
            Input DataFrame containing a "tumor_fraction" column

    Function:
    ---------
        - Categorizes patients into tumor burden groups based on thresholds:
            * <= 0.03 → "low"
            * 0.03 - 0.3 → "medium"
            * >= 0.3 → "high"
        - Adds a new column "tumor_group" to the DataFrame

    Returns:
    --------
        pandas DataFrame
            DataFrame with an additional "tumor_group" column
    """

    # Conditions and values for conditions
    conditions = [
        (df["tumor_fraction"] <= 0.03),
        (df["tumor_fraction"] > 0.03) & (df["tumor_fraction"] < 0.3),
        (df["tumor_fraction"] >= 0.3)
    ]
    values = ['low', 'medium', 'high']

    df["tumor_group"] = np.select(conditions, values, default='unknown')

    return df


def group_df_by_progression_cycle(df_filtered_patients, df_normalized_filtered):
    """
    Parameters:
    -----------
        df_filtered_patients: pandas DataFrame
            Patient metadata DataFrame containing "patient_id" and "progression_group"
        df_normalized_filtered: pandas DataFrame
            Normalized feature matrix with patient IDs as index

    Function:
    ---------
        - Restricts metadata to patients present in normalized data
        - Removes duplicate patient IDs if present
        - Returns DataFrame with mapping of patient_id to progression cycle

    Returns:
    --------
        pandas DataFrame
            Index = patient_id, single column "progression_cycle"
    """

    # Now restrict to patients in normalized filtered index
    df_progression_cycle = (
        df_filtered_patients[df_filtered_patients["patient_id"].isin(df_normalized_filtered.index)]
        .drop_duplicates(subset="patient_id")
        [["patient_id", "progression_group"]]
        .set_index("patient_id")
        .rename(columns={"progression_group": "progression_cycle"})
    )

    return df_progression_cycle


def group_df_by_tumor_fraction(df_filtered_patients, df_normalized_filtered):
    """
    Parameters:
    -----------
        df_filtered_patients: pandas DataFrame
            Patient metadata DataFrame containing "patient_id" and "tumor_fraction"
        df_normalized_filtered: pandas DataFrame
            Normalized feature matrix with patient IDs as index

    Function:
    ---------
        - Restricts metadata to patients present in normalized data
        - Removes duplicate patient IDs if present
        - Returns DataFrame with mapping of patient_id to tumor fraction

    Returns:
    --------
        pandas DataFrame
            Index = patient_id, single column "tumor_fraction"
    """

    # Creating a data frame with tumor fraction and patient id based on the filtered data frame
    df_tumor_fraction = (
        df_filtered_patients[df_filtered_patients["patient_id"].isin(df_normalized_filtered.index)]
        .drop_duplicates(subset="patient_id")
        [["patient_id", "tumor_fraction"]]
        .set_index("patient_id")
    )

    return df_tumor_fraction

# central_depth_differential_activity function
# Take in output path
# For each TFBS it does a man whitney u test on responders vs non responders
# Then it does a fdr correction on the calculated p-values and prints out the significant TFBS
def central_depth_differential_activity(output_path, pvalue_cutoff, normalization_level, tfbs_amount, progression_grouping_dict, get_dif_exp_genes=False):
    # Save correct heatmap_data_table.tsv

    # Creates name of combined data table
    combined_data_table_path = os.path.join(output_path, 'data-tables/combined_data_table.tsv')

    # Selects the data frame
    # Reads in data table into data frame
    df = pd.read_table(combined_data_table_path)

    # Selects for just C1 patients
    df = df[df["time_point"] == 1].reset_index(drop=True)

    print(f"df: \n{df}")

    heatmap_data_table = _create_heatmap_data_frame_(df, "None", tfbs_amount, progression_grouping_dict)

    print(f"heatmap_data_table: \n{heatmap_data_table}")

    _differential_activity_(output_path, heatmap_data_table, pvalue_cutoff, normalization_level, tfbs_amount, get_dif_exp_genes=False)

def _differential_activity_(output_path, df, pvalue_cutoff, normalization_level, tfbs_amount, get_dif_exp_genes=False):
    differential_activity_dict = {}
    responder_list = []
    non_responder_list = []

    print(f"df_dif_act: \n{df}")

    # For every TFBS in the heatmap df calculate the pvalue between the expression values of responders vs nonresponders
    # Using a Mann Whiteny U test, then add it to a dictionary, then using a fdr correction correct all pvalues
    # Print out signficant TFBS
    for column in df.columns:
        if column == "patient_id" or column == "time_point" or column == "progression_cycle" or column == "tumor_fraction" or column == "genomic_instability" or column == "progression_group" or column == "tumor_group":
            continue
        responder_list = []
        non_responder_list = []
        for i, row in df.iterrows():
            value = row[column]
            if row["progression_group"] == "responder":
                responder_list.append(value)
            elif row["progression_group"] == "non-responder":
                non_responder_list.append(value)
            else:
                print("didnt work")

        # Do Man Whitney U Test on lists, do a fdr correction, and add it to the differential activity dictionary 
        mann_statistic, mann_pvalue = stats.mannwhitneyu(responder_list, non_responder_list)
        differential_activity_dict[column] = {"pvalue": mann_pvalue, "responder_median": statistics.median(responder_list), "non_responder_median": statistics.median(non_responder_list)}

    # Need to make map into list
    dictionary_to_list = list(differential_activity_dict.values())

    pvalues_list = []
    for i, x in enumerate(dictionary_to_list):
        pvalues_list.append(x["pvalue"])

    _, adj_differential_activity_pvalues, _, _ = multipletests(pvalues_list, method='fdr_bh')

    # Reassigns dictionary values to corrected pvalues
    keys = list(differential_activity_dict.keys())
    for i, value in enumerate(adj_differential_activity_pvalues):
        if i < len(keys):
            differential_activity_dict[keys[i]]["pvalue"] = float(value)    

    differential_activity_dict_significant = {key: value for key, value in differential_activity_dict.items() if value["pvalue"] <= pvalue_cutoff}

    # Prints out significant tfbs and adds them to a gene list
    dif_exp_gene_list = []
    print("Genes that are expressed at significantly different levels:")
    for i, (key, value) in enumerate(differential_activity_dict_significant.items()):
        significance_side = "responder" if value["responder_median"] > value["non_responder_median"] else "non-responder"
        print(f"{i} {key} {value["pvalue"]} -> expression is higher in {significance_side}")
        dif_exp_gene_list.append(key)

    if get_dif_exp_genes:
        return dif_exp_gene_list

    # Calls the heatmap creator on just the gene list
    df_filtered_patients_no_meta = heatmap_data_frame_strip_meta_data(df)

    print(f"df_filtered_patients_no_meta: \n{df_filtered_patients_no_meta}")

    # Shape title of plot and name of file based on tfbs_amount
    if tfbs_amount > -1:
        tfbs_amount_title = tfbs_amount
        number_of_sites_title = f"Significantly Differentially Expressed Sites from Top {tfbs_amount} by Variance"

    else:
        tfbs_amount_title = "all"
        number_of_sites_title = f"All Significantly Differentially Expressed Sites"

    # Shape title of plot and name of file based on  normalization_level
    normalization_level_title = f"Normalized at {normalization_level}"
        
    plot_title = f"Central Depth Clustered Heatmap - {number_of_sites_title} - {normalization_level_title}"
    file_name = f"sig_diff_exp_sites_normalized_{normalization_level.replace(" ", "_")}_from_{tfbs_amount_title}_by_variance_heatmap.png"

    print(f"dif_exp_gene_list: {dif_exp_gene_list}")

    _heatmap_(
        df_filtered_patients_no_meta, 
        output_path, 
        plot_title, 
        file_name, 
        dif_exp_gene_list
    )
    

def central_depth_vs_tumor_fraction_density_plot(output_path):
    """
    Parameters:
    -----------
        output_path: str
            Directory path where the heatmap image will be saved

    Function:
    ---------
        - Creates a density plot from combinding FGA_data_table and central_depth_data_table
        - Collects a set of prostate specific oncogenes and appies them as a filter to the data frame
        - Plots the oncogenes and only those above 3% TFx on a density plot 

    Returns:
    --------
        None
            Saves plt in output path
    """

    # Selects data table
    FGA_data_file = os.path.join(output_path, "data-tables/FGA_data_table.tsv")
    CD_data_file = os.path.join(output_path, "data-tables/central_depth_data_table.tsv")

    FGA_data_table = pd.read_table(FGA_data_file)
    CD_data_table = pd.read_table(CD_data_file)

    CD_data_table_with_TFx = pd.merge(
        CD_data_table,
        FGA_data_table[["patient_id", "time_point", "tumor_fraction"]],
        on=["patient_id", "time_point"],
        how="left"
    )

    # Find prostate specific oncogenes
    oncogenes = pd.read_table("/fh/fast/ha_g/projects/ProstateTAN/analysis_CN-SV/JASMINE/Plots/data/20230823_Composite_PC_GeneList.txt", header=None)

    # Takes only the names of the oncogenes
    oncogene_set = set(oncogenes[0])

    # Filteres out the table for the genes that are only in the gene set
    CD_data_table_with_TFx = CD_data_table_with_TFx[CD_data_table_with_TFx["site"].isin(oncogene_set)]

    x_axis = CD_data_table_with_TFx["tumor_fraction"]
    y_axis = CD_data_table_with_TFx["central_depth"]

    # Removes all tumor fraction less then 3% and removes any central depth that is 0
    nonzero_mask = (x_axis > 0.03) & (y_axis > 0)
    x_axis = x_axis[nonzero_mask]
    y_axis = y_axis[nonzero_mask]

    # Creates the scatter plot with tumor% as x-axis and genome instability as y-axis
    fig, ax = plt.subplots(subplot_kw={'projection': 'scatter_density'})
    scatter = ax.scatter_density(x_axis, y_axis, cmap='viridis')
    fig.colorbar(scatter, label='Density')
    ax.set_xlabel("Tumor Fraction (%)") 
    ax.set_ylabel("Central Depth")
    ax.set_title("Central Depth vs Tumor Fraction")
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
    output_plot = os.path.join(output_path, "tfx_vs_CD_scatter_plot.png")
    plt.savefig(output_plot, dpi=200) 
    plt.close() 
    

def _regression_(x_data, y_data):
    """
    Parameters:
    -----------
        x_data : np.ndarray
            Tumor fraction values, shaped as (-1, 1).
        y_data : np.ndarray
            Central depth values, shaped as (-1, 1).

    Function:
    ---------
        - Fits a linear regression model to predict central depth from tumor fraction.
        - Calculates residuals as the difference between predicted and actual values.

    Returns:
    --------
        residuals : np.ndarray
            Residual values from the regression fit.
        model : sklearn.linear_model.LinearRegression
            Trained linear regression model.
    """

     # Linear Regression
    model = LinearRegression()

    # Fitting a line to the data
    model.fit(x_data, y_data)
    # Using that line to predict central depth from tumor fraction
    predicted = model.predict(x_data)
    # Calculate residuals
    residuals = predicted - y_data

    return residuals, model


def central_depth_linear_regression(output_path, pvalue_cutoff, normalization_level, tfbs_amount, progression_grouping_dict):
    """
    Parameters:
    -----------
        output_path : str
            Directory path where outputs will be saved.
        pvalue_cutoff : float
            P-value cutoff for differential activity analysis.
        normalization_level : str
            Normalization method or level applied to data.
        tfbs_amount : int
            Number of TFBSs (transcription factor binding sites) to consider (-1 for all).
        progression_grouping_dict : dict
            Dictionary mapping progression groups for clustering and heatmap creation.

    Function:
    ---------
        - Reads combined data table and filters for time point 1.
        - Fits a linear regression model of central depth vs. tumor fraction.
        - Replaces central depth values with residuals for downstream analysis.
        - Identifies significantly expressed TFBSs via differential activity.
        - Plots:
            * Tumor fraction vs. central depth (with regression fit).
            * Residual plots of significant TFBSs.
            * Clustered heatmap of residuals.
        - Saves outputs (plots, heatmap) in the specified output path.

    Returns:
    --------
        None
            Outputs are saved to `output_path`.
    """

    # Selects data table
    combined_data_file = os.path.join(output_path, "data-tables/combined_data_table.tsv")

    combined_data_table = pd.read_table(combined_data_file)

    print(f"combined_data_table: \n{combined_data_table}")

    # Drop all rows that are not time point 1
    combined_data_table = combined_data_table[combined_data_table["time_point"] == 1].reset_index(drop=True)

    print(f"combined_data_table: \n{combined_data_table}")

    tfx_data = np.array(combined_data_table["tumor_fraction"]).reshape((-1,1))
    cd_data = np.array(combined_data_table["central_depth"]).reshape((-1,1))

    residuals, model = _regression_(tfx_data, cd_data)

    print(f"residuals: \n{len(residuals)}")
    
    # Replace central depth column with residuals
    combined_data_table["central_depth"] = residuals

    # Save correct heatmap_data_table.tsv
    print(f"combined_data_table: \n{combined_data_table}")

    df_filtered_patients_with_meta = _create_heatmap_data_frame_(combined_data_table, normalization_level, tfbs_amount, progression_grouping_dict)

    # Get most significatly expressed genes between responders and non-responders
    gene_list = _differential_activity_(output_path, df_filtered_patients_with_meta, pvalue_cutoff, normalization_level, tfbs_amount, True)

    print(f"gene list: \n {gene_list}")
 
    # Creates the scatter plot with tfx as x-axis and central depth as y-axis for specific genes
    filtered_CD_data_table_with_TFx = combined_data_table[combined_data_table["site"].isin(gene_list)]
    tfx_data_specific = np.array(filtered_CD_data_table_with_TFx["tumor_fraction"]).reshape((-1,1))
    cd_data_specific = np.array(filtered_CD_data_table_with_TFx["central_depth"]).reshape((-1,1))

    # Fitting a line to the data
    model.fit(tfx_data_specific, cd_data_specific)
    # Using that line to predict central depth from tumor fraction
    specific_predicted = model.predict(tfx_data_specific)
    # Creates a 
    res = stats.linregress(tfx_data_specific.flatten(), cd_data_specific.flatten())
    print(f"R-squared: {res.rvalue**2:.6f}")
    fig, ax = plt.subplots()
    scatter = ax.scatter(tfx_data_specific, cd_data_specific)
    # Plot the line
    ax.plot(tfx_data_specific, specific_predicted, color="red", linewidth=2, label="Linear Fit")
    ax.text(
        0.05, 0.95, f"$R^2$ = {res.rvalue**2:.6f}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7)
    )
    # Optional: show legend
    ax.legend()
    ax.set_xlabel("Tumor Fraction (%)") 
    ax.set_ylabel("Central Depth")
    ax.set_title("Central Depth vs Tumor Fraction for Specific Genes")
    ax.grid(True)

    # Save plot
    output_plot = os.path.join(output_path, "tfx_vs_CD_scatter_plot_specific_genes.png")
    plt.savefig(output_plot, dpi=200) 
    plt.close() 

    # Call create_central_depth_heatmap on residual table - running it twice to get proper plot formatting
    df_filtered_patients_no_meta = heatmap_data_frame_strip_meta_data(df_filtered_patients_with_meta)

    print(f"df_filtered_patients_no_meta: {df_filtered_patients_no_meta}")

    # Shape title of plot and name of file based on tfbs_amount
    if tfbs_amount > -1:
        tfbs_amount_title = tfbs_amount
        number_of_sites_title = f"Top {tfbs_amount} Significantly Differentially Expressed Residuals by Variance"

    else:
        tfbs_amount_title = "all"
        number_of_sites_title = f"All Significantly Differentially Expressed Residuals at CI {pvalue_cutoff}"

    # Shape title of plot and name of file based on  normalization_level
    normalization_level_title = f"Normalized at {normalization_level}"
        
    plot_title = f"Central Depth Clustered Heatmap - {number_of_sites_title} - {normalization_level_title}"
    file_name = f"diff_sig_exp_residuals_normalized_{normalization_level.replace(" ", "_")}_{tfbs_amount_title}_tfbs_CI_{pvalue_cutoff}_heatmap.png"

    _heatmap_(
        df_filtered_patients_no_meta, 
        output_path, 
        plot_title, 
        file_name, 
        gene_list
    )

    # Plots residules of the significant differencially expressed tfbs
    sns.residplot(
        x=tfx_data_specific, 
        y=cd_data_specific, 
        lowess=True, 
        line_kws=dict(color="r")
    )

    plt.xlabel("Tumor Fraction")
    plt.ylabel("Residuals")
    plt.title("Residual Plot for Tumor Fraction vs. Central Depth")
    plt.tight_layout()

    # Creates output file for plot
    output_file = os.path.join(output_path, f"tfx_central_depth_residual_plot.png")
    plt.savefig(output_file)
    plt.close()


# spearman_rho_of_central_depth_and_tumor_fraction() function
def central_depth_and_tumor_fraction_spearman_rho(output_path, pvalue_cutoff, normalization_level, tfbs_amount, progression_grouping_dict):
    """
    Parameters:
    -----------
        output_path : str
            Directory path where outputs will be saved.
        pvalue_cutoff : float
            Minimum absolute Spearman rho cutoff for significance.
        normalization_level : str
            Normalization method or level applied to data.
        tfbs_amount : int
            Number of TFBSs (transcription factor binding sites) to consider (-1 for all).
        progression_grouping_dict : dict
            Dictionary mapping progression groups for clustering and heatmap creation.

    Function:
    ---------
        - Reads combined data table and creates a heatmap-ready dataframe.
        - Computes Spearman correlation (rho) between tumor fraction and each TFBS.
        - Saves scatter plots of significantly correlated TFBSs.
        - Plots histogram of Spearman rho values across all TFBSs.
        - Identifies prostate-specific oncogenes with significant positive/negative correlations.

    Returns:
    --------
        None
            Outputs (scatter plots, histogram, printed oncogene correlations) are saved to `output_path`.
    """

    combined_data_file = os.path.join(output_path, "data-tables/combined_data_table.tsv")

    combined_data_table = pd.read_table(combined_data_file)

    heatmap_data_table = _create_heatmap_data_frame_(combined_data_table, normalization_level, tfbs_amount, progression_grouping_dict)

    # Making directores for the scatter plots
    try:
        os.mkdir(f"{output_path}/rho-stat-high")
        os.mkdir(f"{output_path}/rho-stat-low")
    except FileExistsError:
        print("directories could not be made")

    spearman_rho_dict = {}
    tfx_column = np.array(heatmap_data_table["tumor_fraction"])

    for column in heatmap_data_table.columns:
        if column == "patient_id" or column == "time_point" or column == "progression_cycle" or column == "tumor_fraction" or column == "genomic_instability" or column == "tumor_group":
            continue
        tfbs_column = np.array(heatmap_data_table[column])
        res = stats.spearmanr(tfx_column, tfbs_column, axis=0, nan_policy='raise', alternative='two-sided')
        spearman_rho_dict[column] = {"stat": res.statistic, "pvalue": res.pvalue}

        if res.pvalue <= 0.05 and abs(res.statistic) >= pvalue_cutoff:
            print(f"{column}: {res}")
            plt.scatter(tfx_column, heatmap_data_table[column])
            plt.xlabel("TFx")
            plt.ylabel("Central Depth")
            plt.title(column)
            if res.statistic >= 0.5:
                scatter_plot_file = os.path.join(output_path, f"rho-stat-high/temp_scatter_{column}.png")
            else:
                scatter_plot_file = os.path.join(output_path, f"rho-stat-low/temp_scatter_{column}.png")
            plt.savefig(scatter_plot_file)
            plt.close()  

    pvalues = [values["stat"] for values in spearman_rho_dict.values()]
    plt.hist(pvalues)
    plt.xlabel("spearman rho bins")
    plt.ylabel("Amount of TFBSs")
    plt.title("Histogram of spearman rho values")
    scatter_plot_file = os.path.join(output_path, f"temp_histogram.png")
    plt.savefig(scatter_plot_file)
    plt.close()  

    # Find prostate specific oncogenes
    oncogenes = pd.read_table("/fh/fast/ha_g/projects/ProstateTAN/analysis_CN-SV/JASMINE/Plots/data/20230823_Composite_PC_GeneList.txt", header=None)
    # Takes only the names of the oncogenes
    oncogene_set = set(oncogenes[0])

    print("positive correlation")
    for key, item in spearman_rho_dict.items():
        if item["pvalue"] <= 0.05 and item["stat"] >= 0.5 and key in oncogene_set:
            print(f"{key}")
    print()
    print("negative correlation")
    for key, item in spearman_rho_dict.items():
        if item["pvalue"] <= 0.05 and item["stat"] <= -0.5 and key in oncogene_set:
            print(f"{key}")
    
# _cox_forest_plot_ function
def _cox_forest_plot_():
    print("hello")
    return

def tfbs_list_vs_clinical_data_cox_forest_plots():
    # Selects data table
    combined_data_file = os.path.join(output_path, "data-tables/combined_data_table.tsv")

    combined_data_table = pd.read_table(combined_data_file)

    print(f"combined_data_table: \n{combined_data_table}")

    # Drop all rows that are not time point 1
    combined_data_table = combined_data_table[combined_data_table["time_point"] == 1].reset_index(drop=True)

    print(f"combined_data_table: \n{combined_data_table}")

    tfx_data = np.array(combined_data_table["tumor_fraction"]).reshape((-1,1))
    cd_data = np.array(combined_data_table["central_depth"]).reshape((-1,1))

    residuals, model = _regression_(tfx_data, cd_data)

    print(f"residuals: \n{len(residuals)}")
    
    # Replace central depth column with residuals
    combined_data_table["central_depth"] = residuals

    # Call create_central_depth_heatmap on residual table - running it twice to get proper plot formatting
    df_filtered_patients_no_meta = heatmap_data_frame_strip_meta_data(df_filtered_patients_with_meta)

    print(f"df_filtered_patients_no_meta: {df_filtered_patients_no_meta}")

    _calculate_all_hr_data()
    return


def _calculate_all_hr_data(df, duration_col, event_col, variables_for_hr_calculation, per_doubling=False):
   
    cph = CoxPHFitter()
    all_hr_data = []
    
    # Choose transformation based on per_doubling flag
    transform_type = 'log2p1' if per_doubling else 'zscore'

    for var in variables_for_hr_calculation:
        subset_cols = [var, duration_col, event_col]
        df_subset = df[subset_cols].dropna().copy()

        if df_subset.shape[0] < 10 or df_subset[var].nunique() < 2:
            print(f"Skipping HR calculation for {var} due to insufficient data (<10 obs) or variance (<2 unique values).")
            continue

        if transform_type == 'log2p1':
            if (df_subset[var] < -1).any():
                print(f"Warning: Variable '{var}' contains values < -1. Log2(x+1) transform might produce NaNs or complex numbers. Affected values will be NaNs.")
            transformed_var = np.log2(df_subset[var] + 1)
            transformed_var.replace([np.inf, -np.inf], np.nan, inplace=True)
            df_subset[var] = transformed_var
            df_subset.dropna(subset=[var], inplace=True)

            if df_subset.shape[0] < 10 or df_subset[var].nunique() < 2:
                print(f"Skipping HR for {var} after log2p1 transform due to insufficient data or variance.")
                continue
        elif transform_type == 'zscore':
            mean_val = df_subset[var].mean()
            std_val = df_subset[var].std()
            if std_val == 0 or pd.isna(std_val):
                print(f"Skipping HR calculation for {var} due to zero or NaN standard deviation (cannot z-score).")
                continue
            df_subset[var] = (df_subset[var] - mean_val) / std_val
        else:
            raise ValueError(f"Unknown transform_type: {transform_type}")

        cph.fit(df_subset, duration_col=duration_col, event_col=event_col, formula=f"`{var}`")
        summary = cph.summary
        hr = summary.loc[var, 'exp(coef)']
        ci_lower = summary.loc[var, 'exp(coef) lower 95%']
        ci_upper = summary.loc[var, 'exp(coef) upper 95%']
        p_val = summary.loc[var, 'p']
        all_hr_data.append({'Variable': var, 'HR': hr, 'CI_Lower': ci_lower, 'CI_Upper': ci_upper, 'p': p_val})

    if not all_hr_data:
        return pd.DataFrame()
    return pd.DataFrame(all_hr_data)