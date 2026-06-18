# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 1/27/26
# Purpose: Calculates fraction genome (FGA) altered per chromosome and outputs a tsv listing patients in decending order by FGA per chromosome.

from python_scripts.imports import *

def normalize_ranges(range_val):
    """Convert input into a clean list of range strings."""
    if pd.isna(range_val):
        return []
    
    # Case 1: already a list
    if isinstance(range_val, (list, tuple)):
        # handle nested stringified list inside list
        if len(range_val) == 1 and isinstance(range_val[0], str) and range_val[0].startswith('['):
            try:
                return ast.literal_eval(range_val[0])
            except Exception:
                return range_val
        return range_val
    
    # Case 2: string that might be a list
    if isinstance(range_val, str):
        range_val = range_val.strip()
        
        if range_val.startswith('[') and range_val.endswith(']'):
            try:
                return ast.literal_eval(range_val)
            except Exception:
                pass
        
        return range_val.split(',')
    
    return []

def row_has_overlap(range_val, intervals):
    """
    range_str: string like '12:63110111-63110111,5:132200733-132204564'
    intervals: list of tuples [(start, end), ...]
    """
    parts = normalize_ranges(range_val)
    
    for part in parts:
        try:
            # remove chromosome part
            _, coords = part.split(':')
            start, end = map(int, coords.split('-'))
        except ValueError:
            continue  # skip malformed entries

        # print(f'{start} {end}')
        # sys.exit(0)
        
        # check overlap with any interval
        for i_start, i_end in intervals:
            if start <= i_end and i_start <= end:
                return True
    
    return False

def filter_by_gene_location(df, location_col, location_array):
    '''
    Parameter:
    ----------
        df: pandas DataFrame
            Dataframe of patient ids, sv caller, sv type, and locations (default integer index).
        
        location_col: String
            Column with location data.

        location_array: Array of tuples.
            List of gene location tuples (idx 0 is start, idx 1 is end).

    Function:
    ---------
        - Finds gene locations that overlap any of the gene locations in location_array.

    Return:
    -------
        df_filtered: pandas DataFrame
            DataFrame filtered to rows that have overlap in location_array.
    '''

    return df[df[location_col].apply(lambda x: row_has_overlap(x, location_array))]

def filter_to_location_in_complex_sv_df(df, location_column, chr_filter):
    '''
    Parameter:
    ----------
        df: pandas DataFrame
            Dataframe of with location column.

        location_column: String
            Name of location column.

        chr_filter: String
            Chromosome to filter down to.

    Function:
    ---------
        - Filters dataframe to only rows that contain chr_filter.

    Return:
    -------
        df_filtered: pandas DataFrame
            Dataframe filtered down to rows containing chr_filter.
    '''
    rows = []
    for idx, row in df.iterrows():
        val = row[location_column]
        # val = "['8:92257042-92265814','9:92257042-92265814','10:92257042-92265814']"
        # val = 'chr8:92257042-92265814,chr9:92257042-92265814,chr10:92257042-92265814'

        # Convert stringified list → actual list
        if isinstance(val, str) and val.startswith('['):
            try:
                val = ast.literal_eval(val)
            except (ValueError, SyntaxError):
                pass  # fallback if it's not valid
            
        # If it's already a list/array-like, iterate it
        if isinstance(val, (list, tuple)):
            items = val
        else:
            # Otherwise treat as a single string
            items = [val]

        split_parts = []
        for item in items:
            split_parts.extend(str(item).split(','))

        cleaned = [p.strip() for p in split_parts if p.strip()]

        chr_list = []
        for item in cleaned:
            split = item.find(':')
            chrom = item[:split]
            chr_list.append(chrom.replace('chr', ''))

        

        if chr_filter in chr_list:
            rows.append(row)

    df_filtered = pd.DataFrame(rows)

    return df_filtered

def create_binary_dataframe(complex_sv_by_patient_df):
    '''
    Parameter:
    ----------
        complex_sv_by_patient_df: pandas DataFrame
            Dataframe of patient ids, sv caller, sv type, and locations (default integer index).

    Function:
    ---------
        - For every sv type in complex_sv_by_patient_df, call binarize_by_patient (this creates a binary present/not present dataframe of patient_id by sv_type).
        - If there are non-unique svs, combine columns.
        - Drop duplicates and return. 

    Return:
    -------
        binary_df: pandas DataFrame
            A binary present/not present dataframe of patient_id by sv_type.
    '''
    binary_df = pd.DataFrame()

    # Add patient_ids to new dataframe
    binary_df['patient_id'] = complex_sv_by_patient_df['sample_id']

    binary_df = binary_df.drop_duplicates()

    # Make a list of all 'unique' svs (svs like 'bfb' and 'BFB' will be counded seperately until later)
    unique_complex_sv_list = complex_sv_by_patient_df['SV_type'].unique().tolist()

    # For every sv create a binary column in new df
    for sv in unique_complex_sv_list:
        binary_df = binarize_by_patient(complex_sv_by_patient_df, binary_df, sv)

    # Make a list of non-unique svs
    unique_complex_sv_list_lower = [item.lower() for item in unique_complex_sv_list]
    counter = Counter(unique_complex_sv_list_lower)
    non_unique_svs = [item for item, count in counter.items() if count > 1]

    # Find non-unique svs and combine columns
    for sv in non_unique_svs:
        binary_df[sv] = (binary_df[sv] | binary_df[sv.upper()]).astype(int)
        binary_df.drop(columns=sv.upper(), inplace=True)

    # Drop duplicates
    binary_df = binary_df.drop_duplicates()

    # Replace any '-' with'_'
    binary_df.columns = binary_df.columns.str.replace('-', '_')

    return binary_df
    

def binarize_by_patient(complex_sv_by_patient_df, binary_df, sv):
    '''
    Parameters:
    -----------
        complex_sv_by_patient_df: pandas DataFrame
            Dataframe of patient ids, sv caller, sv type, and locations (default integer index).

        binary_df: pandas DataFrame
            A binary present/not present dataframe of patient_id by sv_type (default integer index). * Column width is variable *

        sv: String
            Type of sv to filter for in complex_sv_by_patient_df.

    Function:
    ---------
        - Create mask for sv, add a binary column based on mask to binary_df, and return.

    Return:
    -------
        binary_df: pandas DataFrame
            A binary present/not present dataframe of patient_id by sv_type.
    '''
    patients_with_sv = complex_sv_by_patient_df[complex_sv_by_patient_df['SV_type'] == sv]['sample_id'].unique()
    binary_df[sv] = binary_df['patient_id'].isin(patients_with_sv).astype(int)
    return binary_df

def prepare_sv_df_for_survival_analysis(binary_df, pluvicto_master_sheet_path):
    '''
    Parameters:
    -----------
        binary_df: pandas DataFrame
            A binary present/not present dataframe of patient_id by sv_type (default integer index).

        pluvicto_master_sheet_path: String
            Path to pluvicto master sheet. Patient id column should be labeled patient_id

    Function:
    ---------
        - Removes tags from patient ids using clean_df_ids function from entropy_by_chromosome.py.
        - Read in pluvicto master sheet and adds new patients to binary_df, sets their rows to 0 (no sv present), and returns it.

    Return:
    -------
        combined_df: pandas DataFrame
            binary_df but with all patients from pluvicto master sheet. 
    '''

    # Clean off tag from patient_ids
    binary_df_cleaned = clean_df_ids(binary_df, 'patient_id', "C1")

    # Add patients that are not already in binary_df
    pluvicto_df = pd.read_csv(pluvicto_master_sheet_path)
    combined_df = combine_columns(binary_df_cleaned, 'patient_id', pluvicto_df, 'Sample_ID')
    
    # Replace NaNs with 0s (sets all new patient rows to 0 indicating svs are not present)
    combined_df.fillna(0, inplace=True)

    return combined_df

def combine_columns(df_1, column_name_1, df_2, column_name_2):
    '''
    Parameters: 
    -----------
        df_1: pandas DataFrame
            DataFrame with at least one column named column_name_1.

        column_name_1: String
            Name of column that will be added onto.

        df_2: pandas DataFrame
            DataFrame with at least one column named column_name_2.

        column_name_2: String
            Name of column that will be used to add onto df_1.

    Function:
    ---------
        - * Comments say it all *

    Return:
    -------
        df_comnbined: pandas DataFrame
            df_1 with df_2's unique rows added on. NaNs fill in all other columns where new rows were added. 
    '''
    # Take just the column_name_2 column from df_2 and rename it to column_1_name
    df_2_slice = df_2[[column_name_2]]
    df_2_slice = df_2_slice.rename(columns={column_name_2: column_name_1})
    
    # Add rows from df_2 not already in df_1
    new_rows = df_2_slice[~df_2_slice[column_name_1].isin(df_1[column_name_1])]
    df_comnbined = pd.concat([df_1, new_rows], ignore_index=True)
    return df_comnbined

def filter_binary_df(output_path, entropy_csv_path, cycle, chr_to_pick, binary_df, tfx_df, high_entropy):  
    '''
    Parameters:
    -----------
        output_path: String
            Path to output directory.

        entropy_csv_path: String
            Path to entropy csv file.

        cycle: String
            Cycle to filter down to.

        chr_to_pick: String
            Name of entropy column to select for (specifict chromosome). 

        binary_df: pandas DataFrame
            Patient-by-SV binary matrix where rows are patients and columns are structural variants. Values are 1 (SV present) or 0 (SV absent). Must contain a patient identifier column.

        tfx_df: pandas DataFrame
            DataFrame of a column of sample id and TFx. 

        high_entropy: String
            Tag to filter dataframe down to high entropy group or to keep the whole cohort. 

    Function:
    ---------
        - Merge binary_df with tfx_df to filter down to patients with TFx present.
        - If high_entropy is not None, further filter to patients in top 50 percentile of entropy for specified chromosome.
        - Return the filtered binary dataframe.

    Return:
    -------
        binary_df_filtered: pandas DataFrame
            Patient-by-SV binary matrix filtered to patients with TFx present, and optionally to those with high entropy.
    '''
    # Sort down to patients with TFx present (or TFx ≥ threshold set before function call) - drop all tfx_df columns
    how = 'left' if binary_df.shape[0] < tfx_df.shape[0] else 'right'
    binary_df = pd.merge(binary_df, tfx_df, left_on='patient_id', right_on='Sample_ID', how=how).drop(columns=tfx_df.columns)

    # Set binary filtered 
    binary_df_filtered = binary_df

    # if high_entropy is not None filter to patients in top 50 percentile of entropy
    if high_entropy is not None:
        entropy_csv_file = os.path.join(output_path, entropy_csv_path)
        entropy_df = pd.read_csv(entropy_csv_file)
        # Filteres entropy df down to cycle then to specific chromosome column
        entropy_df = entropy_df[entropy_df["cycle"] == cycle]
        entropy_column = entropy_df[["patient_id", chr_to_pick]].copy()
    
        median = entropy_column['chr8'].median()
        entropy_column['temp'] = (entropy_column['chr8'] >= median).astype(int)
        entropy_column = entropy_column[entropy_column['temp'] == 1]
        entropy_column.drop('temp', axis=1, inplace=True)
        binary_df_filtered = binary_df[binary_df['patient_id'].isin(entropy_column['patient_id'])]

    return binary_df_filtered

def sv_entorpy_correlation(output_path, sub_dir, binary_df, entropy_column_name, high_entropy, tfx_thresholding_tag, filter_step):
    '''
    Parameters:
    -----------
        output_path: String
            Path to output directory.

        sub_dir: String
            Subdirectory within output_path for saving plots.

        binary_df: pandas DataFrame
            Patient-by-SV binary matrix where rows are patients and columns are structural variants.

        entropy_column_name: String
            Name of entropy column in binary_df.

        high_entropy: String
            Tag to filter dataframe down to high entropy group or to keep the whole cohort.

        tfx_thresholding_tag: String
            Tag describing TFx threshold applied.

        filter_step: String
            Filter criteria (e.g., 'high_vs_low_entropy').

    Function:
    ---------
        - Box-Cox normalize entropy column.
        - Calculate point-biserial correlations between each SV and entropy.
        - Apply FDR correction to p-values.
        - Generate and save boxplots for each SV-entropy correlation.

    Return:
    -------
        None
    '''
    # Normalize entropy column (add one to it incase there are 0s)
    shift = 1
    transformed_data, lambda_opt = stats.boxcox(binary_df[entropy_column_name] + shift)
    binary_df[entropy_column_name] = transformed_data

    corr_map = {}
    y = binary_df[entropy_column_name]

    for column in binary_df:
        if column == entropy_column_name or column == 'patient_id':
            continue

        binary_df_temp = binary_df.copy()

        # Remove patients that have high entropy, no sv and low entropy, with sv
        if filter_step == 'high_vs_low_entropy':
            # Median chr8 entropy value
            median = binary_df[entropy_column_name].median()

            binary_df_temp = binary_df[~(
                ((binary_df[column] == 0) & (binary_df[entropy_column_name] > median)) | 
                ((binary_df[column] == 1) & (binary_df[entropy_column_name] < median))
                )]
            
            y = binary_df_temp[entropy_column_name]
            
        x = binary_df_temp[column]
        corr = stats.pointbiserialr(list(x), list(y))
        corr_map[column] = corr

    
    # FDR Correction
    keys = list(corr_map.keys())
    pvalues = [corr_map[k].pvalue for k in keys]
    _, adj_pvalue, _, _ = multipletests(pvalues, method='fdr_bh')

    fdr_map = {
        k : {
            'statistic': corr_map[k].statistic,
            'adj_pvalue': adj_pval
        }
        for k, adj_pval in zip(keys, adj_pvalue)
        }

    for column in binary_df:
        if column == entropy_column_name or column == 'patient_id':
            continue
        # Plot correlation
        temp_df = binary_df[[column, entropy_column_name]]
        stat = fdr_map[column]['statistic']
        adj_pval = fdr_map[column]['adj_pvalue']

        # Remove patients that have high entropy, no sv and low entropy, with sv
        if filter_step == 'high_vs_low_entropy':
            # Median chr8 entropy value
            median = binary_df[entropy_column_name].median()

            temp_df = binary_df[~(
                ((binary_df[column] == 0) & (binary_df[entropy_column_name] > median)) | 
                ((binary_df[column] == 1) & (binary_df[entropy_column_name] < median))
                )]

        plot_title = f'{column} vs {entropy_column_name.title()} ({high_entropy.replace('_', ' ').title()}) {f'({filter_step.replace('_', ' ').title()} - {tfx_thresholding_tag})' if filter_step  else ''} Entropy Boxplot'
        x_axis_title = f'{column} Presence (0 = not present, 1 = present)'
        y_axis_title = f'{entropy_column_name.title()} Entropy (Box-Cox Normalized)'
        plot_filename = f'{column}_vs_{entropy_column_name}_entropy_correlation_{high_entropy}{tfx_thresholding_tag}{f'_{filter_step}'if filter_step  else ''}_box_plot.png'

        box_plot(output_path, sub_dir, temp_df, column, entropy_column_name, plot_title, x_axis_title, y_axis_title, plot_filename, stat, adj_pval)

def box_plot(output_path, sub_dir, df, x_column, y_column, title, x_label, y_label, file_name, corr_stat=None, p_val=None, show_points=True):
    '''
    Parameters:
    -----------
        output_path: String
            Path to output directory.

        sub_dir: String
            Subdirectory within output_path for saving plot.

        df: pandas DataFrame
            Data to plot.

        x_column: String
            Name of column for x-axis.

        y_column: String
            Name of column for y-axis.

        title: String
            Title for the plot.

        x_label: String
            Label for x-axis.

        y_label: String
            Label for y-axis.

        file_name: String
            Name for output file.

        corr_stat: Float, optional
            Correlation statistic to display on plot.

        p_val: Float, optional
            P-value to display on plot.

        show_points: Boolean, optional
            Whether to overlay points on boxplot. Default is True.

    Function:
    ---------
        - Create and style a boxplot with optional overlaid points.
        - Display statistical annotations if provided.
        - Save plot to specified directory.

    Return:
    -------
        None
    '''
    # -------------------------
    # Global style settings
    # -------------------------
    sns.set_theme(style="white", context="paper", font_scale=1.2)

    fig, ax = plt.subplots(figsize=(8, 6))

    # -------------------------
    # Boxplot
    # -------------------------
    sns.boxplot(x=x_column, y=y_column, data=df, ax=ax, width=0.6, fliersize=0, linewidth=1.5, boxprops=dict(alpha=0.6), medianprops=dict(color="black", linewidth=2), whiskerprops=dict(linewidth=1.5), capprops=dict(linewidth=1.5))

    # -------------------------
    # Optional: overlay points
    # -------------------------
    if show_points:
        sns.stripplot(x=x_column, y=y_column, data=df, ax=ax, color="black", alpha=0.6, jitter=True, size=4)

    # -------------------------
    # Labels & title
    # -------------------------
    ax.set_title(title, fontsize=14, weight="bold", pad=10)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)

    # -------------------------
    # Statistical annotation
    # -------------------------
    if corr_stat is not None and p_val is not None:
        stats_text = (
            rf"$\mathrm{{r_{{pb}}}} = {corr_stat:.2f}$" "\n"
            rf"$p = {p_val:.2e}$"
        )
        ax.text(
            0.98,
            0.98,
            stats_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10,
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white",
                edgecolor="black",
                linewidth=0.8
            )
        )

    # -------------------------
    # Clean up spines
    # -------------------------
    sns.despine(ax=ax)

    # -------------------------
    # Save
    # -------------------------
    file_dir = os.path.join(output_path, sub_dir)
    os.makedirs(file_dir, exist_ok=True)
    file_path = os.path.join(file_dir, file_name)

    plt.tight_layout()
    plt.savefig(file_path, dpi=300)   # use dpi=600 if needed
    plt.close()

def create_sv_kaplan_meier_plots(output_path, pluvicto_master_sheet_path, entropy_score_csv, binary_df, sequencing_type, high_entropy, tfx_thresholding_tag, sv_grouped_tag=False):
    '''
    Parameters:
    -----------
        output_path: String
            Path to output directory.

        pluvicto_master_sheet_path: String
            Path to pluvicto master sheet with clinical data including survival information.

        entropy_score_csv: String
            Name of CSV with entorpy values per chromosome across the samples (patient x chromosome).

        binary_df: pandas DataFrame
            Patient-by-SV binary matrix where rows are patients and columns are structural variants.

        sequencing_type: String
            Type of sequencing data (e.g., 'deep', 'ulp_curated').

        high_entropy: String
            Tag to filter dataframe down to high entropy group or to keep the whole cohort.

        tfx_thresholding_tag: String
            Tag describing TFx threshold applied.

        sv_grouped_tag: String
            Tag describing if svs are grouped.
            
    Function:
    ---------
        - For each SV in binary_df, fit Cox proportional hazard model.
        - Extract hazard ratios and confidence intervals.
        - Generate and save Kaplan-Meier survival curves stratified by SV presence.
        - Return the binary dataframe.

    Return:
    -------
        binary_df: pandas DataFrame
            Patient-by-SV binary matrix (unchanged).
    '''

    # Drop 'response' column
    if 'response' in binary_df.columns:
        binary_df = binary_df.drop('response', axis=1)

    sv_type_map ={
        'bfb': 'Breakage-fusion-brige',
        'dm': 'double minute',
        'pyrgo': 'High # junctions, low JCN',
        'tic': 'Templated insertion chains',
        'tyfonas': 'Low # junctions, high JCN',
    }

    merged_forest_df = pd.DataFrame()

    for sv in binary_df.columns:
        if sv == 'patient_id':
            continue

        grouping_df = binary_df[['patient_id', sv]]

        metric_to_use = 'Corrected_Copy_Number'

        if tfx_thresholding_tag:
            tfx_thresholding = '10_perc_tfx_threshold'
        else:
            tfx_thresholding = 'no_tfx_threshold'

        cox_forrest_output_path = f'entropy-analysis-{sequencing_type}-{metric_to_use}/cox_hardard_ratios/complex_sv{sv_grouped_tag}_plots/{entropy_score_csv}/{high_entropy}/{tfx_thresholding}'
        
        _, cox_model_temp, grouping_column, cox_model_summary = cox_proportional_hazard_model_from_csvs(
            output_path,
            cox_forrest_output_path, 
            f"data-tables/entropy-tables-{sequencing_type}-Corrected_Copy_Number/{entropy_score_csv}.csv", 
            pluvicto_master_sheet_path, 
            "chr8", 
            "C1", 
            "Death", 
            "survival_days", 
            ["TFx_C1"],
            metric_to_use,
            sequencing_type,
            grouping_df,
            high_entropy,
            'TFx_C1',
            0.10
        )

        # Merge all models together 
        merged_forest_df = pd.concat([merged_forest_df, cox_model_temp.summary.loc[[sv]].copy()])

        if sv in sv_type_map.keys():
            sv_explination = sv_type_map[sv]
        else:
            sv_explination = ""

        # Grabs and format HR and CI from cox model 
        row = cox_model_summary.index[cox_model_summary['covariate'] == sv][0]
        row = cox_model_summary.iloc[row]
        cox_model_hr_text = f'HR: {row["HR"]:.2f} (CI 95%: {row["HR_lower"]:.2f}–{row["HR_upper"]:.2f})'

        km_curve_output_path = f'entropy-analysis-{sequencing_type}-{metric_to_use}/kaplan-meier-curves/complex_sv{sv_grouped_tag}_km_curves/{entropy_score_csv}/{high_entropy}/{tfx_thresholding}'

        master_sheet = pluvicto_master_sheet_path
        master_sheet_df = pd.read_csv(master_sheet, index_col=0)

        print(f"sv: {sv}")
        print(f"grouping_column: \n{grouping_column}")

        # Only create km curve if grouping column is binary
        if grouping_column[sv].isin([0,1]).all():
            kaplan_meier_plot(
                    output_path,
                    km_curve_output_path, 
                    master_sheet_df,
                    "Death", 
                    "survival_days", 
                    sv,
                    grouping_column, 
                    "TFx_C1",
                    "median",
                    "Corrected_Copy_Number",
                    sequencing_type,
                    cox_model_hr_text,
                    None,
                    f'Kaplan-Meier Survival Curves based on presence of {sv} ({sv_explination}){f' - {high_entropy}' if high_entropy is not None else ''}',
                    None,
                    None,
                    sv, 
                    high_entropy
                )
            
    pvalues = merged_forest_df['p'].values
    _, pvalues_corrected, _, _ = multipletests(pvalues, method='fdr_bh')

    # Add corrected pvalues back to dataframe
    merged_forest_df['pvalues_corrected'] = pvalues_corrected

    merged_forest_df.to_csv(f'{output_path}/entropy-analysis-{sequencing_type}-{metric_to_use}/cox_hardard_ratios/complex_sv{sv_grouped_tag}_plots/{entropy_score_csv}/{high_entropy}/{tfx_thresholding}/combined_complex_sv_cox_proportional_hazard_model.csv')
            
    filename=f'combined_complex_sv_hazard_ratio_forest_plots'

    out_path, summary = cox_forest_plot_pub(
        merged_forest_df, 
        f'{output_path}/entropy-analysis-{sequencing_type}-{metric_to_use}/cox_hardard_ratios/complex_sv{sv_grouped_tag}_plots/{entropy_score_csv}/{high_entropy}/{tfx_thresholding}', 
        pval_col='pvalues_corrected',
        filename=filename,
        title='Cox Proportional Hazards Model (FDR corrected p-values)'
    )
        
    return binary_df
        
def subset_df_based_on_median(df, id_column, interest_column):
    '''
    Parameters:
    -----------
        df: pandas DataFrame
            Dataframe to subset.

        id_column: String
            Name of column containing identifiers.

        interest_column: String
            Name of column to calculate median from and filter by.

    Function:
    ---------
        - Calculate median of interest_column.
        - Filter dataframe to rows where interest_column >= median.
        - Reset index and return.

    Return:
    -------
        df_subset: pandas DataFrame
            Dataframe subset to rows with interest_column values above or equal to median.
    '''
    median = df[interest_column].median()
    ids_above_median = df[df[interest_column] >= median][id_column]
    df_subset = df[df[id_column].isin(ids_above_median)]
    df_subset.reset_index(inplace=True, drop=True)
    return df_subset

def scatter_plot_two_columns_one_dataframe(df, col_1_name, col_2_name, title=None, xlabel=None, ylabel=None, out_dir=".", file_name="scatter_plot", jitter=False):
    '''
     Parameters:
    -----------
        df_1: pandas DataFrame
            Dataframe with a columns named col_1_name and col_2_name (features x patients).

        col_1_name: String
            Name of column for x-axis.

        col_2_name: String
            Name of column for x-axis.

        title: String
            Title for the plot.

        xlabel: String
            Label for x-axis.

        ylabel: String
            Label for y-axis.

        out_dir: String
            Output directory for saving the figure.

        file_name: String
            Base name for the output file.

        jitter: boolean
            Add small vertical jitter to y values to unstack points.

    Function:
    ---------
        - Creates a scatter plot.
        - Saves plot.

    Returns:
    --------
        Void: saves figure.
    '''
    os.makedirs(out_dir, exist_ok=True)

    # Clean paired data
    plot_df = df[[col_1_name, col_2_name]].dropna()
    if plot_df.shape[0] == 0:
        raise ValueError(f"No valid paired data found for '{col_1_name}' and '{col_2_name}'")

    x = plot_df[col_1_name].values
    y = plot_df[col_2_name].values

    # ----------------------------
    # JITTER (for ordinal/count y)
    # ----------------------------
    if jitter:
        rng = np.random.default_rng(42)
        y_plot = y + rng.normal(0, 0.08, size=len(y))  # small vertical jitter only
    else:
        y_plot = y
        
    # ----------------------------
    # STATS
    # ----------------------------
    slope, intercept, r_value, p_value, _ = linregress(x, y)
    rho, rho_p = spearmanr(x, y)

    # ----------------------------
    # STYLE
    # ----------------------------
    plt.rcParams.update({
        "font.size": 11,
        "font.family": "sans-serif",
        "axes.linewidth": 1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42
    })

    fig, ax = plt.subplots(figsize=(4.5, 4.5), dpi=600)

    # ----------------------------
    # SCATTER
    # ----------------------------
    ax.scatter(
        x,
        y_plot,
        s=18,
        alpha=0.75,
        edgecolor="black",
        linewidth=0.3,
        zorder=2
    )

    # ----------------------------
    # REGRESSION LINE (de-emphasized)
    # ----------------------------
    x_fit = np.linspace(x.min(), x.max(), 200)
    y_fit = slope * x_fit + intercept

    ax.plot(
        x_fit,
        y_fit,
        linestyle="--",
        linewidth=1.2,
        alpha=0.6,
        zorder=3
    )

    # ----------------------------
    # HELPER FUNCTION 
    # ----------------------------
    def is_discrete(arr):
        arr = np.asarray(arr)
        return np.allclose(arr, np.round(arr))
    
    # ----------------------------
    # AXES LIMITS (safe padding)
    # ----------------------------
    def get_limits(arr):
        span = arr.max() - arr.min()
        pad = span * 0.05 if span > 0 else 1
        return arr.min() - pad, arr.max() + pad

    ax.set_xlim(*get_limits(x))
    ax.set_ylim(*get_limits(y))

    # ----------------------------
    # SMART TICK HANDLING (FIXED)
    # ----------------------------
    x_discrete = is_discrete(x)
    y_discrete = is_discrete(y)

    if x_discrete:
        ax.set_xticks(np.unique(x))

    if y_discrete:
        ax.set_yticks(np.unique(y))

    # ----------------------------
    # LABELS
    # ----------------------------
    ax.set_xlabel(xlabel if xlabel else col_1_name, labelpad=6)
    ax.set_ylabel(ylabel if ylabel else col_2_name, labelpad=6)

    if title:
        ax.set_title(title, pad=10)

    # ----------------------------
    # CLEAN LOOK
    # ----------------------------
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.35, zorder=1)

    # ----------------------------
    # ANNOTATION (smart placement)
    # ----------------------------
    x_mid = np.median(x)
    y_mid = np.median(y)

    if np.mean((x < x_mid) & (y > y_mid)) < 0.25:
        text_x, text_y = 0.05, 0.95
        ha, va = "left", "top"
    else:
        text_x, text_y = 0.95, 0.05
        ha, va = "right", "bottom"

    annotation = (
        f"Spearman ρ = {rho:.2f}\n"
        f"p = {rho_p:.2e}\n"
    )

    ax.text(
        text_x,
        text_y,
        annotation,
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=10
    )

    # Regression label (small, unobtrusive)
    ax.text(
        0.98,
        0.02,
        f"y = {slope:.2f}x + {intercept:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        alpha=0.8
    )

    plt.tight_layout()

    # ----------------------------
    # SAVE (journal standard: PNG + PDF)
    # ----------------------------
    png_path = os.path.join(out_dir, f"{file_name}.png")
    pdf_path = os.path.join(out_dir, f"{file_name}.pdf")

    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

def box_plot_two_columns_one_dataframe(df, col_1_name, col_2_name, title=None, xlabel=None, ylabel=None, out_dir=".", file_name="box_plot", jitter=False, compare_groups=None):
    '''
     Parameters:
    -----------
        df_1: pandas DataFrame
            Dataframe with a columns named col_1_name and col_2_name (features x patients).

        col_1_name: String
            Name of column for x-axis.

        col_2_name: String
            Name of column for x-axis.

        title: String
            Title for the plot.

        xlabel: String
            Label for x-axis.

        ylabel: String
            Label for y-axis.

        out_dir: String
            Output directory for saving the figure.

        file_name: String
            Base name for the output file.

        jitter: boolean
            Add small vertical jitter to y values to unstack points.

        compare_pairs: List of a tuple of lists (e.g. [([0,1],[2,3,4,5])]
            Mode to compare different groups of bars together using mann-whitney u test.

    Function:
    ---------
        - Creates a scatter plot.
        - Saves plot.

    Returns:
    --------
        Void: saves figure.
    '''
    os.makedirs(out_dir, exist_ok=True)

    # Clean data
    plot_df = df[[col_1_name, col_2_name]].dropna()
    if plot_df.shape[0] == 0:
        raise ValueError(f"No valid data for '{col_1_name}' and '{col_2_name}'")

    groups = plot_df.groupby(col_1_name)[col_2_name]
    unique_x = sorted(plot_df[col_1_name].unique())
    data = [groups.get_group(g).values for g in unique_x]

    # ----------------------------
    # STYLE
    # ----------------------------
    plt.rcParams.update({
        "font.size": 11,
        "font.family": "sans-serif",
        "axes.linewidth": 1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42
    })

    fig, ax = plt.subplots(figsize=(5, 4.5), dpi=600)

    # ----------------------------
    # BOX PLOT
    # ----------------------------
    box = ax.boxplot(
        data,
        patch_artist=True,
        widths=0.6,
        showfliers=False,
        medianprops=dict(linewidth=1.5)
    )

    for patch in box["boxes"]:
        patch.set_facecolor("white")
        patch.set_edgecolor("black")
        patch.set_linewidth(1)

    for element in ["whiskers", "caps", "medians"]:
        for item in box[element]:
            item.set_color("black")
            item.set_linewidth(1)

    # ----------------------------
    # OPTIONAL JITTER
    # ----------------------------
    if jitter:
        rng = np.random.default_rng(42)

        for i, g in enumerate(unique_x):
            y_vals = groups.get_group(g).values
            x_vals = rng.normal(i + 1, 0.06, size=len(y_vals))

            ax.scatter(
                x_vals,
                y_vals,
                s=15,
                alpha=0.6,
                edgecolor="black",
                linewidth=0.3,
                zorder=2
            )
    # ----------------------------
    # HELPER FUNCTIONS
    # ----------------------------
    def format_group(groups):
        return "+".join(map(str, groups))
    
    # ----------------------------
    # AXES
    # ----------------------------
    ax.set_xticks(range(1, len(unique_x) + 1))
    ax.set_xticklabels(unique_x)

    ax.set_xlabel(xlabel if xlabel else col_1_name, labelpad=6)
    ax.set_ylabel(ylabel if ylabel else col_2_name, labelpad=6)

    # ----------------------------
    # DYNAMIC TITLE
    # ----------------------------

    # Title
    if compare_groups is None:
        plot_title = title
    else:
        comp_parts = []
        for a, b in compare_groups:
            comp_parts.append(f"{format_group(a)}_vs_{format_group(b)}")
        comp_str = "__".join(comp_parts)
        plot_title = f"{title} (comp_str)"

    if title:
        ax.set_title(plot_title, pad=10)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.35)

    # ----------------------------
    # SAMPLE SIZE LABELS
    # ----------------------------
    for i, g in enumerate(unique_x):
        n = len(groups.get_group(g))
        ax.text(
            i + 1,
            ax.get_ylim()[0],
            f"n={n}",
            ha="center",
            va="top",
            fontsize=9
        )

    # ----------------------------
    # STATISTICAL TESTING (NEW)
    # ----------------------------
    if compare_groups is not None:

        y_min, y_max = ax.get_ylim()
        y_range = y_max - y_min
        y_offset = y_range * 0.12

        for i, (group_a, group_b) in enumerate(compare_groups):

            # collapse data
            data_a = plot_df[plot_df[col_1_name].isin(group_a)][col_2_name].values
            data_b = plot_df[plot_df[col_1_name].isin(group_b)][col_2_name].values

            stat, p = mannwhitneyu(data_a, data_b, alternative="two-sided")

            # x positions (center of groups)
            x_positions = sorted(plot_df[col_1_name].unique())

            def center(groups):
                return np.mean([x_positions.index(g) + 1 for g in groups])

            x1 = center(group_a)
            x2 = center(group_b)

            y = y_max + (i + 1) * y_offset

            # significance label
            if p < 0.001:
                sig = "***"
            elif p < 0.01:
                sig = "**"
            elif p < 0.05:
                sig = "*"
            else:
                sig = "ns"

            # bracket
            ax.plot([x1, x1, x2, x2],
                    [y, y + y_offset*0.2, y + y_offset*0.2, y],
                    lw=1, c="black")

            ax.text((x1 + x2) / 2, y + y_offset*0.25,
                    f"{sig} (p={p:.2e})",
                    ha="center", va="bottom", fontsize=9)

        ax.set_ylim(y_min, y_max + (len(compare_groups) + 1) * y_offset)

    plt.tight_layout()

    # ----------------------------
    # DYNAMIC FILE NAME
    # ----------------------------

    if compare_groups is None:
        group_str = f"{col_1_name}_all"
        comp_str = "no_comparison"
    else:
        comp_parts = []
        for a, b in compare_groups:
            comp_parts.append(f"{format_group(a)}_vs_{format_group(b)}")
        comp_str = "__".join(comp_parts)

        group_str = f"{col_1_name}_{comp_str}"

    final_file_name = file_name + "_" + group_str

    # ----------------------------
    # SAVE
    # ----------------------------
    png_path = os.path.join(out_dir, f"{final_file_name}.png")
    pdf_path = os.path.join(out_dir, f"{final_file_name}.pdf")

    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()