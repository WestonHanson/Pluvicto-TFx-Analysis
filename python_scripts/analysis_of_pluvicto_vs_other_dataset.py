# Author: Weston Hanson
# Place: Fred Hutch Cancer Center, Seattle, WA
# Date Created: 04/17/26
# Purpose: Analyzes and visualizes 

from python_scripts.imports import *

def prepare_distribution_data(csv_path, id_col, list_of_cols, cycle_col, cycle_filter, col_tag):
    '''
    Parameters:
    -----------
        csv_path: String
            Path to CSV file.

        id_col: String
            Name of id column in df.

        list_of_cols: List
            List of column names to keep in df.

        cycle_col: String
            Name of cycle column, 

        cycle_filter: String
            Cycle to filter down to.

        col_tag: String
            Tag to add to the columns to keep to identify which dataset they come from.

    Function:
    ---------
        - Read in CSV.
        - Subset to cycle_filter.
        - Subset the df to the list of columns passed in.
        - Adds tag to columns. 
        - Returns the df.

    Returns:
    --------
        df: pandas DataFrame
            DataFrame with only the columns passed in. 

    '''
    # Read in CSVs
    df = pd.read_csv(csv_path)

    # Filter df
    if cycle_col is not None and cycle_filter is not None:
        df = df[df[cycle_col] == cycle_filter]
    df = df[[id_col] + list_of_cols]

    # Add tag to columns besides id_col
    df.columns = [f'{col}_{col_tag}' if col != id_col else col for col in df.columns]

    return df, col_tag



def compare_column_violin_between_dfs(df_1, col_name_1, df_2, col_name_2, df_3, col_name_3, title=None, xlabel=None, ylabel=None, out_dir=".", file_name="violin_plot"):
    '''
    Parameters:
    -----------
        df_1: pandas DataFrame
            Dataframe with a column named col_name_1 (features x patients).

        col_name_1: String
            Name of column in df_1 to compare.

        df_2: pandas DataFrame
            Dataframe with a column named col_name_2 (features x patients).

        col_name_2: String
            Name of column in df_2 to compare.

        df_3: pandas DataFrame
            Dataframe with a column named col_name_3 (features x patients).

        col_name_3: String
            Name of column in df_3 to compare.

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

    Function:
    ---------
        - Create a publication-ready violin plot comparing a column from two dataframes.
        - Includes:
            - Median
            - Quartiles (25th, 75th percentiles)
            - Whiskers (1.5 * IQR)

    Returns:
    --------
        Void: saves figure.
    '''

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Drop NA values
    data1 = df_1[col_name_1].dropna().values
    data2 = df_2[col_name_2].dropna().values
    data3 = df_3[col_name_3].dropna().values

    if len(data1) == 0:
        raise ValueError(f"No valid data found in column '{col_name_1}' (df_1)")
    if len(data2) == 0:
        raise ValueError(f"No valid data found in column '{col_name_2}' (df_2)")
    if len(data3) == 0:
        raise ValueError(f"No valid data found in column '{col_name_3}' (df_3)")

    # Helper to compute stats
    def compute_stats(data):
        q1 = np.percentile(data, 25)
        median = np.percentile(data, 50)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1

        lower_whisker = np.min(data[data >= q1 - 1.5 * iqr])
        upper_whisker = np.max(data[data <= q3 + 1.5 * iqr])

        return q1, median, q3, lower_whisker, upper_whisker

    stats1 = compute_stats(data1)
    stats2 = compute_stats(data2)
    stats3 = compute_stats(data3)

    # Figure setup (journal style)
    plt.rcParams.update({
        "font.size": 12,
        "font.family": "sans-serif",
        "axes.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42
    })

    fig, ax = plt.subplots(figsize=(5, 6))

    # Violin plot (two groups)
    parts = ax.violinplot(
        [data1, data2, data3],
        positions=[1, 2, 3],
        showmeans=False,
        showmedians=False,
        showextrema=False
    )

    # Style violins
    for pc in parts['bodies']:
        pc.set_facecolor("#82CAFF")
        pc.set_edgecolor("black")
        pc.set_alpha(0.8)
        pc.set_linewidth(1.2)

    # Function to draw box + whiskers
    def draw_box(ax, x, stats):
        q1, median, q3, lower, upper = stats

        # Quartile box
        ax.vlines(x, q1, q3, color='black', linewidth=6)

        # Median
        ax.scatter(x, median, color='white', edgecolor='black', zorder=3, s=60)

        # Whiskers
        ax.vlines(x, lower, upper, color='black', linewidth=1.5)

        # Caps
        ax.hlines(lower, x - 0.05, x + 0.05, color='black', linewidth=1.5)
        ax.hlines(upper, x - 0.05, x + 0.05, color='black', linewidth=1.5)

    draw_box(ax, 1, stats1)
    draw_box(ax, 2, stats2)
    draw_box(ax, 3, stats3)

    # Formatting
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels([col_name_1, col_name_2, col_name_3])

    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)

    if title:
        ax.set_title(title)

    # Remove top/right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    # Save
    png_path = os.path.join(out_dir, f"{file_name}_violin_plot.png")
    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.close()

def compare_column_density_distribution_between_dfs(df_1, col_name_1, df_2, col_name_2, df_3, col_name_3, title=None, xlabel=None, ylabel=None, out_dir=".", file_name="density_distribution"):
    '''
     Parameters:
    -----------
        df_1: pandas DataFrame
            Dataframe with a column named col_name_1 (features x patients).

        col_name_1: String
            Name of column in df_1 to compare.

        df_2: pandas DataFrame
            Dataframe with a column named col_name_2 (features x patients).

        col_name_2: String
            Name of column in df_2 to compare.

        df_3: pandas DataFrame
            Dataframe with a column named col_name_3 (features x patients).

        col_name_3: String
            Name of column in df_3 to compare.

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

    Function:
    ---------
        - Creates a denstiy distribution plot of both dataframes overlayed.
        - Saves plot.

    Returns:
    --------
        Void: saves figure.
    '''

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Extract and clean data
    data1 = df_1[col_name_1].dropna().values
    data2 = df_2[col_name_2].dropna().values
    data3 = df_3[col_name_3].dropna().values

    if len(data1) == 0:
        raise ValueError(f"No valid data found in column '{col_name_1}' (df_1)")
    if len(data2) == 0:
        raise ValueError(f"No valid data found in column '{col_name_2}' (df_2)")
    if len(data3) == 0:
        raise ValueError(f"No valid data found in column '{col_name_3}' (df_3)")

    # Plot style (publication-ready)
    plt.rcParams.update({
        "font.size": 12,
        "axes.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42
    })

    fig, ax = plt.subplots(figsize=(6, 5), dpi=600)

    # Determine shared x range
    xmin = min(data1.min(), data2.min(), data3.min())
    xmax = max(data1.max(), data2.max(), data3.max())
    x = np.linspace(xmin, xmax, 500)

    # Colors (same palette you used)
    colors = [
        "#0044FF",  # black
        "#E69F00",  # orange
        "#56B4E9",
    ]

    # KDE computation
    kde1 = gaussian_kde(data1)
    kde2 = gaussian_kde(data2)
    kde3 = gaussian_kde(data3)

    y1 = kde1(x)
    y2 = kde2(x)
    y3 = kde3(x)

    # Plot
    ax.plot(x, y1, linewidth=2, label=col_name_1, color=colors[0])
    ax.plot(x, y2, linewidth=2, label=col_name_2, color=colors[1])
    ax.plot(x, y3, linewidth=2, label=col_name_3, color=colors[2])

    # Labels
    ax.set_xlabel(xlabel if xlabel else col_name_1)
    ax.set_ylabel(ylabel if ylabel else "Density")

    if title:
        ax.set_title(title)

    # Style cleanup
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(frameon=False)

    plt.tight_layout()

    # Save
    png_path = os.path.join(out_dir, f"{file_name}.png")
    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.close()

def add_column_to_df(df_1, id_column_1, df_2, id_column_2, column_to_add):
    '''
    Parameters:
    -----------
        df_1: pandas DataFrame
            DataFrame 1.

        id_column_1: String
            Name of id column for df_1.

        df_2: pandas DataFrame
            DataFrame 2.

        id_column_2: String
            Name of id column for df_2.

        column_to_add: String
            Name of column from df_2 to add to df_1.

    Function:
    ---------
        - Add column_to_add to df_1 and return.

    Return:
    -------
        merged_df: pandas DataFrame
            A dataframe with new column.
    '''
    df_2 = df_2[[id_column_2, column_to_add]]
    merged_df = pd.merge(df_1, df_2, left_on=id_column_1, right_on=id_column_2, how='inner').drop(columns=id_column_2)

    return merged_df

def final_gene_matrix_modified(annotation_df, annotation_gene_column, CN_matrix_path, CN_matrix_id_column, cycle, meta_data_path, meta_data_id_column, ploidy_column, filter_ids=True):
    '''
    Parameters:
    -----------
        annotation_df: pandas DataFrame
            Filtered annotation gene list.

        annotation_gene_column: String
            Name of gene id column. 
        
        CN_matrix_path: List
            List of file paths to CN gene matrix.

        CN_matrix_id_column: String
            Name of patient id column for CN_matrix_path.

        cycle: String
            Cycle to filter patient ids down to.

        meta_data_path: List
            List of file paths to Titan purity ploidy list.

        meta_data_id_column: String
            Name of patient id column for meta_data_path.

        ploidy_column: String
            Name of column to get ploidy values in ploidy file.

        filter_ids: boolean
            True/False to remove ids.

        clean_ids: boolean
            True/False to pass matrix through clean ids function.

        * Note: Make sure titan_matrix_path_list and meta_data_path are respective to each other (eg. idexes of paths should align) *

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
    # Filter matrix to gene list
    matrix = pd.read_table(CN_matrix_path, sep='\t')
    mask = matrix.columns.isin(gene_list)
    mask[0] = True # Make first index True to keep sample id
    matrix_filtered = matrix.loc[:, mask]

    if filter_ids:
        matrix_filtered_cleaned = matrix_filtered[matrix_filtered[CN_matrix_id_column].str.contains(cycle, na=False)]
    else:
        matrix_filtered_cleaned = matrix_filtered

    # Add ploidy as a column to matrix
    if meta_data_path is not None and ploidy_column is not None:
        ploidy_df = pd.read_csv(meta_data_path)
        
        if ploidy_column not in ploidy_df.columns:
            joined_matrix = pd.concat([joined_matrix, matrix_filtered_cleaned], ignore_index=True)
        else:
            ploidy_df_filtered = ploidy_df[ploidy_df[meta_data_id_column].str.contains(cycle, na=False)]
            ploidy_df_filtered = ploidy_df_filtered[[meta_data_id_column, ploidy_column]]

            matrix_and_ploidy = pd.merge(matrix_filtered_cleaned, ploidy_df_filtered, right_on=CN_matrix_id_column, left_on=meta_data_id_column, how='left')

            joined_matrix = pd.concat([joined_matrix, matrix_and_ploidy], ignore_index=True)
    else:
        joined_matrix = pd.concat([joined_matrix, matrix_filtered_cleaned], ignore_index=True)

    return joined_matrix