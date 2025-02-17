import pandas as pd
import streamlit as st
from io import BytesIO
import altair as alt
import os
from scipy.stats import fisher_exact
from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_object_dtype,
)

from statsmodels.stats.proportion import proportions_ztest
from scipy.stats import chi2_contingency


# Streamlit app title
st.title("Combine SafeSeq Results into a Single Output File")

# File upload widgets
raw_summary_file = st.file_uploader("Upload the raw summary file (.tab format)", type="tab")
run_summary_file = st.file_uploader("Upload the run summary file (.xlsm format)", type="xlsm")
sample_list_file = st.file_uploader("Upload the sample list file (.xlsm format)", type="xlsm")
#CHIP data file would not be always available, if not, the whole function would not be changed. 
chipdatafile = st.file_uploader("Upload the CHIP data file (.tab format)", type="tab")

# Text input to specify custom output file path and name
output_file_path = st.text_input("Enter the full path for the output file (e.g., /Users/username/Downloads/Combined_Output.xlsx)", "/Users/dafualt_path/Combined_Output.xlsx")

# Function to display DataFrame
def display_dataframe(df, title):
    st.subheader(title)
    st.dataframe(df)

def display_dataframe_with_filter(df, title):
    st.subheader(title)
    df = filter_dataframe(df)
    st.dataframe(df)

def read_raw_data(tabfile):
    df = pd.read_csv(tabfile, sep='\t')
    df = df[~df['SampleId'].str.contains('PC')]
    df = df[~df['SampleId'].str.contains('NC')]
    return df

def calculate_p_value_and_odds_ratio_MM(row):
    tumor_var = row['MM_Safeseq_t']
    tumor_ref = row['MM_Safeseq_ref']
    normal_var = row['MM_BC_t']
    normal_ref = row['MM_BC_ref']
    # Add continuity correction to avoid zeros
    tumor_var += 0.1
    tumor_ref += 0.1
    normal_var += 0.1
    normal_ref += 0.1
    # Construct contingency table
    contingency_table = [
        [tumor_var, tumor_ref],
        [normal_var, normal_ref]
        ]
    # Perform Fisher's Exact Test
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
    return pd.Series({'p_value': p_value, 'odds_ratio': odds_ratio})

def two_proportion_ztest(row):
    tumor_var = row['MM_Safeseq_t']
    tumor_ref = row['MM_Safeseq_ref']
    normal_var = row['MM_BC_t']
    normal_ref = row['MM_BC_ref']    
    # Number of "successes" in each group
    count = np.array([tumor_var, normal_var])
    # Total observations in each group
    nobs = np.array([tumor_var + tumor_ref, normal_var + normal_ref])
    stat, p_value = proportions_ztest(count, nobs, alternative='larger')
    return pd.Series({'z_stat': stat, 'p_value': p_value})

def chi_squared_test(row):
    tumor_var = row['MM_Safeseq_t']
    tumor_ref = row['MM_Safeseq_ref']
    normal_var = row['MM_BC_t']
    normal_ref = row['MM_BC_ref']
    contingency_table = [[tumor_var, tumor_ref],
                         [normal_var, normal_ref]]
    chi2, p_value, dof, expected = chi2_contingency(contingency_table)
    return pd.Series({'chi2': chi2, 'p_value': p_value, 'dof': dof})


def update_result_review(df):
    """
    Updates the dataframe based on specific conditions for 'fisher_p_value', 'fisher_odds_ratio', and 'Call' columns.
    Parameters:
        df (pd.DataFrame): The input dataframe.
    Returns:
        pd.DataFrame: The updated dataframe.
    """
    # Define conditions
    condition = (
        (df['fisher_p_value'] > 0.01) | (df['fisher_odds_ratio'] < 2)
    ) & (df['Call'] == 'MD')
    # Identify rows to be changed
    changed_rows = df[condition]
    # Print rows that will be changed
    print("Rows that have been changed:")
    print(changed_rows)
    # Apply updates based on conditions
    df.loc[condition, 'Call'] = 'NMD'
    df.loc[condition, 'Comment Call (Internal)'] = (
        df['Comment Call (Internal)'].fillna('') +
        " BC experiment did not support this mutant being detected"
    )
    return df, changed_rows

def chip_data_process(chipdatafile, df_ResultReview_import):
    df_chip = read_raw_data(chipdatafile)
    df_chip['SampleId']= df_chip['SampleId'].astype(str).str.replace(r'BC', '', regex=True)
    df_chip_filter=df_chip[df_chip['Gene Name'].isin(['TP53', 'KRAS'])]
    df_chip_filter_selected = df_chip_filter[['SampleId','CDSChange', 'AAChange', '#UIDs/Amplicon','#Supermutants', 'Gene Name', 'Call', 'MAF[%]', 'MutantMolecules', 'GE']]
    df_chip_filter_selected = df_chip_filter_selected.rename(columns={"SampleId":"SampleID"})
    df_ResultReview_import_selected = df_ResultReview_import[((df_ResultReview_import['Call']=='MD') & (df_ResultReview_import['Gene Name'].isin(['TP53', 'KRAS'])))]
    merged_df = pd.merge(
        df_ResultReview_import_selected,
        df_chip_filter_selected,
        on=['SampleID', 'CDSChange', 'AAChange'],
        suffixes=('_tumor', '_normal'),
        how='left'
        )
    merged_df_wBC = merged_df[~merged_df['#UIDs/Amplicon_normal'].isna()]
    merged_df_nBC = merged_df[merged_df['#UIDs/Amplicon_normal'].isna()]
    merged_df_wBC['MM_Safeseq_t'] = merged_df_wBC['#Supermutants_tumor']*merged_df_wBC['GE_tumor']/merged_df_wBC['#UIDs/Amplicon_tumor'] 
    merged_df_wBC['MM_Safeseq_ref'] = (merged_df_wBC['#UIDs/Amplicon_tumor'] - merged_df_wBC['#Supermutants_tumor'])*merged_df_wBC['GE_tumor']/merged_df_wBC['#UIDs/Amplicon_tumor']
    merged_df_wBC['MM_BC_t'] = merged_df_wBC['#Supermutants_normal']*merged_df_wBC['GE_normal']/merged_df_wBC['#UIDs/Amplicon_normal'] 
    merged_df_wBC['MM_BC_ref'] = (merged_df_wBC['#UIDs/Amplicon_normal'] - merged_df_wBC['#Supermutants_normal'])*merged_df_wBC['GE_normal']/merged_df_wBC['#UIDs/Amplicon_normal']
    # Combine the columns into one list
    columns_to_fill = ['MM_Safeseq_t', 'MM_Safeseq_ref', 'MM_BC_t', 'MM_BC_ref']
    merged_df_wBC[['fisher_p_value', 'fisher_odds_ratio']] = merged_df_wBC.apply(calculate_p_value_and_odds_ratio_MM, axis=1)
    rows_to_update = merged_df_wBC[(merged_df_wBC['fisher_p_value'] > 0.01) | (merged_df_wBC['fisher_odds_ratio'] < 2)]
    merged_df_wBC['key'] = merged_df_wBC['SampleID'] + '_' + merged_df_wBC['CDSChange'] + '_' + merged_df_wBC['AAChange']
    df_ResultReview_import['key'] = df_ResultReview_import['SampleID'] + '_' + df_ResultReview_import['CDSChange'] + '_' + df_ResultReview_import['AAChange']
    df_ResultReview_import_wBC = df_ResultReview_import.merge(
        merged_df_wBC[['key', 'fisher_p_value', 'fisher_odds_ratio']],
        on='key',
        how='left'
        )
    df_ResultReview_import_wBC_upaded, changed_rows = update_result_review(df_ResultReview_import_wBC)
    if changed_rows.shape == 0:  # Check if the number of rows is 0
        print("No data need to update in this Run")
    df_ResultReview_import.drop(columns=['key'], inplace=True)
    return df_ResultReview_import_wBC_upaded, merged_df_wBC, rows_to_update, merged_df_nBC


def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a UI on top of a dataframe to let viewers filter columns.
    Args:
        df (pd.DataFrame): Original dataframe
    Returns:
        pd.DataFrame: Filtered dataframe
    """
    modify = st.checkbox("Add filters")
    if not modify:
        return df
    df = df.copy()
    # Attempt to convert datetimes into a standard format (datetime without timezone)
    for col in df.columns:
        if is_object_dtype(df[col]):
            try:
                df[col] = pd.to_datetime(df[col])
            except Exception:
                pass
        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)
    modification_container = st.container()
    with modification_container:
        to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
        for column in to_filter_columns:
            left, right = st.columns((1, 20))
            left.write("↳")
            # If fewer than 10 unique values, treat as categorical
            if is_categorical_dtype(df[column]) or df[column].nunique() < 10:
                user_cat_input = right.multiselect(
                    f"Values for {column}",
                    df[column].unique(),
                    default=list(df[column].unique()),
                )
                df = df[df[column].isin(user_cat_input)]
            elif is_numeric_dtype(df[column]):
                _min = float(df[column].min())
                _max = float(df[column].max())
                # Ensure step is not zero
                step = (_max - _min) / 100 if _max != _min else 1.0
                user_num_input = right.slider(
                    f"Values for {column}",
                    _min,
                    _max,
                    (_min, _max),
                    step=step,
                )
                df = df[df[column].between(*user_num_input)]
            elif is_datetime64_any_dtype(df[column]):
                user_date_input = right.date_input(
                    f"Values for {column}",
                    value=(df[column].min(), df[column].max()),
                )
                if len(user_date_input) == 2:
                    start_date, end_date = map(pd.to_datetime, user_date_input)
                    df = df.loc[df[column].between(start_date, end_date)]
            else:
                user_text_input = right.text_input(f"Substring or regex in {column}")
                if user_text_input:
                    # Use na=False to avoid errors with missing data
                    df = df[df[column].astype(str).str.contains(user_text_input, na=False)]
    return df


def plot_chip_data(df):
    """
    Create a grouped bar chart to visualize tumor vs. normal counts and show
    fisher_p_value & fisher_odds_ratio in tooltips.
    """
    required_cols = [
        'SampleID',
        '#UIDs/Amplicon_tumor', '#Supermutants_tumor',
        '#UIDs/Amplicon_normal', '#Supermutants_normal',
        'fisher_p_value', 'fisher_odds_ratio'
    ]
    # Check that the dataframe has all required columns
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        st.warning(f"Missing columns in rows_to_update: {missing_cols}")
        return
    # Melt the dataframe so that each count type is in its own row
    # We'll keep SampleID, fisher_p_value, and fisher_odds_ratio as "id_vars"
    df_long = df.melt(
        id_vars=['SampleID', 'fisher_p_value', 'fisher_odds_ratio'],
        value_vars=[
            'MutantMolecules_tumor',
            'MutantMolecules_normal'
        ],
        var_name='metric',
        value_name='counts'
    )
    # Build an Altair bar chart
    chart = (
        alt.Chart(df_long)
        .mark_bar()
        .encode(
            x=alt.X('SampleID:N', title='SampleID'),
            y=alt.Y('counts:Q', title='Counts'),
            color=alt.Color('metric:N', title='Metric'),
            tooltip=['SampleID', 'metric', 'counts', 'fisher_p_value', 'fisher_odds_ratio']
        )
        .properties(width=600, height=400)
        .interactive()
    )
    st.altair_chart(chart, use_container_width=True)

if "df1_summary_tab" not in st.session_state:
    st.session_state.df1_summary_tab = None

if "df_RunSumm" not in st.session_state:
    st.session_state.df_RunSumm = None

if "sample_mapping_df" not in st.session_state:
    st.session_state.sample_mapping_df = None

if "df_ResultReview_import" not in st.session_state:
    st.session_state.df_ResultReview_import = None

if "merged_df_wBC" not in st.session_state:
    st.session_state.merged_df_wBC = None

if "rows_to_update" not in st.session_state:
    st.session_state.rows_to_update = None

if "merged_df_nBC" not in st.session_state:
    st.session_state.merged_df_nBC = None

# Button to process the files
if st.button("Combine Files") and raw_summary_file and run_summary_file and sample_list_file:
    # Read and preprocess the raw summary file
    df1_summary_tab = pd.read_csv(raw_summary_file, sep='\t')
    df1_summary_tab = df1_summary_tab[~df1_summary_tab['SampleId'].str.contains('PC')]
    df1_summary_tab = df1_summary_tab[~df1_summary_tab['SampleId'].str.contains('NTC')]
    st.session_state.df1_summary_tab = df1_summary_tab
    Samplelist = list(set(df1_summary_tab['SampleId'].to_list()))
        
    # Read the run summary file
    df_RunSumm = pd.read_excel(run_summary_file, sheet_name='7 Run Summary').dropna()
    df_AE_input = pd.read_excel(run_summary_file, sheet_name='9 SafeSEQ AE input')
    FlowCellID = df_AE_input.loc[df_AE_input['[Header]'] == 'Description', 'Unnamed: 1'].values[0]
    df_RunSumm.loc[:, 'FlowcellID'] = FlowCellID
    df_RunSumm = df_RunSumm.drop(columns=['#'])
    st.session_state.df_RunSumm = df_RunSumm

    # Read the sample list file
    dfs_slf_samplelist = pd.read_excel(sample_list_file, sheet_name='Sample List 2.1', skiprows=4)
    filtered_df = dfs_slf_samplelist[dfs_slf_samplelist['Inostics ID'].isin(Samplelist)]
    filtered_df = filtered_df[['Inostics ID', 'External ID1\n(Patient ID-Visit)',
                               'External ID2\n(Collection datetime)', 'Scan External Barcode ']]
    filtered_df['External ID2\n(Collection datetime)'] = filtered_df['External ID2\n(Collection datetime)'].astype(str).str.replace(r'\.0$', '', regex=True)
    Sample_mapping = filtered_df.groupby("Inostics ID", as_index=False).first()
    sample_mapping_df = Sample_mapping.rename(columns={
        'Inostics ID': 'Inostics ID',
        'External ID1\n(Patient ID-Visit)': 'External ID1',
        'External ID2\n(Collection datetime)': 'External ID2',
        'Scan External Barcode ': 'External ID3'
    })
    st.session_state.sample_mapping_df = sample_mapping_df    

    # Display sample_mapping_df
#    display_dataframe(st.session_state.sample_mapping_df, "Sample List File Data")
    
    df_ResultsReview = df1_summary_tab.rename(columns={'Sample ID':'SampleID', 'Total DNA Amount (GE)':'GE', 'Amplicon ID':'AmpliconID', 'CDS Change':'CDSChange', 'AA Change':'AAChange', 'MAF [%]':'MAF[%]', 'Mutant Molecules':'MutantMolecules', 'COSMIC ID':'COSMICID', 'Base specific Cut-off':'Cutoff', 'Comment Call':'UpdateComment'})
    merged_df = df_ResultsReview.merge(sample_mapping_df, left_on='SampleId', right_on='Inostics ID', how='right')
    columns_to_drop = ['Raw Call', 'Reference Transcript ID', 'RUNID', 'SWVersion', 'Config File Name', 'userid', 'machineid', 'ChangeType', 'AvgMutantBaseQuality', 'CoverageStatus', '#positiveWells', 'LoB', 'LoQ', 'MutationHash', 'Mutation classification', 'timestamp', 'indexPlate', 'ampliconPosition', 'ML Plasma', 'MM/ML Plasma', 'MM corrected']
    merged_df = merged_df.drop(columns=columns_to_drop)
    merged_df2 = merged_df.merge(df_RunSumm, left_on='SampleId', right_on='Sample_ID', how='left')
    merged_df2 = merged_df2.rename(columns={'SampleId':'SampleID'})
    columns_to_drop2 = ['Inostics ID','External ID2', 'External ID3', 'Sample_ID', 'Plasma Vol. [mL]', 'FlowcellID']
    merged_df2 = merged_df2.drop(columns=columns_to_drop2)
    Order_list = ['FlowCellID', 'SampleID', 'Comment Call (Internal)', 'Comment Call (External)', 'Call', 'AmpliconID', 'GE', 'Gene Name', 'CDSChange', 'AAChange', 'MAF[%]', 'MutantMolecules', 'CosmicID', 'dbSNP', 'ClinVar', 'Raw Call', 'OOS', 'hg19Pos', '#UIDs/Amplicon', '#Supermutants', 'Cutoff', 'UpdateComment', 'Assay variant', 'Qubit Run ID', 'UID-PCR input (ng/116µl)', 'UID-PCR ID', 'UID-PCR wells', 'Index-PCR ID', 'NextSeq ID']
    df_ResultReview_import = merged_df2.reindex(columns=Order_list)    
#    st.session_state.df_ResultReview_import = df_ResultReview_import    
    st.success("Data combined successfully!")
    # Display ResultsReview imported sheet
    ###Add integrate code after this step if CHIP data is available
    if chipdatafile is not None:
        # Process the uploaded file
        df_ResultReview_import, merged_df_wBC ,rows_to_update, merged_df_nBC = chip_data_process(chipdatafile, df_ResultReview_import) 
        st.session_state.df_ResultReview_import = df_ResultReview_import
        st.session_state.rows_to_update = rows_to_update
        st.session_state.merged_df_nBC = merged_df_nBC
        st.session_state.merged_df_wBC = merged_df_wBC
#        display_dataframe(rows_to_update, "CHIP data Review")
        # Provide download link with the exact custom file name provided
        with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
            df1_summary_tab.to_excel(writer, sheet_name="Import_Raw Data", index=False)
            df_RunSumm.to_excel(writer, sheet_name="Import_Run Summary", index=False)
            sample_mapping_df.to_excel(writer, sheet_name="Sample Mapping", index=False)
            df_ResultReview_import.to_excel(writer, sheet_name="Results Review", index=False)
            # Add more sheets if needed
        st.success(f"File saved successfully to: {output_file_path}")
        
        # Optionally, provide the output in BytesIO for download
        output = BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            df1_summary_tab.to_excel(writer, sheet_name="Import_Raw Data", index=False)
            df_RunSumm.to_excel(writer, sheet_name="Import_Run Summary", index=False)
            sample_mapping_df.to_excel(writer, sheet_name="Sample Mapping", index=False)
            df_ResultReview_import.to_excel(writer, sheet_name="Results Review", index=False)
        output.seek(0)
        download_file_name = os.path.basename(output_file_path)
        st.download_button(
            label="Download Combined Output File with CHIP Data",
            data=output,
            file_name=download_file_name,
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    else:
        st.session_state.df_ResultReview_import = df_ResultReview_import
        with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
            df1_summary_tab.to_excel(writer, sheet_name="Import_Raw Data", index=False)
            df_RunSumm.to_excel(writer, sheet_name="Import_Run Summary", index=False)
            sample_mapping_df.to_excel(writer, sheet_name="Sample Mapping", index=False)
            df_ResultReview_import.to_excel(writer, sheet_name="Results Review", index=False)
            # Add more sheets if needed

        # Notify user of the saved file
        st.success(f"File saved successfully to: {output_file_path}")

        # Optionally, provide the output in BytesIO for download
        output = BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            df1_summary_tab.to_excel(writer, sheet_name="Import_Raw Data", index=False)
            df_RunSumm.to_excel(writer, sheet_name="Import_Run Summary", index=False)
            sample_mapping_df.to_excel(writer, sheet_name="Sample Mapping", index=False)
            df_ResultReview_import.to_excel(writer, sheet_name="Results Review", index=False)
        output.seek(0)
        # Provide download link with the exact custom file name provided
        download_file_name = os.path.basename(output_file_path)
        st.download_button(
            label="Download Combined Output File",
            data=output,
            file_name=download_file_name,
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )


if st.session_state.df1_summary_tab is not None:
    display_dataframe(
        st.session_state.df1_summary_tab,
        "SafeSeq Raw data Review"
    )

if st.session_state.df_RunSumm is not None:
    display_dataframe(
        st.session_state.df_RunSumm,
        "Running Summary data Review"
    )

if st.session_state.sample_mapping_df is not None:
    display_dataframe(
        st.session_state.sample_mapping_df,
        "Sample list data Review"
    )

if st.session_state.merged_df_wBC is not None:
    st.subheader("CHIP data Review Table which mutants have been detected in BC Experiment")
    st.dataframe(st.session_state.merged_df_wBC)


if st.session_state.rows_to_update is not None:
    st.subheader("CHIP data Review Table")
    st.dataframe(st.session_state.rows_to_update)

    st.subheader("CHIP data Visualization")
    plot_chip_data(st.session_state.rows_to_update
    )

if st.session_state.merged_df_nBC is not None:
    st.subheader("CHIP data Review Table which No detected in BC Experiment")
    st.dataframe(st.session_state.merged_df_nBC)

if st.session_state.df_ResultReview_import is not None:
    display_dataframe_with_filter(
        st.session_state.df_ResultReview_import,
        "Combined Results Review data"
    )
