import io
import os
import pandas as pd
import streamlit as st
import warnings

warnings.filterwarnings("ignore")  # Suppress warnings

st.title("SafeSaq Data Processor")

sample_list_file = st.file_uploader("Upload Sample List File", type=["xlsx"])
result_review_file = st.file_uploader("Upload Result Review File", type=["xlsx"])

def display_dataframe(df, title):
    st.subheader(title)
    st.dataframe(df)

###Input data table
if "df_SampleData" not in st.session_state:
    st.session_state.df_SampleData = None

if "df_RR_RawData" not in st.session_state:
    st.session_state.df_RR_RawData = None

###Output data table
if "df_final_sample_inf" not in st.session_state:
    st.session_state.df_final_sample_inf = None

if "df_final_variants" not in st.session_state:
    st.session_state.df_final_variants = None

if "df_final_genes" not in st.session_state:
    st.session_state.df_final_genes = None

if "df_final_mutant_summary" not in st.session_state:
    st.session_state.df_final_mutant_summary = None


if sample_list_file and result_review_file:
    st.success("Files uploaded successfully!")

    # Display file names
    st.write("**Sample List File:**", sample_list_file.name)
    st.write("**Result Review File:**", result_review_file.name)

    # Read data
    df_SampleData = pd.read_excel(sample_list_file, sheet_name='SampleDataFile')
    df_RR_Sample_Information = pd.read_excel(result_review_file, sheet_name='Sample Information')
    df_RR_RawData = pd.read_excel(result_review_file, sheet_name='RawData')
    df_SampleData['ReportDate'] = df_SampleData['ReportDate'].astype(str)
    tmp_merged_table = df_RR_Sample_Information

    st.session_state.df_SampleData = df_SampleData
    st.session_state.df_RR_RawData = df_RR_RawData

    output_file_path_sample = st.text_input("Enter the full path for the Sample Information File", "/Users/user_defined_path/Sample_output.xlsx")

    # Button to combine and process
    if st.button("Combine and Process Sample Information"):
        # Prepare final sample info (as you did previously)
        df_final_sample_inf = df_RR_Sample_Information.rename(columns={
            'Sample ID': 'SAMPID',
            'External ID1': 'SUBJID',
            'External ID2': 'SPECID',
            'External ID3': 'SPECID2',
            'Volume (mL)': 'VOLUME'
        })

#        st.subheader("Final Sample Information")
#        st.dataframe(df_final_sample_inf)
        st.session_state.df_final_sample_inf = df_final_sample_inf

        # Try to save to the specified local path exactly as in the working code
        if not df_final_sample_inf.empty:
            try:
                # Write to the specified output file
                with pd.ExcelWriter(output_file_path_sample, engine='openpyxl') as writer:
                    df_final_sample_inf.to_excel(writer, sheet_name="Sample_Info", index=False)
            
                st.success(f"File saved successfully to: {output_file_path_sample}")

            except Exception as e:
                st.error(f"Failed to save the file locally: {e}")

####
    output_file_path_variant = st.text_input("Enter the full path for the Variants Summary Information File", "/Users/user_defined_path/Variant_output.xlsx")
    if st.button("Combine and Process Variants Information"):
        df_RR_Variants = df_RR_RawData[df_RR_RawData['Call'] == 'MD']
        df_merged = pd.merge(
            df_SampleData,
            df_RR_Variants,
            left_on='InosticsID',   # key in df_SampleData
            right_on='Sample ID', # key in df_RR_Variants
            how='right'           # or 'left', 'right', 'outer' depending on your needs
            )
        df_merged = df_merged.rename(columns={
            'Study': 'PROTOCOL',
            'Visit': 'VISIT',
            'Collection Date': 'CTDNADT',
            'Collection Time': 'CTDNATM',
            'SampleID': 'SUBJID',
            'InosticsID': 'SAMPID',
            'Total DNA Amount (GE)': 'TOTDNAMT',
            'Gene Name': 'GNNAME',
            'COSMIC ID': 'Cosmic ID',
            'MAF [%]': 'MAF[%]',
            'Mutant Molecules': 'MUTMOL',
            'Call': 'CALL',
            'Base specific Cut-off': 'Base Specific Cut-off',
            '#UIDs/Amplicon': 'UIDAMP',
            '#Supermutants': 'SUPMUT',
            'Comment Call': 'COMMCALL'
             # If you want to standardize Amplicon ID, CDS Change, AA Change, ClinVar, dbSNP
                # just leave them as is or rename them accordingly
        })
        final_columns = [
            'PROTOCOL', 'VISIT', 'CTDNADT', 'CTDNATM', 'SAMPID', 'SUBJID',
            'TOTDNAMT', 'GNNAME', 'Amplicon ID', 'CDS Change', 'AA Change',
            'Cosmic ID', 'ClinVar', 'dbSNP', 'MAF[%]', 'MUTMOL', 'CALL',
            'Base Specific Cut-off', 'UIDAMP', 'SUPMUT', 'COMMCALL'
            ]
        df_final_variants = df_merged[final_columns]
        def variant_data_format_modification(df):
            df['PROTOCOL'] = df['PROTOCOL'].apply(lambda x: x.split()[0])
            df['CTDNADT'] = pd.to_datetime(df['CTDNADT'], errors='coerce').dt.strftime('%d %b %Y')
            df['CTDNATM'] = df['CTDNATM'].astype(str).str.replace(r'(\d{2}:\d{2}):\d{2}', r'\1', regex=True)
            df['SUPMUT'] = df['SUPMUT'].astype(int).astype(str)
            df = df.astype(str)
            return df
        df_final_variants = variant_data_format_modification(df_final_variants)
        st.session_state.df_final_variants = df_final_variants      

        if not df_final_variants.empty:
            try:
                # Write to the specified output file
                with pd.ExcelWriter(output_file_path_variant, engine='openpyxl') as writer:
                    df_final_variants.to_excel(writer, sheet_name="Sample_Info", index=False)

                st.success(f"File saved successfully to: {output_file_path_variant}")

            except Exception as e:
                st.error(f"Failed to save the file locally: {e}")

####Process Gene Information
    output_file_path_gene = st.text_input("Enter the full path for the Gene Summary Information File", "/Users/user_defined_path/Gene_output.xlsx")
    if st.button("Combine and Process Genes Information"):
        gene_list = list(set(df_RR_RawData['Gene Name'].to_list()))
        gene_subtables = {}
        for gene in gene_list:
            # Filter for 'Gene Name' == gene and 'Call' == "MD"
            gene_subtables[gene] = df_RR_RawData[
            (df_RR_RawData['Gene Name'] == gene) & (df_RR_RawData['Call'] == "MD")
            ]
            # Select specific columns
            # Check if the table is empty
            if gene_subtables[gene].empty:
                gene_subtables[gene][gene] = "NMD"
            else:
                # Add a column 'Status' with the gene name for non-empty tables
                gene_subtables[gene][gene] = 'MD'
                gene_subtables[gene] = gene_subtables[gene][['Sample ID', gene]]

            tmp_merged_table = pd.merge(tmp_merged_table,
                    gene_subtables[gene][['Sample ID', gene]],
                    on='Sample ID',
                    how='left')
            tmp_merged_table[gene] = tmp_merged_table[gene].fillna('NMD')

        df_total_ge = df_RR_RawData[["Sample ID", 'Total DNA Amount (GE)']].groupby('Sample ID').first().reset_index()
        tmp_merged_table['Overall Status'] = tmp_merged_table[gene_list].apply(
            lambda row: 'MD' if 'MD' in row.values else 'NMD',
            axis=1
            )
        tmp_merged_table = pd.merge(tmp_merged_table,
            df_total_ge,
            on='Sample ID',
            how='left')
        column_order = ['Sample ID', 'External ID1', 'External ID2', 'External ID3',
           'Total DNA Amount (GE)', 'Overall Status'] + gene_list
        df_gene_summary = tmp_merged_table[column_order]
        df_merged_gene_summary = pd.merge(
            df_SampleData,
            df_gene_summary,
            left_on='InosticsID',   # key in df_SampleData
            right_on='Sample ID', # key in df_RR_Variants
            how='right'           # or 'left', 'right', 'outer' depending on your needs
            )
        df_merged_gene_summary = df_merged_gene_summary.rename(columns={
            'Study': 'PROTOCOL',
            'Visit': 'VISIT',
            'Collection Date': 'CTDNADT',
            'Collection Time': 'CTDNATM',
            'SampleID': 'SUBJID',
            'InosticsID': 'SAMPID',
            'Total DNA Amount (GE)': 'TOTDNAMT',
            'Gene Name': 'GNNAME',
            'COSMIC ID': 'Cosmic ID',
            'MAF [%]': 'MAF[%]',
            'Mutant Molecules': 'MUTMOL',
            'Call': 'CALL',
            'Base specific Cut-off': 'Base Specific Cut-off',
            '#UIDs/Amplicon': 'UIDAMP',
            '#Supermutants': 'SUPMUT',
            'External ID3': 'SPECID',
            'External ID2': 'SPECID2',
            'Sample Comment': 'COMMENT',
            'Overall Status':'STATUS',
            # If you want to standardize Amplicon ID, CDS Change, AA Change, ClinVar, dbSNP
            # just leave them as is or rename them accordingly
             })
        final_columns_gene = ['PROTOCOL', 'VISIT', 'CTDNADT', 'CTDNATM', 'SAMPID', 'SUBJID', 'SPECID',
                 'SPECID2', 'TOTDNAMT', 'STATUS' ] + gene_list + ['COMMENT']
        df_final_genes = df_merged_gene_summary[final_columns_gene]
        def gene_data_format_modification(df):
            df['PROTOCOL'] = df['PROTOCOL'].apply(lambda x: x.split()[0])
            df['CTDNATM'] = df['CTDNATM'].astype(str).str.replace(r'(\d{2}:\d{2}):\d{2}', r'\1', regex=True)
            df['CTDNADT'] = pd.to_datetime(df['CTDNADT'], errors='coerce').dt.strftime('%d %b %Y')
            df = df.astype(str)
        df_final_genes = gene_data_format_modification(df_final_genes)
        df_final_genes = df_merged_gene_summary[final_columns_gene]
        df_final_genes['CTDNADT'] = df_final_genes['CTDNADT'].astype("string")
        df_final_genes['CTDNATM'] = df_final_genes['CTDNATM'].astype("string")

        st.session_state.df_final_genes = df_final_genes
         
        if not df_final_genes.empty:
            try:
                # Write to the specified output file
                with pd.ExcelWriter(output_file_path_gene, engine='openpyxl') as writer:
                    df_final_genes.to_excel(writer, sheet_name="Sample_Info", index=False)
                st.success(f"File saved successfully to: {output_file_path_gene}")
            except Exception as e:
                st.error(f"Failed to save the file locally: {e}")
   
####Process Mutants Summary Information
    
    output_file_path_mutant = st.text_input("Enter the full path for the Mutants Summary Information File", "/Users/user_defined_path/mutant_output.xlsx")
    if st.button("Combine and Process Mutants Information"):
        gene_list = list(set(df_RR_RawData['Gene Name'].to_list()))
        gene_subtables = {}
        for gene in gene_list:
            # Filter for 'Gene Name' == gene and 'Call' == "MD"
            gene_subtables[gene] = df_RR_RawData[
                (df_RR_RawData['Gene Name'] == gene) & (df_RR_RawData['Call'] == "MD")
            ]
            # Select specific columns
            # Check if the table is empty
            if gene_subtables[gene].empty:
                gene_subtables[gene][gene] = "NMD"
            else:
                # Add a column 'Status' with the gene name for non-empty tables
                gene_subtables[gene][gene] = 'MD'
                gene_subtables[gene] = gene_subtables[gene][['Sample ID', gene]]

            tmp_merged_table = pd.merge(tmp_merged_table,
                    gene_subtables[gene][['Sample ID', gene]],
                    on='Sample ID',
                    how='left')
            tmp_merged_table[gene] = tmp_merged_table[gene].fillna('NMD')

        df_total_ge = df_RR_RawData[["Sample ID", 'Total DNA Amount (GE)']].groupby('Sample ID').first().reset_index()
        tmp_merged_table['Overall Status'] = tmp_merged_table[gene_list].apply(
            lambda row: 'MD' if 'MD' in row.values else 'NMD',
            axis=1
            )
        tmp_merged_table = pd.merge(tmp_merged_table,
            df_total_ge,
            on='Sample ID',
            how='left')
        def get_description(row, df_RR_RawData):
            # Genes to check
            gene_list = list(set(df_RR_RawData['Gene Name'].to_list()))
            genes = gene_list
            # List to store mutation descriptions
            descriptions = []
            # Check each gene with MD status
            for gene in genes:
                try:  
                    if row[gene] == 'MD':
                    # Filter the raw data based on Sample ID, Call, and Gene Name
                        mutant_df = df_RR_RawData[
                            (df_RR_RawData['Sample ID'] == row['Sample ID']) &
                            (df_RR_RawData['Call'] == 'MD') &
                            (df_RR_RawData['Gene Name'] == gene)
                            ]
                        # If mutations found, format each mutation
                        if not mutant_df.empty:
                            for _, mut_row in mutant_df.iterrows():
                                description = f"{mut_row['Gene Name']} | {mut_row['CDS Change']} | {mut_row['MAF [%]']} | {mut_row['Mutant Molecules']}"
                                descriptions.append(description)
                        # Join multiple descriptions with newline if multiple exist
                except Exception as e:
                    pass
            return '\n'.join(descriptions) if descriptions else 'NMD'
                # Apply the function to create the Description column
        tmp_merged_table['Description   | Gene Name | CDS Change | MAF | MM'] = tmp_merged_table.apply(
                lambda row: get_description(row, df_RR_RawData),
                axis=1
                )
        df_mutant_summary = tmp_merged_table
        df_merged_mutant_summary = pd.merge(
            df_SampleData,
            df_mutant_summary,
            left_on='InosticsID',   # key in df_SampleData
            right_on='Sample ID', # key in df_RR_Variants
            how='right'           # or 'left', 'right', 'outer' depending on your needs
            )
        df_merged_mutant_summary_rename = df_merged_mutant_summary.rename(columns={
            'Study': 'PROTOCOL',
            'Visit': 'VISIT',
            'Collection Date': 'CTDNADT',
            'Collection Time': 'CTDNATM',
            'InosticsID': 'SAMPID',
            'External ID1':'SUBJID',
            'External ID2':'SPECID',
            'External ID3':'SPECID2',
            'Overall Status':'STATUS',
            'Total DNA Amount (GE)':'TOTDNAMT',
            'Description   | Gene Name | CDS Change | MAF | MM':'DESCRP',
            })
        final_order = ['PROTOCOL', 'VISIT', 'CTDNADT', 'CTDNATM', 'SAMPID', 'SUBJID', 'SPECID', 'SPECID2', 'TOTDNAMT', 'STATUS', 'DESCRP']
        df_mutant_summary_rename_reorder = df_merged_mutant_summary_rename[final_order]
        def mutant_data_format_modification(df):
            df['PROTOCOL'] = df['PROTOCOL'].apply(lambda x: x.split()[0])
            df['CTDNATM'] = df['CTDNATM'].astype(str).str.replace(r'(\d{2}:\d{2}):\d{2}', r'\1', regex=True)
            df['CTDNADT'] = pd.to_datetime(df['CTDNADT'], errors='coerce').dt.strftime('%d %b %Y')
            return df
        df_final_mutant_summary = mutant_data_format_modification(df_mutant_summary_rename_reorder)

        st.session_state.df_final_mutant_summary = df_final_mutant_summary
        
        if not df_final_mutant_summary.empty:
            try:
                # Write to the specified output file
                with pd.ExcelWriter(output_file_path_mutant, engine='openpyxl') as writer:
                    df_final_mutant_summary.to_excel(writer, sheet_name="Mutant_Info", index=False)
                st.success(f"File saved successfully to: {output_file_path_mutant}")
            except Exception as e:
                st.error(f"Failed to save the file locally: {e}")

else:
    st.warning("Please upload both files to proceed.")

###Input
if st.session_state.df_SampleData is not None:
    display_dataframe(
        st.session_state.df_SampleData,
        "Input Sample List Data Review"
    )

if st.session_state.df_RR_RawData is not None:
    display_dataframe(
        st.session_state.df_RR_RawData,
        "Input Results Review Data Review"
    )

###Output
if st.session_state.df_final_sample_inf is not None:
    display_dataframe(
        st.session_state.df_final_sample_inf,
        "Output Sample List Data Information Review"
    )

if st.session_state.df_final_variants is not None:
    display_dataframe(
        st.session_state.df_final_variants,
        "Output Variants Summary Data Information Review"
    )

if st.session_state.df_final_genes is not None:
    display_dataframe(
        st.session_state.df_final_genes,
        "Output Gene Summary Data Information Review"
    )

if st.session_state.df_final_mutant_summary is not None:
    display_dataframe(
        st.session_state.df_final_mutant_summary,
        "Output Mutants Summary Data Information Review"
    )
