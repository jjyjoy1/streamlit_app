import streamlit as st
import pandas as pd
import numpy as np
import os
import tempfile
import subprocess
import json

st.title("R ExpressionSet Object Viewer")

# File uploader
uploaded_file = st.file_uploader("Choose an RDS file", type="rds")

if uploaded_file is not None:
    # Save the uploaded file temporarily
    with tempfile.NamedTemporaryFile(suffix='.rds', delete=False) as f:
        temp_filename = f.name
        f.write(uploaded_file.getvalue())
    
    # Create an R script to extract information from the RDS file
    r_script = """
    library(Biobase)
    library(jsonlite)
    
    args <- commandArgs(trailingOnly = TRUE)
    rds_file <- args[1]
    output_dir <- args[2]
    
    # Read the RDS file
    eset <- readRDS(rds_file)
    
    # Extract expression data
    expr_data <- t(exprs(eset))
    write.csv(expr_data, file.path(output_dir, "expression_data.csv"))
    
    # Extract sample metadata
    sample_meta <- pData(eset)
    write.csv(sample_meta, file.path(output_dir, "sample_metadata.csv"))
    
    # Extract feature metadata if available
    if (exists("fData") && is.function(fData)) {
        feature_meta <- fData(eset)
        write.csv(feature_meta, file.path(output_dir, "feature_metadata.csv"))
    } else {
        # Create an empty file to indicate no feature metadata
        file.create(file.path(output_dir, "no_feature_metadata.txt"))
    }
    
    # Extract basic information
    basic_info <- list(
        dataset_name = basename(rds_file),
        num_samples = ncol(exprs(eset)),
        num_features = nrow(exprs(eset)),
        sample_names = sampleNames(eset),
        feature_names = featureNames(eset)
    )
    
    # Save as JSON
    write_json(basic_info, file.path(output_dir, "basic_info.json"))
    
    # If 'Class' column exists in sample metadata, extract it
    if ("Class" %in% colnames(sample_meta)) {
        classes <- sample_meta$Class
        write.csv(data.frame(Class=classes), file.path(output_dir, "classes.csv"))
    }
    
    # If 'Status' and 'Survival_in_days' columns exist in sample metadata, extract them
    if (all(c("Status", "Survival_in_days") %in% colnames(sample_meta))) {
        surv_data <- sample_meta[, c("Status", "Survival_in_days")]
        write.csv(surv_data, file.path(output_dir, "survival_data.csv"))
    }
    """
    
    # Create a temporary directory to store output files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Save the R script to a file
        r_script_path = os.path.join(temp_dir, "extract_data.R")
        with open(r_script_path, "w") as f:
            f.write(r_script)
        
        try:
            # Run the R script
            subprocess.run(["Rscript", r_script_path, temp_filename, temp_dir], check=True)
            
            # Load the extracted data
            basic_info_path = os.path.join(temp_dir, "basic_info.json")
            if os.path.exists(basic_info_path):
                with open(basic_info_path, "r") as f:
                    basic_info = json.load(f)
            else:
                st.error("Failed to extract basic information from the RDS file.")
                if os.path.exists(temp_filename):
                    os.unlink(temp_filename)
                st.stop()
            
            # Load expression data
            expr_data_path = os.path.join(temp_dir, "expression_data.csv")
            if os.path.exists(expr_data_path):
                X = pd.read_csv(expr_data_path, index_col=0)
            else:
                st.error("Failed to extract expression data from the RDS file.")
                if os.path.exists(temp_filename):
                    os.unlink(temp_filename)
                st.stop()
            
            # Load sample metadata
            sample_meta_path = os.path.join(temp_dir, "sample_metadata.csv")
            if os.path.exists(sample_meta_path):
                sample_meta = pd.read_csv(sample_meta_path, index_col=0)
            else:
                st.warning("No sample metadata found in the RDS file.")
                sample_meta = pd.DataFrame(index=X.index)
            
            # Load feature metadata if available
            feature_meta_path = os.path.join(temp_dir, "feature_metadata.csv")
            if os.path.exists(feature_meta_path):
                feature_meta = pd.read_csv(feature_meta_path, index_col=0)
            else:
                st.info("No feature metadata found in the RDS file.")
                feature_meta = pd.DataFrame(index=X.columns)
            
            # Load class information if available
            classes_path = os.path.join(temp_dir, "classes.csv")
            if os.path.exists(classes_path):
                classes = pd.read_csv(classes_path)["Class"].values
            else:
                classes = None
            
            # Load survival data if available
            survival_path = os.path.join(temp_dir, "survival_data.csv")
            if os.path.exists(survival_path):
                survival_data = pd.read_csv(survival_path)
            else:
                survival_data = None
            
            # Display the extracted information
            st.header(f"Dataset: {basic_info['dataset_name']}")
            st.write(f"Number of samples: {basic_info['num_samples']}")
            st.write(f"Number of features: {basic_info['num_features']}")
            
            # Create tabs for different components
            tab1, tab2, tab3, tab4 = st.tabs(["Expression Data", "Sample Metadata", 
                                            "Feature Metadata", "Classes/Survival"])
            
            with tab1:
                st.subheader("Expression Data")
                st.write(f"Shape: {X.shape} - {X.shape[0]} samples × {X.shape[1]} features")
                st.dataframe(X.head())
                
                # Summary statistics
                st.subheader("Summary Statistics")
                st.dataframe(X.describe())
            
            with tab2:
                st.subheader("Sample Metadata")
                st.dataframe(sample_meta)
                
                # Display column info
                if not sample_meta.empty and len(sample_meta.columns) > 0:
                    st.subheader("Sample Metadata Columns")
                    for col in sample_meta.columns:
                        st.write(f"**{col}**: {sample_meta[col].dtype}")
                        if sample_meta[col].dtype == 'object' or pd.api.types.is_categorical_dtype(sample_meta[col]):
                            st.write(sample_meta[col].value_counts())
                        else:
                            st.write(sample_meta[col].describe())
            
            with tab3:
                st.subheader("Feature Metadata")
                st.dataframe(feature_meta)
            
            with tab4:
                if classes is not None:
                    st.subheader("Classes")
                    classes_df = pd.DataFrame({"Class": classes})
                    st.write("Class Distribution:")
                    st.write(classes_df["Class"].value_counts())
                    st.bar_chart(classes_df["Class"].value_counts())
                
                if survival_data is not None:
                    st.subheader("Survival Data")
                    st.dataframe(survival_data.head(20))
                    st.write("Summary statistics of survival time:")
                    st.write(survival_data["Survival_in_days"].describe())
            
        except subprocess.CalledProcessError as e:
            st.error(f"Error running R script: {str(e)}")
        except Exception as e:
            st.error(f"Error processing data: {str(e)}")
    
    # Clean up the temporary file
    if os.path.exists(temp_filename):
        os.unlink(temp_filename)
else:
    st.info("Please upload an RDS file to view its contents.")



