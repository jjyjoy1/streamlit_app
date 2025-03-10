import streamlit as st
import pandas as pd
import os
import subprocess
import tempfile
import re
import threading
import time

st.set_page_config(page_title="BCL2FASTQ Demultiplexer", layout="wide")

def validate_sample_sheet(content):
    """
    Validates the sample sheet format.
    Returns (is_valid, error_message)
    """
    # Check if content is empty
    if not content:
        return False, "Sample sheet is empty"
    
    # Basic format checks
    lines = content.strip().split('\n')
    
    # Check for header section
    if not any(line.startswith('[Header]') for line in lines):
        return False, "Missing [Header] section"
    
    # Check for Data section
    if not any(line.startswith('[Data]') for line in lines):
        return False, "Missing [Data] section"
    
    # Extract the data section for validation
    in_data_section = False
    data_lines = []
    
    for line in lines:
        if line.startswith('[Data]'):
            in_data_section = True
            continue
        if in_data_section and line.strip() and not line.startswith('['):
            data_lines.append(line)
    
    # Check if we have data rows
    if len(data_lines) < 2:  # At least one header line and one data line
        return False, "No sample data found in [Data] section"
    
    # Check for required columns in the data section
    header = data_lines[0].split(',')
    required_columns = ['Sample_ID', 'Sample_Name', 'index']
    
    for col in required_columns:
        if col not in header:
            return False, f"Missing required column: {col}"
    
    # Check for duplicate Sample_IDs
    try:
        # Parse data section into a dataframe for validation
        df = pd.read_csv(pd.StringIO('\n'.join(data_lines)))
        if df['Sample_ID'].duplicated().any():
            duplicate_ids = df[df['Sample_ID'].duplicated()]['Sample_ID'].unique()
            return False, f"Duplicate Sample_IDs found: {', '.join(duplicate_ids)}"
    except Exception as e:
        return False, f"Error parsing data section: {str(e)}"
    
    # All checks passed
    return True, "Sample sheet is valid"

def parse_sample_sheet_to_dataframe(content):
    """
    Parses a sample sheet content into a pandas DataFrame for display
    """
    # Locate the [Data] section
    lines = content.strip().split('\n')
    data_start = None
    
    for i, line in enumerate(lines):
        if line.startswith('[Data]'):
            data_start = i + 1
            break
    
    if data_start is None or data_start >= len(lines):
        return pd.DataFrame()
    
    # Extract data section and parse as CSV
    data_content = '\n'.join(lines[data_start:])
    return pd.read_csv(pd.StringIO(data_content))

def run_bcl2fastq(input_dir, output_dir, sample_sheet_path, threads=4):
    """
    Runs bcl2fastq with the provided parameters
    Returns (success, output_message)
    """
    cmd = [
        "bcl2fastq",
        "--input-dir", input_dir,
        "--output-dir", output_dir,
        "--sample-sheet", sample_sheet_path,
        "--no-lane-splitting",
        "--processing-threads", str(threads)
    ]
    
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            return False, f"bcl2fastq failed with error code {process.returncode}:\n{stderr}"
        
        return True, f"bcl2fastq completed successfully.\nOutput:\n{stdout}"
    except Exception as e:
        return False, f"Failed to run bcl2fastq: {str(e)}"

def background_process(input_dir, output_dir, sample_sheet_path, status_key):
    """Background thread for running bcl2fastq"""
    try:
        success, message = run_bcl2fastq(input_dir, output_dir, sample_sheet_path)
        st.session_state[status_key] = {"complete": True, "success": success, "message": message}
    except Exception as e:
        st.session_state[status_key] = {"complete": True, "success": False, "message": str(e)}

# Initialize session state variables if they don't exist
if 'job_status' not in st.session_state:
    st.session_state.job_status = {"running": False, "complete": False, "success": False, "message": ""}

# App title and description
st.title("BCL2FASTQ Demultiplexer")
st.markdown("""
This app allows you to review NGS sample sheets and submit demultiplexing jobs using bcl2fastq.
""")

# Layout with two columns
col1, col2 = st.columns([1, 1])

with col1:
    # Input directory selection
    st.header("1. Select NGS Data")
    input_dir = st.text_input("BCL Input Directory", 
                             placeholder="Path to the NGS data folder containing BaseCalls")
    
    if input_dir and os.path.exists(input_dir):
        st.success(f"✅ Directory exists: {input_dir}")
        # Show directory contents
        files = os.listdir(input_dir)
        st.write(f"Directory contains {len(files)} files/folders")
        with st.expander("View directory contents"):
            st.write(files)
    elif input_dir:
        st.error(f"❌ Directory does not exist: {input_dir}")
    
    # Output directory
    st.header("2. Select Output Location")
    output_dir = st.text_input("Output Directory", 
                              placeholder="Path where FASTQ files will be saved")
    
    if output_dir:
        if not os.path.exists(output_dir):
            if st.button("Create Output Directory"):
                try:
                    os.makedirs(output_dir, exist_ok=True)
                    st.success(f"✅ Created output directory: {output_dir}")
                except Exception as e:
                    st.error(f"❌ Failed to create directory: {str(e)}")
        else:
            st.success(f"✅ Output directory exists: {output_dir}")

    # Sample sheet upload or input
    st.header("3. Sample Sheet")
    sample_sheet_option = st.radio(
        "How would you like to provide the sample sheet?",
        ["Upload file", "Paste content"]
    )
    
    sample_sheet_content = None
    
    if sample_sheet_option == "Upload file":
        uploaded_file = st.file_uploader("Upload Sample Sheet", type=["csv"])
        if uploaded_file is not None:
            sample_sheet_content = uploaded_file.getvalue().decode("utf-8")
            st.success("✅ Sample sheet uploaded")
    else:
        sample_sheet_content = st.text_area(
            "Paste Sample Sheet Content",
            height=300,
            placeholder="Paste the contents of your sample sheet here...",
            help="The sample sheet should include [Header] and [Data] sections with required columns"
        )

with col2:
    # Sample sheet validation and preview
    st.header("4. Sample Sheet Validation")
    
    if sample_sheet_content:
        is_valid, message = validate_sample_sheet(sample_sheet_content)
        
        if is_valid:
            st.success(f"✅ {message}")
            
            # Display sample data in a table
            df = parse_sample_sheet_to_dataframe(sample_sheet_content)
            if not df.empty:
                st.subheader("Sample Data Preview")
                st.dataframe(df)
                
                # Sample stats
                st.info(f"Total samples: {len(df)}")
                
                # Show unique indexes
                with st.expander("Index Distribution"):
                    if 'index' in df.columns:
                        index_counts = df['index'].value_counts()
                        st.write(index_counts)
                        
                        # Check for potential index collisions
                        if len(index_counts) < len(df):
                            st.warning("⚠️ Some indexes are used multiple times!")
        else:
            st.error(f"❌ {message}")
    
    # Job submission section
    st.header("5. Job Submission")
    
    # Requirements check
    requirements_met = input_dir and os.path.exists(input_dir) and \
                      output_dir and \
                      sample_sheet_content and validate_sample_sheet(sample_sheet_content)[0]
    
    if not requirements_met:
        st.warning("Please complete all previous steps before submitting the job.")
    
    # Threads selection
    threads = st.slider("Number of processing threads", min_value=1, max_value=32, value=4)
    
    # Submit button
    submit_button = st.button("Submit Demultiplexing Job", disabled=not requirements_met)
    
    if submit_button:
        # Create a temporary file for the sample sheet
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as temp_file:
            temp_file.write(sample_sheet_content)
            sample_sheet_path = temp_file.name
        
        # Set status and start background process
        st.session_state.job_status = {
            "running": True, 
            "complete": False, 
            "success": False,
            "message": "Job started...",
            "sample_sheet_path": sample_sheet_path
        }
        
        # Start background thread
        bg_thread = threading.Thread(
            target=background_process,
            args=(input_dir, output_dir, sample_sheet_path, 'job_status')
        )
        bg_thread.start()
        
        st.experimental_rerun()

# Display job status
if st.session_state.job_status["running"]:
    st.header("Job Status")
    
    if not st.session_state.job_status["complete"]:
        st.info("⏳ Job is running... This might take several minutes.")
        
        # Add a progress bar that doesn't actually track real progress
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # A placeholder for moving progress 
        if not hasattr(st.session_state, 'progress_value'):
            st.session_state.progress_value = 0
            
        # Increment progress for visual feedback
        st.session_state.progress_value = (st.session_state.progress_value + 1) % 101
        progress_bar.progress(st.session_state.progress_value)
        status_text.text(f"Processing... {st.session_state.progress_value}%")
        
        # Check again in 0.5 seconds
        time.sleep(0.5)
        st.experimental_rerun()
    else:
        # Job completed
        if st.session_state.job_status["success"]:
            st.success("✅ Job completed successfully!")
        else:
            st.error("❌ Job failed!")
        
        # Show output message
        with st.expander("View job output"):
            st.code(st.session_state.job_status["message"], language="bash")
            
        # Cleanup the temporary sample sheet file
        if "sample_sheet_path" in st.session_state.job_status:
            try:
                os.unlink(st.session_state.job_status["sample_sheet_path"])
            except:
                pass
                
        # Reset button
        if st.button("Start New Job"):
            st.session_state.job_status = {"running": False, "complete": False, "success": False, "message": ""}
            st.experimental_rerun()

