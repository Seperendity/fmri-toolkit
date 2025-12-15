# -*- coding: utf-8 -*-
import numpy as np
import scipy
import scipy.io as sio
import matplotlib.pyplot as plt
import nibabel as nib
import os
import sys
import re
import argparse
from os.path import join, exists, split
import time
import urllib.request
import warnings
from tqdm import tqdm
from pprint import pprint
warnings.filterwarnings('ignore')
import matplotlib 
import logging
from glmsingle.glmsingle import GLM_single

def create_ses_run_dict(base_dir):
    """
    Create a dictionary mapping sessions to runs
    
    Args:
        base_dir: Base path to the sub-01 directory
        
    Returns:
        Dictionary in format {ses-{id}: [run-{id}, run-{id}, ...]}
    """
    # Initialize empty dictionary
    ses_run_dict = {}
    
    # Iterate through all items in base_dir
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        
        # Check if it's a directory starting with 'ses-'
        if os.path.isdir(item_path) and item.startswith('ses-'):
            ses_id = item
            run_list = []
            
            # Iterate through all items in current session directory
            for subitem in os.listdir(item_path):
                subitem_path = os.path.join(item_path, subitem)
                
                # Check if it's a directory starting with 'run-'
                if os.path.isdir(subitem_path) and subitem.startswith('run-'):
                    run_list.append(subitem)
            
            # Sort runs by number
            run_list.sort(key=lambda x: int(x.split('-')[1]))
            
            # Add to dictionary
            ses_run_dict[ses_id] = run_list
    
    # Sort sessions by number in dictionary
    ses_run_dict = dict(sorted(ses_run_dict.items(), key=lambda x: int(x[0].split('-')[1])))
    return ses_run_dict


def filter_sessions(ses_run_dict, target_sessions=None):
    """
    Filter sessions based on specified session IDs
    
    Args:
        ses_run_dict: Dictionary containing all sessions and runs
        target_sessions: List of session IDs to process (e.g., ['ses-01', 'ses-03'])
                        If None, process all sessions
    
    Returns:
        Filtered dictionary containing only specified sessions
    """
    if target_sessions is None:
        return ses_run_dict
    
    # Validate session IDs format
    valid_sessions = []
    for ses in target_sessions:
        if not ses.startswith('ses-'):
            # Try to convert numeric input to proper format
            if isinstance(ses, (int, str)) and str(ses).isdigit():
                ses = f'ses-{int(ses):02d}'
            else:
                print(f"Warning: Invalid session format '{ses}'. Expected format: 'ses-XX' or numeric")
                continue
        valid_sessions.append(ses)
    
    # Filter dictionary
    filtered_dict = {ses: runs for ses, runs in ses_run_dict.items() if ses in valid_sessions}
    
    # Check for missing sessions
    missing_sessions = set(valid_sessions) - set(filtered_dict.keys())
    if missing_sessions:
        print(f"Warning: The following sessions were not found: {list(missing_sessions)}")
    
    return filtered_dict


def parse_session_input(session_str):
    """
    Parse session input string into list of session IDs
    
    Args:
        session_str: String containing session IDs, e.g., "1,3,5" or "ses-01,ses-03"
        
    Returns:
        List of session IDs
    """
    if not session_str:
        return None
    
    # Split by comma and strip whitespace
    session_list = [s.strip() for s in session_str.split(',')]
    return session_list


def setup_logging(subject_id, log_dir='.'):
    """
    Setup logging configuration
    
    Args:
        subject_id: Subject identifier for log filename
        log_dir: Directory to save log files (fixed to current directory)
    """
    log_file = join(log_dir, f'glmsingle_{subject_id}_log.log')
    
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Also log to console
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


def run_glmsingle(design, data, ses_id, glmsingle_obj, outputdir_glmsingle, freeSurfer_subID, stimdur, tr_target):
    """
    Run GLMsingle analysis for a session
    """
    try:
        print(len(design), len(data), ses_id, stimdur, tr_target)
        # Run GLMsingle
        results_glmsingle = glmsingle_obj.fit(
            design,
            data,
            stimdur,
            tr_target,
            outputdir=join(outputdir_glmsingle, ses_id),
            figuredir=f'./glm_run_figures_fixed/{freeSurfer_subID}_ses{ses_id}/'
        )
        # Success log
        logging.info(f"ses-{ses_id} processing completed successfully")
        print("GLMsingle processing completed successfully")
        return results_glmsingle, True
    except Exception as e:
        # Failure log with error details
        logging.error(f"ses-{ses_id} processing failed: {str(e)}", exc_info=True)
        print(f"ses-{ses_id} processing failed: {str(e)}")
        return None, False


def main():
    """
    Main function with command line argument parsing
    """
    parser = argparse.ArgumentParser(
        description='Run GLMsingle analysis on fMRI data with session selection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Process all sessions
  python %(prog)s --subject sub-01
  
  # Process specific sessions (numeric format)
  python %(prog)s --subject sub-01 --sessions 1,3,5
  
  # Process specific sessions (ses- format)  
  python %(prog)s --subject sub-01 --sessions ses-01,ses-03,ses-05
  
  # Process single session
  python %(prog)s --subject sub-01 --sessions 1
        '''
    )
    
    # Required arguments
    parser.add_argument('--subject', '-s', default='sub-04',
                        help='Subject ID (e.g., sub-01)')
    
    # Optional arguments
    parser.add_argument('--sessions', '-ses', default='01,03',
                        help='Comma-separated list of session IDs to process (e.g., "1,3,5" or "ses-01,ses-03"). If not specified, all sessions will be processed.')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.subject, log_dir='.')
    logging.info(f"Starting GLMsingle processing for {args.subject}")
    logging.info(f"Command line arguments: {vars(args)}")
    
    # Parse session input
    target_sessions = parse_session_input(args.sessions)
    
    # Set fixed configuration (previously configurable parameters)
    datatype = 'our'
    stimdur = 3
    datalabel = 'func1pt8mm'
    
    # Set paths based on datatype
    preproc_dir = '/cy/data/preprocdata/'
    designmatrix_dir = '/cy/data/'
    tr_target = 1.5
    Effective_echo_spacing = 0.00027500
    origvoxelSize = np.array([2.0, 2.0, 2.0])
    
    # Set output directory
    Results_dir = "/cy/data/GLMsingle_outputs/"
    outputdir_glmsingle = join(Results_dir, args.subject, datalabel)
    
    # Get all sessions and filter if specified
    freeSurfer_subID = args.subject
    ses_run_dict_all = create_ses_run_dict(join(preproc_dir, freeSurfer_subID))
    ses_run_dict = filter_sessions(ses_run_dict_all, target_sessions)

    # CHECK FOR EXISTING OUTPUT FILES
    print("\nChecking for already completed sessions...")
    logging.info("Checking for already completed sessions to skip.")

    sessions_to_remove = []
    for ses_id in ses_run_dict:
        # Construct the path to the file that indicates a session is complete
        output_check_file = join(outputdir_glmsingle, ses_id, 'TYPED_FITHRF_GLMDENOISE_RR.npy')
        
        if os.path.exists(output_check_file):
            print(f"  - Found existing output for {ses_id}. This session will be skipped.")
            logging.info(f"Skipping {ses_id}: Output file found at {output_check_file}")
            sessions_to_remove.append(ses_id)

    # Remove the completed sessions from the dictionary
    for ses_id in sessions_to_remove:
        del ses_run_dict[ses_id]

    # Display processing information
    if target_sessions is None:
        print(f"Processing ALL sessions for {freeSurfer_subID}: {list(ses_run_dict.keys())}")
        logging.info(f"Processing ALL sessions: {list(ses_run_dict.keys())}")
    else:
        print(f"Processing SELECTED sessions for {freeSurfer_subID}: {list(ses_run_dict.keys())}")
        logging.info(f"Processing SELECTED sessions: {list(ses_run_dict.keys())}")
    
    if not ses_run_dict:
        print("No sessions to process. Exiting.")
        logging.error("No sessions to process")
        return
    
    # Setup GLMsingle options (fixed configuration)
    opt = dict()
    opt['wantlibrary'] = 1
    opt['wantglmdenoise'] = 1
    opt['wantfracridge'] = 1
    
    # For the purpose of this example we will keep the relevant outputs in memory
    # and also save them to the disk
    opt['wantfileoutputs'] = [0,1,0,1]
    opt['wantmemoryoutputs'] = [0,0,0,1]
    
    freesurfer_sub_id_reconname = args.subject  # For 'our' datatype, use subject ID as-is
    
    # Process each session
    total_sessions = len(ses_run_dict)
    processed_sessions = 0
    
    for ses_id, run_id_list in ses_run_dict.items():
        print(f"\nProcessing {freeSurfer_subID}, {ses_id} ({processed_sessions + 1}/{total_sessions})")
        logging.info(f"Starting processing for {freeSurfer_subID}, {ses_id}")
        
        # --- 1. Load current session's design matrix ---
        try:
            # Use simple split method to extract session number
            ses_id_num = ses_id.split('-')[1]
            design_file = join(designmatrix_dir, "design_matrix", freeSurfer_subID, f"session_{ses_id_num}.npy")
            
            # Load design matrix file containing all runs
            design_all_runs = np.load(design_file)  # Shape: (1, 12, 300, 583)
            
            # Extract design matrix for each run in current session
            design_per_ses = [design_all_runs[:, r].squeeze() for r in range(len(run_id_list))]
            print(f"Design matrix loaded: {len(design_per_ses)} runs, shape: {design_per_ses[0].shape}")
            logging.info(f"Successfully loaded design matrix for {ses_id}")
        except FileNotFoundError:
            print(f"Error: Design matrix file not found for {ses_id} at {design_file}")
            logging.error(f"Design matrix file not found for {ses_id} at {design_file}")
            continue # Skip this session
        except Exception as e:
            print(f"An error occurred while loading design matrix for {ses_id}: {e}")
            logging.error(f"Error loading design matrix for {ses_id}: {e}", exc_info=True)
            continue

        # --- 2. Load current session's functional data (fMRI Data) ---
        data_per_ses = []
        total_size_gb = 0
        try:
            for run_id_func in run_id_list:
                # Construct data file path
                func_path = join(preproc_dir, freesurfer_sub_id_reconname, ses_id, run_id_func, 'f_downsampled.nii.gz')
                
                print(f"  Loading data from: {func_path}")
                data_tmp = nib.load(func_path).get_fdata()
                
                # Special preprocessing for 'our' data type
                data_tmp = data_tmp[..., :-2]
                data_per_ses.append(data_tmp)
                
                # Calculate memory usage (optional)
                total_size_gb += data_tmp.nbytes / (1024**3)

            print(f"Successfully loaded {len(data_per_ses)} runs for {ses_id}.")
            print(f"Estimated memory for this session's data: {total_size_gb:.2f} GB")
            logging.info(f"Loaded data for {ses_id}, estimated size: {total_size_gb:.2f} GB")

        except FileNotFoundError as e:
            print(f"Error: Data file not found during processing of {ses_id}. Missing file: {e.filename}")
            logging.error(f"Data file not found for {ses_id}: {e.filename}", exc_info=True)
            continue # Skip this session
        except Exception as e:
            print(f"An error occurred while loading data for {ses_id}: {e}")
            logging.error(f"Error loading data for {ses_id}: {e}", exc_info=True)
            continue

        # --- 3. Run GLMsingle analysis ---
        print(f"Starting GLMsingle fit for {ses_id}...")
        start_time = time.time()

        # Create GLMsingle object
        glmsingle_obj = GLM_single(opt)
        
        # Log GLMsingle parameters
        logging.info(f"GLMsingle parameters: {glmsingle_obj.params}")
        print("GLMsingle parameters:")
        pprint(glmsingle_obj.params)
        
        run_glmsingle(design_per_ses, data_per_ses, ses_id, glmsingle_obj, 
                     outputdir_glmsingle, freeSurfer_subID, stimdur, tr_target)
        
        end_time = time.time()
        print(f"Finished GLMsingle fit for {ses_id} in {end_time - start_time:.2f} seconds.")
        logging.info(f"Finished GLMsingle fit for {ses_id} in {end_time - start_time:.2f} seconds.")

        # Free memory
        del data_per_ses
        del design_per_ses
        
        processed_sessions += 1

    print(f"\nAll sessions have been processed. ({processed_sessions}/{total_sessions} completed)")
    logging.info(f"All sessions have been processed. ({processed_sessions}/{total_sessions} completed)")


if __name__ == "__main__":
    main()