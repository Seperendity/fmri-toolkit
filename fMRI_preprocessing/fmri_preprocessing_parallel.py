"""
Parallel fMRI Preprocessing Pipeline
Author: wrx: Core module implementation, logic design
        cy: Algorithm acceleration, resource optimization
        lq: Code style compliance, security validation
Date: 2025-07-27
====================================

This script performs automated fMRI preprocessing in parallel batches, including slice timing correction,
motion correction, fieldmap-based distortion correction, and registration to anatomical space.

INPUT FILES:
-----------
Anatomical (per subject, ses-01 only):
  - {subject}_ses-01_T1w.nii.gz                    # Raw T1-weighted anatomical image
  - T1_brain.nii.gz                                # Skull-stripped T1 brain (from anatomical preprocessing)
  - T1_wmseg.nii.gz                                # White matter segmentation mask

Fieldmap (per session):
  - {subject}_{session}_magnitude1.nii.gz          # Fieldmap magnitude image  
  - fmap_mag_brain.nii.gz                          # Skull-stripped fieldmap magnitude
  - fmap_rads.nii.gz                               # Fieldmap in radians for distortion correction

Functional (per run):
  - {subject}_{session}_run-{XX}_bold.nii.gz       # Raw functional EPI timeseries

Timing Files:
  - tcustom_epi.txt                                 # Custom slice timing file

OUTPUT FILES:
------------
Motion Correction Output (per run):
  - f_skip_stc.nii.gz                              # Slice-timing corrected data
  - f_skip_stc_mc.nii.gz                           # Motion-corrected data
  - f_skip_stc_mc.nii.gz.par                       # Motion parameters file
  - f_skip_stc_mc.nii.gz_mean_reg.nii.gz          # Mean volume for registration

Distortion Correction Output (per run):
  - f_skip_stc_mc_fm_warp.nii.gz                   # Warp field for distortion correction
  - f_skip_stc_mc_fm.nii.gz                        # Final preprocessed data (registered to T1 space)
  - f_downsampled.nii.gz                           # Downsampled version at original resolution

Registration Output (shared per subject):
  - register.dof6.dat                              # BBR registration matrix (functional to T1)
  - init.register.dof6.dat                         # Initial registration matrix
  - register.downsampledfunc.dof6.dat              # Registration for downsampled data
  - T1_brain.registered.nii.gz                     # T1 brain registered to functional space
  - T1_global_mask.nii.gz                          # Global brain mask in T1 space
  - global_mask_func.nii.gz                        # Global brain mask in functional space

Process Control:
  - .preproc_success                                # Success flag indicating session completion

PROCESSING STEPS:
----------------
1. Slice timing correction using custom timing file
2. Motion correction to mean volume using MCFLIRT
3. Fieldmap-based distortion correction using epi_reg
4. Application of distortion correction to full timeseries
5. Boundary-based registration (BBR) to FreeSurfer anatomical surface
6. Generation of brain masks and registration matrices
7. Downsampling to original voxel resolution

PARALLEL PROCESSING:
-------------------
- Supports multi-worker parallel execution across subjects
- Processes multiple runs concurrently within each worker (MAX_CONCURRENT_RUNS = 6)
- Implements resume functionality using success flags to skip completed sessions
- Batched processing to manage memory and computational resources

USAGE:
------
# Single worker processing all subjects
python fmri_preprocessing_parallel.py

# Multi-worker parallel processing (2 workers total, this is worker 0)
python fmri_preprocessing_parallel.py --total-workers 2 --worker-id 0

# Multi-worker parallel processing (2 workers total, this is worker 1)  
python fmri_preprocessing_parallel.py --total-workers 2 --worker-id 1
"""

import os
import datetime
import random
import subprocess
import argparse
import multiprocessing
import glob
import numpy as np
import time

# Add FSL bin path to environment variables for FSL command execution
os.environ["FSLDIR"] = "/share/home/lq/software/fsl"
os.environ["PATH"] += os.pathsep + "/share/home/lq/software/fsl/bin"
# Set FreeSurfer installation directory for FreeSurfer command execution
os.environ["FREESURFER_HOME"] = "/share/home/lq/software/freesurfer"
os.environ["PATH"] += os.pathsep + os.path.join(os.environ["FREESURFER_HOME"], "bin")
# Set correct SUBJECTS_DIR (location where FreeSurfer looks for reconstruction data)
datatype = 'our' # 'our' here for 42
if datatype == 'nsd':
    os.environ["SUBJECTS_DIR"] = "/share/home/lq/data/NSD/nsddata/freesurfer/subjects"
else:
    os.environ["SUBJECTS_DIR"] = "/share/home/lq/software/freesurfer/subjects"

###-----------------------------------------
### T1.nii.gz is RAS, f.nii.gz-f_skip_stc_mc.nii.gz is LAS.
# after epi_reg, f_skip_stc_mc_fm.nii.gz is RAS, registered to T1.
###-----------------------------------------

def get_run_numbers_from_bids_structure(func_input_dir, freesurfer_sub_id, ses_id):
    pattern = os.path.join(func_input_dir, "*run-*_bold.nii.gz")
    matching_files = glob.glob(pattern)
    
    run_numbers = []
    for file_path in matching_files:
        filename = os.path.basename(file_path)
        import re
        match = re.search(r'run-(\d+)', filename)
        if match:
            # Save the matched string directly, preserving original format
            run_number_str = match.group(1)
            # Ensure two-digit format
            run_numbers.append(f"{int(run_number_str):02d}")

    if not run_numbers:
        print(f"No runs detected for {freesurfer_sub_id} {ses_id}")
        return []

    run_numbers = sorted(list(set(run_numbers)))
    print(f"Detected runs for {freesurfer_sub_id} {ses_id}: {run_numbers}")    
    return run_numbers

def create_batch_log_file(worker_id, total_workers):
    func_dir = os.path.dirname(os.path.realpath(__file__))
    func_filename = os.path.basename(__file__).split('.')[0]

    log_dir = os.path.join(func_dir, 'log')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    current_time = datetime.datetime.now()

    log_file_name = f"{func_filename}_batch_run_worker-{worker_id}-of-{total_workers}_{current_time.strftime('%Y%m%d-%H%M%S')}.log"
    log_file = os.path.join(log_dir, log_file_name)

    log_fid = open(log_file, 'w')
    log_fid.write(f"--- Batch Preprocessing Log (Worker {worker_id+1}/{total_workers}) ---\n")
    log_fid.write(f'Batch process started: {current_time.strftime("%Y-%m-%d %H:%M:%S")}\n\n')
    
    return log_fid

def process_run(run_args):
    # Unpack arguments
    run_id_func, freesurfer_sub_id, ses_id, base_dir, preproc_dir, \
    anat_T1_dir, anat_wmseg_dir, anat_brain_dir, field_input_dir, func_input_dir, \
    freesurfer_sub_id_reconname, TR_s, Effective_echo_spacing, origvoxelSize, tcustom_dir, \
    recon_anat_brain_dir, reg_output_dir, subspec_ROI_path, is_first_run = run_args

    log_messages = []

    def run_command_parallel(cmd):
        try:
            # Use capture_output to get stdout/stderr for better logging
            result = subprocess.run(cmd, check=True, shell=True, text=True, capture_output=True)
            log_messages.append(f"Command executed: {''.join(cmd)}\n")
            log_messages.append(f"STDOUT:\n{result.stdout}\n")
            if result.stderr:
                log_messages.append(f"STDERR:\n{result.stderr}\n")

        except subprocess.CalledProcessError as e:
            log_messages.append(f"---!!! ERROR executing command: {''.join(cmd)}\n")
            log_messages.append(f"Return Code: {e.returncode}\n")
            log_messages.append(f"STDOUT:\n{e.stdout}\n")
            log_messages.append(f"STDERR:\n{e.stderr}\n---!!!\n")

    filetopreproc_dir = os.path.join(func_input_dir, f"{freesurfer_sub_id}_{ses_id}_{run_id_func}_bold.nii.gz")
    field_mag_dir = os.path.join(field_input_dir, f"{freesurfer_sub_id}_{ses_id}_magnitude1.nii.gz")
    field_brain_dir = os.path.join(field_input_dir, 'field_output',  'fmap_mag_brain.nii.gz')
    field_fmap_dir = os.path.join(field_input_dir, 'field_output',  'fmap_rads.nii.gz')
    log_messages.append('\n Preprocessing ' + filetopreproc_dir + '\n')

    run_func_output_dir_name = os.path.join(preproc_dir, freesurfer_sub_id, ses_id, run_id_func)
    os.makedirs(run_func_output_dir_name, exist_ok=True)

    original_cwd = os.getcwd()
    os.chdir(run_func_output_dir_name)
    log_messages.append('\nChange to ' + run_func_output_dir_name + '.\n')

    mc_data_dir = os.path.join(run_func_output_dir_name, 'f_skip_stc_mc.nii.gz.par')
    func_fm_dir = os.path.join(run_func_output_dir_name, 'f_skip_stc_mc_fm.nii.gz')
    func_stc_mc_fm_1vol_dir = os.path.join(run_func_output_dir_name, 'f_skip_stc_mc_fm_1vol.nii.gz')
    func_fm_downsample_dir = os.path.join(run_func_output_dir_name, 'f_downsampled.nii.gz')
    
    anat_mask_global_dir = os.path.join(os.path.dirname(anat_brain_dir),  'T1_global_mask.nii.gz')
    registered_anat_brain_dir = os.path.join(os.path.dirname(anat_brain_dir),  'T1_brain.registered.nii.gz') ## output in RAS
    mask_global_dir_func = os.path.join(subspec_ROI_path, 'global_mask_func.nii.gz')
    register_dir = os.path.join(reg_output_dir, 'register.dof6.dat')
    init_register_dir = os.path.join(reg_output_dir, 'init.register.dof6.dat')
    register_dir_downsampledfunc = os.path.join(reg_output_dir, 'register.downsampledfunc.dof6.dat')
    init_register_dir_downsampledfunc = os.path.join(reg_output_dir, 'init.register.downsampledfunc.dof6.dat')

    if not os.path.exists('f_skip_stc_mc.nii.gz_mean_reg.nii.gz'):
        cmd = f"slicetimer -i {filetopreproc_dir} -o f_skip_stc.nii.gz -r {TR_s} --tcustom={tcustom_dir}"
        run_command_parallel(cmd)


        # Register timeseries to mean volume
        cmd = "mcflirt -in f_skip_stc.nii.gz -out f_skip_stc_mc.nii.gz -cost normmi -meanvol -plots"
        run_command_parallel(cmd)

    if not os.path.exists(func_fm_dir):
        if os.path.exists(field_fmap_dir) and not os.path.exists('f_skip_stc_mc_fm_warp.nii.gz'):
            cmd = f"epi_reg --echospacing={Effective_echo_spacing} --wmseg={anat_wmseg_dir} \
                --fmap={field_fmap_dir} \
                --fmapmag={field_mag_dir} \
                --fmapmagbrain={field_brain_dir} \
                --pedir=-y \
                --epi=f_skip_stc_mc.nii.gz_mean_reg.nii.gz \
                --t1={anat_T1_dir} \
                --t1brain={anat_brain_dir} \
                --out=f_skip_stc_mc_fm"
            run_command_parallel(cmd)


        if os.path.exists('f_skip_stc_mc_fm_warp.nii.gz'):
            cmd = f"applywarp -i f_skip_stc_mc.nii.gz \
                    -r {anat_T1_dir} \
                    -w f_skip_stc_mc_fm_warp.nii.gz \
                    --interp=spline \
                    -o f_skip_stc_mc_fm.nii.gz"
            run_command_parallel(cmd)
        log_messages.append('\n 1st level preprocessing for task fMRI data finished \n')
    else:
        log_messages.append('\n ' + func_fm_dir + ' exists...\n')

    if not os.path.exists(func_fm_downsample_dir):
        cmd = f"mri_convert --voxsize {origvoxelSize[0]} {origvoxelSize[1]} {origvoxelSize[2]} f_skip_stc_mc_fm.nii.gz {func_fm_downsample_dir}"
        run_command_parallel(cmd)

    if is_first_run:
        log_messages.append(f'\n[{run_id_func}] This is the first run of the session, generating shared files...\n')

        if not os.path.exists(register_dir):
            cmd = f"bbregister --s {freesurfer_sub_id_reconname} --init-fsl --6 --bold --mov {func_stc_mc_fm_1vol_dir} --reg {register_dir} --init-reg-out {init_register_dir} --frame 0"
            run_command_parallel(cmd)
        else:
            log_messages.append('\n ' + register_dir + 'for ' + func_stc_mc_fm_1vol_dir + ' exists...\n')

        if not os.path.exists(mask_global_dir_func):
            cmd = f"mri_vol2vol --mov {recon_anat_brain_dir} --regheader \
                    --targ {func_stc_mc_fm_1vol_dir} \
                    --o {registered_anat_brain_dir} --no-save-reg"
            run_command_parallel(cmd)
            cmd = f"mri_binarize --i {registered_anat_brain_dir} --o {anat_mask_global_dir} --min 0.1"
            run_command_parallel(cmd)
            cmd = f"mri_convert --voxsize {origvoxelSize[0]} {origvoxelSize[1]} {origvoxelSize[2]} {anat_mask_global_dir} {mask_global_dir_func}"
            run_command_parallel(cmd)
        else:
            log_messages.append(f'\n[{run_id_func}] ' + mask_global_dir_func + ' exists...\n')

        if not os.path.exists(register_dir_downsampledfunc):
            cmd = f"bbregister --s {freesurfer_sub_id_reconname} --init-fsl --6 --bold --mov {func_fm_downsample_dir} --reg {register_dir_downsampledfunc} --init-reg-out {init_register_dir_downsampledfunc} --frame 0"
            run_command_parallel(cmd)
        else:
            log_messages.append(f'\n[{run_id_func}] ' + register_dir_downsampledfunc + 'for downsampled func:' + func_fm_downsample_dir + ' exists...\n')
    else:
        log_messages.append(f'\n[{run_id_func}] Not the first run, skipping shared files generation.\n')

    if os.path.exists('f_skip_stc_mc.nii.gz'):
        try:
            os.remove('f_skip_stc_mc.nii.gz')
            os.remove('f_skip_stc.nii.gz')
        except OSError as e:
            log_messages.append(f"Error removing temp files: {e}\n")
    
    os.chdir(original_cwd)
    return log_messages


# --- SCRIPT CONFIGURATION ---
freesurfer_recon_dir = os.environ.get('SUBJECTS_DIR')

if datatype == 'nsd':
    base_dir = '/share/home/lq/data/NSD/nsddata_rawdata'
    preproc_dir = '/share/home/lq/data/NSD/nsddata_preprocdata_multiprocessing/'

    TR_s = 1.6 # nsd: 1.6; our data: 1.5
    Effective_echo_spacing = 0.000329994 # nsd: 0.000329994; our data: 0.000275003
    origvoxelSize = np.array([1.8, 1.8, 1.8])
elif datatype == 'our':
    base_dir = '/share/home/lq/data/our_nsd/bids'
    preproc_dir = '/share/home/lq/data/our_nsd/preprocdata/'

    TR_s = 1.5
    Effective_echo_spacing = 0.000275003
    origvoxelSize = np.array([2.0, 2.0, 2.0])

tcustom_dir = os.path.join(base_dir, 'tcustom_epi.txt')

MAX_CONCURRENT_RUNS = 10

# Define success flag filename
SUCCESS_FLAG_FILENAME = ".preproc_success"

# -------------------------------------------------------------------------------------
# --- Main Processing Logic (with resume functionality) ---
# -------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Batch fMRI preprocessing with parallel worker support.")
parser.add_argument('--total-workers', type=int, default=2, help='The total number of parallel worker scripts being run.')
parser.add_argument('--worker-id', type=int, default=0, help='The 0-indexed ID of this specific worker instance.')
args = parser.parse_args()

total_workers = args.total_workers
worker_id = args.worker_id

if worker_id >= total_workers:
    print(f"Error: worker-id ({worker_id}) must be less than total-workers ({total_workers}).")
    exit(1)

start_time = datetime.datetime.now()

print(f"--- Starting Batch Preprocessing [Worker {worker_id+1} of {total_workers}] ---")
print("Scanning for all subjects to determine workload...")

all_subject_dirs = sorted(glob.glob(os.path.join(base_dir, 'sub-*')))
if not all_subject_dirs:
    print(f"Error: No subject directories found in '{base_dir}'. Please check the path.")
    exit()

# Divide subject data among workers
subject_splits = np.array_split(all_subject_dirs, total_workers)

subject_dirs_for_this_worker = subject_splits[worker_id]

if len(subject_dirs_for_this_worker) == 0:
    print(f"Worker {worker_id+1} has no subjects to process. Exiting.")
    exit(0)

print(f"\nWorker {worker_id+1} is assigned the following {len(subject_dirs_for_this_worker)} subjects:")
for sub_path in subject_dirs_for_this_worker:
    print(f"- {os.path.basename(sub_path)}")

all_sessions_to_process = []
# Build work queue
for sub_path in subject_dirs_for_this_worker:
    subject_id = os.path.basename(sub_path)
    session_dirs = sorted(glob.glob(os.path.join(sub_path, 'ses-*')))
    for ses_path in session_dirs:
        session_id = os.path.basename(ses_path)

        # Resume functionality check
        session_preproc_dir = os.path.join(preproc_dir, subject_id, session_id)
        success_flag_path = os.path.join(session_preproc_dir, SUCCESS_FLAG_FILENAME)

        if os.path.exists(success_flag_path):
            print(f"--> Skipping [{subject_id}/{session_id}]: Already processed.")
            continue

        all_sessions_to_process.append((subject_id, session_id))

if not all_sessions_to_process:
    print("\nAll assigned sessions have already been processed. Nothing to do. Exiting.")
    exit()

print(f"\nFound {len(all_sessions_to_process)} unprocessed sessions for this worker.")
print("--- Work Queue for this Worker ---")
for sub, ses in all_sessions_to_process:
    print(f"- {sub}, {ses}")
print("----------------------------------\n")

log_fid = create_batch_log_file(worker_id, total_workers)
log_fid.write("Discovered Work Queue for this worker:\n")
for sub, ses in all_sessions_to_process:
    log_fid.write(f"- {sub}, {ses}\n")
log_fid.write("----------------------------------\n\n")
log_fid.flush()

# Collect all runs from all sessions into batches
print("Collecting all runs from unprocessed sessions...")
all_runs_to_process = []

for freesurfer_sub_id, ses_id in all_sessions_to_process:
    freesurfer_sub_id_reconname = freesurfer_sub_id.replace("sub-", "subj") if datatype == 'nsd' else freesurfer_sub_id.replace("sub-", "sub")
    
    # Set up directories for this session
    anat_input_dir = os.path.join(base_dir, freesurfer_sub_id, 'ses-01', 'anat')
    anat_output_dir = os.path.join(anat_input_dir, 'anat_output')
    reg_output_dir = os.path.join(preproc_dir, freesurfer_sub_id, 'reg_output')    
    subspec_ROI_path = os.path.join(preproc_dir, freesurfer_sub_id, 'ROIs')
    
    # Create necessary directories
    folders = [anat_output_dir, reg_output_dir, subspec_ROI_path]
    for folder in folders:
        os.makedirs(folder, exist_ok=True)
    
    # Set up file paths
    anat_T1_dir = os.path.join(anat_input_dir, f"{freesurfer_sub_id}_ses-01_T1w.nii.gz")
    anat_brain_dir = os.path.join(anat_output_dir, 'T1_brain.nii.gz')
    anat_wmseg_dir = os.path.join(anat_output_dir, 'T1_wmseg.nii.gz')
    recon_anat_brain_dir = os.path.join(freesurfer_recon_dir, freesurfer_sub_id_reconname, 'mri', 'brainmask.mgz')
    field_input_dir = os.path.join(base_dir, freesurfer_sub_id, ses_id, 'fmap')
    func_input_dir = os.path.join(base_dir, freesurfer_sub_id, ses_id, 'func')

    run_numbers = get_run_numbers_from_bids_structure(func_input_dir, freesurfer_sub_id, ses_id)

    if not run_numbers:
        print(f"No runs found for {freesurfer_sub_id} {ses_id}, skipping...")
        continue

    # Sort run numbers to ensure consistent ordering
    run_numbers = sorted(run_numbers)

    for run_number in run_numbers:
        run_id_func = f'run-{run_number}'

        # is_first_run is True only for the first run of the first session per subject
        is_first_run = True if ses_id == 'ses-01' and run_number == '01' else False
        
        filetopreproc_dir = os.path.join(func_input_dir, f"{freesurfer_sub_id}_{ses_id}_{run_id_func}_bold.nii.gz")
        if not os.path.exists(filetopreproc_dir):
            print(f"Warning: File {filetopreproc_dir} does not exist, skipping run {run_id_func}")
            continue
        
        run_args = (
            run_id_func, freesurfer_sub_id, ses_id, base_dir, preproc_dir,
            anat_T1_dir, anat_wmseg_dir, anat_brain_dir, field_input_dir, func_input_dir,
            freesurfer_sub_id_reconname, TR_s, Effective_echo_spacing, origvoxelSize, tcustom_dir,
            recon_anat_brain_dir, reg_output_dir, subspec_ROI_path, is_first_run
        )
        all_runs_to_process.append(run_args)

if not all_runs_to_process:
    print("No runs to process. Exiting.")
    log_fid.write("No runs to process. Exiting.\n")
    log_fid.close()
    exit()

print(f"Total runs collected: {len(all_runs_to_process)}")
log_fid.write(f"Total runs collected: {len(all_runs_to_process)}\n")

# Process runs in batches of MAX_CONCURRENT_RUNS
batch_size = MAX_CONCURRENT_RUNS
total_batches = (len(all_runs_to_process) + batch_size - 1) // batch_size
batch_success_status = []

print(f"Processing {len(all_runs_to_process)} runs in {total_batches} batches of up to {batch_size} runs each")
log_fid.write(f"Processing {len(all_runs_to_process)} runs in {total_batches} batches of up to {batch_size} runs each\n")
log_fid.flush()

for batch_idx in range(total_batches):
    batch_start_time = datetime.datetime.now()
    start_idx = batch_idx * batch_size
    end_idx = min(start_idx + batch_size, len(all_runs_to_process))
    batch_runs = all_runs_to_process[start_idx:end_idx]
    
    batch_message = f"--- Processing Batch {batch_idx + 1}/{total_batches}: {len(batch_runs)} runs ---"
    print(batch_message)
    log_fid.write(f"\n{batch_message}\n")
    log_fid.write(f"Batch start time: {batch_start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Log which runs are in this batch
    for run_args in batch_runs:
        run_id_func, freesurfer_sub_id, ses_id = run_args[0], run_args[1], run_args[2]
        log_fid.write(f"  - {freesurfer_sub_id}/{ses_id}/{run_id_func}\n")
    log_fid.flush()
    
    batch_successful = True
    try:
        with multiprocessing.Pool(processes=len(batch_runs)) as pool:
            results = pool.map(process_run, batch_runs)
        
        # Write all log messages from this batch
        for log_messages in results:
            for message in log_messages:
                log_fid.write(message)
    
    except Exception as e:
        batch_successful = False
        error_message = f"---!!! CRITICAL ERROR during batch {batch_idx + 1}. !!!---\nError: {e}\n"
        print(error_message)
        log_fid.write(error_message)

    # Store batch success status
    batch_success_status.append(batch_successful)
    
    batch_end_time = datetime.datetime.now()
    batch_duration = batch_end_time - batch_start_time
    
    if batch_successful:
        success_message = f"--- Finished batch {batch_idx + 1}/{total_batches} in {batch_duration} ---"
        print(success_message)
        log_fid.write(f"\n{success_message}\n")
    else:
        failure_message = f"--- FAILED batch {batch_idx + 1}/{total_batches}. See logs for details. ---"
        print(failure_message)
        log_fid.write(f"\n{failure_message}\n")
    
    log_fid.flush()
    
    # Add small delay between batches to prevent resource conflicts
    time.sleep(2)

# After all batches, create success flags for completed sessions
print("Creating success flags for completed sessions...")
log_fid.write("\nCreating success flags for completed sessions...\n")

# Track which runs were successfully processed
successfully_processed_runs = set()
# Use the batch_success_status that was already populated during batch processing
for batch_idx, batch_successful in enumerate(batch_success_status):
    if batch_successful:
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, len(all_runs_to_process))
        batch_runs = all_runs_to_process[start_idx:end_idx]
        
        for run_args in batch_runs:
            run_id_func, freesurfer_sub_id, ses_id = run_args[0], run_args[1], run_args[2]
            successfully_processed_runs.add((freesurfer_sub_id, ses_id, run_id_func))

# Group runs by session and check completion
session_run_counts = {}
session_successful_counts = {}

for run_args in all_runs_to_process:
    run_id_func, freesurfer_sub_id, ses_id = run_args[0], run_args[1], run_args[2]
    session_key = (freesurfer_sub_id, ses_id)
    
    session_run_counts[session_key] = session_run_counts.get(session_key, 0) + 1
    
    if (freesurfer_sub_id, ses_id, run_id_func) in successfully_processed_runs:
        session_successful_counts[session_key] = session_successful_counts.get(session_key, 0) + 1

# Create success flags only for sessions where all runs were successful
for session_key, total_runs in session_run_counts.items():
    freesurfer_sub_id, ses_id = session_key
    successful_runs = session_successful_counts.get(session_key, 0)
    
    if successful_runs == total_runs:
        session_preproc_dir = os.path.join(preproc_dir, freesurfer_sub_id, ses_id)
        os.makedirs(session_preproc_dir, exist_ok=True)
        success_flag_path = os.path.join(session_preproc_dir, SUCCESS_FLAG_FILENAME)
        
        with open(success_flag_path, 'w') as f:
            f.write(f"Completed on: {datetime.datetime.now().isoformat()}\n")
            f.write(f"Total runs processed: {successful_runs}/{total_runs}\n")
        
        log_fid.write(f"Created success flag: {success_flag_path} ({successful_runs}/{total_runs} runs)\n")
        print(f"Session {freesurfer_sub_id}/{ses_id}: All {total_runs} runs completed successfully")
    else:
        log_fid.write(f"Session {freesurfer_sub_id}/{ses_id}: Only {successful_runs}/{total_runs} runs completed - no success flag created\n")
        print(f"Session {freesurfer_sub_id}/{ses_id}: Only {successful_runs}/{total_runs} runs completed - will retry on next run")

end_time = datetime.datetime.now()
total_execution_time = end_time - start_time
final_message = f"\n--- Batch job finished or no more sessions to process. ---"
time_message = f"Total execution time (including wait times): {total_execution_time}"

print(final_message)
print(time_message)
log_fid.write(f"\n{final_message}\n")
log_fid.write(f"{time_message}\n")
log_fid.write(f"Batch process finished: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
log_fid.close()

# Usage: Run script specifying total workers as 2, this is worker 0 (each processes 2 subjects)
# python your_script_name.py --total-workers 2 --worker-id 0
