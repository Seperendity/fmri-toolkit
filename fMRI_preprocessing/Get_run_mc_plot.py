#!/usr/bin/env python3
import numpy as np
import os
import matplotlib.pyplot as plt
from os.path import join, exists
import logging
import sys
import glob
from datetime import datetime


def setup_logging(log_dir):
    """Set up logging to both file and console"""
    os.makedirs(log_dir, exist_ok=True)

    # Create log filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = join(log_dir, f'run_motion_QC_{timestamp}.log')

    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

    return log_file

def log_print(message):
    """Print and log message"""
    logging.info(message)

def find_sessions_and_runs(base_dir, subID):
    """
    Automatically search for all sessions and runs for a specified subject.

    Args:
        base_dir: Base data directory
        subID: Subject ID, e.g., 'sub-01'

    Returns:
        dict: {session_id: [run_list]}
    """
    subject_dir = join(base_dir, subID)

    if not exists(subject_dir):
        log_print(f"Subject directory not found: {subject_dir}")
        return {}

    # Search for all session directories
    session_pattern = join(subject_dir, "ses-*")
    session_dirs = glob.glob(session_pattern)

    sessions_runs = {}

    for session_dir in sorted(session_dirs):
        session_id = os.path.basename(session_dir)

        # Search for all run directories under this session
        run_pattern = join(session_dir, "run-*")
        run_dirs = glob.glob(run_pattern)

        runs = []
        for run_dir in sorted(run_dirs):
            run_id = os.path.basename(run_dir)

            # Check if motion parameter file exists
            motion_file = join(run_dir, 'f_skip_stc_mc.nii.gz.par')
            if exists(motion_file):
                runs.append(run_id)
            else:
                log_print(f"Motion file not found: {motion_file}")

        if runs:
            sessions_runs[session_id] = runs
            log_print(f"Found {subID} {session_id}: {len(runs)} runs - {runs}")

    return sessions_runs

def create_single_run_plot(mc_data, subID, session_id, run_id, run_index, output_dir):
    """Create motion parameter plot for a single run."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot rotations
    motion_labels_rot = ['rot_x', 'rot_y', 'rot_z']
    for i, label in enumerate(motion_labels_rot):
        axes[0].plot(mc_data[:, i], label=label, alpha=0.7)
    axes[0].set_title(f"{subID} {session_id} {run_id} - Rotations (rad)")
    axes[0].set_xlabel('Time (TR)')
    axes[0].set_ylabel('Rotation (rad)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot translations
    motion_labels_trans = ['trans_x', 'trans_y', 'trans_z']
    for i, label in enumerate(motion_labels_trans):
        axes[1].plot(mc_data[:, i+3], label=label, alpha=0.7)
    axes[1].set_title(f"{subID} {session_id} {run_id} - Translations (mm)")
    axes[1].set_xlabel('Time (TR)')
    axes[1].set_ylabel('Translation (mm)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(join(output_dir, f'{subID}_{session_id}_{run_id}_motion_index{run_index}.png'),
               dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    """Main analysis function."""
    # Configuration parameters
    base_dir = '/cy/data/preprocdata/'
    subID_list = [f'sub-{i:02d}' for i in range(1, 5)]  # sub-01 to sub-04
    hrf_latency_nTR = 4  # 6s

    # Set up logging
    log_dir = join(base_dir,  'mc_QC', 'logs')
    log_file = setup_logging(log_dir)

    log_print("Starting automatic run motion QC analysis...")
    log_print(f"Log file: {log_file}")
    log_print(f"Base directory: {base_dir}")
    log_print(f"Subjects to process: {subID_list}")

    # Store data for all subjects
    all_subjects_data = {}

    # Process each subject
    for subID in subID_list:
        log_print(f"\n{'='*60}")
        log_print(f"Processing subject: {subID}")
        log_print(f"{'='*60}")

        # Automatically search for sessions and runs
        sessions_runs = find_sessions_and_runs(base_dir, subID)

        if not sessions_runs:
            log_print(f"No sessions/runs found for {subID}")
            continue

        # Create output directory
        subject_output_dir = join(base_dir, 'mc_QC', subID)
        os.makedirs(subject_output_dir, exist_ok=True)

        global_run_index = 0

        log_print(f"Processing subject {subID}")

        for session_id in sorted(sessions_runs.keys()):
            runs = sessions_runs[session_id]
            log_print(f"  Processing {session_id} with {len(runs)} runs")

            for run_id in runs:
                motion_file = join(base_dir, subID, session_id, run_id, 'f_skip_stc_mc.nii.gz.par')

                try:
                    # Load motion parameters
                    mc_data = np.loadtxt(motion_file)
                    mc_data_currentrun = mc_data[hrf_latency_nTR:]

                    # Create motion parameter plot for single run
                    create_single_run_plot(mc_data_currentrun,  subID, session_id, run_id,
                                        global_run_index, subject_output_dir)

                    global_run_index += 1

                except Exception as e:
                    log_print(f"    Error processing {session_id}/{run_id}: {e}")
                    run_key = f"{session_id}_{run_id}"

                    global_run_index += 1

        log_print(f"Plots saved to: {subject_output_dir}")

    # Execute QC analysis
    log_print(f"\n{'='*60}")
    log_print("Starting QC Analysis")
    log_print(f"{'='*60}")


    log_print("\n" + "="*50)
    log_print("ANALYSIS COMPLETE")
    log_print("="*50)


if __name__ == "__main__":
    results = main()
