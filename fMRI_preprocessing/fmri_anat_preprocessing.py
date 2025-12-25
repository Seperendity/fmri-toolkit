#!/usr/bin/env python
# coding: utf-8
"""
Anatomical and fieldmap preprocessing script.

Performs automated preprocessing including anatomical skull-stripping,
tissue segmentation, and fieldmap correction for multiple sessions.

See README.md for detailed documentation on input/output files and usage.
"""

import os
import datetime
import numpy as np
import nibabel as nib
import scipy
import scipy.signal as ss
import scipy.io as io
from matplotlib import pyplot as plt
import argparse
import sys
import re

from glob import glob
from nilearn.plotting import plot_anat, plot_roi
from nilearn import plotting

# Add FSL bin path to environment variables for calling FSL commands
os.environ["FSLDIR"] = "/share/home/lq/software/fsl/"
os.environ["PATH"] += os.pathsep + "/share/home/lq/software/fsl/bin"

# Set base directory
base_dir = '/share/home/lq/data/our_nsd/bids'


def get_available_sessions(subject_id):
    """Automatically detect all available sessions for a subject."""
    subject_dir = os.path.join(base_dir, subject_id)

    if not os.path.exists(subject_dir):
        print(f"Error: Subject directory does not exist: {subject_dir}")
        return []

    # Find all directories matching ses-xx format
    ses_pattern = os.path.join(subject_dir, 'ses-*')
    ses_dirs = glob(ses_pattern)

    # Extract session numbers and sort
    sessions = []
    for ses_dir in ses_dirs:
        if os.path.isdir(ses_dir):
            ses_name = os.path.basename(ses_dir)
            # Verify if it matches ses-xx format
            if re.match(r'ses-\d+', ses_name):
                sessions.append(ses_name)

    sessions.sort()  # Sort alphabetically, so ses-01 comes before ses-02

    print(f"Detected {len(sessions)} sessions: {sessions}")
    return sessions


def setup_anat_directories(freesurfer_sub_id):
    """Set up directories and file paths for anatomical images."""
    ## set up input folders
    anat_input_dir = os.path.join(base_dir, freesurfer_sub_id, 'ses-01', 'anat')
    anat_output_dir = os.path.join(anat_input_dir, 'anat_output')
    reg_output_dir = os.path.join(anat_input_dir, 'reg_output')

    ## check and create folders
    folders = [anat_input_dir, anat_output_dir, reg_output_dir]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder, exist_ok=True)

    ## set up input and output filenames
    anat_T1_dir = os.path.join(anat_input_dir, f"{freesurfer_sub_id}_ses-01_T1w.nii.gz")
    anat_brain_dir = os.path.join(anat_output_dir, 'T1_brain.nii.gz')
    anat_wmseg_dir = os.path.join(anat_output_dir, 'T1_wmseg.nii.gz')
    anat_pve0_dir = os.path.join(anat_output_dir, 'T1_brain_pve_0.nii.gz')
    anat_pve1_dir = os.path.join(anat_output_dir, 'T1_brain_pve_1.nii.gz')
    anat_pve2_dir = os.path.join(anat_output_dir, 'T1_brain_pve_2.nii.gz')
    anat_gm_mask_dir = os.path.join(anat_output_dir, "T1_gm_mask.nii.gz")

    tcustom_dir = os.path.join(base_dir, 'tcustom_epi.txt')

    return {
        'anat_input_dir': anat_input_dir,
        'anat_output_dir': anat_output_dir,
        'reg_output_dir': reg_output_dir,
        'anat_T1_dir': anat_T1_dir,
        'anat_brain_dir': anat_brain_dir,
        'anat_wmseg_dir': anat_wmseg_dir,
        'anat_pve0_dir': anat_pve0_dir,
        'anat_pve1_dir': anat_pve1_dir,
        'anat_pve2_dir': anat_pve2_dir,
        'anat_gm_mask_dir': anat_gm_mask_dir,
        'tcustom_dir': tcustom_dir,
    }

def setup_field_directories(freesurfer_sub_id, ses_id):
    """Set up directories and file paths for field map processing."""
    ## set up input folders
    field_input_dir = os.path.join(base_dir, freesurfer_sub_id, ses_id, 'fmap')

    ## set up output folders
    field_output_dir = os.path.join(field_input_dir, 'field_output')
    field_QC_fig_dir = os.path.join(field_input_dir, 'field_output', 'QC')

    ## check and create folders
    folders = [field_input_dir, field_output_dir, field_QC_fig_dir]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder, exist_ok=True)

    ## set up input and output filenames
    field_mag_dir = os.path.join(field_input_dir, f"{freesurfer_sub_id}_{ses_id}_magnitude1.nii.gz")
    field_pha_dir = os.path.join(field_input_dir, f"{freesurfer_sub_id}_{ses_id}_phasediff.nii.gz")
    field_brain1_dir = os.path.join(field_output_dir, 'fmap_mag_brain1.nii.gz')
    field_brain_dir = os.path.join(field_output_dir, 'fmap_mag_brain.nii.gz')
    field_fmap_dir = os.path.join(field_output_dir, 'fmap_rads.nii.gz')

    return {
        'field_input_dir': field_input_dir,
        'field_output_dir': field_output_dir,
        'field_QC_fig_dir': field_QC_fig_dir,
        'field_mag_dir': field_mag_dir,
        'field_pha_dir': field_pha_dir,
        'field_brain1_dir': field_brain1_dir,
        'field_brain_dir': field_brain_dir,
        'field_fmap_dir': field_fmap_dir
    }


def plot_brain_slices(img_path, output_dir, prefix):
    """Plot brain image slices and save to file."""
    # Define slice coordinates
    x_slices = [-40, -20, 0, 20, 40, 60]
    y_slices = [-60, -40, -20, 0, 20, 40]
    z_slices = [10, 30, 40, 50, 60, 70]

    plot_img = nib.load(img_path)

    # X-axis slices
    display_x = plot_anat(
        plot_img,
        display_mode='x',
        cut_coords=x_slices,
        title='X Slices',
        draw_cross=False,
        annotate=True,
        bg_img=None,
        black_bg=False
    )
    display_x.savefig(os.path.join(output_dir, f"{prefix}_x_slices.png"), dpi=500)
    display_x.close()
    print(f"Saved {prefix} x_slices.png")

    # Y-axis slices
    display_y = plot_anat(
        plot_img,
        display_mode='y',
        cut_coords=y_slices,
        title='Y Slices',
        draw_cross=False,
        annotate=True,
        bg_img=None,
        black_bg=False
    )
    display_y.savefig(os.path.join(output_dir, f"{prefix}_y_slices.png"), dpi=500)
    display_y.close()
    print(f"Saved {prefix} y_slices.png")

    # Z-axis slices
    display_z = plot_anat(
        plot_img,
        display_mode='z',
        cut_coords=z_slices,
        title='Z Slices',
        draw_cross=False,
        annotate=True,
        bg_img=None,
        black_bg=False
    )
    display_z.savefig(os.path.join(output_dir, f"{prefix}_z_slices.png"), dpi=500)
    display_z.close()
    print(f"Saved {prefix} z_slices.png")

def check_outputs_exist(paths):
    """Check if key output files already exist."""
    key_outputs = [
        paths['field_brain_dir'],
        paths['field_fmap_dir']
    ]
    return all(os.path.exists(f) for f in key_outputs)

def process_fieldmap(freesurfer_sub_id, ses_id):
    """Process fieldmap for a single session."""
    print(f"\n=== Processing fieldmap for {ses_id} ===")

    # Set up directories and file paths
    paths = setup_field_directories(freesurfer_sub_id, ses_id)

    print(f"Field magnitude: {paths['field_mag_dir']}")
    print(f"Field phase: {paths['field_pha_dir']}")

    # Check if input files exist
    if not os.path.exists(paths['field_mag_dir']):
        print(f"Warning: {paths['field_mag_dir']} does not exist, skipping this session")
        return False
    if not os.path.exists(paths['field_pha_dir']):
        print(f"Warning: {paths['field_pha_dir']} does not exist, skipping this session")
        return False

    # Check if output files already exist
    if check_outputs_exist(paths):
        print(f"Output files already exist, skipping {ses_id}")
        return True

    ## 1. Skull-strip the magnitude nifti image
    print("1. Performing skull-stripping on magnitude image...")
    cmd = f"bet2 {paths['field_mag_dir']} {paths['field_brain1_dir']} -f 0.5"
    os.system(cmd)

    cmd = f"fslmaths {paths['field_brain1_dir']} -ero {paths['field_brain_dir']}"
    os.system(cmd)

    ## 2. Plot brain slices for QC
    print("2. Generating QC images for brain...")
    plot_brain_slices(paths['field_brain_dir'], paths['field_QC_fig_dir'], "field_brain")

    ## 3. Create fmap_rads.nii.gz file
    print("3. Creating fieldmap...")
    cmd = f"fsl_prepare_fieldmap SIEMENS {paths['field_pha_dir']} {paths['field_brain_dir']} {paths['field_fmap_dir']} 2.46"
    os.system(cmd)

    ## 4. Plot fieldmap slices for QC
    print("4. Generating QC images for fieldmap...")
    plot_brain_slices(paths['field_fmap_dir'], paths['field_QC_fig_dir'], "field_fmap")

    print(f"=== {ses_id} fieldmap processing complete ===\n")
    return True

def process_anatomical(freesurfer_sub_id):
    """Process anatomical images (only runs once)."""
    print("\n=== Processing anatomical images ===")

    # Set up directories and file paths
    paths = setup_anat_directories(freesurfer_sub_id)

    # Check if input file exists
    if not os.path.exists(paths['anat_T1_dir']):
        print(f"Error: {paths['anat_T1_dir']} does not exist")
        return False

    ## 1. Skull-strip the T1 image
    print("1. Performing skull-stripping on T1 image...")
    cmd = f"bet2 {paths['anat_T1_dir']} {paths['anat_brain_dir']}"
    os.system(cmd)

    print("2. Performing tissue segmentation...")
    cmd = f"fast -B -I 10 -l 10 {paths['anat_brain_dir']}"
    os.system(cmd)

    ## 2. Create a binary white-matter map
    print("3. Creating binary white matter map...")
    cmd = f"fslmaths {paths['anat_pve2_dir']} -thr 0.5 -bin {paths['anat_wmseg_dir']}"
    os.system(cmd)

    ## 3. Create a binary gray-matter map
    print("4. Creating binary gray matter map...")
    cmd = f"fslmaths {paths['anat_pve1_dir']} -thr 0.4 -bin {paths['anat_gm_mask_dir']}"
    os.system(cmd)

    print("=== Anatomical image processing complete ===\n")
    return True

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='fMRI data preprocessing script',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('subject_id',
                       help='Subject ID (e.g., sub-02)')

    parser.add_argument('--skip-anatomical',
                       action='store_true',
                       help='Skip anatomical image processing')

    parser.add_argument('--sessions',
                       nargs='+',
                       help='Specify sessions to process (e.g., ses-01 ses-02)')

    parser.add_argument('--base-dir',
                       default='/share/home/lq/data/our_nsd/bids',
                       help='Data root directory (default: /share/home/lq/data/our_nsd/bids)')

    return parser.parse_args()

def main():
    """Main function: coordinates the entire processing workflow."""
    # Parse command-line arguments
    args = parse_arguments()

    # Update global base_dir
    global base_dir
    base_dir = args.base_dir

    freesurfer_sub_id = args.subject_id

    print(f"Starting processing for {freesurfer_sub_id}")
    print(f"Data directory: {base_dir}")

    # Check if subject directory exists
    subject_dir = os.path.join(base_dir, freesurfer_sub_id)
    if not os.path.exists(subject_dir):
        print(f"Error: Subject directory does not exist: {subject_dir}")
        sys.exit(1)

    # Determine which sessions to process
    if args.sessions:
        sessions = args.sessions
        print(f"User-specified sessions: {sessions}")
    else:
        sessions = get_available_sessions(freesurfer_sub_id)
        if not sessions:
            print("Error: No session directories found")
            sys.exit(1)

    # 1. Process anatomical images (unless user chooses to skip)
    if not args.skip_anatomical:
        if not process_anatomical(freesurfer_sub_id):
            print("Anatomical image processing failed, exiting program")
            sys.exit(1)
    else:
        print("Skipping anatomical image processing")

    # 2. Loop through and process fieldmap for each session
    successful_sessions = 0
    failed_sessions = []

    for ses_id in sessions:
        success = process_fieldmap(freesurfer_sub_id, ses_id)
        if success:
            successful_sessions += 1
        else:
            failed_sessions.append(ses_id)

    # Output processing results summary
    print("\n" + "="*50)
    print("Processing Results Summary:")
    print(f"Total sessions attempted: {len(sessions)}")
    print(f"Successfully processed: {successful_sessions} sessions")

    if failed_sessions:
        print(f"Failed sessions: {failed_sessions}")

    print("All processing complete!")

if __name__ == "__main__":
    main()