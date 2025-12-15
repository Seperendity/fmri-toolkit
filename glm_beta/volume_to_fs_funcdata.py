import os
import nibabel as nib
from os.path import join, exists
import numpy as np
import glob

# Add FSL bin path to environment variables for calling FSL commands
os.environ["FSLDIR"] = "/usr/local/fsl"
os.environ["PATH"] += os.pathsep + "/usr/local/fsl/bin"
# Set FreeSurfer installation directory for calling FreeSurfer commands
os.environ["FREESURFER_HOME"] = "/usr/local/freesurfer"
os.environ["PATH"] += os.pathsep + os.path.join(os.environ["FREESURFER_HOME"], "bin")
os.environ["SUBJECTS_DIR"] = "/usr/local/freesurfer/subjects/"

# Data type setting ('nsd' or 'our')
datatype = 'our'

# Set core paths based on data type (only keeping actually used paths)
if datatype == 'nsd':
    preproc_dir = '/opt/data/private/lq/data/NSD/nsddata_preprocdata/'
else:  # datatype == 'our'
    preproc_dir = '/cy/data/preprocdata/'

# Results save directory
results_dir = "/cy/data/betas_test/"

# Experiment parameters (only keeping necessary configuration)
subjects = ['sub-03']  # Subject list
#sessions = ['ses-01']  # Session list
glm_types = ['betas_assumehrf', 'betas_fithrf', 'betas_fithrf_GLMdenoise_RR']
used_glm = glm_types[1]  # Directly specify the GLM type to use (originally glm_idx=2)
surf_fwhm = 4  # Surface smoothing parameter

def get_all_sessions(subject_path):
    session_pattern = join(subject_path, 'ses-*')
    session_dirs = sorted(glob.glob(session_pattern))
    return [os.path.basename(ses_dir) for ses_dir in session_dirs]

def get_all_runs(session_path):
    """Get all run folders under the session directory."""
    run_pattern = join(session_path, 'run-*')
    run_dirs = sorted(glob.glob(run_pattern))
    return [os.path.basename(run_dir) for run_dir in run_dirs]

# Iterate over subjects (currently only 1 subject, keeping loop structure for extensibility)
for sub in subjects:
    # Subject name processing (NSD data needs special conversion, other data keeps original name)
    sub_recon_name = sub.replace("sub-", "sub")

    # Core path definitions (only keeping actually used paths)
    roi_path = join(preproc_dir, sub, 'ROIs')  # ROI-related path
    reg_output_path = join(preproc_dir, sub, 'reg_output')  # Registration output path

    # Registration file paths
    reg_file = join(reg_output_path, 'register.dof6.dat')
    reg_file_downsampled = join(reg_output_path, 'register.downsampledfunc.dof6.dat')

    sessions = get_all_sessions(join(preproc_dir, sub))
    # Iterate over each session
    for session in sessions:
        session_path = join(preproc_dir, sub, session)

        ses_id_num = session.split('-')[1]
        if int(ses_id_num) < 20:
            continue

        # Get all runs in the current session
        runs = get_all_runs(session_path)
        print(f"Processing {sub} {session}, found runs: {runs}")

        # Generate surface masks (only on first run, subsequent runs share the same masks)
        mask_lh = join(roi_path, 'brain_mask.self.lh.nii.gz')
        mask_rh = join(roi_path, 'brain_mask.self.rh.nii.gz')

        if not exists(mask_lh) and runs:  # If mask doesn't exist and there are run folders
            first_run_path = join(session_path, runs[0])
            first_vol = join(first_run_path, 'f_skip_stc_mc_fm_1vol.nii.gz')

            if exists(first_vol):
                print(f"Generating surface masks using {runs[0]}...")

                # Left hemisphere mask generation (calling FreeSurfer command)
                cmd_lh = (f"mri_vol2surf --mov {first_vol} "
                          f"--reg {reg_file} --trgsubject {sub_recon_name} "
                          f"--interp nearest --projfrac 0.5 --hemi lh --o {mask_lh} --noreshape --cortex")
                os.system(cmd_lh)
                # Binarize mask
                os.system(f"mri_binarize --i {mask_lh} --min 0.00001 --o {mask_lh}")

                # Right hemisphere mask generation
                cmd_rh = (f"mri_vol2surf --mov {first_vol} "
                          f"--reg {reg_file} --trgsubject {sub_recon_name} "
                          f"--interp nearest --projfrac 0.5 --hemi rh --o {mask_rh} --noreshape --cortex")
                os.system(cmd_rh)
                os.system(f"mri_binarize --i {mask_rh} --min 0.00001 --o {mask_rh}")

        # Iterate over each run
        for run in runs:
            run_func_path = join(session_path, run)
            print(f"Processing {sub} {session} {run}...")

            # Functional image path and loading
            func_downsample = join(run_func_path, 'f_downsampled.nii.gz')

            # Check if input file exists
            if not exists(func_downsample):
                print(f"Warning: {func_downsample} not found, skipping {run}")
                continue

            func_img = nib.load(func_downsample)
            func_affine = func_img.affine  # Affine matrix of functional image

            # Surface data paths
            surf_lh = join(run_func_path, 'lh.func.nativesurface.nii.gz')
            surf_rh = join(run_func_path, 'rh.func.nativesurface.nii.gz')
            fs_lh = join(run_func_path, 'lh.func.fs5.nii.gz')
            fs_rh = join(run_func_path, 'rh.func.fs5.nii.gz')
            fs_lh_sm = join(run_func_path, 'lh.func.fs5_sm.nii.gz')
            fs_rh_sm = join(run_func_path, 'rh.func.fs5_sm.nii.gz')

            # Volume to Surface conversion (left hemisphere)
            print(f"  Converting volume to surface (left hemisphere)...")
            cmd_vol2surf_lh = (f"mri_vol2surf --mov {func_downsample} "
                              f"--reg {reg_file_downsampled} --trgsubject {sub_recon_name} "
                              f"--interp trilin --projfrac 0.5 --hemi lh --o {surf_lh} --noreshape --cortex")
            os.system(cmd_vol2surf_lh)

            # Volume to Surface conversion (right hemisphere)
            print(f"  Converting volume to surface (right hemisphere)...")
            cmd_vol2surf_rh = (f"mri_vol2surf --mov {func_downsample} "
                              f"--reg {reg_file_downsampled} --trgsubject {sub_recon_name} "
                              f"--interp trilin --projfrac 0.5 --hemi rh --o {surf_rh} --noreshape --cortex")
            os.system(cmd_vol2surf_rh)

            # Transform surface data to fsaverage5 template (left hemisphere)
            print(f"  Converting to fsaverage5 template (left hemisphere)...")
            cmd_surf2fs_lh = (f"mri_surf2surf --sval {surf_lh} --srcsubject {sub_recon_name} "
                             f"--trgsubject fsaverage5 --hemi lh --tval {fs_lh} --cortex")
            os.system(cmd_surf2fs_lh)

            # Transform surface data to fsaverage5 template (right hemisphere)
            print(f"  Converting to fsaverage5 template (right hemisphere)...")
            cmd_surf2fs_rh = (f"mri_surf2surf --sval {surf_rh} --srcsubject {sub_recon_name} "
                             f"--trgsubject fsaverage5 --hemi rh --tval {fs_rh} --cortex")
            os.system(cmd_surf2fs_rh)


            cmd_smooth_lh = (f"mris_fwhm --s fsaverage5 --hemi lh --smooth-only "
                        f"--i {fs_lh} --fwhm 4.0 --o {fs_lh_sm} --cortex")
            os.system(cmd_smooth_lh)

            cmd_smooth_rh = (f"mris_fwhm --s fsaverage5 --hemi rh --smooth-only "
                                    f"--i {fs_rh} --fwhm 4.0 --o {fs_rh_sm} --cortex")
            os.system(cmd_smooth_rh)
            
            print(f"  Completed processing {run}")

print("All runs processed successfully!")