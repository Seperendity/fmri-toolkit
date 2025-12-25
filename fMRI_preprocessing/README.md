# fMRI Preprocessing Pipeline

Parallel batch preprocessing pipeline for functional MRI data with automated slice timing correction, motion correction, fieldmap-based distortion correction, and registration to anatomical space.

**Authors:**

- wrx: Core module implementation, logic design
- cy: Algorithm acceleration, resource optimization
- lq: Code style compliance, security validation

**Date:** 2025-07-27

---

## Table of Contents

- [Overview](#overview)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Processing Steps](#processing-steps)
- [Parallel Processing](#parallel-processing)
- [Usage](#usage)
- [Command-Line Arguments](#command-line-arguments)
- [Configuration](#configuration)
- [Resume Functionality](#resume-functionality)
- [Requirements](#requirements)
- [Anatomical and Fieldmap Preprocessing](#anatomical-and-fieldmap-preprocessing)

---

## Overview

The `fmri_preprocessing_parallel.py` script performs automated fMRI preprocessing in parallel batches with the following capabilities:

- **Slice timing correction** using custom timing files
- **Motion correction** to mean volume using MCFLIRT
- **Fieldmap-based distortion correction** using FSL's epi_reg
- **Boundary-based registration (BBR)** to FreeSurfer anatomical surfaces
- **Parallel processing** across multiple subjects and runs
- **Resume functionality** to skip already-processed sessions
- **Comprehensive logging** with detailed execution tracking

---

## Input Files

### Anatomical Data (per subject, ses-01 only)

```
{subject}_ses-01_T1w.nii.gz                    # Raw T1-weighted anatomical image
T1_brain.nii.gz                                # Skull-stripped T1 brain (from anatomical preprocessing)
T1_wmseg.nii.gz                                # White matter segmentation mask
```

### Fieldmap Data (per session)

```
{subject}_{session}_magnitude1.nii.gz          # Fieldmap magnitude image  
fmap_mag_brain.nii.gz                          # Skull-stripped fieldmap magnitude
fmap_rads.nii.gz                               # Fieldmap in radians for distortion correction
```

### Functional Data (per run)

```
{subject}_{session}_run-{XX}_bold.nii.gz       # Raw functional EPI timeseries
```

### Timing Files

```
tcustom_epi.txt                                # Custom slice timing file
```

---

## Output Files

### Motion Correction Output (per run)

```
f_skip_stc.nii.gz                              # Slice-timing corrected data
f_skip_stc_mc.nii.gz                           # Motion-corrected data
f_skip_stc_mc.nii.gz.par                       # Motion parameters file
f_skip_stc_mc.nii.gz_mean_reg.nii.gz          # Mean volume for registration
```

### Distortion Correction Output (per run)

```
f_skip_stc_mc_fm_warp.nii.gz                   # Warp field for distortion correction
f_skip_stc_mc_fm.nii.gz                        # Final preprocessed data (registered to T1 space)
f_downsampled.nii.gz                           # Downsampled version at original resolution
```

### Registration Output (shared per subject)

```
register.dof6.dat                              # BBR registration matrix (functional to T1)
init.register.dof6.dat                         # Initial registration matrix
register.downsampledfunc.dof6.dat              # Registration for downsampled data
T1_brain.registered.nii.gz                     # T1 brain registered to functional space
T1_global_mask.nii.gz                          # Global brain mask in T1 space
global_mask_func.nii.gz                        # Global brain mask in functional space
```

### Process Control

```
.preproc_success                               # Success flag indicating session completion
```

---

## Processing Steps

1. **Slice timing correction** using custom timing file
2. **Motion correction** to mean volume using MCFLIRT
3. **Fieldmap-based distortion correction** using epi_reg
4. **Application of distortion correction** to full timeseries
5. **Boundary-based registration (BBR)** to FreeSurfer anatomical surface
6. **Generation of brain masks** and registration matrices
7. **Downsampling** to original voxel resolution

---

## Parallel Processing

The pipeline supports efficient parallel processing with the following features:

- **Multi-worker parallel execution** across subjects
- **Concurrent run processing** within each worker (configurable via `MAX_CONCURRENT_RUNS`)
- **Resume functionality** using success flags (`.preproc_success`) to skip completed sessions
- **Batched processing** to manage memory and computational resources
- **Automatic workload distribution** across workers

### Parallel Processing Architecture

- Subjects are automatically divided among workers
- Each worker processes multiple runs concurrently (default: 10 concurrent runs)
- Sessions are processed in batches to optimize resource usage
- Success flags prevent reprocessing of completed sessions

---

## Usage

### Single Worker Processing

Process all subjects with a single worker:

```bash
python fmri_preprocessing_parallel.py
```

### Multi-Worker Parallel Processing

Run multiple workers in parallel to process different subjects simultaneously:

**Worker 0 (Terminal 1):**

```bash
python fmri_preprocessing_parallel.py --total-workers 2 --worker-id 0
```

**Worker 1 (Terminal 2):**

```bash
python fmri_preprocessing_parallel.py --total-workers 2 --worker-id 1
```

### Example with 4 Workers

```bash
# Terminal 1
python fmri_preprocessing_parallel.py --total-workers 4 --worker-id 0

# Terminal 2
python fmri_preprocessing_parallel.py --total-workers 4 --worker-id 1

# Terminal 3
python fmri_preprocessing_parallel.py --total-workers 4 --worker-id 2

# Terminal 4
python fmri_preprocessing_parallel.py --total-workers 4 --worker-id 3
```

---

## Command-Line Arguments

- `--total-workers` - Total number of parallel worker scripts being run (default: 2)
- `--worker-id` - Zero-indexed ID of this specific worker instance (default: 0)

**Note:** `worker-id` must be less than `total-workers`

---

## Configuration

### Script Configuration Variables

Edit the following variables in the script to match your dataset:

```python
# Data type selection
datatype = 'our'  # 'nsd' or 'our'

# Directory paths
base_dir = '/path/to/bids/data'
preproc_dir = '/path/to/preprocessed/output'

# Acquisition parameters
TR_s = 1.5                          # Repetition time in seconds
Effective_echo_spacing = 0.000275003  # Echo spacing for distortion correction
origvoxelSize = [2.0, 2.0, 2.0]     # Original voxel size in mm

# Parallel processing
MAX_CONCURRENT_RUNS = 10            # Maximum concurrent runs per worker
```

### Environment Variables

The script automatically sets up FSL and FreeSurfer environments:

```python
os.environ["FSLDIR"] = "/path/to/fsl"
os.environ["FREESURFER_HOME"] = "/path/to/freesurfer"
os.environ["SUBJECTS_DIR"] = "/path/to/freesurfer/subjects"
```

---

## Resume Functionality

The pipeline implements automatic resume functionality:

- **Success flags** (`.preproc_success`) are created after successful session completion
- **Automatic skip** of already-processed sessions on subsequent runs
- **Partial completion handling** - sessions with incomplete runs will be reprocessed
- **Per-session tracking** - each session is tracked independently

This allows you to:

- Resume processing after interruptions
- Add new sessions without reprocessing existing data
- Safely re-run the script without duplicating work

---

## Requirements

### Software Dependencies

- **FSL** (FMRIB Software Library, v6.0.5) - for motion correction and distortion correction
- **FreeSurfer** (v6.0.0)- for anatomical surface reconstruction and registration

### FSL Tools Used

- `slicetimer` - Slice timing correction
- `mcflirt` - Motion correction
- `epi_reg` - EPI registration with fieldmap correction
- `applywarp` - Apply warp fields
- `mri_convert` - Format conversion and resampling
- `bbregister` - Boundary-based registration
- `mri_vol2vol` - Volume-to-volume registration
- `mri_binarize` - Binary mask creation

---

## Logging

The script generates detailed log files in the `log/` directory:

```
log/fmri_preprocessing_parallel_batch_run_worker-{worker_id}-of-{total_workers}_{timestamp}.log
```

Log files include:

- Work queue for each worker
- Batch processing progress
- Command execution details (stdout/stderr)
- Error messages and warnings
- Session completion status
- Total execution time

---

## Important Notes

### First Run Special Processing

The first run of the first session (`ses-01`, `run-01`) generates shared files:

- Registration matrices
- Brain masks
- Anatomical-to-functional transformations

These files are reused for all subsequent runs to ensure consistency.

### Temporary File Cleanup

Intermediate files (`f_skip_stc.nii.gz`, `f_skip_stc_mc.nii.gz`) are automatically removed after processing to save disk space.

---

## Anatomical and Fieldmap Preprocessing

The `fmri_anat_preprocessing.py` script performs automated preprocessing including anatomical skull-stripping, tissue segmentation, and fieldmap correction for multiple sessions of a specified subject.

**Author:** wrx & cy
**Date:** 2025-07-26

### Input Files

#### Anatomical (per subject, ses-01 only)

```
{subject}_ses-01_T1w.nii.gz                    # Raw T1-weighted anatomical image
```

#### Fieldmap (per session)

```
{subject}_{session}_magnitude1.nii.gz          # Fieldmap magnitude image
{subject}_{session}_phasediff.nii.gz           # Fieldmap phase difference image
```

### Output Files

#### Anatomical Processing Output

```
T1_brain.nii.gz                                # Skull-stripped T1 brain
T1_brain_pve_0.nii.gz                          # CSF probability map
T1_brain_pve_1.nii.gz                          # Gray matter probability map
T1_brain_pve_2.nii.gz                          # White matter probability map
T1_wmseg.nii.gz                                # Binary white matter mask (threshold: 0.5)
T1_gm_mask.nii.gz                              # Binary gray matter mask (threshold: 0.4)
```

#### Fieldmap Processing Output (per session)

```
fmap_mag_brain1.nii.gz                         # Initial skull-stripped magnitude image
fmap_mag_brain.nii.gz                          # Final processed magnitude brain (eroded)
fmap_rads.nii.gz                               # Fieldmap in radians for distortion correction
```

#### Quality Control Images (per session)

```
field_brain_x_slices.png                       # Magnitude brain sagittal view QC
field_brain_y_slices.png                       # Magnitude brain coronal view QC
field_brain_z_slices.png                       # Magnitude brain axial view QC
field_fmap_x_slices.png                        # Fieldmap sagittal view QC
field_fmap_y_slices.png                        # Fieldmap coronal view QC
field_fmap_z_slices.png                        # Fieldmap axial view QC
```

### Anatomical Preprocessing Usage

```bash
# Basic usage - process all sessions for sub-02
python fmri_anat_preprocessing.py sub-02

# Skip anatomical image processing
python fmri_anat_preprocessing.py sub-02 --skip-anatomical

# Process only specific sessions
python fmri_anat_preprocessing.py sub-02 --sessions ses-01 ses-02 ses-05

# Use custom data directory
python fmri_anat_preprocessing.py sub-02 --base-dir /path/to/your/data
```

---

## License

Please refer to the main project license.
