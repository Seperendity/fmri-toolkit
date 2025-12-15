# GLM Beta Processing Scripts

This directory contains three Python scripts for processing fMRI data using GLM (General Linear Model) analysis. These scripts handle the complete pipeline from volume-to-surface conversion through GLM analysis to beta coefficient extraction.

---

## Table of Contents

1. [glmsingle_session.py](#2-glmsingle_sessionpy)
2. [volume_to_fs_ours_funcdata.py](#1-volume_to_fs_ours_funcdatapy)
3. [save_beta_from_fs5.py](#3-save_beta_from_fs5py)

## 1. glmsingle_session.py

### Purpose

Runs GLMsingle analysis on preprocessed fMRI data to estimate beta coefficients for individual stimuli. GLMsingle improves upon standard GLM by incorporating denoising and regularization techniques.

### Functionality

- Loads design matrices and functional data for specified sessions
- Applies GLMdenoise for noise reduction
- Uses fractional ridge regression for improved beta estimates
- Generates diagnostic figures for quality control
- Supports selective session processing
- Automatically skips already-completed sessions
- Comprehensive logging of processing status

### Prerequisites

- **Python packages**: `glmsingle`, `nibabel`, `numpy`, `scipy`, `matplotlib`
- **Input data**:
  - Preprocessed functional data: `f_downsampled.nii.gz`
  - Design matrices: `session_{XX}.npy` files
- **Directory structure**: Data organized by subject/session/run

### Command-Line Arguments

- `--subject`, `-s`: Subject ID (e.g., `sub-01`). Default: `sub-04`
- `--sessions`, `-ses`: Comma-separated session IDs. Default: `01,03`
  - Numeric format: `1,3,5`
  - Full format: `ses-01,ses-03,ses-05`
  - If not specified, processes all available sessions

### Usage Examples

```bash
# Process all sessions for a subject
python glmsingle_session.py --subject sub-01

# Process specific sessions (numeric format)
python glmsingle_session.py --subject sub-01 --sessions 1,3,5

# Process specific sessions (ses- format)
python glmsingle_session.py --subject sub-01 --sessions ses-01,ses-03,ses-05

# Process single session
python glmsingle_session.py --subject sub-01 --sessions 1
```

### Configuration

Fixed parameters (edit script to modify):

- `stimdur = 3`: Stimulus duration in seconds
- `tr_target = 1.5`: Target TR (repetition time) in seconds
- `preproc_dir = '/xx/data/preprocdata/'`: Input data directory
- `designmatrix_dir = '/xx/data/'`: Design matrix directory
- `Results_dir = '/xx/data/GLMsingle_outputs/'`: Output directory

### Output Files

For each session, creates:

- `TYPED_FITHRF_GLMDENOISE_RR.npy` - Final beta estimates with all optimizations
- Diagnostic figures in `./glm_run_figures_fixed/{subject}_{session}/`
- Log file: `glmsingle_{subject}_log.log`

### Special Features

- **Automatic skip**: Detects and skips sessions with existing output files
- **Memory management**: Frees memory after each session to handle large datasets
- **Error handling**: Continues processing remaining sessions if one fails
- **Logging**: Detailed logs for debugging and tracking progress

## 2. volume_to_fs_funcdata.py

### Purpose

Converts volumetric fMRI functional data to FreeSurfer surface space. This script performs volume-to-surface mapping, transforms data to the fsaverage5 template, and applies surface-based smoothing.

### Functionality

- Generates surface masks for left and right hemispheres
- Converts volumetric functional data to native surface space using FreeSurfer's `mri_vol2surf`
- Transforms surface data to fsaverage5 template space using `mri_surf2surf`
- Applies surface-based smoothing (4mm FWHM) using `mris_fwhm`
- Processes all sessions and runs for specified subjects

### Prerequisites

- **FreeSurfer**: Must be installed at `/usr/local/freesurfer`
- **FSL**: Must be installed at `/usr/local/fsl`
- **Python packages**: `nibabel`, `numpy`
- **Input data**: Preprocessed functional data in `f_downsampled.nii.gz` format
- **Registration files**: `register.dof6.dat` and `register.downsampledfunc.dof6.dat`

### Configuration

Edit the script to modify these parameters:

```python
subjects = ['sub-03']  # List of subjects to process
surf_fwhm = 4          # Surface smoothing FWHM in mm
preproc_dir = '/xx/data/preprocdata/'  # Input directory
```

### Usage

```bash
python volume_to_fs_ours_funcdata.py
```

**Note**: This script processes sessions starting from session 20 onwards (hardcoded filter: `if int(ses_id_num) < 20: continue`).

### Output Files

For each run, generates:

- `lh.func.nativesurface.nii.gz` - Left hemisphere native surface data
- `rh.func.nativesurface.nii.gz` - Right hemisphere native surface data
- `lh.func.fs5.nii.gz` - Left hemisphere fsaverage5 template data
- `rh.func.fs5.nii.gz` - Right hemisphere fsaverage5 template data
- `lh.func.fs5_sm.nii.gz` - Left hemisphere smoothed data
- `rh.func.fs5_sm.nii.gz` - Right hemisphere smoothed data

---

## 3. save_beta_from_fs5.py

### Purpose

Processes fMRI beta coefficient data from GLMdenoise analysis. This script loads beta values, applies spatial masking using the NSD general mask, removes specific stimuli at regular intervals, and saves the processed data for further analysis.

### Functionality

1. Loads left and right hemisphere beta coefficient files (`.nii.gz` format)
2. Concatenates hemispheric data into a unified brain representation
3. Applies NSD general cortical mask to focus on visually responsive regions
4. Removes stimuli at specified intervals (default: every 20th stimulus)
5. Saves processed data and metadata for each session

### Prerequisites

- **Python packages**: `nibabel`, `numpy`
- **Input files**:
  - Beta coefficient files: `{data_dir}/sub-{subject}/lh.betas_session{session}.{fs_suffix}_sm.nii.gz`
  - NSD general mask files: `/usr/local/freesurfer/subjects/{space}/label/{lh|rh}.nsdgeneral.mgz`

### Command-Line Arguments

- `--sub`: Subject ID (default: `04`)
- `--ses`: Session IDs, supports ranges and multiple values (default: `['01']`)
  - Single session: `--ses 01`
  - Multiple sessions: `--ses 01 02 03`
  - Range format: `--ses 01-05` (processes sessions 01, 02, 03, 04, 05)
  - Mixed format: `--ses 01 03-05 08`
- `--space`: Surface space, choices: `fsaverage` or `fsaverage5` (default: `fsaverage5`)
- `--data_dir`: Root directory containing beta data (default: `/xx/data/fs`)
- `--output_dir`: Output directory (default: `/xx/data/masked_betas`)
- `--stimuli_per_run`: Number of stimuli per run (default: `105`)
- `--remove_every_n`: Remove every nth stimulus (default: `21`)

### Usage Examples

```bash
# Process single subject, multiple sessions
python save_beta_from_fs5.py --sub 01 --ses 01-05 --space fsaverage5

# Process specific sessions with custom parameters
python save_beta_from_fs5.py --sub 04 --ses 01 03 05 --stimuli_per_run 105

# Process with range and individual sessions
python save_beta_from_fs5.py --sub 02 --ses 01-03 07 09-12 --space fsaverage
```

### Input Files

- **Beta coefficient files**:
  - `{data_dir}/sub-{subject}/lh.betas_session{session}.fs5_sm.nii.gz` (for fsaverage5)
  - `{data_dir}/sub-{subject}/rh.betas_session{session}.fs5_sm.nii.gz` (for fsaverage5)
  - `{data_dir}/sub-{subject}/lh.betas_session{session}.fs_sm.nii.gz` (for fsaverage)
  - `{data_dir}/sub-{subject}/rh.betas_session{session}.fs_sm.nii.gz` (for fsaverage)
- **NSD general mask files**:
  - `/usr/local/freesurfer/subjects/{space}/label/lh.nsdgeneral.mgz`
  - `/usr/local/freesurfer/subjects/{space}/label/rh.nsdgeneral.mgz`

### Output Files

For each session, creates:

- `{output_dir}/sub-{subject}/{space}/session{session}_betas.npy` - Processed beta data (masked vertices × stimuli)
- `{output_dir}/sub-{subject}/{space}/session{session}_metadata.npy` - Session metadata including:
  - Session ID
  - Number of stimuli
  - Start and end stimulus indices
  - Data shape
  - Number of vertices

### Key Features

- **Session range parsing**: Supports flexible session specification with ranges (e.g., `01-05`)
- **Stimulus removal**: Automatically removes every 20th stimulus (button press images) from each run
- **NSD masking**: Applies NSD general mask to focus on visually responsive cortical regions
- **Metadata tracking**: Saves comprehensive metadata for each session for downstream analysis
- **Memory efficient**: Processes sessions individually to manage memory usage

### Loading Saved Data

The script includes a helper function to load processed data:

```python
from save_beta_from_fs5 import load_session_betas

# Load beta data and metadata for a specific session
betas, metadata = load_session_betas(
    save_dir='/cy/data/masked_betas/full',
    sub='01',
    ses='01',
    space='fsaverage5'
)

print(f"Beta shape: {betas.shape}")
print(f"Metadata: {metadata}")
```

---

## Pipeline Overview

The typical processing pipeline using these scripts:

1. **Volume-to-Surface Conversion** (`volume_to_fs_ours_funcdata.py`)

   - Converts volumetric fMRI data to surface space
   - Transforms to fsaverage5 template
   - Applies surface smoothing
2. **GLM Analysis** (`glmsingle_session.py`)

   - Estimates beta coefficients for each stimulus
   - Applies denoising and regularization
   - Generates quality control figures
3. **Beta Extraction and Masking** (`save_beta_from_fs5.py`)

   - Loads GLM beta coefficients
   - Applies NSD cortical mask
   - Removes non-stimulus trials
   - Saves processed data for analysis

---

## Directory Structure

Expected input/output directory structure:

```
/xx/data/
├── preprocdata/
│   └── sub-XX/
│       ├── ses-XX/
│       │   └── run-XX/
│       │       └── f_downsampled.nii.gz
│       ├── ROIs/
│       └── reg_output/
├── design_matrix/
│   └── sub-XX/
│       └── session_XX.npy
├── GLMsingle_outputs/
│   └── sub-XX/
│       └── func1pt8mm/
│           └── ses-XX/
│               └── TYPED_FITHRF_GLMDENOISE_RR.npy
└── masked_betas/
    └── full/
        └── sub-XX/
            └── fsaverage5/
                ├── sessionXX_betas.npy
                └── sessionXX_metadata.npy
```

---

## Notes

- All scripts assume data is organized following BIDS-like conventions
- FreeSurfer and FSL must be properly installed and configured
- Processing can be memory-intensive; monitor system resources
- Check log files for detailed processing information and errors
- Adjust hardcoded paths in scripts to match your local environment
