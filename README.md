# fMRI Toolkit

[![Python 3.9+](https://img.shields.io/badge/python-3.9.21-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![FreeSurfer](https://img.shields.io/badge/FreeSurfer-6.0.0-green.svg)](https://surfer.nmr.mgh.harvard.edu/)
[![FSL](https://img.shields.io/badge/FSL-6.0.5-orange.svg)](https://fsl.fmrib.ox.ac.uk/fsl/)

A comprehensive toolkit for processing functional MRI (fMRI) data, from DICOM conversion through preprocessing to GLM-based beta coefficient extraction. This toolkit provides an end-to-end pipeline for neuroimaging research with support for parallel processing and surface-based analysis.

---

## Table of Contents

- [fMRI Toolkit](#fmri-toolkit)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [Features](#features)
  - [Directory Structure](#directory-structure)
  - [Installation \& Prerequisites](#installation--prerequisites)
    - [Software Requirements](#software-requirements)
    - [Python Dependencies](#python-dependencies)
    - [Environment Setup](#environment-setup)
  - [Pipeline Workflow](#pipeline-workflow)
  - [Quick Start](#quick-start)
    - [Step 1: Convert DICOM to BIDS](#step-1-convert-dicom-to-bids)
    - [Step 2: Preprocess fMRI Data](#step-2-preprocess-fmri-data)
    - [Step 3: Run GLM Analysis](#step-3-run-glm-analysis)
  - [Module Documentation](#module-documentation)
    - [1. BIDS Conversion](#1-bids-conversion)
    - [2. fMRI Preprocessing](#2-fmri-preprocessing)
    - [3. GLM Beta Analysis](#3-glm-beta-analysis)
  - [Data Organization](#data-organization)
    - [Expected Input Structure (BIDS)](#expected-input-structure-bids)
    - [Output Structure (Preprocessed)](#output-structure-preprocessed)
    - [Output Structure (GLM Results)](#output-structure-glm-results)
  - [Contributing](#contributing)
    - [Code Style](#code-style)
  - [License](#license)
  - [Acknowledgments](#acknowledgments)

---

## Overview

The **fMRI Toolkit** is designed for researchers working with functional neuroimaging data. It provides a modular, reproducible pipeline that handles:

1. **DICOM to BIDS Conversion** — Automated batch conversion of raw DICOM files to the standardized Brain Imaging Data Structure (BIDS) format
2. **fMRI Preprocessing** — Parallel preprocessing including slice timing correction, motion correction, distortion correction, and registration
3. **GLM Analysis** — General Linear Model analysis using GLMsingle for single-trial beta estimation in volume space, followed by surface conversion and visual cortex extraction

The toolkit is optimized for large-scale studies with multiple subjects and sessions, featuring parallel processing capabilities and resume functionality.

---

## Features

- 🔄 **Automated BIDS Conversion** — Batch processing of DICOM files with validation
- ⚡ **Parallel Processing** — Multi-worker architecture for efficient preprocessing
- 🧠 **Surface-Based Analysis** — FreeSurfer integration for cortical surface mapping
- 📊 **GLMsingle Integration** — State-of-the-art single-trial beta estimation
- 🔧 **Distortion Correction** — Fieldmap-based EPI distortion correction
- 📁 **Resume Functionality** — Automatic detection and skipping of completed sessions
- 📝 **Comprehensive Logging** — Detailed logs for debugging and quality control
- ✅ **BIDS Validation** — Built-in tools to verify conversion accuracy

---

## Directory Structure

```
fmri-toolkit/
├── README.md                          # This file
├── bids_conversion/                   # DICOM to BIDS conversion tools
│   ├── README.md                      # Detailed module documentation
│   ├── dcm2bids_batch.py              # Batch DICOM conversion script
│   ├── check_bids.py                  # BIDS validation tool
│   ├── config.json                    # dcm2bids configuration
│   └── run_scripts/                   # Helper scripts
├── fMRI_preprocessing/                # Preprocessing pipeline
│   ├── README.md                      # Detailed module documentation
│   └── fmri_preprocessing_parallel.py # Parallel preprocessing script
└── glm_beta/                          # GLM analysis and beta extraction
    ├── README.md                      # Detailed module documentation
    ├── glmsingle_session.py           # GLMsingle analysis script
    ├── volume_to_fs_funcdata.py       # Volume-to-surface conversion
    └── save_beta_from_fs5.py          # Beta coefficient extraction
```

---

## Installation & Prerequisites

### Software Requirements

| Software             | Version       | Purpose                                     |
| -------------------- | ------------- | ------------------------------------------- |
| **Python**     | 3.9.21        | Core scripting language                     |
| **FreeSurfer** | 6.0.0         | Surface reconstruction, registration        |
| **FSL**        | 6.0.5         | Motion/distortion correction, preprocessing |
| **dcm2bids**   | 3.2.0         | DICOM to BIDS conversion                    |
| **dcm2niix**   | v1.0.20250505 | NIfTI conversion (dcm2bids dependency)      |

### Python Dependencies

```bash
# Core dependencies
pip install numpy scipy nibabel matplotlib

# For GLM analysis
pip install git+https://github.com/cvnlab/GLMsingle.git

# For BIDS validation
pip install pydicom
```

### Environment Setup

Ensure FreeSurfer and FSL are properly configured:

```bash
# FreeSurfer setup
export FREESURFER_HOME=/usr/local/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=/path/to/your/subjects

# FSL setup
export FSLDIR=/usr/local/fsl
source $FSLDIR/etc/fslconf/fsl.sh
```

---

## Pipeline Workflow

The toolkit implements a three-stage pipeline for complete fMRI data processing:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           fMRI PROCESSING PIPELINE                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────────────┐      ┌──────────────────┐     ┌──────────────────┐    │
│  │  STAGE 1         │      │  STAGE 2         │     │  STAGE 3         │    │
│  │  BIDS Conversion │──▶  │ Preprocessing    |───▶ │  GLM Analysis    │    │  
│  └──────────────────┘      └──────────────────┘     └──────────────────┘    │
│         │                        │                       │                  │
│         ▼                        ▼                       ▼                  │
│  • DICOM → NIfTI          • Slice Timing          • GLMsingle (volume)      │
│  • BIDS organization      • Motion Correction     • Beta → Surface          │
│  • Metadata extraction    • Distortion Correction • Visual cortex extract   │
│  • Validation             • Registration (BBR)    • fsaverage5 mapping      │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Quick Start

### Step 1: Convert DICOM to BIDS

```bash
cd bids_conversion

# Convert sessions 1-10 for subject 01
python dcm2bids_batch.py -s 01 --start-session 1 --end-session 10 \
    -b /path/to/dicom/data \
    -c config.json \
    -o /path/to/bids/output

# Validate the conversion
python check_bids.py
```

### Step 2: Preprocess fMRI Data

```bash
cd fMRI_preprocessing

# Single worker processing
python fmri_preprocessing_parallel.py

# Multi-worker parallel processing (run in separate terminals)
python fmri_preprocessing_parallel.py --total-workers 4 --worker-id 0
python fmri_preprocessing_parallel.py --total-workers 4 --worker-id 1
python fmri_preprocessing_parallel.py --total-workers 4 --worker-id 2
python fmri_preprocessing_parallel.py --total-workers 4 --worker-id 3
```

### Step 3: Run GLM Analysis

```bash
cd glm_beta

# Run GLMsingle analysis on volume space data
python glmsingle_session.py --subject sub-01 --sessions 1,2,3,4,5

# Convert beta coefficients from volume to surface space
python volume_to_fs_funcdata.py

# Extract visual cortex vertices and save final betas
python save_beta_from_fs5.py --sub 01 --ses 01-05 --space fsaverage5
```

---

## Module Documentation

### 1. BIDS Conversion

**Location:** `bids_conversion/`

Automated batch conversion of DICOM files to BIDS format with validation.

| Script                | Description                                                    |
| --------------------- | -------------------------------------------------------------- |
| `dcm2bids_batch.py` | Batch DICOM to BIDS conversion with auto session detection     |
| `check_bids.py`     | Validates conversion accuracy by comparing DICOM and BIDS data |
| `config.json`       | Configuration file defining conversion rules                   |

**Key Features:**

- Automatic session directory detection
- Comprehensive logging with timestamps
- Support for custom session naming patterns
- Validation of DICOM file presence before conversion

📖 **Detailed documentation:** [bids_conversion/README.md](bids_conversion/README.md)

---

### 2. fMRI Preprocessing

**Location:** `fMRI_preprocessing/`

Parallel batch preprocessing pipeline with full distortion correction and registration.

| Script                             | Description                                         |
| ---------------------------------- | --------------------------------------------------- |
| `fmri_preprocessing_parallel.py` | Main preprocessing pipeline with parallel execution |

**Processing Steps:**

1. Slice timing correction
2. Motion correction (MCFLIRT)
3. Fieldmap-based distortion correction
4. Boundary-based registration (BBR) to FreeSurfer surfaces
5. Brain mask generation
6. Downsampling to original resolution

**Key Features:**

- Multi-worker parallel execution
- Resume functionality (skips completed sessions)
- Concurrent run processing within each worker
- Automatic cleanup of intermediate files

📖 **Detailed documentation:** [fMRI_preprocessing/README.md](fMRI_preprocessing/README.md)

---

### 3. GLM Beta Analysis

**Location:** `glm_beta/`

GLM analysis pipeline for single-trial beta coefficient estimation in volume space, followed by surface projection and visual cortex extraction.

| Script                       | Description                                                     |
| ---------------------------- | --------------------------------------------------------------- |
| `glmsingle_session.py`     | GLMsingle analysis for beta estimation (volume space)           |
| `volume_to_fs_funcdata.py` | Convert volume betas to FreeSurfer surface space with smoothing |
| `save_beta_from_fs5.py`    | Extract visual cortex vertices and save final betas             |

**Processing Steps:**

1. GLMsingle analysis on volume space data with denoising
2. Fractional ridge regression for beta estimation
3. Volume-to-surface conversion of beta coefficients
4. Transformation to fsaverage5 template with smoothing (4mm FWHM)
5. Visual cortex extraction using NSD cortical masks

**Key Features:**

- GLMdenoise for noise reduction in volume space
- Fractional ridge regression for robust beta estimation
- Post-hoc surface projection preserving beta accuracy
- Visual cortex vertex extraction for downstream analysis
- Flexible session specification with ranges

📖 **Detailed documentation:** [glm_beta/README.md](glm_beta/README.md)

---

## Data Organization

### Expected Input Structure (BIDS)

```
bids_dataset/
├── sub-01/
│   ├── ses-01/
│   │   ├── anat/
│   │   │   └── sub-01_ses-01_T1w.nii.gz
│   │   ├── func/
│   │   │   ├── sub-01_ses-01_task-visual_run-01_bold.nii.gz
│   │   │   ├── sub-01_ses-01_task-visual_run-02_bold.nii.gz
│   │   │   └── ...
│   │   └── fmap/
│   │       ├── sub-01_ses-01_magnitude1.nii.gz
│   │       └── sub-01_ses-01_phasediff.nii.gz
│   ├── ses-02/
│   │   └── ...
│   └── ...
└── sub-02/
    └── ...
```

### Output Structure (Preprocessed)

```
preprocessed_data/
├── sub-01/
│   ├── ses-01/
│   │   └── run-01/
│   │       ├── f_downsampled.nii.gz          # Preprocessed functional
│   │       └── f_skip_stc_mc_fm.nii.gz       # Full resolution output
│   ├── ROIs/
│   │   └── brain_mask.*.nii.gz               # Brain masks
│   └── reg_output/
│       └── register.dof6.dat                 # Registration matrices
└── ...
```

### Output Structure (GLM Results)

```
glm_outputs/
├── GLMsingle_outputs/
│   └── sub-01/
│       └── func1pt8mm/
│           └── ses-01/
│               └── TYPED_FITHRF_GLMDENOISE_RR.npy
└── masked_betas/
    └── full/
        └── sub-01/
            └── fsaverage5/
                ├── session01_betas.npy
                └── session01_metadata.npy
```

---

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Code Style

- Follow PEP 8 for Python code
- Include docstrings for all functions
- Add comments for complex logic
- Update documentation for new features

---

## License

This project is licensed under the MIT License - see below for details.

```
MIT License

Copyright (c) 2025 fMRI Toolkit Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## Acknowledgments

This toolkit builds upon several excellent open-source projects:

- [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) — Cortical surface reconstruction
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/) — FMRI preprocessing tools
- [dcm2bids](https://unfmontreal.github.io/Dcm2Bids/) — DICOM to BIDS conversion
- [GLMsingle](https://github.com/cvnlab/GLMsingle) — Single-trial beta estimation
- [NiBabel](https://nipy.org/nibabel/) — Neuroimaging file I/O

---

<p align="center">
  Made with ❤️ for the neuroimaging community
</p>
