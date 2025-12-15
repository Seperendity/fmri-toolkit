import argparse
from os.path import join, exists
import nibabel as nib
import numpy as np
import os

def parse_session_ranges(session_args):
    """
    Parse session arguments with support for range input.

    Args:
        session_args: List of session arguments from command line

    Returns:
        sessions: Parsed list of sessions (string format, preserving zero-padding)

    Example:
        ['01', '02', '05-08'] -> ['01', '02', '05', '06', '07', '08']
    """
    sessions = []

    for arg in session_args:
        if '-' in arg and len(arg.split('-')) == 2:
            # Handle range input, e.g., '01-05'
            start_str, end_str = arg.split('-')
            try:
                start = int(start_str)
                end = int(end_str)

                # Ensure the range is valid
                if start > end:
                    raise ValueError(f"Invalid range: {arg}. Start ({start}) must be <= end ({end})")

                # Generate all sessions in the range, preserving zero-padding format
                width = len(start_str)  # Preserve original digit width
                for i in range(start, end + 1):
                    sessions.append(str(i).zfill(width))

            except ValueError as e:
                if "invalid literal" in str(e):
                    raise ValueError(f"Invalid session range format: {arg}. Use format like '01-05'")
                else:
                    raise e
        else:
            # Handle single session
            sessions.append(arg)

    # Remove duplicates and sort by numeric value
    unique_sessions = list(set(sessions))
    unique_sessions.sort(key=int)

    return unique_sessions

def create_removal_mask(num_runs, stimuli_per_run=105, remove_every_n=21):
    """
    Create a mask to remove button-press images after every 20 stimuli.

    Args:
        num_runs: Number of runs
        stimuli_per_run: Number of stimuli per run, default 105
        remove_every_n: Remove one stimulus every n stimuli, default 21 (i.e., remove button-press image after every 20th)

    Returns:
        mask: Boolean array, True means keep, False means remove
    """
    total_stimuli = num_runs * stimuli_per_run
    mask = np.ones(total_stimuli, dtype=bool)

    # Calculate indices to be removed
    indices_to_remove = []
    for run_idx in range(num_runs):
        run_start = run_idx * stimuli_per_run
        # Remove button-press images after every 20th stimulus in each run (indices: 20, 41, 62, 83, 104)
        run_remove_indices = np.arange(run_start + 20, run_start + stimuli_per_run, remove_every_n)
        indices_to_remove.extend(run_remove_indices)

    indices_to_remove = np.array(indices_to_remove)
    print(f"Indices to be removed: {indices_to_remove}")
    print(f"Removing {len(indices_to_remove) // num_runs} stimuli per run")

    mask[indices_to_remove] = False
    return mask


def load_and_process_session_data(data_dir, sub, ses, space, stimuli_per_run=105, remove_every_n=21):
    """
    Load and process data for a single session.

    Args:
        data_dir: Root data directory
        sub: Subject ID
        ses: Session ID
        stimuli_per_run: Number of stimuli per run
        remove_every_n: Removal interval

    Returns:
        processed_data: Processed data (vertices, stimuli)
        num_runs: Number of runs in this session
    """
    # Define surface space suffix

    if space == 'fsaverage5':
        fs_suffix = 'fs5'
    else:
        fs_suffix = 'fs'

    # Build file paths
    lh_file = join(data_dir, f'sub-{sub}',
                   f'lh.betas_session{ses}.{fs_suffix}_sm.nii.gz')
    rh_file = join(data_dir, f'sub-{sub}',
                   f'rh.betas_session{ses}.{fs_suffix}_sm.nii.gz')


    print(f"Loading data for sub-{sub} session{ses}...")
    print(f"  -> Found data in: {os.path.dirname(lh_file)}")

    # Load data
    lh = nib.load(lh_file).get_fdata().squeeze()
    rh = nib.load(rh_file).get_fdata().squeeze()

    # Concatenate left and right hemisphere data
    data = np.concatenate((lh, rh), axis=0)

    print(f"Original data shape: {data.shape}")

    # Calculate number of runs
    num_cols = data.shape[1]
    num_runs = num_cols // stimuli_per_run

    if num_cols != num_runs * stimuli_per_run:
        print(f"Warning: Total stimuli ({num_cols}) not divisible by stimuli_per_run ({stimuli_per_run})")
        print(f"Assuming {num_runs} complete runs")

    print(f"Detected {num_runs} runs with {stimuli_per_run} stimuli each")

    # Create removal mask
    removal_mask = create_removal_mask(num_runs, stimuli_per_run, remove_every_n)

    # Apply mask to remove button-press images
    processed_data = data[:, removal_mask]

    print(f"After removal, data shape: {processed_data.shape}")
    print(f"Removed {num_cols - processed_data.shape[1]} stimuli")

    return processed_data, num_runs


def load_nsd_mask(space):
    """
    Load NSD general cortical mask.

    Returns:
        loc_verts_index_lh: Left hemisphere mask indices
        loc_verts_index_rh: Right hemisphere mask indices
        loc_verts_index_combined: Combined mask indices
        n_verts_lh: Number of left hemisphere vertices
        n_verts_rh: Number of right hemisphere vertices
    """
    NSD_dir = join(f'/usr/local/freesurfer/subjects/{space}/label/')
    lh_mask_file = join(NSD_dir, "lh.nsdgeneral.mgz")
    rh_mask_file = join(NSD_dir, "rh.nsdgeneral.mgz")

    if not exists(lh_mask_file) or not exists(rh_mask_file):
        raise FileNotFoundError(f"NSD mask files not found in {NSD_dir}")

    lh_mask_data = nib.load(lh_mask_file).get_fdata().squeeze()
    rh_mask_data = nib.load(rh_mask_file).get_fdata().squeeze()

    n_verts_lh = len(lh_mask_data)
    n_verts_rh = len(rh_mask_data)

    # Get vertex indices within the nsdgeneral mask
    loc_verts_index_lh = np.where(lh_mask_data != 0)[0]
    loc_verts_index_rh = np.where(rh_mask_data != 0)[0]

    # Combine left and right hemisphere mask vertex indices
    loc_verts_index_combined = np.concatenate([
        loc_verts_index_lh,  # Left hemisphere indices
        loc_verts_index_rh + n_verts_lh  # Right hemisphere indices (with offset)
    ])

    print(f"NSD mask info({space} space):")
    print(f"  Left hemisphere vertices: {n_verts_lh}, masked: {len(loc_verts_index_lh)}")
    print(f"  Right hemisphere vertices: {n_verts_rh}, masked: {len(loc_verts_index_rh)}")
    print(f"  Total masked vertices: {len(loc_verts_index_combined)}")

    return (loc_verts_index_lh, loc_verts_index_rh, loc_verts_index_combined,
            n_verts_lh, n_verts_rh)


def save_masked_betas(data, loc_verts_index_combined, save_dir, sub, ses, space, start_stim_idx=0):
    """
    Save masked beta values (entire session saved as one file).

    Args:
        data: Original data (vertices, stimuli)
        loc_verts_index_combined: Mask indices
        save_dir: Save directory
        sub: Subject ID
        ses: Session ID
        start_stim_idx: Starting stimulus index (for continuous numbering across sessions)

    Returns:
        next_stim_idx: Next stimulus index
    """
    # Apply nsdgeneral mask
    masked_data = data[loc_verts_index_combined, :]
    print(f"Masked data shape: {masked_data.shape}")

    # Create subject save directory
    space_save_dir = join(save_dir, f'sub-{sub}')
    subject_save_dir = join(space_save_dir, space)
    os.makedirs(subject_save_dir, exist_ok=True)

    # Save entire session's beta values to one file
    num_stimuli = masked_data.shape[1]

    filename = f'session{ses}_betas.npy'
    filepath = join(subject_save_dir, filename)
    np.save(filepath, masked_data)

    # Save metadata information to a separate file
    metadata = {
        'session': ses,
        'num_stimuli': num_stimuli,
        'start_stim_idx': start_stim_idx,
        'end_stim_idx': start_stim_idx + num_stimuli - 1,
        'data_shape': masked_data.shape,
        'num_vertices': len(loc_verts_index_combined)
    }

    metadata_filename = f'session{ses}_metadata.npy'
    metadata_filepath = join(subject_save_dir, metadata_filename)
    np.save(metadata_filepath, metadata)

    print(f"Session{ses}: Saved {num_stimuli} stimuli to '{filepath}'")
    print(f"Session{ses}: Saved metadata to '{metadata_filepath}'")

    return start_stim_idx + num_stimuli

# Function to load saved data
def load_session_betas(save_dir, sub, ses, space):
    """
    Load beta data for a specified session.

    Args:
        save_dir: Save directory
        sub: Subject ID
        ses: Session ID

    Returns:
        betas: Beta data (vertices, stimuli)
        metadata: Metadata dictionary
    """
    space_save_dir = join(save_dir, space)
    subject_save_dir = join(space_save_dir, f'sub-{sub}')

    # Load beta data
    betas_file = join(subject_save_dir, f'session{ses}_betas.npy')
    betas = np.load(betas_file)

    # Load metadata
    metadata_file = join(subject_save_dir, f'session{ses}_metadata.npy')
    metadata = np.load(metadata_file, allow_pickle=True).item()

    return betas, metadata

def main():
    parser = argparse.ArgumentParser(description='Save masked beta values for multiple subjects and sessions')
    parser.add_argument('--sub', type=str, default='04', 
                        help='Subject ID (default: 01)')
    parser.add_argument('--ses', type=str, nargs='+', default=['01'], 
                        help='Session IDs (default: 01). Can specify multiple sessions like --ses 01 02 03')
    parser.add_argument('--space', type=str, default='fsaverage5', choices=['fsaverage', 'fsaverage5'],
                        help='Surface space (default: fsaverage)')    
    parser.add_argument('--data_dir', type=str, default='/cy/data/fs', 
                        help='Root directory containing beta data (default: /cy/data/betas)')
    parser.add_argument('--output_dir', type=str, default='/cy/data/masked_betas', 
                        help='Output directory (default: masked_betas)')
    parser.add_argument('--stimuli_per_run', type=int, default=105,
                        help='Number of stimuli per run (default: 105)')
    parser.add_argument('--remove_every_n', type=int, default=21,
                        help='Remove every nth stimulus (default: 21, i.e., remove 21st, 42nd, etc.)')
    
    args = parser.parse_args()

    # Parse session arguments
    try:
        sessions = parse_session_ranges(args.ses)
    except ValueError as e:
        print(f"Error parsing session arguments: {e}")
        return

    data_type = 'full'
    args.data_dir = join(args.data_dir, data_type)
    args.output_dir = join(args.output_dir, data_type)

    print("=== Multi-Subject Multi-Session Beta Processing ===")
    print(f"Subject: {args.sub}")
    print(f"Sessions to process: {sessions} (total: {len(sessions)} sessions)")
    print(f"Data directory: {args.data_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Stimuli per run: {args.stimuli_per_run}")
    print(f"Remove every {args.remove_every_n} stimuli")

    # Load NSD mask (shared by all subjects and sessions)
    try:
        (loc_verts_index_lh, loc_verts_index_rh, loc_verts_index_combined,
         n_verts_lh, n_verts_rh) = load_nsd_mask(args.space)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # Create main output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Process all sessions for the specified subject
    total_stimuli_processed = 0
    session_info = {}
    processed_sessions = []
    skipped_sessions = []

    for ses in sessions:
        print(f"\n--- Processing Subject {args.sub}, Session {ses} ---")

        # Load and process session data
        processed_data, num_runs = load_and_process_session_data(
            args.data_dir, args.sub, ses, args.space, args.stimuli_per_run, args.remove_every_n
        )

        if processed_data is None:
            print(f"Skipping session {ses} due to missing files")
            skipped_sessions.append(ses)
            continue

        # Save masked beta values
        next_stim_idx = save_masked_betas(
            processed_data, loc_verts_index_combined, args.output_dir,
            args.sub, ses, args.space, total_stimuli_processed
        )

        # Record session information
        num_stimuli_this_session = next_stim_idx - total_stimuli_processed
        session_info[f'session{ses}'] = {
            'num_runs': num_runs,
            'num_stimuli': num_stimuli_this_session,
            'start_stim_idx': total_stimuli_processed,
            'end_stim_idx': next_stim_idx - 1
        }
        processed_sessions.append(ses)
        total_stimuli_processed = next_stim_idx
    
    print(f"\n=== Processing Complete ===")
    print(f"Successfully processed sessions: {processed_sessions} ({len(processed_sessions)} sessions)")
    if skipped_sessions:
        print(f"Skipped sessions (missing files): {skipped_sessions} ({len(skipped_sessions)} sessions)")
    print(f"Total stimuli processed: {total_stimuli_processed}")
    
    if session_info:
        print(f"\nSession Summary:")
        for ses_name, info in session_info.items():
            print(f"  {ses_name}: {info['num_runs']} runs, {info['num_stimuli']} stimuli "
                  f"(indices {info['start_stim_idx']}-{info['end_stim_idx']})")
        
        print(f"\nData saved to: {args.output_dir}/sub-{args.sub}/")
        print("Each session file contains {} voxel values (nsdgeneral masked vertices)".format(
            len(loc_verts_index_combined)))
    else:
        print("No sessions were successfully processed.")


if __name__ == "__main__":
    main()