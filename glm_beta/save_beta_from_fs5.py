"""
Beta Value Processing Script

DESCRIPTION:
This script processes fMRI beta coefficient data.
It loads beta values from GLMdenoise analysis, applies spatial masking using the NSD 
general mask, removes specific stimuli at regular intervals, 
and saves the processed data for further analysis.

FUNCTIONALITY:
1. Loads left and right hemisphere beta coefficient files (.nii.gz format)
2. Concatenates hemispheric data into a unified brain representation
3. Applies NSD general cortical mask to focus on visually responsive regions
4. Removes stimuli at specified intervals (default: every 20 stimulus)
5. Saves processed data and metadata for each session

INPUT FILES:
- Beta coefficient files: 
  {data_dir}/sub-{subject}/fsaverage/betas_fithrf_GLMdenoise_RR/lh.betas_session{session}.fs.nii.gz
  {data_dir}/sub-{subject}/fsaverage/betas_fithrf_GLMdenoise_RR/rh.betas_session{session}.fs.nii.gz
- NSD general mask files:
  /usr/local/freesurfer/subjects/fsaverage/label/lh.nsdgeneral.mgz
  /usr/local/freesurfer/subjects/fsaverage/label/rh.nsdgeneral.mgz

OUTPUT FILES:
- Processed beta data: {output_dir}/sub-{subject}/session{session}_betas.npy
- Session metadata: {output_dir}/sub-{subject}/session{session}_metadata.npy

USAGE EXAMPLES:
# Process single subject, multiple sessions
python save_beta.py --sub 01 --ses 01-05 --space fsaverage
"""

import argparse
from os.path import join, exists
import nibabel as nib
import numpy as np
import os

def parse_session_ranges(session_args):
    """
    解析session参数，支持范围输入
    
    参数:
    session_args: 从命令行传入的session参数列表
    
    返回:
    sessions: 解析后的session列表（字符串格式，保持原有的零填充）
    
    示例:
    ['01', '02', '05-08'] -> ['01', '02', '05', '06', '07', '08']
    """
    sessions = []
    
    for arg in session_args:
        if '-' in arg and len(arg.split('-')) == 2:
            # 处理范围输入，如 '01-05'
            start_str, end_str = arg.split('-')
            try:
                start = int(start_str)
                end = int(end_str)
                
                # 确保范围有效
                if start > end:
                    raise ValueError(f"Invalid range: {arg}. Start ({start}) must be <= end ({end})")
                
                # 生成范围内的所有session，保持原有的零填充格式
                width = len(start_str)  # 保持原有的位数
                for i in range(start, end + 1):
                    sessions.append(str(i).zfill(width))
                    
            except ValueError as e:
                if "invalid literal" in str(e):
                    raise ValueError(f"Invalid session range format: {arg}. Use format like '01-05'")
                else:
                    raise e
        else:
            # 处理单个session
            sessions.append(arg)
    
    # 去重并排序（按数值排序）
    unique_sessions = list(set(sessions))
    unique_sessions.sort(key=int)
    
    return unique_sessions

def create_removal_mask(num_runs, stimuli_per_run=105, remove_every_n=21):
    """
    创建删除每20个后按键图片的掩码
    
    参数:
    num_runs: run的数量
    stimuli_per_run: 每个run的刺激数量，默认105
    remove_every_n: 每多少个删除一个，默认21（即删除第20个后的按键图片）
    
    返回:
    mask: 布尔数组，True表示保留，False表示删除
    """
    total_stimuli = num_runs * stimuli_per_run
    mask = np.ones(total_stimuli, dtype=bool)
    
    # 计算需要删除的索引位置
    indices_to_remove = []
    for run_idx in range(num_runs):
        run_start = run_idx * stimuli_per_run
        # 每个run内删除每第20个后的按键图片（索引为20, 41, 62, 83, 104）
        run_remove_indices = np.arange(run_start + 20, run_start + stimuli_per_run, remove_every_n)
        indices_to_remove.extend(run_remove_indices)
    
    indices_to_remove = np.array(indices_to_remove)
    print(f"将要被删除的索引位置: {indices_to_remove}")
    print(f"每个run删除 {len(indices_to_remove) // num_runs} 个刺激")
    
    mask[indices_to_remove] = False
    return mask


def load_and_process_session_data(data_dir, sub, ses, space, stimuli_per_run=105, remove_every_n=21):
    """
    加载并处理单个session的数据
    
    参数:
    data_dir: 数据根目录
    sub: 被试编号
    ses: session编号
    stimuli_per_run: 每个run的刺激数量
    remove_every_n: 删除间隔
    
    返回:
    processed_data: 处理后的数据 (vertices, stimuli)
    num_runs: 该session的run数量
    """
    # 定义首选和备用目录名

    if space == 'fsaverage5':
        fs_suffix = 'fs5'
    else:
        fs_suffix = 'fs'    

    # 1. 尝试使用首选路径构建文件路径
    lh_file = join(data_dir, f'sub-{sub}',
                   f'lh.betas_session{ses}.{fs_suffix}_sm.nii.gz')
    rh_file = join(data_dir, f'sub-{sub}', 
                   f'rh.betas_session{ses}.{fs_suffix}_sm.nii.gz')
    
            
    print(f"Loading data for sub-{sub} session{ses}...")
    print(f"  -> Found data in: {os.path.dirname(lh_file)}")
    
    # 加载数据
    lh = nib.load(lh_file).get_fdata().squeeze()
    rh = nib.load(rh_file).get_fdata().squeeze()
    
    # 合并左右半球数据
    data = np.concatenate((lh, rh), axis=0)
    
    print(f"Original data shape: {data.shape}")
    
    # 计算run数量
    num_cols = data.shape[1]
    num_runs = num_cols // stimuli_per_run
    
    if num_cols != num_runs * stimuli_per_run:
        print(f"Warning: Total stimuli ({num_cols}) not divisible by stimuli_per_run ({stimuli_per_run})")
        print(f"Assuming {num_runs} complete runs")
    
    print(f"Detected {num_runs} runs with {stimuli_per_run} stimuli each")
    
    # 创建删除掩码
    removal_mask = create_removal_mask(num_runs, stimuli_per_run, remove_every_n)
    
    # 应用掩码删除按键图片
    processed_data = data[:, removal_mask]
    
    print(f"After removal, data shape: {processed_data.shape}")
    print(f"Removed {num_cols - processed_data.shape[1]} stimuli")
    
    return processed_data, num_runs


def load_nsd_mask(space):
    """
    加载NSD general掩码
    
    返回:
    loc_verts_index_lh: 左半球掩码索引
    loc_verts_index_rh: 右半球掩码索引
    loc_verts_index_combined: 合并的掩码索引
    n_verts_lh: 左半球顶点数
    n_verts_rh: 右半球顶点数
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
    
    # 获取nsdgeneral掩码内的顶点索引
    loc_verts_index_lh = np.where(lh_mask_data != 0)[0]
    loc_verts_index_rh = np.where(rh_mask_data != 0)[0]
    
    # 合并左右半球的掩码顶点索引
    loc_verts_index_combined = np.concatenate([
        loc_verts_index_lh,  # 左半球索引
        loc_verts_index_rh + n_verts_lh  # 右半球索引（加上偏移量）
    ])
    
    print(f"NSD mask info({space} space):")
    print(f"  Left hemisphere vertices: {n_verts_lh}, masked: {len(loc_verts_index_lh)}")
    print(f"  Right hemisphere vertices: {n_verts_rh}, masked: {len(loc_verts_index_rh)}")
    print(f"  Total masked vertices: {len(loc_verts_index_combined)}")
    
    return (loc_verts_index_lh, loc_verts_index_rh, loc_verts_index_combined, 
            n_verts_lh, n_verts_rh)


def save_masked_betas(data, loc_verts_index_combined, save_dir, sub, ses, space, start_stim_idx=0):
    """
    保存经过掩码处理的beta值（整个session保存为一个文件）
    
    参数:
    data: 原始数据 (vertices, stimuli)
    loc_verts_index_combined: 掩码索引
    save_dir: 保存目录
    sub: 被试编号
    ses: session编号
    start_stim_idx: 起始刺激索引（用于多session连续编号）
    
    返回:
    next_stim_idx: 下一个刺激索引
    """
    # 应用nsdgeneral掩码
    masked_data = data[loc_verts_index_combined, :]
    print(f"Masked data shape: {masked_data.shape}")
    
    # 创建被试的保存目录
    space_save_dir = join(save_dir, f'sub-{sub}')
    subject_save_dir = join(space_save_dir, space)
    os.makedirs(subject_save_dir, exist_ok=True)
    
    # 保存整个session的beta值到一个文件
    num_stimuli = masked_data.shape[1]
    
    filename = f'session{ses}_betas.npy'
    filepath = join(subject_save_dir, filename)
    np.save(filepath, masked_data)
    
    # 保存元数据信息到单独的文件
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

# 后续读取数据
def load_session_betas(save_dir, sub, ses, space):
    """
    加载指定session的beta数据
    
    参数:
    save_dir: 保存目录
    sub: 被试编号
    ses: session编号
    
    返回:
    betas: beta数据 (vertices, stimuli)
    metadata: 元数据字典
    """
    space_save_dir = join(save_dir, space)
    subject_save_dir = join(space_save_dir, f'sub-{sub}')
    
    # 加载beta数据
    betas_file = join(subject_save_dir, f'session{ses}_betas.npy')
    betas = np.load(betas_file)
    
    # 加载元数据
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

    # 解析session参数
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
    
    # 加载NSD掩码（所有被试和session共用）
    try:
        (loc_verts_index_lh, loc_verts_index_rh, loc_verts_index_combined, 
         n_verts_lh, n_verts_rh) = load_nsd_mask(args.space)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    
    # 创建主输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 处理指定被试的所有session
    total_stimuli_processed = 0
    session_info = {}
    processed_sessions = []
    skipped_sessions = []

    for ses in sessions:
        print(f"\n--- Processing Subject {args.sub}, Session {ses} ---")
        
        # 加载和处理session数据
        processed_data, num_runs = load_and_process_session_data(
            args.data_dir, args.sub, ses, args.space, args.stimuli_per_run, args.remove_every_n
        )
        
        if processed_data is None:
            print(f"Skipping session {ses} due to missing files")
            skipped_sessions.append(ses)
            continue
        
        # 保存masked beta值
        next_stim_idx = save_masked_betas(
            processed_data, loc_verts_index_combined, args.output_dir, 
            args.sub, ses, args.space, total_stimuli_processed
        )
        
        # 记录session信息
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