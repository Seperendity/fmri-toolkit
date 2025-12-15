import os
import nibabel as nib
from os.path import join, exists
import numpy as np
import glob

# 添加 FSL bin 路径到环境变量，便于调用fsl命令
os.environ["FSLDIR"] = "/usr/local/fsl"
os.environ["PATH"] += os.pathsep + "/usr/local/fsl/bin"
# 设置 FreeSurfer 安装目录，便于调用fsl命令
os.environ["FREESURFER_HOME"] = "/usr/local/freesurfer"
os.environ["PATH"] += os.pathsep + os.path.join(os.environ["FREESURFER_HOME"], "bin")
os.environ["SUBJECTS_DIR"] = "/usr/local/freesurfer/subjects/"

# 数据类型设置（'nsd'或'our'）
datatype = 'our'

# 根据数据类型设置核心路径（仅保留实际使用的路径）
if datatype == 'nsd':
    preproc_dir = '/opt/data/private/lq/data/NSD/nsddata_preprocdata/'
else:  # datatype == 'our'
    preproc_dir = '/cy/data/preprocdata/'

# 结果保存目录
results_dir = "/cy/data/betas_test/"

# 实验参数（仅保留必要配置）
subjects = ['sub-03']  # 被试列表
#sessions = ['ses-01']  # 会话列表
glm_types = ['betas_assumehrf', 'betas_fithrf', 'betas_fithrf_GLMdenoise_RR']
used_glm = glm_types[1]  # 直接指定当前使用的GLM类型（原glm_idx=2）
surf_fwhm = 4  # 表面平滑参数

def get_all_sessions(subject_path):
    session_pattern = join(subject_path, 'ses-*')
    session_dirs = sorted(glob.glob(session_pattern))
    return [os.path.basename(ses_dir) for ses_dir in session_dirs]

def get_all_runs(session_path):
    """获取session目录下的所有run文件夹"""
    run_pattern = join(session_path, 'run-*')
    run_dirs = sorted(glob.glob(run_pattern))
    return [os.path.basename(run_dir) for run_dir in run_dirs]

# 遍历被试（当前仅1个被试，保留循环结构以便扩展）
for sub in subjects:
    # 被试名称处理（NSD数据需特殊转换，其他数据保持原名）
    sub_recon_name = sub.replace("sub-", "sub")
    
    # 核心路径定义（仅保留实际使用的路径）
    roi_path = join(preproc_dir, sub, 'ROIs')  # ROI相关路径
    reg_output_path = join(preproc_dir, sub, 'reg_output')  # 配准结果路径
    
    # 配准文件路径
    reg_file = join(reg_output_path, 'register.dof6.dat')
    reg_file_downsampled = join(reg_output_path, 'register.downsampledfunc.dof6.dat')

    sessions = get_all_sessions(join(preproc_dir, sub))
    # 遍历每个session
    for session in sessions:
        session_path = join(preproc_dir, sub, session)

        ses_id_num = session.split('-')[1]
        if int(ses_id_num) < 20:
            continue
        
        # 获取当前session下所有的run
        runs = get_all_runs(session_path)
        print(f"Processing {sub} {session}, found runs: {runs}")
        
        # 生成表面mask（仅在第一个run时生成，后续run共用）
        mask_lh = join(roi_path, 'brain_mask.self.lh.nii.gz')
        mask_rh = join(roi_path, 'brain_mask.self.rh.nii.gz')
        
        if not exists(mask_lh) and runs:  # 如果mask不存在且有run文件夹
            first_run_path = join(session_path, runs[0])
            first_vol = join(first_run_path, 'f_skip_stc_mc_fm_1vol.nii.gz')
            
            if exists(first_vol):
                print(f"Generating surface masks using {runs[0]}...")
                
                # 左半球mask生成（调用FreeSurfer命令）
                cmd_lh = (f"mri_vol2surf --mov {first_vol} "
                          f"--reg {reg_file} --trgsubject {sub_recon_name} "
                          f"--interp nearest --projfrac 0.5 --hemi lh --o {mask_lh} --noreshape --cortex")
                os.system(cmd_lh)
                # 二值化mask
                os.system(f"mri_binarize --i {mask_lh} --min 0.00001 --o {mask_lh}")
                
                # 右半球mask生成
                cmd_rh = (f"mri_vol2surf --mov {first_vol} "
                          f"--reg {reg_file} --trgsubject {sub_recon_name} "
                          f"--interp nearest --projfrac 0.5 --hemi rh --o {mask_rh} --noreshape --cortex")
                os.system(cmd_rh)
                os.system(f"mri_binarize --i {mask_rh} --min 0.00001 --o {mask_rh}")
        
        # 遍历每个run
        for run in runs:
            run_func_path = join(session_path, run)
            print(f"Processing {sub} {session} {run}...")
            
            # 功能影像路径及加载
            func_downsample = join(run_func_path, 'f_downsampled.nii.gz')
            
            # 检查输入文件是否存在
            if not exists(func_downsample):
                print(f"Warning: {func_downsample} not found, skipping {run}")
                continue
                
            func_img = nib.load(func_downsample)
            func_affine = func_img.affine  # 功能影像的仿射矩阵
            
            # 表面数据路径
            surf_lh = join(run_func_path, 'lh.func.nativesurface.nii.gz')
            surf_rh = join(run_func_path, 'rh.func.nativesurface.nii.gz')
            fs_lh = join(run_func_path, 'lh.func.fs5.nii.gz')
            fs_rh = join(run_func_path, 'rh.func.fs5.nii.gz')
            fs_lh_sm = join(run_func_path, 'lh.func.fs5_sm.nii.gz')
            fs_rh_sm = join(run_func_path, 'rh.func.fs5_sm.nii.gz')
            
            # Volume转Surface（左半球）
            print(f"  Converting volume to surface (left hemisphere)...")
            cmd_vol2surf_lh = (f"mri_vol2surf --mov {func_downsample} "
                              f"--reg {reg_file_downsampled} --trgsubject {sub_recon_name} "
                              f"--interp trilin --projfrac 0.5 --hemi lh --o {surf_lh} --noreshape --cortex")
            os.system(cmd_vol2surf_lh)
            
            # Volume转Surface（右半球）
            print(f"  Converting volume to surface (right hemisphere)...")
            cmd_vol2surf_rh = (f"mri_vol2surf --mov {func_downsample} "
                              f"--reg {reg_file_downsampled} --trgsubject {sub_recon_name} "
                              f"--interp trilin --projfrac 0.5 --hemi rh --o {surf_rh} --noreshape --cortex")
            os.system(cmd_vol2surf_rh)
            
            # 表面数据转换到fsaverage5模板（左半球）
            print(f"  Converting to fsaverage5 template (left hemisphere)...")
            cmd_surf2fs_lh = (f"mri_surf2surf --sval {surf_lh} --srcsubject {sub_recon_name} "
                             f"--trgsubject fsaverage5 --hemi lh --tval {fs_lh} --cortex")
            os.system(cmd_surf2fs_lh)
            
            # 表面数据转换到fsaverage5模板（右半球）
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