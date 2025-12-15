import os
import json
import pydicom
import re
from datetime import datetime

def find_dicom_root_dir(base_dir, ses_id):
    """Find root folder starting with xxx-xx_xxx-{ses_id}"""
    prefix = f"xxx-xx_xxx-{ses_id}"
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path) and item.startswith(prefix):
            return item_path
    return None

def find_task_folders(root_dir):
    """Recursively find all subfolders containing task/TASK under root"""
    task_folders = []
    # Recursively traverse all subdirectories
    for dirpath, dirnames, _ in os.walk(root_dir):
        for dirname in dirnames:
            if "task" in dirname.lower():  # Match folders containing 'task' or 'TASK'
                # Extract task number (e.g., TASK1/task_02/Task-3 -> 1/2/3)
                num_match = re.search(r'(\d+)', dirname)
                task_num = int(num_match.group(1)) if num_match else None
                task_folders.append({
                    "path": os.path.join(dirpath, dirname),
                    "name": dirname,
                    "task_num": task_num
                })
    return task_folders

def sort_task_folders(task_folders):
    """Sort task folders by task number (numeric first, else name)"""
    # Group by whether number exists; sort numbered by value, others by name
    return sorted(
        task_folders,
        key=lambda x: (x["task_num"] is None, x["task_num"] if x["task_num"] is not None else x["name"])
    )

def get_dicom_acquisition_time(dicom_path):
    """
    Safely get DICOM acquisition time, handling missing fields or naming variants
    Tries multiple possible tag names (case variants and common alternatives)
    """
    try:
        try:
            # 1. Read DICOM file
            dicom_file = pydicom.dcmread(dicom_path)
            print(f"Successfully read file: {dicom_path}")

            # 2. Get SeriesDate and SeriesTime tag values
            series_date_str = dicom_file.get("SeriesDate")  # Tag (0008, 0021)
            series_time_str = dicom_file.get("SeriesTime")  # Tag (0008, 0031)
            # print(f"\nRaw values read from file:")
            # print(f"  - SeriesDate: {series_date_str}")
            # print(f"  - SeriesTime: {series_time_str}")

            # 3. Validate and parse timestamp
            if series_date_str and series_time_str:
                # Combine date and time strings
                # e.g., '20240820' + '144753.756000' -> '20240820144753.756000'
                datetime_str = f"{series_date_str}{series_time_str}"
                
                # Define datetime format for DICOM
                # %Y%m%d -> YYYYMMDD
                # %H%M%S -> HHMMSS
                # .%f -> fractional seconds
                # Note: strptime handles fractional digits
                dicom_format = "%Y%m%d%H%M%S.%f"

                try:
                    # 4. Parse the timestamp
                    parsed_datetime = datetime.strptime(datetime_str, dicom_format)
                    
                    # 5. Format result

                    return parsed_datetime.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]

                except ValueError:
                    print("\n[Error] Date or time string not in expected format; cannot parse.")

            else:
                # If any tag is missing, notify user
                print("\n[Info] Missing 'SeriesDate' or 'SeriesTime' tag; cannot build full timestamp.")


        except FileNotFoundError:
            print(f"[Error] File '{dicom_path}' not found. Please verify the path.")
        except Exception as e:
            print(f"[Error] Unexpected exception while reading/processing file: {e}")
    except Exception as e:
        print(f"⚠️ Failed to read DICOM file {os.path.basename(dicom_path)}: {str(e)}")
        return None

def extract_dicom_from_task_folder(task_folder):
    """Extract DICOM acquisition time from one task folder (handles missing tags)"""
    # Find DICOM files in folder (.dcm or tag-less DICOM)
    dicom_files = []
    for f in os.listdir(task_folder["path"]):
        f_path = os.path.join(task_folder["path"], f)
        if os.path.isfile(f_path): #and (f.endswith(".dcm") or "DICOM" in f.upper()):
            dicom_files.append(f_path)
    
    if not dicom_files:
        print(f"⚠️ No DICOM files found in task folder {task_folder['name']}")
        return None
    
    # Get acquisition time of the first DICOM file (try multiple fields)
    acq_time = get_dicom_acquisition_time(dicom_files[0])
    if acq_time is None:
        print(f"⚠️ Task folder {task_folder['name']} acquisition time unavailable")
        return None
    
    return {
        "task_num": task_folder["task_num"],
        "task_name": task_folder["name"],
        "acquisition_time": acq_time,
        "folder_path": task_folder["path"],
        "time_source": "DICOM tag" if not acq_time.isdigit() else "File modification time"
    }

def extract_bids_run_info(bids_func_dir):
    """Extract run info (run number, acquisition time) from BIDS func folder"""
    bids_run_info = []
    for json_filename in os.listdir(bids_func_dir):
        if json_filename.endswith("_bold.json") and "run-" in json_filename:
            json_path = os.path.join(bids_func_dir, json_filename)
            try:
                with open(json_path, "r") as f:
                    json_data = json.load(f)
                acq_time = json_data.get("AcquisitionTime")
                run_num = int(json_filename.split("run-")[1].split("_")[0])
                bids_run_info.append({
                    "run_num": run_num,
                    "acquisition_time": acq_time,
                    "filename": json_filename
                })
            except Exception as e:
                print(f"⚠️ Failed to read BIDS file {json_filename}: {str(e)}")
    
    # Sort by acquisition time
    return sorted(bids_run_info, key=lambda x: x["acquisition_time"])

def compare_single_session(dicom_root_dir, bids_func_dir, ses_id):
    """Compare task folder order with BIDS run order for one session"""
    print(f"\n{'='*80}")
    print(f"Checking session-{ses_id}")
    print(f"Root DICOM folder: {os.path.basename(dicom_root_dir)}")
    print(f"BIDS func path: {bids_func_dir}")
    print(f"{'-'*80}")

    # 1. Find and sort task folders
    task_folders = find_task_folders(dicom_root_dir)
    if not task_folders:
        print("❌ No subfolders containing task/TASK found")
        return False
    sorted_tasks = sort_task_folders(task_folders)
    print(f"Found {len(sorted_tasks)} task folders, sorted:")
    for i, task in enumerate(sorted_tasks, 1):
        print(f"  {i}. {task['name']} (number: {task['task_num'] or 'none'})")

    # 2. Extract DICOM acquisition time for task folders
    task_dicom_info = [extract_dicom_from_task_folder(t) for t in sorted_tasks]
    task_dicom_info = [t for t in task_dicom_info if t is not None]  # Filter out invalid items
    if not task_dicom_info:
        print("❌ Failed to extract any valid DICOM acquisition time")
        return False

    # 3. Extract BIDS run info
    bids_runs = extract_bids_run_info(bids_func_dir)
    if not bids_runs:
        print("❌ No valid BIDS run files found")
        return False

    # # 4. Compare order (match task number with run number)
    # print(f"\n{'-'*80}")
    # print(f"{'Original Order':<10} {'Task Folder':<25} {'Task Num':<10} {'BIDS Run Num':<15} Match Status")
    # print(f"{'-'*80}")
    
    all_match = True
    min_count = min(len(task_dicom_info), len(bids_runs))
    for i in range(min_count):
        task = task_dicom_info[i]
        run = bids_runs[i]
        # If task has no number, only check if the order is consistent (do not verify number match)
        match = (task["task_num"] == run["run_num"]) if task["task_num"] is not None else True
        print(task['task_num'], task['acquisition_time'])
        print(run["run_num"], run['acquisition_time'])
        print('----------------------------------------')

def batch_compare_sessions(base_dicom_dir, base_bids_dir, sub_id, ses_ids):
    """Batch-compare multiple sessions"""
    print(f"Starting batch check: subject {sub_id}, total {len(ses_ids)} sessions")
    results = []
    
    for ses_id in ses_ids:
        # Locate root DICOM folder
        dicom_root = find_dicom_root_dir(base_dicom_dir, ses_id)
        if not dicom_root:
            print(f"❌ Skip session-{ses_id}: no root folder starting with YXMLQ-FAN_HUA_001-{ses_id}")
            continue
        
        # Build BIDS path
        bids_func_dir = os.path.join(base_bids_dir, f"sub-{sub_id}", f"ses-{ses_id:02d}", "func")
        if not os.path.exists(bids_func_dir):
            print(f"❌ Skip session-{ses_id}: BIDS path not found {bids_func_dir}")
            continue
        
        # Compare current session
        result = compare_single_session(dicom_root, bids_func_dir, ses_id)
        results.append({
            "session": ses_id,
            "dicom_root": os.path.basename(dicom_root),
            "status": "match" if result else "mismatch"
        })
    
    # # Summary of results
    # print(f"\n{'='*80}")
    # print("Batch Check Summary")
    # print(f"{'-'*80}")
    # print(f"{'session':<10} {'Main DICOM Folder':<30} Status")
    # print(f"{'-'*80}")
    # for res in results:
    #     print(f"{res['session']:<10} {res['dicom_root']:<30} {res['status']}")
    # print(f"{'='*80}")

# ------------------- Config -------------------
if __name__ == "__main__":
    BASE_DICOM_DIR = "/share/home/xx/data/our_nsd/subj_xx"  # DICOM root directory
    BASE_BIDS_DIR = "/share/home/xx/data/our_nsd/subj_xx/bids"  # BIDS root directory
    SUBJECT_ID = "01"  # Subject ID (01 in sub-01)
    SESSION_IDS = [i for i in range(1, 10)]  # List of sessions to check
    
    batch_compare_sessions(
        base_dicom_dir=BASE_DICOM_DIR,
        base_bids_dir=BASE_BIDS_DIR,
        sub_id=SUBJECT_ID,
        ses_ids=SESSION_IDS
    )