#!/usr/bin/env python
# coding: utf-8
import os
import subprocess
import logging
import argparse
from pathlib import Path
from datetime import datetime

def setup_logging(output_dir, subject_id):
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(f'{output_dir}/batch_{subject_id}_conversion.log'),
            logging.StreamHandler()
        ]
    )

def check_dicom_files(data_dir):
    """Check if directory contains DICOM files"""
    dicom_extensions = ['.dcm', '.IMA', '.ima', '']
    for root, dirs, files in os.walk(data_dir):
        for file in files:
            if any(file.lower().endswith(ext) for ext in dicom_extensions):
                return True
    return False

def run_dcm2bids(data_dir, subject_id, session_id, config_file, output_dir):
    """Run dcm2bids conversion"""
    cmd = [
        'dcm2bids',
        '-d', data_dir,
        '-p', subject_id,
        '-s', session_id,
        '-c', config_file,
        '-o', output_dir,
        '--auto_extract_entities',
        '--clobber'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.info(f"Session {session_id} conversion successful")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Session {session_id} conversion failed: {e}")
        logging.error(f"Error output: {e.stderr}")
        return False

def main():
    """Main function to process all sessions"""
    parser = argparse.ArgumentParser(description='Batch convert DICOM to BIDS format')
    parser.add_argument('-s', '--subject', default='04', help='Subject ID')

    parser.add_argument('-b', '--base-dir', default='/share/home/user/data/our_nsd/subj_xx', 
                       help='Base directory containing session folders')

    parser.add_argument('-c', '--config', default='/share/home/user/data/our_nsd/config.json', 
                       help='Path to dcm2bids config file')
    parser.add_argument('-o', '--output', default='/share/home/user/data/our_nsd/bids/', 
                       help='Output BIDS directory')
    parser.add_argument('--start-session', type=int, default=23, help='Start session number')
    parser.add_argument('--end-session', type=int, default=40, help='End session number')

    parser.add_argument('--session-pattern', default='YXMLQ-WANG_XIN_LIANG004-{}',
                help='Session directory pattern')
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output).mkdir(parents=True, exist_ok=True)

    # Setup logging
    setup_logging(args.output, args.subject)
    
    logging.info(f"Starting batch conversion for subject {args.subject}")
    logging.info(f"Base directory: {args.base_dir}")
    logging.info(f"Config file: {args.config}")
    logging.info(f"Output directory: {args.output}")
    
    successful_sessions = []
    failed_sessions = []
    
    # Process all sessions
    for session in range(args.start_session, args.end_session):
        session_id = f"{session:02d}"  # Format as 02d
        
        # Dynamically find session directory
        session_pattern = args.session_pattern.format(session)
        session_dir = None
        
        # Search for matching directory in base_dir
        for dir_name in os.listdir(args.base_dir):
            if (dir_name.startswith(session_pattern) and 
                not dir_name.endswith('.zip') and 
                os.path.isdir(os.path.join(args.base_dir, dir_name))):
                session_dir = os.path.join(args.base_dir, dir_name)
                break
        
        # Check if session directory was found
        if session_dir is None or not os.path.exists(session_dir):
            logging.warning(f"Session {session_id} directory not found, pattern: {session_pattern}")
            continue
        
        logging.info(f"Session {session_id} found directory: {session_dir}")
        
        # Automatically find subdirectories (usually only one)
        subdirs = [d for d in os.listdir(session_dir) if os.path.isdir(os.path.join(session_dir, d))]
        if not subdirs:
            logging.warning(f"Session {session_id} no subdirectories found")
            continue
        
        # Use first subdirectory as data_dir
        data_dir = os.path.join(session_dir, subdirs[0])
        logging.info(f"Session {session_id} using data directory: {data_dir}")
        
        # Check if data directory exists
        if not os.path.exists(data_dir):
            logging.warning(f"Session {session_id} data directory does not exist: {data_dir}")
            continue
        
        # Check for DICOM files
        if not check_dicom_files(data_dir):
            logging.warning(f"Session {session_id} no DICOM files found")
            continue
        
        logging.info(f"Processing Session {session_id}...")
        
        # Run dcm2bids
        if run_dcm2bids(data_dir, args.subject, session_id, args.config, args.output):
            successful_sessions.append(session_id)
        else:
            failed_sessions.append(session_id)
    
    # Output summary
    logging.info("=" * 50)
    logging.info("Batch processing completed")
    logging.info(f"Successfully converted sessions: {successful_sessions}")
    logging.info(f"Failed sessions: {failed_sessions}")

if __name__ == "__main__":
    main()