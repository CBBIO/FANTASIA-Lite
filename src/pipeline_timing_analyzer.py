#!/usr/bin/env python3
"""Benchmark the Lite pipeline across one or more FASTA files."""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def check_csv_compatibility(csv_path: Path, expected_fieldnames: List[str]) -> Tuple[bool, List[str]]:
    """
    Check if existing CSV has compatible columns.
    Returns (is_compatible, existing_fieldnames)
    """
    if not csv_path.exists():
        return True, []
    
    try:
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            existing_fieldnames = reader.fieldnames or []
            
        # Check if all expected fields are present
        missing_fields = set(expected_fieldnames) - set(existing_fieldnames)
        extra_fields = set(existing_fieldnames) - set(expected_fieldnames)
        
        if missing_fields or extra_fields:
            return False, existing_fieldnames
        
        return True, existing_fieldnames
    except Exception as e:
        print(f"Warning: Could not read existing CSV file: {e}")
        return False, []


def handle_existing_report_file(report_csv: Path, fieldnames: List[str]) -> str:
    """
    Handle existing report file. Returns mode: 'w' for overwrite, 'a' for append, 'abort' to stop.
    """
    if not report_csv.exists():
        return 'w'
    
    print(f"Report file {report_csv} already exists.")
    
    # Check CSV compatibility
    is_compatible, existing_fieldnames = check_csv_compatibility(report_csv, fieldnames)
    
    if not is_compatible:
        print(f"ERROR: Existing CSV file has incompatible columns!")
        print(f"Expected columns: {', '.join(fieldnames)}")
        print(f"Existing columns: {', '.join(existing_fieldnames)}")
        print(f"Please choose a different report filename or fix the existing file.")
        return 'abort'
    
    # CSV is compatible, ask user what to do
    while True:
        choice = input("What would you like to do?\n"
                      "  [o] Overwrite the existing file\n"
                      "  [a] Append to the existing file\n"
                      "  [c] Cancel and exit\n"
                      "Choice (o/a/c): ").strip().lower()
        
        if choice in ['o', 'overwrite']:
            return 'w'
        elif choice in ['a', 'append']:
            return 'a'
        elif choice in ['c', 'cancel']:
            return 'abort'
        else:
            print("Please enter 'o', 'a', or 'c'")


def get_gpu_info() -> Dict[str, str]:
    """
    Get GPU information using nvidia-smi if available.
    Returns dictionary with GPU details or empty dict if no GPU.
    """
    gpu_info = {
        'device_type': 'cpu',
        'gpu_model': '',
        'gpu_memory_total': '',
        'gpu_memory_used': '',
        'gpu_memory_free': ''
    }
    
    # Check if nvidia-smi is available
    nvidia_smi = shutil.which('nvidia-smi')
    if not nvidia_smi:
        return gpu_info
    
    try:
        # Get GPU info in JSON format
        cmd = [
            'nvidia-smi', 
            '--query-gpu=name,memory.total,memory.used,memory.free',
            '--format=csv,noheader,nounits'
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            # Use first GPU if multiple
            first_gpu = lines[0].split(', ')
            if len(first_gpu) >= 4:
                gpu_info.update({
                    'device_type': 'gpu',
                    'gpu_model': first_gpu[0].strip(),
                    'gpu_memory_total': f"{first_gpu[1].strip()} MB",
                    'gpu_memory_used': f"{first_gpu[2].strip()} MB", 
                    'gpu_memory_free': f"{first_gpu[3].strip()} MB"
                })
    
    except (subprocess.CalledProcessError, IndexError, ValueError) as e:
        print(f"Warning: Could not get GPU info: {e}")
    
    return gpu_info


def count_csv_rows(csv_path: Path) -> int:
    """Count rows in CSV file (excluding header)."""
    if not csv_path.exists():
        return 0
    
    try:
        with open(csv_path, 'r', encoding='utf-8') as f:
            # Count lines and subtract header
            line_count = sum(1 for line in f)
            return max(0, line_count - 1)
    except Exception as e:
        print(f"Warning: Could not count rows in {csv_path}: {e}")
        return 0


def count_unique_sequences_in_results(csv_path: Path) -> int:
    """Count unique protein sequences processed (by unique query_accession)."""
    if not csv_path.exists():
        return 0
    
    try:
        import csv
        unique_sequences = set()
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'query_accession' in row and row['query_accession']:
                    unique_sequences.add(row['query_accession'])
        return len(unique_sequences)
    except Exception as e:
        print(f"Warning: Could not count unique sequences in {csv_path}: {e}")
        return 0


def read_run_metadata(metadata_path: Path) -> Dict[str, object]:
    """Read pipeline run metadata written by fantasia_pipeline.py."""
    if not metadata_path.exists():
        return {}

    try:
        import yaml

        with metadata_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
        return data if isinstance(data, dict) else {}
    except Exception:
        try:
            return json.loads(metadata_path.read_text(encoding="utf-8"))
        except Exception:
            return {}


def find_fasta_files(directory: Path) -> List[Path]:
    """Find all FASTA files in directory with common extensions."""
    fasta_extensions = {'.fa', '.faa', '.fasta'}
    fasta_files = []
    
    if not directory.exists():
        print(f"Error: Directory {directory} does not exist")
        return fasta_files
    
    for file_path in directory.iterdir():
        if file_path.is_file() and file_path.suffix.lower() in fasta_extensions:
            fasta_files.append(file_path)
    
    return sorted(fasta_files)


def count_fasta_sequences(fasta_path: Path) -> int:
    """Count number of sequences in FASTA file."""
    if not fasta_path.exists():
        return 0
    
    try:
        count = 0
        with open(fasta_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    except Exception as e:
        print(f"Warning: Could not count sequences in {fasta_path}: {e}")
        return 0


def run_fantasia_pipeline(
    fasta_path: Path,
    pipeline_script: Path,
    model_name: str
) -> Tuple[bool, float, str, Optional[Path], Dict[str, object]]:
    """
    Run FANTASIA pipeline on a single FASTA file.
    
    Returns:
        - Success status (bool)
        - Runtime in seconds (float)  
        - Error message if failed (str)
        - Path to the created pipeline output directory (Path or None)
    """
    start_time = time.time()
    
    # Record the current timestamp to find the pipeline output directory later
    pipeline_start_time = datetime.now()
    
    # Prepare command - let pipeline use its default timestamped outputs
    python_bin = shutil.which("python3") or sys.executable or "python3"
    cmd = [
        python_bin,
        str(pipeline_script),
        '--serial-models',
        '--embed-models', model_name,
        str(fasta_path)
    ]
    
    try:
        print(f"  Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            check=True,
            cwd=pipeline_script.parent
        )
        
        runtime = time.time() - start_time
        
        # Find the pipeline output directory created by fantasia_pipeline.py
        # It should be named outputs_YYYYMMDD_HHMMSS
        pipeline_output_dir = None
        base_dir = pipeline_script.parent
        
        # Look for the most recent directory that matches our pattern
        # and was created during our pipeline run (with some buffer time)
        pipeline_end_time = datetime.now()
        buffer_minutes = 2  # Allow 2 minutes buffer for timing differences
        
        best_match = None
        for item in base_dir.iterdir():
            if item.is_dir() and item.name.startswith('outputs_'):
                try:
                    # Extract timestamp from directory name (outputs_YYYYMMDD_HHMMSS)
                    timestamp_str = item.name[8:]  # Remove 'outputs_' prefix
                    dir_time = datetime.strptime(timestamp_str, '%Y%m%d_%H%M%S')
                    
                    # Check if this directory was created during our pipeline run
                    # (between start time minus buffer and end time plus buffer)
                    time_before = pipeline_start_time - timedelta(minutes=buffer_minutes)
                    time_after = pipeline_end_time + timedelta(minutes=buffer_minutes)
                    
                    if time_before <= dir_time <= time_after:
                        # If we haven't found a match yet, or this one is more recent
                        if best_match is None or dir_time > best_match[1]:
                            best_match = (item, dir_time)
                            
                except ValueError:
                    # Skip directories that don't match the timestamp format
                    continue
        
        if best_match:
            pipeline_output_dir = best_match[0]
            
        if pipeline_output_dir and pipeline_output_dir.exists():
            print(f"  Pipeline output saved to: {pipeline_output_dir.name}")
            metadata = read_run_metadata(pipeline_output_dir / "run_metadata.yaml")
            return True, runtime, "", pipeline_output_dir, metadata
        else:
            # Try to find any recent outputs directory as fallback
            recent_dirs = []
            for item in base_dir.iterdir():
                if item.is_dir() and item.name.startswith('outputs_'):
                    try:
                        timestamp_str = item.name[8:]
                        dir_time = datetime.strptime(timestamp_str, '%Y%m%d_%H%M%S')
                        recent_dirs.append((item, dir_time))
                    except ValueError:
                        continue
            
            if recent_dirs:
                # Get the most recent one
                recent_dirs.sort(key=lambda x: x[1], reverse=True)
                pipeline_output_dir = recent_dirs[0][0]
                print(f"  Using most recent output directory: {pipeline_output_dir.name}")
                metadata = read_run_metadata(pipeline_output_dir / "run_metadata.yaml")
                return True, runtime, "", pipeline_output_dir, metadata
            else:
                return True, runtime, "Warning: Could not find pipeline output directory", None, {}
        
    except subprocess.CalledProcessError as e:
        runtime = time.time() - start_time
        error_msg = f"Exit code {e.returncode}: {e.stderr[:200] if e.stderr else e.stdout[:200]}"
        return False, runtime, error_msg, None, {}


def process_all_fasta_files(
    fasta_dir: Path,
    pipeline_script: Path,
    model_names: List[str],
    report_csv: Path,
    specific_files: Optional[List[Path]] = None
) -> None:
    """Process all FASTA files and generate comprehensive report."""
    
    # Get system info once
    gpu_info = get_gpu_info()
    
    print(f"System info:")
    print(f"  Device type: {gpu_info['device_type']}")
    if gpu_info['device_type'] == 'gpu':
        print(f"  GPU model: {gpu_info['gpu_model']}")
        print(f"  GPU memory: {gpu_info['gpu_memory_total']}")
    print(f"  Models: {', '.join(model_names)}")
    print()
    
    # Get FASTA files to process
    if specific_files:
        fasta_files = specific_files
        print(f"Processing {len(fasta_files)} specified files:")
        for f in fasta_files:
            print(f"  - {f.name}")
    else:
        fasta_files = find_fasta_files(fasta_dir)
        if not fasta_files:
            print(f"No FASTA files found in {fasta_dir}")
            return
        print(f"Found {len(fasta_files)} FASTA files:")
        for f in fasta_files:
            print(f"  - {f.name}")
    print()
    
    # Prepare CSV report
    fieldnames = [
        'gpu_model',
        'gpu_memory_total',
        'model_name',
        'runtime_seconds',
        'embedding_time_seconds',
        'lookup_time_seconds',
        'postprocess_time_seconds',
        'measured_stage_time_seconds',
        'input_sequences',
        'processing_rate_seq_per_sec',
        'sequences_processed',
        'sequences_failed',
        'total_results',
        'total_annotations',
        'success',
        'device_type',
        'gpu_memory_used_before',
        'gpu_memory_free_before',
        'fasta_file',
        'pipeline_output_dir',
        'timestamp',
        'error_message'
    ]
    
    # Handle existing report file
    file_mode = handle_existing_report_file(report_csv, fieldnames)
    if file_mode == 'abort':
        print("Operation cancelled.")
        return
    
    # Create or append to report CSV
    with open(report_csv, file_mode, newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        # Write header only if we're overwriting or creating new file
        if file_mode == 'w':
            writer.writeheader()
        
        # Process each FASTA file with each model
        for i, fasta_path in enumerate(fasta_files, 1):
            print(f"Processing {i}/{len(fasta_files)}: {fasta_path.name}")
            
            # Count input sequences
            input_seqs = count_fasta_sequences(fasta_path)
            print(f"  Input sequences: {input_seqs}")
            
            # Process with each model
            for model_name in model_names:
                print(f"  Running with model: {model_name}")
                
                # Get GPU info before processing (in case it changes)
                current_gpu_info = get_gpu_info()
                
                # Run pipeline - let it create its own timestamped directory
                success, runtime, error_msg, pipeline_output_path, run_metadata = run_fantasia_pipeline(
                    fasta_path, pipeline_script, model_name
                )
                stage_timing = run_metadata.get('stage_timing_seconds', {}) if isinstance(run_metadata, dict) else {}
                embedding_time = float(stage_timing.get('embedding', 0.0) or 0.0)
                lookup_time = float(stage_timing.get('lookup', 0.0) or 0.0)
                postprocess_time = float(stage_timing.get('postprocess', 0.0) or 0.0)
                measured_stage_total = float(stage_timing.get('measured_total', 0.0) or 0.0)
                
                # Count results from pipeline output directory
                if pipeline_output_path and success:
                    results_file = pipeline_output_path / 'results.csv'
                    failures_file = pipeline_output_path / 'failed_sequences.csv'
                    
                    # Count unique protein sequences processed, not GO annotations
                    sequences_processed = count_unique_sequences_in_results(results_file)
                    sequences_failed = count_csv_rows(failures_file)
                    
                    # Also track total annotations found
                    total_annotations = count_csv_rows(results_file)
                else:
                    sequences_processed = 0
                    sequences_failed = 0
                    total_annotations = 0 
                
                # Calculate processing rate
                processing_rate = input_seqs / runtime if runtime > 0 else 0
                
                print(f"    Success: {success}")
                print(f"    Runtime: {runtime:.2f}s")
                if measured_stage_total > 0:
                    print(
                        "    Stage timing: "
                        f"embedding={embedding_time:.2f}s, "
                        f"lookup={lookup_time:.2f}s, "
                        f"postprocess={postprocess_time:.2f}s"
                    )
                print(f"    Processed: {sequences_processed} sequences")
                print(f"    Failed: {sequences_failed} sequences")
                print(f"    Total annotations: {total_annotations}")
                print(f"    Rate: {processing_rate:.2f} seq/s")
                if pipeline_output_path:
                    print(f"    Pipeline output: {pipeline_output_path.name}")
                if error_msg:
                    print(f"    Error: {error_msg}")
                
                # Write to CSV
                row = {
                    'timestamp': datetime.now().isoformat(),
                    'fasta_file': fasta_path.name,
                    'input_sequences': input_seqs,
                    'device_type': current_gpu_info['device_type'],
                    'gpu_model': current_gpu_info['gpu_model'],
                    'gpu_memory_total': current_gpu_info['gpu_memory_total'],
                    'gpu_memory_used_before': current_gpu_info['gpu_memory_used'],
                    'gpu_memory_free_before': current_gpu_info['gpu_memory_free'],
                    'model_name': model_name,
                    'success': success,
                    'runtime_seconds': f"{runtime:.2f}",
                    'embedding_time_seconds': f"{embedding_time:.2f}",
                    'lookup_time_seconds': f"{lookup_time:.2f}",
                    'postprocess_time_seconds': f"{postprocess_time:.2f}",
                    'measured_stage_time_seconds': f"{measured_stage_total:.2f}",
                    'sequences_processed': sequences_processed,
                    'sequences_failed': sequences_failed,
                    'total_results': sequences_processed + sequences_failed,
                    'total_annotations': total_annotations,
                    'processing_rate_seq_per_sec': f"{processing_rate:.2f}",
                    'pipeline_output_dir': pipeline_output_path.name if pipeline_output_path else "",
                    'error_message': error_msg
                }
                
                writer.writerow(row)
                csvfile.flush()  # Ensure data is written immediately
            
            print()  # Extra line break between files
    
    print(f"Timing analysis complete! Results saved to {report_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze FANTASIA pipeline timing and performance across multiple FASTA files"
    )
    parser.add_argument(
        "--fasta-dir",
        default="fasta_test",
        help="Directory containing FASTA files to process (default: fasta_test)"
    )
    parser.add_argument(
        "--pipeline-script",
        default="src/fantasia_pipeline.py",
        help="Path to fantasia_pipeline.py script (default: src/fantasia_pipeline.py)"
    )
    parser.add_argument(
        "--model",
        default=['prot_t5'],
        choices=['prot_t5', 'ankh3'],
        nargs='*',
        help="Embedding model(s) to use (default: prot_t5)"
    )
    parser.add_argument(
        "--files",
        nargs='*',
        help="Specific FASTA files to process (if not provided, processes all files in --fasta-dir)"
    )
    parser.add_argument(
        "--report-csv",
        default="pipeline_timing_analysis.csv",
        help="CSV file for timing analysis report (default: pipeline_timing_analysis.csv)"
    )
    
    args = parser.parse_args()
    
    # Convert paths
    fasta_dir = Path(args.fasta_dir)
    pipeline_script = Path(args.pipeline_script)
    report_csv = Path(args.report_csv)
    
    # Handle file specification
    if args.files:
        # Specific files provided
        fasta_files = [Path(f).resolve() for f in args.files]
        # Validate that all files exist
        missing_files = [f for f in fasta_files if not f.exists()]
        if missing_files:
            print(f"Error: The following files do not exist:")
            for f in missing_files:
                print(f"  - {f}")
            return 1
        print(f"FANTASIA Pipeline Timing Analyzer")
        print(f"==================================")
        print(f"Files to process: {len(fasta_files)} specified files")
        for f in fasta_files:
            print(f"  - {f.name}")
    else:
        # Use directory mode (original behavior)
        if not fasta_dir.exists():
            print(f"Error: FASTA directory {fasta_dir} does not exist")
            return 1
        fasta_files = None  # Will be determined in process_all_fasta_files
        print(f"FANTASIA Pipeline Timing Analyzer")
        print(f"==================================")
        print(f"FASTA directory: {fasta_dir.absolute()}")
    
    # Validate pipeline script
    if not pipeline_script.exists():
        print(f"Error: Pipeline script {pipeline_script} does not exist")
        return 1
    
    print(f"Pipeline script: {pipeline_script.absolute()}")
    print(f"Models: {', '.join(args.model)}")
    print(f"Analysis report: {report_csv.absolute()}")
    print()
    
    # Run timing analysis
    process_all_fasta_files(
        fasta_dir, pipeline_script, args.model, report_csv, specific_files=fasta_files
    )
    
    return 0


if __name__ == "__main__":
    exit(main())
