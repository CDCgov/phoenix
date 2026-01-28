#!/usr/bin/env python3

import argparse
import os
import glob
from pathlib import Path

def rename_griphin_files(directory):
    """
    Rename files ending in *_GRiPHin_Summary.tsv or *_GRiPHin_Summary.xlsx
    to {prefix}_GRiPHin.tsv or {prefix}_GRiPHin.xlsx
    """
    # Search for both .tsv and .xlsx files
    patterns = ['*_GRiPHin_Summary.tsv', '*_GRiPHin_Summary.xlsx']
    
    renamed_count = 0
    
    for pattern in patterns:
        search_path = os.path.join(directory, pattern)
        files = glob.glob(search_path)
        
        for file_path in files:
            # Get the file components
            dir_name = os.path.dirname(file_path)
            file_name = os.path.basename(file_path)
            
            # Determine the extension
            if file_name.endswith('_GRiPHin_Summary.tsv'):
                prefix = file_name.replace('_GRiPHin_Summary.tsv', '')
                new_name = f"{prefix}_GRiPHin.tsv"
            elif file_name.endswith('_GRiPHin_Summary.xlsx'):
                prefix = file_name.replace('_GRiPHin_Summary.xlsx', '')
                new_name = f"{prefix}_GRiPHin.xlsx"
            else:
                continue
            
            # Create the new file path
            new_path = os.path.join(dir_name, new_name)
            
            # Rename the file
            try:
                os.rename(file_path, new_path)
                print(f"Renamed: {file_name} -> {new_name}")
                renamed_count += 1
            except Exception as e:
                print(f"Error renaming {file_name}: {e}")
    
    print(f"\nTotal files renamed: {renamed_count}")

def main():
    parser = argparse.ArgumentParser(description='Rename GRiPHin_Summary files to shorter format')
    parser.add_argument( '-d', '--directory', type=str, default='.', help='Directory containing files to rename (default: current directory)' )
    parser.add_argument('--dry-run', action='store_true', help='Show what would be renamed without actually renaming' )
    args = parser.parse_args()
    
    # Check if directory exists
    if not os.path.isdir(args.directory):
        print(f"Error: Directory '{args.directory}' does not exist")
        return
    
    if args.dry_run:
        print("DRY RUN MODE - No files will be renamed\n")
        patterns = ['*_GRiPHin_Summary.tsv', '*_GRiPHin_Summary.xlsx']
        for pattern in patterns:
            search_path = os.path.join(args.directory, pattern)
            files = glob.glob(search_path)
            for file_path in files:
                file_name = os.path.basename(file_path)
                if file_name.endswith('_GRiPHin_Summary.tsv'):
                    prefix = file_name.replace('_GRiPHin_Summary.tsv', '')
                    new_name = f"{prefix}_GRiPHin.tsv"
                elif file_name.endswith('_GRiPHin_Summary.xlsx'):
                    prefix = file_name.replace('_GRiPHin_Summary.xlsx', '')
                    new_name = f"{prefix}_GRiPHin.xlsx"
                print(f"Would rename: {file_name} -> {new_name}")
    else:
        rename_griphin_files(args.directory)

if __name__ == '__main__':
    main()