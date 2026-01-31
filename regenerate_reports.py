#!/usr/bin/env python3
"""
Standalone script to regenerate MetaPilot reports from existing data.
Use this to test report fixes without re-running the full workflow.

Prerequisites: pip install pandas plotly numpy scikit-learn

Usage:
  python regenerate_reports.py <result_directory>
  
Example:
  python regenerate_reports.py "z:\\PRIDE_dataset\\PXD069517\\MetaPilot\\pfind_MAG_search"
  
This will regenerate report_function_summary.html and report_taxonomy_summary.html
in the report/reports/ subdirectory, using metainfo.tsv for experiment names and
metadata classification.

Alternatively: Rebuild MetaPilot JAR (pack_jar.bat or mvn clean package) and re-run
the workflow - the packaged Python scripts will then include these fixes.
"""

import sys
import os

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nError: Please provide the result directory path (e.g. pfind_MAG_search folder)")
        sys.exit(1)
    
    result_dir = os.path.abspath(sys.argv[1])
    if not os.path.isdir(result_dir):
        print(f"Error: Directory not found: {result_dir}")
        sys.exit(1)
    
    # Use canonical report scripts from metapilot-windows-gui (single source for GUI and future HPC)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    python_scripts = os.path.join(script_dir, "metapilot-windows-gui", "src", "main", "resources", "python_scripts")
    if not os.path.exists(python_scripts):
        print(f"Error: Python scripts not found. Expected at: {python_scripts}")
        print("Run this script from the MetaPilot project root (parent of metapilot-windows-gui).")
        sys.exit(1)
    sys.path.insert(0, python_scripts)
    
    from report_generation.metapilot_report_generator import MetaPilotReportGenerator
    import pandas as pd
    
    # Paths
    metainfo = os.path.join(result_dir, "metainfo.tsv")
    taxa_file = os.path.join(result_dir, "taxonomy_analysis", "Taxa_report.tsv")
    func_dir = os.path.join(result_dir, "functional_annotation")
    report_dir = os.path.join(result_dir, "report", "reports")
    
    os.makedirs(report_dir, exist_ok=True)
    
    # Load metadata
    metadata = None
    if os.path.exists(metainfo):
        try:
            rg = MetaPilotReportGenerator(report_dir)
            metadata = rg.load_metadata(metainfo)
            print(f"Loaded metadata from {metainfo} ({len(metadata)} rows)")
        except Exception as e:
            print(f"Warning: Could not load metadata: {e}")
    
    if metadata is None:
        print("Warning: No metadata - reports will use raw file names")
    
    rg = MetaPilotReportGenerator(report_dir)
    
    # Regenerate taxonomy report
    if os.path.exists(taxa_file):
        try:
            taxa_df = pd.read_csv(taxa_file, sep='\t')
            fig = rg.create_taxonomy_summary(taxa_df, metadata)
            out_path = os.path.join(report_dir, "report_taxonomy_summary.html")
            rg.generate_html_report(fig, "Taxonomy Summary Report", out_path)
            print(f"Regenerated: {out_path}")
        except Exception as e:
            print(f"Error regenerating taxonomy report: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Skipped taxonomy: {taxa_file} not found")
    
    # Regenerate function report
    if os.path.isdir(func_dir):
        try:
            fig = rg.create_function_summary(metadata, input_file_path=func_dir)
            out_path = os.path.join(report_dir, "report_function_summary.html")
            rg.generate_html_report(fig, "Function Summary Report", out_path)
            print(f"Regenerated: {out_path}")
        except Exception as e:
            print(f"Error regenerating function report: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Skipped function: {func_dir} not found")
    
    print("\nDone. Open the HTML files in your browser to verify.")

if __name__ == "__main__":
    main()
