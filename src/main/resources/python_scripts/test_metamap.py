#!/usr/bin/env python3
"""
Test script for MetaMap generator using the tree.js data format.
"""

import sys
import os
import json
from pathlib import Path

# Use canonical report scripts from metapilot-windows-gui (single source)
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.normpath(os.path.join(_script_dir, "..", "..", "..", ".."))
_canonical_scripts = os.path.join(_project_root, "metapilot-windows-gui", "src", "main", "resources", "python_scripts")
if os.path.exists(_canonical_scripts):
    sys.path.insert(0, _canonical_scripts)
else:
    sys.path.insert(0, _script_dir)  # fallback when run from extracted JAR (report_generation is alongside)

from report_generation.metapilot_report_generator import MetaPilotReportGenerator

def test_metamap_generation():
    """Test MetaMap generation with sample tree data."""
    
    # Path to the tree.js file
    tree_file = Path(__file__).parent.parent.parent / "tests" / "tree.js"
    
    if not tree_file.exists():
        print(f"Error: Tree file not found at {tree_file}")
        return False
    
    # Read the tree data
    with open(tree_file, 'r') as f:
        content = f.read()
    
    # Extract the JSON data from the JavaScript file
    # The data is in the format: let data = '[{...}]';
    start_idx = content.find("'[")
    end_idx = content.rfind("]'") + 1
    
    if start_idx == -1 or end_idx == -1:
        print("Error: Could not find JSON data in tree.js file")
        return False
    
    json_str = content[start_idx + 1:end_idx]
    
    try:
        tree_data = json.loads(json_str)
        print(f"Loaded {len(tree_data)} tree entries")
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON: {e}")
        return False
    
    # Create output directory
    output_dir = Path(__file__).parent / "test_output"
    output_dir.mkdir(exist_ok=True)
    
    # Initialize the report generator
    generator = MetaPilotReportGenerator(str(output_dir))
    
    # Generate bubble chart report
    print("Generating bubble chart report...")
    fig = generator.create_bubble_chart_summary(
        tree_data=tree_data,
        metadata=None,
        max_depth=4,
        min_abundance=0.001
    )
    
    if fig:
        # Generate HTML report
        output_file = output_dir / "metamap_report.html"
        success = generator.generate_html_report(
            fig=fig,
        title="Bubble Chart: Hierarchical Taxonomic Abundance",
        output_file=str(output_file),
        summary="Interactive bubble chart showing hierarchical taxonomic abundance with treemap, sunburst, and distribution visualizations.",
        tags=["bubble-chart", "taxonomy", "hierarchical", "treemap"]
        )
        
        if success:
            print(f"✅ Bubble chart report generated successfully: {output_file}")
            return True
        else:
            print("❌ Failed to generate HTML report")
            return False
    else:
        print("❌ Failed to create bubble chart visualization")
        return False

if __name__ == "__main__":
    print("Testing Bubble Chart Generator...")
    success = test_metamap_generation()
    if success:
        print("🎉 Test completed successfully!")
    else:
        print("💥 Test failed!")
        sys.exit(1)
