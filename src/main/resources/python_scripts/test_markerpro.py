import csv
import re
import sys
import os
import pandas as pd
import plotly.express as px

# Check if file path was provided as argument
if len(sys.argv) < 2:
    print("Usage: python test_markerpro.py <path_to_functions.tsv>")
    sys.exit(1)

tsv_file_path = sys.argv[1]

# Check if file exists
if not os.path.exists(tsv_file_path):
    print(f"Error: File '{tsv_file_path}' not found.")
    sys.exit(1)

# 1. Define Regex Patterns for your 4 Categories
patterns = {
    "Anchor": r"(?i).*ribosomal\s+protein.*",
    "High Resolution": r"(?i).*(polymerase|gyrase).*",
    "Abundance": r"(?i).*elongation\s+factor.*",
    "Stability": r"(?i).*tRNA\s+(ligase|synthetase).*"
}

data = []

# 2. Parse the File
with open(tsv_file_path, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader) 
    
    # Map columns by name
    try:
        idx_name = header.index("Name")
        idx_prot = header.index("Protein name")
        # Grab intensity columns dynamically
        int_indices = [i for i, col in enumerate(header) if "Intensity" in col]
    except ValueError:
        print("Error: Could not find required columns in header.")
        exit()

    for row in reader:
        # Safety check for short rows
        if len(row) <= max(idx_prot, max(int_indices, default=0)): 
            continue
            
        name = row[idx_name]
        prot_name = row[idx_prot]
        
        # Calculate Mean Intensity
        try:
            intensities = [float(row[i]) for i in int_indices]
            if not intensities: continue
            mean_intensity = sum(intensities) / len(intensities)
        except ValueError:
            continue
            
        # 3. Categorize Protein
        category = None
        for cat, pat in patterns.items():
            if re.match(pat, prot_name):
                category = cat
                break 
        
        if category:
            # Extract Genome ID (Part before the first underscore)
            genome = name.split('_')[0]
            data.append({
                "Genome": genome,
                "Protein": prot_name,
                "Category": category,
                "Intensity": mean_intensity
            })

# 4. Create DataFrame and Plot
df_res = pd.DataFrame(data)

if not df_res.empty:
    # Determine output directory (same as input TSV file)
    output_dir = os.path.dirname(os.path.abspath(tsv_file_path))
    html_output = os.path.join(output_dir, 'intensity_distribution.html')
    png_output = os.path.join(output_dir, 'intensity_distribution.png')
    csv_output = os.path.join(output_dir, 'marker_analysis_results.csv')
    
    # Save the data table
    df_res.to_csv(csv_output, index=False)
    print(f"Data table saved as: {csv_output}")
    
    # Create Interactive Boxplot with Plotly
    fig = px.box(
        df_res, 
        x='Genome', 
        y='Intensity', 
        color='Category',
        log_y=True,  # Log scale for intensity
        title='Intensity Distribution of Functional Markers per Genome',
        labels={
            'Genome': 'Genome ID',
            'Intensity': 'Intensity (Log Scale)',
            'Category': 'Marker Category'
        },
        color_discrete_sequence=px.colors.qualitative.Set2
    )
    
    # Update layout for better readability
    fig.update_layout(
        xaxis_tickangle=-45,
        width=1200,
        height=800,
        hovermode='closest',
        legend=dict(
            title='Marker Category',
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        )
    )
    
    # Save as HTML (interactive) and PNG (static)
    fig.write_html(html_output)
    fig.write_image(png_output)
    
    # Show interactive plot in browser
    fig.show()
    
    print(f"\nAnalysis Complete. Found {len(df_res)} matching proteins across {df_res['Genome'].nunique()} genomes.")
    print(f"Interactive plot saved as: {html_output}")
    print(f"Static plot saved as: {png_output}")
else:
    print("No matching proteins found.")