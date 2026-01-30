"""
Violin plot visualization for protein and peptide distribution per genome across 4 marker categories.

Input: TSV with columns: GenomeId, Category, ProteinCount, PeptideCount

Output: Two violin plots (HTML) showing:
1. Distribution of protein counts per genome for each category
2. Distribution of peptide counts per genome for each category

Usage:
  python protein_peptide_dist_viz.py input.tsv output.html
"""

import sys
import csv
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def read_table(path: Path):
    """Read TSV file and return as DataFrame."""
    with path.open("r", newline="") as f:
        # sniff delimiter (tab vs comma)
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.reader(f, dialect)
        header = next(reader)
        rows = list(reader)
    
    df = pd.DataFrame(rows, columns=header)
    # Convert numeric columns
    if 'ProteinCount' in df.columns:
        df['ProteinCount'] = pd.to_numeric(df['ProteinCount'], errors='coerce')
    if 'PeptideCount' in df.columns:
        df['PeptideCount'] = pd.to_numeric(df['PeptideCount'], errors='coerce')
    
    return df


def create_violin_plots(df, out_path: Path):
    """
    Create two violin plots: one for protein counts, one for peptide counts.
    X-axis: 4 marker categories, Y-axis: count distribution per genome.
    """
    categories = ['Anchor', 'HighRes', 'Abundance', 'Stability']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    # Create subplots: 2 rows, 1 column
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('Protein Count Distribution per Genome', 
                       'Peptide Count Distribution per Genome'),
        vertical_spacing=0.15
    )
    
    # Plot 1: Protein counts
    for idx, cat in enumerate(categories):
        cat_data = df[df['Category'] == cat]['ProteinCount'].dropna()
        if len(cat_data) > 0:
            fig.add_trace(
                go.Violin(
                    y=cat_data,
                    name=cat,
                    box_visible=True,
                    meanline_visible=True,
                    fillcolor=colors[idx],
                    line_color=colors[idx],
                    opacity=0.7
                ),
                row=1, col=1
            )
    
    # Plot 2: Peptide counts
    for idx, cat in enumerate(categories):
        cat_data = df[df['Category'] == cat]['PeptideCount'].dropna()
        if len(cat_data) > 0:
            fig.add_trace(
                go.Violin(
                    y=cat_data,
                    name=cat,
                    box_visible=True,
                    meanline_visible=True,
                    fillcolor=colors[idx],
                    line_color=colors[idx],
                    opacity=0.7,
                    showlegend=False  # Only show legend once
                ),
                row=2, col=1
            )
    
    # Update layout
    fig.update_layout(
        title_text='Distribution of Protein and Peptide Counts per Genome Across Marker Categories',
        title_x=0.5,
        height=800,
        xaxis_title='Category',
        yaxis_title='Count per Genome',
        yaxis2_title='Count per Genome',
        hovermode='closest'
    )
    
    # Update x-axis labels for both subplots
    fig.update_xaxes(title_text='Category', row=1, col=1)
    fig.update_xaxes(title_text='Category', row=2, col=1)
    
    # Save as HTML
    if out_path.suffix.lower() == '.html':
        fig.write_html(str(out_path))
    else:
        # Try to save as static image if kaleido is available
        try:
            fig.write_image(str(out_path))
        except Exception:
            # Fallback to HTML
            html_path = out_path.with_suffix('.html')
            fig.write_html(str(html_path))
            print(f"Kaleido not available, saved as HTML instead: {html_path}")
    
    print(f"Saved violin plots to {out_path}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python protein_peptide_dist_viz.py <input_tsv> <output_html>")
        sys.exit(1)
    
    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    
    if not in_path.exists():
        print(f"Error: Input file not found: {in_path}")
        sys.exit(1)
    
    df = read_table(in_path)
    print(f"Loaded {len(df)} rows from {in_path}")
    print(f"Categories: {df['Category'].unique()}")
    print(f"Summary statistics:")
    print(df.groupby('Category')[['ProteinCount', 'PeptideCount']].describe())
    
    create_violin_plots(df, out_path)


if __name__ == "__main__":
    main()

