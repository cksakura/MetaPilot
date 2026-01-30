"""
Violin plot visualization for distinct peptide distribution per genome across 4 marker categories.

Input: TSV with columns: GenomeId, Category, DistinctPeptideCount
(peptides that appear ONLY in that genome, not shared with others)

Output: Violin plot (HTML) showing distribution of distinct peptide counts per genome
for each marker category.

Usage:
  python distinct_peptide_dist_viz.py input.tsv output.html
"""

import sys
import csv
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go


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
    if 'DistinctPeptideCount' in df.columns:
        df['DistinctPeptideCount'] = pd.to_numeric(df['DistinctPeptideCount'], errors='coerce')
    
    return df


def create_violin_plot(df, out_path: Path):
    """
    Create violin plot showing distribution of distinct peptide counts per genome
    for each marker category.
    X-axis: 4 marker categories, Y-axis: distinct peptide count distribution per genome.
    """
    categories = ['Anchor', 'HighRes', 'Abundance', 'Stability']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    fig = go.Figure()
    
    for idx, cat in enumerate(categories):
        cat_data = df[df['Category'] == cat]['DistinctPeptideCount'].dropna()
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
                )
            )
    
    fig.update_layout(
        title_text='Distribution of Distinct Peptide Counts per Genome Across Marker Categories<br><sub>Peptides unique to each genome (not shared with others)</sub>',
        title_x=0.5,
        height=600,
        xaxis_title='Category',
        yaxis_title='Distinct Peptide Count per Genome',
        hovermode='closest',
        showlegend=True
    )
    
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
    
    print(f"Saved violin plot to {out_path}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python distinct_peptide_dist_viz.py <input_tsv> <output_html>")
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
    print(df.groupby('Category')['DistinctPeptideCount'].describe())
    
    create_violin_plot(df, out_path)


if __name__ == "__main__":
    main()

