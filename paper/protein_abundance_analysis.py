"""
Statistical analysis of protein abundances for 4 marker categories.

Input: TSV with columns: Sample, GenomeId, Category, Abundance_All, Abundance_Distinct, Rank_All, Rank_Distinct

Output: 
- Statistical summary (console + text file)
- Visualizations (HTML with multiple plots)
- Category comparison analysis

Usage:
  python protein_abundance_analysis.py input.tsv output_prefix
"""

import sys
import csv
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


def read_data(path: Path):
    """Read TSV file and return as DataFrame."""
    df = pd.read_csv(path, sep='\t')
    
    # Convert numeric columns
    numeric_cols = ['Abundance', 'Ratio', 'Total_Genome_Abundance', 
                   'Abundance_All', 'Abundance_Distinct', 'Ratio_All', 'Ratio_Distinct', 
                   'Marker4_Total_All', 'Marker4_Total_Distinct', 'Rank_All', 'Rank_Distinct']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    return df


def calculate_statistics(df, output_prefix: Path):
    """Calculate comprehensive statistics and save to text file."""
    stats_file = output_prefix.with_suffix('.stats.txt')
    
    with open(stats_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PROTEIN ABUNDANCE STATISTICAL ANALYSIS\n")
        f.write("=" * 80 + "\n\n")
        
        # Basic info
        f.write(f"Total samples: {df['Sample'].nunique()}\n")
        f.write(f"Total genomes: {df['GenomeId'].nunique()}\n")
        f.write(f"Total observations: {len(df)}\n\n")
        
        categories = ['Anchor', 'HighRes', 'Abundance', 'Stability']
        
        # Summary statistics per category
        abundance_col = 'Abundance' if 'Abundance' in df.columns else 'Abundance_All'
        f.write("=" * 80 + "\n")
        f.write(f"SUMMARY STATISTICS BY CATEGORY ({abundance_col})\n")
        f.write("=" * 80 + "\n\n")
        
        summary_stats = df.groupby('Category')[abundance_col].agg([
            'count', 'mean', 'median', 'std', 'min', 'max',
            lambda x: np.percentile(x, 25),  # Q1
            lambda x: np.percentile(x, 75)  # Q3
        ])
        summary_stats.columns = ['Count', 'Mean', 'Median', 'Std', 'Min', 'Max', 'Q1', 'Q3']
        f.write(summary_stats.to_string())
        f.write("\n\n")
        
        # Ratio analysis
        ratio_col = 'Ratio' if 'Ratio' in df.columns else 'Ratio_All'
        f.write("=" * 80 + "\n")
        f.write(f"ABUNDANCE RATIO ANALYSIS ({ratio_col} - within each genome)\n")
        f.write("=" * 80 + "\n\n")
        f.write("Ratio = Category Abundance / Total Genome Abundance (all proteins)\n")
        f.write("This shows what percentage of total protein abundance comes from each marker category.\n\n")
        
        if ratio_col in df.columns:
            ratio_stats = df[df[ratio_col] > 0].groupby('Category')[ratio_col].agg([
                'count', 'mean', 'median', 'std', 'min', 'max'
            ])
            ratio_stats.columns = ['Count', 'Mean', 'Median', 'Std', 'Min', 'Max']
            f.write("Ratio statistics (category abundance / total genome abundance per genome):\n")
            f.write(ratio_stats.to_string())
            f.write("\n\n")
            f.write("Note: Ratios are shown as decimals. Multiply by 100 to get percentages.\n")
            f.write("Example: Ratio = 0.15 means 15%% of total genome abundance comes from this category.\n\n")
            
            # Dominant category (highest ratio) per genome
            f.write("=" * 80 + "\n")
            f.write(f"DOMINANT CATEGORY PER GENOME (Highest {ratio_col})\n")
            f.write("=" * 80 + "\n\n")
            
            # For each sample-genome combination, find category with highest ratio
            dominant_categories = []
            for (sample, genome), group in df.groupby(['Sample', 'GenomeId']):
                if ratio_col in group.columns:
                    max_idx = group[ratio_col].idxmax()
                    dominant_cat = group.loc[max_idx, 'Category']
                    dominant_categories.append(dominant_cat)
            
            if dominant_categories:
                from collections import Counter
                dominant_counts = Counter(dominant_categories)
                f.write("Frequency of dominant category:\n")
                for cat, count in dominant_counts.most_common():
                    f.write(f"  {cat}: {count} ({100*count/len(dominant_categories):.2f}%)\n")
                f.write("\n")
        
        # Category dominance (which category has highest total abundance)
        f.write("=" * 80 + "\n")
        f.write("TOTAL ABUNDANCE BY CATEGORY (across all samples and genomes)\n")
        f.write("=" * 80 + "\n\n")
        
        total_by_category = df.groupby('Category')[abundance_col].sum().sort_values(ascending=False)
        f.write(total_by_category.to_string())
        f.write("\n\n")
        f.write("Percentage of total:\n")
        total_sum = total_by_category.sum()
        for cat, val in total_by_category.items():
            f.write(f"  {cat}: {100*val/total_sum:.2f}%\n")
        f.write("\n")
        
        # Sample-level analysis
        f.write("=" * 80 + "\n")
        f.write("SAMPLE-LEVEL SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        sample_summary = df.groupby('Sample')['Abundance_All'].agg(['sum', 'mean', 'count'])
        sample_summary.columns = ['Total_Abundance', 'Mean_Abundance', 'Genome_Count']
        f.write(f"Number of samples: {len(sample_summary)}\n")
        f.write(f"Mean total abundance per sample: {sample_summary['Total_Abundance'].mean():.2f}\n")
        f.write(f"Mean genomes detected per sample: {sample_summary['Genome_Count'].mean():.2f}\n\n")
    
    print(f"Statistics saved to: {stats_file}")
    return stats_file


def create_visualizations(df, output_prefix: Path):
    """Create comprehensive visualizations."""
    categories = ['Anchor', 'HighRes', 'Abundance', 'Stability']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    # Create subplots
    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=(
            'Abundance Distribution by Category (All Peptides)',
            'Abundance Distribution by Category (Log Scale)',
            'Abundance Ratio Distribution by Category',
            'Total Abundance by Category',
            'Category Dominance (Highest Ratio Frequency)',
            'Sample-Level Total Abundance Distribution'
        ),
        specs=[[{"type": "violin"}, {"type": "violin"}],
               [{"type": "violin"}, {"type": "bar"}],
               [{"type": "bar"}, {"type": "histogram"}]],
        vertical_spacing=0.12,
        horizontal_spacing=0.1
    )
    
    # Determine column names (support both old and new format)
    abundance_col = 'Abundance' if 'Abundance' in df.columns else 'Abundance_All'
    ratio_col = 'Ratio' if 'Ratio' in df.columns else 'Ratio_All'
    
    # 1. Violin plot: Abundance distribution (linear scale)
    for idx, cat in enumerate(categories):
        cat_data = df[df['Category'] == cat][abundance_col].dropna()
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
    
    # 2. Violin plot: Abundance distribution (log scale)
    df_log = df.copy()
    df_log['Abundance_Log'] = np.log10(df_log[abundance_col] + 1)
    for idx, cat in enumerate(categories):
        cat_data = df_log[df_log['Category'] == cat]['Abundance_Log'].dropna()
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
                    showlegend=False
                ),
                row=1, col=2
            )
    
    # 3. Ratio distribution
    if ratio_col in df.columns:
        for idx, cat in enumerate(categories):
            cat_data = df[(df['Category'] == cat) & (df[ratio_col] > 0)][ratio_col].dropna()
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
                        showlegend=False
                    ),
                    row=2, col=1
                )
    
    # 4. Total abundance by category
    total_by_category = df.groupby('Category')[abundance_col].sum().sort_values(ascending=False)
    fig.add_trace(
        go.Bar(
            x=total_by_category.index,
            y=total_by_category.values,
            marker_color=[colors[categories.index(c)] for c in total_by_category.index],
            opacity=0.7,
            showlegend=False
        ),
        row=2, col=2
    )
    
    # 5. Category dominance (highest ratio frequency)
    if ratio_col in df.columns:
        dominant_categories = []
        for (sample, genome), group in df.groupby(['Sample', 'GenomeId']):
            if len(group) > 0 and ratio_col in group.columns:
                max_idx = group[ratio_col].idxmax()
                dominant_cat = group.loc[max_idx, 'Category']
                dominant_categories.append(dominant_cat)
        
        if dominant_categories:
            from collections import Counter
            dominant_counts = Counter(dominant_categories)
            dominant_df = pd.DataFrame(list(dominant_counts.items()), columns=['Category', 'Count'])
            dominant_df = dominant_df.sort_values('Count', ascending=False)
            
            fig.add_trace(
                go.Bar(
                    x=dominant_df['Category'],
                    y=dominant_df['Count'],
                    marker_color=[colors[categories.index(c)] for c in dominant_df['Category']],
                    opacity=0.7,
                    showlegend=False
                ),
                row=3, col=1
            )
    
    # 6. Sample-level total abundance
    sample_totals = df.groupby('Sample')[abundance_col].sum()
    fig.add_trace(
        go.Histogram(
            x=sample_totals.values,
            nbinsx=30,
            marker_color='#9467bd',
            opacity=0.7,
            showlegend=False
        ),
        row=3, col=2
    )
    
    # Update layout
    fig.update_layout(
        title_text='Protein Abundance Statistical Analysis Dashboard',
        title_x=0.5,
        height=1400,
        showlegend=True
    )
    
    # Update axes labels
    fig.update_yaxes(title_text="Abundance", row=1, col=1)
    fig.update_yaxes(title_text="Log10(Abundance + 1)", row=1, col=2)
    fig.update_yaxes(title_text="Ratio (Category/Total Marker-4)", row=2, col=1)
    fig.update_yaxes(title_text="Total Abundance", row=2, col=2)
    fig.update_yaxes(title_text="Frequency (Dominant Category)", row=3, col=1)
    fig.update_yaxes(title_text="Frequency", row=3, col=2)
    
    fig.update_xaxes(title_text="Category", row=1, col=1)
    fig.update_xaxes(title_text="Category", row=1, col=2)
    fig.update_xaxes(title_text="Category", row=2, col=1)
    fig.update_xaxes(title_text="Category", row=2, col=2)
    fig.update_xaxes(title_text="Category", row=3, col=1)
    fig.update_xaxes(title_text="Total Abundance per Sample", row=3, col=2)
    
    # Save
    html_file = output_prefix.with_suffix('.html')
    fig.write_html(str(html_file))
    print(f"Visualizations saved to: {html_file}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python protein_abundance_analysis.py <input_tsv> [output_prefix]")
        sys.exit(1)
    
    in_path = Path(sys.argv[1])
    if len(sys.argv) >= 3:
        out_prefix = Path(sys.argv[2])
    else:
        out_prefix = in_path.with_suffix('')  # Remove extension
    
    if not in_path.exists():
        print(f"Error: Input file not found: {in_path}")
        sys.exit(1)
    
    print(f"Loading data from {in_path}...")
    df = read_data(in_path)
    print(f"Loaded {len(df)} rows")
    print(f"Categories: {df['Category'].unique()}")
    print(f"Samples: {df['Sample'].nunique()}")
    print(f"Genomes: {df['GenomeId'].nunique()}\n")
    
    print("Calculating statistics...")
    calculate_statistics(df, out_prefix)
    
    print("Creating visualizations...")
    create_visualizations(df, out_prefix)
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()

