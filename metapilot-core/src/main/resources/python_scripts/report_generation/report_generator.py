import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from pathlib import Path
import json
from datetime import datetime

class ReportGenerator:
    def __init__(self, output_dir="reports"):
        """Initialize the report generator with output directory."""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Define color schemes
        self.colors = {
            'primary': '#1f77b4',
            'secondary': '#ff7f0e',
            'tertiary': '#2ca02c',
            'background': '#f8f9fa',
            'text': '#2c3e50'
        }
        
        # Define plot templates
        self.template = {
            'layout': {
                'paper_bgcolor': self.colors['background'],
                'plot_bgcolor': self.colors['background'],
                'font': {'color': self.colors['text']},
                'margin': {'t': 30, 'l': 30, 'r': 30, 'b': 30}
            }
        }

    def create_peptides_summary(self, peptides_df):
        """Generate peptides summary report with visualizations."""
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Peptide Length Distribution',
                'Peptide Charge Distribution',
                'Peptide Mass Distribution',
                'Peptide Score Distribution'
            )
        )

        # Length distribution
        fig.add_trace(
            go.Histogram(
                x=peptides_df['Length'],
                name='Length',
                marker_color=self.colors['primary']
            ),
            row=1, col=1
        )

        # Charge distribution
        fig.add_trace(
            go.Histogram(
                x=peptides_df['Charge'],
                name='Charge',
                marker_color=self.colors['secondary']
            ),
            row=1, col=2
        )

        # Mass distribution
        fig.add_trace(
            go.Histogram(
                x=peptides_df['Mass'],
                name='Mass',
                marker_color=self.colors['tertiary']
            ),
            row=2, col=1
        )

        # Score distribution
        fig.add_trace(
            go.Histogram(
                x=peptides_df['Score'],
                name='Score',
                marker_color=self.colors['primary']
            ),
            row=2, col=2
        )

        # Update layout
        fig.update_layout(
            height=800,
            showlegend=False,
            template=self.template
        )

        return fig

    def create_proteins_summary(self, proteins_df):
        """Generate proteins summary report with visualizations."""
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Protein Length Distribution',
                'Protein Coverage Distribution',
                'Number of Peptides per Protein',
                'Protein Score Distribution'
            )
        )

        # Length distribution
        fig.add_trace(
            go.Histogram(
                x=proteins_df['Length'],
                name='Length',
                marker_color=self.colors['primary']
            ),
            row=1, col=1
        )

        # Coverage distribution
        fig.add_trace(
            go.Histogram(
                x=proteins_df['Coverage'],
                name='Coverage',
                marker_color=self.colors['secondary']
            ),
            row=1, col=2
        )

        # Peptides per protein
        fig.add_trace(
            go.Histogram(
                x=proteins_df['Peptides'],
                name='Peptides',
                marker_color=self.colors['tertiary']
            ),
            row=2, col=1
        )

        # Score distribution
        fig.add_trace(
            go.Histogram(
                x=proteins_df['Score'],
                name='Score',
                marker_color=self.colors['primary']
            ),
            row=2, col=2
        )

        # Update layout
        fig.update_layout(
            height=800,
            showlegend=False,
            template=self.template
        )

        return fig

    def create_id_summary(self, summary_df):
        """Generate ID summary report with visualizations."""
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'PSM Distribution',
                'Peptide Distribution',
                'Protein Distribution',
                'Score Distribution'
            )
        )

        # PSM distribution
        fig.add_trace(
            go.Bar(
                x=summary_df['Sample'],
                y=summary_df['PSMs'],
                name='PSMs',
                marker_color=self.colors['primary']
            ),
            row=1, col=1
        )

        # Peptide distribution
        fig.add_trace(
            go.Bar(
                x=summary_df['Sample'],
                y=summary_df['Peptides'],
                name='Peptides',
                marker_color=self.colors['secondary']
            ),
            row=1, col=2
        )

        # Protein distribution
        fig.add_trace(
            go.Bar(
                x=summary_df['Sample'],
                y=summary_df['Proteins'],
                name='Proteins',
                marker_color=self.colors['tertiary']
            ),
            row=2, col=1
        )

        # Score distribution
        fig.add_trace(
            go.Box(
                y=summary_df['Score'],
                name='Score',
                marker_color=self.colors['primary']
            ),
            row=2, col=2
        )

        # Update layout
        fig.update_layout(
            height=800,
            showlegend=False,
            template=self.template
        )

        return fig

    def create_taxonomy_summary(self, taxa_df):
        """Create taxonomy summary visualization with individual plots for each taxonomic level"""
        try:
            # Extract sample names from metadata
            intensity_cols = [col for col in taxa_df.columns if col.startswith('Intensity_')]
            
            if len(intensity_cols) == 0:
                return self.create_empty_data_plot("Taxonomy Summary Report", "No intensity data available")

            # Use simple names from column headers
            simple_names = []
            for col in intensity_cols:
                # Extract sample name from column header like "Intensity_Haonan_20210924_isotope_mix_Stool1"
                sample_name = col.replace('Intensity_', '')
                simple_names.append(sample_name)

            # Ensure intensity columns are numeric
            for c in intensity_cols:
                if c in taxa_df.columns:
                    taxa_df[c] = pd.to_numeric(taxa_df[c], errors="coerce").fillna(0)

            # Define taxonomic ranks to process
            taxonomic_ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
            
            # Create a single figure with all taxonomic levels
            # Calculate total rows: number of available ranks × 3 plots per rank
            available_ranks = [rank for rank in taxonomic_ranks if rank in taxa_df.columns and len(taxa_df[taxa_df[rank].notna() & (taxa_df[rank] != "")]) > 0]
            total_rows = len(available_ranks) * 3
            
            if total_rows == 0:
                return self.create_empty_data_plot("Taxonomy Summary Report", "No taxonomic data available")
            
            # Create subplot titles and specs
            subplot_titles = []
            specs = []
            
            for rank in available_ranks:
                subplot_titles.extend([
                    f"{rank.capitalize()} Intensity Distribution (Log Intensity)",
                    f"{rank.capitalize()} Log Intensity Heatmap (Clustered)",
                    f"{rank.capitalize()} Composition by Sample"
                ])
                specs.extend([
                    [{"type": "violin"}],
                    [{"type": "heatmap"}],
                    [{"type": "bar"}]
                ])
            
            # Create the main figure
            fig = make_subplots(
                rows=total_rows,
                cols=1,
                subplot_titles=subplot_titles,
                specs=specs,
                vertical_spacing=0.05
            )
            
            current_row = 1
            
            # Process each taxonomic rank
            for rank in available_ranks:
                # Get data for this rank
                rank_data = taxa_df[taxa_df[rank].notna() & (taxa_df[rank] != "")]
                
                # 1. Violin Plot - Intensity Distribution with LOG intensity
                intensity_values = []
                sample_labels = []
                for _, row in rank_data.iterrows():
                    for col in intensity_cols:
                        if pd.notna(row[col]) and row[col] > 0:
                            # Convert to log intensity for better visualization
                            log_intensity = np.log10(row[col])
                            intensity_values.append(log_intensity)
                            # Get sample name for x-axis
                            sample_idx = intensity_cols.index(col)
                            sample_name = simple_names[sample_idx] if sample_idx < len(simple_names) else f"Sample_{sample_idx+1}"
                            sample_labels.append(sample_name)

                if intensity_values:
                    fig.add_trace(
                        go.Violin(
                            x=sample_labels,
                            y=intensity_values,
                            name=f"{rank.capitalize()} Intensity Distribution",
                            box_visible=True,
                            meanline_visible=True,
                            points='all',
                            pointpos=0,
                            jitter=0.3,
                            marker_color=self.colors["primary"]
                        ),
                        row=current_row,
                        col=1
                    )

                current_row += 1

                # 2. Heatmap - All taxa vs experiments
                # Prepare heatmap data
                heatmap_data = []
                taxon_labels = []

                # Get all unique taxa for this rank
                unique_taxa = rank_data[rank].unique()
                unique_taxa = [t for t in unique_taxa if pd.notna(t) and t != ""]

                for taxon in unique_taxa:
                    taxon_data = rank_data[rank_data[rank] == taxon]
                    sample_intensities = []
                    for col in intensity_cols:
                        intensity_val = taxon_data[col].sum() if len(taxon_data) > 0 else 0
                        sample_intensities.append(intensity_val)
                    heatmap_data.append(sample_intensities)
                    taxon_labels.append(taxon)

                if heatmap_data:
                    # Apply hierarchical clustering to heatmap data
                    from scipy.cluster.hierarchy import linkage, dendrogram
                    from scipy.spatial.distance import pdist
                    
                    # Convert to numpy array for clustering
                    heatmap_array = np.array(heatmap_data)
                    
                    # Cluster rows (taxa) and columns (samples) separately
                    if heatmap_array.shape[0] > 1:  # Need at least 2 taxa for clustering
                        # Cluster taxa (rows)
                        row_linkage = linkage(pdist(heatmap_array), method='ward')
                        row_order = dendrogram(row_linkage, no_plot=True)['leaves']
                        heatmap_array = heatmap_array[row_order]
                        taxon_labels_clustered = [taxon_labels[i] for i in row_order]
                    else:
                        taxon_labels_clustered = taxon_labels
                    
                    if heatmap_array.shape[1] > 1:  # Need at least 2 samples for clustering
                        # Cluster samples (columns)
                        col_linkage = linkage(pdist(heatmap_array.T), method='ward')
                        col_order = dendrogram(col_linkage, no_plot=True)['leaves']
                        heatmap_array = heatmap_array[:, col_order]
                        sample_names_clustered = [simple_names[i] if i < len(simple_names) else f"Sample_{i+1}" for i in col_order]
                    else:
                        sample_names_clustered = simple_names[:len(intensity_cols)] if len(simple_names) > 0 else [f"Sample_{i+1}" for i in range(len(intensity_cols))]
                    
                    # Apply log transformation to intensity values
                    log_heatmap_data = []
                    for row in heatmap_array:
                        log_row = []
                        for val in row:
                            if val > 0:
                                log_val = np.log10(val)
                            else:
                                log_val = 0
                            log_row.append(log_val)
                        log_heatmap_data.append(log_row)
                    
                    fig.add_trace(
                        go.Heatmap(
                            z=log_heatmap_data,
                            x=sample_names_clustered,
                            y=taxon_labels_clustered,
                            colorscale='Viridis',
                            showscale=False,  # Remove color bar
                            name=f"{rank.capitalize()} Log Intensity Heatmap (Clustered)",
                            hovertemplate='<b>%{y}</b><br><b>Sample:</b> %{x}<br><b>Log Intensity:</b> %{z:.2f}<br><extra></extra>'
                        ),
                        row=current_row,
                        col=1
                    )

                current_row += 1

                # 3. Stacked Bar Chart - Composition by Sample
                # Calculate relative abundance for each sample
                sample_compositions = {}
                for i, col in enumerate(intensity_cols):
                    sample_data = rank_data[col]
                    total_sample = sample_data.sum()
                    if total_sample > 0:
                        # Get top 9 taxa for this sample (plus "Others")
                        sample_taxa_abundance = {}
                        for taxon in unique_taxa:
                            taxon_data = rank_data[rank_data[rank] == taxon]
                            taxon_abundance = taxon_data[col].sum()
                            if taxon_abundance > 0:
                                sample_taxa_abundance[taxon] = taxon_abundance
                        
                        top_sample_taxa = sorted(sample_taxa_abundance.items(), key=lambda x: x[1], reverse=True)[:9]
                        
                        # Calculate relative abundances
                        rel_abundances = []
                        taxon_names = []
                        for taxon_name, abundance in top_sample_taxa:
                            rel_abundances.append(abundance / total_sample * 100)
                            taxon_names.append(taxon_name)
                        
                        # Add "Others" category for remaining abundance
                        top_total = sum([ab for _, ab in top_sample_taxa])
                        others_abundance = (total_sample - top_total) / total_sample * 100
                        if others_abundance > 0:
                            rel_abundances.append(others_abundance)
                            taxon_names.append("Others")
                        
                        # Use File ID instead of long sample name
                        file_id = simple_names[i] if i < len(simple_names) else f"File_{i+1}"
                        sample_compositions[file_id] = {
                            'taxa': taxon_names,
                            'abundances': rel_abundances
                        }

                # Create stacked bar chart with distinct colors for each component
                if sample_compositions:
                    sample_names_list = list(sample_compositions.keys())
                    all_taxa = set()
                    for comp in sample_compositions.values():
                        all_taxa.update(comp['taxa'])
                    all_taxa = list(all_taxa)
                    
                    # Use distinct colors for each component (up to 10 components)
                    color_palette = [
                        '#B5BD00', '#8ED1FC', '#FF6B6B', '#4ECDC4', '#45B7D1', 
                        '#96CEB4', '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F'
                    ]
                    
                    # Create traces for each taxon
                    for i, taxon in enumerate(all_taxa):
                        abundances = []
                        for sample_name in sample_names_list:
                            if taxon in sample_compositions[sample_name]['taxa']:
                                idx = sample_compositions[sample_name]['taxa'].index(taxon)
                                abundances.append(sample_compositions[sample_name]['abundances'][idx])
                            else:
                                abundances.append(0)
                        
                        fig.add_trace(
                            go.Bar(
                                x=sample_names_list,
                                y=abundances,
                                name=taxon,
                                marker_color=color_palette[i % len(color_palette)],  # Each component gets distinct color
                                showlegend=False  # Hide legend for cleaner stacked bars
                            ),
                            row=current_row,
                            col=1
                        )

                current_row += 1

            # Update layout for the entire figure
            fig.update_layout(
                height=400 * len(available_ranks),  # Dynamic height based on number of ranks
                showlegend=False,
                title_text="Taxonomy Summary Report",
                barmode='stack'  # Enable stacked bars for composition charts
            )
            
            # Update axis labels for all rows
            for i in range(1, total_rows + 1, 3):
                fig.update_xaxes(title_text="Samples", row=i, col=1)
                fig.update_yaxes(title_text="Log Intensity", row=i, col=1)
                if i + 1 <= total_rows:
                    fig.update_xaxes(title_text="Samples", row=i+1, col=1)
                    fig.update_yaxes(title_text="Taxa", row=i+1, col=1)
                if i + 2 <= total_rows:
                    fig.update_xaxes(title_text="Samples", row=i+2, col=1)
                    fig.update_yaxes(title_text="Relative Abundance (%)", row=i+2, col=1)

            return fig
            
        except Exception as e:
            print(f"Error creating taxonomy summary: {e}")
            return self.create_empty_data_plot("Taxonomy Summary Report", f"Error generating report: {str(e)}")

    def create_function_summary(self, func_df):
        """Generate function summary report with visualizations."""
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Top 10 GO Terms',
                'Top 10 KEGG Pathways',
                'Top 10 COG Categories',
                'Top 10 EC Numbers'
            )
        )

        # GO terms distribution
        go_counts = func_df['GO_Term'].value_counts().head(10)
        fig.add_trace(
            go.Bar(
                x=go_counts.index,
                y=go_counts.values,
                name='GO Terms',
                marker_color=self.colors['primary']
            ),
            row=1, col=1
        )

        # KEGG pathways distribution
        kegg_counts = func_df['KEGG_Pathway'].value_counts().head(10)
        fig.add_trace(
            go.Bar(
                x=kegg_counts.index,
                y=kegg_counts.values,
                name='KEGG Pathways',
                marker_color=self.colors['secondary']
            ),
            row=1, col=2
        )

        # COG categories distribution
        cog_counts = func_df['COG_Category'].value_counts().head(10)
        fig.add_trace(
            go.Bar(
                x=cog_counts.index,
                y=cog_counts.values,
                name='COG Categories',
                marker_color=self.colors['tertiary']
            ),
            row=2, col=1
        )

        # EC numbers distribution
        ec_counts = func_df['EC_Number'].value_counts().head(10)
        fig.add_trace(
            go.Bar(
                x=ec_counts.index,
                y=ec_counts.values,
                name='EC Numbers',
                marker_color=self.colors['primary']
            ),
            row=2, col=2
        )

        # Update layout
        fig.update_layout(
            height=800,
            showlegend=False,
            template=self.template
        )

        return fig

    def generate_html_report(self, fig, title, output_file):
        """Generate HTML report with the given figure."""
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>{title}</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
            <style>
                body {{
                    background-color: {self.colors['background']};
                    color: {self.colors['text']};
                    padding: 20px;
                }}
                .report-header {{
                    text-align: center;
                    margin-bottom: 30px;
                }}
                .plot-container {{
                    background-color: white;
                    border-radius: 10px;
                    padding: 20px;
                    box-shadow: 0 0 10px rgba(0,0,0,0.1);
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <div class="report-header">
                    <h1>{title}</h1>
                    <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                </div>
                <div class="plot-container">
                    {fig.to_html(full_html=False, include_plotlyjs='cdn')}
                </div>
            </div>
        </body>
        </html>
        """
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

    def generate_all_reports(self, data_dir):
        """Generate all reports from the data directory."""
        # Read data files
        peptides_df = pd.read_csv(os.path.join(data_dir, 'final_peptides.tsv'), sep='\t')
        proteins_path = os.path.join(data_dir, 'razor_proteins.tsv')
        if not os.path.exists(proteins_path):
            proteins_path = os.path.join(data_dir, 'final_proteins.tsv')
        proteins_df = pd.read_csv(proteins_path, sep='\t')
        summary_df = pd.read_csv(os.path.join(data_dir, 'final_summary.tsv'), sep='\t')
        taxa_df = pd.read_csv(os.path.join(data_dir, 'Taxa_report.tsv'), sep='\t')
        func_df = pd.read_csv(os.path.join(data_dir, 'functions_report.tsv'), sep='\t')

        # Generate reports
        self.generate_html_report(
            self.create_peptides_summary(peptides_df),
            'Peptides Summary Report',
            os.path.join(self.output_dir, 'report_peptides_summary.html'),
            summary='Peptide identification coverage, quality metrics, and quantitation trends.',
            tags=['peptides', 'quantitation']
        )

        self.generate_html_report(
            self.create_proteins_summary(proteins_df),
            'Proteins Summary Report',
            os.path.join(self.output_dir, 'report_proteins_summary.html'),
            summary='Protein group abundance, coverage, and differential analysis highlights.',
            tags=['proteins']
        )

        self.generate_html_report(
            self.create_id_summary(summary_df),
            'ID Summary Report',
            os.path.join(self.output_dir, 'report_ID_summary.html'),
            summary='Run-level identification overview with search engine performance metrics.',
            tags=['identifications']
        )

        self.generate_html_report(
            self.create_taxonomy_summary(taxa_df),
            'Taxonomy Summary Report',
            os.path.join(self.output_dir, 'report_taxonomy_summary.html'),
            summary='Abundance and diversity of taxa across ranks with multivariate exploration.',
            tags=['taxonomy']
        )

        self.generate_html_report(
            self.create_function_summary(func_df),
            'Function Summary Report',
            os.path.join(self.output_dir, 'report_function_summary.html'),
            summary='Functional annotation, pathway enrichment, and ontology-level breakdowns.',
            tags=['function', 'ontology']
        )

if __name__ == '__main__':
    # Example usage
    generator = ReportGenerator()
    generator.generate_all_reports('src/test/java/org/metapilot/data') 
