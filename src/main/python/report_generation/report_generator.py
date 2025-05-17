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
        """Generate taxonomy summary report with visualizations."""
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Top 10 Phyla',
                'Top 10 Classes',
                'Top 10 Orders',
                'Top 10 Families'
            )
        )

        # Phyla distribution
        phyla_counts = taxa_df['Phylum'].value_counts().head(10)
        fig.add_trace(
            go.Bar(
                x=phyla_counts.index,
                y=phyla_counts.values,
                name='Phyla',
                marker_color=self.colors['primary']
            ),
            row=1, col=1
        )

        # Class distribution
        class_counts = taxa_df['Class'].value_counts().head(10)
        fig.add_trace(
            go.Bar(
                x=class_counts.index,
                y=class_counts.values,
                name='Classes',
                marker_color=self.colors['secondary']
            ),
            row=1, col=2
        )

        # Order distribution
        order_counts = taxa_df['Order'].value_counts().head(10)
        fig.add_trace(
            go.Bar(
                x=order_counts.index,
                y=order_counts.values,
                name='Orders',
                marker_color=self.colors['tertiary']
            ),
            row=2, col=1
        )

        # Family distribution
        family_counts = taxa_df['Family'].value_counts().head(10)
        fig.add_trace(
            go.Bar(
                x=family_counts.index,
                y=family_counts.values,
                name='Families',
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
        proteins_df = pd.read_csv(os.path.join(data_dir, 'final_proteins.tsv'), sep='\t')
        summary_df = pd.read_csv(os.path.join(data_dir, 'final_summary.tsv'), sep='\t')
        taxa_df = pd.read_csv(os.path.join(data_dir, 'Taxa_report.tsv'), sep='\t')
        func_df = pd.read_csv(os.path.join(data_dir, 'functions_report.tsv'), sep='\t')

        # Generate reports
        self.generate_html_report(
            self.create_peptides_summary(peptides_df),
            'Peptides Summary Report',
            os.path.join(self.output_dir, 'report_peptides_summary.html')
        )

        self.generate_html_report(
            self.create_proteins_summary(proteins_df),
            'Proteins Summary Report',
            os.path.join(self.output_dir, 'report_proteins_summary.html')
        )

        self.generate_html_report(
            self.create_id_summary(summary_df),
            'ID Summary Report',
            os.path.join(self.output_dir, 'report_ID_summary.html')
        )

        self.generate_html_report(
            self.create_taxonomy_summary(taxa_df),
            'Taxonomy Summary Report',
            os.path.join(self.output_dir, 'report_taxonomy_summary.html')
        )

        self.generate_html_report(
            self.create_function_summary(func_df),
            'Function Summary Report',
            os.path.join(self.output_dir, 'report_function_summary.html')
        )

if __name__ == '__main__':
    # Example usage
    generator = ReportGenerator()
    generator.generate_all_reports('src/test/java/org/metapilot/data') 