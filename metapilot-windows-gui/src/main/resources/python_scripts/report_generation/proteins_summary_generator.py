import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from .base_generator import BaseReportGenerator

try:
    from umap import UMAP
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False

class ProteinsSummaryGenerator(BaseReportGenerator):
    """Generator for proteins summary reports."""
    
    def create_proteins_summary(self, proteins_df, metadata=None):
        """Enhanced proteins summary with advanced dimensionality reduction."""
        # Validate data before processing
        if not self.validate_data(proteins_df, "proteins"):
            return self.create_empty_data_plot(
                "Proteins Summary Report",
                "No proteins data available for analysis.<br>Please check if the protein identification completed successfully."
            )
        
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=(
                'Peptide Counts (All) Distribution',
                'Peptide Counts (Razor) Distribution',
                'Number of PSMs Distribution',
                'Intensity Distribution',
                'Protein Intensity Heatmap (All Proteins by Total Intensity)',
                ''
            ),
            specs=[[{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 2}, None]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(proteins_df, metadata)

        # Use experiment names directly (no file ID mapping needed)
        simple_names = sample_names

        # 1. Peptide Counts (All) Distribution by Sample (row 1, col 1)
        if 'Peptide counts (all)' in proteins_df.columns and len(intensity_cols) > 0:
            all_peptide_data = []
            sample_labels = []
            for i, col in enumerate(intensity_cols):
                # Use File ID instead of long sample name
                file_id = simple_names[i] if i < len(simple_names) else f"File_{i+1}"
                # Get peptide counts for proteins with intensity > 0 in this sample
                sample_proteins = proteins_df[proteins_df[col] > 0]
                if len(sample_proteins) > 0:
                    peptide_counts = sample_proteins['Peptide counts (all)'].tolist()
                    all_peptide_data.extend(peptide_counts)
                    sample_labels.extend([file_id] * len(peptide_counts))
            
            if all_peptide_data:
                fig.add_trace(go.Violin(x=sample_labels, y=all_peptide_data, name='Peptide Counts (All)',
                                      box_visible=True, meanline_visible=True,
                                      points='all', pointpos=0, jitter=0.3,
                                      marker_color=self.colors['primary']), row=1, col=1)

        # 2. Peptide Counts (Razor) Distribution by Sample (row 1, col 2)
        if 'Peptide counts (razor)' in proteins_df.columns and len(intensity_cols) > 0:
            razor_peptide_data = []
            sample_labels = []
            for i, col in enumerate(intensity_cols):
                # Use File ID instead of long sample name
                file_id = simple_names[i] if i < len(simple_names) else f"File_{i+1}"
                # Get peptide counts for proteins with intensity > 0 in this sample
                sample_proteins = proteins_df[proteins_df[col] > 0]
                if len(sample_proteins) > 0:
                    peptide_counts = sample_proteins['Peptide counts (razor)'].tolist()
                    razor_peptide_data.extend(peptide_counts)
                    sample_labels.extend([file_id] * len(peptide_counts))
            
            if razor_peptide_data:
                fig.add_trace(go.Violin(x=sample_labels, y=razor_peptide_data, name='Peptide Counts (Razor)',
                                      box_visible=True, meanline_visible=True,
                                      points='all', pointpos=0, jitter=0.3,
                                      marker_color=self.colors['secondary']), row=1, col=2)

        # 3. Number of PSMs Distribution by Sample (row 2, col 1)
        if 'Number of PSMs' in proteins_df.columns and len(intensity_cols) > 0:
            psm_data = []
            sample_labels = []
            for i, col in enumerate(intensity_cols):
                # Use File ID instead of long sample name
                file_id = simple_names[i] if i < len(simple_names) else f"File_{i+1}"
                # Get PSM counts for proteins with intensity > 0 in this sample
                sample_proteins = proteins_df[proteins_df[col] > 0]
                if len(sample_proteins) > 0:
                    psm_counts = sample_proteins['Number of PSMs'].tolist()
                    psm_data.extend(psm_counts)
                    sample_labels.extend([file_id] * len(psm_counts))
            
            if psm_data:
                fig.add_trace(go.Violin(x=sample_labels, y=psm_data, name='Number of PSMs',
                                      box_visible=True, meanline_visible=True,
                                      points='all', pointpos=0, jitter=0.3,
                                      marker_color=self.colors['tertiary']), row=2, col=1)

        # 4. Intensity Distribution by Sample (row 2, col 2)
        if len(intensity_cols) > 0:
            intensity_data = []
            sample_labels = []
            for i, col in enumerate(intensity_cols):
                # Use File ID instead of long sample name
                file_id = simple_names[i] if i < len(simple_names) else f"File_{i+1}"
                # Get intensity values > 0 for this sample
                non_zero_intensities = proteins_df[proteins_df[col] > 0][col].tolist()
                if non_zero_intensities:
                    intensity_data.extend(non_zero_intensities)
                    sample_labels.extend([file_id] * len(non_zero_intensities))
            
            if intensity_data:
                fig.add_trace(go.Violin(x=sample_labels, y=intensity_data, name='Intensity',
                                      box_visible=True, meanline_visible=True,
                                      points='all', pointpos=0, jitter=0.3,
                                      marker_color=self.colors['success']), row=2, col=2)

        # 5. Protein Intensity Heatmap (row 3, spanning both columns)
        if 'Intensity' in proteins_df.columns and len(intensity_cols) > 0:
            # Use all proteins ordered by total intensity
            all_proteins = proteins_df.sort_values('Intensity', ascending=False)
            
            # Prepare data for heatmap with log transformation
            heatmap_data = []
            
            for _, row in all_proteins.iterrows():
                # Get intensity values for each sample
                sample_intensities = []
                for col in intensity_cols:
                    intensity_val = row[col] if pd.notna(row[col]) else 0
                    # Apply log transformation but keep 0 as 0
                    if intensity_val > 0:
                        log_intensity = np.log10(intensity_val)
                    else:
                        log_intensity = 0
                    sample_intensities.append(log_intensity)
                heatmap_data.append(sample_intensities)
            
            # Create heatmap with File IDs
            heatmap_sample_names = simple_names[:len(intensity_cols)] if len(simple_names) > 0 else [f"Sample_{i+1}" for i in range(len(intensity_cols))]
            
            fig.add_trace(go.Heatmap(
                z=heatmap_data,
                x=heatmap_sample_names,
                y=None,  # Don't show protein labels - too messy with all proteins
                colorscale='Viridis',
                showscale=False,  # Remove color bar
                name='Log Protein Intensity Heatmap',
                hovertemplate='<b>Sample:</b> %{x}<br><b>Log Intensity:</b> %{z:.2f}<br><extra></extra>'
            ), row=3, col=1)

        fig.update_layout(height=1200, showlegend=False, title_text="MetaPilot Proteins Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Samples", row=1, col=1)
        fig.update_xaxes(title_text="Samples", row=1, col=2)
        fig.update_xaxes(title_text="Samples", row=2, col=1)
        fig.update_xaxes(title_text="Samples", row=2, col=2)
        fig.update_xaxes(title_text="Samples", row=3, col=1)
        
        fig.update_yaxes(title_text="Peptide Counts (All)", row=1, col=1)
        fig.update_yaxes(title_text="Peptide Counts (Razor)", row=1, col=2)
        fig.update_yaxes(title_text="Number of PSMs", row=2, col=1)
        fig.update_yaxes(title_text="Intensity", row=2, col=2)
        fig.update_yaxes(title_text="All Proteins (by Total Intensity)", row=3, col=1)

        return fig
