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

        print(f"INFO: proteins report: {len(proteins_df)} proteins, {len(intensity_cols)} sample columns",
              flush=True)

        sample_labels = [simple_names[i] if i < len(simple_names) else f"File_{i+1}"
                         for i in range(len(intensity_cols))]
        # Presence mask per sample, computed once from numpy arrays. The previous
        # code filtered the whole DataFrame once per sample and plot (copying every
        # column each time) and drew every value as a violin point, which with
        # thousands of TMT channels does not finish.
        presence = [np.nan_to_num(proteins_df[col].to_numpy(dtype='float64'), nan=0.0) > 0
                    for col in intensity_cols]

        def per_sample(values_column):
            values = proteins_df[values_column].to_numpy()
            return [values[mask] for mask in presence]

        # 1. Peptide Counts (All) Distribution by Sample (row 1, col 1)
        if 'Peptide counts (all)' in proteins_df.columns and len(intensity_cols) > 0:
            self.add_sample_distribution(fig, sample_labels, per_sample('Peptide counts (all)'),
                                         'Peptide Counts (All)', self.colors['primary'], 1, 1)

        # 2. Peptide Counts (Razor) Distribution by Sample (row 1, col 2)
        if 'Peptide counts (razor)' in proteins_df.columns and len(intensity_cols) > 0:
            self.add_sample_distribution(fig, sample_labels, per_sample('Peptide counts (razor)'),
                                         'Peptide Counts (Razor)', self.colors['secondary'], 1, 2)

        # 3. Number of PSMs Distribution by Sample (row 2, col 1)
        if 'Number of PSMs' in proteins_df.columns and len(intensity_cols) > 0:
            self.add_sample_distribution(fig, sample_labels, per_sample('Number of PSMs'),
                                         'Number of PSMs', self.colors['tertiary'], 2, 1)

        # 4. Intensity Distribution by Sample (row 2, col 2)
        if len(intensity_cols) > 0:
            intensity_values = [proteins_df[col].to_numpy()[mask]
                                for col, mask in zip(intensity_cols, presence)]
            self.add_sample_distribution(fig, sample_labels, intensity_values,
                                         'Intensity', self.colors['success'], 2, 2)

        # 5. Protein Intensity Heatmap (row 3, spanning both columns): top proteins
        # by total intensity, limited to HEATMAP_CELL_BUDGET cells (all proteins
        # when they fit).
        if 'Intensity' in proteins_df.columns and len(intensity_cols) > 0:
            heatmap_data, rows_shown, rows_total = self.top_rows_log_matrix(proteins_df, intensity_cols,
                                                                            'Intensity')
            if rows_shown < rows_total:
                self.set_subplot_title(fig, 'Protein Intensity Heatmap',
                                       f'Protein Intensity Heatmap (Top {rows_shown:,} of {rows_total:,} '
                                       f'Proteins by Total Intensity)')
            print(f"INFO: protein heatmap: {rows_shown} of {rows_total} proteins x {len(intensity_cols)} samples",
                  flush=True)

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
        fig.update_yaxes(title_text="Proteins (by Total Intensity)", row=3, col=1)

        return fig
