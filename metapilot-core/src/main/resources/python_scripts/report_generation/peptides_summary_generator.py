import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from .base_generator import BaseReportGenerator

class PeptidesSummaryGenerator(BaseReportGenerator):
    """Generator for peptides summary reports."""
    
    def create_peptides_summary(self, peptides_df, metadata=None):
        """Enhanced peptides summary with metadata support."""
        # Validate data before processing
        if not self.validate_data(peptides_df, "peptides"):
            return self.create_empty_data_plot(
                "Peptides Summary Report",
                "No peptides data available for analysis.<br>Please check if the peptide identification completed successfully."
            )
        
        fig = make_subplots(
            rows=5, cols=2,
            subplot_titles=(
                'Peptide Length Distribution',
                'Missed Cleavages Distribution', 
                'Charge State Distribution',
                'Mass Distribution (Da)',
                'Peptide Counts per Sample',
                'Unique Sequence Counts per Sample',
                'MS/MS Count Distribution by Sample',
                'Intensity Distribution by Sample',
                'Peptide Intensity Heatmap (All Peptides by Total Intensity)',
                ''
            ),
            specs=[[{"type": "bar"}, {"type": "bar"}],  # Pre-binned histograms (rendered as bar charts)
                   [{"type": "bar"}, {"type": "bar"}],  # Pre-binned histograms (rendered as bar charts)
                   [{"type": "bar"}, {"type": "bar"}],
                   [{"type": "bar"}, {"type": "bar"}],
                   [{"type": "heatmap", "colspan": 2}, None]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(peptides_df, metadata)

        # Use experiment names directly (no file ID mapping needed)
        simple_names = sample_names

        print(f"INFO: peptides report: {len(peptides_df)} peptides, {len(intensity_cols)} sample columns",
              flush=True)

        # Peptide counts per sample (column-wise on numpy arrays; no DataFrame
        # copies per sample, which cost minutes with thousands of TMT channels)
        if len(intensity_cols) > 0:
            peptide_counts = []
            unique_counts = []
            sequences = peptides_df['Sequence'].to_numpy() if 'Sequence' in peptides_df.columns else None
            sequences_unique = sequences is not None and peptides_df['Sequence'].is_unique

            for col in intensity_cols:
                present = np.nan_to_num(peptides_df[col].to_numpy(dtype='float64'), nan=0.0) > 0
                peptide_count = int(present.sum())
                peptide_counts.append(peptide_count)

                if sequences is None or sequences_unique:
                    unique_counts.append(peptide_count)
                else:
                    unique_counts.append(int(pd.unique(sequences[present]).size))
            
            # Color by groups if available
            colors = [group_colors[groups[i]] if groups and group_colors else self.colors['primary'] 
                     for i in range(len(sample_names))]
            
            fig.add_trace(go.Bar(x=simple_names, y=peptide_counts, name='Total Peptides',
                               marker_color=colors), row=3, col=1)
            
            fig.add_trace(go.Bar(x=simple_names, y=unique_counts, name='Unique Sequences',
                               marker_color=colors), row=3, col=2)

        # Store mapping in figure for potential use (if mapping exists)
        if len(simple_names) > 0:
            fig.update_layout(
                annotations=[
                    dict(
                        text=f"Sample Names: {', '.join(simple_names[:5])}{'...' if len(simple_names) > 5 else ''}",
                        xref="paper", yref="paper",
                        x=0, y=-0.1, xanchor='left', yanchor='top',
                        showarrow=False, font=dict(size=10),
                        visible=False  # Hidden but stored
                    )
                ]
            )

        # 1. Peptide Length Distribution (row 1, col 1)
        length_series = None
        if 'Length' in peptides_df.columns:
            length_series = peptides_df['Length'].dropna()
        elif 'Sequence' in peptides_df.columns:
            length_series = peptides_df['Sequence'].str.len().dropna()

        if length_series is not None and len(length_series) > 0:
            length_counts = length_series.astype(int).value_counts().sort_index()
            fig.add_trace(
                go.Bar(
                    x=length_counts.index.astype(int).tolist(),
                    y=length_counts.values.tolist(),
                    name='Length Distribution',
                    marker_color=self.colors['tertiary']
                ),
                row=1,
                col=1
            )
        else:
            fig.add_annotation(text="No length or sequence data", xref="paper", yref="paper",
                             x=0.25, y=0.9, showarrow=False, row=1, col=1)

        # 2. Missed Cleavages Distribution (row 1, col 2)
        # For discrete data like missed cleavages, use actual unique values as bins
        mc_cols = ['Missed cleavages', 'Missed Cleavages', 'missed_cleavages']
        mc_col = next((col for col in mc_cols if col in peptides_df.columns), None)
        
        if mc_col:
            mc_data = peptides_df[mc_col].dropna()
            if len(mc_data) > 0:
                # For discrete integer data, count occurrences of each unique value
                value_counts = mc_data.value_counts().sort_index()
                fig.add_trace(go.Bar(x=value_counts.index.tolist(), y=value_counts.values.tolist(), 
                                   name='Missed Cleavages',
                                   marker_color=self.colors['secondary']), row=1, col=2)
            else:
                fig.add_annotation(text="No missed cleavages data", xref="paper", yref="paper", 
                                 x=0.75, y=0.9, showarrow=False, row=1, col=2)
        else:
            fig.add_annotation(text="Missed cleavages column not found", xref="paper", yref="paper", 
                             x=0.75, y=0.9, showarrow=False, row=1, col=2)

        # 3. Charge State Distribution (row 2, col 1)
        # Use pre-binned data to avoid binary encoding issues with large datasets
        charge_cols = ['Charges', 'Charge', 'charge']
        charge_col = next((col for col in charge_cols if col in peptides_df.columns), None)
        
        if charge_col:
            if peptides_df[charge_col].dtype == 'object':
                # Handle semicolon-separated charge states like "2;3"
                all_charges = []
                for charge_str in peptides_df[charge_col].astype(str):
                    if ';' in charge_str:
                        # Split by semicolon and add all charge states
                        individual_charges = charge_str.split(';')
                        for ch in individual_charges:
                            try:
                                all_charges.append(int(ch.strip()))
                            except (ValueError, AttributeError):
                                continue
                    else:
                        # Single charge state
                        try:
                            all_charges.append(int(charge_str.strip()))
                        except (ValueError, AttributeError):
                            continue
                charges = pd.Series(all_charges).dropna()
            else:
                charges = peptides_df[charge_col].dropna()
            
            if len(charges) > 0:
                # For discrete integer data, count occurrences of each unique value
                value_counts = charges.value_counts().sort_index()
                fig.add_trace(go.Bar(x=value_counts.index.tolist(), y=value_counts.values.tolist(),
                                   name='Charge States',
                                   marker_color=self.colors['success']), row=2, col=1)
            else:
                fig.add_annotation(text="No charge data available", xref="paper", yref="paper", 
                                 x=0.25, y=0.65, showarrow=False, row=2, col=1)
        else:
            fig.add_annotation(text="Charge column not found", xref="paper", yref="paper", 
                             x=0.25, y=0.65, showarrow=False, row=2, col=1)

        # 4. Mass Distribution (row 2, col 2)
        # Use pre-binned data to avoid binary encoding issues with large datasets
        mass_cols = ['Mass', 'Molecular Weight', 'mass']
        mass_col = next((col for col in mass_cols if col in peptides_df.columns), None)
        
        if mass_col:
            mass_data = peptides_df[mass_col].dropna()
            if len(mass_data) > 0:
                # Determine optimal number of bins based on data range
                data_range = mass_data.max() - mass_data.min()
                num_bins = min(40, max(15, int(data_range / 100)))  # Between 15-40 bins
                counts, bins = np.histogram(mass_data, bins=num_bins)
                bin_centers = ((bins[:-1] + bins[1:]) / 2).tolist()
                counts_list = counts.tolist()
                fig.add_trace(go.Bar(x=bin_centers, y=counts_list, name='Mass Distribution',
                                   marker_color=self.colors['warning']), row=2, col=2)
            else:
                fig.add_annotation(text="No mass data available", xref="paper", yref="paper", 
                                 x=0.75, y=0.65, showarrow=False, row=2, col=2)
        else:
            fig.add_annotation(text="Mass column not found", xref="paper", yref="paper", 
                             x=0.75, y=0.65, showarrow=False, row=2, col=2)

        # 5. MS/MS Count Distribution by Sample (row 4, col 1) - Bar chart
        # MS/MS counts are per raw file. Label them with the sample names only when
        # the sample columns are also per raw file (label-free); with isobaric
        # labels there are several channels per raw file, so use the file names.
        msms_count_cols = [col for col in peptides_df.columns if 'MS/MS count' in col]
        if msms_count_cols and len(simple_names) > 0:
            msms_totals = [float(peptides_df[col].sum()) for col in msms_count_cols]
            if len(msms_count_cols) == len(intensity_cols):
                msms_labels = simple_names[:len(msms_totals)]
            else:
                msms_labels = [col.replace('MS/MS count', '', 1).strip() for col in msms_count_cols]

            if msms_totals:
                fig.add_trace(go.Bar(x=msms_labels, y=msms_totals, name='MS/MS Count Total',
                                      marker_color=self.colors['primary']), row=4, col=1)

        # 6. Intensity Distribution by Sample (row 4, col 2) - Bar chart
        intensity_sample_cols = intensity_cols
        if intensity_sample_cols and len(simple_names) > 0:
            intensity_totals = [float(peptides_df[col].to_numpy().sum(dtype='float64'))
                                for col in intensity_sample_cols]
            
            if intensity_totals:
                fig.add_trace(go.Bar(x=simple_names[:len(intensity_totals)], y=intensity_totals, name='Intensity Total',
                                      marker_color=self.colors['secondary']), row=4, col=2)

        # 9. Peptide Intensity Heatmap (row 5, spanning both columns): top peptides
        # by total intensity, limited to HEATMAP_CELL_BUDGET cells (all peptides
        # when they fit). Writing every peptide x sample cell made TMT reports
        # (hundreds of millions of cells) run for hours.
        heatmap_rows_shown = None
        if 'Intensity' in peptides_df.columns and len(intensity_sample_cols) > 0:
            heatmap_data, heatmap_rows_shown, heatmap_rows_total = self.top_rows_log_matrix(
                peptides_df, intensity_sample_cols, 'Intensity')
            if heatmap_rows_shown < heatmap_rows_total:
                self.set_subplot_title(fig, 'Peptide Intensity Heatmap',
                                       f'Peptide Intensity Heatmap (Top {heatmap_rows_shown:,} of '
                                       f'{heatmap_rows_total:,} Peptides by Total Intensity)')
            print(f"INFO: peptide heatmap: {heatmap_rows_shown} of {heatmap_rows_total} peptides x "
                  f"{len(intensity_sample_cols)} samples", flush=True)

            # Create heatmap with file IDs
            heatmap_sample_names = simple_names[:len(intensity_sample_cols)] if len(simple_names) > 0 else [f"Sample_{i+1}" for i in range(len(intensity_sample_cols))]
            
            fig.add_trace(go.Heatmap(
                z=heatmap_data,
                x=heatmap_sample_names,
                y=None,  # Don't show peptide labels - too messy with all peptides
                colorscale='Viridis',
                showscale=False,  # Remove color bar
                name='Log Intensity Heatmap',
                hovertemplate='<b>Sample:</b> %{x}<br><b>Log Intensity:</b> %{z:.2f}<br><extra></extra>'
            ), row=5, col=1)

        fig.update_layout(height=1800, showlegend=False, title_text="MetaPilot Peptides Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Peptide Length (aa)", row=1, col=1)
        fig.update_xaxes(title_text="Missed Cleavages", row=1, col=2)
        fig.update_xaxes(title_text="Charge State", row=2, col=1)
        fig.update_xaxes(title_text="Mass (Da)", row=2, col=2)
        fig.update_xaxes(title_text="Samples", row=3, col=1)
        fig.update_xaxes(title_text="Samples", row=3, col=2)
        fig.update_xaxes(title_text="Samples", row=4, col=1)
        fig.update_xaxes(title_text="Samples", row=4, col=2)
        fig.update_xaxes(title_text="Samples", row=5, col=1)
        
        fig.update_yaxes(title_text="Count", row=1, col=1)
        fig.update_yaxes(title_text="Count", row=1, col=2)
        fig.update_yaxes(title_text="Count", row=2, col=1)
        fig.update_yaxes(title_text="Count", row=2, col=2)
        fig.update_yaxes(title_text="Peptide Count", row=3, col=1)
        fig.update_yaxes(title_text="Unique Sequences", row=3, col=2)
        fig.update_yaxes(title_text="Total MS/MS Count", row=4, col=1)
        fig.update_yaxes(title_text="Total Intensity", row=4, col=2)
        fig.update_yaxes(title_text="Peptides (by Total Intensity)", row=5, col=1)

        return fig

