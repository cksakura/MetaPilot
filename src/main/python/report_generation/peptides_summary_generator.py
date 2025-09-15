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
            rows=3, cols=2,
            subplot_titles=(
                'Peptide Counts per Raw File',
                'Unique Sequence Counts per Raw File',
                'Peptide Length Distribution',
                'Charge State Distribution',
                'Mass Distribution (Da)',
                'Missed Cleavages Distribution'
            ),
            specs=[[{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(peptides_df, metadata)

        # Peptide counts per raw file
        if len(intensity_cols) > 0:
            peptide_counts = []
            unique_counts = []
            
            for col in intensity_cols:
                peptide_count = (peptides_df[col] > 0).sum()
                peptide_counts.append(peptide_count)
                
                if 'Sequence' in peptides_df.columns:
                    unique_seqs = peptides_df[peptides_df[col] > 0]['Sequence'].nunique()
                    unique_counts.append(unique_seqs)
                else:
                    unique_counts.append(peptide_count)
            
            # Color by groups if available
            colors = [group_colors[groups[i]] if groups and group_colors else self.colors['primary'] 
                     for i in range(len(sample_names))]
            
            fig.add_trace(go.Bar(x=sample_names, y=peptide_counts, name='Total Peptides',
                               marker_color=colors), row=1, col=1)
            
            fig.add_trace(go.Bar(x=sample_names, y=unique_counts, name='Unique Sequences',
                               marker_color=colors), row=1, col=2)

        # Statistical distributions
        if 'Sequence' in peptides_df.columns:
            peptides_df['Length'] = peptides_df['Sequence'].str.len()
            fig.add_trace(go.Histogram(x=peptides_df['Length'], name='Length Distribution',
                                     marker_color=self.colors['tertiary'], nbinsx=30), row=2, col=1)

        charge_cols = ['Charges', 'Charge', 'charge']
        charge_col = next((col for col in charge_cols if col in peptides_df.columns), None)
        
        if charge_col:
            if peptides_df[charge_col].dtype == 'object':
                charges = peptides_df[charge_col].astype(str).str.split(';').str[0]
                charges = pd.to_numeric(charges, errors='coerce')
            else:
                charges = peptides_df[charge_col]
                
            fig.add_trace(go.Histogram(x=charges, name='Charge States',
                                     marker_color=self.colors['success'], nbinsx=8), row=2, col=2)

        mass_cols = ['Mass', 'Molecular Weight', 'mass']
        mass_col = next((col for col in mass_cols if col in peptides_df.columns), None)
        
        if mass_col:
            fig.add_trace(go.Histogram(x=peptides_df[mass_col], name='Mass Distribution',
                                     marker_color=self.colors['light_blue'], nbinsx=40), row=3, col=1)

        mc_cols = ['Missed cleavages', 'Missed Cleavages', 'missed_cleavages']
        mc_col = next((col for col in mc_cols if col in peptides_df.columns), None)
        
        if mc_col:
            fig.add_trace(go.Histogram(x=peptides_df[mc_col], name='Missed Cleavages',
                                     marker_color=self.colors['light_purple'], nbinsx=5), row=3, col=2)

        fig.update_layout(height=1200, showlegend=False, title_text="MetaPilot Peptides Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Raw Files", row=1, col=1)
        fig.update_xaxes(title_text="Raw Files", row=1, col=2)
        fig.update_yaxes(title_text="Count", row=1, col=1)
        fig.update_yaxes(title_text="Count", row=1, col=2)

        return fig
