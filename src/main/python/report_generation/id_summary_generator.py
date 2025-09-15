import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from .base_generator import BaseReportGenerator

class IDSummaryGenerator(BaseReportGenerator):
    """Generator for ID summary reports."""
    
    def create_id_summary(self, summary_df, metadata=None):
        """Generate ID summary report with metadata support."""
        # Validate data before processing
        if not self.validate_data(summary_df, "ID summary"):
            return self.create_empty_data_plot(
                "ID Summary Report",
                "No ID data available for analysis.<br>Please check if the analysis completed successfully."
            )
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Total Identifications by Sample',
                'PSM Distribution', 
                'Peptide Distribution',
                'Protein Distribution'
            )
        )

        sample_cols = [col for col in summary_df.columns if not col.startswith('Unnamed')]
        sample_names = sample_cols[:10]
        
        # Extract group information
        groups = None
        group_colors = None
        if metadata is not None:
            try:
                sample_col = 'Sample' if 'Sample' in metadata.columns else metadata.columns[0]
                group_col = 'Group' if 'Group' in metadata.columns else None
                
                if group_col:
                    sample_groups = dict(zip(metadata[sample_col], metadata[group_col]))
                    groups = [sample_groups.get(sample, 'Unknown') for sample in sample_names]
                    unique_groups = list(set(groups))
                    group_colors = {group: self.group_colors[i % len(self.group_colors)] 
                                  for i, group in enumerate(unique_groups)}
            except:
                pass
        
        if len(sample_names) > 0:
            # Color by groups if available
            colors = [group_colors[groups[i]] if groups and group_colors else self.colors['primary'] 
                     for i in range(len(sample_names))]
            
            # 1. Total Identifications by Sample
            total_ids = [summary_df[sample].sum() if pd.api.types.is_numeric_dtype(summary_df[sample]) else 0 
                        for sample in sample_names]
            fig.add_trace(go.Bar(x=sample_names, y=total_ids, name='Total IDs', 
                               marker_color=colors), row=1, col=1)
            
            # 2. PSM Distribution - try to find PSM-related columns
            psm_data = []
            for sample in sample_names:
                # Look for PSM count columns
                psm_cols = [col for col in summary_df.columns if 'psm' in col.lower() and sample in col.lower()]
                if psm_cols:
                    psm_count = summary_df[psm_cols[0]].sum() if pd.api.types.is_numeric_dtype(summary_df[psm_cols[0]]) else 0
                else:
                    # Fallback: use total IDs as proxy for PSMs
                    psm_count = total_ids[sample_names.index(sample)]
                psm_data.append(psm_count)
            
            fig.add_trace(go.Bar(x=sample_names, y=psm_data, name='PSMs', 
                               marker_color=colors), row=1, col=2)
            
            # 3. Peptide Distribution - try to find peptide-related columns
            peptide_data = []
            for sample in sample_names:
                # Look for peptide count columns
                peptide_cols = [col for col in summary_df.columns if 'peptide' in col.lower() and sample in col.lower()]
                if peptide_cols:
                    peptide_count = summary_df[peptide_cols[0]].sum() if pd.api.types.is_numeric_dtype(summary_df[peptide_cols[0]]) else 0
                else:
                    # Fallback: use total IDs as proxy for peptides
                    peptide_count = total_ids[sample_names.index(sample)]
                peptide_data.append(peptide_count)
            
            fig.add_trace(go.Bar(x=sample_names, y=peptide_data, name='Peptides', 
                               marker_color=colors), row=2, col=1)
            
            # 4. Protein Distribution - try to find protein-related columns
            protein_data = []
            for sample in sample_names:
                # Look for protein count columns
                protein_cols = [col for col in summary_df.columns if 'protein' in col.lower() and sample in col.lower()]
                if protein_cols:
                    protein_count = summary_df[protein_cols[0]].sum() if pd.api.types.is_numeric_dtype(summary_df[protein_cols[0]]) else 0
                else:
                    # Fallback: use total IDs as proxy for proteins
                    protein_count = total_ids[sample_names.index(sample)]
                protein_data.append(protein_count)
            
            fig.add_trace(go.Bar(x=sample_names, y=protein_data, name='Proteins', 
                               marker_color=colors), row=2, col=2)

        fig.update_layout(height=800, showlegend=False, title_text="ID Summary Report")
        return fig
