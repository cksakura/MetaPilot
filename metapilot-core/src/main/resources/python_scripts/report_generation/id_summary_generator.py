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
            rows=4, cols=2,
            subplot_titles=(
                'MS/MS Scans by Sample',
                'MS/MS Identified (PSM Count)', 
                'Identification Rate (%)',
                'Peptide Sequences Count',
                'Distribution of MS/MS Scans',
                'Distribution of PSM Counts',
                'Distribution of Identification Rates',
                'Distribution of Peptide Counts'
            )
        )

        # Check if this is a MetaPilot summary format (with Raw file, Experiment columns)
        if 'Raw file' in summary_df.columns and 'Experiment' in summary_df.columns:
            # Filter out any rows with missing or invalid data (like "Total" rows)
            valid_rows = summary_df[
                (summary_df['Raw file'].notna()) & 
                (summary_df['Raw file'] != 'Total') &
                (summary_df['Experiment'].notna()) &
                (summary_df['Experiment'] != '')
            ].copy()
            
            if len(valid_rows) == 0:
                return self.create_empty_data_plot(
                    "ID Summary Report",
                    "No valid sample data found in the summary file."
                )
            
            # MetaPilot format: each row is a sample
            raw_sample_names = valid_rows['Raw file'].tolist()
            
            # Map raw file names to experiment names using metadata
            sample_names = raw_sample_names
            if metadata is not None and len(raw_sample_names) > 0:
                try:
                    # Create mapping from raw file names to experiment (source) names
                    raw_to_exp = {}
                    if 'comment[data file]' in metadata.columns and 'source name' in metadata.columns:
                        for _, row in metadata.iterrows():
                            raw_name = row['comment[data file]']
                            exp_name = row['source name']
                            raw_to_exp[raw_name] = exp_name
                    
                    # Map sample names to experiment names
                    sample_names = []
                    for raw_name in raw_sample_names:
                        # Try exact match first
                        if raw_name in raw_to_exp:
                            sample_names.append(raw_to_exp[raw_name])
                        else:
                            # Try partial match (in case of slight differences)
                            matched = False
                            for raw_key, exp_name in raw_to_exp.items():
                                if raw_name in raw_key or raw_key in raw_name:
                                    sample_names.append(exp_name)
                                    matched = True
                                    break
                            if not matched:
                                sample_names.append(raw_name)  # Fallback to original name
                                
                except Exception as e:
                    print(f"Warning: Could not map to experiment names: {e}")
                    sample_names = raw_sample_names
            
            # Use experiment names directly
            simple_names = sample_names
            
            # Extract group information from experiment names if available
            groups = None
            group_colors = None
            if 'Experiment' in valid_rows.columns:
                # Try to extract group from experiment name (e.g., "DIA_10507_19" -> "DIA_10507")
                # Handle Experiment column - convert to string and handle NaN values
                experiment_col = valid_rows['Experiment'].astype(str).fillna('Unknown')
                groups = []
                for exp in experiment_col:
                    if pd.isna(exp) or exp == 'nan' or exp == 'Unknown' or not isinstance(exp, str):
                        groups.append(['Unknown'])
                    elif '_' in exp:
                        groups.append(exp.split('_')[:-1])
                    else:
                        groups.append([exp])
                groups = ['_'.join(g) if len(g) > 1 else g[0] for g in groups]
                unique_groups = list(set(groups))
                # Define group colors if not available
                default_group_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
                group_colors = {group: default_group_colors[i % len(default_group_colors)] 
                              for i, group in enumerate(unique_groups)}
            
            # Color by groups if available
            colors = [group_colors[groups[i]] if groups and group_colors else self.colors['primary'] 
                     for i in range(len(sample_names))]
            
            # 1. MS/MS Scans by Sample (spectra count)
            if 'MS/MS' in valid_rows.columns:
                msms_data = valid_rows['MS/MS'].tolist()
            else:
                msms_data = [0] * len(sample_names)
            fig.add_trace(go.Bar(x=simple_names, y=msms_data, name='MS/MS Scans', 
                               marker_color=colors), row=1, col=1)
            
            # 2. MS/MS Identified (PSM count)
            if 'MS/MS Identified' in valid_rows.columns:
                psm_data = valid_rows['MS/MS Identified'].tolist()
            else:
                psm_data = [0] * len(sample_names)
            fig.add_trace(go.Bar(x=simple_names, y=psm_data, name='PSM Count', 
                               marker_color=colors), row=1, col=2)
            
            # 3. Identification Rate (%)
            if 'MS/MS Identified [%]' in valid_rows.columns:
                id_rate_data = valid_rows['MS/MS Identified [%]'].tolist()
            else:
                id_rate_data = [0] * len(sample_names)
            fig.add_trace(go.Bar(x=simple_names, y=id_rate_data, name='ID Rate (%)', 
                               marker_color=colors), row=2, col=1)
            
            # 4. Peptide Sequences Count
            if 'Peptide Sequences Identified' in valid_rows.columns:
                peptide_data = valid_rows['Peptide Sequences Identified'].tolist()
            else:
                peptide_data = [0] * len(sample_names)
            fig.add_trace(go.Bar(x=simple_names, y=peptide_data, name='Peptide Sequences', 
                               marker_color=colors), row=2, col=2)
            
            # 5. Distribution of MS/MS Scans (Violin plot)
            fig.add_trace(go.Violin(y=msms_data, name='MS/MS Distribution', 
                                  box_visible=True, meanline_visible=True,
                                  points='all', pointpos=0, jitter=0.3,
                                  marker_color=colors[0]), row=3, col=1)
            
            # 6. Distribution of PSM Counts (Violin plot)
            fig.add_trace(go.Violin(y=psm_data, name='PSM Distribution', 
                                  box_visible=True, meanline_visible=True,
                                  points='all', pointpos=0, jitter=0.3,
                                  marker_color=colors[0]), row=3, col=2)
            
            # 7. Distribution of Identification Rates (Violin plot)
            fig.add_trace(go.Violin(y=id_rate_data, name='ID Rate Distribution', 
                                  box_visible=True, meanline_visible=True,
                                  points='all', pointpos=0, jitter=0.3,
                                  marker_color=colors[0]), row=4, col=1)
            
            # 8. Distribution of Peptide Counts (Violin plot)
            fig.add_trace(go.Violin(y=peptide_data, name='Peptide Distribution', 
                                  box_visible=True, meanline_visible=True,
                                  points='all', pointpos=0, jitter=0.3,
                                  marker_color=colors[0]), row=4, col=2)
        
        else:
            # Original format: sample columns
            sample_cols = [col for col in summary_df.columns if not col.startswith('Unnamed')]
            sample_names = sample_cols[:10]
            
            # Extract group information
            groups = None
            group_colors = None
            if metadata is not None:
                try:
                    sample_col = 'source name' if 'source name' in metadata.columns else metadata.columns[0]
                    if 'factor value[condition]' in metadata.columns:
                        group_col = 'factor value[condition]'
                    elif 'characteristics[condition]' in metadata.columns:
                        group_col = 'characteristics[condition]'
                    else:
                        group_col = None

                    if group_col:
                        sample_groups = dict(zip(metadata[sample_col], metadata[group_col]))
                        groups = [sample_groups.get(sample, 'Unknown') for sample in sample_names]
                        unique_groups = list(set(groups))
                        group_colors = self._assign_group_colors(unique_groups)
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

        fig.update_layout(height=1200, showlegend=False, title_text="ID Summary Report")
        return fig
