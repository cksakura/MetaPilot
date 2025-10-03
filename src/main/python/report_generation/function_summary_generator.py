import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from .base_generator import BaseReportGenerator

class FunctionSummaryGenerator(BaseReportGenerator):
    """Generator for function summary reports."""
    
    def create_function_summary(self, func_df, metadata=None):
        """Create function summary with data-driven functional category analysis."""
        # Validate data before processing
        if not self.validate_data(func_df, "function"):
            return self.create_empty_data_plot(
                "Function Summary Report",
                "No function data available for analysis.<br>Please check if the functional annotation completed successfully."
            )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(func_df, metadata)
        
        # Use experiment names directly (no file ID mapping needed)
        simple_names = sample_names
        
        # Ensure all intensity columns are numeric to prevent dtype errors
        for col in intensity_cols:
            if col in func_df.columns:
                func_df[col] = pd.to_numeric(func_df[col], errors='coerce').fillna(0)
        
        # Check for available functional categories and files
        available_plots = self._check_available_functional_plots(func_df)
        
        if not available_plots:
            return self.create_empty_data_plot(
                "Function Summary Report",
                "No functional categories detected in the data.<br>Please check if the functional annotation completed successfully."
            )
        
        # Create dynamic subplot layout based on available plots
        num_rows = len(available_plots)
        subplot_titles = []
        specs = []
        
        for plot_info in available_plots:
            plot_name = plot_info['name']
            plot_type = plot_info['type']
            subplot_titles.append(f"{plot_name} {plot_type}")
            specs.append([{"type": "bar"}])
        
        fig = make_subplots(
            rows=num_rows, cols=1,
            subplot_titles=subplot_titles,
            specs=specs,
            vertical_spacing=0.1
        )
        
        # Create visualizations for each available plot
        for row_idx, plot_info in enumerate(available_plots, 1):
            plot_name = plot_info['name']
            plot_type = plot_info['type']
            data_source = plot_info['data_source']
            column = plot_info['column']
            
            if plot_type == "Stacked Bar Chart":
                self._create_stacked_bar_chart(fig, data_source, column, intensity_cols, simple_names, row_idx, plot_name)
            elif plot_type == "PCA":
                self._create_pca_plot(fig, data_source, column, intensity_cols, simple_names, row_idx, plot_name)
            elif plot_type == "Heatmap":
                self._create_heatmap(fig, data_source, column, intensity_cols, simple_names, row_idx, plot_name)
        
        # Update layout with fixed height per plot
        height_per_plot = 400
        total_height = num_rows * height_per_plot
        fig.update_layout(height=total_height, showlegend=True, title_text="Function Summary Report")
        
        return fig
    
    def _detect_functional_categories(self, func_df):
        """Detect available functional categories from column names."""
        available_categories = []
        
        # Only detect COG and NOG categories
        category_patterns = {
            'COG Category': ['COG category', 'COG_Category', 'COG_category', 'cog_category', 'COG'],
            'NOG Category': ['NOG category', 'NOG_Category', 'NOG_category', 'nog_category', 'NOG']
        }
        
        # Check for each category pattern
        for category_name, patterns in category_patterns.items():
            for pattern in patterns:
                if pattern in func_df.columns:
                    available_categories.append({
                        'name': category_name,
                        'column': pattern
                    })
                    break  # Found the category, move to next
        
        return available_categories
    
    def _check_available_functional_plots(self, func_df):
        """
        Check for available functional categories and files, return list of plots to create.
        
        Returns:
            list: List of plot information dictionaries
        """
        import os
        
        plots = []
        
        # Get the directory of the current file
        current_file_path = func_df.attrs.get('file_path', '') if hasattr(func_df, 'attrs') else ''
        func_dir = os.path.dirname(current_file_path) if current_file_path else ''
        
        # Check for stacked bar charts (COG category, NOG category from main file)
        if 'COG category' in func_df.columns:
            cog_data = func_df[func_df['COG category'].notna() & (func_df['COG category'] != "")]
            if len(cog_data) > 0:
                plots.append({
                    'name': 'COG Category',
                    'type': 'Stacked Bar Chart',
                    'data_source': func_df,
                    'column': 'COG category'
                })
        
        if 'NOG category' in func_df.columns:
            nog_data = func_df[func_df['NOG category'].notna() & (func_df['NOG category'] != "")]
            if len(nog_data) > 0:
                plots.append({
                    'name': 'NOG Category',
                    'type': 'Stacked Bar Chart',
                    'data_source': func_df,
                    'column': 'NOG category'
                })
        
        # Check for PCA + Heatmap from main file (COG accession, NOG accession)
        if 'COG accession' in func_df.columns:
            cog_acc_data = func_df[func_df['COG accession'].notna() & (func_df['COG accession'] != "")]
            if len(cog_acc_data) > 0:
                plots.append({
                    'name': 'COG Accession',
                    'type': 'PCA',
                    'data_source': func_df,
                    'column': 'COG accession'
                })
                plots.append({
                    'name': 'COG Accession',
                    'type': 'Heatmap',
                    'data_source': func_df,
                    'column': 'COG accession'
                })
        
        if 'NOG accession' in func_df.columns:
            nog_acc_data = func_df[func_df['NOG accession'].notna() & (func_df['NOG accession'] != "")]
            if len(nog_acc_data) > 0:
                plots.append({
                    'name': 'NOG Accession',
                    'type': 'PCA',
                    'data_source': func_df,
                    'column': 'NOG accession'
                })
                plots.append({
                    'name': 'NOG Accession',
                    'type': 'Heatmap',
                    'data_source': func_df,
                    'column': 'NOG accession'
                })
        
        # Check for separate files (Gene Ontology, EC, KEGG, BRITE)
        if func_dir:
            separate_files = {
                'Gene Ontology': 'Gene_Ontology_id_report.tsv',
                'EC': 'EC_id_report.tsv',
                'KEGG': 'KEGG_ko_report.tsv',
                'BRITE': 'BRITE_report.tsv'
            }
            
            for category_name, filename in separate_files.items():
                file_path = os.path.join(func_dir, filename)
                if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                    try:
                        separate_df = pd.read_csv(file_path, sep='\t')
                        if len(separate_df) > 0 and len(separate_df.columns) > 1:
                            plots.append({
                                'name': category_name,
                                'type': 'PCA',
                                'data_source': separate_df,
                                'column': separate_df.columns[0]
                            })
                            plots.append({
                                'name': category_name,
                                'type': 'Heatmap',
                                'data_source': separate_df,
                                'column': separate_df.columns[0]
                            })
                    except Exception:
                        continue
        
        return plots
    
    def _create_stacked_bar_chart(self, fig, data_source, category_col, intensity_cols, simple_names, row, plot_name):
        """Create a stacked bar chart for COG/NOG categories with relative abundance."""
        try:
            
            # Get data for this category
            category_data = data_source[data_source[category_col].notna() & (data_source[category_col] != "")]
            if len(category_data) == 0:
                fig.add_annotation(text=f"No {plot_name} data available", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Debug: Check if category_col exists
            if category_col not in data_source.columns:
                fig.add_annotation(text=f"Column '{category_col}' not found in data", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Calculate relative abundance for each sample
            sample_compositions = {}
            all_categories = set()
            
            for i, col in enumerate(intensity_cols):
                sample_name = simple_names[i] if i < len(simple_names) else f"Sample_{i+1}"
                
                # Get abundance for each category in this sample
                sample_abundances = {}
                for category in category_data[category_col].unique():
                    if pd.notna(category) and category != "":
                        category_mask = category_data[category_col] == category
                        # Handle NaN values in intensity column
                        abundance_values = category_data[category_mask][col].fillna(0)
                        abundance = abundance_values.sum()
                        if abundance > 0 and not pd.isna(abundance):
                            sample_abundances[category] = abundance
                            all_categories.add(category)
                
                # Calculate relative abundance
                total_abundance = sum(sample_abundances.values())
                if total_abundance > 0 and not pd.isna(total_abundance):
                    sample_compositions[sample_name] = {
                        category: (abundance / total_abundance) * 100 
                        for category, abundance in sample_abundances.items()
                    }
            
            if not sample_compositions:
                fig.add_annotation(text=f"No {plot_name} abundance data available", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Create stacked bar chart
            sample_names_list = list(sample_compositions.keys())
            all_categories = sorted(list(all_categories))
            
            # Color palette for categories
            color_palette = [
                self.colors["primary"], self.colors["secondary"], self.colors["tertiary"], 
                self.colors["success"], self.colors["warning"], self.colors["danger"], 
                self.colors["light_blue"], self.colors["light_purple"], "#FF6B6B", "#4ECDC4",
                "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD", "#98D8C8"
            ]
            
            # Create stacked bar chart - one trace per category, all samples in each trace
            for i, category in enumerate(all_categories):
                abundances = []
                for sample_name in sample_names_list:
                    abundance = sample_compositions[sample_name].get(category, 0)
                    # Ensure abundance is a valid number
                    if pd.isna(abundance):
                        abundance = 0
                    abundances.append(abundance)
                
                fig.add_trace(
                    go.Bar(
                        x=sample_names_list,
                        y=abundances,
                        name=category,
                        marker_color=color_palette[i % len(color_palette)],
                        showlegend=True,
                        hovertemplate=f'<b>{category}</b><br>Sample: %{{x}}<br>Abundance: %{{y:.1f}}%<extra></extra>'
                    ),
                    row=row, col=1
                )
            
            # Update axes labels
            fig.update_xaxes(title_text="Samples", row=row, col=1)
            fig.update_yaxes(title_text="Relative Abundance (%)", row=row, col=1)
            
        except Exception as e:
            fig.add_annotation(text=f"{plot_name} chart creation failed: {str(e)[:30]}...", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
    
    def _create_pca_plot(self, fig, data_source, category_col, intensity_cols, simple_names, row, plot_name):
        """Create PCA plot for functional categories."""
        try:
            import numpy as np
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            
            # Get data for this category
            category_data = data_source[data_source[category_col].notna() & (data_source[category_col] != "")]
            if len(category_data) == 0:
                fig.add_annotation(text=f"No {plot_name} data available", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Check if category_col exists
            if category_col not in data_source.columns:
                fig.add_annotation(text=f"Column '{category_col}' not found in data", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Prepare data matrix: rows = functional IDs, columns = samples
            functional_ids = category_data[category_col].unique()
            functional_ids = [fid for fid in functional_ids if pd.notna(fid) and fid != ""]
            
            if len(functional_ids) == 0:
                fig.add_annotation(text=f"No valid {plot_name} IDs found", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Create abundance matrix
            abundance_matrix = []
            for fid in functional_ids:
                fid_data = category_data[category_data[category_col] == fid]
                abundances = []
                for col in intensity_cols:
                    # Handle NaN values
                    abundance_values = fid_data[col].fillna(0)
                    total_abundance = abundance_values.sum()
                    abundances.append(total_abundance if not pd.isna(total_abundance) else 0)
                abundance_matrix.append(abundances)
            
            abundance_matrix = np.array(abundance_matrix)
            
            # Check if we have enough data for PCA (at least 2 samples and 2 features)
            if abundance_matrix.shape[0] < 2 or abundance_matrix.shape[1] < 2:
                fig.add_annotation(text=f"Insufficient data for PCA: {abundance_matrix.shape[0]} features, {abundance_matrix.shape[1]} samples", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Standardize the data
            scaler = StandardScaler()
            abundance_matrix_scaled = scaler.fit_transform(abundance_matrix)
            
            # Perform PCA
            pca = PCA(n_components=min(2, abundance_matrix.shape[1]))
            pca_result = pca.fit_transform(abundance_matrix_scaled)
            
            # Create scatter plot
            fig.add_trace(
                go.Scatter(
                    x=pca_result[:, 0],
                    y=pca_result[:, 1],
                    mode='markers',
                    marker=dict(
                        size=8,
                        color=self.colors["primary"],
                        opacity=0.7
                    ),
                    text=functional_ids,
                    hovertemplate=f'<b>%{{text}}</b><br>PC1: %{{x:.2f}}<br>PC2: %{{y:.2f}}<extra></extra>',
                    name=f"{plot_name} IDs",
                    showlegend=False
                ),
                row=row, col=1
            )
            
            # Update axes labels
            pc1_var = pca.explained_variance_ratio_[0] * 100
            pc2_var = pca.explained_variance_ratio_[1] * 100 if len(pca.explained_variance_ratio_) > 1 else 0
            
            fig.update_xaxes(title_text=f"PC1 ({pc1_var:.1f}%)", row=row, col=1)
            fig.update_yaxes(title_text=f"PC2 ({pc2_var:.1f}%)", row=row, col=1)
            
        except Exception as e:
            fig.add_annotation(text=f"{plot_name} PCA creation failed: {str(e)[:50]}...", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
    
    def _create_heatmap(self, fig, data_source, category_col, intensity_cols, simple_names, row, plot_name):
        """Create heatmap for functional categories."""
        try:
            import numpy as np
            from scipy.cluster.hierarchy import linkage, dendrogram
            from scipy.spatial.distance import pdist
            
            # Get data for this category
            category_data = data_source[data_source[category_col].notna() & (data_source[category_col] != "")]
            if len(category_data) == 0:
                fig.add_annotation(text=f"No {plot_name} data available", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Check if category_col exists
            if category_col not in data_source.columns:
                fig.add_annotation(text=f"Column '{category_col}' not found in data", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Prepare data matrix: rows = functional IDs, columns = samples
            functional_ids = category_data[category_col].unique()
            functional_ids = [fid for fid in functional_ids if pd.notna(fid) and fid != ""]
            
            if len(functional_ids) == 0:
                fig.add_annotation(text=f"No valid {plot_name} IDs found", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Create abundance matrix
            abundance_matrix = []
            for fid in functional_ids:
                fid_data = category_data[category_data[category_col] == fid]
                abundances = []
                for col in intensity_cols:
                    # Handle NaN values
                    abundance_values = fid_data[col].fillna(0)
                    total_abundance = abundance_values.sum()
                    abundances.append(total_abundance if not pd.isna(total_abundance) else 0)
                abundance_matrix.append(abundances)
            
            abundance_matrix = np.array(abundance_matrix)
            
            # Check if we have enough data for heatmap
            if abundance_matrix.shape[0] < 2 or abundance_matrix.shape[1] < 2:
                fig.add_annotation(text=f"Insufficient data for heatmap: {abundance_matrix.shape[0]} features, {abundance_matrix.shape[1]} samples", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            # Limit the number of functional IDs for better visualization (top 50 by total abundance)
            if len(functional_ids) > 50:
                total_abundances = abundance_matrix.sum(axis=1)
                top_indices = np.argsort(total_abundances)[-50:]
                abundance_matrix = abundance_matrix[top_indices]
                functional_ids = [functional_ids[i] for i in top_indices]
            
            # Create heatmap
            fig.add_trace(
                go.Heatmap(
                    z=abundance_matrix,
                    x=simple_names,
                    y=functional_ids,
                    colorscale='Viridis',
                    hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Abundance: %{z:.2f}<extra></extra>',
                    name=f"{plot_name} Heatmap",
                    showscale=True
                ),
                row=row, col=1
            )
            
            # Update axes labels
            fig.update_xaxes(title_text="Samples", row=row, col=1)
            fig.update_yaxes(title_text=f"{plot_name} IDs", row=row, col=1)
            
        except Exception as e:
            fig.add_annotation(text=f"{plot_name} heatmap creation failed: {str(e)[:50]}...", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
    
    
    