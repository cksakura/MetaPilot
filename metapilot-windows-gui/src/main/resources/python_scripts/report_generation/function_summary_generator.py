import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from .base_generator import BaseReportGenerator

class FunctionSummaryGenerator(BaseReportGenerator):
    """Generator for function summary reports."""
    
    def create_function_summary(self, metadata=None, input_file_path=None):
        """Create function summary with data-driven functional category analysis."""
        # All data comes from separate files in the input_file_path folder
        # We'll get sample info from the first separate file we find
        
        # Check for available functional categories and files
        available_plots = self._check_available_functional_plots(input_file_path)
        
        if not available_plots:
            return self.create_empty_data_plot(
                "Function Summary Report",
                "No functional categories detected in the data.<br>Please check if the functional annotation completed successfully."
            )
        
        # Get sample info from the first available separate file
        sample_df = available_plots[0]['data_source']
        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(sample_df, metadata)
        
        # Use experiment names directly (no file ID mapping needed)
        simple_names = sample_names
        
        # Create dynamic subplot layout based on available plots
        num_rows = len(available_plots)
        subplot_titles = []
        specs = []
        
        for plot_info in available_plots:
            plot_name = plot_info['name']
            plot_type = plot_info['type']
            subplot_titles.append(f"{plot_name} {plot_type}")
            
            # Set appropriate subplot type based on plot type
            if plot_type == "Stacked Bar Chart":
                specs.append([{"type": "bar"}])
            elif plot_type == "PCA":
                specs.append([{"type": "xy"}])
            elif plot_type == "Heatmap":
                specs.append([{"type": "xy"}])  # Use xy for heatmap as it's more flexible
            else:
                specs.append([{"type": "bar"}])  # Default fallback
        
        # Calculate vertical spacing based on number of rows
        # Maximum spacing is 1/(rows-1), use 25% of that to be safe
        if num_rows > 1:
            max_spacing = 1.0 / (num_rows - 1)
            vertical_spacing = min(0.1, max_spacing * 0.25)
        else:
            vertical_spacing = 0.1
        
        fig = make_subplots(
            rows=num_rows, cols=1,
            subplot_titles=subplot_titles,
            specs=specs,
            vertical_spacing=vertical_spacing
        )
        
        # Pre-build data matrices for PCA and heatmap plots to avoid rebuilding
        for plot_info in available_plots:
            if plot_info['type'] in ['PCA', 'Heatmap']:
                data_matrix_info = self._build_data_matrix(
                    plot_info['data_source'], 
                    plot_info['column'], 
                    intensity_cols
                )
                plot_info['data_matrix_info'] = data_matrix_info
        
        # Create visualizations for each available plot
        for row_idx, plot_info in enumerate(available_plots, 1):
            plot_name = plot_info['name']
            plot_type = plot_info['type']
            data_source = plot_info['data_source']
            column = plot_info['column']
            
            if plot_type == "Stacked Bar Chart":
                self._create_stacked_bar_chart(fig, data_source, column, intensity_cols, simple_names, row_idx, plot_name)
            elif plot_type == "PCA":
                self._create_pca_plot(fig, plot_info.get('data_matrix_info'), simple_names, row_idx, plot_name)
            elif plot_type == "Heatmap":
                self._create_heatmap(fig, plot_info.get('data_matrix_info'), simple_names, row_idx, plot_name)
        
        # Update layout with fixed height per plot
        height_per_plot = 400
        total_height = num_rows * height_per_plot
        
        # Configure legend to prevent overlap with subplots
        fig.update_layout(
            height=total_height, 
            showlegend=True, 
            title_text="Function Summary Report", 
            barmode="stack",
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1.0,
                xanchor="left", 
                x=1.02,
                bgcolor="rgba(255,255,255,0.8)",
                bordercolor="rgba(0,0,0,0.2)",
                borderwidth=1
            )
        )
        
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
    
    def _check_available_functional_plots(self, input_file_path=None):
        """
        Check for available functional categories and files, return list of plots to create.
        
        Returns:
            list: List of plot information dictionaries
        """
        import os
        
        plots = []
        
        # For function reports, input_file_path is now the folder path containing all separate files
        # Get the directory from the input_file_path parameter (which is now the folder path)
        func_dir = input_file_path if input_file_path else ''
        
        print(f"DEBUG: input_file_path (folder): '{input_file_path}'")
        print(f"DEBUG: func_dir: '{func_dir}'")
        
        # All function data now comes from separate files (COG Category, NOG Category, COG Accession, NOG Accession, GO, EC, KEGG, BRITE)
        # The main functions_report.tsv file is only used for validation and basic info
        
        # Check for separate files (COG, NOG, Gene Ontology, EC, KEGG, BRITE)
        print(f"DEBUG: Checking for separate files. func_dir exists: {bool(func_dir)}")
        if func_dir:
            separate_files = {
                'COG Category': 'COG_category_report.tsv',
                'NOG Category': 'NOG_category_report.tsv',
                'COG Accession': 'COG_accession_report.tsv',
                'NOG Accession': 'NOG_accession_report.tsv',
                'Gene Ontology': 'Gene_Ontology_id_report.tsv',
                'EC': 'EC_id_report.tsv',
                'KEGG': 'KEGG_ko_report.tsv',
                'BRITE': 'BRITE_report.tsv'
            }
            
            for category_name, filename in separate_files.items():
                file_path = os.path.join(func_dir, filename)
                print(f"DEBUG: Looking for {category_name} file at: {file_path}")
                print(f"DEBUG:   File exists: {os.path.exists(file_path)}")
                if os.path.exists(file_path):
                    print(f"DEBUG:   File size: {os.path.getsize(file_path)} bytes")
                if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                    try:
                        print(f"DEBUG:   Reading {category_name} file...")
                        separate_df = pd.read_csv(file_path, sep='\t')
                        print(f"DEBUG:   Loaded {len(separate_df)} rows, {len(separate_df.columns)} columns")
                        if len(separate_df) > 0 and len(separate_df.columns) > 1:
                            # Ensure intensity columns are numeric (same as main file)
                            separate_intensity_cols = self._detect_intensity_columns(separate_df)
                            for col in separate_intensity_cols:
                                if col in separate_df.columns:
                                    separate_df[col] = pd.to_numeric(separate_df[col], errors='coerce').fillna(0)
                            
                            # Determine plot type based on category name
                            if "Category" in category_name:
                                plot_type = "Stacked Bar Chart"
                                print(f"DEBUG:   Adding {plot_type} plot for {category_name}")
                                plots.append({
                                    'name': category_name,
                                    'type': plot_type,
                                    'data_source': separate_df,
                                    'column': separate_df.columns[0]
                                })
                            else:
                                plot_type = "PCA"
                                print(f"DEBUG:   Adding {plot_type} and Heatmap plots for {category_name}")
                                plots.append({
                                    'name': category_name,
                                    'type': plot_type,
                                    'data_source': separate_df,
                                    'column': separate_df.columns[0]
                                })
                                plots.append({
                                    'name': category_name,
                                    'type': 'Heatmap',
                                    'data_source': separate_df,
                                    'column': separate_df.columns[0]
                                })
                        else:
                            print(f"DEBUG:   Skipping {category_name} - insufficient data")
                    except Exception as e:
                        print(f"DEBUG:   Error reading {category_name} file: {e}")
                        continue
        else:
            print(f"DEBUG: No func_dir specified, cannot check for separate files")
        
        print(f"DEBUG: Total plots to create: {len(plots)}")
        for i, plot in enumerate(plots):
            print(f"DEBUG:   Plot {i+1}: {plot['name']} - {plot['type']}")
        
        return plots
    
    def _build_data_matrix(self, data_source, category_col, intensity_cols):
        """Build abundance matrix once for reuse in PCA and heatmap.
        Returns a dictionary with matrix, labels, and metadata."""
        
        try:
            # Filter data
            category_data = data_source[data_source[category_col].notna() & (data_source[category_col] != "") & (data_source[category_col] != "-")]
            
            if len(category_data) == 0 or category_col not in data_source.columns:
                return None
            
            # Get intensity columns that exist in this data source
            data_source_intensity_cols = self._detect_intensity_columns(data_source)
            
            # Prepare data matrix: rows = functional IDs, columns = samples
            functional_ids = category_data[category_col].unique()
            functional_ids = [fid for fid in functional_ids if pd.notna(fid) and fid != ""]
            
            if len(functional_ids) == 0:
                return None
            
            # Create abundance matrix
            abundance_matrix = []
            for fid in functional_ids:
                fid_data = category_data[category_data[category_col] == fid]
                abundances = []
                for col in data_source_intensity_cols:
                    # Handle NaN values
                    abundance_values = fid_data[col].fillna(0)
                    total_abundance = abundance_values.sum()
                    abundances.append(total_abundance if not pd.isna(total_abundance) else 0)
                abundance_matrix.append(abundances)
            
            abundance_matrix = np.array(abundance_matrix)
            
            # Filter out zero/near-zero variance rows
            if abundance_matrix.shape[0] > 0:
                row_variances = np.var(abundance_matrix, axis=1)
                non_zero_var_mask = row_variances > 0
                abundance_matrix = abundance_matrix[non_zero_var_mask]
                functional_ids = [functional_ids[i] for i in range(len(functional_ids)) if non_zero_var_mask[i]]
            
            # Check if we have enough data
            if abundance_matrix.shape[0] < 2 or abundance_matrix.shape[1] < 2:
                return None
            
            return {
                'matrix': abundance_matrix,
                'labels': functional_ids,
                'shape': abundance_matrix.shape
            }
            
        except Exception as e:
            print(f"DEBUG: Error building data matrix: {e}")
            return None
    
    def _create_stacked_bar_chart(self, fig, data_source, category_col, intensity_cols, simple_names, row, plot_name):
        """Create a stacked bar chart for COG/NOG categories with relative abundance."""
        try:
            
            # Get data for this category
            category_data = data_source[data_source[category_col].notna() & (data_source[category_col] != "") & (data_source[category_col] != "-")]
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
    
    def _create_pca_plot(self, fig, data_matrix_info, simple_names, row, plot_name):
        """Create PCA plot using pre-built data matrix."""
        try:
            import numpy as np
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            
            # Check if we have valid data
            if data_matrix_info is None:
                fig.add_annotation(text=f"Insufficient data for PCA", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            abundance_matrix = data_matrix_info['matrix']
            functional_ids = data_matrix_info['labels']
            
            # Standardize the data
            scaler = StandardScaler()
            abundance_matrix_scaled = scaler.fit_transform(abundance_matrix)
            
            # Transpose for sample-level PCA (samples as points, functional IDs as features)
            abundance_matrix_scaled = abundance_matrix_scaled.T
            pca = PCA(n_components=min(2, abundance_matrix_scaled.shape[1]))
            pca_result = pca.fit_transform(abundance_matrix_scaled)
            point_labels = simple_names  # Samples as points
            
            print(f"DEBUG PCA: {plot_name}")
            print(f"  PCA result shape: {pca_result.shape}")
            print(f"  PC1 range: [{pca_result[:, 0].min():.2f}, {pca_result[:, 0].max():.2f}]")
            print(f"  PC2 range: [{pca_result[:, 1].min():.2f}, {pca_result[:, 1].max():.2f}]")
            print(f"  Number of points: {len(point_labels)}")
            print(f"  Sample names: {point_labels}")
            print(f"  Adding trace to subplot row {row}, col 1")
            
            # Convert to lists explicitly
            x_coords = pca_result[:, 0].tolist()
            y_coords = pca_result[:, 1].tolist()
            print(f"  X coordinates: {x_coords}")
            print(f"  Y coordinates: {y_coords}")
            
            # Create scatter plot
            fig.add_trace(
                go.Scatter(
                    x=x_coords,
                    y=y_coords,
                    mode='markers+text',  # Add text mode to see if points exist
                    marker=dict(
                        size=12,  # Larger markers
                        color=self.colors["primary"],
                        opacity=0.8,
                        line=dict(width=2, color='#333333')  # Add border to markers
                    ),
                    text=point_labels,
                    textposition='top center',
                    hovertemplate=f'<b>%{{text}}</b><br>PC1: %{{x:.2f}}<br>PC2: %{{y:.2f}}<extra></extra>',
                    name=f"{plot_name} Sample PCA",
                    showlegend=False
                ),
                row=row, col=1
            )
            print(f"  Trace added successfully")
            
            # Update axes labels
            pc1_var = pca.explained_variance_ratio_[0] * 100
            pc2_var = pca.explained_variance_ratio_[1] * 100 if len(pca.explained_variance_ratio_) > 1 else 0
            
            fig.update_xaxes(title_text=f"PC1 ({pc1_var:.1f}%)", row=row, col=1)
            fig.update_yaxes(title_text=f"PC2 ({pc2_var:.1f}%)", row=row, col=1)
            
        except Exception as e:
            fig.add_annotation(text=f"{plot_name} PCA creation failed: {str(e)[:50]}...", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
    
    def _create_heatmap(self, fig, data_matrix_info, simple_names, row, plot_name):
        """Create heatmap using pre-built data matrix."""
        try:
            import numpy as np
            from scipy.cluster.hierarchy import linkage, dendrogram
            from scipy.spatial.distance import pdist
            
            # Check if we have valid data
            if data_matrix_info is None:
                fig.add_annotation(text=f"Insufficient data for heatmap", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
                return
            
            abundance_matrix = data_matrix_info['matrix'].copy()
            functional_ids = data_matrix_info['labels'].copy()
            
            # Limit the number of functional IDs for better visualization (top 50 by total abundance)
            if len(functional_ids) > 50:
                total_abundances = abundance_matrix.sum(axis=1)
                top_indices = np.argsort(total_abundances)[-50:]
                abundance_matrix = abundance_matrix[top_indices]
                functional_ids = [functional_ids[i] for i in top_indices]
            
            # Hierarchical clustering on rows (functional IDs)
            if abundance_matrix.shape[0] > 1 and abundance_matrix.size > 0:
                row_linkage = linkage(pdist(abundance_matrix), method='ward')
                row_order = dendrogram(row_linkage, no_plot=True)['leaves']
                abundance_matrix = abundance_matrix[row_order]
                functional_ids = [functional_ids[i] for i in row_order]
            
            # Hierarchical clustering on columns (samples)
            if abundance_matrix.shape[1] > 1:
                col_linkage = linkage(pdist(abundance_matrix.T), method='ward')
                col_order = dendrogram(col_linkage, no_plot=True)['leaves']
                abundance_matrix = abundance_matrix[:, col_order]
                sample_names_clustered = [simple_names[i] if i < len(simple_names) else f"Sample_{i+1}" for i in col_order]
            else:
                sample_names_clustered = simple_names
            
            # Create heatmap
            fig.add_trace(
                go.Heatmap(
                    z=abundance_matrix.tolist(),  # Convert to list like taxonomy
                    x=sample_names_clustered,
                    y=functional_ids,
                    colorscale='Viridis',
                    hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Abundance: %{z:.2f}<extra></extra>',
                    name=f"{plot_name} Heatmap",
                    showscale=False
                ),
                row=row, col=1
            )
            
            # Update axes labels
            fig.update_xaxes(title_text="Samples", row=row, col=1)
            fig.update_yaxes(title_text=f"{plot_name} IDs", row=row, col=1)
            
        except Exception as e:
            fig.add_annotation(text=f"{plot_name} heatmap creation failed: {str(e)[:50]}...", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False, row=row, col=1)
    
    
    