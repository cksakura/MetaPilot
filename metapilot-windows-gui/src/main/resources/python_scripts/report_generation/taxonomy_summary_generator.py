import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from .base_generator import BaseReportGenerator

class TaxonomySummaryGenerator(BaseReportGenerator):
    """Generator for taxonomy summary reports."""
    
    def create_taxonomy_summary(self, taxa_df, metadata=None):
        """Create taxonomy summary visualization with individual plots for each taxonomic level"""
        # Validate data before processing
        if not self.validate_data(taxa_df, "taxonomy"):
            return self.create_empty_data_plot(
                "Taxonomy Summary Report",
                "No taxonomy data available for analysis.<br>Please check if the taxonomic analysis completed successfully."
            )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(taxa_df, metadata)
        simple_names = sample_names

        for c in intensity_cols:
            if c in taxa_df.columns:
                taxa_df[c] = pd.to_numeric(taxa_df[c], errors="coerce").fillna(0)

        taxonomic_ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
        available_ranks = [
            rank for rank in taxonomic_ranks
            if rank in taxa_df.columns
            and len(taxa_df[taxa_df[rank].notna() & (taxa_df[rank] != "")]) > 0
        ]
        
        if not available_ranks:
            return self.create_empty_data_plot("Taxonomy Summary Report", "No taxonomic data available")
        
        total_rows = len(available_ranks) * 4
        subplot_titles = []
        specs = []
        
        for rank in available_ranks:
            subplot_titles.extend([
                f"{rank.capitalize()} Intensity Distribution (Log Intensity)",
                f"{rank.capitalize()} PCA (Taxa in Sample Space)",
                f"{rank.capitalize()} Log Intensity Heatmap (Clustered)",
                f"{rank.capitalize()} Composition by Sample"
            ])
            specs.extend([
                [{"type": "violin"}],
                [{"type": "scatter"}],
                [{"type": "heatmap"}],
                [{"type": "bar"}]
            ])
        
        subplot_spacing = min(0.06, 0.25 / max(total_rows - 1, 1)) if total_rows > 1 else 0
        fig = make_subplots(
            rows=total_rows,
            cols=1,
            subplot_titles=subplot_titles,
            specs=specs,
            vertical_spacing=subplot_spacing
        )
        
        current_row = 1
        
        for rank in available_ranks:
            rank_data = taxa_df[taxa_df[rank].notna() & (taxa_df[rank] != "")]
            if rank_data.empty:
                for _ in range(4):
                    fig.add_annotation(
                        text="No data",
                        x=0.5,
                        y=0.5,
                        showarrow=False,
                        row=current_row,
                        col=1
                    )
                    current_row += 1
                continue

            rank_summary = rank_data.groupby(rank)[intensity_cols].sum().reset_index()
            taxon_labels = rank_summary[rank].tolist()

            if len(intensity_cols) > 0:
                intensity_matrix = rank_summary[intensity_cols].to_numpy(dtype=float)
                log_intensity_matrix = np.where(intensity_matrix > 0, np.log10(intensity_matrix), 0.0)
            else:
                intensity_matrix = np.empty((len(rank_summary), 0))
                log_intensity_matrix = intensity_matrix

            # 1. Violin plot
            sample_data = {}
            for i, col in enumerate(intensity_cols):
                sample_name = simple_names[i] if i < len(simple_names) else f"Sample_{i+1}"
                sample_intensities = []
                for val in rank_summary[col]:
                    if pd.notna(val) and val > 0:
                        sample_intensities.append(np.log10(val))
                if sample_intensities:
                    sample_data[sample_name] = sample_intensities

            if sample_data:
                for sample_name, intensities in sample_data.items():
                    fig.add_trace(
                        go.Violin(
                            y=intensities,
                            x=[sample_name] * len(intensities),
                            name=sample_name,
                            box_visible=True,
                            meanline_visible=True,
                            points='all',
                            pointpos=0,
                            jitter=0.3,
                            marker_color=self.colors["primary"],
                            showlegend=False
                        ),
                        row=current_row,
                        col=1
                    )
            else:
                fig.add_annotation(
                    text="No intensity data",
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    row=current_row,
                    col=1
                )

            fig.update_xaxes(title_text="Samples", row=current_row, col=1)
            fig.update_yaxes(title_text="Log Intensity", row=current_row, col=1)
            current_row += 1

            # 2. PCA scatter (with group classification if metadata available)
            if log_intensity_matrix.shape[0] >= 2 and log_intensity_matrix.shape[1] >= 2:
                try:
                    # Filter out zero/near-zero variance rows (all-zero taxa)
                    row_variances = np.var(log_intensity_matrix, axis=1)
                    non_zero_var_mask = row_variances > 0
                    filtered_matrix = log_intensity_matrix[non_zero_var_mask]
                    filtered_taxon_labels = [taxon_labels[i] for i in range(len(taxon_labels)) if non_zero_var_mask[i]]
                    
                    # Filter out zero-variance sample columns
                    if filtered_matrix.shape[0] > 0:
                        col_variances = np.var(filtered_matrix, axis=0)
                        non_zero_col_mask = col_variances > 0
                        filtered_matrix = filtered_matrix[:, non_zero_col_mask]
                    
                    # Check if we still have enough data after filtering
                    if filtered_matrix.shape[0] < 2 or filtered_matrix.shape[1] < 2:
                        raise ValueError("Insufficient data after filtering zero variance")
                    
                    scaled = StandardScaler().fit_transform(filtered_matrix)
                    
                    # Transpose for sample-level PCA (samples as points, taxa as features)
                    scaled = scaled.T
                    pca = PCA(n_components=min(2, scaled.shape[1]))
                    components = pca.fit_transform(scaled)
                    point_labels = simple_names  # Samples as points
                    explained = pca.explained_variance_ratio_

                    # Convert to lists explicitly for better Plotly compatibility
                    x_coords = components[:, 0].tolist()
                    y_coords = components[:, 1].tolist()
                    
                    # Plot with group classification if available
                    if groups is not None and group_colors is not None:
                        # Group samples by their classification
                        unique_groups = list(dict.fromkeys(groups))
                        for group_name in unique_groups:
                            group_indices = [i for i, g in enumerate(groups) if g == group_name]
                            group_x = [x_coords[i] for i in group_indices]
                            group_y = [y_coords[i] for i in group_indices]
                            group_labels = [point_labels[i] for i in group_indices]
                            group_color = group_colors.get(group_name, self.colors['secondary'])
                            
                            # Add scatter points for this group
                            fig.add_trace(
                                go.Scatter(
                                    x=group_x,
                                    y=group_y,
                                    mode='markers+text',
                                    marker=dict(
                                        color=group_color,
                                        size=12,
                                        opacity=0.8,
                                        line=dict(width=2, color='#333333')
                                    ),
                                    text=group_labels,
                                    textposition='top center',
                                    name=group_name,
                                    legendgroup=group_name,
                                    hovertemplate=f'<b>%{{text}}</b><br>Group: {group_name}<br>PC1: %{{x:.2f}}<br>PC2: %{{y:.2f}}<extra></extra>',
                                    showlegend=False
                                ),
                                row=current_row,
                                col=1
                            )
                            
                            # Add confidence ellipse if group has 3+ samples
                            if len(group_x) >= 3:
                                self._add_confidence_ellipse(fig, group_x, group_y, group_color, group_name, current_row)
                    else:
                        # No grouping - single color for all points
                        fig.add_trace(
                            go.Scatter(
                                x=x_coords,
                                y=y_coords,
                                mode='markers+text',
                                marker=dict(
                                    color=self.colors['secondary'], 
                                    size=12,
                                    opacity=0.8,
                                    line=dict(width=2, color='#333333')
                                ),
                                text=point_labels,
                                textposition='top center',
                                hovertemplate='<b>%{text}</b><br>PC1: %{x:.2f}<br>PC2: %{y:.2f}<br><extra></extra>',
                                showlegend=False
                            ),
                            row=current_row,
                            col=1
                        )

                    pc1_label = f"PC1 ({explained[0]*100:.1f}% var)" if len(explained) > 0 else "PC1"
                    pc2_label = f"PC2 ({explained[1]*100:.1f}% var)" if len(explained) > 1 else "PC2"
                    fig.update_xaxes(title_text=pc1_label, row=current_row, col=1)
                    fig.update_yaxes(title_text=pc2_label, row=current_row, col=1)

                    x_vals = components[:, 0]
                    y_vals = components[:, 1]
                    x_min, x_max = float(x_vals.min()), float(x_vals.max())
                    y_min, y_max = float(y_vals.min()), float(y_vals.max())
                    x_pad = (x_max - x_min) * 0.15 if x_max > x_min else 0.5
                    y_pad = (y_max - y_min) * 0.15 if y_max > y_min else 0.5
                    fig.update_xaxes(range=[x_min - x_pad, x_max + x_pad], row=current_row, col=1)
                    fig.update_yaxes(range=[y_min - y_pad, y_max + y_pad], row=current_row, col=1)

                    # Small legend inside this PCA subplot (group names only)
                    if groups is not None and group_colors is not None:
                        unique_groups = list(dict.fromkeys(groups))
                        legend_text = "<br>".join(unique_groups)
                        xref = "x domain" if current_row == 1 else f"x{current_row} domain"
                        yref = "y domain" if current_row == 1 else f"y{current_row} domain"
                        fig.add_annotation(
                            text=legend_text,
                            xref=xref, yref=yref,
                            x=0.99, y=0.99,
                            xanchor="right", yanchor="top",
                            showarrow=False,
                            font=dict(size=10),
                            bgcolor="rgba(255,255,255,0.85)",
                            bordercolor="rgba(0,0,0,0.2)",
                            borderwidth=1,
                            borderpad=4
                        )

                except Exception:
                    fig.add_annotation(
                        text="PCA failed",
                        x=0.5,
                        y=0.5,
                        showarrow=False,
                        row=current_row,
                        col=1
                    )
                    fig.update_xaxes(title_text="PC1", row=current_row, col=1)
                    fig.update_yaxes(title_text="PC2", row=current_row, col=1)
            else:
                fig.add_annotation(
                    text="Insufficient data for PCA",
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    row=current_row,
                    col=1
                )
                fig.update_xaxes(title_text="PC1", row=current_row, col=1)
                fig.update_yaxes(title_text="PC2", row=current_row, col=1)

            current_row += 1

            # 3. Heatmap
            if log_intensity_matrix.size > 0 and log_intensity_matrix.shape[1] > 0:
                # Filter out zero/near-zero variance rows (all-zero taxa) before clustering
                row_variances = np.var(log_intensity_matrix, axis=1)
                non_zero_var_mask = row_variances > 0
                heatmap_values = log_intensity_matrix[non_zero_var_mask].copy()
                heatmap_labels = [taxon_labels[i] for i in range(len(taxon_labels)) if non_zero_var_mask[i]]

                if heatmap_values.shape[0] > 1 and heatmap_values.size > 0:
                    from scipy.cluster.hierarchy import linkage, dendrogram
                    from scipy.spatial.distance import pdist
                    row_linkage = linkage(pdist(heatmap_values), method='ward')
                    row_order = dendrogram(row_linkage, no_plot=True)['leaves']
                    heatmap_values = heatmap_values[row_order]
                    heatmap_labels = [heatmap_labels[i] for i in row_order]

                col_order = list(range(heatmap_values.shape[1]))
                if heatmap_values.shape[1] > 1:
                    from scipy.cluster.hierarchy import linkage, dendrogram
                    from scipy.spatial.distance import pdist
                    col_linkage = linkage(pdist(heatmap_values.T), method='ward')
                    col_order = dendrogram(col_linkage, no_plot=True)['leaves']
                    heatmap_values = heatmap_values[:, col_order]
                
                sample_names_clustered = [
                    simple_names[i] if i < len(simple_names) else f"Sample_{i+1}"
                    for i in col_order
                ]
                # X-axis: experiment name only (no group/Condition appended)
                heatmap_x_labels = sample_names_clustered

                fig.add_trace(
                    go.Heatmap(
                        z=heatmap_values.tolist(),
                        x=heatmap_x_labels,
                        y=heatmap_labels,
                        colorscale='Viridis',
                        showscale=False,
                        name=f"{rank.capitalize()} Log Intensity Heatmap (Clustered)",
                        hovertemplate='<b>%{y}</b><br><b>Sample:</b> %{x}<br><b>Log Intensity:</b> %{z:.2f}<br><extra></extra>'
                    ),
                    row=current_row,
                    col=1
                )
            else:
                fig.add_annotation(
                    text="Insufficient data for heatmap",
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    row=current_row,
                    col=1
                )

            fig.update_xaxes(title_text="Samples", row=current_row, col=1, tickangle=45)
            fig.update_yaxes(title_text="Taxa", row=current_row, col=1)
            current_row += 1

            # 4. Composition bar chart
            sample_compositions = {}
            taxon_list = rank_summary[rank].tolist()
            for i, col in enumerate(intensity_cols):
                sample_values = rank_summary[col]
                total_sample = sample_values.sum()
                if total_sample > 0:
                    sample_taxa_abundance = {}
                    for idx, taxon in enumerate(taxon_list):
                        taxon_abundance = sample_values.iloc[idx]
                        if taxon_abundance > 0:
                            sample_taxa_abundance[taxon] = taxon_abundance
                    
                    top_sample_taxa = sorted(sample_taxa_abundance.items(), key=lambda x: x[1], reverse=True)[:9]
                    
                    rel_abundances = []
                    taxon_names = []
                    for taxon_name, abundance in top_sample_taxa:
                        rel_abundances.append(abundance / total_sample * 100)
                        taxon_names.append(taxon_name)
                    
                    top_total = sum(ab for _, ab in top_sample_taxa)
                    others_abundance = (total_sample - top_total) / total_sample * 100
                    if others_abundance > 0:
                        rel_abundances.append(others_abundance)
                        taxon_names.append("Others")
                    
                    file_id = simple_names[i] if i < len(simple_names) else f"File_{i+1}"
                    sample_compositions[file_id] = {
                        'taxa': taxon_names,
                        'abundances': rel_abundances
                    }

            if sample_compositions:
                sample_names_list = list(sample_compositions.keys())
                all_taxa = set()
                for comp in sample_compositions.values():
                    all_taxa.update(comp['taxa'])
                all_taxa = list(all_taxa)

                color_palette = [
                    '#B5BD00', '#8ED1FC', '#FF6B6B', '#4ECDC4', '#45B7D1',
                    '#96CEB4', '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F'
                ]

                for idx, taxon in enumerate(all_taxa):
                    abundances = []
                    for sample_name in sample_names_list:
                        if taxon in sample_compositions[sample_name]['taxa']:
                            taxon_idx = sample_compositions[sample_name]['taxa'].index(taxon)
                            abundances.append(sample_compositions[sample_name]['abundances'][taxon_idx])
                        else:
                            abundances.append(0)
                    
                    fig.add_trace(
                        go.Bar(
                            x=sample_names_list,
                            y=abundances,
                            name=taxon,
                            marker_color=color_palette[idx % len(color_palette)],
                            showlegend=False,
                            hovertemplate=f'<b>{taxon}</b><br>Sample: %{{x}}<br>Abundance: %{{y:.1f}}%<extra></extra>'
                        ),
                        row=current_row,
                        col=1
                    )
                
                # Constrain y-axis to [0, 100] for percentage display
                fig.update_yaxes(range=[0, 100], row=current_row, col=1)
            else:
                fig.add_annotation(
                    text="Insufficient data",
                    x=0.5,
                    y=0.5,
                    showarrow=False,
                    row=current_row,
                    col=1
                )

            fig.update_xaxes(title_text="Samples", row=current_row, col=1)
            fig.update_yaxes(title_text="Relative Abundance (%)", row=current_row, col=1)
            current_row += 1

        height_per_subplot = 900
        # No global Sample Groups legend; each PCA has a small in-figure legend
        fig.update_layout(
            height=height_per_subplot * total_rows,
            showlegend=False,
            title_text="Taxonomy Summary Report",
            barmode='stack'
        )
        
        return fig

    def _add_confidence_ellipse(self, fig, x_coords, y_coords, color, group_name, row):
        """Add a 95% confidence ellipse around group points in PCA plot."""
        try:
            if len(x_coords) < 3:
                return
            
            x_arr = np.array(x_coords)
            y_arr = np.array(y_coords)
            
            # Calculate mean and covariance
            mean_x = np.mean(x_arr)
            mean_y = np.mean(y_arr)
            
            # Calculate covariance matrix
            cov_matrix = np.cov(x_arr, y_arr)
            
            # Get eigenvalues and eigenvectors
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
            
            # Sort by eigenvalue (largest first)
            order = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[order]
            eigenvectors = eigenvectors[:, order]
            
            # Calculate ellipse parameters (95% confidence = 2.447 standard deviations)
            chi2_val = 2.447  # 95% confidence for 2D
            width = 2 * chi2_val * np.sqrt(eigenvalues[0])
            height = 2 * chi2_val * np.sqrt(eigenvalues[1])
            
            # Calculate rotation angle
            angle = np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])
            
            # Generate ellipse points
            t = np.linspace(0, 2 * np.pi, 50)
            ellipse_x = width/2 * np.cos(t)
            ellipse_y = height/2 * np.sin(t)
            
            # Rotate ellipse
            cos_angle = np.cos(angle)
            sin_angle = np.sin(angle)
            rotated_x = cos_angle * ellipse_x - sin_angle * ellipse_y + mean_x
            rotated_y = sin_angle * ellipse_x + cos_angle * ellipse_y + mean_y
            
            # Add ellipse trace
            fig.add_trace(
                go.Scatter(
                    x=rotated_x.tolist(),
                    y=rotated_y.tolist(),
                    mode='lines',
                    line=dict(color=color, width=2, dash='dash'),
                    fill='toself',
                    fillcolor=color.replace(')', ', 0.1)').replace('rgb', 'rgba') if 'rgb' in color else f"rgba{tuple(list(int(color.lstrip('#')[i:i+2], 16) for i in (0, 2, 4)) + [0.1])}",
                    name=f"{group_name} (95% CI)",
                    legendgroup=group_name,
                    showlegend=False,
                    hoverinfo='skip'
                ),
                row=row,
                col=1
            )
        except Exception as e:
            # Silently skip ellipse if calculation fails
            print(f"Warning: Could not add confidence ellipse for {group_name}: {e}")

