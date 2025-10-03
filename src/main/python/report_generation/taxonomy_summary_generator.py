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
        
        subplot_spacing = min(0.06, 0.8 / max(total_rows - 1, 1)) if total_rows > 1 else 0
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
                            name=sample_name,
                            box_visible=True,
                            meanline_visible=True,
                            points='all',
                            pointpos=0,
                            jitter=0.3,
                            marker_color=self.colors["primary"],
                            showlegend=False,
                            x0=sample_name
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

            # 2. PCA scatter
            if log_intensity_matrix.shape[0] >= 2 and log_intensity_matrix.shape[1] >= 2:
                try:
                    scaled = StandardScaler().fit_transform(log_intensity_matrix)
                    pca = PCA(n_components=2)
                    components = pca.fit_transform(scaled)
                    explained = pca.explained_variance_ratio_

                    fig.add_trace(
                        go.Scatter(
                            x=components[:, 0],
                            y=components[:, 1],
                            mode='markers',
                            marker=dict(color=self.colors['secondary'], size=9, line=dict(width=1, color='#333333')),
                            text=taxon_labels,
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
                    x_pad = (x_max - x_min) * 0.1 if x_max > x_min else 0.5
                    y_pad = (y_max - y_min) * 0.1 if y_max > y_min else 0.5
                    fig.update_xaxes(range=[x_min - x_pad, x_max + x_pad], row=current_row, col=1)
                    fig.update_yaxes(range=[y_min - y_pad, y_max + y_pad], row=current_row, col=1)

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
                heatmap_values = log_intensity_matrix.copy()
                heatmap_labels = taxon_labels.copy()

                if heatmap_values.shape[0] > 1:
                    from scipy.cluster.hierarchy import linkage, dendrogram
                    from scipy.spatial.distance import pdist
                    row_linkage = linkage(pdist(heatmap_values), method='ward')
                    row_order = dendrogram(row_linkage, no_plot=True)['leaves']
                    heatmap_values = heatmap_values[row_order]
                    heatmap_labels = [heatmap_labels[i] for i in row_order]

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
                else:
                    sample_names_clustered = (
                        simple_names[:len(intensity_cols)]
                        if len(simple_names) > 0
                        else [f"Sample_{i+1}" for i in range(len(intensity_cols))]
                    )

                fig.add_trace(
                    go.Heatmap(
                        z=heatmap_values.tolist(),
                        x=sample_names_clustered,
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

            fig.update_xaxes(title_text="Samples", row=current_row, col=1)
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
                            showlegend=False
                        ),
                        row=current_row,
                        col=1
                    )
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
        fig.update_layout(
            height=height_per_subplot * total_rows,
            showlegend=False,
            title_text="Taxonomy Summary Report",
            barmode='stack'
        )
        
        return fig

