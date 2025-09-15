import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
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
            rows=4, cols=2,
            subplot_titles=(
                'Protein Counts per Raw File',
                'Peptides per Protein Distribution',
                'PCA Analysis of Samples',
                'UMAP Analysis of Samples' if UMAP_AVAILABLE else 'Sample Clustering',
                'Sample Clustering Dendrogram',
                'PSMs per Protein',
                'Protein Intensity Distribution',
                'Group Comparison (if metadata available)'
            ),
            specs=[[{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(proteins_df, metadata)

        # Protein counts per raw file
        if len(intensity_cols) > 0:
            protein_counts = [((proteins_df[col] > 0).sum()) for col in intensity_cols]
            
            colors = [group_colors[groups[i]] if groups and group_colors else self.colors['primary'] 
                     for i in range(len(sample_names))]
            
            fig.add_trace(go.Bar(x=sample_names, y=protein_counts, name='Protein Counts',
                               marker_color=colors), row=1, col=1)

        # Peptides per protein distribution
        all_pep_cols = ['Peptide counts (all)', 'Peptides', 'Peptide count']
        all_pep_col = next((col for col in all_pep_cols if col in proteins_df.columns), None)
        
        if all_pep_col:
            fig.add_trace(go.Histogram(x=proteins_df[all_pep_col], name='Peptides per Protein',
                                     marker_color=self.colors['secondary'], nbinsx=25), row=1, col=2)

        # Advanced dimensionality reduction analysis
        if len(intensity_cols) >= 2:
            try:
                # Prepare data
                intensity_data = proteins_df[intensity_cols].fillna(0)
                intensity_data = intensity_data.loc[intensity_data.sum(axis=1) > 0]
                
                if len(intensity_data) > 0:
                    log_data = np.log10(intensity_data + 1)
                    scaler = StandardScaler()
                    scaled_data = scaler.fit_transform(log_data.T)
                    
                    # PCA Analysis
                    pca = PCA(n_components=min(2, len(intensity_cols)))
                    pca_result = pca.fit_transform(scaled_data)
                    
                    # Color points by groups if available
                    if groups and group_colors:
                        colors = [group_colors[group] for group in groups]
                        for group in set(groups):
                            group_indices = [i for i, g in enumerate(groups) if g == group]
                            fig.add_trace(go.Scatter(
                                x=[pca_result[i, 0] for i in group_indices],
                                y=[pca_result[i, 1] for i in group_indices] if pca_result.shape[1] > 1 else [0] * len(group_indices),
                                mode='markers+text',
                                text=[sample_names[i] for i in group_indices],
                                textposition="top center",
                                marker=dict(size=12, color=group_colors[group], line=dict(width=2, color='white')),
                                name=f'Group {group}',
                                showlegend=True
                            ), row=2, col=1)
                    else:
                        fig.add_trace(go.Scatter(
                            x=pca_result[:, 0],
                            y=pca_result[:, 1] if pca_result.shape[1] > 1 else np.zeros(len(pca_result)),
                            mode='markers+text',
                            text=sample_names,
                            textposition="top center",
                            marker=dict(size=12, color=self.colors['tertiary'], line=dict(width=2, color='white')),
                            name='Samples'
                        ), row=2, col=1)
                    
                    # UMAP Analysis (if available)
                    if UMAP_AVAILABLE and len(scaled_data) > 3:
                        try:
                            umap_reducer = UMAP(n_components=2, random_state=42)
                            umap_result = umap_reducer.fit_transform(scaled_data)
                            
                            if groups and group_colors:
                                for group in set(groups):
                                    group_indices = [i for i, g in enumerate(groups) if g == group]
                                    fig.add_trace(go.Scatter(
                                        x=[umap_result[i, 0] for i in group_indices],
                                        y=[umap_result[i, 1] for i in group_indices],
                                        mode='markers+text',
                                        text=[sample_names[i] for i in group_indices],
                                        textposition="top center",
                                        marker=dict(size=12, color=group_colors[group], line=dict(width=2, color='white')),
                                        name=f'Group {group}',
                                        showlegend=False
                                    ), row=2, col=2)
                            else:
                                fig.add_trace(go.Scatter(
                                    x=umap_result[:, 0], y=umap_result[:, 1],
                                    mode='markers+text', text=sample_names, textposition="top center",
                                    marker=dict(size=12, color=self.colors['success'], line=dict(width=2, color='white')),
                                    name='UMAP'
                                ), row=2, col=2)
                        except Exception as e:
                            fig.add_annotation(text=f"UMAP failed: {str(e)[:30]}...", 
                                             xref="x4", yref="y4", x=0.5, y=0.5, showarrow=False, row=2, col=2)
                    
                    # Hierarchical Clustering
                    distance_matrix = pdist(log_data.T, metric='euclidean')
                    linkage_matrix = linkage(distance_matrix, method='ward')
                    dendro = dendrogram(linkage_matrix, labels=sample_names, no_plot=True)
                    
                    icoord = np.array(dendro['icoord'])
                    dcoord = np.array(dendro['dcoord'])
                    
                    for i in range(len(icoord)):
                        fig.add_trace(go.Scatter(
                            x=[icoord[i][0], icoord[i][1], icoord[i][2], icoord[i][3]],
                            y=[dcoord[i][0], dcoord[i][1], dcoord[i][2], dcoord[i][3]],
                            mode='lines', line=dict(color=self.colors['success'], width=2),
                            showlegend=False
                        ), row=3, col=1)
                        
            except Exception as e:
                fig.add_annotation(text=f"Analysis failed: {str(e)[:50]}...", 
                                 xref="x3", yref="y3", x=0.5, y=0.5, showarrow=False, row=2, col=1)

        # PSM and Intensity distributions
        psm_cols = ['Number of PSMs', 'PSMs', 'MS/MS count']
        psm_col = next((col for col in psm_cols if col in proteins_df.columns), None)
        
        if psm_col:
            fig.add_trace(go.Histogram(x=proteins_df[psm_col], name='PSMs per Protein',
                                     marker_color=self.colors['light_blue'], nbinsx=25), row=3, col=2)

        intensity_col = next((col for col in ['Intensity', 'Total Intensity'] if col in proteins_df.columns), None)
        
        if intensity_col:
            intensities = proteins_df[intensity_col][proteins_df[intensity_col] > 0]
            fig.add_trace(go.Histogram(x=np.log10(intensities), name='Log10(Intensity)',
                                     marker_color=self.colors['light_purple'], nbinsx=40), row=4, col=1)

        # Group comparison box plot
        if groups and group_colors and intensity_col:
            for group in set(groups):
                group_samples = [sample_names[i] for i, g in enumerate(groups) if g == group]
                group_cols = [col for col in intensity_cols if any(sample in col for sample in group_samples)]
                
                if group_cols:
                    group_intensities = proteins_df[group_cols].sum(axis=1)
                    group_intensities = group_intensities[group_intensities > 0]
                    
                    fig.add_trace(go.Box(y=np.log10(group_intensities), name=f'Group {group}',
                                       marker_color=group_colors[group]), row=4, col=2)

        fig.update_layout(height=1600, title_text="MetaPilot Proteins Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Raw Files", row=1, col=1)
        fig.update_xaxes(title_text="PC1", row=2, col=1)
        fig.update_xaxes(title_text="UMAP1" if UMAP_AVAILABLE else "Samples", row=2, col=2)
        fig.update_yaxes(title_text="PC2", row=2, col=1)
        fig.update_yaxes(title_text="UMAP2" if UMAP_AVAILABLE else "Distance", row=2, col=2)

        return fig
