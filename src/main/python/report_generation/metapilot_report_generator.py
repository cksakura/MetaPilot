import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from pathlib import Path
import json
from datetime import datetime
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
try:
    from umap import UMAP
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False
try:
    from sklearn.manifold import TSNE
    TSNE_AVAILABLE = True
except ImportError:
    TSNE_AVAILABLE = False

class MetaPilotReportGenerator:
    """Enhanced report generator for MetaPilot TSV file formats with metadata support."""
    
    def __init__(self, output_dir="reports"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # MetaPilot color scheme
        self.colors = {
            'primary': '#2E86AB',
            'secondary': '#A23B72', 
            'tertiary': '#F18F01',
            'success': '#C73E1D',
            'background': '#f8f9fa',
            'text': '#2c3e50',
            'light_blue': '#87CEEB',
            'light_purple': '#DDA0DD',
            'light_orange': '#FFE4B5'
        }
        
        # Color palette for groups
        self.group_colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ]

    def load_metadata(self, metadata_file):
        """Load metadata file if available."""
        try:
            if metadata_file and os.path.exists(metadata_file):
                metadata = pd.read_csv(metadata_file, sep='\t')
                return metadata
        except Exception as e:
            print(f"Warning: Could not load metadata: {e}")
        return None

    def extract_sample_info(self, df, metadata=None):
        """Extract sample names and group information."""
        # Find intensity columns (representing raw files/samples)
        intensity_cols = []
        for col in df.columns:
            if 'Intensity' in col and col != 'Intensity':
                intensity_cols.append(col)
        
        if len(intensity_cols) == 0:
            for col in df.columns:
                if col.startswith('Intensity '):
                    intensity_cols.append(col)
        
        sample_names = [col.replace('Intensity ', '').replace('Intensity', '') for col in intensity_cols]
        
        # Extract group information from metadata
        groups = None
        group_colors = None
        if metadata is not None and len(sample_names) > 0:
            try:
                # Match samples with metadata
                sample_col = 'Sample' if 'Sample' in metadata.columns else metadata.columns[0]
                group_col = 'Group' if 'Group' in metadata.columns else None
                
                if group_col and group_col in metadata.columns:
                    sample_groups = {}
                    for _, row in metadata.iterrows():
                        sample_groups[row[sample_col]] = row[group_col]
                    
                    groups = [sample_groups.get(sample, 'Unknown') for sample in sample_names]
                    unique_groups = list(set(groups))
                    group_colors = {group: self.group_colors[i % len(self.group_colors)] 
                                  for i, group in enumerate(unique_groups)}
            except Exception as e:
                print(f"Warning: Could not extract group information: {e}")
        
        return intensity_cols, sample_names, groups, group_colors

    def create_id_summary(self, summary_df, metadata=None):
        """Generate ID summary report with metadata support."""
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
            total_ids = [summary_df[sample].sum() if pd.api.types.is_numeric_dtype(summary_df[sample]) else 0 
                        for sample in sample_names]
            
            # Color by groups if available
            colors = [group_colors[groups[i]] if groups and group_colors else self.colors['primary'] 
                     for i in range(len(sample_names))]
            
            fig.add_trace(go.Bar(x=sample_names, y=total_ids, name='Total IDs', 
                               marker_color=colors), row=1, col=1)

        fig.update_layout(height=800, showlegend=False, title_text="ID Summary Report")
        return fig

    def create_peptides_summary(self, peptides_df, metadata=None):
        """Enhanced peptides summary with metadata support."""
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

        # Statistical distributions (unchanged)
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

    def create_proteins_summary(self, proteins_df, metadata=None):
        """Enhanced proteins summary with advanced dimensionality reduction."""
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

    def create_taxonomy_summary(self, taxa_df, metadata=None):
        """Enhanced taxonomy summary with hierarchical visualizations."""
        # Create complex subplot layout
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=(
                'Taxonomic Composition (Treemap)',
                'Relative Abundance by Rank',
                'Top 20 Taxa by Abundance',
                'PCA Analysis of Taxonomic Profiles',
                'Sample Clustering (Taxonomy)',
                'Diversity Indices by Sample'
            ),
            specs=[[{"type": "treemap"}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(taxa_df, metadata)
        
        # Find taxonomic columns
        taxon_cols = ['Taxon', 'taxon', 'Name', 'Taxa', 'Species']
        taxon_col = next((col for col in taxon_cols if col in taxa_df.columns), None)
        
        rank_cols = ['Rank', 'rank', 'Level', 'Taxonomic_Rank']
        rank_col = next((col for col in rank_cols if col in taxa_df.columns), None)

        if len(intensity_cols) > 0 and taxon_col:
            # Calculate total abundance
            taxa_df['Total_Abundance'] = taxa_df[intensity_cols].sum(axis=1)
            
            # Calculate relative abundance
            total_sum = taxa_df['Total_Abundance'].sum()
            taxa_df['Relative_Abundance'] = (taxa_df['Total_Abundance'] / total_sum) * 100
            
            # Treemap visualization
            top_taxa = taxa_df.nlargest(30, 'Total_Abundance')  # Top 30 for treemap
            
            fig.add_trace(go.Treemap(
                labels=top_taxa[taxon_col],
                values=top_taxa['Relative_Abundance'],
                parents=[""] * len(top_taxa),  # All at root level
                textinfo="label+percent entry",
                marker=dict(colorscale="Viridis"),
                name="Taxonomic Composition"
            ), row=1, col=1)
            
            # Relative abundance by rank
            if rank_col:
                rank_abundance = taxa_df.groupby(rank_col)['Relative_Abundance'].sum().sort_values(ascending=False)
                
                fig.add_trace(go.Bar(
                    x=rank_abundance.index,
                    y=rank_abundance.values,
                    name='Abundance by Rank',
                    marker_color=self.colors['secondary']
                ), row=1, col=2)
            
            # Top taxa bar chart
            top_20_taxa = taxa_df.nlargest(20, 'Total_Abundance')
            
            fig.add_trace(go.Bar(
                x=top_20_taxa['Relative_Abundance'],
                y=top_20_taxa[taxon_col],
                orientation='h',
                name='Top Taxa',
                marker_color=self.colors['tertiary']
            ), row=2, col=1)
            
            # PCA analysis of taxonomic profiles
            if len(intensity_cols) >= 2:
                try:
                    # Prepare taxonomic profile data
                    profile_data = taxa_df[intensity_cols].fillna(0)
                    profile_data = profile_data.loc[profile_data.sum(axis=1) > 0]
                    
                    if len(profile_data) > 0:
                        # Relative abundance transformation
                        rel_abundance = profile_data.div(profile_data.sum(axis=0), axis=1) * 100
                        
                        # Log transform and standardize
                        log_data = np.log10(rel_abundance + 0.01)  # Small pseudocount
                        scaler = StandardScaler()
                        scaled_data = scaler.fit_transform(log_data.T)
                        
                        # PCA
                        pca = PCA(n_components=min(2, len(intensity_cols)))
                        pca_result = pca.fit_transform(scaled_data)
                        
                        if groups and group_colors:
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
                                ), row=2, col=2)
                        else:
                            fig.add_trace(go.Scatter(
                                x=pca_result[:, 0],
                                y=pca_result[:, 1] if pca_result.shape[1] > 1 else np.zeros(len(pca_result)),
                                mode='markers+text', text=sample_names, textposition="top center",
                                marker=dict(size=12, color=self.colors['primary'], line=dict(width=2, color='white')),
                                name='Taxonomic PCA'
                            ), row=2, col=2)
                            
                        # Hierarchical clustering
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
                    fig.add_annotation(text=f"Taxonomic analysis failed: {str(e)[:40]}...",
                                     xref="x4", yref="y4", x=0.5, y=0.5, showarrow=False, row=2, col=2)
            
            # Diversity indices
            if len(intensity_cols) > 0:
                # Calculate Shannon diversity for each sample
                shannon_diversity = []
                for col in intensity_cols:
                    sample_data = taxa_df[col][taxa_df[col] > 0]
                    if len(sample_data) > 0:
                        rel_abund = sample_data / sample_data.sum()
                        shannon = -np.sum(rel_abund * np.log(rel_abund))
                        shannon_diversity.append(shannon)
                    else:
                        shannon_diversity.append(0)
                
                colors = [group_colors[groups[i]] if groups and group_colors else self.colors['light_blue'] 
                         for i in range(len(sample_names))]
                
                fig.add_trace(go.Bar(
                    x=sample_names,
                    y=shannon_diversity,
                    name='Shannon Diversity',
                    marker_color=colors
                ), row=3, col=2)

        fig.update_layout(height=1400, title_text="MetaPilot Taxonomy Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Taxonomic Rank", row=1, col=2)
        fig.update_xaxes(title_text="Relative Abundance (%)", row=2, col=1)
        fig.update_xaxes(title_text="PC1", row=2, col=2)
        fig.update_xaxes(title_text="Samples", row=3, col=1)
        fig.update_xaxes(title_text="Samples", row=3, col=2)
        
        fig.update_yaxes(title_text="Relative Abundance (%)", row=1, col=2)
        fig.update_yaxes(title_text="Taxa", row=2, col=1)
        fig.update_yaxes(title_text="PC2", row=2, col=2)
        fig.update_yaxes(title_text="Distance", row=3, col=1)
        fig.update_yaxes(title_text="Shannon Index", row=3, col=2)

        return fig

    def create_function_summary(self, func_df, metadata=None):
        """Enhanced function summary with category analysis and dimensionality reduction."""
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=(
                'Functional Categories Distribution',
                'Top 20 Functions by Abundance',
                'PCA Analysis of Functional Profiles',
                'UMAP Analysis of Functions' if UMAP_AVAILABLE else 'Function Clustering',
                'Sample Clustering (Functions)',
                'Pathway Enrichment Analysis'
            ),
            specs=[[{"type": "pie"}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(func_df, metadata)
        
        # Find function columns
        func_cols = ['Function', 'function', 'Description', 'Name', 'Gene_Name']
        func_col = next((col for col in func_cols if col in func_df.columns), None)
        
        cat_cols = ['Category', 'category', 'Class', 'Type', 'COG_Category', 'KEGG_Category']
        cat_col = next((col for col in cat_cols if col in func_df.columns), None)

        if len(intensity_cols) > 0 and func_col:
            # Calculate total abundance
            func_df['Total_Abundance'] = func_df[intensity_cols].sum(axis=1)
            
            # Functional categories pie chart
            if cat_col:
                cat_abundance = func_df.groupby(cat_col)['Total_Abundance'].sum()
                
                fig.add_trace(go.Pie(
                    labels=cat_abundance.index,
                    values=cat_abundance.values,
                    name="Functional Categories",
                    marker=dict(colors=px.colors.qualitative.Set3)
                ), row=1, col=1)
            
            # Top functions bar chart
            top_20_funcs = func_df.nlargest(20, 'Total_Abundance')
            
            fig.add_trace(go.Bar(
                x=top_20_funcs['Total_Abundance'],
                y=top_20_funcs[func_col],
                orientation='h',
                name='Top Functions',
                marker_color=self.colors['tertiary']
            ), row=1, col=2)
            
            # PCA analysis of functional profiles
            if len(intensity_cols) >= 2:
                try:
                    # Prepare functional profile data
                    profile_data = func_df[intensity_cols].fillna(0)
                    profile_data = profile_data.loc[profile_data.sum(axis=1) > 0]
                    
                    if len(profile_data) > 0:
                        # Log transform and standardize
                        log_data = np.log10(profile_data + 1)
                        scaler = StandardScaler()
                        scaled_data = scaler.fit_transform(log_data.T)
                        
                        # PCA
                        pca = PCA(n_components=min(2, len(intensity_cols)))
                        pca_result = pca.fit_transform(scaled_data)
                        
                        if groups and group_colors:
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
                                mode='markers+text', text=sample_names, textposition="top center",
                                marker=dict(size=12, color=self.colors['primary'], line=dict(width=2, color='white')),
                                name='Functional PCA'
                            ), row=2, col=1)
                        
                        # UMAP for functions (if available)
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
                                        name='Functional UMAP'
                                    ), row=2, col=2)
                            except Exception as e:
                                fig.add_annotation(text=f"UMAP failed: {str(e)[:30]}...",
                                                 xref="x4", yref="y4", x=0.5, y=0.5, showarrow=False, row=2, col=2)
                        
                        # Hierarchical clustering
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
                    fig.add_annotation(text=f"Functional analysis failed: {str(e)[:40]}...",
                                     xref="x3", yref="y3", x=0.5, y=0.5, showarrow=False, row=2, col=1)
            
            # Pathway enrichment (simplified version)
            if cat_col and groups:
                # Compare functional categories between groups
                try:
                    group_comparison = {}
                    for group in set(groups):
                        group_samples = [sample_names[i] for i, g in enumerate(groups) if g == group]
                        group_cols = [col for col in intensity_cols if any(sample in col for sample in group_samples)]
                        
                        if group_cols:
                            group_func_data = func_df[group_cols].sum(axis=1)
                            group_cat_abundance = func_df.groupby(cat_col)[group_cols].sum().sum(axis=1)
                            group_comparison[group] = group_cat_abundance
                    
                    if len(group_comparison) > 1:
                        comparison_df = pd.DataFrame(group_comparison).fillna(0)
                        
                        # Create heatmap-like visualization
                        for i, (category, row) in enumerate(comparison_df.iterrows()):
                            fig.add_trace(go.Bar(
                                x=list(row.index),
                                y=[row[group] for group in row.index],
                                name=category,
                                offsetgroup=i
                            ), row=3, col=2)
                            
                except Exception as e:
                    fig.add_annotation(text=f"Group comparison failed: {str(e)[:30]}...",
                                     xref="x6", yref="y6", x=0.5, y=0.5, showarrow=False, row=3, col=2)

        fig.update_layout(height=1400, title_text="MetaPilot Function Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Abundance", row=1, col=2)
        fig.update_xaxes(title_text="PC1", row=2, col=1)
        fig.update_xaxes(title_text="UMAP1" if UMAP_AVAILABLE else "Samples", row=2, col=2)
        fig.update_xaxes(title_text="Groups", row=3, col=2)
        
        fig.update_yaxes(title_text="Functions", row=1, col=2)
        fig.update_yaxes(title_text="PC2", row=2, col=1)
        fig.update_yaxes(title_text="UMAP2" if UMAP_AVAILABLE else "Distance", row=2, col=2)
        fig.update_yaxes(title_text="Distance", row=3, col=1)
        fig.update_yaxes(title_text="Abundance", row=3, col=2)

        return fig

    def create_taxonomy_summary(self, taxa_df, metadata=None):
        """Enhanced taxonomy summary with hierarchical visualizations."""
        # Create complex subplot layout
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=(
                'Taxonomic Composition (Treemap)',
                'Relative Abundance by Rank',
                'Top 20 Taxa by Abundance',
                'PCA Analysis of Taxonomic Profiles',
                'Sample Clustering (Taxonomy)',
                'Diversity Indices by Sample'
            ),
            specs=[[{"type": "treemap"}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(taxa_df, metadata)
        
        # Find taxonomic columns
        taxon_cols = ['Taxon', 'taxon', 'Name', 'Taxa', 'Species']
        taxon_col = next((col for col in taxon_cols if col in taxa_df.columns), None)
        
        rank_cols = ['Rank', 'rank', 'Level', 'Taxonomic_Rank']
        rank_col = next((col for col in rank_cols if col in taxa_df.columns), None)

        if len(intensity_cols) > 0 and taxon_col:
            # Calculate total abundance
            taxa_df['Total_Abundance'] = taxa_df[intensity_cols].sum(axis=1)
            
            # Calculate relative abundance
            total_sum = taxa_df['Total_Abundance'].sum()
            taxa_df['Relative_Abundance'] = (taxa_df['Total_Abundance'] / total_sum) * 100
            
            # Treemap visualization
            top_taxa = taxa_df.nlargest(30, 'Total_Abundance')  # Top 30 for treemap
            
            fig.add_trace(go.Treemap(
                labels=top_taxa[taxon_col],
                values=top_taxa['Relative_Abundance'],
                parents=[""] * len(top_taxa),  # All at root level
                textinfo="label+percent entry",
                marker=dict(colorscale="Viridis"),
                name="Taxonomic Composition"
            ), row=1, col=1)
            
            # Relative abundance by rank
            if rank_col:
                rank_abundance = taxa_df.groupby(rank_col)['Relative_Abundance'].sum().sort_values(ascending=False)
                
                fig.add_trace(go.Bar(
                    x=rank_abundance.index,
                    y=rank_abundance.values,
                    name='Abundance by Rank',
                    marker_color=self.colors['secondary']
                ), row=1, col=2)
            
            # Top taxa bar chart
            top_20_taxa = taxa_df.nlargest(20, 'Total_Abundance')
            
            fig.add_trace(go.Bar(
                x=top_20_taxa['Relative_Abundance'],
                y=top_20_taxa[taxon_col],
                orientation='h',
                name='Top Taxa',
                marker_color=self.colors['tertiary']
            ), row=2, col=1)
            
            # PCA analysis of taxonomic profiles
            if len(intensity_cols) >= 2:
                try:
                    # Prepare taxonomic profile data
                    profile_data = taxa_df[intensity_cols].fillna(0)
                    profile_data = profile_data.loc[profile_data.sum(axis=1) > 0]
                    
                    if len(profile_data) > 0:
                        # Relative abundance transformation
                        rel_abundance = profile_data.div(profile_data.sum(axis=0), axis=1) * 100
                        
                        # Log transform and standardize
                        log_data = np.log10(rel_abundance + 0.01)  # Small pseudocount
                        scaler = StandardScaler()
                        scaled_data = scaler.fit_transform(log_data.T)
                        
                        # PCA
                        pca = PCA(n_components=min(2, len(intensity_cols)))
                        pca_result = pca.fit_transform(scaled_data)
                        
                        if groups and group_colors:
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
                                ), row=2, col=2)
                        else:
                            fig.add_trace(go.Scatter(
                                x=pca_result[:, 0],
                                y=pca_result[:, 1] if pca_result.shape[1] > 1 else np.zeros(len(pca_result)),
                                mode='markers+text', text=sample_names, textposition="top center",
                                marker=dict(size=12, color=self.colors['primary'], line=dict(width=2, color='white')),
                                name='Taxonomic PCA'
                            ), row=2, col=2)
                            
                        # Hierarchical clustering
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
                    fig.add_annotation(text=f"Taxonomic analysis failed: {str(e)[:40]}...",
                                     xref="x4", yref="y4", x=0.5, y=0.5, showarrow=False, row=2, col=2)
            
            # Diversity indices
            if len(intensity_cols) > 0:
                # Calculate Shannon diversity for each sample
                shannon_diversity = []
                for col in intensity_cols:
                    sample_data = taxa_df[col][taxa_df[col] > 0]
                    if len(sample_data) > 0:
                        rel_abund = sample_data / sample_data.sum()
                        shannon = -np.sum(rel_abund * np.log(rel_abund))
                        shannon_diversity.append(shannon)
                    else:
                        shannon_diversity.append(0)
                
                colors = [group_colors[groups[i]] if groups and group_colors else self.colors['light_blue'] 
                         for i in range(len(sample_names))]
                
                fig.add_trace(go.Bar(
                    x=sample_names,
                    y=shannon_diversity,
                    name='Shannon Diversity',
                    marker_color=colors
                ), row=3, col=2)

        fig.update_layout(height=1400, title_text="MetaPilot Taxonomy Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Taxonomic Rank", row=1, col=2)
        fig.update_xaxes(title_text="Relative Abundance (%)", row=2, col=1)
        fig.update_xaxes(title_text="PC1", row=2, col=2)
        fig.update_xaxes(title_text="Samples", row=3, col=1)
        fig.update_xaxes(title_text="Samples", row=3, col=2)
        
        fig.update_yaxes(title_text="Relative Abundance (%)", row=1, col=2)
        fig.update_yaxes(title_text="Taxa", row=2, col=1)
        fig.update_yaxes(title_text="PC2", row=2, col=2)
        fig.update_yaxes(title_text="Distance", row=3, col=1)
        fig.update_yaxes(title_text="Shannon Index", row=3, col=2)

        return fig

    def create_function_summary(self, func_df, metadata=None):
        """Enhanced function summary with category analysis and dimensionality reduction."""
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=(
                'Functional Categories Distribution',
                'Top 20 Functions by Abundance',
                'PCA Analysis of Functional Profiles',
                'UMAP Analysis of Functions' if UMAP_AVAILABLE else 'Function Clustering',
                'Sample Clustering (Functions)',
                'Pathway Enrichment Analysis'
            ),
            specs=[[{"type": "pie"}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}],
                   [{"colspan": 1}, {"colspan": 1}]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(func_df, metadata)
        
        # Find function columns
        func_cols = ['Function', 'function', 'Description', 'Name', 'Gene_Name']
        func_col = next((col for col in func_cols if col in func_df.columns), None)
        
        cat_cols = ['Category', 'category', 'Class', 'Type', 'COG_Category', 'KEGG_Category']
        cat_col = next((col for col in cat_cols if col in func_df.columns), None)

        if len(intensity_cols) > 0 and func_col:
            # Calculate total abundance
            func_df['Total_Abundance'] = func_df[intensity_cols].sum(axis=1)
            
            # Functional categories pie chart
            if cat_col:
                cat_abundance = func_df.groupby(cat_col)['Total_Abundance'].sum()
                
                fig.add_trace(go.Pie(
                    labels=cat_abundance.index,
                    values=cat_abundance.values,
                    name="Functional Categories",
                    marker=dict(colors=px.colors.qualitative.Set3)
                ), row=1, col=1)
            
            # Top functions bar chart
            top_20_funcs = func_df.nlargest(20, 'Total_Abundance')
            
            fig.add_trace(go.Bar(
                x=top_20_funcs['Total_Abundance'],
                y=top_20_funcs[func_col],
                orientation='h',
                name='Top Functions',
                marker_color=self.colors['tertiary']
            ), row=1, col=2)
            
            # PCA analysis of functional profiles
            if len(intensity_cols) >= 2:
                try:
                    # Prepare functional profile data
                    profile_data = func_df[intensity_cols].fillna(0)
                    profile_data = profile_data.loc[profile_data.sum(axis=1) > 0]
                    
                    if len(profile_data) > 0:
                        # Log transform and standardize
                        log_data = np.log10(profile_data + 1)
                        scaler = StandardScaler()
                        scaled_data = scaler.fit_transform(log_data.T)
                        
                        # PCA
                        pca = PCA(n_components=min(2, len(intensity_cols)))
                        pca_result = pca.fit_transform(scaled_data)
                        
                        if groups and group_colors:
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
                                mode='markers+text', text=sample_names, textposition="top center",
                                marker=dict(size=12, color=self.colors['primary'], line=dict(width=2, color='white')),
                                name='Functional PCA'
                            ), row=2, col=1)
                        
                        # UMAP for functions (if available)
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
                                        name='Functional UMAP'
                                    ), row=2, col=2)
                            except Exception as e:
                                fig.add_annotation(text=f"UMAP failed: {str(e)[:30]}...",
                                                 xref="x4", yref="y4", x=0.5, y=0.5, showarrow=False, row=2, col=2)
                        
                        # Hierarchical clustering
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
                    fig.add_annotation(text=f"Functional analysis failed: {str(e)[:40]}...",
                                     xref="x3", yref="y3", x=0.5, y=0.5, showarrow=False, row=2, col=1)
            
            # Pathway enrichment (simplified version)
            if cat_col and groups:
                # Compare functional categories between groups
                try:
                    group_comparison = {}
                    for group in set(groups):
                        group_samples = [sample_names[i] for i, g in enumerate(groups) if g == group]
                        group_cols = [col for col in intensity_cols if any(sample in col for sample in group_samples)]
                        
                        if group_cols:
                            group_func_data = func_df[group_cols].sum(axis=1)
                            group_cat_abundance = func_df.groupby(cat_col)[group_cols].sum().sum(axis=1)
                            group_comparison[group] = group_cat_abundance
                    
                    if len(group_comparison) > 1:
                        comparison_df = pd.DataFrame(group_comparison).fillna(0)
                        
                        # Create heatmap-like visualization
                        for i, (category, row) in enumerate(comparison_df.iterrows()):
                            fig.add_trace(go.Bar(
                                x=list(row.index),
                                y=[row[group] for group in row.index],
                                name=category,
                                offsetgroup=i
                            ), row=3, col=2)
                            
                except Exception as e:
                    fig.add_annotation(text=f"Group comparison failed: {str(e)[:30]}...",
                                     xref="x6", yref="y6", x=0.5, y=0.5, showarrow=False, row=3, col=2)

        fig.update_layout(height=1400, title_text="MetaPilot Function Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Abundance", row=1, col=2)
        fig.update_xaxes(title_text="PC1", row=2, col=1)
        fig.update_xaxes(title_text="UMAP1" if UMAP_AVAILABLE else "Samples", row=2, col=2)
        fig.update_xaxes(title_text="Groups", row=3, col=2)
        
        fig.update_yaxes(title_text="Functions", row=1, col=2)
        fig.update_yaxes(title_text="PC2", row=2, col=1)
        fig.update_yaxes(title_text="UMAP2" if UMAP_AVAILABLE else "Distance", row=2, col=2)
        fig.update_yaxes(title_text="Distance", row=3, col=1)
        fig.update_yaxes(title_text="Abundance", row=3, col=2)

        return fig

    def generate_html_report(self, fig, title, output_file):
        """Generate HTML report with MetaPilot styling."""
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>{title}</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>
                body {{
                    background-color: {self.colors['background']};
                    color: {self.colors['text']};
                    font-family: Arial, sans-serif;
                    padding: 20px;
                }}
                .report-header {{
                    text-align: center;
                    margin-bottom: 30px;
                    padding: 20px;
                    background: linear-gradient(135deg, {self.colors['primary']}, {self.colors['secondary']});
                    color: white;
                    border-radius: 10px;
                }}
                .plot-container {{
                    background-color: white;
                    border-radius: 10px;
                    padding: 20px;
                    box-shadow: 0 4px 6px rgba(0,0,0,0.1);
                }}
                .info-box {{
                    background-color: #e3f2fd;
                    border-left: 4px solid {self.colors['primary']};
                    padding: 15px;
                    margin: 20px 0;
                    border-radius: 5px;
                }}
            </style>
        </head>
        <body>
            <div class="report-header">
                <h1>MetaPilot</h1>
                <h2>{title}</h2>
                <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            <div class="info-box">
                <strong>Analysis Features:</strong>
                <ul>
                    <li>✓ Metadata-aware grouping and coloring</li>
                    <li>✓ Advanced dimensionality reduction (PCA, UMAP)</li>
                    <li>✓ Hierarchical clustering analysis</li>
                    <li>✓ Interactive visualizations with zoom and hover</li>
                </ul>
            </div>
            <div class="plot-container">
                {fig.to_html(full_html=False, include_plotlyjs='cdn')}
            </div>
        </body>
        </html>
        """
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

# Backward compatibility
ReportGenerator = MetaPilotReportGenerator 