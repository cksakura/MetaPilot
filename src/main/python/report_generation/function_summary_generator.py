import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from .base_generator import BaseReportGenerator

try:
    from umap import UMAP
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False

class FunctionSummaryGenerator(BaseReportGenerator):
    """Generator for function summary reports."""
    
    def create_function_summary(self, func_df, metadata=None):
        """Enhanced function summary with category analysis and dimensionality reduction."""
        # Validate data before processing
        if not self.validate_data(func_df, "function"):
            return self.create_empty_data_plot(
                "Function Summary Report",
                "No function data available for analysis.<br>Please check if the functional annotation completed successfully."
            )

        fig = make_subplots(
            rows=7, cols=1,
            subplot_titles=(
                'Functional Categories Distribution (Top 5)',
                'All Functional Categories',
                'Top 20 Functions by Abundance',
                'PCA Analysis of Functional Profiles',
                'UMAP Analysis of Functions' if UMAP_AVAILABLE else 'Function Clustering',
                'Sample Clustering (Functions)',
                'Pathway Enrichment Analysis'
            ),
            specs=[[{"type": "pie"}],
                   [{"type": "table"}],
                   [{"type": "bar"}],
                   [{"type": "scatter"}],
                   [{"type": "scatter"}],
                   [{"type": "scatter"}],
                   [{"type": "bar"}]]
        )

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(func_df, metadata)
        
        # Ensure all intensity columns are numeric to prevent dtype errors
        for col in intensity_cols:
            if col in func_df.columns:
                func_df[col] = pd.to_numeric(func_df[col], errors='coerce').fillna(0)
        
        # Find function columns
        func_cols = ['Function', 'function', 'Description', 'Name', 'Gene_Name']
        func_col = next((col for col in func_cols if col in func_df.columns), None)
        
        cat_cols = ['NOG category', 'NOG_Category', 'NOG_category', 'nog_category', 'Category', 'category', 'Class', 'Type', 'COG category', 'COG_Category', 'KEGG_Category']
        cat_col = next((col for col in cat_cols if col in func_df.columns), None)

        if len(intensity_cols) > 0 and func_col:
            # Calculate total abundance
            func_df['Total_Abundance'] = func_df[intensity_cols].sum(axis=1)
            
            # NOG Category composition pie chart
            if cat_col:
                # Calculate total abundance per category
                cat_abundance = func_df.groupby(cat_col)['Total_Abundance'].sum()
                
                # Calculate percentage composition
                total_abundance = cat_abundance.sum()
                cat_percentages = (cat_abundance / total_abundance) * 100
                
                # Sort by abundance for better visualization
                cat_percentages = cat_percentages.sort_values(ascending=False)
                
                # Create pie chart with top 5 categories + Others
                top_5_categories = cat_percentages.head(5)
                others_percentage = cat_percentages.iloc[5:].sum() if len(cat_percentages) > 5 else 0
                
                # Prepare pie chart data
                if others_percentage > 0:
                    pie_labels = list(top_5_categories.index) + ['Others']
                    pie_values = list(top_5_categories.values) + [others_percentage]
                else:
                    pie_labels = list(top_5_categories.index)
                    pie_values = list(top_5_categories.values)
                
                # Use Quadram Institute colors for the pie chart
                colors = [self.colors['primary'], self.colors['secondary'], self.colors['tertiary'], 
                         self.colors['light_blue'], self.colors['light_purple'], '#cccccc']  # Grey for Others
                
                fig.add_trace(go.Pie(
                    labels=pie_labels,
                    values=pie_values,
                    name="Top Categories",
                    marker=dict(colors=colors[:len(pie_labels)]),
                    textinfo='label+percent',
                    textposition='auto',
                    pull=[0.1 if label != 'Others' else 0 for label in pie_labels],  # Slightly separate main categories
                    hovertemplate='<b>%{label}</b><br>Percentage: %{percent}<br>Value: %{value:.1f}<extra></extra>'
                ), row=1, col=1)
                
                # Add comprehensive table with ALL categories
                cat_table_data = pd.DataFrame({
                    'Category': cat_percentages.index,
                    'Abundance': cat_abundance[cat_percentages.index].values,
                    'Percentage': cat_percentages.values
                }).round(2)
                
                fig.add_trace(go.Table(
                    header=dict(
                        values=['<b>Functional Category</b>', '<b>Total Abundance</b>', '<b>Percentage (%)</b>'],
                        fill_color=self.colors['primary'],
                        font=dict(color='white', size=12),
                        align='left'
                    ),
                    cells=dict(
                        values=[cat_table_data['Category'], 
                               cat_table_data['Abundance'], 
                               cat_table_data['Percentage']],
                        fill_color=[['white', '#f0f8ff'] * len(cat_table_data)],
                        font=dict(color='black', size=11),
                        align='left',
                        height=25
                    )
                ), row=2, col=1)
            
            # Top functions bar chart
            top_20_funcs = func_df.nlargest(20, 'Total_Abundance')
            
            fig.add_trace(go.Bar(
                x=top_20_funcs['Total_Abundance'],
                y=top_20_funcs[func_col],
                orientation='h',
                name='Top Functions',
                marker_color=self.colors['tertiary']
            ), row=3, col=1)
            
            # PCA analysis of functional profiles
            if len(intensity_cols) >= 2:
                try:
                    # Prepare matrices
                    scaled_data, sample_feature = self.prepare_sample_feature_matrices(
                        func_df, intensity_cols, log_transform=True, pseudocount=1.0, standardize=True
                    )
                    if scaled_data is not None:
                        # PCA
                        _, pca_result = self.compute_pca(scaled_data, n_components=2)
                        if pca_result is not None:
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
                                    ), row=4, col=1)
                            else:
                                fig.add_trace(go.Scatter(
                                    x=pca_result[:, 0],
                                    y=pca_result[:, 1] if pca_result.shape[1] > 1 else np.zeros(len(pca_result)),
                                    mode='markers+text', text=sample_names, textposition="top center",
                                    marker=dict(size=12, color=self.colors['primary'], line=dict(width=2, color='white')),
                                    name='Functional PCA'
                                ), row=4, col=1)

                        # UMAP for functions (if available)
                        if scaled_data.shape[0] > 3:
                            _, umap_result = self.compute_umap(scaled_data, n_components=2, random_state=42)
                            if umap_result is not None:
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
                                        ), row=5, col=1)
                                else:
                                    fig.add_trace(go.Scatter(
                                        x=umap_result[:, 0], y=umap_result[:, 1],
                                        mode='markers+text', text=sample_names, textposition="top center",
                                        marker=dict(size=12, color=self.colors['success'], line=dict(width=2, color='white')),
                                        name='Functional UMAP'
                                    ), row=5, col=1)

                        # Hierarchical clustering
                        icoord, dcoord = self.compute_dendrogram_coords(scaled_data)
                        if icoord is not None and dcoord is not None:
                            for i in range(len(icoord)):
                                fig.add_trace(go.Scatter(
                                    x=[icoord[i][0], icoord[i][1], icoord[i][2], icoord[i][3]],
                                    y=[dcoord[i][0], dcoord[i][1], dcoord[i][2], dcoord[i][3]],
                                    mode='lines',                                 line=dict(color=self.colors['success'], width=2),
                                showlegend=False
                            ), row=6, col=1)
                            
                except Exception as e:
                    fig.add_annotation(text=f"Functional analysis failed: {str(e)[:40]}...",
                                     xref="x3", yref="y3", x=0.5, y=0.5, showarrow=False, row=3, col=1)
            
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
                            ), row=7, col=1)
                            
                except Exception as e:
                    fig.add_annotation(text=f"Group comparison failed: {str(e)[:30]}...",
                                     xref="x6", yref="y6", x=0.5, y=0.5, showarrow=False, row=4, col=2)

        fig.update_layout(height=2800, title_text="MetaPilot Function Summary Report")
        
        # Update axis labels
        fig.update_xaxes(title_text="Abundance", row=3, col=1)
        fig.update_xaxes(title_text="PC1", row=4, col=1)
        fig.update_xaxes(title_text="UMAP1" if UMAP_AVAILABLE else "Samples", row=5, col=1)
        fig.update_xaxes(title_text="Groups", row=7, col=1)
        
        fig.update_yaxes(title_text="Functions", row=3, col=1)
        fig.update_yaxes(title_text="PC2", row=4, col=1)
        fig.update_yaxes(title_text="UMAP2" if UMAP_AVAILABLE else "Distance", row=5, col=1)
        fig.update_yaxes(title_text="Distance", row=6, col=1)
        fig.update_yaxes(title_text="Abundance", row=7, col=1)

        return fig
