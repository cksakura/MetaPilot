import json
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
from .base_generator import BaseReportGenerator
from typing import Dict, List, Optional, Tuple

class BubbleChartGenerator(BaseReportGenerator):
    """Generator for bubble chart reports showing hierarchical taxonomic abundance."""
    
    def create_bubble_chart_summary(self, tree_data, metadata=None, max_depth=4, min_abundance=0.001):
        """
        Create a hierarchical bubble chart from tree data.
        
        Args:
            tree_data: JSON string or list of dictionaries containing hierarchical tree data
            metadata: Optional metadata for sample grouping
            max_depth: Maximum depth to display (default: 4)
            min_abundance: Minimum relative abundance to display (default: 0.001)
        """
        try:
            # Parse tree data
            if isinstance(tree_data, str):
                tree_data = json.loads(tree_data)
            
            if not tree_data or len(tree_data) == 0:
                return self.create_empty_data_plot(
                    "Bubble Chart Report",
                    "No tree data available for bubble chart visualization."
                )
            
            # Convert to DataFrame for easier processing
            df = pd.DataFrame(tree_data)
            
            # Extract sample information
            intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(df, metadata)
            
            if len(intensity_cols) == 0:
                return self.create_empty_data_plot(
                    "Bubble Chart Report",
                    "No intensity data found in tree data."
                )
            
            # Process tree data for visualization
            processed_data = self._process_tree_data(df, intensity_cols, max_depth, min_abundance)
            
            if not processed_data:
                return self.create_empty_data_plot(
                    "Bubble Chart Report",
                    "No data meets the minimum abundance threshold."
                )
            
            # Create the MetaMap visualization
            fig = self._create_metamap_visualization(processed_data, sample_names, groups, group_colors)
            
            return fig
            
        except Exception as e:
            print(f"Error creating MetaMap: {e}")
            return self.create_empty_data_plot(
                "Bubble Chart Report",
                f"Error processing tree data: {str(e)}"
            )
    
    def _process_tree_data(self, df, intensity_cols, max_depth, min_abundance):
        """Process tree data for bubble chart visualization."""
        processed_data = []
        
        for _, row in df.iterrows():
            # Parse hierarchical ID
            hierarchy = row['id'].split('@')
            depth = len(hierarchy) - 1
            
            # Skip if too deep
            if depth > max_depth:
                continue
            
            # Get abundance values
            total_value = float(row['value'])
            
            # Calculate relative abundance
            if total_value <= 0:
                continue
            
            # Get sample intensities
            sample_intensities = []
            for col in intensity_cols:
                if col in row:
                    # Handle semicolon-separated values
                    intensity_str = str(row[col])
                    if ';' in intensity_str:
                        intensities = [float(x) for x in intensity_str.split(';')]
                        sample_intensities.extend(intensities)
                    else:
                        sample_intensities.append(float(intensity_str))
                else:
                    sample_intensities.append(0.0)
            
            # Calculate relative abundance for each sample
            total_sample_intensity = sum(sample_intensities)
            if total_sample_intensity > 0:
                relative_abundances = [intensity / total_sample_intensity for intensity in sample_intensities]
                max_relative_abundance = max(relative_abundances)
                
                # Skip if below minimum abundance threshold
                if max_relative_abundance < min_abundance:
                    continue
                
                # Create data entry
                entry = {
                    'id': row['id'],
                    'name': hierarchy[-1],  # Last element is the taxon name
                    'parent': '@'.join(hierarchy[:-1]) if len(hierarchy) > 1 else '',
                    'depth': depth,
                    'total_value': total_value,
                    'sample_intensities': sample_intensities,
                    'relative_abundances': relative_abundances,
                    'max_relative_abundance': max_relative_abundance,
                    'hierarchy': hierarchy
                }
                processed_data.append(entry)
        
        # Sort by total value (largest first)
        processed_data.sort(key=lambda x: x['total_value'], reverse=True)
        
        return processed_data
    
    def _create_metamap_visualization(self, processed_data, sample_names, groups, group_colors):
        """Create the bubble chart visualization."""
        
        # Create subplots for different views
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Bubble Chart - Total Abundance',
                'Bubble Chart - Sample Comparison',
                'Taxonomic Hierarchy Tree',
                'Abundance Distribution'
            ),
            specs=[[{"type": "treemap"}, {"type": "treemap"}],
                   [{"type": "sunburst"}, {"type": "bar"}]]
        )
        
        # 1. Total Abundance Treemap
        self._add_total_abundance_treemap(fig, processed_data, row=1, col=1)
        
        # 2. Sample Comparison Treemap (if multiple samples)
        if len(sample_names) > 1:
            self._add_sample_comparison_treemap(fig, processed_data, sample_names, groups, group_colors, row=1, col=2)
        else:
            # If only one sample, show the same as total abundance
            self._add_total_abundance_treemap(fig, processed_data, row=1, col=2)
        
        # 3. Taxonomic Hierarchy Sunburst
        self._add_hierarchy_sunburst(fig, processed_data, row=2, col=1)
        
        # 4. Abundance Distribution Bar Chart
        self._add_abundance_distribution(fig, processed_data, row=2, col=2)
        
        # Update layout
        fig.update_layout(
            title_text="Bubble Chart: Hierarchical Taxonomic Abundance Visualization",
            title_x=0.5,
            height=800,
            showlegend=True,
            font=dict(size=12)
        )
        
        return fig
    
    def _add_total_abundance_treemap(self, fig, processed_data, row, col):
        """Add total abundance treemap."""
        # Prepare data for treemap
        treemap_data = []
        for item in processed_data[:50]:  # Limit to top 50 for performance
            treemap_data.append({
                'ids': item['id'],
                'labels': item['name'],
                'parents': item['parent'],
                'values': item['total_value'],
                'text': f"{item['name']}<br>Total: {item['total_value']:,.0f}<br>Max Rel. Abundance: {item['max_relative_abundance']:.3f}",
                'max_relative_abundance': item['max_relative_abundance'],
                'textinfo': 'label+value'
            })
        
        if treemap_data:
            fig.add_trace(
                go.Treemap(
                    ids=[d['ids'] for d in treemap_data],
                    labels=[d['labels'] for d in treemap_data],
                    parents=[d['parents'] for d in treemap_data],
                    values=[d['values'] for d in treemap_data],
                    text=[d['text'] for d in treemap_data],
                    textinfo="label+value",
                    textfont=dict(size=10),
                    marker=dict(
                        colorscale='Viridis',
                        colorbar=dict(title="Total Abundance")
                    ),
                    hovertemplate="<b>%{label}</b><br>" +
                                "Total Abundance: %{value:,.0f}<br>" +
                                "Max Relative Abundance: %{customdata:.3f}<extra></extra>",
                    customdata=[d['max_relative_abundance'] for d in treemap_data]
                ),
                row=row, col=col
            )
    
    def _add_sample_comparison_treemap(self, fig, processed_data, sample_names, groups, group_colors, row, col):
        """Add sample comparison treemap."""
        # For sample comparison, we'll show the most abundant sample for each taxon
        treemap_data = []
        for item in processed_data[:30]:  # Limit for performance
            # Find the sample with highest abundance for this taxon
            max_idx = item['relative_abundances'].index(item['max_relative_abundance'])
            sample_name = sample_names[max_idx] if max_idx < len(sample_names) else f"Sample {max_idx + 1}"
            
            treemap_data.append({
                'ids': item['id'],
                'labels': item['name'],
                'parents': item['parent'],
                'values': item['sample_intensities'][max_idx],
                'text': f"{item['name']}<br>Sample: {sample_name}<br>Intensity: {item['sample_intensities'][max_idx]:,.0f}",
                'sample_name': sample_name
            })
        
        if treemap_data:
            fig.add_trace(
                go.Treemap(
                    ids=[d['ids'] for d in treemap_data],
                    labels=[d['labels'] for d in treemap_data],
                    parents=[d['parents'] for d in treemap_data],
                    values=[d['values'] for d in treemap_data],
                    text=[d['text'] for d in treemap_data],
                    textinfo="label+value",
                    textfont=dict(size=10),
                    marker=dict(
                        colorscale='Plasma',
                        colorbar=dict(title="Sample Intensity")
                    ),
                    hovertemplate="<b>%{label}</b><br>" +
                                "Sample: %{customdata}<br>" +
                                "Intensity: %{value:,.0f}<extra></extra>",
                    customdata=[d['sample_name'] for d in treemap_data]
                ),
                row=row, col=col
            )
    
    def _add_hierarchy_sunburst(self, fig, processed_data, row, col):
        """Add taxonomic hierarchy sunburst chart."""
        # Prepare data for sunburst
        sunburst_data = []
        for item in processed_data[:100]:  # Limit for performance
            sunburst_data.append({
                'ids': item['id'],
                'labels': item['name'],
                'parents': item['parent'],
                'values': item['total_value']
            })
        
        if sunburst_data:
            fig.add_trace(
                go.Sunburst(
                    ids=[d['ids'] for d in sunburst_data],
                    labels=[d['labels'] for d in sunburst_data],
                    parents=[d['parents'] for d in sunburst_data],
                    values=[d['values'] for d in sunburst_data],
                    branchvalues="total",
                    textinfo="label+value",
                    textfont=dict(size=10),
                    marker=dict(
                        colorscale='Rainbow',
                        colorbar=dict(title="Abundance")
                    ),
                    hovertemplate="<b>%{label}</b><br>" +
                                "Abundance: %{value:,.0f}<br>" +
                                "Path: %{id}<extra></extra>"
                ),
                row=row, col=col
            )
    
    def _add_abundance_distribution(self, fig, processed_data, row, col):
        """Add abundance distribution bar chart."""
        # Show top 20 most abundant taxa
        top_taxa = processed_data[:20]
        
        if top_taxa:
            fig.add_trace(
                go.Bar(
                    x=[item['name'] for item in top_taxa],
                    y=[item['total_value'] for item in top_taxa],
                    text=[f"{item['total_value']:,.0f}" for item in top_taxa],
                    textposition='auto',
                    marker=dict(
                        color='lightblue',
                        line=dict(color='darkblue', width=1)
                    ),
                    hovertemplate="<b>%{x}</b><br>" +
                                "Total Abundance: %{y:,.0f}<br>" +
                                "Max Rel. Abundance: %{customdata:.3f}<extra></extra>",
                    customdata=[item['max_relative_abundance'] for item in top_taxa]
                ),
                row=row, col=col
            )
            
            # Update x-axis for better readability
            fig.update_xaxes(tickangle=45, row=row, col=col)
            fig.update_yaxes(title_text="Total Abundance", row=row, col=col)

