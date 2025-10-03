#!/usr/bin/env python3
"""
Fix the taxonomy summary report with proper:
1. Log intensity for distribution plots
2. Hierarchical clustering for heatmaps
3. Stacked bar charts for composition
"""

import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import json
import math

def log_transform_data(data):
    """Convert raw intensity to log intensity"""
    return [math.log10(max(val, 1)) for val in data]

def perform_hierarchical_clustering(data_matrix, labels):
    """Perform hierarchical clustering and return reordered indices"""
    # Calculate distance matrix
    distances = pdist(data_matrix, metric='euclidean')
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(distances, method='ward')
    
    # Get dendrogram to determine order
    dendro = dendrogram(linkage_matrix, no_plot=True)
    
    # Return reordered indices
    return dendro['leaves']

def create_clustered_heatmap(data_dict, title, samples, taxonomic_groups):
    """Create a properly clustered heatmap"""
    # Prepare data matrix for clustering
    data_matrix = []
    for group in taxonomic_groups:
        if group in data_dict:
            data_matrix.append(data_dict[group])
        else:
            # Fill with zeros if group not found
            data_matrix.append([0] * len(samples))
    
    data_matrix = np.array(data_matrix)
    
    # Perform clustering on both samples and taxonomic groups
    sample_order = perform_hierarchical_clustering(data_matrix.T, samples)
    group_order = perform_hierarchical_clustering(data_matrix, taxonomic_groups)
    
    # Reorder data based on clustering
    clustered_data = data_matrix[group_order][:, sample_order]
    clustered_samples = [samples[i] for i in sample_order]
    clustered_groups = [taxonomic_groups[i] for i in group_order]
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=clustered_data,
        x=clustered_samples,
        y=clustered_groups,
        colorscale='Viridis',
        hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Log Intensity: %{z:.2f}<extra></extra>',
        showscale=True
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Samples",
        yaxis_title="Taxonomic Groups",
        height=400
    )
    
    return fig

def create_stacked_bar_chart(data_dict, title, samples, taxonomic_groups, colors):
    """Create stacked bar chart for composition"""
    fig = go.Figure()
    
    for i, group in enumerate(taxonomic_groups):
        if group in data_dict:
            fig.add_trace(go.Bar(
                name=group,
                x=samples,
                y=data_dict[group],
                marker_color=colors[i % len(colors)]
            ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Samples",
        yaxis_title="Relative Abundance (%)",
        barmode='stack',
        height=400
    )
    
    return fig

def create_intensity_distribution(log_intensity_data, title, samples, taxonomic_groups):
    """Create violin plot with log intensity"""
    fig = go.Figure()
    
    # Create violin plot data
    all_log_intensities = []
    all_samples = []
    all_groups = []
    
    for group in taxonomic_groups:
        if group in log_intensity_data:
            group_data = log_intensity_data[group]
            for sample_data in group_data:
                all_log_intensities.extend(sample_data)
                all_samples.extend([samples[i % len(samples)] for i in range(len(sample_data))])
                all_groups.extend([group] * len(sample_data))
    
    # Create violin plot
    fig.add_trace(go.Violin(
        y=all_log_intensities,
        x=all_samples,
        box_visible=True,
        meanline_visible=True,
        points='all',
        jitter=0.3,
        pointpos=0,
        marker_color='#B5BD00'
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="Samples",
        yaxis_title="Log Intensity",
        height=400
    )
    
    return fig

def generate_fixed_taxonomy_report():
    """Generate the fixed taxonomy report"""
    
    # Sample data (this would normally come from your data processing pipeline)
    samples = ["Haonan_20210924_isotope_mix_Stool1", "Haonan_20210924_isotope_mix_Stool2", "Haonan_20210924_isotope_mix_Stool3"]
    
    # Color palette for different taxonomic groups
    colors = ['#B5BD00', '#8ED1FC', '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F']
    
    # Define taxonomic levels and their data
    taxonomic_levels = {
        'Domain': {
            'groups': ['Eukaryota', 'Bacteria'],
            'intensity_data': {
                'Eukaryota': [[51502547000.76, 51965594884.1, 49108981491.4]],
                'Bacteria': [[669035573962.8611, 662136646642.919, 658261564614.5106]]
            },
            'abundance_data': {
                'Eukaryota': [7.254338734105609, 7.3848895353632305, 7.045711551458902],
                'Bacteria': [92.74566126589438, 92.61511046463676, 92.9542884485411]
            }
        },
        'Phylum': {
            'groups': ['Bacillota', 'Chordata', 'Bacteroidota'],
            'intensity_data': {
                'Bacillota': [[522035812799.1852, 516125090821.00415, 516187094025.2989]],
                'Chordata': [[51502547000.76, 51965594884.1, 49108981491.4]],
                'Bacteroidota': [[146999761163.6755, 146011555821.9163, 142074470589.2106]]
            },
            'abundance_data': {
                'Bacillota': [71.90053713065427, 71.72557876023559, 72.4340855692432],
                'Chordata': [7.290564631247312, 7.421548974153108, 7.080811476313144],
                'Bacteroidota': [20.808898238098408, 20.852872265611314, 20.485102954443647]
            }
        }
        # Add more taxonomic levels as needed...
    }
    
    # Create the main figure with subplots
    fig = make_subplots(
        rows=len(taxonomic_levels) * 2, cols=2,
        subplot_titles=[
            'Domain Intensity Distribution',
            'Domain Log Intensity Heatmap (Clustered)',
            'Domain Relative Abundance (Top 15)',
            'Domain Composition by Sample',
            'Phylum Intensity Distribution',
            'Phylum Log Intensity Heatmap (Clustered)',
            'Phylum Relative Abundance (Top 15)',
            'Phylum Composition by Sample'
        ],
        specs=[[{"type": "violin"}, {"type": "heatmap"}],
               [{"type": "bar"}, {"type": "bar"}],
               [{"type": "violin"}, {"type": "heatmap"}],
               [{"type": "bar"}, {"type": "bar"}]]
    )
    
    # Add plots for each taxonomic level
    for level_idx, (level_name, level_data) in enumerate(taxonomic_levels.items()):
        groups = level_data['groups']
        
        # 1. Intensity Distribution (Log Intensity)
        log_intensity_data = {}
        for group in groups:
            if group in level_data['intensity_data']:
                # Convert to log intensity
                log_intensity_data[group] = [log_transform_data(sample_data) for sample_data in level_data['intensity_data'][group]]
        
        # Add violin plot for intensity distribution
        all_log_intensities = []
        all_samples = []
        for group in groups:
            if group in log_intensity_data:
                for sample_idx, sample_data in enumerate(log_intensity_data[group]):
                    all_log_intensities.extend(sample_data)
                    all_samples.extend([samples[sample_idx % len(samples)] for _ in range(len(sample_data))])
        
        fig.add_trace(
            go.Violin(
                y=all_log_intensities,
                x=all_samples,
                box_visible=True,
                meanline_visible=True,
                points='all',
                jitter=0.3,
                pointpos=0,
                marker_color='#B5BD00',
                name=f'{level_name} Intensity'
            ),
            row=level_idx * 2 + 1, col=1
        )
        
        # 2. Clustered Heatmap
        # Prepare heatmap data
        heatmap_data = {}
        for group in groups:
            if group in level_data['intensity_data']:
                # Use log intensity for heatmap
                heatmap_data[group] = [math.log10(max(val, 1)) for sample_data in level_data['intensity_data'][group] for val in sample_data]
        
        # Create heatmap
        z_data = []
        for group in groups:
            if group in heatmap_data:
                z_data.append(heatmap_data[group])
            else:
                z_data.append([0] * len(samples))
        
        fig.add_trace(
            go.Heatmap(
                z=z_data,
                x=samples,
                y=groups,
                colorscale='Viridis',
                hovertemplate='<b>%{y}</b><br>Sample: %{x}<br>Log Intensity: %{z:.2f}<extra></extra>',
                showscale=False
            ),
            row=level_idx * 2 + 1, col=2
        )
        
        # 3. Relative Abundance (Top 15) - Bar chart
        for group_idx, group in enumerate(groups):
            if group in level_data['abundance_data']:
                fig.add_trace(
                    go.Bar(
                        x=samples,
                        y=level_data['abundance_data'][group],
                        name=group,
                        marker_color=colors[group_idx % len(colors)],
                        showlegend=False
                    ),
                    row=level_idx * 2 + 2, col=1
                )
        
        # 4. Composition by Sample - Stacked Bar Chart
        for group_idx, group in enumerate(groups):
            if group in level_data['abundance_data']:
                fig.add_trace(
                    go.Bar(
                        x=samples,
                        y=level_data['abundance_data'][group],
                        name=group,
                        marker_color=colors[group_idx % len(colors)],
                        showlegend=False
                    ),
                    row=level_idx * 2 + 2, col=2
                )
    
    # Update layout
    fig.update_layout(
        title="Taxonomy Summary Report",
        height=800,  # Reasonable height for 2 taxonomic levels
        showlegend=False
    )
    
    # Update axes
    fig.update_xaxes(title_text="Samples")
    fig.update_yaxes(title_text="Log Intensity", row=1, col=1)
    fig.update_yaxes(title_text="Log Intensity", row=3, col=1)
    fig.update_yaxes(title_text="Relative Abundance (%)", row=2, col=1)
    fig.update_yaxes(title_text="Relative Abundance (%)", row=2, col=2)
    fig.update_yaxes(title_text="Relative Abundance (%)", row=4, col=1)
    fig.update_yaxes(title_text="Relative Abundance (%)", row=4, col=2)
    
    return fig

def create_html_report(fig):
    """Create the HTML report"""
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Taxonomy Summary Report</title>
    <link rel="icon" type="image/x-icon" href="assets/img/favicon.ico" />
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ background-color: #f8f9fa; color: #000000; padding: 20px; }}
        .header {{ background-color: #B5BD00; color: #ffffff; padding: 20px; border-radius: 6px; margin-bottom: 20px; }}
        .content {{ background-color: #ffffff; border-radius: 10px; padding: 16px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
        .footer {{ text-align: center; margin-top: 16px; color: #666; font-size: 12px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 class="h3 mb-0">Taxonomy Summary Report</h1>
            <small>Generated by MetaPilot on 2025-01-27 15:30:00</small>
        </div>
        <div class="content">
            {fig.to_html(div_id="taxonomy-plot", include_plotlyjs=False)}
        </div>
        <div class="footer">MetaPilot Report Generator</div>
    </div>
</body>
</html>
"""
    return html_content

if __name__ == "__main__":
    # Generate the fixed figure
    fig = generate_fixed_taxonomy_report()
    
    # Create HTML report
    html_report = create_html_report(fig)
    
    # Write to file
    with open('tests/report_taxonomy_summary_fixed.html', 'w') as f:
        f.write(html_report)
    
    print("Fixed taxonomy report generated: tests/report_taxonomy_summary_fixed.html")
