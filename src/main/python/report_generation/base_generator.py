import os
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from datetime import datetime

class BaseReportGenerator:
    """Base class for MetaPilot report generators with common functionality."""
    
    def __init__(self, output_dir="reports"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Quadram Institute color scheme
        self.colors = {
            'primary': '#B5BD00',      # Buddha Gold
            'secondary': '#8ED1FC',    # Malibu
            'tertiary': '#000000',     # Black
            'success': '#B5BD00',      # Buddha Gold for success
            'background': '#f8f9fa',
            'text': '#000000',         # Black text
            'light_blue': '#8ED1FC',   # Malibu
            'light_purple': '#B5BD00', # Buddha Gold
            'light_orange': '#8ED1FC'  # Malibu
        }
        
        # Color palette for groups using Quadram Institute colors
        self.group_colors = [
            '#B5BD00', '#8ED1FC', '#000000', '#B5BD00', '#8ED1FC',
            '#000000', '#B5BD00', '#8ED1FC', '#000000', '#B5BD00'
        ]

    def validate_data(self, df, data_type="data"):
        """Validate that the dataframe contains actual data (not just headers)."""
        if df is None:
            print(f"WARNING: {data_type} dataframe is None")
            return False
        
        if df.empty:
            print(f"WARNING: {data_type} dataframe is empty (no rows)")
            return False
        
        # Check if all rows are empty (all NaN or empty strings)
        if df.dropna(how='all').empty:
            print(f"WARNING: {data_type} dataframe contains only empty rows")
            return False
        
        # Check if we have meaningful data (at least one non-empty cell)
        has_data = False
        for col in df.columns:
            if not df[col].dropna().empty:
                has_data = True
                break
        
        if not has_data:
            print(f"WARNING: {data_type} dataframe contains no meaningful data")
            return False
        
        print(f"INFO: {data_type} validation passed - {len(df)} rows, {len(df.columns)} columns")
        return True

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

    def create_empty_data_plot(self, title, message):
        """Create a placeholder plot for empty data."""
        fig = make_subplots(rows=1, cols=1, subplot_titles=(title,))
        fig.add_annotation(
            x=0.5, y=0.5,
            xref="paper", yref="paper",
            text=message,
            showarrow=False,
            font=dict(size=16, color="red"),
            align="center"
        )
        fig.update_layout(
            title_text=f"{title} - No Data Available",
            title_x=0.5,
            height=400,
            showlegend=False
        )
        return fig

    def generate_html_report(self, fig, title, output_file):
        """Generate HTML report from a Plotly figure."""
        try:
            # Update the figure layout with the title
            fig.update_layout(title_text=title)
            
            # Generate HTML content (just the div and script, not full document)
            html_content = fig.to_html(
                include_plotlyjs=True,
                div_id="plotly-div",
                full_html=False,  # This prevents the full HTML document structure
                config={
                    'displayModeBar': True,
                    'displaylogo': False,
                    'modeBarButtonsToRemove': ['pan2d', 'lasso2d', 'select2d']
                }
            )
            
            # Create a more complete HTML document
            full_html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f8f9fa;
        }}
        .header {{
            background-color: #B5BD00;
            color: white;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .content {{
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .footer {{
            text-align: center;
            margin-top: 20px;
            color: #666;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>{title}</h1>
        <p>Generated by MetaPilot on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    <div class="content">
        {html_content}
    </div>
    <div class="footer">
        <p>MetaPilot Report Generator - Quadram Institute</p>
    </div>
</body>
</html>
"""
            
            # Write to file
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(full_html)
            
            print(f"HTML report generated successfully: {output_file}")
            return True
            
        except Exception as e:
            print(f"Error generating HTML report: {e}")
            return False
