import os
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from datetime import datetime
from typing import Dict, List, Optional, Tuple

class BaseReportGenerator:
    """Base class for MetaPilot report generators with common functionality."""
    
    def __init__(self, output_dir="reports"):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Base color scheme (kept for backward compatibility)
        self.colors = {
            'primary': '#B5BD00',
            'secondary': '#8ED1FC',
            'tertiary': '#000000',
            'success': '#B5BD00',
            'background': '#f8f9fa',
            'text': '#000000',
            'light_blue': '#8ED1FC',
            'light_purple': '#B5BD00',
            'light_orange': '#8ED1FC'
        }

        # Color-blind friendly palette (Okabe–Ito)
        self.palette_okabe_ito = [
            '#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442',
            '#0072B2', '#D55E00', '#CC79A7', '#999999'
        ]

        # Extended group palette prioritizing Okabe–Ito, falling back to Quadram colors
        self.group_palette = (
            self.palette_okabe_ito
            + ['#B5BD00', '#8ED1FC', '#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499']
        )

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
        # Detect intensity/sample columns across common tools
        intensity_cols = self._detect_intensity_columns(df)

        # Derive sample names by stripping known prefixes
        sample_names = [
            (
                col.replace('LFQ intensity ', '')
                   .replace('iBAQ ', '')
                   .replace('Intensity ', '')
                   .replace('Intensity', '')
            )
            for col in intensity_cols
        ]
        
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
                    unique_groups = list(dict.fromkeys(groups))
                    group_colors = self._assign_group_colors(unique_groups)
            except Exception as e:
                print(f"Warning: Could not extract group information: {e}")
        
        return intensity_cols, sample_names, groups, group_colors

    def _detect_intensity_columns(self, df) -> List[str]:
        """Detect sample intensity columns across common proteomics outputs."""
        candidates: List[str] = []
        patterns = [
            'LFQ intensity ', 'iBAQ ', 'Intensity ',
        ]
        # DIA-NN/other tools may produce columns like "Intensity RawFile"; include generic contains rule
        for col in df.columns:
            if col == 'Intensity':
                continue
            if any(col.startswith(p) for p in patterns) or (('Intensity' in col) and col != 'Intensity'):
                candidates.append(col)
        return candidates

    def _assign_group_colors(self, unique_groups: List[str]) -> Dict[str, str]:
        """Assign consistent, distinct colors to groups using a shared palette."""
        color_map: Dict[str, str] = {}
        for i, group in enumerate(unique_groups):
            color_map[group] = self.group_palette[i % len(self.group_palette)]
        return color_map

    def prepare_sample_feature_matrices(
        self,
        df: pd.DataFrame,
        intensity_cols: List[str],
        log_transform: bool = True,
        pseudocount: float = 1.0,
        standardize: bool = True,
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """Prepare (samples x features) matrices: scaled_data and log_data.

        - Filters out features with zero total intensity across samples
        - Applies log transform with pseudo-count if requested
        - Standardizes features per sample matrix if requested
        """
        try:
            if len(intensity_cols) == 0:
                return None, None
            intensity_data = df[intensity_cols].fillna(0)
            intensity_data = intensity_data.loc[intensity_data.sum(axis=1) > 0]
            if len(intensity_data) == 0:
                return None, None
            if log_transform:
                log_data = np.log10(intensity_data + pseudocount)
            else:
                log_data = intensity_data.to_numpy(dtype=float)
            # Samples x features
            sample_feature = log_data.T
            if standardize:
                from sklearn.preprocessing import StandardScaler
                scaler = StandardScaler()
                scaled = scaler.fit_transform(sample_feature)
            else:
                scaled = sample_feature
            return scaled, sample_feature
        except Exception as e:
            print(f"Warning: Matrix preparation failed: {e}")
            return None, None

    def compute_pca(self, data: np.ndarray, n_components: int = 2):
        """Compute PCA embedding for samples x features matrix."""
        try:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=min(n_components, data.shape[1]))
            embedding = pca.fit_transform(data)
            return pca, embedding
        except Exception as e:
            print(f"Warning: PCA failed: {e}")
            return None, None

    def compute_umap(self, data: np.ndarray, n_components: int = 2, random_state: int = 42):
        """Compute UMAP embedding if umap-learn is available; else return (None, None)."""
        try:
            from umap import UMAP  # type: ignore
            reducer = UMAP(n_components=n_components, random_state=random_state)
            embedding = reducer.fit_transform(data)
            return reducer, embedding
        except Exception:
            return None, None

    def compute_dendrogram_coords(self, data: np.ndarray):
        """Compute dendrogram coordinates for samples using Ward linkage on Euclidean distances."""
        try:
            from scipy.cluster.hierarchy import dendrogram, linkage
            from scipy.spatial.distance import pdist
            distance_matrix = pdist(data, metric='euclidean')
            linkage_matrix = linkage(distance_matrix, method='ward')
            dendro = dendrogram(linkage_matrix, no_plot=True)
            return np.array(dendro['icoord']), np.array(dendro['dcoord'])
        except Exception as e:
            print(f"Warning: Dendrogram failed: {e}")
            return None, None

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
        """Generate a consistent HTML report using Plotly CDN and a light Bootstrap style."""
        try:
            fig.update_layout(title_text=title)
            fig_html = fig.to_html(full_html=False, include_plotlyjs='cdn', config={
                'displayModeBar': True,
                'displaylogo': False,
                'modeBarButtonsToRemove': ['pan2d', 'lasso2d', 'select2d']
            })

            full_html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{title}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ background-color: {self.colors['background']}; color: {self.colors['text']}; padding: 20px; }}
        .header {{ background-color: {self.colors['primary']}; color: #ffffff; padding: 20px; border-radius: 6px; margin-bottom: 20px; }}
        .content {{ background-color: #ffffff; border-radius: 10px; padding: 16px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
        .footer {{ text-align: center; margin-top: 16px; color: #666; font-size: 12px; }}
    </style>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    
    <!-- MetaPilot Report -->
</head>
<body>
    <div class="container">
        <div class="header">
            <h1 class="h3 mb-0">{title}</h1>
            <small>Generated by MetaPilot on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</small>
        </div>
        <div class="content">
            {fig_html}
        </div>
        <div class="footer">MetaPilot Report Generator</div>
    </div>
</body>
</html>
"""

            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(full_html)
            print(f"HTML report generated successfully: {output_file}")
            return True
        except Exception as e:
            print(f"Error generating HTML report: {e}")
            return False
