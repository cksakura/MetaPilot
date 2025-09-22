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
        """Enhanced taxonomy summary with rank-specific analysis for all taxonomic ranks."""
        # Validate data before processing
        if not self.validate_data(taxa_df, "taxonomy"):
            return self.create_empty_data_plot(
                "Taxonomy Summary Report",
                "No taxonomy data available for analysis.<br>Please check if the taxonomic analysis completed successfully."
            )

        # Identify taxonomic columns
        rank_cols = ["Rank", "rank", "Level", "Taxonomic_Rank"]
        rank_col = next((c for c in rank_cols if c in taxa_df.columns), None)
        taxon_cols = ["Taxon", "taxon", "Name", "Taxa", "Species"]
        taxon_col = next((c for c in taxon_cols if c in taxa_df.columns), None)

        intensity_cols, sample_names, groups, group_colors = self.extract_sample_info(taxa_df, metadata)

        # Ensure intensity columns are numeric
        for c in intensity_cols:
            if c in taxa_df.columns:
                taxa_df[c] = pd.to_numeric(taxa_df[c], errors="coerce").fillna(0)

        # Fallback simple layout if key columns are missing
        if rank_col is None or taxon_col is None or len(intensity_cols) == 0:
            fig = make_subplots(
                rows=1,
                cols=2,
                subplot_titles=("Taxonomic Composition (Treemap)", "Top 20 Taxa by Abundance"),
                specs=[[{"type": "treemap"}, {"type": "bar"}]],
            )
            if len(intensity_cols) > 0 and taxon_col:
                tmp = taxa_df.copy()
                tmp["Total_Abundance"] = tmp[intensity_cols].sum(axis=1)
                top_20 = tmp.nlargest(20, "Total_Abundance")
                fig.add_trace(
                    go.Treemap(
                        labels=top_20[taxon_col].head(10),
                        values=top_20["Total_Abundance"].head(10),
                        parents=[""] * min(10, len(top_20)),
                        textinfo="label+percent entry",
                        marker=dict(colorscale="Viridis"),
                        name="Taxonomic Composition",
                    ),
                    row=1,
                    col=1,
                )
                fig.add_trace(
                    go.Bar(
                        x=top_20["Total_Abundance"],
                        y=top_20[taxon_col],
                        orientation="h",
                        name="Top Taxa",
                        marker_color=self.colors["tertiary"],
                    ),
                    row=1,
                    col=2,
                )
            fig.update_layout(height=1600, title_text="MetaPilot Taxonomy Summary Report")
            return fig

        # Compute total abundance and sort ranks by hierarchy
        taxa_df = taxa_df.copy()
        taxa_df["Total_Abundance"] = taxa_df[intensity_cols].sum(axis=1)
        available_ranks = pd.Series(taxa_df[rank_col].dropna().unique())
        rank_hierarchy = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
        rank_order = {}
        for i, r in enumerate(rank_hierarchy):
            for ar in available_ranks:
                if r in str(ar).lower():
                    rank_order[ar] = i
        sorted_ranks = sorted(available_ranks, key=lambda x: (rank_order.get(x, 999), str(x)))

        # Build subplot grid: Single column layout with double height
        num_ranks = min(len(sorted_ranks), 4)
        total_rows = 1 + num_ranks + 2
        subplot_titles = [
            "Taxonomic Hierarchy (Multiple Ring)",
            "Top 20 Taxa by Abundance",
        ]
        specs = [[{"type": "sunburst"}], [{"type": "bar"}]]

        for r in sorted_ranks[:num_ranks]:
            subplot_titles.extend([f"{r} Composition (Treemap)", f"Top 20 {r} by Abundance"])
            specs.append([{"type": "treemap"}])
            specs.append([{"type": "bar"}])

        subplot_titles.extend(["Taxa Count by Rank", "PCA of Taxonomic Profiles"])
        specs.append([{"type": "bar"}])
        specs.append([{"type": "scatter"}])

        subplot_titles.extend(["Diversity Indices by Sample"])
        specs.append([{"type": "bar"}])

        fig = make_subplots(rows=total_rows, cols=1, subplot_titles=subplot_titles, specs=specs)

        current_row = 1
        # Sunburst
        sunburst_data = self._create_sunburst_data(taxa_df, rank_col, taxon_col, sorted_ranks)
        if sunburst_data:
            fig.add_trace(
                go.Sunburst(
                    labels=sunburst_data["labels"],
                    parents=sunburst_data["parents"],
                    values=sunburst_data["values"],
                    branchvalues="total",
                    textinfo="label+percent entry",
                    textfont=dict(size=10),
                    marker=dict(colorscale="Viridis", colors=sunburst_data["colors"], line=dict(width=1, color="white")),
                    name="Taxonomic Hierarchy",
                ),
                row=current_row,
                col=1,
            )

        # Top-20 overall
        top_20_taxa = taxa_df.nlargest(20, "Total_Abundance")
        fig.add_trace(
            go.Bar(
                x=top_20_taxa["Total_Abundance"],
                y=top_20_taxa[taxon_col],
                orientation="h",
                name="Top Taxa Overall",
                marker_color=self.colors["tertiary"],
            ),
            row=current_row,
            col=2,
        )
        current_row += 1

        # Rank-specific rows
        for r in sorted_ranks[:num_ranks]:
            rdata = taxa_df[taxa_df[rank_col] == r]
            if len(rdata) == 0:
                current_row += 1
                continue
            rtotal = rdata["Total_Abundance"].sum()
            if rtotal <= 0:
                current_row += 1
                continue
            rdata = rdata.copy()
            rdata["Rank_Relative_Abundance"] = (rdata["Total_Abundance"] / rtotal) * 100
            top_20_rank = rdata.nlargest(20, "Total_Abundance")
            top_10_rank = top_20_rank.head(10)
            fig.add_trace(
                go.Treemap(
                    labels=top_10_rank[taxon_col],
                    values=top_10_rank["Rank_Relative_Abundance"],
                    parents=[""] * len(top_10_rank),
                    textinfo="label+percent entry",
                    marker=dict(colorscale="Viridis"),
                    name=f"{r} Composition",
                ),
                row=current_row,
                col=1,
            )
            fig.add_trace(
                go.Bar(
                    x=top_20_rank["Rank_Relative_Abundance"],
                    y=top_20_rank[taxon_col],
                    orientation="h",
                    name=f"Top {r}",
                    marker_color=self.colors["tertiary"],
                ),
                row=current_row,
                col=2,
            )
            current_row += 1

        # Summary: taxa count
        rank_counts = taxa_df[rank_col].value_counts()
        rank_counts_sorted = rank_counts.reindex(sorted_ranks)
        fig.add_trace(
            go.Bar(
                x=list(rank_counts_sorted.index.astype(str)),
                y=rank_counts_sorted.values,
                name="Taxa Count by Rank",
                marker_color=self.colors["secondary"],
            ),
            row=current_row,
            col=1,
        )

        # Summary: PCA of profiles
        if len(intensity_cols) >= 2:
            try:
                profile_data = taxa_df[intensity_cols].fillna(0)
                profile_data = profile_data.loc[profile_data.sum(axis=1) > 0]
                if len(profile_data) > 0:
                    rel_abundance = profile_data.div(profile_data.sum(axis=0), axis=1) * 100
                    log_data = np.log10(rel_abundance + 0.01)
                    scaler = StandardScaler()
                    scaled = scaler.fit_transform(log_data.T)
                    pca = PCA(n_components=min(2, len(intensity_cols)))
                    pcr = pca.fit_transform(scaled)
                    if groups and group_colors:
                        for g in sorted(set(groups)):
                            idx = [i for i, gg in enumerate(groups) if gg == g]
                            fig.add_trace(
                                go.Scatter(
                                    x=[pcr[i, 0] for i in idx],
                                    y=[pcr[i, 1] for i in idx] if pcr.shape[1] > 1 else [0] * len(idx),
                                    mode="markers+text",
                                    text=[sample_names[i] for i in idx],
                                    textposition="top center",
                                    marker=dict(size=12, color=group_colors[g], line=dict(width=2, color="white")),
                                    name=f"Group {g}",
                                    showlegend=True,
                                ),
                                row=current_row,
                                col=2,
                            )
                    else:
                        fig.add_trace(
                            go.Scatter(
                                x=pcr[:, 0],
                                y=pcr[:, 1] if pcr.shape[1] > 1 else np.zeros(len(pcr)),
                                mode="markers+text",
                                text=sample_names,
                                textposition="top center",
                                marker=dict(size=12, color=self.colors["primary"], line=dict(width=2, color="white")),
                                name="Taxonomic PCA",
                            ),
                            row=current_row,
                            col=2,
                        )
            except Exception as e:
                fig.add_annotation(
                    text=f"Taxonomic PCA failed: {str(e)[:60]}...",
                    xref="paper",
                    yref="paper",
                    x=0.75,
                    y=1 - (current_row - 1) / total_rows,
                    showarrow=False,
                )

        current_row += 1

        # Summary: Diversity indices (Shannon)
        if len(intensity_cols) > 0:
            shannon = []
            for c in intensity_cols:
                s = taxa_df[c]
                s = s[s > 0]
                if len(s) > 0:
                    rel = s / s.sum()
                    shannon.append(float(-np.sum(rel * np.log(rel))))
                else:
                    shannon.append(0.0)
            colors = [group_colors[groups[i]] if groups and group_colors else self.colors["light_blue"] for i in range(len(sample_names))]
            fig.add_trace(
                go.Bar(x=sample_names, y=shannon, name="Shannon Diversity", marker_color=colors),
                row=current_row,
                col=1,
            )

        # Layout sizing and labels
        base_height = 500
        rank_height = 560  # Double the height for single column layout
        total_height = base_height + num_ranks * rank_height
        fig.update_layout(height=total_height, title_text="MetaPilot Taxonomy Summary Report - All Taxonomic Ranks")
        fig.update_xaxes(title_text="Taxonomic Rank", row=(1 + num_ranks), col=1)
        fig.update_xaxes(title_text="PC1", row=(1 + num_ranks + 1), col=1)
        fig.update_yaxes(title_text="Number of Taxa", row=(1 + num_ranks), col=1)
        return fig

    def _create_sunburst_data(self, taxa_df, rank_col, taxon_col, sorted_ranks):
        """Create hierarchical data structure for sunburst chart showing taxonomic relationships."""
        try:
            # Define taxonomic hierarchy order (from highest to lowest rank)
            rank_hierarchy = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            
            # Filter to only include ranks that exist in our data
            available_hierarchy = []
            for rank in rank_hierarchy:
                for available_rank in sorted_ranks:
                    if rank in str(available_rank).lower():
                        available_hierarchy.append(available_rank)
                        break
            
            if len(available_hierarchy) < 2:
                return None
            
            # Create hierarchical structure
            labels = []
            parents = []
            values = []
            colors = []
            
            # Color palette for different taxonomic levels
            color_palette = [
                '#B5BD00', '#8ED1FC', '#000000', '#B5BD00', '#8ED1FC',
                '#000000', '#B5BD00', '#8ED1FC', '#000000', '#B5BD00'
            ]
            
            # Track parent-child relationships
            parent_map = {}  # Maps taxon name to its parent
            rank_data = {}   # Stores data for each rank
            
            # Process each rank level
            for i, rank in enumerate(available_hierarchy):
                rank_data[rank] = taxa_df[taxa_df[rank_col] == rank].copy()
                if len(rank_data[rank]) == 0:
                    continue
                
                # Sort by abundance and take top taxa for this rank
                rank_data[rank] = rank_data[rank].sort_values('Total_Abundance', ascending=False)
                
                # Limit number of taxa per rank to avoid overcrowding
                max_taxa_per_rank = 15 if i == 0 else 10  # More for higher ranks
                rank_data[rank] = rank_data[rank].head(max_taxa_per_rank)
                
                # Add taxa for this rank
                for _, row in rank_data[rank].iterrows():
                    taxon_name = str(row[taxon_col])
                    abundance = row['Total_Abundance']
                    
                    # Skip if abundance is too low
                    if abundance < taxa_df['Total_Abundance'].max() * 0.001:  # 0.1% threshold
                        continue
                    
                    labels.append(taxon_name)
                    values.append(abundance)
                    
                    # Determine parent
                    if i == 0:  # Highest rank (root level)
                        parents.append("")
                        parent_map[taxon_name] = ""
                    else:
                        # Find parent from previous rank
                        parent_rank = available_hierarchy[i-1]
                        parent_found = False
                        
                        # Try to find a logical parent based on taxonomic relationships
                        for parent_taxon, parent_abundance in parent_map.items():
                            if parent_abundance != "" and parent_abundance in rank_data[parent_rank][taxon_col].values:
                                # Check if this taxon might belong to this parent
                                # This is a simplified approach - in reality, you'd need proper taxonomic lineage
                                if abundance < parent_abundance * 0.5:  # Child should be smaller than parent
                                    parents.append(parent_taxon)
                                    parent_map[taxon_name] = parent_taxon
                                    parent_found = True
                                    break
                
                        if not parent_found:
                            # Default to root if no parent found
                            parents.append("")
                            parent_map[taxon_name] = ""
                    
                    # Assign color based on rank level
                    color_index = i % len(color_palette)
                    colors.append(color_palette[color_index])
            
            # If we don't have enough data, return None
            if len(labels) < 3:
                return None
            
            return {
                'labels': labels,
                'parents': parents,
                'values': values,
                'colors': colors
            }
                            
        except Exception as e:
            print(f"Error creating sunburst data: {e}")
            return None
