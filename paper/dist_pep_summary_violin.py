"""
Violin plot visualization for category-level peptide sharing (Plotly).

Input: TSV/CSV with columns:
Category,TotalDistinctPeptides,Peptides_in_1_genomes,Peptides_in_2_genomes,...

Output: an interactive violin plot (HTML). If kaleido is installed, can also write PNG.

Usage:
  python dist_pep_summary_violin.py dist_pep_summary_humangut.tsv out.html
  # or out.png (requires kaleido)
"""

import sys
import csv
from pathlib import Path
import pandas as pd
import plotly.express as px


def read_table(path: Path):
    with path.open("r", newline="") as f:
        # sniff delimiter (tab vs comma)
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.reader(f, dialect)
        header = next(reader)
        rows = list(reader)
    return header, rows


def load_data(path: Path):
    header, rows = read_table(path)
    if len(header) < 3:
        raise ValueError("Unexpected header format; need Category, TotalDistinctPeptides, share columns.")

    share_labels = header[2:]
    # Extract share number from column name
    share_bins = []
    for col in share_labels:
        parts = col.split("_")
        # expect ..._in_<N>_genomes
        try:
            share_bins.append(int(parts[-2]))
        except (IndexError, ValueError):
            raise ValueError(f"Cannot parse share bin from column name: {col}")

    categories = []
    data = []
    for row in rows:
        if not row or len(row) < len(header):
            continue
        cat = row[0]
        categories.append(cat)
        counts = [int(row[i + 2]) for i in range(len(share_bins))]
        data.append(counts)
    return categories, share_bins, data


def prepare_violin_data(categories, share_bins, data, max_share=15):
    """
    Prepare data for violin plot: expand counts into individual observations.
    For each category and sharing level, create that many data points.
    """
    plot_data = []
    plot_labels = []

    for cat_idx, (cat, counts) in enumerate(zip(categories, data)):
        for share_idx, share_level in enumerate(share_bins):
            if share_level > max_share:
                continue
            count = counts[share_idx]
            if count > 0:
                # Add 'count' number of observations at this sharing level
                plot_data.extend([share_level] * count)
                plot_labels.extend([cat] * count)

    return plot_data, plot_labels


def plot_violin(categories, share_bins, data, out_path: Path, max_share=15):
    """
    Create violin plot showing distribution of peptides across sharing levels for each category (Plotly).
    """
    # Filter to max_share and drop zero columns
    keep_idx = [i for i, s in enumerate(share_bins) if s <= max_share]
    share_bins_filtered = [share_bins[i] for i in keep_idx]
    data_filtered = [[row[i] for i in keep_idx] for row in data]

    # Drop zero-only columns
    keep_cols = []
    for idx in range(len(share_bins_filtered)):
        if any(d[idx] != 0 for d in data_filtered):
            keep_cols.append(idx)

    if not keep_cols:
        raise ValueError("All share columns are zero; nothing to plot.")

    share_bins_filtered = [share_bins_filtered[i] for i in keep_cols]
    data_filtered = [[row[i] for i in keep_cols] for row in data_filtered]

    # Prepare data for violin plot
    plot_data, plot_labels = prepare_violin_data(categories, share_bins_filtered, data_filtered, max_share)

    if not plot_data:
        raise ValueError("No data to plot after filtering.")

    df = pd.DataFrame({
        "Sharing_Level": plot_data,
        "Category": plot_labels
    })

    fig = px.violin(
        df,
        x="Category",
        y="Sharing_Level",
        color="Category",
        box=True,
        points=False,
        template="plotly_white",
        title="Distribution of peptide sharing levels across marker categories",
    )
    fig.update_yaxes(title="Found in N genomes")
    fig.update_xaxes(title="Category")

    ext = out_path.suffix.lower()
    if ext == ".html":
        fig.write_html(out_path)
    else:
        # requires kaleido for static export
        fig.write_image(out_path)


def main():
    if len(sys.argv) != 3:
        print("Usage: python dist_pep_summary_violin.py <input_tsv> <output_png>")
        sys.exit(1)

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    categories, share_bins, data = load_data(in_path)
    plot_violin(categories, share_bins, data, out_path, max_share=15)
    print(f"Saved violin plot to {out_path}")


if __name__ == "__main__":
    main()

