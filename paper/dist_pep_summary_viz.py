"""
Simple visualization for category-level peptide sharing.

Input: TSV/CSV with columns:
Category,TotalDistinctPeptides,Peptides_in_1_genomes,Peptides_in_2_genomes,...

Output: a stacked bar plot (PNG) showing how many distinct peptides per category
are observed in N genomes. Zero-only columns are dropped automatically.

Usage:
  python dist_pep_summary_viz.py dist_pep_summary_humangut.tsv out.png
"""

import sys
import csv
from pathlib import Path

import matplotlib.pyplot as plt


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


def drop_zero_columns(share_bins, data):
    # Remove columns where all categories are zero
    keep = []
    for idx in range(len(share_bins)):
        if any(d[idx] != 0 for d in data):
            keep.append(idx)
    share_bins = [share_bins[i] for i in keep]
    data = [[row[i] for i in keep] for row in data]
    return share_bins, data


def plot(categories, share_bins, data, out_path: Path):
    """
    Bar chart with x-axis = "found in N genomes", grouped by category.
    """
    share_bins, data = drop_zero_columns(share_bins, data)
    if not share_bins:
        raise ValueError("All share columns are zero; nothing to plot.")

    # Drop very sparse bins (total peptides across all categories <= threshold)
    MIN_TOTAL_PER_BIN = 2
    MAX_BINS = 15  # cap at "found in 15 genomes" maximum

    totals = [sum(row[idx] for row in data) for idx in range(len(share_bins))]
    keep_idx = [i for i, t in enumerate(totals) if t >= MIN_TOTAL_PER_BIN]

    # Trim to last bin that passes the threshold
    if keep_idx:
        last = max(keep_idx)
        keep_idx = [i for i in range(last + 1) if totals[i] >= MIN_TOTAL_PER_BIN]

    # Apply max bin cap if still too wide
    if len(keep_idx) > MAX_BINS:
        keep_idx = keep_idx[:MAX_BINS]

    share_bins = [share_bins[i] for i in keep_idx]
    data = [[row[i] for i in keep_idx] for row in data]
    if not share_bins:
        raise ValueError("All bins pruned as too sparse; nothing to plot.")

    fig, ax = plt.subplots(figsize=(8, 5))

    x = range(len(share_bins))
    width = 0.15
    num_cats = len(categories)
    colors = plt.cm.tab10.colors

    for idx, (cat, counts) in enumerate(zip(categories, data)):
        offset = (idx - (num_cats - 1) / 2) * width
        ax.bar([p + offset for p in x], counts, width=width,
               label=cat, color=colors[idx % len(colors)])

    ax.set_xticks(list(x))
    ax.set_xticklabels([f"found in {s} genomes" for s in share_bins], rotation=15, ha="right")
    ax.set_ylabel("Distinct peptides")
    ax.set_title("Distinct peptides per marker category by sharing across genomes")
    ax.legend(title="Category")
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        print("Usage: python dist_pep_summary_viz.py <input_tsv> <output_png>")
        sys.exit(1)

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    categories, share_bins, data = load_data(in_path)
    plot(categories, share_bins, data, out_path)
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()

