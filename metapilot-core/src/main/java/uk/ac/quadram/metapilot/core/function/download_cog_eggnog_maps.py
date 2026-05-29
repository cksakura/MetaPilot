#!/usr/bin/env python3
"""
download_cog_eggnog_maps.py

Download / generate lookup tables for:

1. NCBI COG 2024:
   - COG ID -> name / category / gene / pathway
   - COG functional category letter -> description

2. eggNOG_OGs helper table:
   - This script does NOT download the full eggNOG hierarchy by default because
     eggNOG contains millions of OGs and some files are very large.
   - Instead, it creates:
       a) a parser-ready note explaining how to handle eggNOG_OGs values
       b) a COG-compatible lookup using NCBI COG IDs, useful for entries like:
          COG1197@1|root, COG1197@2|Bacteria

Example:
    python download_cog_eggnog_maps.py --outdir annotation_maps

Outputs:
    annotation_maps/cog/COG.tsv
    annotation_maps/cog/COG_category.tsv
    annotation_maps/cog/COG_long.tsv
    annotation_maps/eggnog_ogs/README_eggNOG_OGs.txt
"""

from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


COG_DEF_URL = "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab"
COG_FUN_URL = "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab"


def log(msg: str) -> None:
    print(msg, flush=True)


def download_text(url: str, retries: int = 3, sleep_seconds: int = 5) -> str:
    last_error = None

    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "cog-eggnog-map-downloader/1.0"},
            )
            with urllib.request.urlopen(req, timeout=180) as response:
                raw = response.read()
            return raw.decode("utf-8", errors="replace")

        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_error = e
            log(f"  attempt {attempt}/{retries} failed: {e}")
            if attempt < retries:
                time.sleep(sleep_seconds)

    raise RuntimeError(f"Failed to download {url}: {last_error}")


def parse_cog_def(text: str) -> list[dict[str, str]]:
    """
    Parse NCBI COG2024 cog-24.def.tab.

    Expected columns:
        1. COG ID
        2. COG functional category
        3. COG name
        4. Gene name, optional
        5. Functional pathway, optional
        6. PubMed IDs, optional
        7. PDB IDs, optional
    """
    rows = []

    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        parts = parts + [""] * (7 - len(parts))

        cog_id, category, name, gene, pathway, pubmed, pdb = parts[:7]

        if not cog_id.startswith("COG"):
            continue

        rows.append(
            {
                "id": cog_id,
                "category": category,
                "name": name,
                "gene": gene,
                "pathway": pathway,
                "pubmed": pubmed,
                "pdb": pdb,
            }
        )

    return rows


def parse_cog_fun(text: str) -> list[dict[str, str]]:
    """
    Parse NCBI COG2024 cog-24.fun.tab.

    Expected columns:
        1. Functional category ID, one letter
        2. Functional group number
        3. RGB color
        4. Functional category description
    """
    rows = []

    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        parts = parts + [""] * (4 - len(parts))

        category_id, group, color, description = parts[:4]

        if not category_id:
            continue

        rows.append(
            {
                "id": category_id,
                "group": group,
                "color": color,
                "description": description,
            }
        )

    return rows


def write_cog_outputs(outdir: Path, cog_rows: list[dict[str, str]], cat_rows: list[dict[str, str]]) -> None:
    cog_dir = outdir / "cog"
    cog_dir.mkdir(parents=True, exist_ok=True)

    # Simple ID -> description table
    cog_simple = cog_dir / "COG.tsv"
    with cog_simple.open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tdescription\n")
        for r in cog_rows:
            desc = r["name"]
            if r["gene"]:
                desc = f'{r["gene"]}; {desc}'
            f.write(f'{r["id"]}\t{desc}\n')

    # COG category letters
    cog_category = cog_dir / "COG_category.tsv"
    with cog_category.open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tdescription\tgroup\tcolor\n")
        for r in cat_rows:
            f.write(f'{r["id"]}\t{r["description"]}\t{r["group"]}\t{r["color"]}\n')

    # More detailed COG table
    cog_long = cog_dir / "COG_long.tsv"
    with cog_long.open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tcategory\tname\tgene\tpathway\tpubmed\tpdb\n")
        for r in cog_rows:
            f.write(
                f'{r["id"]}\t{r["category"]}\t{r["name"]}\t{r["gene"]}\t'
                f'{r["pathway"]}\t{r["pubmed"]}\t{r["pdb"]}\n'
            )

    log(f"Wrote {len(cog_rows):,} COG rows -> {cog_simple}")
    log(f"Wrote {len(cat_rows):,} COG category rows -> {cog_category}")
    log(f"Wrote detailed COG table -> {cog_long}")


def write_eggnog_ogs_readme(outdir: Path) -> None:
    eggnog_dir = outdir / "eggnog_ogs"
    eggnog_dir.mkdir(parents=True, exist_ok=True)

    readme = eggnog_dir / "README_eggNOG_OGs.txt"
    readme.write_text(
        """How to interpret eggNOG_OGs from eggNOG-mapper

Example value:
    COG1197@1|root,COG1197@2|Bacteria,1TPF1@1239|Firmicutes,248D8@186801|Clostridia

Each item usually has this form:
    OG_ID@taxid|taxonomic_level_name

Examples:
    COG1197@2|Bacteria
        OG_ID = COG1197
        taxid = 2
        taxonomic level = Bacteria

    1TPF1@1239|Firmicutes
        OG_ID = 1TPF1
        taxid = 1239
        taxonomic level = Firmicutes

Practical recommendation:
    1. For OG_ID values starting with COG, use cog/COG.tsv or cog/COG_long.tsv
       to map COG1197 -> description.

    2. For non-COG eggNOG OGs, such as 1TPF1 or 248D8, the ID represents an
       eggNOG orthologous group at a specific taxonomic level. These are useful
       for orthology/provenance, but they are usually not the best function label
       for reports.

    3. For biological interpretation, use:
           Description
           Preferred_name
           GOs
           EC
           KEGG_ko
           KEGG_Pathway
           PFAMs
           CAZy
       rather than trying to force every eggNOG OG ID into a readable function.

Suggested parser:
    Split eggNOG_OGs by comma.
    For each item, parse with:
        ^([^@]+)@(\\d+)\\|(.+)$

This gives:
    og_id, taxid, taxonomic_level
""",
        encoding="utf-8",
    )

    log(f"Wrote eggNOG_OGs note -> {readme}")


def extract_eggnog_ogs(value: str) -> list[dict[str, str]]:
    """
    Utility function you can copy into downstream scripts.

    Input:
        COG1197@1|root,COG1197@2|Bacteria,1TPF1@1239|Firmicutes

    Output:
        [
          {"og_id": "COG1197", "taxid": "1", "tax_level": "root"},
          {"og_id": "COG1197", "taxid": "2", "tax_level": "Bacteria"},
          {"og_id": "1TPF1", "taxid": "1239", "tax_level": "Firmicutes"},
        ]
    """
    if value is None or value == "-" or not str(value).strip():
        return []

    out = []
    for item in str(value).split(","):
        item = item.strip()
        if not item:
            continue

        m = re.match(r"^([^@]+)@(\d+)\|(.+)$", item)
        if m:
            og_id, taxid, tax_level = m.groups()
            out.append(
                {
                    "og_id": og_id,
                    "taxid": taxid,
                    "tax_level": tax_level,
                }
            )
        else:
            out.append(
                {
                    "og_id": item,
                    "taxid": "",
                    "tax_level": "",
                }
            )

    return out


def append_to_merged_lookup(outdir: Path) -> None:
    """
    If annotation_lookup_merged.tsv already exists from the previous downloader,
    append COG and COG_category rows. Existing rows are preserved.
    """
    merged = outdir / "annotation_lookup_merged.tsv"
    cog_simple = outdir / "cog" / "COG.tsv"
    cog_category = outdir / "cog" / "COG_category.tsv"

    if not merged.exists():
        with merged.open("w", encoding="utf-8", newline="\n") as f:
            f.write("source\tid\tdescription\n")

    existing = set()
    with merged.open("r", encoding="utf-8", errors="replace") as f:
        next(f, None)
        for line in f:
            parts = line.rstrip("\n").split("\t", 2)
            if len(parts) >= 2:
                existing.add((parts[0], parts[1]))

    added = 0
    with merged.open("a", encoding="utf-8", newline="\n") as out:
        if cog_simple.exists():
            with cog_simple.open("r", encoding="utf-8", errors="replace") as f:
                next(f, None)
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    cog_id, desc = line.split("\t", 1)
                    key = ("COG", cog_id)
                    if key not in existing:
                        out.write(f"COG\t{cog_id}\t{desc}\n")
                        added += 1

        if cog_category.exists():
            with cog_category.open("r", encoding="utf-8", errors="replace") as f:
                next(f, None)
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    cat_id = parts[0]
                    desc = parts[1] if len(parts) > 1 else ""
                    key = ("COG_category", cat_id)
                    if key not in existing:
                        out.write(f"COG_category\t{cat_id}\t{desc}\n")
                        added += 1

    log(f"Appended {added:,} COG-related rows -> {merged}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Download COG and eggNOG_OG helper lookup tables.")
    parser.add_argument("--outdir", required=True, help="Output folder, e.g. annotation_maps")
    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    log(f"Output folder: {outdir}")

    log("Downloading NCBI COG definitions...")
    cog_def_text = download_text(COG_DEF_URL)

    log("Downloading NCBI COG functional categories...")
    cog_fun_text = download_text(COG_FUN_URL)

    cog_rows = parse_cog_def(cog_def_text)
    cat_rows = parse_cog_fun(cog_fun_text)

    # Save raw files too
    raw_dir = outdir / "cog" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    (raw_dir / "cog-24.def.tab").write_text(cog_def_text, encoding="utf-8")
    (raw_dir / "cog-24.fun.tab").write_text(cog_fun_text, encoding="utf-8")

    write_cog_outputs(outdir, cog_rows, cat_rows)
    write_eggnog_ogs_readme(outdir)
    append_to_merged_lookup(outdir)

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
