#!/usr/bin/env python3
"""
download_annotation_maps.py

Download local ID -> description lookup tables for common eggNOG-mapper
annotation columns:

- GO
- EC
- KEGG_ko
- KEGG_Pathway
- KEGG_Module
- KEGG_Reaction
- KEGG_rclass
- BRITE
- KEGG_TC
- CAZy
- BiGG_Reaction
- PFAMs

Example:
    python download_annotation_maps.py --outdir annotation_maps

Notes:
    - KEGG REST is used for KEGG-related tables.
    - GO is parsed from go-basic.obo.
    - Pfam is parsed from Pfam-A.hmm.dat.gz.
    - BiGG reactions are downloaded from the BiGG API.
    - CAZy is generated as a simple family-level table, e.g. GH13 -> Glycoside Hydrolase family 13.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Iterable


KEGG_TABLES = {
    "KEGG_ko": {
        "url": "https://rest.kegg.jp/list/ko",
        "prefix_to_remove": "ko:",
    },
    "KEGG_Pathway": {
        "url": "https://rest.kegg.jp/list/pathway",
        "prefix_to_remove": "path:",
    },
    "KEGG_Module": {
        "url": "https://rest.kegg.jp/list/module",
        "prefix_to_remove": "md:",
    },
    "KEGG_Reaction": {
        "url": "https://rest.kegg.jp/list/reaction",
        "prefix_to_remove": "rn:",
    },
    "KEGG_rclass": {
        "url": "https://rest.kegg.jp/list/rclass",
        "prefix_to_remove": "rc:",
    },
    "EC": {
        "url": "https://rest.kegg.jp/list/enzyme",
        "prefix_to_remove": "ec:",
    },
    "BRITE": {
        "url": "https://rest.kegg.jp/list/brite",
        "prefix_to_remove": "br:",
    },
}

GO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"
BIGG_REACTIONS_URL = "http://bigg.ucsd.edu/api/v2/universal/reactions"


def log(message: str) -> None:
    print(message, flush=True)


def download_text(url: str, retries: int = 3, sleep_seconds: int = 5) -> str:
    last_error = None

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "annotation-map-downloader/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                raw = response.read()
            return raw.decode("utf-8", errors="replace")

        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_error = e
            log(f"  attempt {attempt}/{retries} failed: {e}")
            if attempt < retries:
                time.sleep(sleep_seconds)

    raise RuntimeError(f"Failed to download {url}: {last_error}")


def download_binary(url: str, retries: int = 3, sleep_seconds: int = 5) -> bytes:
    last_error = None

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "annotation-map-downloader/1.0"},
            )
            with urllib.request.urlopen(request, timeout=300) as response:
                return response.read()

        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_error = e
            log(f"  attempt {attempt}/{retries} failed: {e}")
            if attempt < retries:
                time.sleep(sleep_seconds)

    raise RuntimeError(f"Failed to download {url}: {last_error}")


def write_tsv(path: Path, rows: Iterable[tuple[str, str]]) -> int:
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tdescription\n")
        for key, value in rows:
            key = key.strip()
            value = value.strip()
            if not key or not value:
                continue
            f.write(f"{key}\t{value}\n")
            count += 1
    return count


def download_kegg_tables(outdir: Path) -> list[Path]:
    kegg_dir = outdir / "kegg"
    kegg_dir.mkdir(parents=True, exist_ok=True)

    written = []

    for table_name, spec in KEGG_TABLES.items():
        url = spec["url"]
        prefix = spec["prefix_to_remove"]
        outfile = kegg_dir / f"{table_name}.tsv"

        log(f"Downloading {table_name} from KEGG...")
        text = download_text(url)

        rows = []
        for line in text.splitlines():
            if not line.strip():
                continue

            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue

            entry_id, description = parts
            if entry_id.startswith(prefix):
                entry_id = entry_id[len(prefix):]

            rows.append((entry_id, description))

        n = write_tsv(outfile, rows)
        log(f"  wrote {n:,} rows -> {outfile}")
        written.append(outfile)

    return written


def parse_go_obo(obo_text: str) -> list[tuple[str, str]]:
    rows = []

    current_id = None
    current_name = None
    is_obsolete = False
    in_term = False

    def flush_current() -> None:
        nonlocal current_id, current_name, is_obsolete, in_term
        if in_term and current_id and current_name and not is_obsolete:
            rows.append((current_id, current_name))
        current_id = None
        current_name = None
        is_obsolete = False
        in_term = False

    for raw_line in obo_text.splitlines():
        line = raw_line.strip()

        if line == "[Term]":
            flush_current()
            in_term = True

        elif line.startswith("[") and line != "[Term]":
            flush_current()

        elif in_term and line.startswith("id: GO:"):
            current_id = line.replace("id: ", "", 1)

        elif in_term and line.startswith("name: "):
            current_name = line.replace("name: ", "", 1)

        elif in_term and line == "is_obsolete: true":
            is_obsolete = True

        elif line == "":
            flush_current()

    flush_current()
    return rows


def download_go(outdir: Path) -> Path:
    go_dir = outdir / "go"
    go_dir.mkdir(parents=True, exist_ok=True)

    raw_file = go_dir / "go-basic.obo"
    table_file = go_dir / "GO.tsv"

    log("Downloading GO ontology...")
    text = download_text(GO_URL)

    raw_file.write_text(text, encoding="utf-8")
    rows = parse_go_obo(text)
    n = write_tsv(table_file, rows)

    log(f"  wrote raw OBO -> {raw_file}")
    log(f"  wrote {n:,} GO rows -> {table_file}")

    return table_file


def parse_pfam_hmm_dat(text: str) -> list[tuple[str, str]]:
    rows = []

    pfam_id = None
    pfam_acc = None
    description = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")

        if line.startswith("#=GF ID"):
            parts = line.split(maxsplit=2)
            if len(parts) >= 3:
                pfam_id = parts[2].strip()

        elif line.startswith("#=GF AC"):
            parts = line.split(maxsplit=2)
            if len(parts) >= 3:
                # Remove version suffix, e.g. PF00001.27 -> PF00001
                pfam_acc = parts[2].strip().split(".")[0]

        elif line.startswith("#=GF DE"):
            parts = line.split(maxsplit=2)
            if len(parts) >= 3:
                description = parts[2].strip()

        elif line == "//":
            if description:
                if pfam_id:
                    rows.append((pfam_id, description))
                if pfam_acc:
                    rows.append((pfam_acc, description))
            pfam_id = None
            pfam_acc = None
            description = None

    return rows


def download_pfam(outdir: Path) -> Path:
    pfam_dir = outdir / "pfam"
    pfam_dir.mkdir(parents=True, exist_ok=True)

    gz_file = pfam_dir / "Pfam-A.hmm.dat.gz"
    dat_file = pfam_dir / "Pfam-A.hmm.dat"
    table_file = pfam_dir / "PFAM.tsv"

    log("Downloading Pfam metadata...")
    binary = download_binary(PFAM_URL)
    gz_file.write_bytes(binary)

    log("  decompressing Pfam metadata...")
    text = gzip.decompress(binary).decode("utf-8", errors="replace")
    dat_file.write_text(text, encoding="utf-8")

    rows = parse_pfam_hmm_dat(text)
    n = write_tsv(table_file, rows)

    log(f"  wrote compressed raw file -> {gz_file}")
    log(f"  wrote decompressed raw file -> {dat_file}")
    log(f"  wrote {n:,} PFAM rows -> {table_file}")

    return table_file


def download_bigg_reactions(outdir: Path) -> Path:
    bigg_dir = outdir / "bigg"
    bigg_dir.mkdir(parents=True, exist_ok=True)

    json_file = bigg_dir / "bigg_universal_reactions.json"
    table_file = bigg_dir / "BiGG_Reaction.tsv"

    log("Downloading BiGG universal reactions...")
    text = download_text(BIGG_REACTIONS_URL)
    json_file.write_text(text, encoding="utf-8")

    data = json.loads(text)
    results = data.get("results", [])

    rows = []
    for item in results:
        bigg_id = item.get("bigg_id")
        name = item.get("name")
        if bigg_id and name:
            rows.append((bigg_id, name))

    n = write_tsv(table_file, rows)

    log(f"  wrote raw JSON -> {json_file}")
    log(f"  wrote {n:,} BiGG reaction rows -> {table_file}")

    return table_file


def generate_cazy_table(outdir: Path, max_family_number: int = 500) -> Path:
    """
    Generate a simple CAZy family-level lookup.

    This does not download the curated CAZy family pages. It creates useful
    generic family descriptions such as:
        GH13  Glycoside Hydrolase family 13
        CBM50 Carbohydrate-Binding Module family 50

    This is often sufficient for annotating eggNOG-mapper CAZy columns.
    """

    cazy_dir = outdir / "cazy"
    cazy_dir.mkdir(parents=True, exist_ok=True)

    table_file = cazy_dir / "CAZy.tsv"

    prefix_names = {
        "GH": "Glycoside Hydrolase",
        "GT": "GlycosylTransferase",
        "PL": "Polysaccharide Lyase",
        "CE": "Carbohydrate Esterase",
        "AA": "Auxiliary Activity",
        "CBM": "Carbohydrate-Binding Module",
    }

    rows = []
    for prefix, name in prefix_names.items():
        for i in range(1, max_family_number + 1):
            rows.append((f"{prefix}{i}", f"{name} family {i}"))

    n = write_tsv(table_file, rows)
    log(f"Generated {n:,} generic CAZy family rows -> {table_file}")

    return table_file


def download_kegg_tc(outdir: Path) -> Path:
    """
    KEGG_TC values in eggNOG output are usually TCDB-like IDs.

    KEGG REST does not provide a simple 'list/tc' database equivalent in the
    same way as KO/pathway/module/reaction. This function creates an empty
    placeholder table and documents the situation.

    If you later want true TCDB descriptions, they should be downloaded from
    TCDB directly and parsed separately.
    """
    tc_dir = outdir / "kegg_tc"
    tc_dir.mkdir(parents=True, exist_ok=True)

    table_file = tc_dir / "KEGG_TC.tsv"
    note_file = tc_dir / "README.txt"

    table_file.write_text("id\tdescription\n", encoding="utf-8")

    note_file.write_text(
        "KEGG_TC values from eggNOG-mapper are transport-classification IDs. "
        "KEGG REST does not provide a simple list/tc endpoint equivalent to "
        "list/ko or list/pathway. This script therefore creates an empty "
        "placeholder table. For detailed transport-classification descriptions, "
        "download and parse TCDB data separately.\n",
        encoding="utf-8",
    )

    log(f"Created empty KEGG_TC placeholder -> {table_file}")
    log(f"  note -> {note_file}")

    return table_file


def create_merged_table(outdir: Path) -> Path:
    """
    Merge all per-database tables into one long-format table:

        source    id    description

    This makes downstream lookup easy.
    """

    merged_file = outdir / "annotation_lookup_merged.tsv"

    table_specs = [
        ("KEGG_ko", outdir / "kegg" / "KEGG_ko.tsv"),
        ("KEGG_Pathway", outdir / "kegg" / "KEGG_Pathway.tsv"),
        ("KEGG_Module", outdir / "kegg" / "KEGG_Module.tsv"),
        ("KEGG_Reaction", outdir / "kegg" / "KEGG_Reaction.tsv"),
        ("KEGG_rclass", outdir / "kegg" / "KEGG_rclass.tsv"),
        ("EC", outdir / "kegg" / "EC.tsv"),
        ("BRITE", outdir / "kegg" / "BRITE.tsv"),
        ("GO", outdir / "go" / "GO.tsv"),
        ("PFAM", outdir / "pfam" / "PFAM.tsv"),
        ("BiGG_Reaction", outdir / "bigg" / "BiGG_Reaction.tsv"),
        ("CAZy", outdir / "cazy" / "CAZy.tsv"),
        ("KEGG_TC", outdir / "kegg_tc" / "KEGG_TC.tsv"),
    ]

    n = 0
    with merged_file.open("w", encoding="utf-8", newline="\n") as out:
        out.write("source\tid\tdescription\n")

        for source, table in table_specs:
            if not table.exists():
                continue

            with table.open("r", encoding="utf-8", errors="replace") as f:
                header = next(f, None)
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t", 1)
                    if len(parts) != 2:
                        continue
                    entry_id, description = parts
                    out.write(f"{source}\t{entry_id}\t{description}\n")
                    n += 1

    log(f"Created merged lookup table with {n:,} rows -> {merged_file}")
    return merged_file


def write_readme(outdir: Path) -> Path:
    readme = outdir / "README_annotation_maps.txt"
    readme.write_text(
        """Annotation lookup tables generated by download_annotation_maps.py

Main output:
    annotation_lookup_merged.tsv

Format:
    source    id    description

Per-database tables:
    kegg/KEGG_ko.tsv
    kegg/KEGG_Pathway.tsv
    kegg/KEGG_Module.tsv
    kegg/KEGG_Reaction.tsv
    kegg/KEGG_rclass.tsv
    kegg/EC.tsv
    kegg/BRITE.tsv
    go/GO.tsv
    pfam/PFAM.tsv
    bigg/BiGG_Reaction.tsv
    cazy/CAZy.tsv
    kegg_tc/KEGG_TC.tsv

Important notes:
    1. KEGG identifiers from eggNOG may contain prefixes such as ko:K01810,
       map00010, ko00010, rn:R00001, md:M00001, etc. You may need to normalise
       prefixes before joining.

    2. GO values may contain comma-separated terms. Split by comma first.

    3. PFAM values from eggNOG may be short names such as DEAD or Helicase_C,
       not only PFxxxxx accessions. The PFAM table includes both short names
       and accessions when available.

    4. CAZy descriptions here are generic family-level names, not curated
       substrate/activity descriptions.

    5. KEGG_TC is created as an empty placeholder because KEGG REST does not
       provide a simple list/tc endpoint equivalent to list/ko.
""",
        encoding="utf-8",
    )
    log(f"Wrote README -> {readme}")
    return readme


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download ID-to-description lookup tables for eggNOG-mapper annotations."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output folder for downloaded annotation maps.",
    )
    parser.add_argument(
        "--skip-pfam",
        action="store_true",
        help="Skip Pfam download. Useful if the EBI FTP/HTTPS server is slow.",
    )
    parser.add_argument(
        "--skip-bigg",
        action="store_true",
        help="Skip BiGG reaction download.",
    )
    parser.add_argument(
        "--skip-go",
        action="store_true",
        help="Skip GO ontology download.",
    )
    parser.add_argument(
        "--skip-kegg",
        action="store_true",
        help="Skip KEGG downloads.",
    )
    parser.add_argument(
        "--cazy-max-family-number",
        type=int,
        default=500,
        help="Maximum family number to generate for each CAZy prefix. Default: 500.",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    log(f"Output folder: {outdir}")

    failed = []

    if not args.skip_kegg:
        try:
            download_kegg_tables(outdir)
        except Exception as e:
            failed.append(("KEGG", e))
            log(f"ERROR: KEGG download failed: {e}")

    if not args.skip_go:
        try:
            download_go(outdir)
        except Exception as e:
            failed.append(("GO", e))
            log(f"ERROR: GO download failed: {e}")

    if not args.skip_pfam:
        try:
            download_pfam(outdir)
        except Exception as e:
            failed.append(("PFAM", e))
            log(f"ERROR: PFAM download failed: {e}")

    if not args.skip_bigg:
        try:
            download_bigg_reactions(outdir)
        except Exception as e:
            failed.append(("BiGG", e))
            log(f"ERROR: BiGG download failed: {e}")

    try:
        generate_cazy_table(outdir, max_family_number=args.cazy_max_family_number)
    except Exception as e:
        failed.append(("CAZy", e))
        log(f"ERROR: CAZy table generation failed: {e}")

    try:
        download_kegg_tc(outdir)
    except Exception as e:
        failed.append(("KEGG_TC", e))
        log(f"ERROR: KEGG_TC placeholder generation failed: {e}")

    try:
        create_merged_table(outdir)
    except Exception as e:
        failed.append(("merged table", e))
        log(f"ERROR: merged table generation failed: {e}")

    write_readme(outdir)

    if failed:
        log("\nCompleted with some failures:")
        for name, error in failed:
            log(f"  - {name}: {error}")
        return 1

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
