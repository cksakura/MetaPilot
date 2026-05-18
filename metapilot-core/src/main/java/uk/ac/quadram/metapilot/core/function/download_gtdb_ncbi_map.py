#!/usr/bin/env python3
"""
download_gtdb_ncbi_map.py

Download GTDB metadata files and build genome-level mapping tables between
GTDB taxonomy and NCBI taxonomy.

Main outputs:
    gtdb_ncbi_mapping_all.tsv
    gtdb_ncbi_mapping_bacteria.tsv
    gtdb_ncbi_mapping_archaea.tsv

Useful columns:
    domain
    accession
    accession_no_prefix
    ncbi_genbank_assembly_accession
    ncbi_refseq_assembly_accession
    ncbi_taxid
    ncbi_organism_name
    gtdb_taxonomy
    ncbi_taxonomy
    gtdb_domain ... gtdb_species
    ncbi_domain ... ncbi_species
    gtdb_release

Example:
    python download_gtdb_ncbi_map.py --outdir gtdb_ncbi_map

Use a fixed GTDB release instead of latest:
    python download_gtdb_ncbi_map.py \
        --outdir gtdb_ncbi_map \
        --release-url https://data.gtdb.aau.ecogenomic.org/releases/release220/220.0/

Notes:
    - By default this uses https://data.gtdb.aau.ecogenomic.org/releases/latest/
    - GTDB and NCBI taxonomy are not always one-to-one. This script creates a
      genome/assembly-level cross-reference, which is usually the safest mapping.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional


DEFAULT_RELEASE_URL = "https://data.gtdb.aau.ecogenomic.org/releases/latest/"

METADATA_FILES = {
    "bacteria": "bac120_metadata.tsv.gz",
    "archaea": "ar53_metadata.tsv.gz",
}

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def log(message: str) -> None:
    print(message, flush=True)


def ensure_trailing_slash(url: str) -> str:
    return url if url.endswith("/") else url + "/"


def download_binary(url: str, retries: int = 3, sleep_seconds: int = 5) -> bytes:
    last_error = None

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "gtdb-ncbi-map-downloader/1.0"},
            )
            with urllib.request.urlopen(request, timeout=600) as response:
                return response.read()

        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_error = e
            log(f"  attempt {attempt}/{retries} failed: {e}")
            if attempt < retries:
                time.sleep(sleep_seconds)

    raise RuntimeError(f"Failed to download {url}: {last_error}")


def download_text(url: str, retries: int = 3, sleep_seconds: int = 5) -> str:
    raw = download_binary(url, retries=retries, sleep_seconds=sleep_seconds)
    return raw.decode("utf-8", errors="replace")


def maybe_download(url: str, outfile: Path, force: bool = False) -> Path:
    if outfile.exists() and outfile.stat().st_size > 0 and not force:
        log(f"Using existing file -> {outfile}")
        return outfile

    log(f"Downloading {url}")
    raw = download_binary(url)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfile.write_bytes(raw)
    log(f"  wrote {outfile} ({outfile.stat().st_size:,} bytes)")
    return outfile


def read_gzip_tsv(path: Path) -> Iterable[Dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield {k: (v if v is not None else "") for k, v in row.items()}


def normalise_accession(accession: str) -> str:
    """
    Convert GTDB-style accessions to assembly accessions where possible.

    Examples:
        RS_GCF_000005845.2 -> GCF_000005845.2
        GB_GCA_000006155.2 -> GCA_000006155.2
    """
    if accession is None:
        return ""

    accession = accession.strip()
    accession = re.sub(r"^(RS_|GB_)", "", accession)
    return accession


def split_taxonomy(taxonomy: str, prefix: str) -> Dict[str, str]:
    """
    Split GTDB/NCBI taxonomy string into rank columns.

    Input example:
        d__Bacteria;p__Bacillota;c__Clostridia;o__...

    Output:
        gtdb_domain, gtdb_phylum, ...
    """
    out = {f"{prefix}_{rank}": "" for rank in RANKS}

    if not taxonomy or taxonomy == "none":
        return out

    rank_prefix_to_name = {
        "d__": "domain",
        "p__": "phylum",
        "c__": "class",
        "o__": "order",
        "f__": "family",
        "g__": "genus",
        "s__": "species",
    }

    for item in taxonomy.split(";"):
        item = item.strip()
        if not item:
            continue

        matched = False
        for rank_prefix, rank_name in rank_prefix_to_name.items():
            if item.startswith(rank_prefix):
                out[f"{prefix}_{rank_name}"] = item[len(rank_prefix):]
                matched = True
                break

        # Some metadata fields may contain unprefixed taxonomy. Keep a safe fallback.
        if not matched:
            continue

    return out


def pick_first(row: Dict[str, str], candidates: List[str]) -> str:
    for c in candidates:
        if c in row and row[c] not in (None, ""):
            return str(row[c])
    return ""


def build_mapping_rows(
    metadata_file: Path,
    domain: str,
    release_name: str,
) -> List[Dict[str, str]]:
    rows = []

    for row in read_gzip_tsv(metadata_file):
        accession = pick_first(row, ["accession", "genome_id", "gid"])
        accession_no_prefix = normalise_accession(accession)

        gtdb_taxonomy = pick_first(row, ["gtdb_taxonomy", "gtdb_tax"])
        ncbi_taxonomy = pick_first(row, ["ncbi_taxonomy", "ncbi_tax"])

        out = {
            "domain": domain,
            "accession": accession,
            "accession_no_prefix": accession_no_prefix,
            "ncbi_genbank_assembly_accession": pick_first(
                row,
                [
                    "ncbi_genbank_assembly_accession",
                    "genbank_assembly_accession",
                    "ncbi_genbank_accession",
                ],
            ),
            "ncbi_refseq_assembly_accession": pick_first(
                row,
                [
                    "ncbi_refseq_assembly_accession",
                    "refseq_assembly_accession",
                    "ncbi_refseq_accession",
                ],
            ),
            "ncbi_taxid": pick_first(row, ["ncbi_taxid", "taxid", "ncbi_species_taxid"]),
            "ncbi_organism_name": pick_first(
                row,
                [
                    "ncbi_organism_name",
                    "organism_name",
                    "ncbi_organism",
                ],
            ),
            "ncbi_species_taxid": pick_first(row, ["ncbi_species_taxid"]),
            "ncbi_strain_identifiers": pick_first(row, ["ncbi_strain_identifiers"]),
            "gtdb_representative": pick_first(row, ["gtdb_representative"]),
            "gtdb_type_designation": pick_first(row, ["gtdb_type_designation"]),
            "gtdb_type_designation_sources": pick_first(row, ["gtdb_type_designation_sources"]),
            "ncbi_type_material_designation": pick_first(row, ["ncbi_type_material_designation"]),
            "gtdb_taxonomy": gtdb_taxonomy,
            "ncbi_taxonomy": ncbi_taxonomy,
            "gtdb_release": release_name,
        }

        out.update(split_taxonomy(gtdb_taxonomy, "gtdb"))
        out.update(split_taxonomy(ncbi_taxonomy, "ncbi"))

        rows.append(out)

    return rows


def write_tsv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    log(f"Wrote {len(rows):,} rows -> {path}")


def get_release_name(release_url: str, outdir: Path, force: bool = False) -> str:
    """
    Try to download the GTDB VERSION file. If unavailable, use release URL as label.
    """
    version_file = outdir / "VERSION"

    try:
        if version_file.exists() and version_file.stat().st_size > 0 and not force:
            text = version_file.read_text(encoding="utf-8", errors="replace")
        else:
            text = download_text(ensure_trailing_slash(release_url) + "VERSION")
            version_file.write_text(text, encoding="utf-8")

        version = text.strip().replace("\n", " ")
        return version if version else release_url

    except Exception as e:
        log(f"Could not read VERSION file: {e}")
        return release_url


def download_correspondence_excels(release_url: str, outdir: Path, force: bool = False) -> None:
    """
    Download optional GTDB-vs-NCBI / NCBI-vs-GTDB Excel files.

    The exact filenames historically contain release placeholders, so this
    function tries several common names. If they do not exist, this step is
    skipped without failing the main genome-level mapping.
    """
    aux_dir = outdir / "auxillary_files"
    aux_dir.mkdir(parents=True, exist_ok=True)

    # Common names used in recent GTDB release directories.
    candidates = [
        "gtdb_vs_ncbi_bacteria.xlsx",
        "gtdb_vs_ncbi_archaea.xlsx",
        "ncbi_vs_gtdb_bacteria.xlsx",
        "ncbi_vs_gtdb_archaea.xlsx",
    ]

    base = ensure_trailing_slash(release_url) + "auxillary_files/"

    for name in candidates:
        url = base + name
        outfile = aux_dir / name
        try:
            maybe_download(url, outfile, force=force)
        except Exception as e:
            log(f"Optional Excel file not downloaded ({name}): {e}")


def write_readme(outdir: Path, release_url: str, release_name: str) -> None:
    readme = outdir / "README_gtdb_ncbi_mapping.txt"
    readme.write_text(
        f"""GTDB-NCBI mapping files

Release URL:
    {release_url}

Release label:
    {release_name}

Main tables:
    gtdb_ncbi_mapping_all.tsv
    gtdb_ncbi_mapping_bacteria.tsv
    gtdb_ncbi_mapping_archaea.tsv

Interpretation:
    These are genome/assembly-level mapping tables. Each row is a GTDB genome
    record with associated GTDB taxonomy and NCBI taxonomy/TaxID fields from
    the GTDB metadata file.

Important:
    GTDB and NCBI taxonomy are not perfectly one-to-one. The safest stable
    bridge is genome/assembly accession, not taxon name alone.

Useful join keys:
    accession
        GTDB accession, often prefixed with RS_ or GB_.

    accession_no_prefix
        Same accession with RS_ or GB_ removed, e.g. GCF_... or GCA_....

    ncbi_genbank_assembly_accession
    ncbi_refseq_assembly_accession
        NCBI assembly accession fields when present in the GTDB metadata.

    ncbi_taxid
        NCBI taxonomic identifier from GTDB metadata.

Taxonomy columns:
    gtdb_taxonomy, ncbi_taxonomy
        Full semicolon-separated taxonomy strings.

    gtdb_domain ... gtdb_species
    ncbi_domain ... ncbi_species
        Parsed rank columns.
""",
        encoding="utf-8",
    )
    log(f"Wrote README -> {readme}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download GTDB metadata and create NCBI-GTDB mapping tables."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output folder, e.g. gtdb_ncbi_map",
    )
    parser.add_argument(
        "--release-url",
        default=DEFAULT_RELEASE_URL,
        help=f"GTDB release URL. Default: {DEFAULT_RELEASE_URL}",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even if they already exist.",
    )
    parser.add_argument(
        "--skip-excel-correspondence",
        action="store_true",
        help="Do not try to download optional GTDB-vs-NCBI Excel correspondence files.",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    release_url = ensure_trailing_slash(args.release_url)

    log(f"Output folder: {outdir}")
    log(f"GTDB release URL: {release_url}")

    release_name = get_release_name(release_url, outdir, force=args.force)
    log(f"GTDB release label: {release_name}")

    raw_dir = outdir / "raw_metadata"
    raw_dir.mkdir(parents=True, exist_ok=True)

    all_rows: List[Dict[str, str]] = []
    rows_by_domain: Dict[str, List[Dict[str, str]]] = {}

    for domain, filename in METADATA_FILES.items():
        url = release_url + filename
        outfile = raw_dir / filename

        try:
            maybe_download(url, outfile, force=args.force)
            rows = build_mapping_rows(outfile, domain, release_name)
            rows_by_domain[domain] = rows
            all_rows.extend(rows)
        except Exception as e:
            log(f"ERROR: failed to process {domain} metadata ({filename}): {e}")
            rows_by_domain[domain] = []

    columns = [
        "domain",
        "accession",
        "accession_no_prefix",
        "ncbi_genbank_assembly_accession",
        "ncbi_refseq_assembly_accession",
        "ncbi_taxid",
        "ncbi_species_taxid",
        "ncbi_organism_name",
        "ncbi_strain_identifiers",
        "gtdb_representative",
        "gtdb_type_designation",
        "gtdb_type_designation_sources",
        "ncbi_type_material_designation",
        "gtdb_taxonomy",
        "ncbi_taxonomy",
        "gtdb_domain",
        "gtdb_phylum",
        "gtdb_class",
        "gtdb_order",
        "gtdb_family",
        "gtdb_genus",
        "gtdb_species",
        "ncbi_domain",
        "ncbi_phylum",
        "ncbi_class",
        "ncbi_order",
        "ncbi_family",
        "ncbi_genus",
        "ncbi_species",
        "gtdb_release",
    ]

    write_tsv(outdir / "gtdb_ncbi_mapping_all.tsv", all_rows, columns)

    for domain in ["bacteria", "archaea"]:
        write_tsv(
            outdir / f"gtdb_ncbi_mapping_{domain}.tsv",
            rows_by_domain.get(domain, []),
            columns,
        )

    if not args.skip_excel_correspondence:
        download_correspondence_excels(release_url, outdir, force=args.force)

    write_readme(outdir, release_url, release_name)

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
