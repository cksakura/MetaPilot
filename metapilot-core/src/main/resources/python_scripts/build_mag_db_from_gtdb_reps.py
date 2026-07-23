#!/usr/bin/env python3
"""
Build a MetaPilot-compatible MAG database from a folder of downloaded GTDB
representative genomes (proteins + eggNOG), e.g. GTDB_rep_db.

Expected input layout:
  proteins/<accession>.faa                 Prodigal protein FASTA per genome
  eggnog/<accession>.emapper.annotations   eggNOG-mapper annotations per genome

Where <accession> is an NCBI assembly accession such as GCA_000154345.1.

This differs from build_mag_db_from_binning.py in two ways:
  1. Genome IDs are accessions that contain '_' and '.'. MetaPilot parses protein
     accessions as genomeID_proteinID using the FIRST underscore, so the genome ID
     itself must contain no underscore. Each accession is remapped to a clean,
     underscore-free ID (GCA_000154345.1 -> GCA000154345v1) and recorded in
     genome_id_map.tsv.
  2. There is no gtdbtk_retry/ folder. Lineage comes from a taxonomy table keyed by
     accession (the non-redundant GTDB species list produced earlier), not GTDB-Tk
     summaries.

Output layout (MetaPilot local MAG DB version folder):
  original_db/<cleanID>.faa
  eggNOG/<cleanID>_eggNOG.tsv
  genomes-all_metadata.tsv
  genome_id_map.tsv
  protein_id_map.tsv

FASTA headers follow the MGnify / MetaPilot convention:
  >cleanID_00001 30S ribosomal protein S2
  >cleanID_00002 hypothetical protein

Only genomes with BOTH a protein FASTA and an eggNOG annotation are written
(MetaPilot pairing validation). Marker2/Marker4 are generated later by MetaPilot
during database registration.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


# Reuse the binning converter's helpers (eggNOG parsing, name resolution, FASTA
# rewrite, eggNOG rewrite, taxonomy consensus) so behaviour stays consistent.
_SIBLING = Path(__file__).resolve().parent / "build_mag_db_from_binning.py"
_spec = importlib.util.spec_from_file_location("build_mag_db_from_binning", _SIBLING)
if _spec is None or _spec.loader is None:
    raise ImportError(f"Cannot load sibling module: {_SIBLING}")
_binning = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_binning)

parse_eggnog_annotation = _binning.parse_eggnog_annotation
summarize_eggnog_taxonomy = _binning.summarize_eggnog_taxonomy
convert_fasta = _binning.convert_fasta
write_eggnog = _binning.write_eggnog
natural_key = _binning.natural_key
ensure_output_ready = _binning.ensure_output_ready

EGGNOG_ANNOT_SUFFIX = _binning.EGGNOG_ANNOT_SUFFIX


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a GTDB representative genome folder (proteins + eggNOG) "
        "into a MetaPilot MAG database."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input root (contains proteins/ and eggnog/).",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output MAG DB version folder (original_db/, eggNOG/, genomes-all_metadata.tsv).",
    )
    parser.add_argument(
        "--taxonomy-table",
        type=Path,
        required=True,
        help="TSV mapping accession -> GTDB species/lineage. Accepts the "
        "gtdb_retry_nonredundant_species_list.tsv format (accession + lineage columns).",
    )
    parser.add_argument("--proteins-dir", type=Path, default=None,
                        help="Override proteins directory (default: <input>/proteins).")
    parser.add_argument("--eggnog-dir", type=Path, default=None,
                        help="Override eggNOG directory (default: <input>/eggnog).")
    parser.add_argument("--db-name", default="original_db",
                        help="Protein FASTA folder name (default: original_db).")
    parser.add_argument("--func-name", default="eggNOG",
                        help="Functional annotation folder name (default: eggNOG).")
    parser.add_argument("--taxa-name", default="genomes-all_metadata.tsv",
                        help="Taxonomy metadata filename (default: genomes-all_metadata.tsv).")
    parser.add_argument("--force", action="store_true", help="Overwrite existing output files.")
    return parser.parse_args(argv)


def find_protein_files(proteins_dir: Path) -> Dict[str, Path]:
    files: Dict[str, Path] = {}
    for path in sorted(proteins_dir.glob("*.faa")):
        files[path.name[: -len(".faa")]] = path
    for path in sorted(proteins_dir.glob("*.fasta")):
        files.setdefault(path.name[: -len(".fasta")], path)
    return files


def find_eggnog_files(eggnog_dir: Path) -> Dict[str, Path]:
    files: Dict[str, Path] = {}
    for path in sorted(eggnog_dir.glob(f"*{EGGNOG_ANNOT_SUFFIX}")):
        files[path.name[: -len(EGGNOG_ANNOT_SUFFIX)]] = path
    return files


def clean_genome_id(accession: str) -> str:
    """
    Turn an assembly accession into an underscore-free, dot-free MetaPilot genome ID.
      GCA_000154345.1 -> GCA000154345v1
    Any remaining non-alphanumeric characters are stripped.
    """
    cid = accession.strip().replace("_", "").replace(".", "v")
    cid = re.sub(r"[^A-Za-z0-9]", "", cid)
    if not cid:
        raise ValueError(f"Could not derive a clean genome ID from accession: {accession!r}")
    return cid


def _first_present(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    lower = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def load_taxonomy_table(path: Path) -> Dict[str, Dict[str, str]]:
    """
    Load accession -> {species, genus, family, lineage} from a flexible TSV.

    Recognised accession columns: GTDB_representative_accession, accession, Accession,
    closest_genome_reference.
    Recognised lineage columns: Example_lineage, Lineage, classification, gtdb_taxonomy.
    Recognised species columns: Species, gtdb_species.
    """
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Empty taxonomy table: {path}")
        acc_col = _first_present(
            reader.fieldnames,
            ["GTDB_representative_accession", "accession", "Accession", "closest_genome_reference"],
        )
        lin_col = _first_present(
            reader.fieldnames,
            ["Example_lineage", "Lineage", "classification", "gtdb_taxonomy"],
        )
        sp_col = _first_present(reader.fieldnames, ["Species", "gtdb_species"])
        gen_col = _first_present(reader.fieldnames, ["Genus", "gtdb_genus"])
        fam_col = _first_present(reader.fieldnames, ["Family", "gtdb_family"])
        if acc_col is None:
            raise ValueError(
                f"No accession column found in {path} (looked for "
                "GTDB_representative_accession/accession/closest_genome_reference)."
            )

        out: Dict[str, Dict[str, str]] = {}
        for row in reader:
            acc = _normalize_accession((row.get(acc_col) or "").strip())
            if not acc:
                continue
            out[acc] = {
                "species": (row.get(sp_col) or "").strip() if sp_col else "",
                "genus": (row.get(gen_col) or "").strip() if gen_col else "",
                "family": (row.get(fam_col) or "").strip() if fam_col else "",
                "lineage": (row.get(lin_col) or "").strip() if lin_col else "",
            }
    if not out:
        raise ValueError(f"No usable rows in taxonomy table: {path}")
    return out


def _normalize_accession(acc: str) -> str:
    acc = acc.strip()
    if acc.startswith("GB_") or acc.startswith("RS_"):
        acc = acc[3:]
    return acc


def _species_from_lineage(lineage: str) -> str:
    for part in (lineage or "").split(";"):
        part = part.strip()
        if part.startswith("s__"):
            return part[3:].strip()
    return ""


METADATA_HEADERS = [
    "Genome",
    "Lineage",
    "TaxonSource",
    "AssignmentConfidence",
    "Accession",
    "GTDB_species",
    "GTDB_genus",
    "GTDB_family",
    "eggNOG_proteins_annotated",
    "eggNOG_top_max_annot_lvl",
    "eggNOG_top_max_annot_lvl_fraction",
    "eggNOG_top_specific_taxon",
    "eggNOG_top_specific_taxon_fraction",
    "eggNOG_consensus_lineage",
    "AssignmentNote",
]


def build(args: argparse.Namespace) -> int:
    input_dir = args.input.resolve()
    output_dir = args.output.resolve()
    proteins_dir = (args.proteins_dir or (input_dir / "proteins")).resolve()
    eggnog_dir = (args.eggnog_dir or (input_dir / "eggnog")).resolve()

    if not proteins_dir.is_dir():
        raise FileNotFoundError(f"Proteins directory not found: {proteins_dir}")
    if not eggnog_dir.is_dir():
        raise FileNotFoundError(f"eggNOG directory not found: {eggnog_dir}")

    protein_files = find_protein_files(proteins_dir)
    eggnog_files = find_eggnog_files(eggnog_dir)
    taxonomy = load_taxonomy_table(args.taxonomy_table.resolve())

    paired = sorted(set(protein_files) & set(eggnog_files), key=natural_key)
    missing_eggnog = sorted(set(protein_files) - set(eggnog_files), key=natural_key)
    missing_proteins = sorted(set(eggnog_files) - set(protein_files), key=natural_key)
    if not paired:
        raise RuntimeError("No accessions paired between proteins/ and eggnog/.")

    # Remap accessions to clean genome IDs; guard against collisions.
    clean_ids: Dict[str, str] = {}
    seen: Dict[str, str] = {}
    for accession in paired:
        cid = clean_genome_id(accession)
        if cid in seen:
            raise ValueError(
                f"Clean genome ID collision: {accession} and {seen[cid]} both map to {cid}"
            )
        seen[cid] = accession
        clean_ids[accession] = cid

    ensure_output_ready(output_dir, args.force)
    db_dir = output_dir / args.db_name
    func_dir = output_dir / args.func_name
    taxa_path = output_dir / args.taxa_name
    db_dir.mkdir(parents=True, exist_ok=True)
    func_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input:            {input_dir}")
    print(f"Output:           {output_dir}")
    print(f"Protein files:    {len(protein_files)}")
    print(f"eggNOG files:     {len(eggnog_files)}")
    print(f"Paired genomes:   {len(paired)}")
    print(f"Taxonomy entries: {len(taxonomy)}")
    if missing_eggnog:
        print(f"Skipped (no eggNOG): {len(missing_eggnog)}")
    if missing_proteins:
        print(f"Ignored eggNOG without protein FASTA: {len(missing_proteins)}")

    total_proteins = 0
    total_annots = 0
    n_with_lineage = 0
    n_without_lineage = 0

    id_map_path = output_dir / "genome_id_map.tsv"
    protein_map_path = output_dir / "protein_id_map.tsv"
    metadata_rows: List[List[object]] = []

    with id_map_path.open("w", encoding="utf-8", newline="") as id_handle, \
            protein_map_path.open("w", encoding="utf-8", newline="") as pm_handle:
        id_writer = csv.writer(id_handle, delimiter="\t", lineterminator="\n")
        id_writer.writerow(["Accession", "Genome", "GTDB_species"])
        pm_writer = csv.writer(pm_handle, delimiter="\t", lineterminator="\n")
        pm_writer.writerow(["Genome", "Accession", "OldProteinID", "NewProteinID", "ProteinName"])

        for accession in paired:
            genome_id = clean_ids[accession]
            fasta_in = protein_files[accession]
            eggnog_in = eggnog_files[accession]
            fasta_out = db_dir / f"{genome_id}.faa"
            eggnog_out = func_dir / f"{genome_id}_eggNOG.tsv"

            header, rows = parse_eggnog_annotation(eggnog_in)
            consensus = summarize_eggnog_taxonomy(header, rows)
            n_prots, id_map, map_rows = convert_fasta(
                fasta_in, fasta_out, genome_id, header, rows
            )
            n_ann = write_eggnog(eggnog_out, header, rows, id_map)

            for old_id, new_id, protein_name in map_rows:
                pm_writer.writerow([genome_id, accession, old_id, new_id, protein_name])

            tax = taxonomy.get(_normalize_accession(accession), {})
            lineage = tax.get("lineage", "")
            species = tax.get("species", "") or _species_from_lineage(lineage)
            genus = tax.get("genus", "")
            family = tax.get("family", "")
            egg_lineage = str(consensus.get("eggnog_lineage") or "")

            if lineage:
                source = "GTDB_rep"
                confidence = "high"
                note = "GTDB species representative genome; lineage from GTDB species table"
                n_with_lineage += 1
            elif egg_lineage:
                lineage = egg_lineage
                source = "eggNOG"
                confidence = "low"
                note = "Accession not in taxonomy table; lineage from eggNOG protein consensus"
                n_without_lineage += 1
            else:
                lineage = "Unclassified"
                source = "none"
                confidence = "none"
                note = "Accession not in taxonomy table and insufficient eggNOG consensus"
                n_without_lineage += 1

            id_writer.writerow([accession, genome_id, species])
            metadata_rows.append([
                genome_id,
                lineage,
                source,
                confidence,
                accession,
                species,
                genus,
                family,
                consensus.get("n_proteins", 0),
                consensus.get("top_max_annot_lvl", ""),
                consensus.get("top_max_annot_lvl_fraction", 0.0),
                consensus.get("top_specific_taxon", ""),
                consensus.get("top_specific_taxon_fraction", 0.0),
                egg_lineage,
                note,
            ])

            total_proteins += n_prots
            total_annots += n_ann
            print(f"  {accession} -> {genome_id}: proteins={n_prots}, annotated={n_ann}")

    with taxa_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(METADATA_HEADERS)
        for row in sorted(metadata_rows, key=lambda r: natural_key(str(r[0]))):
            writer.writerow(row)

    print("\nDone.")
    print(f"  original_db genomes: {len(paired)}")
    print(f"  eggNOG genomes:      {len(paired)}")
    print(f"  total proteins:      {total_proteins}")
    print(f"  eggNOG rows written: {total_annots}")
    print(f"  lineage from GTDB:   {n_with_lineage}")
    print(f"  lineage fallback:    {n_without_lineage}")
    print(f"  taxonomy file:       {taxa_path}")
    print(f"  genome ID map:       {id_map_path}")
    print(f"  protein ID map:      {protein_map_path}")
    print(
        "\nNext: register this folder in MetaPilot MAG DB settings "
        f"(db={args.db_name}, func={args.func_name}, taxa={args.taxa_name})."
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        return build(args)
    except Exception as exc:  # noqa: BLE001 - CLI boundary
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
