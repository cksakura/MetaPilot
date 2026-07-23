#!/usr/bin/env python3
"""
Build a MetaPilot-compatible MAG database from a binning result folder.

Expected input layout (e.g. SRR*_result):
  proteins/*.faa              Prodigal protein FASTA per bin
  eggnog/*.emapper.annotations
  gtdbtk_retry/classify/gtdbtk.bac120.summary.tsv
  gtdbtk_retry/classify/gtdbtk.ar53.summary.tsv   (optional)

Output layout (MetaPilot local MAG DB version folder):
  original_db/<genomeID>.faa
  eggNOG/<genomeID>_eggNOG.tsv
  genomes-all_metadata.tsv
  protein_id_map.tsv

FASTA headers follow the MGnify / MetaPilot convention:
  >genomeID_00001 30S ribosomal protein S2
  >genomeID_00002 hypothetical protein

Protein accessions are renumbered sequentially (genomeID_00001, ...), not Prodigal
contig IDs. The protein/gene name prefers eggNOG Description, then Preferred_name,
otherwise "hypothetical protein". eggNOG query IDs are rewritten to match.

genomes-all_metadata.tsv always starts with MetaPilot-required columns:
  Genome, Lineage
and adds provenance/confidence columns (TaxonSource, GTDB_*, eggNOG_*, CheckM2_*).
Lineage prefers GTDB-Tk when classified; otherwise an eggNOG protein-consensus
lineage is used so every genome has a taxon assignment.

Only genomes with both a protein FASTA and an eggNOG annotation file are written
(required for MetaPilot pairing validation). Marker2/Marker4 are generated later
by MetaPilot during database registration.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


PRODIGAL_HEADER_RE = re.compile(r"^>(\S+)(?:\s+(.*))?$")
EGGNOG_ANNOT_SUFFIX = ".emapper.annotations"
DEFAULT_PROTEIN_NAME = "hypothetical protein"


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a binning/eggNOG/GTDB-Tk result folder into a MetaPilot MAG database."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Binning result root (contains proteins/, eggnog/, gtdbtk_retry/).",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output MAG DB version folder (will contain original_db/, eggNOG/, genomes-all_metadata.tsv).",
    )
    parser.add_argument(
        "--proteins-dir",
        type=Path,
        default=None,
        help="Override proteins directory (default: <input>/proteins).",
    )
    parser.add_argument(
        "--eggnog-dir",
        type=Path,
        default=None,
        help="Override eggNOG directory (default: <input>/eggnog).",
    )
    parser.add_argument(
        "--gtdbtk-dir",
        type=Path,
        default=None,
        help="Override GTDB-Tk directory (default: <input>/gtdbtk_retry, else <input>/gtdbtk).",
    )
    parser.add_argument(
        "--db-name",
        default="original_db",
        help="Protein FASTA folder name (default: original_db).",
    )
    parser.add_argument(
        "--func-name",
        default="eggNOG",
        help="Functional annotation folder name (default: eggNOG).",
    )
    parser.add_argument(
        "--taxa-name",
        default="genomes-all_metadata.tsv",
        help="Taxonomy metadata filename (default: genomes-all_metadata.tsv).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args(argv)


def find_protein_files(proteins_dir: Path) -> Dict[str, Path]:
    files: Dict[str, Path] = {}
    for path in sorted(proteins_dir.glob("*.faa")):
        files[path.stem] = path
    for path in sorted(proteins_dir.glob("*.fasta")):
        files.setdefault(path.stem, path)
    return files


def find_eggnog_files(eggnog_dir: Path) -> Dict[str, Path]:
    files: Dict[str, Path] = {}
    for path in sorted(eggnog_dir.glob(f"*{EGGNOG_ANNOT_SUFFIX}")):
        genome_id = path.name[: -len(EGGNOG_ANNOT_SUFFIX)]
        files[genome_id] = path
    return files


def resolve_gtdbtk_dir(input_dir: Path, override: Optional[Path]) -> Path:
    if override is not None:
        return override
    retry = input_dir / "gtdbtk_retry"
    if retry.is_dir():
        return retry
    plain = input_dir / "gtdbtk"
    if plain.is_dir():
        return plain
    raise FileNotFoundError(
        f"No GTDB-Tk folder found under {input_dir} (expected gtdbtk_retry/ or gtdbtk/)."
    )


def find_summary_file(gtdbtk_dir: Path, name: str) -> Optional[Path]:
    candidates = [
        gtdbtk_dir / "classify" / name,
        gtdbtk_dir / name,
    ]
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            return path
    return None


def load_gtdb_rows(gtdbtk_dir: Path) -> Dict[str, Dict[str, str]]:
    """Load full GTDB-Tk summary rows keyed by genome ID."""
    rows: Dict[str, Dict[str, str]] = {}
    for summary_name in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
        path = find_summary_file(gtdbtk_dir, summary_name)
        if path is None:
            continue
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "user_genome" not in reader.fieldnames:
                raise ValueError(f"Unexpected GTDB-Tk summary format: {path}")
            if "classification" not in reader.fieldnames:
                raise ValueError(f"Missing classification column in {path}")
            for row in reader:
                genome = (row.get("user_genome") or "").strip()
                if genome:
                    rows[genome] = {k: (v if v is not None else "") for k, v in row.items()}
    if not rows:
        raise FileNotFoundError(
            f"No usable GTDB-Tk summary TSV found under {gtdbtk_dir}/classify/"
        )
    return rows


def load_quality_report(input_dir: Path) -> Dict[str, Dict[str, str]]:
    path = input_dir / "quality_report.tsv"
    if not path.is_file():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        out: Dict[str, Dict[str, str]] = {}
        for row in reader:
            name = (row.get("Name") or "").strip()
            if name:
                out[name] = row
        return out


# NCBI-ish eggNOG taxon names -> GTDB-style rank tokens used in MetaPilot Lineage parsing.
EGGNOG_NAME_TO_GTDB = {
    "Bacteria": "d__Bacteria",
    "Archaea": "d__Archaea",
    "Eukaryota": "d__Eukaryota",
    "Firmicutes": "p__Bacillota",
    "Bacillota": "p__Bacillota",
    "Bacteroidetes": "p__Bacteroidota",
    "Bacteroidota": "p__Bacteroidota",
    "Actinobacteria": "p__Actinomycetota",
    "Actinomycetota": "p__Actinomycetota",
    "Proteobacteria": "p__Pseudomonadota",
    "Pseudomonadota": "p__Pseudomonadota",
    "Fusobacteria": "p__Fusobacteriota",
    "Verrucomicrobia": "p__Verrucomicrobiota",
    "Clostridia": "c__Clostridia",
    "Bacilli": "c__Bacilli",
    "Bacteroidia": "c__Bacteroidia",
    "Coriobacteriia": "c__Coriobacteriia",
    "Erysipelotrichia": "c__Erysipelotrichia",
    "Negativicutes": "c__Negativicutes",
    "Alphaproteobacteria": "c__Alphaproteobacteria",
    "Gammaproteobacteria": "c__Gammaproteobacteria",
    "Clostridiales": "o__Clostridiales",
    "Oscillospirales": "o__Oscillospirales",
    "Lachnospirales": "o__Lachnospirales",
    "Eubacteriales": "o__Eubacteriales",
    "Bacteroidales": "o__Bacteroidales",
    "Coriobacteriales": "o__Coriobacteriales",
    "Erysipelotrichales": "o__Erysipelotrichales",
    "Ruminococcaceae": "f__Ruminococcaceae",
    "Oscillospiraceae": "f__Oscillospiraceae",
    "Lachnospiraceae": "f__Lachnospiraceae",
    "Clostridiaceae": "f__Clostridiaceae",
    "Eubacteriaceae": "f__Eubacteriaceae",
    "Peptostreptococcaceae": "f__Peptostreptococcaceae",
    "Enterobacteriaceae": "f__Enterobacteriaceae",
    "Bacteroidaceae": "f__Bacteroidaceae",
    "Prevotellaceae": "f__Prevotellaceae",
    "Eggerthellaceae": "f__Eggerthellaceae",
    "Coriobacteriaceae": "f__Coriobacteriaceae",
    "Blautia": "g__Blautia",
    "Faecalibacterium": "g__Faecalibacterium",
    "Bacteroides": "g__Bacteroides",
    "Prevotella": "g__Prevotella",
    "Collinsella": "g__Collinsella",
}


def _col_index(header: List[str], name: str) -> int:
    try:
        return header.index(name)
    except ValueError:
        return -1


def _field(fields: List[str], idx: int) -> str:
    if idx < 0 or idx >= len(fields):
        return ""
    value = fields[idx].strip()
    return "" if value in {"", "-"} else value


def parse_eggnog_og_taxa(ogs: str) -> List[Tuple[str, str]]:
    """Return [(taxid, taxname), ...] from broad to specific."""
    out: List[Tuple[str, str]] = []
    if not ogs:
        return out
    for tok in ogs.split(","):
        tok = tok.strip()
        if "@" not in tok or "|" not in tok:
            continue
        left, taxname = tok.split("|", 1)
        _nog, taxid = left.split("@", 1)
        taxname = taxname.strip()
        if not taxname or taxname.lower() == "root":
            continue
        # Skip uninformative "unclassified X" as lineage anchors when better options exist,
        # but still keep them for reporting.
        out.append((taxid, taxname))
    return out


def summarize_eggnog_taxonomy(
    eggnog_header: List[str], eggnog_rows: Dict[str, List[str]]
) -> Dict[str, object]:
    """Protein-level eggNOG taxonomy consensus for one genome."""
    idx_ogs = _col_index(eggnog_header, "eggNOG_OGs")
    idx_max = _col_index(eggnog_header, "max_annot_lvl")
    idx_seed = _col_index(eggnog_header, "seed_ortholog")

    n_proteins = len(eggnog_rows)
    max_lvl: Counter = Counter()
    specific: Counter = Counter()
    seed_tax: Counter = Counter()
    domain_votes: Counter = Counter()

    junk_taxnames = {
        "delta/epsilon subdivisions",
        "environmental samples",
        "cellular organisms",
        "unclassified sequences",
    }

    for fields in eggnog_rows.values():
        mx = _field(fields, idx_max)
        if mx:
            max_lvl[mx] += 1
        ogs = _field(fields, idx_ogs)
        taxa = parse_eggnog_og_taxa(ogs)
        if taxa:
            for _taxid, taxname in taxa:
                if taxname in {"Bacteria", "Archaea", "Eukaryota"}:
                    domain_votes[EGGNOG_NAME_TO_GTDB[taxname]] += 1
            sp_taxid, sp_name = taxa[-1]
            if (
                sp_name not in junk_taxnames
                and not sp_name.startswith("unclassified ")
            ):
                specific[f"{sp_taxid}|{sp_name}"] += 1
        seed = _field(fields, idx_seed)
        if seed and "." in seed:
            taxid = seed.split(".", 1)[0]
            if taxid.isdigit():
                seed_tax[taxid] += 1

    def top_with_frac(counter: Counter) -> Tuple[str, int, float]:
        if not counter or n_proteins == 0:
            return "", 0, 0.0
        label, count = counter.most_common(1)[0]
        return label, count, count / float(n_proteins)

    top_max, top_max_n, top_max_frac = top_with_frac(max_lvl)
    top_sp, top_sp_n, top_sp_frac = top_with_frac(specific)
    top_seed, top_seed_n, top_seed_frac = top_with_frac(seed_tax)
    top2 = ""
    if len(max_lvl) >= 2:
        top2 = ";".join(f"{lab}:{cnt}" for lab, cnt in max_lvl.most_common(2))

    # Conservative lineage: domain + mapped top max_annot_lvl (+ optional family/genus
    # from top specific taxon only when support is decent). Do NOT mix independent
    # rank majority votes — that creates chimeric lineages on LQ bins.
    lineage_parts: List[str] = []
    if domain_votes:
        lineage_parts.append(domain_votes.most_common(1)[0][0])
    elif top_max:
        lineage_parts.append("d__Bacteria")

    top_max_token = ""
    if top_max and "|" in top_max:
        max_name = top_max.split("|", 1)[1]
        top_max_token = EGGNOG_NAME_TO_GTDB.get(max_name, "")
        if top_max_token and top_max_token not in lineage_parts:
            lineage_parts.append(top_max_token)

    if top_sp and "|" in top_sp and top_sp_frac >= 0.20 and top_max_frac >= 0.50:
        sp_name = top_sp.split("|", 1)[1]
        sp_token = EGGNOG_NAME_TO_GTDB.get(sp_name, "")
        if sp_token and sp_token not in lineage_parts:
            # Only deepen with family/genus tokens; ignore noisy species strings.
            if sp_token.startswith(("f__", "g__")):
                lineage_parts.append(sp_token)

    return {
        "n_proteins": n_proteins,
        "top_max_annot_lvl": top_max,
        "top_max_annot_lvl_count": top_max_n,
        "top_max_annot_lvl_fraction": round(top_max_frac, 4),
        "top2_max_annot_lvl": top2,
        "top_specific_taxon": top_sp,
        "top_specific_taxon_count": top_sp_n,
        "top_specific_taxon_fraction": round(top_sp_frac, 4),
        "top_seed_taxid": top_seed,
        "top_seed_taxid_count": top_seed_n,
        "top_seed_taxid_fraction": round(top_seed_frac, 4),
        "eggnog_lineage": ";".join(lineage_parts),
    }


def is_gtdb_classified(classification: str) -> bool:
    c = (classification or "").strip()
    return bool(c) and not c.startswith("Unclassified")


def assignment_confidence(
    source: str,
    gtdb_method: str,
    gtdb_warnings: str,
    max_frac: float,
    specific_frac: float,
    completeness: Optional[float],
) -> str:
    if source == "GTDB":
        if gtdb_warnings and gtdb_warnings not in {"", "N/A", "nan"}:
            return "medium"
        if gtdb_method and gtdb_method not in {"", "N/A"}:
            return "high"
        return "high"
    if source == "eggNOG":
        if max_frac >= 0.7 and specific_frac >= 0.25:
            return "medium"
        if max_frac >= 0.5:
            return "low"
        return "low"
    if completeness is not None and completeness < 50:
        return "low"
    return "none"


def write_taxonomy_metadata(
    out_path: Path,
    genome_ids: Iterable[str],
    gtdb_rows: Dict[str, Dict[str, str]],
    eggnog_by_genome: Dict[str, Dict[str, object]],
    quality_by_genome: Dict[str, Dict[str, str]],
) -> Tuple[int, int, int]:
    """
    Write genomes-all_metadata.tsv.

    First two columns remain MetaPilot-required:
      Genome, Lineage
    Extra columns document provenance and confidence.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "Genome",
        "Lineage",
        "TaxonSource",
        "AssignmentConfidence",
        "GTDB_classification",
        "GTDB_method",
        "GTDB_warnings",
        "GTDB_closest_reference",
        "GTDB_closest_taxonomy",
        "GTDB_closest_ani",
        "GTDB_closest_af",
        "GTDB_msa_percent",
        "eggNOG_proteins_annotated",
        "eggNOG_top_max_annot_lvl",
        "eggNOG_top_max_annot_lvl_count",
        "eggNOG_top_max_annot_lvl_fraction",
        "eggNOG_top2_max_annot_lvl",
        "eggNOG_top_specific_taxon",
        "eggNOG_top_specific_taxon_count",
        "eggNOG_top_specific_taxon_fraction",
        "eggNOG_top_seed_taxid",
        "eggNOG_top_seed_taxid_count",
        "eggNOG_top_seed_taxid_fraction",
        "eggNOG_consensus_lineage",
        "CheckM2_completeness",
        "CheckM2_contamination",
        "CheckM2_genome_size",
        "AssignmentNote",
    ]

    n_gtdb = 0
    n_eggnog = 0
    n_written = 0

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)

        for genome_id in sorted(genome_ids, key=natural_key):
            gtdb = gtdb_rows.get(genome_id, {})
            egg = eggnog_by_genome.get(genome_id, {})
            qual = quality_by_genome.get(genome_id, {})

            gtdb_class = (gtdb.get("classification") or "").strip() or "Unclassified"
            gtdb_method = (gtdb.get("classification_method") or "").strip()
            gtdb_warnings = (gtdb.get("warnings") or "").strip()
            closest_ref = (gtdb.get("closest_genome_reference") or "").strip()
            closest_tax = (gtdb.get("closest_genome_taxonomy") or "").strip()
            closest_ani = (gtdb.get("closest_genome_ani") or "").strip()
            closest_af = (gtdb.get("closest_genome_af") or "").strip()
            msa = (gtdb.get("msa_percent") or "").strip()

            max_frac = float(egg.get("top_max_annot_lvl_fraction") or 0.0)
            specific_frac = float(egg.get("top_specific_taxon_fraction") or 0.0)
            egg_lineage = str(egg.get("eggnog_lineage") or "")

            completeness: Optional[float] = None
            contam = ""
            gsize = ""
            if qual:
                try:
                    completeness = float(qual.get("Completeness") or "")
                except ValueError:
                    completeness = None
                contam = str(qual.get("Contamination") or "")
                gsize = str(qual.get("Genome_Size") or "")

            note = ""
            if is_gtdb_classified(gtdb_class):
                lineage = gtdb_class
                source = "GTDB"
                n_gtdb += 1
                note = "Assigned from GTDB-Tk classification"
            elif egg_lineage:
                lineage = egg_lineage
                source = "eggNOG"
                n_eggnog += 1
                note = (
                    "GTDB unclassified/failed; Lineage inferred from eggNOG protein consensus "
                    f"(top max_annot_lvl fraction={max_frac:.3f})"
                )
            else:
                lineage = gtdb_class if gtdb_class else "Unclassified"
                source = "none"
                note = "No GTDB classification and insufficient eggNOG consensus"

            conf = assignment_confidence(
                source, gtdb_method, gtdb_warnings, max_frac, specific_frac, completeness
            )

            writer.writerow(
                [
                    genome_id,
                    lineage,
                    source,
                    conf,
                    gtdb_class,
                    gtdb_method,
                    gtdb_warnings,
                    closest_ref,
                    closest_tax,
                    closest_ani,
                    closest_af,
                    msa,
                    egg.get("n_proteins", 0),
                    egg.get("top_max_annot_lvl", ""),
                    egg.get("top_max_annot_lvl_count", 0),
                    egg.get("top_max_annot_lvl_fraction", 0.0),
                    egg.get("top2_max_annot_lvl", ""),
                    egg.get("top_specific_taxon", ""),
                    egg.get("top_specific_taxon_count", 0),
                    egg.get("top_specific_taxon_fraction", 0.0),
                    egg.get("top_seed_taxid", ""),
                    egg.get("top_seed_taxid_count", 0),
                    egg.get("top_seed_taxid_fraction", 0.0),
                    egg_lineage,
                    "" if completeness is None else completeness,
                    contam,
                    gsize,
                    note,
                ]
            )
            n_written += 1

    return n_written, n_gtdb, n_eggnog


def parse_eggnog_annotation(path: Path) -> Tuple[List[str], Dict[str, List[str]]]:
    """
    Return (header_columns, query -> row fields).
    Accepts emapper ## comment lines and a #query header line.
    """
    header: Optional[List[str]] = None
    rows: Dict[str, List[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].split("\t")
                continue
            fields = line.split("\t")
            if header is None:
                # Fallback if no header was present
                header = [f"col{i}" for i in range(len(fields))]
                header[0] = "query"
            if not fields or not fields[0].strip():
                continue
            rows[fields[0].strip()] = fields
    if header is None:
        raise ValueError(f"No header/data found in eggNOG file: {path}")
    # Normalise first column name used by MetaPilot / emapper
    if header and header[0] in {"query", "query_name"}:
        header[0] = "query_name"
    return header, rows


def protein_display_name(eggnog_row: Optional[List[str]], header: List[str]) -> str:
    """
    MGnify-style product name after the accession, e.g.:
      >MGYG000321000_00003 30S ribosomal protein S2
      >MGYG000321000_00001 hypothetical protein
    """
    if eggnog_row is None:
        return DEFAULT_PROTEIN_NAME

    def col(name: str) -> str:
        try:
            idx = header.index(name)
        except ValueError:
            return ""
        if idx >= len(eggnog_row):
            return ""
        value = eggnog_row[idx].strip()
        return "" if value in {"", "-"} else value

    # Prefer functional product description; fall back to gene symbol; else hypothetical.
    description = col("Description")
    preferred = col("Preferred_name")
    name = description or preferred or DEFAULT_PROTEIN_NAME
    name = re.sub(r"\s+", " ", name.replace("\t", " ")).strip()
    return name or DEFAULT_PROTEIN_NAME


def clean_sequence(seq: str) -> str:
    # Search engines generally reject stop codon asterisks in protein FASTA.
    return seq.replace("*", "").replace(" ", "").replace("\t", "")


def convert_fasta(
    fasta_path: Path,
    out_path: Path,
    genome_id: str,
    eggnog_header: List[str],
    eggnog_rows: Dict[str, List[str]],
) -> Tuple[int, Dict[str, str], List[Tuple[str, str, str]]]:
    """
    Rewrite FASTA to MetaPilot/MGnify format.

    Returns:
      protein_count,
      old_query_id -> new_accession (genomeID_00001),
      mapping rows (old_id, new_id, protein_name)
    """
    id_map: Dict[str, str] = {}
    map_rows: List[Tuple[str, str, str]] = []
    protein_count = 0
    current_old_id: Optional[str] = None
    current_name = DEFAULT_PROTEIN_NAME
    seq_chunks: List[str] = []

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with fasta_path.open("r", encoding="utf-8", errors="replace") as in_handle, out_path.open(
        "w", encoding="utf-8", newline="\n"
    ) as out_handle:

        def flush() -> None:
            nonlocal protein_count
            if current_old_id is None:
                return
            sequence = clean_sequence("".join(seq_chunks))
            if not sequence:
                return
            protein_count += 1
            new_id = f"{genome_id}_{protein_count:05d}"
            # MetaPilot rule: >genomeID_proteinID protein/gene name
            out_handle.write(f">{new_id} {current_name}\n")
            for i in range(0, len(sequence), 60):
                out_handle.write(sequence[i : i + 60] + "\n")
            id_map[current_old_id] = new_id
            map_rows.append((current_old_id, new_id, current_name))

        for raw in in_handle:
            line = raw.rstrip("\n\r")
            if line.startswith(">"):
                flush()
                seq_chunks = []
                match = PRODIGAL_HEADER_RE.match(line)
                if not match:
                    raise ValueError(f"Unrecognized FASTA header in {fasta_path}: {line}")
                old_id = match.group(1)
                if not old_id:
                    raise ValueError(f"Empty protein ID in {fasta_path}: {line}")
                current_old_id = old_id
                current_name = protein_display_name(eggnog_rows.get(old_id), eggnog_header)
            elif line:
                seq_chunks.append(line.strip())
        flush()

    if protein_count == 0:
        raise ValueError(f"No proteins written from {fasta_path}")
    return protein_count, id_map, map_rows


def write_eggnog(
    out_path: Path,
    header: List[str],
    eggnog_rows: Dict[str, List[str]],
    id_map: Dict[str, str],
) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        # MetaPilot MagFuncSearcher skips the first line as header.
        writer.writerow(header)
        for old_id, new_id in id_map.items():
            row = eggnog_rows.get(old_id)
            if row is None:
                continue
            fields = list(row)
            # Pad/truncate to header width for stable column indices
            if len(fields) < len(header):
                fields.extend(["-"] * (len(header) - len(fields)))
            elif len(fields) > len(header):
                fields = fields[: len(header)]
            fields[0] = new_id
            writer.writerow(fields)
            written += 1
    return written


def natural_key(value: str) -> Tuple:
    parts = re.split(r"(\d+)", value)
    key: List[object] = []
    for part in parts:
        if part.isdigit():
            key.append(int(part))
        else:
            key.append(part.lower())
    return tuple(key)


def ensure_output_ready(output_dir: Path, force: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()) and not force:
        raise FileExistsError(
            f"Output directory is not empty: {output_dir}\n"
            "Use --force to overwrite, or choose an empty folder."
        )
    output_dir.mkdir(parents=True, exist_ok=True)


def build_mag_db(args: argparse.Namespace) -> int:
    input_dir = args.input.resolve()
    output_dir = args.output.resolve()
    proteins_dir = (args.proteins_dir or (input_dir / "proteins")).resolve()
    eggnog_dir = (args.eggnog_dir or (input_dir / "eggnog")).resolve()
    gtdbtk_dir = resolve_gtdbtk_dir(input_dir, args.gtdbtk_dir).resolve()

    if not proteins_dir.is_dir():
        raise FileNotFoundError(f"Proteins directory not found: {proteins_dir}")
    if not eggnog_dir.is_dir():
        raise FileNotFoundError(f"eggNOG directory not found: {eggnog_dir}")

    protein_files = find_protein_files(proteins_dir)
    eggnog_files = find_eggnog_files(eggnog_dir)
    gtdb_rows = load_gtdb_rows(gtdbtk_dir)
    quality_rows = load_quality_report(input_dir)

    paired = sorted(set(protein_files) & set(eggnog_files), key=natural_key)
    missing_eggnog = sorted(set(protein_files) - set(eggnog_files), key=natural_key)
    missing_proteins = sorted(set(eggnog_files) - set(protein_files), key=natural_key)

    if not paired:
        raise RuntimeError("No genome IDs paired between proteins/ and eggnog/.")

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
    print(f"GTDB entries:     {len(gtdb_rows)}")
    print(f"Quality entries:  {len(quality_rows)}")
    if missing_eggnog:
        print(f"Skipped (no eggNOG): {len(missing_eggnog)}")
        print("  " + ", ".join(missing_eggnog))
    if missing_proteins:
        print(f"Ignored eggNOG without protein FASTA: {len(missing_proteins)}")

    converted: List[str] = []
    total_proteins = 0
    total_annots = 0
    eggnog_consensus: Dict[str, Dict[str, object]] = {}
    map_path = output_dir / "protein_id_map.tsv"

    with map_path.open("w", encoding="utf-8", newline="") as map_handle:
        map_writer = csv.writer(map_handle, delimiter="\t", lineterminator="\n")
        map_writer.writerow(["Genome", "OldProteinID", "NewProteinID", "ProteinName"])

        for genome_id in paired:
            fasta_in = protein_files[genome_id]
            eggnog_in = eggnog_files[genome_id]
            fasta_out = db_dir / f"{genome_id}.faa"
            eggnog_out = func_dir / f"{genome_id}_eggNOG.tsv"

            header, rows = parse_eggnog_annotation(eggnog_in)
            eggnog_consensus[genome_id] = summarize_eggnog_taxonomy(header, rows)
            n_prots, id_map, map_rows = convert_fasta(
                fasta_in, fasta_out, genome_id, header, rows
            )
            n_ann = write_eggnog(eggnog_out, header, rows, id_map)
            for old_id, new_id, protein_name in map_rows:
                map_writer.writerow([genome_id, old_id, new_id, protein_name])

            converted.append(genome_id)
            total_proteins += n_prots
            total_annots += n_ann
            print(f"  {genome_id}: proteins={n_prots}, annotated={n_ann}")

    taxa_written, n_gtdb, n_eggnog = write_taxonomy_metadata(
        taxa_path, converted, gtdb_rows, eggnog_consensus, quality_rows
    )

    print("\nDone.")
    print(f"  original_db genomes: {len(converted)}")
    print(f"  eggNOG genomes:      {len(converted)}")
    print(f"  total proteins:      {total_proteins}")
    print(f"  eggNOG rows written: {total_annots}")
    print(f"  taxonomy rows:       {taxa_written}")
    print(f"    TaxonSource=GTDB:  {n_gtdb}")
    print(f"    TaxonSource=eggNOG:{n_eggnog}")
    print(f"    TaxonSource=none:  {taxa_written - n_gtdb - n_eggnog}")
    print(f"  taxonomy file:       {taxa_path}")
    print(f"  protein ID map:      {map_path}")
    print(
        "\nNext: register this folder in MetaPilot MAG DB settings "
        f"(db={args.db_name}, func={args.func_name}, taxa={args.taxa_name}). "
        "MetaPilot uses Genome + Lineage (columns 1-2); extra columns are provenance."
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        return build_mag_db(args)
    except Exception as exc:  # noqa: BLE001 - CLI boundary
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
