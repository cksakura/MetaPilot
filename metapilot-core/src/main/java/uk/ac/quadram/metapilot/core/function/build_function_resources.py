#!/usr/bin/env python3
"""
Build MetaPilot functional annotation resources.

This script is the maintained replacement for the older one-off download
helpers in this directory. It downloads/parses all supported public functional
definition resources except the full eggNOG OG annotation table, which is very
large and should be supplied manually with --eggnog-og-annotations when needed.

Main outputs under --outdir:
    annotation_lookup_merged.tsv
    kegg/*.tsv
    go/GO.tsv
    go/GO_full.tsv
    ec/EC.tsv
    ec/EC_full.tsv
    cog/COG.tsv
    cog/COG_long.tsv
    cog/COG_category.tsv
    pfam/PFAM.tsv
    bigg/BiGG_Reaction.tsv
    cazy/CAZy.tsv
    kegg_tc/KEGG_TC.tsv

Optional output:
    func_def.db, when --sqlite-out is provided.

Example:
    python build_function_resources.py --outdir annotation_maps --sqlite-out Resources/function/func_def.db

With manually downloaded eggNOG OG annotations:
    python build_function_resources.py --outdir annotation_maps --sqlite-out Resources/function/func_def.db \
        --eggnog-og-annotations C:/path/to/e5.og_annotations.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


KEGG_TABLES = {
    "KEGG_ko": ("https://rest.kegg.jp/list/ko", "ko:"),
    "KEGG_Pathway": ("https://rest.kegg.jp/list/pathway", "path:"),
    "KEGG_Module": ("https://rest.kegg.jp/list/module", "md:"),
    "KEGG_Reaction": ("https://rest.kegg.jp/list/reaction", "rn:"),
    "KEGG_rclass": ("https://rest.kegg.jp/list/rclass", "rc:"),
    "KEGG_EC": ("https://rest.kegg.jp/list/enzyme", "ec:"),
    "BRITE": ("https://rest.kegg.jp/list/brite", "br:"),
}

GO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"
BIGG_REACTIONS_URL = "http://bigg.ucsd.edu/api/v2/universal/reactions"
COG_DEF_URL = "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab"
COG_FUN_URL = "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.fun.tab"
ENZYME_DAT_URL = "https://ftp.expasy.org/databases/enzyme/enzyme.dat"


@dataclass(frozen=True)
class Term:
    source: str
    term_id: str
    name: str = ""
    description: str = ""
    namespace: str = ""
    category: str = ""


def log(message: str) -> None:
    print(message, flush=True)


def download_binary(url: str, retries: int = 3, sleep_seconds: int = 5, timeout: int = 300) -> bytes:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "metapilot-function-resource-builder/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_error = e
            log(f"  attempt {attempt}/{retries} failed: {e}")
            if attempt < retries:
                time.sleep(sleep_seconds)
    raise RuntimeError(f"Failed to download {url}: {last_error}")


def download_text(url: str, retries: int = 3, sleep_seconds: int = 5, timeout: int = 300) -> str:
    return download_binary(url, retries=retries, sleep_seconds=sleep_seconds, timeout=timeout).decode(
        "utf-8", errors="replace"
    )


def write_simple_tsv(path: Path, rows: Iterable[tuple[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tdescription\n")
        for term_id, description in rows:
            term_id = clean_cell(term_id)
            description = clean_cell(description)
            if not term_id:
                continue
            f.write(f"{term_id}\t{description}\n")
            count += 1
    return count


def clean_cell(value: object) -> str:
    return "" if value is None else str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def download_kegg(outdir: Path, force: bool = False) -> None:
    kegg_dir = outdir / "kegg"
    kegg_dir.mkdir(parents=True, exist_ok=True)

    for table_name, (url, prefix) in KEGG_TABLES.items():
        filename = "EC.tsv" if table_name == "KEGG_EC" else f"{table_name}.tsv"
        outfile = kegg_dir / filename
        if outfile.exists() and outfile.stat().st_size > 0 and not force:
            log(f"Using existing KEGG table -> {outfile}")
            continue

        log(f"Downloading {table_name} from KEGG...")
        rows: list[tuple[str, str]] = []
        for line in download_text(url, timeout=180).splitlines():
            if not line.strip():
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            term_id, description = parts
            if term_id.startswith(prefix):
                term_id = term_id[len(prefix):]
            rows.append((term_id, description))
        count = write_simple_tsv(outfile, rows)
        log(f"  wrote {count:,} rows -> {outfile}")


def parse_go_obo(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    current: dict[str, str] = {}
    in_term = False

    def flush() -> None:
        nonlocal current, in_term
        if in_term and current.get("id") and current.get("name") and current.get("is_obsolete") != "true":
            rows.append(
                {
                    "id": current.get("id", ""),
                    "name": current.get("name", ""),
                    "namespace": current.get("namespace", ""),
                }
            )
        current = {}
        in_term = False

    for raw in text.splitlines():
        line = raw.strip()
        if line == "[Term]":
            flush()
            in_term = True
        elif line.startswith("["):
            flush()
        elif not line:
            flush()
        elif in_term and ": " in line:
            key, value = line.split(": ", 1)
            if key in {"id", "name", "namespace", "is_obsolete"}:
                current[key] = value
    flush()
    return rows


def download_go(outdir: Path, force: bool = False) -> None:
    go_dir = outdir / "go"
    go_dir.mkdir(parents=True, exist_ok=True)
    raw_file = go_dir / "go-basic.obo"
    full_file = go_dir / "GO_full.tsv"
    simple_file = go_dir / "GO.tsv"

    if raw_file.exists() and raw_file.stat().st_size > 0 and not force:
        text = raw_file.read_text(encoding="utf-8", errors="replace")
        log(f"Using existing GO OBO -> {raw_file}")
    else:
        log("Downloading GO ontology...")
        text = download_text(GO_URL, timeout=300)
        raw_file.write_text(text, encoding="utf-8")

    rows = parse_go_obo(text)
    with full_file.open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tname\tnamespace\n")
        for row in rows:
            f.write(f"{row['id']}\t{clean_cell(row['name'])}\t{clean_cell(row['namespace'])}\n")
    write_simple_tsv(simple_file, ((row["id"], row["name"]) for row in rows))
    log(f"  wrote {len(rows):,} GO rows -> {simple_file}")


def parse_enzyme_dat(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    current: dict[str, object] | None = None

    def flush() -> None:
        nonlocal current
        if not current or not current.get("id"):
            current = None
            return
        rows.append(
            {
                "id": clean_cell(current.get("id", "")),
                "description": clean_cell(current.get("description", "")),
                "reaction": clean_cell(" ".join(current.get("reaction", []))),
                "synonyms": clean_cell(";".join(current.get("synonyms", []))),
            }
        )
        current = None

    for raw in text.splitlines():
        if raw.startswith("//"):
            flush()
            continue
        if len(raw) < 5:
            continue
        key = raw[:2]
        value = raw[5:].strip() if raw[2:5].strip() == "" else raw[2:].strip()
        if key == "ID":
            flush()
            current = {"id": value, "reaction": [], "synonyms": []}
        elif current is not None and key == "DE":
            current["description"] = value
        elif current is not None and key == "CA":
            current.setdefault("reaction", []).append(value)
        elif current is not None and key == "AN":
            current.setdefault("synonyms", []).append(value.rstrip("."))
    flush()
    return rows


def download_ec(outdir: Path, force: bool = False) -> None:
    ec_dir = outdir / "ec"
    ec_dir.mkdir(parents=True, exist_ok=True)
    raw_file = ec_dir / "enzyme.dat"
    full_file = ec_dir / "EC_full.tsv"
    simple_file = ec_dir / "EC.tsv"

    if raw_file.exists() and raw_file.stat().st_size > 0 and not force:
        text = raw_file.read_text(encoding="utf-8", errors="replace")
        log(f"Using existing enzyme.dat -> {raw_file}")
    else:
        log("Downloading ExPASy enzyme.dat...")
        text = download_text(ENZYME_DAT_URL, timeout=300)
        raw_file.write_text(text, encoding="utf-8")

    rows = parse_enzyme_dat(text)
    with full_file.open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tdescription\treaction\tsynonyms\n")
        for row in rows:
            f.write(f"{row['id']}\t{row['description']}\t{row['reaction']}\t{row['synonyms']}\n")
    write_simple_tsv(simple_file, ((row["id"], row["description"]) for row in rows))
    log(f"  wrote {len(rows):,} EC rows -> {full_file}")


def parse_pfam_dat(text: str) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    pfam_id = ""
    pfam_acc = ""
    description = ""
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if line.startswith("#=GF ID"):
            pfam_id = line.split(maxsplit=2)[2].strip()
        elif line.startswith("#=GF AC"):
            pfam_acc = line.split(maxsplit=2)[2].strip().split(".")[0]
        elif line.startswith("#=GF DE"):
            description = line.split(maxsplit=2)[2].strip()
        elif line == "//":
            if description:
                if pfam_id:
                    rows.append((pfam_id, description))
                if pfam_acc:
                    rows.append((pfam_acc, description))
            pfam_id = pfam_acc = description = ""
    return rows


def download_pfam(outdir: Path, force: bool = False) -> None:
    pfam_dir = outdir / "pfam"
    pfam_dir.mkdir(parents=True, exist_ok=True)
    gz_file = pfam_dir / "Pfam-A.hmm.dat.gz"
    dat_file = pfam_dir / "Pfam-A.hmm.dat"
    table_file = pfam_dir / "PFAM.tsv"

    if gz_file.exists() and gz_file.stat().st_size > 0 and not force:
        binary = gz_file.read_bytes()
        log(f"Using existing Pfam gzip -> {gz_file}")
    else:
        log("Downloading Pfam metadata...")
        binary = download_binary(PFAM_URL, timeout=600)
        gz_file.write_bytes(binary)

    text = gzip.decompress(binary).decode("utf-8", errors="replace")
    dat_file.write_text(text, encoding="utf-8")
    rows = parse_pfam_dat(text)
    count = write_simple_tsv(table_file, rows)
    log(f"  wrote {count:,} PFAM rows -> {table_file}")


def download_bigg(outdir: Path, force: bool = False) -> None:
    bigg_dir = outdir / "bigg"
    bigg_dir.mkdir(parents=True, exist_ok=True)
    json_file = bigg_dir / "bigg_universal_reactions.json"
    table_file = bigg_dir / "BiGG_Reaction.tsv"

    if json_file.exists() and json_file.stat().st_size > 0 and not force:
        text = json_file.read_text(encoding="utf-8", errors="replace")
        log(f"Using existing BiGG JSON -> {json_file}")
    else:
        log("Downloading BiGG universal reactions...")
        text = download_text(BIGG_REACTIONS_URL, timeout=300)
        json_file.write_text(text, encoding="utf-8")

    data = json.loads(text)
    rows = [(item.get("bigg_id", ""), item.get("name", "")) for item in data.get("results", [])]
    count = write_simple_tsv(table_file, rows)
    log(f"  wrote {count:,} BiGG rows -> {table_file}")


def download_cog(outdir: Path, force: bool = False) -> None:
    cog_dir = outdir / "cog"
    raw_dir = cog_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    def_file = raw_dir / "cog-24.def.tab"
    fun_file = raw_dir / "cog-24.fun.tab"

    if def_file.exists() and def_file.stat().st_size > 0 and not force:
        def_text = def_file.read_text(encoding="utf-8", errors="replace")
        log(f"Using existing COG definitions -> {def_file}")
    else:
        log("Downloading NCBI COG definitions...")
        def_text = download_text(COG_DEF_URL, timeout=300)
        def_file.write_text(def_text, encoding="utf-8")

    if fun_file.exists() and fun_file.stat().st_size > 0 and not force:
        fun_text = fun_file.read_text(encoding="utf-8", errors="replace")
        log(f"Using existing COG categories -> {fun_file}")
    else:
        log("Downloading NCBI COG categories...")
        fun_text = download_text(COG_FUN_URL, timeout=300)
        fun_file.write_text(fun_text, encoding="utf-8")

    cog_rows: list[list[str]] = []
    for line in def_text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        parts += [""] * (7 - len(parts))
        if parts[0].startswith("COG"):
            cog_rows.append(parts[:7])

    cat_rows: list[list[str]] = []
    for line in fun_text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        parts += [""] * (4 - len(parts))
        cat_rows.append(parts[:4])

    with (cog_dir / "COG.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tdescription\n")
        for cog_id, _category, name, gene, *_rest in cog_rows:
            desc = f"{gene}; {name}" if gene else name
            f.write(f"{cog_id}\t{clean_cell(desc)}\n")

    with (cog_dir / "COG_long.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tcategory\tname\tgene\tpathway\tpubmed\tpdb\n")
        for row in cog_rows:
            f.write("\t".join(clean_cell(v) for v in row) + "\n")

    with (cog_dir / "COG_category.tsv").open("w", encoding="utf-8", newline="\n") as f:
        f.write("id\tdescription\tgroup\tcolor\n")
        for category_id, group, color, description in cat_rows:
            f.write(f"{category_id}\t{clean_cell(description)}\t{clean_cell(group)}\t{clean_cell(color)}\n")

    log(f"  wrote {len(cog_rows):,} COG rows and {len(cat_rows):,} category rows -> {cog_dir}")


def generate_cazy(outdir: Path, max_family_number: int = 500) -> None:
    cazy_dir = outdir / "cazy"
    cazy_dir.mkdir(parents=True, exist_ok=True)
    prefix_names = {
        "GH": "Glycoside Hydrolase",
        "GT": "GlycosylTransferase",
        "PL": "Polysaccharide Lyase",
        "CE": "Carbohydrate Esterase",
        "AA": "Auxiliary Activity",
        "CBM": "Carbohydrate-Binding Module",
    }
    rows = [
        (f"{prefix}{i}", f"{description} family {i}")
        for prefix, description in prefix_names.items()
        for i in range(1, max_family_number + 1)
    ]
    count = write_simple_tsv(cazy_dir / "CAZy.tsv", rows)
    log(f"Generated {count:,} CAZy rows -> {cazy_dir / 'CAZy.tsv'}")


def create_kegg_tc_placeholder(outdir: Path) -> None:
    tc_dir = outdir / "kegg_tc"
    tc_dir.mkdir(parents=True, exist_ok=True)
    (tc_dir / "KEGG_TC.tsv").write_text("id\tdescription\n", encoding="utf-8")
    (tc_dir / "README.txt").write_text(
        "KEGG_TC values are transport classification IDs. KEGG REST does not provide "
        "a simple list/tc endpoint; add a TCDB-derived table here when needed.\n",
        encoding="utf-8",
    )
    log(f"Created KEGG_TC placeholder -> {tc_dir / 'KEGG_TC.tsv'}")


def copy_manual_eggnog(path: Path, outdir: Path) -> Path:
    eggnog_dir = outdir / "eggnog_ogs"
    eggnog_dir.mkdir(parents=True, exist_ok=True)
    target = eggnog_dir / "e5.og_annotations.tsv"
    source = path.expanduser().resolve()
    if source != target.resolve():
        shutil.copy2(source, target)
    log(f"Registered manual eggNOG OG annotations -> {target}")
    return target


def create_merged_lookup(outdir: Path) -> Path:
    ec_simple = outdir / "ec" / "EC.tsv"
    if not ec_simple.exists():
        ec_simple = outdir / "kegg" / "EC.tsv"

    specs = [
        ("KEGG_ko", outdir / "kegg" / "KEGG_ko.tsv"),
        ("KEGG_Pathway", outdir / "kegg" / "KEGG_Pathway.tsv"),
        ("KEGG_Module", outdir / "kegg" / "KEGG_Module.tsv"),
        ("KEGG_Reaction", outdir / "kegg" / "KEGG_Reaction.tsv"),
        ("KEGG_rclass", outdir / "kegg" / "KEGG_rclass.tsv"),
        ("EC", ec_simple),
        ("BRITE", outdir / "kegg" / "BRITE.tsv"),
        ("GO", outdir / "go" / "GO.tsv"),
        ("PFAM", outdir / "pfam" / "PFAM.tsv"),
        ("BiGG_Reaction", outdir / "bigg" / "BiGG_Reaction.tsv"),
        ("CAZy", outdir / "cazy" / "CAZy.tsv"),
        ("KEGG_TC", outdir / "kegg_tc" / "KEGG_TC.tsv"),
        ("COG", outdir / "cog" / "COG.tsv"),
        ("COG_category", outdir / "cog" / "COG_category.tsv"),
    ]
    merged = outdir / "annotation_lookup_merged.tsv"
    count = 0
    with merged.open("w", encoding="utf-8", newline="\n") as out:
        out.write("source\tid\tdescription\n")
        for source, table in specs:
            if not table.exists():
                continue
            with table.open("r", encoding="utf-8", errors="replace") as f:
                next(f, None)
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    out.write(f"{source}\t{parts[0]}\t{parts[1]}\n")
                    count += 1
    log(f"Created merged lookup with {count:,} rows -> {merged}")
    return merged


def iter_simple_terms(source: str, path: Path) -> Iterator[Term]:
    if not path.exists():
        return
    with path.open("r", encoding="utf-8", errors="replace") as f:
        next(f, None)
        for line in f:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2 and parts[0]:
                yield Term(source=source, term_id=parts[0], description=parts[1])


def build_sqlite(outdir: Path, sqlite_out: Path, eggnog_ogs: Path | None = None) -> None:
    sqlite_out.parent.mkdir(parents=True, exist_ok=True)
    if sqlite_out.exists():
        sqlite_out.unlink()

    log(f"Building SQLite database -> {sqlite_out}")
    conn = sqlite3.connect(sqlite_out)
    try:
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.executescript(
            """
            CREATE TABLE metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );

            CREATE TABLE annotation_term (
                source TEXT NOT NULL,
                term_id TEXT NOT NULL,
                name TEXT NOT NULL DEFAULT '',
                description TEXT NOT NULL DEFAULT '',
                namespace TEXT NOT NULL DEFAULT '',
                category TEXT NOT NULL DEFAULT '',
                PRIMARY KEY (source, term_id)
            );

            CREATE TABLE annotation_property (
                source TEXT NOT NULL,
                term_id TEXT NOT NULL,
                property_key TEXT NOT NULL,
                property_value TEXT NOT NULL,
                PRIMARY KEY (source, term_id, property_key)
            );

            CREATE TABLE annotation_alias (
                source TEXT NOT NULL,
                alias TEXT NOT NULL,
                canonical_id TEXT NOT NULL,
                PRIMARY KEY (source, alias)
            );
            """
        )
        conn.executemany(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            [
                ("schema_version", "1"),
                ("builder", Path(__file__).name),
                ("build_time_utc", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())),
                ("eggnog_ogs_included", "true" if eggnog_ogs else "false"),
            ],
        )

        terms: list[Term] = []
        ec_simple = outdir / "ec" / "EC.tsv"
        if not ec_simple.exists():
            ec_simple = outdir / "kegg" / "EC.tsv"

        simple_specs = [
            ("KEGG_ko", outdir / "kegg" / "KEGG_ko.tsv"),
            ("KEGG_Pathway", outdir / "kegg" / "KEGG_Pathway.tsv"),
            ("KEGG_Module", outdir / "kegg" / "KEGG_Module.tsv"),
            ("KEGG_Reaction", outdir / "kegg" / "KEGG_Reaction.tsv"),
            ("KEGG_rclass", outdir / "kegg" / "KEGG_rclass.tsv"),
            ("BRITE", outdir / "kegg" / "BRITE.tsv"),
            ("PFAM", outdir / "pfam" / "PFAM.tsv"),
            ("BiGG_Reaction", outdir / "bigg" / "BiGG_Reaction.tsv"),
            ("CAZy", outdir / "cazy" / "CAZy.tsv"),
            ("KEGG_TC", outdir / "kegg_tc" / "KEGG_TC.tsv"),
            ("COG_category", outdir / "cog" / "COG_category.tsv"),
        ]
        for source, path in simple_specs:
            terms.extend(iter_simple_terms(source, path) or [])

        go_full = outdir / "go" / "GO_full.tsv"
        if go_full.exists():
            with go_full.open("r", encoding="utf-8", errors="replace") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    terms.append(
                        Term(
                            source="GO",
                            term_id=row.get("id", ""),
                            name=row.get("name", ""),
                            description=row.get("name", ""),
                            namespace=row.get("namespace", ""),
                        )
                    )
        else:
            terms.extend(iter_simple_terms("GO", outdir / "go" / "GO.tsv") or [])

        cog_long = outdir / "cog" / "COG_long.tsv"
        if cog_long.exists():
            with cog_long.open("r", encoding="utf-8", errors="replace") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    terms.append(
                        Term(
                            source="COG",
                            term_id=row.get("id", ""),
                            name=row.get("name", ""),
                            description=row.get("name", ""),
                            category=row.get("category", ""),
                        )
                    )

        ec_full = outdir / "ec" / "EC_full.tsv"
        ec_properties: list[tuple[str, str, str, str]] = []
        if ec_full.exists():
            with ec_full.open("r", encoding="utf-8", errors="replace") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    term_id = row.get("id", "")
                    terms.append(Term(source="EC", term_id=term_id, description=row.get("description", "")))
                    ec_properties.append(("EC", term_id, "reaction", row.get("reaction", "")))
                    ec_properties.append(("EC", term_id, "synonyms", row.get("synonyms", "")))
        else:
            terms.extend(iter_simple_terms("EC", ec_simple) or [])

        conn.executemany(
            """
            INSERT OR REPLACE INTO annotation_term(source, term_id, name, description, namespace, category)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            [(t.source, t.term_id, t.name, t.description, t.namespace, t.category) for t in terms if t.term_id],
        )
        conn.executemany(
            """
            INSERT OR REPLACE INTO annotation_property(source, term_id, property_key, property_value)
            VALUES (?, ?, ?, ?)
            """,
            [p for p in ec_properties if p[1]],
        )

        # Useful aliases for KEGG identifiers commonly seen in eggNOG output.
        aliases: list[tuple[str, str, str]] = []
        for source, prefix, table in [
            ("KEGG_ko", "ko:", outdir / "kegg" / "KEGG_ko.tsv"),
            ("KEGG_Reaction", "rn:", outdir / "kegg" / "KEGG_Reaction.tsv"),
            ("KEGG_Module", "md:", outdir / "kegg" / "KEGG_Module.tsv"),
            ("KEGG_rclass", "rc:", outdir / "kegg" / "KEGG_rclass.tsv"),
            ("EC", "ec:", ec_simple),
        ]:
            for term in iter_simple_terms(source, table) or []:
                aliases.append((source, prefix + term.term_id, term.term_id))
        for term in iter_simple_terms("KEGG_Pathway", outdir / "kegg" / "KEGG_Pathway.tsv") or []:
            aliases.append(("KEGG_Pathway", "path:" + term.term_id, term.term_id))
            if term.term_id.startswith("map"):
                aliases.append(("KEGG_Pathway", term.term_id.replace("map", "ko", 1), term.term_id))
        conn.executemany(
            "INSERT OR IGNORE INTO annotation_alias(source, alias, canonical_id) VALUES (?, ?, ?)",
            aliases,
        )

        if eggnog_ogs and eggnog_ogs.exists():
            conn.execute(
                """
                CREATE TABLE eggnog_og_annotation (
                    tax_id TEXT NOT NULL,
                    og_id TEXT NOT NULL,
                    category TEXT NOT NULL DEFAULT '',
                    description TEXT NOT NULL DEFAULT '',
                    id TEXT GENERATED ALWAYS AS (og_id || '@' || tax_id) VIRTUAL,
                    PRIMARY KEY (tax_id, og_id)
                )
                """
            )
            rows = []
            with eggnog_ogs.open("r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    tax_id, og_id, category, description = parts[:4]
                    rows.append((tax_id, og_id, category, description))
                    if len(rows) >= 100_000:
                        conn.executemany(
                            "INSERT OR REPLACE INTO eggnog_og_annotation(tax_id, og_id, category, description) VALUES (?, ?, ?, ?)",
                            rows,
                        )
                        rows.clear()
                if rows:
                    conn.executemany(
                        "INSERT OR REPLACE INTO eggnog_og_annotation(tax_id, og_id, category, description) VALUES (?, ?, ?, ?)",
                        rows,
                    )
            conn.execute("CREATE INDEX idx_eggnog_og_annotation_id ON eggnog_og_annotation(id)")

        conn.executescript(
            """
            CREATE INDEX idx_annotation_term_id ON annotation_term(term_id);
            CREATE INDEX idx_annotation_alias_canonical ON annotation_alias(source, canonical_id);

            CREATE VIEW go_obo AS
            SELECT term_id AS id, name, namespace FROM annotation_term WHERE source = 'GO';

            CREATE VIEW enzyme AS
            SELECT
                t.term_id AS id,
                t.description AS de,
                COALESCE(ca.property_value, '') AS ca,
                COALESCE(an.property_value, '') AS an
            FROM annotation_term t
            LEFT JOIN annotation_property ca
                ON ca.source = t.source AND ca.term_id = t.term_id AND ca.property_key = 'reaction'
            LEFT JOIN annotation_property an
                ON an.source = t.source AND an.term_id = t.term_id AND an.property_key = 'synonyms'
            WHERE t.source = 'EC';

            CREATE VIEW cog_20_def_tab AS
            SELECT term_id AS id, category, description FROM annotation_term WHERE source = 'COG';

            CREATE VIEW kegg_pathway AS
            SELECT
                CASE WHEN term_id LIKE 'map%' THEN substr(term_id, 4) ELSE term_id END AS id,
                term_id AS name,
                description
            FROM annotation_term
            WHERE source = 'KEGG_Pathway';
            """
        )
        if eggnog_ogs and eggnog_ogs.exists():
            conn.execute(
                """
                CREATE VIEW nog_e5_og_annotations AS
                SELECT id, category, description FROM eggnog_og_annotation
                """
            )

        conn.commit()
    finally:
        conn.close()

    log(f"Finished SQLite database -> {sqlite_out}")


def write_readme(outdir: Path) -> None:
    (outdir / "README_annotation_maps.txt").write_text(
        """MetaPilot functional annotation resources

Generated by build_function_resources.py.

The TSV files are the editable/source artifacts. The optional SQLite database
is the runtime lookup index used by report generation.

The full eggNOG OG annotation table is intentionally not downloaded by this
script. Download it manually from eggNOG and pass it with --eggnog-og-annotations
when building func_def.db.

Important notes:
    - EC_full.tsv is parsed from ExPASy enzyme.dat and keeps description,
      reaction, and synonym fields.
    - GO_full.tsv keeps GO namespace from go-basic.obo.
    - KEGG_TC.tsv is a placeholder until a TCDB-derived source is added.
    - annotation_lookup_merged.tsv is a simple source/id/description table for
      lightweight lookups; func_def.db keeps richer fields.
""",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Download/build MetaPilot functional annotation resources.")
    parser.add_argument("--outdir", required=True, help="Output folder for annotation TSV resources.")
    parser.add_argument("--sqlite-out", help="Optional SQLite database path to build, e.g. Resources/function/func_def.db.")
    parser.add_argument("--eggnog-og-annotations", help="Manual path to eggNOG e5.og_annotations.tsv.")
    parser.add_argument("--force", action="store_true", help="Re-download sources even when local files exist.")
    parser.add_argument("--skip-downloads", action="store_true", help="Only rebuild merged TSV/SQLite from existing files.")
    parser.add_argument("--skip-pfam", action="store_true", help="Skip Pfam download/parsing.")
    parser.add_argument("--skip-bigg", action="store_true", help="Skip BiGG download/parsing.")
    parser.add_argument("--cazy-max-family-number", type=int, default=500)
    args = parser.parse_args()

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    failed: list[tuple[str, Exception]] = []
    if not args.skip_downloads:
        steps = [
            ("KEGG", lambda: download_kegg(outdir, force=args.force)),
            ("GO", lambda: download_go(outdir, force=args.force)),
            ("EC", lambda: download_ec(outdir, force=args.force)),
            ("COG", lambda: download_cog(outdir, force=args.force)),
            ("CAZy", lambda: generate_cazy(outdir, max_family_number=args.cazy_max_family_number)),
            ("KEGG_TC", lambda: create_kegg_tc_placeholder(outdir)),
        ]
        if not args.skip_pfam:
            steps.append(("PFAM", lambda: download_pfam(outdir, force=args.force)))
        if not args.skip_bigg:
            steps.append(("BiGG", lambda: download_bigg(outdir, force=args.force)))

        for name, step in steps:
            try:
                step()
            except Exception as e:
                failed.append((name, e))
                log(f"ERROR: {name} failed: {e}")

    eggnog_path: Path | None = None
    if args.eggnog_og_annotations:
        eggnog_path = copy_manual_eggnog(Path(args.eggnog_og_annotations), outdir)
    else:
        existing = outdir / "eggnog_ogs" / "e5.og_annotations.tsv"
        if existing.exists():
            eggnog_path = existing

    create_merged_lookup(outdir)
    write_readme(outdir)

    if args.sqlite_out:
        build_sqlite(outdir, Path(args.sqlite_out).expanduser().resolve(), eggnog_path)

    if failed:
        log("\nCompleted with failures:")
        for name, error in failed:
            log(f"  - {name}: {error}")
        return 1

    log("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
