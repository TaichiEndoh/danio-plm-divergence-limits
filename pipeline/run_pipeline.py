import argparse
import csv
import gzip
import logging
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

GENEID_PATTERN = re.compile(r"GeneID:(\d+)")


def setup_logger(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "pipeline.log"
    logger = logging.getLogger("zebrafish_go_pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers = []
    file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger


def resolve_path(base: Path, path_str: str) -> Path:
    path = Path(path_str)
    if not path.is_absolute():
        path = (base / path_str).resolve()
    return path


def extract_from_zip(zip_path: Path, target_filename: str, dest_dir: Path, logger: logging.Logger) -> Optional[Path]:
    if not zip_path.exists():
        return None
    dest_dir.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(zip_path) as zf:
            matches = [name for name in zf.namelist() if Path(name).name == target_filename]
            if not matches:
                logger.warning("%s does not contain %s", zip_path, target_filename)
                return None
            member = matches[0]
            target_path = dest_dir / member
            target_path.parent.mkdir(parents=True, exist_ok=True)
            zf.extract(member, dest_dir)
            logger.info("Extracted %s from %s to %s", target_filename, zip_path, target_path)
            return target_path
    except zipfile.BadZipFile:
        logger.error("Could not read %s as a zip archive", zip_path)
        return None


def locate_input_file(
    hint: str,
    expected_filename: str,
    base_dir: Path,
    data_raw_dir: Path,
    logger: logging.Logger,
    alternate_patterns: Optional[Sequence[str]] = None,
) -> Optional[Path]:
    candidates: List[Path] = []
    if hint and hint.lower() != "auto":
        path = Path(hint)
        if not path.is_absolute():
            path = (base_dir / hint).resolve()
        candidates.append(path)
    # look for default file placements
    candidates.extend(
        [
            data_raw_dir / expected_filename,
            base_dir / expected_filename,
        ]
    )

    for cand in candidates:
        if cand.is_file():
            if cand.suffix == ".zip":
                extracted = extract_from_zip(cand, expected_filename, data_raw_dir, logger)
                if extracted and extracted.is_file():
                    return extracted
            else:
                logger.info("Using %s for %s", cand, expected_filename)
                return cand

    # check bundled archives (local + default mount)
    archive_candidates = [
        Path("/mnt/data/public-archivedwl-762.zip"),
        base_dir / "public-archivedwl-762.zip",
        data_raw_dir / "public-archivedwl-762.zip",
    ]
    for archive_path in archive_candidates:
        if archive_path.exists():
            extracted = extract_from_zip(archive_path, expected_filename, data_raw_dir / archive_path.stem, logger)
            if extracted and extracted.is_file():
                return extracted

    # fallback: search by basename to allow csv/tsv variants
    stem = Path(expected_filename).stem
    patterns = [stem]
    if alternate_patterns:
        patterns.extend(alternate_patterns)
    for search_dir in [data_raw_dir, base_dir]:
        for pattern in patterns:
            for candidate in sorted(search_dir.glob(f"{pattern}*")):
                if not candidate.is_file():
                    continue
                if candidate.name == expected_filename:
                    return candidate
                logger.info("Using %s as fallback for %s", candidate, expected_filename)
                return candidate

    logger.warning("Could not find %s; expected at %s or provide via --%s", expected_filename, data_raw_dir, expected_filename.split('.')[0])
    return None


def parse_fasta(fasta_path: Optional[Path], output_path: Path, logger: logging.Logger) -> List[Dict[str, str]]:
    entries: List[Dict[str, str]] = []
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not fasta_path or not fasta_path.exists():
        logger.warning("FASTA file not found; creating empty %s", output_path)
        write_tsv(output_path, ["protein_id_full", "protein_id_nover"], [])
        return entries

    with fasta_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    continue
                token = header.split()[0]
                if not token:
                    continue
                nover = token.split(".")[0]
                entries.append({"protein_id_full": token, "protein_id_nover": nover})
    write_tsv(output_path, ["protein_id_full", "protein_id_nover"], entries)
    logger.info("Extracted %d protein IDs from %s", len(entries), fasta_path)
    return entries


def open_textmaybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return path.open("r", encoding="utf-8", errors="ignore")


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Dict[str, str]]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow([row.get(col, "") for col in header])


def map_proteins_to_geneid(
    proteins: List[Dict[str, str]],
    gene2accession_path: Path,
    output_path: Path,
    logger: logging.Logger,
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not gene2accession_path or not gene2accession_path.exists():
        logger.warning("gene2accession file missing at %s; protein to gene mapping skipped", gene2accession_path)
        write_tsv(output_path, ["protein_id_full", "protein_id_nover", "gene_id"], [])
        return defaultdict(set), defaultdict(set), defaultdict(set)

    target_ids = {p["protein_id_nover"] for p in proteins}
    accession_to_gene: Dict[str, Set[str]] = defaultdict(set)
    total_lines = 0
    matched_lines = 0
    symbol_to_gene: Dict[str, Set[str]] = defaultdict(set)
    with open_textmaybe_gzip(gene2accession_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            total_lines += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            tax_id = parts[0]
            if tax_id != "7955":
                continue
            protein_accession = parts[5]
            gene_id = parts[1]
            if not protein_accession or protein_accession == "-" or not gene_id:
                continue
            nover = protein_accession.split(".")[0]
            if len(parts) > 15:
                symbol = parts[15].strip()
                if symbol and symbol != "-":
                    symbol_to_gene[symbol.upper()].add(gene_id)
            if nover in target_ids:
                accession_to_gene[nover].add(gene_id)
                matched_lines += 1
    logger.info("Processed %d Danio rerio rows from gene2accession; matched %d entries to target proteins", total_lines, matched_lines)

    rows = []
    matched_proteins = 0
    protein_to_genes: Dict[str, Set[str]] = defaultdict(set)
    gene_to_proteins: Dict[str, Set[str]] = defaultdict(set)
    for protein in proteins:
        nover = protein["protein_id_nover"]
        genes = sorted(accession_to_gene.get(nover, []))
        gene_str = ";".join(genes)
        if genes:
            matched_proteins += 1
            for gene in genes:
                protein_to_genes[protein["protein_id_full"]].add(gene)
                gene_to_proteins[gene].add(protein["protein_id_full"])
        rows.append({
            "protein_id_full": protein["protein_id_full"],
            "protein_id_nover": nover,
            "gene_id": gene_str,
        })
    write_tsv(output_path, ["protein_id_full", "protein_id_nover", "gene_id"], rows)
    total_proteins = len(proteins)
    ratio = (matched_proteins / total_proteins * 100) if total_proteins else 0
    logger.info(
        "Mapped %d/%d proteins to GeneIDs (%.2f%%)",
        matched_proteins,
        total_proteins,
        ratio,
    )
    return protein_to_genes, gene_to_proteins, symbol_to_gene


def load_zfin_gene_map(mapping_path: Path, logger: logging.Logger) -> Dict[str, str]:
    if not mapping_path.exists():
        logger.warning("Optional ZFIN gene to GeneID mapping file missing at %s", mapping_path)
        return {}
    mapping: Dict[str, str] = {}
    with mapping_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                continue
            zfin_id = row[0].strip()
            gene_id = row[1].strip()
            if zfin_id and gene_id:
                mapping[zfin_id] = gene_id
    logger.info("Loaded %d ZFIN↔GeneID mappings", len(mapping))
    return mapping


def resolve_symbol_to_gene(
    symbol: str,
    synonyms: str,
    symbol_map: Dict[str, Set[str]],
) -> str:
    candidates: List[str] = []
    if symbol:
        candidates.append(symbol)
    if synonyms:
        candidates.extend(synonyms.split("|"))
    for candidate in candidates:
        key = candidate.strip().upper()
        if not key:
            continue
        genes = symbol_map.get(key)
        if not genes:
            continue
        if len(genes) == 1:
            return next(iter(genes))
        return sorted(genes)[0]
    return ""


def detect_gene_id(row: List[str], zfin_id: str, mapping: Dict[str, str]) -> str:
    if zfin_id in mapping:
        return mapping[zfin_id]
    candidate_fields = []
    for idx in (7, 9, 10, 15, 16):
        if idx < len(row):
            candidate_fields.append(row[idx])
    for field in candidate_fields:
        if not field or field == "-":
            continue
        match = GENEID_PATTERN.search(field)
        if match:
            return match.group(1)
    return ""


def parse_gaf(
    gaf_path: Path,
    mapping_path: Path,
    symbol_to_gene: Dict[str, Set[str]],
    logger: logging.Logger,
) -> List[Dict[str, str]]:
    annotations: List[Dict[str, str]] = []
    if not gaf_path.exists():
        logger.warning("ZFIN GAF file not found at %s", gaf_path)
        return annotations
    zfin_map = load_zfin_gene_map(mapping_path, logger)
    with open_textmaybe_gzip(gaf_path) as handle:
        for line in handle:
            if not line or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            zfin_gene_id = parts[1]
            go_id = parts[4]
            evidence = parts[6]
            aspect = parts[8]
            synonyms = parts[10] if len(parts) > 10 else ""
            assigned_by = parts[14] if len(parts) > 14 else ""
            gene_id = detect_gene_id(parts, zfin_gene_id, zfin_map)
            if not gene_id:
                symbol = parts[2] if len(parts) > 2 else ""
                gene_id = resolve_symbol_to_gene(symbol, synonyms, symbol_to_gene)
            annotations.append(
                {
                    "zfin_gene_id": zfin_gene_id,
                    "gene_id": gene_id,
                    "go_id": go_id,
                    "aspect": aspect,
                    "evidence": evidence,
                    "assigned_by": assigned_by,
                }
            )
    logger.info("Parsed %d GO annotations from %s", len(annotations), gaf_path)
    return annotations


def assign_go_to_proteins(
    protein_to_genes: Dict[str, Set[str]],
    gene_to_proteins: Dict[str, Set[str]],
    annotations: List[Dict[str, str]],
    output_path: Path,
    logger: logging.Logger,
) -> List[Dict[str, str]]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not annotations:
        logger.warning("No GO annotations available; writing empty %s", output_path)
        write_tsv(
            output_path,
            ["protein_id_full", "gene_id", "zfin_gene_id", "go_id", "aspect", "evidence", "assigned_by"],
            [],
        )
        return []
    gene_ids = set(gene_to_proteins.keys())
    rows: List[Dict[str, str]] = []
    seen: Set[Tuple[str, str]] = set()
    for annot in annotations:
        candidate_proteins: Set[str] = set()
        gene_id = annot.get("gene_id")
        if gene_id and gene_id in gene_ids:
            candidate_proteins = gene_to_proteins[gene_id]
        if not candidate_proteins:
            continue
        for protein_id in candidate_proteins:
            dedup_key = (protein_id, annot["go_id"])
            if dedup_key in seen:
                continue
            seen.add(dedup_key)
            rows.append(
                {
                    "protein_id_full": protein_id,
                    "gene_id": annot.get("gene_id", ""),
                    "zfin_gene_id": annot.get("zfin_gene_id", ""),
                    "go_id": annot["go_id"],
                    "aspect": annot.get("aspect", ""),
                    "evidence": annot.get("evidence", ""),
                    "assigned_by": annot.get("assigned_by", ""),
                }
            )
    write_tsv(
        output_path,
        ["protein_id_full", "gene_id", "zfin_gene_id", "go_id", "aspect", "evidence", "assigned_by"],
        rows,
    )
    logger.info("Assigned %d unique GO annotations across %d proteins", len(rows), len({r['protein_id_full'] for r in rows}))
    return rows


def parse_obo(obo_path: Path, logger: logging.Logger) -> Dict[str, Dict[str, str]]:
    terms: Dict[str, Dict[str, str]] = {}
    if not obo_path.exists():
        logger.warning("OBO file %s not found", obo_path)
        return terms
    current_id: Optional[str] = None
    current_data: Dict[str, str] = {}
    with obo_path.open("r", encoding="utf-8", errors="ignore") as handle:
        in_term = False
        for line in handle:
            line = line.strip()
            if not line:
                if in_term and current_id:
                    terms[current_id] = current_data
                in_term = False
                current_id = None
                current_data = {}
                continue
            if line == "[Term]":
                in_term = True
                current_id = None
                current_data = {}
                continue
            if not in_term:
                continue
            if line.startswith("id:"):
                current_id = line.split("id:", 1)[1].strip()
            elif line.startswith("name:"):
                current_data["name"] = line.split("name:", 1)[1].strip()
            elif line.startswith("namespace:"):
                current_data["namespace"] = line.split("namespace:", 1)[1].strip()
            elif line.startswith("def:"):
                current_data["def"] = line.split("def:", 1)[1].strip().strip('"')
        if in_term and current_id:
            terms[current_id] = current_data
    logger.info("Parsed %d GO terms from %s", len(terms), obo_path)
    return terms


def annotate_go_rows(
    protein_rows: List[Dict[str, str]],
    go_terms: Dict[str, Dict[str, str]],
    output_path: Path,
):
    output_rows = []
    for row in protein_rows:
        term = go_terms.get(row["go_id"], {})
        output_rows.append(
            {
                **row,
                "go_name": term.get("name", ""),
                "go_namespace": term.get("namespace", ""),
                "go_definition": term.get("def", ""),
            }
        )
    write_tsv(
        output_path,
        [
            "protein_id_full",
            "gene_id",
            "zfin_gene_id",
            "go_id",
            "aspect",
            "evidence",
            "assigned_by",
            "go_name",
            "go_namespace",
            "go_definition",
        ],
        output_rows,
    )


def summarize_by_aspect(rows: List[Dict[str, str]], output_path: Path):
    aspect_to_annotations: Dict[str, int] = Counter()
    aspect_to_proteins: Dict[str, Set[str]] = defaultdict(set)
    for row in rows:
        aspect = row.get("aspect", "")
        aspect_to_annotations[aspect] += 1
        aspect_to_proteins[aspect].add(row.get("protein_id_full", ""))
    summary_rows = []
    for aspect, count in aspect_to_annotations.items():
        summary_rows.append(
            {
                "aspect": aspect,
                "annotation_count": str(count),
                "protein_count": str(len(aspect_to_proteins.get(aspect, set()))),
            }
        )
    write_tsv(output_path, ["aspect", "annotation_count", "protein_count"], summary_rows)


def top_go_terms(rows: List[Dict[str, str]], aspect: str, go_terms: Dict[str, Dict[str, str]], output_path: Path, top_n: int = 20):
    counter: Dict[str, Set[str]] = defaultdict(set)
    for row in rows:
        if row.get("aspect") != aspect:
            continue
        counter[row["go_id"]].add(row.get("protein_id_full", ""))
    counted = sorted(counter.items(), key=lambda item: len(item[1]), reverse=True)[:top_n]
    output_rows = []
    for go_id, proteins in counted:
        term = go_terms.get(go_id, {})
        output_rows.append(
            {
                "go_id": go_id,
                "go_name": term.get("name", ""),
                "protein_count": str(len(proteins)),
            }
        )
    write_tsv(output_path, ["go_id", "go_name", "protein_count"], output_rows)


def load_slim_terms(slim_path: Path, logger: logging.Logger) -> Set[str]:
    if not slim_path.exists():
        logger.info("GO-slim file %s not found; skipping slim summary", slim_path)
        return set()
    terms = parse_obo(slim_path, logger)
    logger.info("Loaded %d slim GO terms", len(terms))
    return set(terms.keys())


def summarize_slim(rows: List[Dict[str, str]], slim_terms: Set[str], output_path: Path):
    if not slim_terms:
        write_tsv(output_path, ["go_id", "protein_count"], [])
        return
    counter: Dict[str, Set[str]] = defaultdict(set)
    for row in rows:
        go_id = row.get("go_id")
        if go_id in slim_terms:
            counter[go_id].add(row.get("protein_id_full", ""))
    output_rows = []
    for go_id, proteins in sorted(counter.items(), key=lambda item: len(item[1]), reverse=True):
        output_rows.append(
            {
                "go_id": go_id,
                "protein_count": str(len(proteins)),
            }
        )
    write_tsv(output_path, ["go_id", "protein_count"], output_rows)


def normalize_header(value: str) -> str:
    return value.strip().lower().replace(" ", "_")


def read_ortholog_pairs(path: Path, logger: logging.Logger) -> List[Tuple[str, str]]:
    if not path.exists():
        logger.warning("Ortholog file %s missing", path)
        return []
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    rows: List[List[str]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        for row in reader:
            cleaned = [cell.strip() for cell in row]
            if any(cleaned):
                rows.append(cleaned)
    if not rows:
        return []

    header = rows[0]
    normalized_header = [normalize_header(col) for col in header]
    header_text = " ".join(normalized_header)
    has_header = any(keyword in header_text for keyword in ["zebrafish", "aescallii", "danio"])
    start_idx = 1 if has_header else 0

    def pick(record: Dict[str, str], candidates: Sequence[str], default: str) -> str:
        for cand in candidates:
            if cand in record and record[cand]:
                return record[cand]
        return default

    pairs: List[Tuple[str, str]] = []
    for row in rows[start_idx:]:
        if len(row) < 2:
            continue
        if has_header:
            record = {normalized_header[i]: row[i] for i in range(min(len(row), len(normalized_header)))}
            zebrafish_id = pick(
                record,
                ["zebrafish_id", "danio_rerio", "zebrafish", "np_id", "source"],
                row[0],
            )
            aescallii_id = pick(
                record,
                ["aescallii_id", "danio_aesculapii", "aescallii", "target", "other_id"],
                row[1],
            )
        else:
            zebrafish_id, aescallii_id = row[0], row[1]
        zebrafish_id = zebrafish_id.strip()
        aescallii_id = aescallii_id.strip()
        if zebrafish_id and aescallii_id:
            pairs.append((zebrafish_id, aescallii_id))
    logger.info("Loaded %d ortholog pairs from %s", len(pairs), path)
    return pairs


def build_protein_lookup(rows: List[Dict[str, str]]) -> Tuple[Dict[str, List[Dict[str, str]]], Dict[str, List[Dict[str, str]]]]:
    by_full: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    by_nover: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        protein_id = row.get("protein_id_full", "")
        by_full[protein_id].append(row)
        nover = protein_id.split(".")[0]
        by_nover[nover].append(row)
    return by_full, by_nover


def transfer_annotations(
    ortholog_path: Optional[Path],
    protein_rows: List[Dict[str, str]],
    output_path: Path,
    logger: logging.Logger,
):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pairs: List[Tuple[str, str]] = []
    if ortholog_path and ortholog_path.exists():
        pairs = read_ortholog_pairs(ortholog_path, logger)
    else:
        logger.warning("Ortholog pairs file missing; cannot transfer annotations")
    if not pairs or not protein_rows:
        write_tsv(
            output_path,
            ["aescallii_id", "zebrafish_id", "go_id", "aspect", "evidence", "note", "assigned_by"],
            [],
        )
        return
    by_full, by_nover = build_protein_lookup(protein_rows)
    rows = []
    annotated_pairs = 0
    for zebrafish_id, aescallii_id in pairs:
        matched_rows = list(by_full.get(zebrafish_id, []))
        if not matched_rows:
            nover = zebrafish_id.split(".")[0]
            matched_rows = list(by_nover.get(nover, []))
        if not matched_rows:
            continue
        annotated_pairs += 1
        for row in matched_rows:
            rows.append(
                {
                    "aescallii_id": aescallii_id,
                    "zebrafish_id": row.get("protein_id_full", zebrafish_id),
                    "go_id": row.get("go_id", ""),
                    "aspect": row.get("aspect", ""),
                    "evidence": row.get("evidence", ""),
                    "assigned_by": row.get("assigned_by", ""),
                    "note": "transferred_from_zebrafish_ortholog",
                }
            )
    write_tsv(
        output_path,
        ["aescallii_id", "zebrafish_id", "go_id", "aspect", "evidence", "assigned_by", "note"],
        rows,
    )
    logger.info(
        "Transferred GO annotations to %d/%d aescallii ortholog entries",
        annotated_pairs,
        len(pairs),
    )


def ensure_dirs(paths: Sequence[Path]):
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def main():
    parser = argparse.ArgumentParser(description="Danio rerio GO annotation pipeline")
    parser.add_argument("--fasta", default="auto", help="Path to zebrafish RefSeq protein FASTA or 'auto'")
    parser.add_argument("--ortholog", default="auto", help="Path to ortholog pairs TSV or 'auto'")
    parser.add_argument("--outdir", default="results", help="Output directory for result tables")
    parser.add_argument("--data-raw-dir", default="data_raw", help="Directory for raw data")
    parser.add_argument("--data-go-dir", default="data_go", help="Directory containing GO resources")
    parser.add_argument("--logs-dir", default="logs", help="Directory for logs")
    parser.add_argument("--zfin-gene-map", default="data_go/zfin_geneid_mappings.tsv", help="Optional ZFIN↔GeneID mapping TSV")
    parser.add_argument("--goslim", default="data_go/goslim_generic.obo", help="Optional GO-slim OBO file")
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[1]
    data_raw_dir = resolve_path(project_root, args.data_raw_dir)
    data_go_dir = resolve_path(project_root, args.data_go_dir)
    out_dir = resolve_path(project_root, args.outdir)
    logs_dir = resolve_path(project_root, args.logs_dir)
    ensure_dirs([data_raw_dir, data_go_dir, out_dir, logs_dir])

    logger = setup_logger(logs_dir)
    logger.info("Starting zebrafish→aescallii GO annotation pipeline")

    fasta_path = locate_input_file(args.fasta, "zebrafish_np_proteins.fasta", project_root, data_raw_dir, logger)
    ortholog_path = locate_input_file(
        args.ortholog,
        "ortholog_pairs.tsv",
        project_root,
        data_raw_dir,
        logger,
        alternate_patterns=["ortholog", "ortholog_table"],
    )

    protein_list_path = out_dir / "protein_list.tsv"
    proteins = parse_fasta(fasta_path, protein_list_path, logger)

    gene2accession_path = data_go_dir / "gene2accession.gz"
    protein_to_genes, gene_to_proteins, symbol_to_gene = map_proteins_to_geneid(
        proteins,
        gene2accession_path,
        out_dir / "protein_to_geneid.tsv",
        logger,
    )

    gaf_path = data_go_dir / "gene_association.zfin.gz"
    zfin_gene_map_path = resolve_path(project_root, args.zfin_gene_map) if not Path(args.zfin_gene_map).is_absolute() else Path(args.zfin_gene_map)
    annotations = parse_gaf(gaf_path, zfin_gene_map_path, symbol_to_gene, logger)

    protein_go_rows = assign_go_to_proteins(
        protein_to_genes,
        gene_to_proteins,
        annotations,
        out_dir / "zebrafish_protein_go.tsv",
        logger,
    )

    go_terms = parse_obo(data_go_dir / "go-basic.obo", logger)
    annotate_go_rows(
        protein_go_rows,
        go_terms,
        out_dir / "zebrafish_protein_go_annotated.tsv",
    )

    summarize_by_aspect(protein_go_rows, out_dir / "go_counts_by_aspect.tsv")
    for aspect, suffix in [("P", "BP"), ("F", "MF"), ("C", "CC")]:
        top_go_terms(
            protein_go_rows,
            aspect,
            go_terms,
            out_dir / f"top_go_terms_{suffix}.tsv",
        )

    slim_terms = load_slim_terms(resolve_path(project_root, args.goslim), logger)
    summarize_slim(protein_go_rows, slim_terms, out_dir / "go_slim_counts.tsv")

    transfer_annotations(
        ortholog_path,
        protein_go_rows,
        out_dir / "aescallii_inferred_go.tsv",
        logger,
    )

    logger.info("Pipeline completed. Outputs are in %s", out_dir)


if __name__ == "__main__":
    main()
