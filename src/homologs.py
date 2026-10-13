#!/usr/bin/env python3

import os
import argparse
import csv
import logging
import subprocess
import tempfile
from pathlib import Path


logger = logging.getLogger(__name__)


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "staxids",
]


RANKS = [
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

def run_taxonkit(taxids, taxdump_dir):
    """
    Resolve taxids and return:

        taxid -> {
            "kingdom": (taxid, name),
            "phylum":  (taxid, name),
            ...
        }

    TaxonKit is supplied with an actual temporary file.
    """

    taxids = sorted(set(taxids))

    if not taxids:
        return {}

    with tempfile.TemporaryDirectory() as tmpdir:

        tmpdir = Path(tmpdir)

        taxid_file = tmpdir / "taxids.txt"
        lineage_file = tmpdir / "lineage.tsv"
        reformatted_file = tmpdir / "reformatted.tsv"

        # ------------------------------------------------------
        # Write taxids
        # ------------------------------------------------------

        with taxid_file.open("w") as f:
            for taxid in taxids:
                f.write(f"{taxid}\n")

        # ------------------------------------------------------
        # taxonkit lineage
        # ------------------------------------------------------

        with lineage_file.open("w") as out:
            subprocess.run(
                ["taxonkit", "lineage", "-t", "--data-dir", taxdump_dir, str(taxid_file)],
                stdout=out,
                check=True
            )

        # ------------------------------------------------------
        # taxonkit reformat
        # ------------------------------------------------------

        with reformatted_file.open("w") as out:
            subprocess.run(
                ["taxonkit", "reformat", "--data-dir", taxdump_dir, str(lineage_file), "-f", "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}"],
                stdout=out,
                check=True
            )

        # ------------------------------------------------------
        # Parse results
        # ------------------------------------------------------

        taxonomy = {}

        with reformatted_file.open() as f:

            for line in f:

                line = line.rstrip("\n")

                if not line:
                    continue

                fields = line.split("\t")

                if len(fields) < 10:
                    continue

                try:
                    original_taxid = int(fields[0])
                except ValueError:
                    continue

                # taxonkit lineage -t: taxid, lineage names, lineage taxids
                lineage_names = fields[1].split(";")
                lineage_taxids = fields[2].split(";")

                rank_names = {
                    "kingdom": fields[3],
                    "phylum": fields[4],
                    "class": fields[5],
                    "order": fields[6],
                    "family": fields[7],
                    "genus": fields[8],
                    "species": fields[9],
                }

                rank_taxids = {}

                # Map lineage name -> lineage taxid.
                name_to_taxid = {}

                for tid, name in zip(
                    lineage_taxids,
                    lineage_names,
                ):
                    if name and name != "-":
                        try:
                            name_to_taxid[name] = int(tid)
                        except ValueError:
                            pass

                for rank, name in rank_names.items():

                    if not name or name == "-":
                        continue

                    rank_taxid = name_to_taxid.get(name)

                    if rank_taxid is not None:
                        rank_taxids[rank] = (
                            rank_taxid,
                            name,
                        )

                taxonomy[original_taxid] = rank_taxids

        return taxonomy


def parse_taxids(value):
    """
    Parse BLAST staxids.

    Examples:

        562
        562;570
        N/A
        -
    """

    if not value or value in {"N/A", "-"}:
        return []

    result = []

    for x in value.split(";"):

        x = x.strip()

        if not x:
            continue

        try:
            result.append(int(x))
        except ValueError:
            pass

    return result


def read_blast(path):

    hits = []

    with open(path) as f:

        reader = csv.reader(
            f,
            delimiter="\t",
        )

        for line_number, row in enumerate(
            reader,
            start=1,
        ):

            if len(row) < len(BLAST_COLUMNS):
                raise ValueError(
                    f"Line {line_number}: expected at least {len(BLAST_COLUMNS)} columns, "
                    f"got {len(row)}"
                )

            hit = dict(
                zip(
                    BLAST_COLUMNS,
                    row,
                )
            )

            hit["evalue_float"] = float(
                hit["evalue"]
            )

            hit["bitscore_float"] = float(
                hit["bitscore"]
            )

            hit["taxids"] = parse_taxids(
                hit["staxids"]
            )

            hits.append(hit)

    return hits


def get_top_hits(hits):
    top_hits = {}
    for hit in hits:
        qseqid = hit["qseqid"]
        if qseqid not in top_hits:
            top_hits[qseqid] = hit
            continue
        current = top_hits[qseqid]
        if (hit["evalue_float"], -hit["bitscore_float"]) < (current["evalue_float"], -current["bitscore_float"]):
            top_hits[qseqid] = hit
    return top_hits


def taxonomic_context(subject_taxonomy, reference_taxonomy):
    """
    Return the most specific rank shared by subject and reference
    (e.g. "species" if same species, "genus" if same genus but different
    species, ...), "none" if no rank is shared, or "unclassified" if the
    subject could not be resolved.
    """

    if not subject_taxonomy:
        return "unclassified"

    for rank in reversed(RANKS):
        subject_rank = subject_taxonomy.get(rank)
        reference_rank = reference_taxonomy.get(rank)
        if subject_rank is not None and reference_rank is not None and subject_rank[0] == reference_rank[0]:
            return rank

    return "none"


def infer_reference_taxid(hits):
    """Taxid of the overall top hit (lowest E-value, highest bitscore)."""
    top_hits = get_top_hits(hits)
    best_hit = min(top_hits.values(), key=lambda hit: (hit["evalue_float"], -hit["bitscore_float"]))
    if not best_hit["taxids"]:
        raise RuntimeError("Top BLAST hit has no staxid; cannot infer reference taxid.")
    logger.info(f"Inferred taxid from top hit ({best_hit['qseqid']} vs {best_hit['sseqid']}): {best_hit['taxids'][0]}")
    return best_hit["taxids"][0]


def reference_taxon(reference_taxid, level, taxdump_dir):
    """Taxid and name of the reference's taxon at `level`."""
    taxonomy = run_taxonkit([reference_taxid], taxdump_dir)
    if reference_taxid not in taxonomy:
        raise RuntimeError(f"TaxonKit could not resolve taxid {reference_taxid}")
    rank = taxonomy[reference_taxid].get(level)
    if rank is None:
        raise RuntimeError(f"TaxonKit could not determine {level} for taxid {reference_taxid}")
    return rank


def run_captured(command):
    """Run a command and return its stdout; failures are raised with the command's stderr."""
    # empty stdin: taxonkit reads its input from stdin if that is a pipe
    result = subprocess.run(command, stdin=subprocess.DEVNULL, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed (exit {result.returncode}): {result.stderr.strip()}")
    return result.stdout


def blastdb_taxids(BLASTdb_path):
    """Taxids of the sequences in the BLAST database."""
    stdout = run_captured(["blastdbcmd", "-db", BLASTdb_path, "-tax_info"])
    return {int(line.split("\t")[0]) for line in stdout.splitlines()
            if line and not line.startswith("#")}


def write_search_taxids(taxon_taxid, BLASTdb_path, taxidlist_path, taxdump_dir):
    """
    Write the taxids of all database sequences that belong to the taxon (the taxon
    and all its descendants in the taxdump) and return how many there are.
    """
    stdout = run_captured(["taxonkit", "list", "--data-dir", taxdump_dir, "--ids", str(taxon_taxid), "--indent", ""])
    members = {int(line) for line in stdout.split() if line}
    taxids = sorted(members & blastdb_taxids(BLASTdb_path))
    with open(taxidlist_path, "w") as f:
        f.writelines(f"{taxid}\n" for taxid in taxids)
    return len(taxids)


def annotate_tax(blast_tsv, output_path, taxdump_dir, reference_taxid=None):
    """
    Annotate each BLAST hit with its taxonomic context relative to the reference
    taxid (extra last column): the most specific shared rank, "none" or
    "unclassified". If `reference_taxid` is None, it is inferred from the overall
    top hit.
    """

    hits = read_blast(blast_tsv)
    if not hits:
        logger.warning(f"No BLAST hits found in {blast_tsv}.")
        return None

    if reference_taxid is None:
        reference_taxid = infer_reference_taxid(hits)

    all_taxids = {reference_taxid}
    for hit in hits:
        all_taxids.update(hit["taxids"])

    logger.debug(f"Resolving {len(all_taxids):,} taxids with TaxonKit...")
    taxonomy = run_taxonkit(all_taxids, taxdump_dir)

    if reference_taxid not in taxonomy:
        raise RuntimeError(f"TaxonKit could not resolve taxid {reference_taxid}")
    reference_taxonomy = taxonomy[reference_taxid]

    with open(output_path, "w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for hit in hits:
            # With multiple staxids, report the closest context.
            contexts = [taxonomic_context(taxonomy.get(taxid, {}), reference_taxonomy) for taxid in hit["taxids"]] or ["unclassified"]
            context = max(contexts, key=lambda c: RANKS.index(c) if c in RANKS else -1)
            writer.writerow([hit[column] for column in BLAST_COLUMNS] + [context])

    logger.info(f"Annotated {len(hits):,} hits with their taxonomic context relative to taxid {reference_taxid}.")
    return output_path


# tblastn only benefits from splitting threads by query (-mt_mode 1) for small databases
MT_MODE_QUERY_MAX_BASES = 280_000_000


def subset_size(BLASTdb_path, taxidlist):
    """Number of sequences and bases of the database restricted to the taxid list."""
    stdout = run_captured(["blastdbcmd", "-db", BLASTdb_path, "-taxidlist", taxidlist,
                           "-no_taxid_expansion", "-outfmt", "%l"])
    lengths = [int(x) for x in stdout.split()]
    return len(lengths), sum(lengths)


def run_tBLASTn(aa_fasta, BLASTdb_path, n_threads, blast_tsv, taxidlist=None, max_target_seqs=None):
    command = ["tblastn", "-db", BLASTdb_path, "-query", aa_fasta, "-evalue", "1e-5",
               "-outfmt", "6 " + " ".join(BLAST_COLUMNS), "-out", blast_tsv]
    if taxidlist is not None:
        # the list already contains all descendants (GTDB taxdump), so BLAST must not expand
        # it with its own (NCBI) taxonomy
        command += ["-taxidlist", taxidlist, "-no_taxid_expansion"]
        n_seqs, n_bases = subset_size(BLASTdb_path, taxidlist)
        if n_bases < MT_MODE_QUERY_MAX_BASES:
            # small subset: split threads by query; splitting it between threads is slow
            command += ["-mt_mode", "1"]
            mode = "threads split by query"
        else:
            # large subset: split threads over the subject sequences; more threads than
            # subject sequences crash tblastn
            n_threads = min(n_threads, n_seqs)
            mode = "threads split over the database"
        logger.info(f"Searching {n_seqs:,} sequences ({n_bases / 1e6:,.1f} Mbp) with {n_threads} threads, {mode}")
    command += ["-num_threads", str(n_threads)]
    if max_target_seqs is not None:
        command += ["-max_target_seqs", str(max_target_seqs), "-max_hsps", "1"]
    logger.debug(f"Running: {' '.join(command)}")
    subprocess.run(command, check=True)
    return blast_tsv


def find_homologs(orfs_aa_fasta, BLASTdb_path, taxdump_dir, n_threads, out_dir, level=None, tax_id=None):
    """
    tBLASTn search of the ORF proteins. With `level`, the search is restricted to the
    database sequences of the reference's taxon at that level; the reference is
    `tax_id` or, if not given, the taxid of the top hit of a preliminary search of the
    whole database. All hits are annotated with their taxonomic context and written to
    homologs/tBLASTn.tsv; intermediate files are removed.
    """
    if level is not None and level not in RANKS:
        raise ValueError(f"Invalid taxonomic level '{level}'; choose from {', '.join(RANKS)}")

    homologs_path = os.path.join(out_dir, 'homologs')
    os.makedirs(homologs_path, exist_ok=True)
    raw_path = os.path.join(homologs_path, 'tBLASTn_raw.tsv')
    reference_path = os.path.join(homologs_path, 'tBLASTn_reference.tsv')
    taxidlist_path = os.path.join(homologs_path, 'search_taxids.txt')
    result_path = os.path.join(homologs_path, 'tBLASTn.tsv')

    def finish(reference_taxid):
        if annotate_tax(raw_path, result_path, taxdump_dir, reference_taxid) is None:
            open(result_path, 'w').close()   # no hits
        for path in (raw_path, reference_path, taxidlist_path):
            if os.path.exists(path):
                os.remove(path)
        return result_path

    if level is None:
        logger.info("No taxonomic level given; no filtering is performed, the whole database is searched "
                    "and hits are only annotated with their taxonomic context.")
        run_tBLASTn(orfs_aa_fasta, BLASTdb_path, n_threads, raw_path)
        return finish(tax_id)

    if tax_id is not None:
        reference_taxid = int(tax_id)
        logger.info(f"Using supplied taxid: {reference_taxid}")
    else:
        logger.info("No taxid given; running a preliminary search of the whole database to find the top hit "
                    "(give --taxid to skip this)...")
        run_tBLASTn(orfs_aa_fasta, BLASTdb_path, n_threads, reference_path, max_target_seqs=5)
        hits = read_blast(reference_path)
        os.remove(reference_path)
        if not hits:
            logger.warning("No BLAST hits found in the preliminary search.")
            open(result_path, 'w').close()
            return result_path
        reference_taxid = infer_reference_taxid(hits)

    taxon_taxid, taxon_name = reference_taxon(reference_taxid, level, taxdump_dir)
    n_taxids = write_search_taxids(taxon_taxid, BLASTdb_path, taxidlist_path, taxdump_dir)
    if n_taxids == 0:
        os.remove(taxidlist_path)
        raise RuntimeError(f"The BLAST database contains no sequences of {level} {taxon_name} (taxid {taxon_taxid})")
    logger.info(f"Restricting the search to {level} {taxon_name} (taxid {taxon_taxid}): {n_taxids:,} taxids in the database")

    run_tBLASTn(orfs_aa_fasta, BLASTdb_path, n_threads, raw_path, taxidlist=taxidlist_path)
    return finish(reference_taxid)
