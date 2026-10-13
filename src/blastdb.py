#!/usr/bin/env python3

"""
Build the BLAST nucleotide database with taxids used by find_homologs and align_codons.

Genome FASTA files (e.g. the GTDB representative genomes) are named by their assembly
accession (GCA_/GCF_), while their sequences have their own accessions
(e.g. NZ_JBHTNH010000001.1). The taxonomy dump's taxid.map assigns taxids to assemblies,
makeblastdb needs them per sequence: every sequence gets the taxid of its assembly.
"""

import os, re, gzip, logging, subprocess
from pathlib import Path


logger = logging.getLogger(__name__)


ASSEMBLY_PATTERN = re.compile(r"(GC[AF]_\d+\.\d+)")


def read_assembly_taxids(taxid_map):
    """Assembly accession -> taxid from the taxonomy dump's taxid.map."""
    taxids = {}
    with open(taxid_map) as f:
        for line_number, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                logger.warning(f"Ignoring malformed line {line_number} of {taxid_map}: {line}")
                continue
            taxids[fields[0]] = fields[1]
    return taxids


def assembly_accession(filename):
    """GCA_/GCF_ accession in a file name, None if there is none."""
    match = ASSEMBLY_PATTERN.search(filename)
    return match.group(1) if match else None


def append_genome(fasta_file, combined, taxid, sequence_taxids):
    """
    Append one genome FASTA to the combined FASTA and record sequence -> taxid.
    Sequences whose accession was already added are skipped (makeblastdb -parse_seqids
    does not allow duplicate IDs). Returns the number of added and skipped sequences.
    """
    opener = gzip.open if fasta_file.name.endswith(".gz") else open
    added = skipped = 0
    keep = False
    with opener(fasta_file, "rt") as f:
        for line in f:
            if line.startswith(">"):
                sequence_id = line[1:].split()[0]
                keep = sequence_id not in sequence_taxids
                if keep:
                    sequence_taxids[sequence_id] = taxid
                    added += 1
                else:
                    skipped += 1
                    previous = sequence_taxids[sequence_id]
                    detail = f"taxid {taxid} vs {previous}" if previous != taxid else f"taxid {taxid}"
                    logger.warning(f"Skipping duplicate sequence {sequence_id} in {fasta_file.name} ({detail})")
            if keep:
                combined.write(line)
    return added, skipped


def combine_genomes(genome_dir, assembly_taxids, combined_fasta, sequence_map):
    """Combine all genome FASTA files (.fna, .fna.gz) and write the sequence -> taxid map."""
    fasta_files = sorted(list(Path(genome_dir).rglob("*.fna.gz")) + list(Path(genome_dir).rglob("*.fna")))
    if not fasta_files:
        raise RuntimeError(f"No .fna or .fna.gz files found in {genome_dir}")
    logger.info(f"Found {len(fasta_files):,} genome FASTA files")

    sequence_taxids = {}
    no_accession, no_taxid = [], []
    skipped_total = 0
    with open(combined_fasta, "w") as combined:
        for i, fasta_file in enumerate(fasta_files, 1):
            assembly = assembly_accession(fasta_file.name)
            if assembly is None:
                no_accession.append(fasta_file.name)
                logger.warning(f"No assembly accession in file name, skipping: {fasta_file.name}")
                continue
            if assembly not in assembly_taxids:
                no_taxid.append(assembly)
                logger.warning(f"No taxid for {assembly} in the taxid map, skipping: {fasta_file.name}")
                continue
            added, skipped = append_genome(fasta_file, combined, assembly_taxids[assembly], sequence_taxids)
            skipped_total += skipped
            logger.debug(f"[{i:,}/{len(fasta_files):,}] {assembly}: {added:,} sequences -> taxid {assembly_taxids[assembly]}")
            if i % 10000 == 0:
                logger.info(f"Combined {i:,} / {len(fasta_files):,} genomes")

    with open(sequence_map, "w") as out:
        for sequence_id, taxid in sequence_taxids.items():
            out.write(f"{sequence_id}\t{taxid}\n")

    logger.info(f"Combined {len(sequence_taxids):,} sequences of {len(fasta_files) - len(no_accession) - len(no_taxid):,} genomes "
                f"({len(no_accession):,} files without accession, {len(no_taxid):,} assemblies without taxid, "
                f"{skipped_total:,} duplicate sequences skipped)")


def build_blastdb(genome_dir, taxdump_dir, blastdb_path, keep_intermediate=False):
    """
    Build a BLAST v5 nucleotide database with taxids at `blastdb_path` (path prefix)
    from the genome FASTA files in `genome_dir`, using `taxdump_dir`/taxid.map.
    """
    taxid_map = os.path.join(taxdump_dir, "taxid.map")
    assembly_taxids = read_assembly_taxids(taxid_map)
    logger.info(f"Loaded {len(assembly_taxids):,} assembly taxids from {taxid_map}")

    os.makedirs(os.path.dirname(os.path.abspath(blastdb_path)), exist_ok=True)
    combined_fasta = blastdb_path + ".combined.fna"
    sequence_map = blastdb_path + ".taxid.map"

    combine_genomes(genome_dir, assembly_taxids, combined_fasta, sequence_map)

    logger.info(f"Building BLAST database {blastdb_path} with makeblastdb...")
    subprocess.run(
        ["makeblastdb", "-in", combined_fasta, "-dbtype", "nucl", "-parse_seqids",
         "-taxid_map", sequence_map, "-blastdb_version", "5",
         "-title", os.path.basename(blastdb_path), "-out", blastdb_path],
        stdout=subprocess.DEVNULL, check=True
    )

    result = subprocess.run(["blastdbcmd", "-db", blastdb_path, "-info"], capture_output=True, text=True, check=True)
    summary = next((line.strip() for line in result.stdout.splitlines() if "sequences;" in line), "")
    logger.info(f"BLAST database {blastdb_path}: {summary}")

    if not keep_intermediate:
        for path in (combined_fasta, sequence_map):
            os.remove(path)
    return blastdb_path
