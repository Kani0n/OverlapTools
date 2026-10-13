#!/usr/bin/env python3

"""Shared ORF/genome handling for all scoring methods."""

import re, logging
from dataclasses import dataclass


logger = logging.getLogger(__name__)


STOP_CODONS = {"TAA", "TAG", "TGA"}
NORMAL_STARTS = ("ATG", "GTG", "TTG")

GENETIC_CODE = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        ((a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}

# ORF headers written by classify_ORFs: <contig>_<start>-<end>_<strand>_<length>_...
# <start>/<end> are the 1-based first bases (in strand orientation) of the start and the
# stop codon, ordered low-high.
HEADER_PATTERN = re.compile(r"^(?P<contig>.+?)_(?P<low>\d+)-(?P<high>\d+)_(?P<strand>[+-])(?:_|$)")


def read_fasta(path):
    records = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line)
    return {name: "".join(seq) for name, seq in records.items()}


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def translate(seq):
    return "".join(GENETIC_CODE.get(seq[i:i+3], "X") for i in range(0, len(seq) - 2, 3))


def normal_start(codon):
    """
    Map a start codon to the normal start codon (ATG/GTG/TTG) whose weight it gets:
    the one with the fewest mismatches, ties broken in the order ATG, GTG, TTG
    (ATA/ATC/ATT -> ATG, CTG -> GTG).
    """
    return min(NORMAL_STARTS, key=lambda s: sum(a != b for a, b in zip(s, codon)))


class Genome:
    """Genome contigs, with both strands available 5'->3'."""

    def __init__(self, fasta):
        self.contigs = {name: seq.upper() for name, seq in read_fasta(fasta).items()}
        self._strands = {}

    def strand(self, contig, strand):
        key = (contig, strand)
        if key not in self._strands:
            seq = self.contigs[contig]
            self._strands[key] = seq if strand == "+" else reverse_complement(seq)
        return self._strands[key]


@dataclass
class ORF:
    header: str
    protein: str                # protein from the ORFs FASTA (without stop)
    contig: str = None
    strand: str = None
    start: int = None           # 0-based first base of the start codon on the strand sequence
    stop: int = None            # 0-based first base of the stop codon on the strand sequence
    seq: str = None             # strand sequence of the contig (5'->3')

    @property
    def start_codon(self):
        return self.seq[self.start:self.start+3] if self.seq else None

    @property
    def initiator_protein(self):
        """Protein as translated in vivo: alternative start codons are read as Met."""
        return "M" + self.protein[1:]


def load_orfs(orfs_fasta, genome=None):
    """Read the ORF proteins and, if a genome is given, locate each ORF on it."""
    orfs = []
    for header, protein in read_fasta(orfs_fasta).items():
        orf = ORF(header=header, protein=protein.rstrip("*"))
        if genome is not None:
            match = HEADER_PATTERN.match(header)
            if not match:
                raise ValueError(f"Cannot parse ORF coordinates from header: {header}")
            contig = match["contig"]
            if contig not in genome.contigs:
                raise ValueError(f"Contig {contig} of ORF {header} not found in genome")
            low, high, strand = int(match["low"]), int(match["high"]), match["strand"]
            n = len(genome.contigs[contig])
            orf.contig, orf.strand, orf.seq = contig, strand, genome.strand(contig, strand)
            orf.start, orf.stop = (low - 1, high - 1) if strand == "+" else (n - high, n - low)

            if orf.seq[orf.stop:orf.stop+3] not in STOP_CODONS:
                logger.warning(f"No stop codon at the end of {header} in the genome")
            if translate(orf.seq[orf.start:orf.stop])[1:] != orf.protein[1:]:
                logger.warning(f"Genome translation does not match the ORF sequence of {header}")
        orfs.append(orf)
    return orfs
