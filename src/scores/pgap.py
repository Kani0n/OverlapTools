#!/usr/bin/env python3

"""
PGAP protein family evidence for arbitrary ORFs.

PGAP (NCBI Prokaryotic Genome Annotation Pipeline) supports and names proteins
with its HMM library (NCBIfams, TIGRFAMs, Pfam; hmm_PGAP.LIB). Every ORF protein
is searched against that library with HMMER (pyhmmer hmmsearch); the score is
the best full-sequence bitscore of a significant hit (E-value <= 0.01, HMMER's
inclusion threshold), 0 if the ORF has no significant hit. Rare start codons are
translated as Met. The library is downloaded once into the database directory.
"""

import os, shutil, logging, urllib.request
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMFile


logger = logging.getLogger(__name__)


LIBRARY_URL = "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB"
INCLUSION_EVALUE = 0.01


def library_path(db_dir):
    pgap_dir = os.path.join(db_dir, "pgap")
    path = os.path.join(pgap_dir, "hmm_PGAP.LIB")
    if not os.path.isfile(path):
        os.makedirs(pgap_dir, exist_ok=True)
        logger.info(f"Downloading PGAP HMM library (~3 GB) to {path}...")
        with urllib.request.urlopen(LIBRARY_URL) as response, open(path + ".part", "wb") as out:
            shutil.copyfileobj(response, out)
        os.rename(path + ".part", path)
    return path


def score(orfs, genome, db_dir, n_threads=1, **kwargs):
    """Best significant PGAP HMM bitscore of every ORF protein."""
    alphabet = Alphabet.amino()
    sequences = [
        TextSequence(name=str(i), sequence=orf.initiator_protein).digitize(alphabet)
        for i, orf in enumerate(orfs)
    ]

    best = [0.0] * len(orfs)
    path = library_path(db_dir)
    logger.info(f"Searching {len(orfs)} ORF proteins against the PGAP HMM library...")
    with HMMFile(path) as hmms:
        for hits in pyhmmer.hmmsearch(hmms, sequences, cpus=n_threads, incE=INCLUSION_EVALUE):
            for hit in hits:
                if hit.included:
                    i = int(hit.name)
                    if hit.score > best[i]:
                        best[i] = hit.score
                        logger.debug(f"pgap {orfs[i].header}: {hits.query.name} bitscore={hit.score:.1f} E={hit.evalue:.2g}")

    return {orf.header: s for orf, s in zip(orfs, best)}
