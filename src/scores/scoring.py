#!/usr/bin/env python3

"""Runs the selected scoring methods and writes the score table."""

import os, shutil, logging, importlib

from .orfs import Genome, load_orfs


logger = logging.getLogger(__name__)


# scoring method -> module in this package; every module provides
# score(orfs, genome, db_dir=..., n_threads=..., tmp_dir=...) -> {ORF header: score}
METHODS = ['prodigal', 'balrog', 'pgap', 'bakta']
NEEDS_GENOME = {'prodigal', 'balrog'}


def score_ORFs(orfs_fasta, methods, out_dir, genome_fasta=None, db_dir=None, n_threads=1):

    scores_path = os.path.join(out_dir, 'scores')
    os.makedirs(scores_path, exist_ok=True)
    scores_file = os.path.join(scores_path, 'ORF_scores.tsv')

    methods = list(dict.fromkeys(methods))
    if NEEDS_GENOME.intersection(methods) and genome_fasta is None:
        raise ValueError(f"--genome is required for {', '.join(sorted(NEEDS_GENOME.intersection(methods)))}")

    genome = Genome(genome_fasta) if genome_fasta else None
    orfs = load_orfs(orfs_fasta, genome)
    logger.info(f"Read {len(orfs)} ORFs" + (f" and {len(genome.contigs)} contigs" if genome else ""))

    tmp_dir = os.path.join(scores_path, 'tmp')
    results = {}
    for method in methods:
        logger.info(f"Scoring ORFs with {method}...")
        module = importlib.import_module(f'.{method}', __package__)
        results[method] = module.score(orfs, genome, db_dir=db_dir, n_threads=n_threads, tmp_dir=tmp_dir)

    if not logger.isEnabledFor(logging.DEBUG):
        shutil.rmtree(tmp_dir, ignore_errors=True)

    with open(scores_file, 'w') as out:
        out.write('\t'.join(['ORF'] + methods) + '\n')
        for orf in orfs:
            values = [results[m].get(orf.header) for m in methods]
            out.write('\t'.join([orf.header] + ['NA' if v is None else f'{v:.3f}' for v in values]) + '\n')

    logger.info(f"Wrote {', '.join(methods)} scores of {len(orfs)} ORFs to {scores_file}")
    return scores_file
