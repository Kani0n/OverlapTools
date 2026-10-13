#!/usr/bin/env python3

"""
Balrog ORF scores for arbitrary ORFs.

Balrog (Sommer & Salzberg 2021, BalrogCPP v0.5) scores every ORF with two
pretrained temporal convolutional networks and fixed, previously optimized
weights, then chooses genes by a maximum-weight path with overlap penalties.
This module reproduces the per-ORF score (GeneFinder::combine_scores):

    score = (w_gene * gene_prob + w_TIS * TIS_prob + w_start - threshold) * length_nt

- gene_prob: protein gene model, mean logit over the ORF's residues between start
  and stop codon, read 3'->5' (the network is causal in that direction)
- TIS_prob:  translation initiation site model on the 16 bp up- and downstream of
  the start codon
- w_start:   start codon weight; rare start codons get the weight of the nearest
  normal start codon

Positive scores mean the ORF is more gene-like than Balrog's threshold, scaled by
its length. The models are loaded with torch.hub from Markusjsommer/balrog_models
(downloaded once into the database directory).
"""

import os, logging, warnings
import torch
import torch.nn.functional as F

from .orfs import translate, normal_start


logger = logging.getLogger(__name__)


MODEL_REPO = "Markusjsommer/balrog_models"

# previously optimized parameters (BalrogCPP include/GeneFinder.h)
WEIGHT_GENE_PROB = 0.9746869839852076
WEIGHT_TIS_PROB = 0.25380288790532707
START_WEIGHTS = {"ATG": 0.84249804151264, "TTG": 0.7512400826652517, "GTG": 0.7083689705744909}
SCORE_THRESHOLD = 0.47256101519707244
PROB_THRESHOLD = SCORE_THRESHOLD * (WEIGHT_GENE_PROB + WEIGHT_TIS_PROB + max(START_WEIGHTS.values()))

AA_TOKENS = {aa: i for i, aa in enumerate("*LVIMCAGSTPFYWEDNQKRH")}   # "*" and "X" -> 0
NUC_TOKENS = {"A": 0, "T": 1, "G": 2, "C": 3}                          # others -> 0
TIS_FLANK = 16
BATCH_SIZE = 64


def load_models(db_dir):
    torch.hub.set_dir(os.path.join(db_dir, "balrog"))
    with warnings.catch_warnings():
        # the Balrog models use the deprecated torch.nn.utils.weight_norm
        warnings.filterwarnings("ignore", message=".*weight_norm.*", category=FutureWarning)
        gene_model = torch.hub.load(MODEL_REPO, "geneTCN", trust_repo=True, verbose=False)
        tis_model = torch.hub.load(MODEL_REPO, "tisTCN", trust_repo=True, verbose=False)
    return gene_model.eval(), tis_model.eval()


def gene_probabilities(proteins, model):
    """Gene model probability per protein (each given 3'->5', without start and stop codon)."""
    probs = [None] * len(proteins)
    order = sorted(range(len(proteins)), key=lambda i: len(proteins[i]))
    with torch.no_grad():
        for b in range(0, len(order), BATCH_SIZE):
            batch = order[b:b + BATCH_SIZE]
            max_len = max(len(proteins[i]) for i in batch)
            tokens = torch.zeros((len(batch), max_len), dtype=torch.long)   # padded at the 5' end
            for row, i in enumerate(batch):
                tokens[row, :len(proteins[i])] = torch.tensor([AA_TOKENS.get(aa, 0) for aa in proteins[i]])
            logits = model(F.one_hot(tokens, 21).permute(0, 2, 1).float())
            for row, i in enumerate(batch):
                if proteins[i]:
                    probs[i] = torch.sigmoid(logits[row, 0, :len(proteins[i])].mean()).item()
                else:
                    probs[i] = 0.0   # ORF without residues between start and stop codon
    return probs


def tis_probabilities(windows, model):
    """TIS model probability per 32 bp window (given 3'->5', start codon removed)."""
    tokens = torch.tensor([[NUC_TOKENS.get(n, 0) for n in w] for w in windows], dtype=torch.long)
    with torch.no_grad():
        logits = model(F.one_hot(tokens, 4).permute(0, 2, 1).float())
    return torch.sigmoid(logits[:, 0]).tolist()


def tis_window(orf):
    if orf.start < TIS_FLANK or len(orf.seq) - orf.start < TIS_FLANK + 3:
        return "A" * (2 * TIS_FLANK)   # Balrog's placeholder for starts at contig edges
    window = orf.seq[orf.start - TIS_FLANK:orf.start] + orf.seq[orf.start + 3:orf.start + 3 + TIS_FLANK]
    return window[::-1]


def score(orfs, genome, db_dir, **kwargs):
    """Balrog ORF score of every ORF."""
    logger.info("Loading Balrog gene and TIS models...")
    gene_model, tis_model = load_models(db_dir)

    proteins = [translate(orf.seq[orf.start + 3:orf.stop])[::-1] for orf in orfs]
    gene_probs = gene_probabilities(proteins, gene_model)
    tis_probs = tis_probabilities([tis_window(orf) for orf in orfs], tis_model)

    scores = {}
    for orf, gene_prob, tis_prob in zip(orfs, gene_probs, tis_probs):
        start_weight = START_WEIGHTS[normal_start(orf.start_codon)]
        per_base = WEIGHT_GENE_PROB * gene_prob + WEIGHT_TIS_PROB * tis_prob + start_weight - PROB_THRESHOLD
        scores[orf.header] = per_base * (orf.stop - orf.start)
        logger.debug(f"balrog {orf.header}: gene_prob={gene_prob:.3f}, TIS_prob={tis_prob:.3f}, start_weight={start_weight:.3f}")
    return scores
