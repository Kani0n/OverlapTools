#!/usr/bin/env python3

"""
Prodigal-like "gene-likeness" scores for arbitrary ORFs.

Prodigal (and pyrodigal) only scores start nodes it creates itself, i.e.
ATG/GTG/TTG starts of ORFs above a minimum length. Many ORFs of interest
(alternative start codons, short ORFs, ORFs nested in other genes) are
therefore never scored. This module trains pyrodigal on the genome and then
re-implements Prodigal's node scoring (node.c: raw_coding_score, rbs_score,
score_upstream_composition, score_nodes) in Python, so that every ORF can be
scored with the trained parameters, exactly as Prodigal would score that
start/stop pair. The score is Prodigal's gene score (coding + start score).

Ported from pyrodigal (Martin Larralde, https://github.com/althonos/pyrodigal) and
Prodigal (Doug Hyatt, https://github.com/hyattpd/Prodigal), both GPL-3.0.
"""

import math, logging
import pyrodigal

from .orfs import STOP_CODONS, normal_start


logger = logging.getLogger(__name__)


TRANSLATION_TABLE = 11
PRODIGAL_START_TYPES = {"ATG": 0, "GTG": 1, "TTG": 2}

# 2-bit nucleotide encoding used by Prodigal/pyrodigal (A=0, G=1, C=2, T=3, N=0b110 -> 2)
DIGITS = {"A": 0, "G": 1, "C": 2, "T": 3}
N_DIGIT = 2


class OrientedContig:
    """A contig strand: the sequence read 5'->3' on that strand plus its digit encoding."""

    def __init__(self, seq):
        self.seq = seq
        self.digits = [DIGITS.get(base, N_DIGIT) for base in seq]
        self.length = len(seq)

    def mer(self, i, length):
        ndx = 0
        for j in range(length):
            ndx |= self.digits[i + j] << (2 * j)
        return ndx

    def base_is(self, i, base):
        return 0 <= i < self.length and self.seq[i] == base


class Parameters:
    """Plain Python copy of the pyrodigal TrainingInfo values needed for scoring."""

    def __init__(self, tinf):
        self.gc = tinf.gc
        self.st_wt = tinf.start_weight
        self.uses_sd = tinf.uses_sd
        self.no_mot = tinf.missing_motif_weight
        self.type_wt = list(tinf.type_weights)
        self.rbs_wt = list(tinf.rbs_weights)
        self.gene_dc = list(tinf.coding_statistics)
        self.ups_comp = [list(row) for row in tinf.upstream_compositions.tolist()]
        self.mot_wt = tinf.motif_weights.tolist() if not self.uses_sd else None

        no_stop = ((1 - self.gc) * (1 - self.gc) * self.gc) / 4.0
        no_stop += ((1 - self.gc) * (1 - self.gc) * (1 - self.gc)) / 8.0
        self.no_stop = 1 - no_stop
        self.lfac_max = math.log((1 - self.no_stop ** 1000.0) / self.no_stop ** 1000.0)
        self.lfac_min = math.log((1 - self.no_stop ** 80) / self.no_stop ** 80)


# ----------------------------------------------------------------------------
# Prodigal scoring functions (ported from pyrodigal lib.pyx / Prodigal node.c)
# ----------------------------------------------------------------------------

def length_factor(start, stop, par):
    """Length-based log-likelihood of an ORF of this size not occurring by chance."""
    gsize = (stop - start + 3.0) / 3.0
    if gsize > 1000.0:
        return (par.lfac_max - par.lfac_min) * (gsize - 80) / 920.0
    tmp = par.no_stop ** gsize
    return math.log((1 - tmp) / tmp) - par.lfac_min


def raw_coding(contig, start, stop, par):
    """Sum of in-frame hexamer (dicodon) log-likelihoods from start up to the stop codon."""
    return sum(par.gene_dc[contig.mer(j, 6)] for j in range(start, stop, 3))


def coding_score(contig, start, stop, upstream_starts, par):
    """
    Prodigal's cscore: raw coding potential, penalized if a longer (upstream) start of the
    same ORF has more coding potential, plus the (penalized) length factor.
    """
    cscore = raw_coding(contig, start, stop, par)
    lfac = length_factor(start, stop, par)

    if upstream_starts:
        max_upstream = max(raw_coding(contig, s, stop, par) for s in upstream_starts)
        if cscore <= max_upstream:
            cscore -= max_upstream - cscore
        max_lfac = max(length_factor(s, stop, par) for s in upstream_starts)
        if lfac <= max_lfac:
            lfac -= max(min(max_lfac - lfac, lfac), 0)

    if lfac > 3.0 and cscore < 0.5 * lfac:
        cscore = 0.5 * lfac
    return cscore + lfac


def shine_dalgarno_exact(contig, pos, start, rbs_wt):
    match = [-10] * 6
    limit = min(6, start - 4 - pos)
    for i in range(limit):
        if 0 <= pos + i < contig.length:
            if i % 3 == 0:
                if contig.base_is(pos + i, "A"):
                    match[i] = 2
            elif contig.base_is(pos + i, "G"):
                match[i] = 3

    max_val = 0
    for i in range(limit, 2, -1):
        for j in range(limit + 1 - i):
            cur_ctr = -2 + sum(match[j:j+i])
            if cur_ctr < 6:
                continue
            rdis = start - (pos + j + i)
            if rdis < 5:
                dis_flag = 2 if i < 5 else 1
            elif rdis < 11:
                dis_flag = 0
            elif rdis < 13:
                dis_flag = 1 if i < 5 else 2
            elif rdis < 16:
                dis_flag = 3
            else:
                continue
            table = {
                6: (13, 6, 1, 2),     # GGA
                8: (15, 12, 11, 3),   # AGGA
                9: (16, 12, 11, 3),   # GGAG
                11: (22, 21, 20, 10), # AGGAG
                12: (24, 23, 20, 10), # GGAGG
                14: (27, 26, 25, 10), # AGGAGG
            }
            cur_val = table[cur_ctr][dis_flag] if cur_ctr in table else 0
            if rbs_wt[cur_val] < rbs_wt[max_val]:
                continue
            if rbs_wt[cur_val] == rbs_wt[max_val] and cur_val < max_val:
                continue
            max_val = cur_val
    return max_val


def shine_dalgarno_mm(contig, pos, start, rbs_wt):
    match = [0] * 6
    limit = min(6, start - 4 - pos)
    for i in range(limit):
        if 0 <= pos + i < contig.length:
            if i % 3 == 0:
                match[i] = 2 if contig.base_is(pos + i, "A") else -3
            else:
                match[i] = 3 if contig.base_is(pos + i, "G") else -2
        else:
            match[i] = -3 if i % 3 == 0 else -2

    max_val = 0
    for i in range(limit, 4, -1):
        for j in range(limit + 1 - i):
            cur_ctr = -2
            mism = 0
            for k in range(j, j + i):
                cur_ctr += match[k]
                if match[k] < 0:
                    mism += 1
                    if k <= j + 1 or k >= j + i - 2:
                        cur_ctr -= 10
            if mism != 1 or cur_ctr < 6:
                continue
            rdis = start - (pos + j + i)
            if rdis < 5:
                dis_flag = 1
            elif rdis < 11:
                dis_flag = 0
            elif rdis < 13:
                dis_flag = 2
            elif rdis < 16:
                dis_flag = 3
            else:
                continue
            table = {
                6: (9, 5, 4, 2),      # AGxAG
                7: (14, 8, 7, 2),     # GGxGG
                9: (19, 18, 17, 3),   # AGGxGG
            }
            cur_val = table[cur_ctr][dis_flag] if cur_ctr in table else 0
            if rbs_wt[cur_val] < rbs_wt[max_val]:
                continue
            if rbs_wt[cur_val] == rbs_wt[max_val] and cur_val < max_val:
                continue
            max_val = cur_val
    return max_val


def best_upstream_motif(contig, start, par):
    """Score of the best trained upstream motif (used instead of SD if the genome doesn't use SD)."""
    max_sc = -100.0
    for i in reversed(range(4)):
        for j in range(start - 18 - i, start - 5 - i):
            if j < 0:
                continue
            if j <= start - 16 - i:
                spacendx = 3
            elif j <= start - 14 - i:
                spacendx = 2
            elif j >= start - 7 - i:
                spacendx = 1
            else:
                spacendx = 0
            score = par.mot_wt[i][spacendx][contig.mer(j, i + 3)]
            if score > max_sc:
                max_sc = score
    if max_sc == -4.0 or max_sc < par.no_mot + 0.69:
        return par.no_mot
    return max_sc


def rbs_score(contig, start, par):
    rbs = [0, 0]
    for j in range(start - 20, start - 5):
        if j < 0:
            continue
        rbs[0] = max(rbs[0], shine_dalgarno_exact(contig, j, start, par.rbs_wt))
        rbs[1] = max(rbs[1], shine_dalgarno_mm(contig, j, start, par.rbs_wt))
    sd_score = max(par.rbs_wt[rbs[0]], par.rbs_wt[rbs[1]]) * par.st_wt
    if par.uses_sd:
        return sd_score
    rscore = par.st_wt * best_upstream_motif(contig, start, par)
    if rscore < sd_score and par.no_mot > -0.5:
        rscore = sd_score
    return rscore


def upstream_score(contig, start, par):
    """Base composition at -1,-2 and -15..-44 relative to the start codon."""
    uscore = 0.0
    count = 0
    for i in list(range(1, 3)) + list(range(15, 45)):
        if i > start:
            break
        uscore += 0.4 * par.st_wt * par.ups_comp[count][contig.mer(start - i, 1)]
        count += 1
    return uscore


def type_score(start_codon, par):
    """Start codon weight; rare start codons get the weight of the nearest normal start codon."""
    return par.type_wt[PRODIGAL_START_TYPES[normal_start(start_codon)]] * par.st_wt


def score_orf(contig, start, stop, par):
    """Score one ORF (start/stop: 0-based first base of start/stop codon on the oriented contig)."""

    # Alternative Prodigal start nodes of the same ORF that lie upstream (between the
    # previous in-frame stop and this start).
    upstream_starts = []
    pos = start - 3
    while pos >= 0 and contig.seq[pos:pos+3] not in STOP_CODONS:
        if contig.seq[pos:pos+3] in PRODIGAL_START_TYPES:
            upstream_starts.append(pos)
        pos -= 3

    cscore = coding_score(contig, start, stop, upstream_starts, par)
    tscore = type_score(contig.seq[start:start+3], par)
    rscore = rbs_score(contig, start, par)
    uscore = upstream_score(contig, start, par)

    # Penalize genes < 250bp: shrink positive and inflate negative start components
    orf_length = stop - start
    if orf_length < 250:
        negf = 250.0 / orf_length
        posf = orf_length / 250.0
        rscore *= negf if rscore < 0 else posf
        uscore *= negf if uscore < 0 else posf
        tscore *= negf if tscore < 0 else posf

    sscore = tscore + rscore + uscore
    if cscore < 0.0:
        sscore -= 0.5

    return {
        "coding": cscore,
        "rbs": rscore,
        "upstream": uscore,
        "type": tscore,
        "start": sscore,
        "total": cscore + sscore,
    }


def score(orfs, genome, **kwargs):
    """Prodigal gene score (coding + start score) of every ORF."""
    logger.info("Training pyrodigal on the genome...")
    finder = pyrodigal.GeneFinder(closed=True)
    tinf = finder.train(*genome.contigs.values(), translation_table=TRANSLATION_TABLE)
    par = Parameters(tinf)
    logger.debug(f"Training: GC={par.gc:.3f}, start weight={par.st_wt:.2f}, uses SD={par.uses_sd}")

    contigs = {}
    scores = {}
    for orf in orfs:
        key = (orf.contig, orf.strand)
        if key not in contigs:
            contigs[key] = OrientedContig(orf.seq)
        result = score_orf(contigs[key], orf.start, orf.stop, par)
        logger.debug(f"prodigal {orf.header}: " + ", ".join(f"{k}={v:.3f}" for k, v in result.items()))
        scores[orf.header] = result["total"]
    return scores
