#!/usr/bin/env python3

"""
Bakta protein homology evidence for arbitrary ORFs.

Bakta calls its CDS with pyrodigal (see the prodigal score) and annotates them by
homology to its protein database. Every ORF protein is annotated with
bakta_proteins (light database: UniRef50 protein clusters, expert systems such as
AMRFinderPlus and curated BlastRules, Pfam); the score is the best alignment
bitscore of the evidence Bakta uses for the annotation (protein cluster / expert
hits), 0 if Bakta calls the protein hypothetical. Rare start codons are translated
as Met. The database is downloaded once into the database directory.
"""

import os, sys, json, logging, subprocess


logger = logging.getLogger(__name__)


DB_TYPE = "light"


def bakta_env():
    """Environment with the conda env's bin directory on PATH (diamond, hmmer, amrfinder)."""
    bin_dir = os.path.dirname(sys.executable)
    return dict(os.environ, PATH=bin_dir + os.pathsep + os.environ.get("PATH", ""))


def database_path(db_dir):
    bakta_dir = os.path.join(db_dir, "bakta")
    path = os.path.join(bakta_dir, f"db-{DB_TYPE}")
    if not os.path.isfile(os.path.join(path, "version.json")):
        os.makedirs(bakta_dir, exist_ok=True)
        logger.info(f"Downloading Bakta {DB_TYPE} database (~1.5 GB) to {path}...")
        subprocess.run(["bakta_db", "download", "--output", bakta_dir, "--type", DB_TYPE],
                       stdout=subprocess.DEVNULL, check=True, env=bakta_env())
    return path


def evidence_score(feature):
    """Best bitscore of the alignments Bakta based the annotation on."""
    if feature.get("hypothetical"):
        return 0.0
    evidence = [feature.get("psc"), feature.get("pscc")] + feature.get("expert", [])
    scores = [e["score"] for e in evidence if e and e.get("score") is not None]
    return max(scores) if scores else None


def score(orfs, genome, db_dir, n_threads=1, tmp_dir="tmp", **kwargs):
    """Bakta annotation evidence bitscore of every ORF protein."""
    db = database_path(db_dir)
    bakta_dir = os.path.join(tmp_dir, "bakta")
    os.makedirs(bakta_dir, exist_ok=True)

    proteins = os.path.join(bakta_dir, "orfs.faa")
    with open(proteins, "w") as f:
        for i, orf in enumerate(orfs):
            f.write(f">{i}\n{orf.initiator_protein}\n")

    logger.info(f"Annotating {len(orfs)} ORF proteins with bakta_proteins...")
    subprocess.run(
        ["bakta_proteins", "--db", db, "--output", bakta_dir, "--prefix", "orfs",
         "--force", "--threads", str(n_threads), proteins],
        stdout=subprocess.DEVNULL, check=True, env=bakta_env()
    )

    with open(os.path.join(bakta_dir, "orfs.json")) as f:
        features = {feature["id"]: feature for feature in json.load(f)["features"]}

    scores = {}
    for i, orf in enumerate(orfs):
        feature = features.get(str(i), {})
        scores[orf.header] = evidence_score(feature)
        logger.debug(f"bakta {orf.header}: {feature.get('product')}, score={scores[orf.header]}")
    return scores
