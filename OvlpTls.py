#!/usr/bin/env python3

# OverlapTools - https://github.com/Kani0n/OverlapTools
# Copyright (C) 2026 Korbinian Pürckhauer
# SPDX-License-Identifier: GPL-3.0-or-later

import argparse, logging, os, subprocess, sys
from src import classify, scores, homologs, align, tools, blastdb


def validate_file(filepath, param_name):
    """Validate that a file exists."""
    if not os.path.isfile(filepath):
        raise ValueError(f"File does not exist for {param_name}: {filepath}")


def validate_dir(dirpath, param_name):
    """Validate that a directory exists."""
    if not os.path.isdir(dirpath):
        raise ValueError(f"Directory does not exist for {param_name}: {dirpath}")


def validate_blastdb(path):
    """Validate that a BLAST nucleotide database exists at the path prefix."""
    if not any(os.path.isfile(path + ext) for ext in ('.ndb', '.nal', '.nin')):
        raise ValueError(f"No BLAST nucleotide database found at {path} (--blastdb); build it with build_blastdb (see SETUP.md)")


def validate_file_or_dir(path, param_name):
    """Validate that a path exists as either a file or directory."""
    if not os.path.isfile(path) and not os.path.isdir(path):
        raise ValueError(f"File/Directory does not exist for {param_name}: {path}")


DEFAULT_BLASTDB = 'GTDBtax_r232'   # in <db-dir>/blast/


def create_parser():
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        description="OverlapTools - Annotate ORFs by overlap with GTF gene ranges.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Positional argument for analysis
    parser.add_argument(
        "analysis",
        choices=['classify_ORFs', 'score_ORFs', 'find_homologs', 'align_codons', 'run_RNAcode', 'run_OLGenie', 'run_all', 'build_blastdb'],
        help="Analysis to run; run_all runs classify_ORFs, score_ORFs (if --score is given), find_homologs, "
             "align_codons, run_RNAcode and run_OLGenie one after the other"
    )

    # Unified input argument
    parser.add_argument(
        "--input", "-i",
        help="Input file or directory"
    )
    
    # Global optional arguments
    parser.add_argument(
        "--out-dir", "-o",
        default="./out",
        help="Output directory (default: ./out)"
    )
    parser.add_argument(
        "--log-level",
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help="Logging level (default: INFO)"
    )
    parser.add_argument(
        "--blastdb", "-b",
        help="Path (prefix) of the BLAST database for find_homologs and align_codons; output path for build_blastdb "
             f"(default: <db-dir>/blast/{DEFAULT_BLASTDB})"
    )
    parser.add_argument(
        "--taxdump", "-t",
        help="TaxonKit taxonomy dump directory (nodes.dmp, names.dmp, ...) matching the taxids of the BLAST database (required for find_homologs and build_blastdb)"
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="build_blastdb: keep the combined genome FASTA and the sequence -> taxid map after building the database"
    )
    
    # classify_ORFs arguments
    parser.add_argument(
        "--mode", "-m",
        choices=['REP', 'DR', 'CUSTOM'],
        help="Classification mode for classify_ORFs"
    )
    parser.add_argument(
        "--gtf-path",
        help="Path to genome annotation .gtf file for classify_ORFs"
    )
    parser.add_argument(
        "--n-orfs",
        type=int,
        default=100,
        help="Number of ORFs to consider for classify_ORFs (default: 100)"
    )

    # score ORFs arguments
    parser.add_argument(
        "--score",
        nargs="+",
        choices=scores.METHODS,
        help="Scoring method(s) for score_ORFs, e.g. --score prodigal balrog"
    )
    parser.add_argument(
        "--genome",
        help="Genome FASTA the ORFs come from (required for --score prodigal and balrog)"
    )
    parser.add_argument(
        "--db-dir",
        default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'db'),
        help="Directory for the scoring models and databases (downloaded on first use) and the default BLAST database (default: OverlapTools/db)"
    )
    
    # find_homologs arguments
    parser.add_argument(
        "--n-threads",
        type=int,
        default=5,
        help="Number of threads for tBLASTn (find_homologs) and the pgap/bakta scores (score_ORFs) (default: 5)"
    )
    parser.add_argument(
        "--level",
        choices=[
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
        default=None,
        help="Restrict the find_homologs search to the reference's taxon at this level (default: whole database).",
    )
    parser.add_argument(
        "--taxid",
        type=int,
        default=None,
        help=(
            "Reference organism taxid for the --level search restriction and the taxonomic context column. "
            "If omitted, infer it from the top BLAST hit."
        ),
    )

    # align_codons arguments
    parser.add_argument(
        "--n-homologs",
        type=int,
        default=100,
        help="Maximum number of homologs per ORF for align_codons (default: 100)"
    )
    
    return parser


def validate_args(args):
    """Validate arguments based on the analysis being run."""

    if args.analysis == 'classify_ORFs':
        if not args.mode:
            raise ValueError("--mode is required for classify_ORFs")
        if not args.input:
            raise ValueError("--input is required once for classify_ORFs")
        if not args.gtf_path:
            raise ValueError("--gtf-path is required for classify_ORFs")
        validate_file(args.input, "--input (ORFs file)")
        validate_file(args.gtf_path, "--gtf-path")
    
    elif args.analysis == 'score_ORFs':
        if not args.input:
            raise ValueError("--input is required once for score_ORFs")
        if not args.score:
            raise ValueError("--score is required for score_ORFs")
        needs_genome = scores.NEEDS_GENOME.intersection(args.score)
        if needs_genome and not args.genome:
            raise ValueError(f"--genome is required for --score {' '.join(sorted(needs_genome))}")
        validate_file(args.input, "--input (classified ORFs)")
        if args.genome:
            validate_file(args.genome, "--genome")

    elif args.analysis == 'find_homologs':
        if not args.taxdump:
            raise ValueError("--taxdump is required for find_homologs")
        if not args.input:
            raise ValueError("--input is required once for find_homologs")
        validate_file(args.input, "--input (classified ORFs)")
        validate_blastdb(args.blastdb)
        validate_dir(args.taxdump, "--taxdump")
        validate_file(os.path.join(args.taxdump, "nodes.dmp"), "--taxdump (nodes.dmp)")
    
    elif args.analysis == 'align_codons':
        if not args.input:
            raise ValueError("--input is required once for align_codons")
        validate_file(args.input, "--input (BLAST results)")
        validate_blastdb(args.blastdb)
    
    elif args.analysis == 'build_blastdb':
        missing = [flag for flag, value in {'--input': args.input, '--taxdump': args.taxdump}.items() if not value]
        if missing:
            raise ValueError(f"{', '.join(missing)} required for build_blastdb")
        validate_dir(args.input, "--input (genome FASTA directory)")
        validate_file(os.path.join(args.taxdump, "taxid.map"), "--taxdump (taxid.map)")

    elif args.analysis == 'run_all':
        required = {'--input': args.input, '--mode': args.mode, '--gtf-path': args.gtf_path,
                    '--taxdump': args.taxdump}
        missing = [flag for flag, value in required.items() if not value]
        if missing:
            raise ValueError(f"{', '.join(missing)} required for run_all")
        if args.score:
            needs_genome = scores.NEEDS_GENOME.intersection(args.score)
            if needs_genome and not args.genome:
                raise ValueError(f"--genome is required for --score {' '.join(sorted(needs_genome))}")
        validate_file(args.input, "--input (ORFs file)")
        validate_file(args.gtf_path, "--gtf-path")
        if args.genome:
            validate_file(args.genome, "--genome")
        validate_blastdb(args.blastdb)
        validate_dir(args.taxdump, "--taxdump")
        validate_file(os.path.join(args.taxdump, "nodes.dmp"), "--taxdump (nodes.dmp)")

    elif args.analysis in ['run_RNAcode', 'run_OLGenie']:
        if not args.input:
            raise ValueError(f"--input is required once for {args.analysis}")
        validate_file_or_dir(args.input, "--input (alignment file/directory)")


def OLGenie_scripts():
    OLGenie_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'src', 'OLGenie')
    return os.path.join(OLGenie_dir, 'OLGenie.pl'), os.path.join(OLGenie_dir, 'OLGenie_bootstrap.R')


def run_all(args, out_dir):
    """Run the whole pipeline, each step on the output of the previous one."""
    logger = logging.getLogger("OverlapTools")
    steps = 6 if args.score else 5
    step = 0

    def next_step(name):
        nonlocal step
        step += 1
        logger.info(f"Step {step}/{steps}: {name}")

    next_step("classify_ORFs")
    classified = classify.classify_ORFs(args.mode, args.input, args.gtf_path, args.n_orfs, out_dir)

    if args.score:
        next_step("score_ORFs")
        scores.score_ORFs(classified, args.score, out_dir, args.genome, args.db_dir, args.n_threads)

    next_step("find_homologs")
    blast_tsv = homologs.find_homologs(classified, args.blastdb, args.taxdump, args.n_threads, out_dir, args.level, args.taxid)

    next_step("align_codons")
    align_dir = align.align_codons(blast_tsv, args.blastdb, args.n_homologs, out_dir)

    next_step("run_RNAcode")
    tools.run_RNAcode(align_dir, out_dir)

    next_step("run_OLGenie")
    tools.run_OLGenie(*OLGenie_scripts(), align_dir, out_dir)


def run_analysis(args, out_dir):
    """Run the requested analysis."""
    if args.analysis == 'run_all':
        run_all(args, out_dir)

    elif args.analysis == 'build_blastdb':
        blastdb.build_blastdb(args.input, args.taxdump, args.blastdb, args.keep_intermediate)

    elif args.analysis == 'classify_ORFs':
        classify.classify_ORFs(args.mode, args.input, args.gtf_path, args.n_orfs, out_dir)

    elif args.analysis == 'score_ORFs':
        scores.score_ORFs(args.input, args.score, out_dir, args.genome, args.db_dir, args.n_threads)

    elif args.analysis == 'find_homologs':
        homologs.find_homologs(args.input, args.blastdb, args.taxdump, args.n_threads, out_dir, args.level, args.taxid)
    
    elif args.analysis == 'align_codons':
        align.align_codons(args.input, args.blastdb, args.n_homologs, out_dir)
    
    elif args.analysis == 'run_RNAcode':
        tools.run_RNAcode(args.input, out_dir)
    
    elif args.analysis == 'run_OLGenie':
        tools.run_OLGenie(*OLGenie_scripts(), args.input, out_dir)


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s [%(levelname)s] %(module)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    logger = logging.getLogger("OverlapTools")
    
    if not args.blastdb:
        args.blastdb = os.path.join(args.db_dir, 'blast', DEFAULT_BLASTDB)

    try:
        validate_args(args)
    except ValueError as e:
        parser.error(str(e))
    
    out_dir = args.out_dir
    if args.analysis == 'build_blastdb':
        logger.info(f"Running build_blastdb, database: {args.blastdb}")
    else:
        os.makedirs(out_dir, exist_ok=True)
        logger.info(f"Running {args.analysis}, output directory: {out_dir}")
    logger.debug(f"Arguments: {vars(args)}")

    try:
        run_analysis(args, out_dir)
    except (subprocess.CalledProcessError, RuntimeError, ValueError, FileNotFoundError) as e:
        logger.error(str(e), exc_info=logger.isEnabledFor(logging.DEBUG))
        sys.exit(1)

    logger.info(f"Finished {args.analysis}")
