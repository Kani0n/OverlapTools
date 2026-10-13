import os, subprocess, logging

logger = logging.getLogger(__name__)


MIN_RNACODE_SEQUENCES = 3


def count_aln_sequences(codons_aln):
    """Number of sequences in a CLUSTAL alignment (conservation lines start with a blank)."""
    names = set()
    with open(codons_aln, 'r') as aln:
        for line in aln:
            if line.strip() and not line.startswith(('CLUSTAL', ' ', '\t')):
                names.add(line.split()[0])
    return len(names)


def run_RNAcode_single(codons_aln, rnacode_out_dir):
    rnacode_outfile = os.path.join(rnacode_out_dir, os.path.basename(codons_aln).replace('.aln', '.txt'))
    n_sequences = count_aln_sequences(codons_aln)
    if n_sequences < MIN_RNACODE_SEQUENCES:
        logger.warning(f'Skipping {os.path.basename(codons_aln)}: RNAcode needs at least {MIN_RNACODE_SEQUENCES} sequences, the alignment has {n_sequences}')
        if os.path.exists(rnacode_outfile):
            os.remove(rnacode_outfile)   # left over from a previous run
        return False
    with open(rnacode_outfile, 'w') as out:
        result = subprocess.run(['RNAcode', codons_aln], stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        # a failing alignment should not stop the whole run
        reason = f'signal {-result.returncode}' if result.returncode < 0 else f'exit {result.returncode}'
        logger.warning(f'RNAcode failed on {os.path.basename(codons_aln)} ({reason}): {result.stderr.strip() or "no error message"}')
        os.remove(rnacode_outfile)
        return None
    return True


def run_RNAcode(codons_aln_path, out_dir):
    rnacode_out_dir = os.path.join(out_dir, 'RNAcode')
    if not os.path.exists(rnacode_out_dir):
        os.mkdir(rnacode_out_dir)
    if os.path.isfile(codons_aln_path):
        run_RNAcode_single(codons_aln_path, rnacode_out_dir)
    else:
        alignments = [a for a in os.listdir(codons_aln_path) if a.endswith('.aln')]
        results = [run_RNAcode_single(os.path.join(codons_aln_path, a), rnacode_out_dir) for a in alignments]
        n_skipped, n_failed = results.count(False), results.count(None)
        logger.info(f'Ran RNAcode on {results.count(True)} / {len(alignments)} alignments '
                    f'({n_skipped} skipped with fewer than {MIN_RNACODE_SEQUENCES} sequences, {n_failed} failed)')


def read_info(info_file):
    with open(info_file, 'r') as f:
        for line in f:
            line = line.strip()
            key, value = line.split(':')
            if key == 'classification':
                if value == 'intergenic':
                    return value
            elif key == 'frameshift':
                return value


def convert_aln_to_fasta(codons_aln, orf_outdir):
    seqs = {}
    with open(codons_aln, 'r') as aln:
        for line in aln:
            line = line.rstrip('\n')
            if not line or line.startswith('CLUSTAL'):
                continue
            elements = line.split()
            header = elements[0]
            seq = elements[1]
            if header not in seqs.keys():
                seqs[header] = ''
            seqs[header] += seq
    codons_fasta = os.path.join(orf_outdir, os.path.basename(codons_aln).replace('.aln', '.fasta'))
    with open(codons_fasta, 'w') as fasta:
        for header, seq in seqs.items():
            fasta.write(f'>{header}\n{seq}\n')
    return codons_fasta


def run_OLGenie_and_bootstrap(OLGenie, OLGenie_bootstrap, codons_fasta, frame_shift, olgenie_outfile):
    subprocess.run([OLGenie, f'--fasta_file={codons_fasta}', f'--frame={frame_shift}', f'--output_file={olgenie_outfile}'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    olgenie_bootstrap_outfile = olgenie_outfile.replace('.tsv', '_bootstrap.tsv')
    with open(olgenie_bootstrap_outfile, 'w') as out:
        subprocess.run(['Rscript', OLGenie_bootstrap, olgenie_outfile, '2', '1000', '4'], stdout=out, check=True)


def run_OLGenie_single(OLGenie, OLGenie_bootstrap, codons_aln, olgenie_out_dir):
    info_file = codons_aln.replace('.aln', '.info')
    if not os.path.isfile(info_file):
        raise FileNotFoundError(f"Failed to find corresponding .info file for alignment file: {codons_aln}")
    frame_shift = read_info(info_file)
    if frame_shift == '=':
        logger.info(f'Skipping same frame ORF: {codons_aln}')
        return
    orf_outdir = os.path.join(olgenie_out_dir, os.path.basename(codons_aln).strip('.aln'))
    if not os.path.exists(orf_outdir):
        os.mkdir(orf_outdir)
    codons_fasta = convert_aln_to_fasta(codons_aln, orf_outdir)
    if frame_shift == 'intergenic':
        for frame_shift in ['ss12', 'ss13', 'sas11', 'sas12', 'sas13']:
            olgenie_outfile = os.path.join(orf_outdir, f'OLGenie_{frame_shift}.tsv')
            run_OLGenie_and_bootstrap(OLGenie, OLGenie_bootstrap, codons_fasta, frame_shift, olgenie_outfile) 
    else:
        olgenie_outfile = os.path.join(orf_outdir, 'OLGenie.tsv')
        run_OLGenie_and_bootstrap(OLGenie, OLGenie_bootstrap, codons_fasta, frame_shift, olgenie_outfile)
    os.remove(codons_fasta)


def run_OLGenie(OLGenie, OLGenie_bootstrap, codons_aln_path, out_dir):
    olgenie_out_dir = os.path.join(out_dir, 'OLGenie')
    if not os.path.exists(olgenie_out_dir):
        os.mkdir(olgenie_out_dir)
    if os.path.isfile(codons_aln_path):
        run_OLGenie_single(OLGenie, OLGenie_bootstrap, codons_aln_path, olgenie_out_dir)
    else:
        for codons_aln in os.listdir(codons_aln_path):
            if codons_aln.endswith('.aln'):
                run_OLGenie_single(OLGenie, OLGenie_bootstrap, os.path.join(codons_aln_path, codons_aln), olgenie_out_dir)
