import os, shutil, subprocess, logging

logger = logging.getLogger(__name__)


def filter_blast_results(blast_csv, n_homologs):
    ORFs = {}
    with open(blast_csv, 'r') as f:
        for line in f:
            if line.startswith('qaccver'):
                continue
            elements = line.rstrip('\n').split('\t')    # tab-separated tBLASTn output (outfmt 6)
            ID = elements[0]
            start = float(elements[8])
            end = float(elements[9])
            if start < end:
                strand = 'plus' 
            else:
                start, end = end, start
                strand = 'minus'
            start = str(start).removesuffix('.0')
            end = str(end).removesuffix('.0')
            if ID in ORFs.keys():
                if len(ORFs[ID]) == n_homologs:
                    ORFs[ID].sort(key=lambda x: x[1])
                    if float(elements[10]) < ORFs[ID][-1][1]:
                        ORFs[ID][-1] = ('\t'.join([elements[1], f'{start}-{end}', strand]), float(elements[10]))
                else:
                    ORFs[ID].append(('\t'.join([elements[1], f'{start}-{end}', strand]), float(elements[10])))
            else:
                ORFs[ID] = [('\t'.join([elements[1], f'{start}-{end}', strand]), float(elements[10]))]
    return ORFs


def write_ORF_tsv(line_list, tmp_dir):
    orf_tsv = os.path.join(tmp_dir, 'homologs.tsv')
    with open(orf_tsv, 'w') as orf:
        for line, _ in line_list:
            orf.write(line+'\n')
    return orf_tsv


def run_blastdbcmd(orf_tsv, BLASTdb_path):
    orf_fasta = orf_tsv.replace('.tsv', '.fasta')
    subprocess.run(['blastdbcmd', '-db', BLASTdb_path, '-dbtype', 'nucl', '-entry_batch', orf_tsv, '-out', orf_fasta, '-outfmt', '%f'], check=True)
    return orf_fasta


def read_single_line_fasta(text):
    # parses fasta written with seqkit -w 0 (header line followed by one sequence line)
    lines = text.splitlines()
    return [(lines[i][1:], lines[i+1]) for i in range(0, len(lines) - 1, 2)]


def filter_orfs(orf_fasta, tmp_dir):
    # removes duplicate sequences | makes fasta onto 1 line | saves to ntd1.fasta
    ntd1_fasta = os.path.join(tmp_dir, 'ntd1.fasta')
    subprocess.run(['seqkit', 'rmdup', '-s', '--quiet', '-w', '0', '-o', ntd1_fasta, orf_fasta], check=True)
    # 1 line per entry aa-seqs
    prots = subprocess.run(['seqkit', 'translate', '-w', '0', ntd1_fasta], capture_output=True, text=True, check=True).stdout
    with open(ntd1_fasta, 'r') as f:
        ntds = read_single_line_fasta(f.read())
    # removes ones with stop codons | keep cleaned nuc entries
    cleaned_ntd_fasta = os.path.join(tmp_dir, 'cleaned_ntd.fasta')
    with open(cleaned_ntd_fasta, 'w') as out:
        for (header, ntd_seq), (_, prot_seq) in zip(ntds, read_single_line_fasta(prots)):
            if '*' not in prot_seq:
                out.write(f'>{header}\n{ntd_seq}\n')
    # translate back to Proteins, the frame should be the correct one from tblastn
    prot_fasta = os.path.join(tmp_dir, 'prots.fasta')
    subprocess.run(['seqkit', 'translate', '--frame', '1', '-o', prot_fasta, cleaned_ntd_fasta], check=True)
    return cleaned_ntd_fasta, prot_fasta


def run_muscle(prot_fasta, tmp_dir):
    aln_prot_fasta = os.path.join(tmp_dir, 'aln_prots.fasta')
    subprocess.run(['muscle', '-align', prot_fasta, '-output', aln_prot_fasta],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return aln_prot_fasta


def remove_duplicate_headers(aln_prot_fasta, cleaned_ntd_fasta, tmp_dir):
    # remove duplicate headers
    aln_prot_fasta2 = os.path.join(tmp_dir, 'aln_prots2.fasta')
    subprocess.run(['seqkit', 'rmdup', '--quiet', '-o', aln_prot_fasta2, aln_prot_fasta], check=True)
    # same for nuc fasta
    cleaned_ntd_fasta2 = os.path.join(tmp_dir, 'cleaned_ntd2.fasta')
    subprocess.run(['seqkit', 'rmdup', '--quiet', '-o', cleaned_ntd_fasta2, cleaned_ntd_fasta], check=True)
    return aln_prot_fasta2, cleaned_ntd_fasta2


def run_pal2nal(aln_prot_fasta2, cleaned_ntd_fasta2, align_path, tmp_dir, ID):
    # Run PAL2NAL, Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options], pep is the aligned file, nuc are the DNA sequences
    codons_aln = subprocess.run(['pal2nal.pl', aln_prot_fasta2, cleaned_ntd_fasta2], capture_output=True, text=True, check=True).stdout
    # : in header cause errors in RNAcode
    codons_aln_final = os.path.join(align_path, f'{ID}.aln')
    with open(codons_aln_final, 'w') as out:
        out.write(codons_aln.replace(':', '_'))
    return codons_aln_final


def get_n_homologs(fasta):
    n = 0
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


def write_info(codons_aln, ID, n_homologs):
    classification = ID.split('_')[5]
    info = f'classification:{classification}\n'
    if classification != 'intergenic':
        frame_shift = ID.split('_')[6]
        info += f'frameshift:{frame_shift}\n'
    info += f'homologs:{n_homologs}\n'
    with open(codons_aln.replace('.aln', '.info'), 'w') as f:
        f.write(info)


def remove_alignment(align_path, ID):
    """Remove the alignment files of an ORF left over from a previous run."""
    for ext in ('.aln', '.info'):
        path = os.path.join(align_path, f'{ID}{ext}')
        if os.path.exists(path):
            os.remove(path)


def process_ORFs(ORFs, align_path, BLASTdb_path):
    n = len(ORFs)
    skipped = 0
    for i, (ID, line_list) in enumerate(ORFs.items()):
        i = i+1
        logger.debug(f'ORF {i}/{n}: {ID}')

        tmp_dir = os.path.join(align_path, 'tmp')
        shutil.rmtree(tmp_dir, ignore_errors=True)   # left over if a previous run failed
        os.mkdir(tmp_dir)

        orf_tsv = write_ORF_tsv(line_list, tmp_dir)

        orf_fasta = run_blastdbcmd(orf_tsv, BLASTdb_path)
        
        cleaned_ntd_fasta, prot_fasta = filter_orfs(orf_fasta, tmp_dir)

        # an alignment needs at least two homologs (MUSCLE crashes with fewer)
        n_cleaned = get_n_homologs(prot_fasta)
        if n_cleaned < 2:
            logger.warning(f'Skipping {ID}: {n_cleaned} homolog(s) left after removing duplicates and sequences with stop codons')
            remove_alignment(align_path, ID)
            shutil.rmtree(tmp_dir)
            skipped += 1
            continue

        aln_prot_fasta = run_muscle(prot_fasta, tmp_dir)

        aln_prot_fasta2, cleaned_ntd_fasta2 = remove_duplicate_headers(aln_prot_fasta, cleaned_ntd_fasta, tmp_dir)

        codons_aln = run_pal2nal(aln_prot_fasta2, cleaned_ntd_fasta2, align_path, tmp_dir, ID)

        write_info(codons_aln, ID, get_n_homologs(aln_prot_fasta2))

        shutil.rmtree(tmp_dir)

    logger.info(f'Aligned {n - skipped} / {n} ORFs ({skipped} skipped with fewer than 2 homologs)')


def align_codons(BLAST_res_path, BLASTdb_path, n_homologs, out_dir):
    align_path = os.path.join(out_dir, 'align')
    if not os.path.exists(align_path):
        os.mkdir(align_path)
    ORFs = filter_blast_results(BLAST_res_path, n_homologs)
    process_ORFs(ORFs, align_path, BLASTdb_path)
    return align_path
