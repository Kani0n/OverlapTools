import os, logging

logger = logging.getLogger(__name__)


def parse_attributes(attribute_text):
    attrs = {}
    for item in attribute_text.split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
        else:
            key, value = item, ""
        value = value.strip().strip('"')
        attrs[key] = value
    return attrs


def parse_gtf(gtf):
    genes_by_id = {}

    with open(gtf, 'r', encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom, _, feature, start, end, _, strand, _, attr_text = parts
            attrs = parse_attributes(attr_text)
            gene_id = attrs.get("gene_id") or attrs.get("Parent") or attrs.get("ID") or attrs.get("locus_tag")
            if not gene_id:
                logger.warning(f'Failed to extract gene_id from line: {line}')
                continue
            if feature == "gene":
                pseudogene = attrs.get("gene_biotype") == "pseudogene"
                gene_record = {
                    "chrom": chrom,
                    "strand": strand,
                    "cds_ranges": [],
                    "pseudogene": pseudogene
                }
                genes_by_id[gene_id] = gene_record
            elif feature == "CDS" and gene_id in genes_by_id:
                if strand == "+":
                    cds_start = int(start)
                    cds_end = int(end)
                elif strand == "-":
                   cds_start = int(start)
                   cds_end = int(end)
                incomplete = "incomplete" in attrs.keys() or attrs.get("note") == "incomplete"
                cds_record = {
                    "start": cds_start,
                    "end": cds_end,
                    "incomplete": incomplete
                }
                genes_by_id[gene_id]["cds_ranges"].append(cds_record)
    return genes_by_id


def swap_frame(frame):
    frame_map = {
        'ss12': 'ss13',
        'ss13': 'ss12',
        'sas11': 'sas11',
        'sas12': 'sas12',
        'sas13': 'sas13',
        '=': '='
    }
    return frame_map.get(frame)


def classify_orf(orf_start, orf_end, orf_strand, seq_id, genes):

    if orf_start > orf_end:
            orf_start, orf_end = orf_end, orf_start
    
    match = 'intergenic'
    gene_id = None
    frame = None
    incomplete = None
    pseudogene = None

    for id, gene in genes.items():
        if gene_id:
            break
        if str(gene["chrom"]) != str(seq_id):
            continue
        for cds_entry in gene["cds_ranges"]:
            cds_start = cds_entry["start"]
            cds_end = cds_entry["end"]
            if orf_end < cds_start or orf_start > cds_end:
                continue
            gene_id = id

            if orf_start == cds_start:
                if orf_end < cds_end:
                    match = 'within'
                else:
                    match = 'overlapping'
            elif orf_end == cds_end:
                if orf_start > cds_start:
                    match = 'within'
                else:
                    match = 'overlapping'
            elif orf_start > cds_start and orf_end < cds_end:
                match = 'within'
            else:
                match = 'overlapping'

            if orf_strand == '+' and gene["strand"] == '+':
                shift = (cds_start - orf_start) % 3
                if shift == 0:
                    frame = '='
                elif shift == 1:
                    frame = f'ss12'
                elif shift == 2:
                    frame = f'ss13'
            elif orf_strand == '-' and gene["strand"] == '-':
                shift = (cds_end - orf_end) % 3
                if shift == 0:
                    frame = '='
                elif shift == 1:
                    frame = f'ss13'
                elif shift == 2:
                    frame = f'ss12'
            elif orf_strand == '+' and gene["strand"] == '-':
                shift = (cds_end - orf_start) % 3
                if shift == 0:
                    frame = 'sas11'
                elif shift == 1:
                    frame = f'sas12'
                elif shift == 2:
                    frame = f'sas13'
            elif orf_strand == '-' and gene["strand"] == '+':
                shift = (cds_start - orf_end) % 3
                if shift == 0:
                    frame = 'sas11'
                elif shift == 1:
                    frame = f'sas13'
                elif shift == 2:
                    frame = f'sas12'
            
            incomplete = cds_entry["incomplete"]
            pseudogene = gene["pseudogene"]
            break
    
    # uncommend if you want to swap Parent-Child relationship
    frame = swap_frame(frame)
    return match, gene_id, frame, incomplete, pseudogene


def parse_fasta(path, genes):
    entries = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                elements = line.split('_')
                locus = elements[0]
                chr = locus.split(':')[0]
                start = int(locus.split(':')[1].split('-')[0])
                end = int(locus.split(':')[1].split('-')[1])
                strand = elements[1]

                if start > end:
                    start, end = end, start
                # give start codon coordinade inclusive, but stop codon exclusive
                if strand == '+':
                    end = end - 1       # exclude stop codon for gtf compatibility
                elif strand == '-':
                    start = start + 1   # exclude stop codon for gtf compatibility
                
                classification, gene_id, frame_shift, incomplete, pseudogene = classify_orf(start, end, strand, chr, genes)
            else:
                entries.append((locus, strand, classification, frame_shift, gene_id, incomplete, pseudogene, line))
    return entries


def parse_reparation_fasta(path, genes):
    entries = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                elements = line.split('|')
                locus = elements[1]
                chr = locus.split(':')[0]
                start = int(locus.split(':')[1].split('-')[0])
                end = int(locus.split(':')[1].split('-')[1])
                strand = elements[2].split(' ')[2].split(':')[1]

                if start > end:
                    start, end = end, start
                # REP gives start codon coordinade inclusive, but stop codon exclusive
                if strand == '+':
                    end = end - 1       # exclude stop codon for gtf compatibility
                elif strand == '-':
                    start = start + 1   # exclude stop codon for gtf compatibility
                
                classification, gene_id, frame_shift, incomplete, pseudogene = classify_orf(start, end, strand, chr, genes)
            else:
                entries.append((locus, strand, classification, frame_shift, gene_id, incomplete, pseudogene, line))
    return entries


def parse_deepribo_csv(path, genes, n):
    entries = []
    with open(path, 'r') as csv:
        for line in csv:
            line = line.strip()
            if line.startswith(',filename'):
                continue
            elements = line.split(',')

            # discard annotated genes
            label = elements[3]
            if label == 'True':
                continue
            # discard non-SS ORFs
            SS_pred_rank = elements[21].strip()
            if SS_pred_rank == '999999':
                continue

            # used for ranking the top n predictions
            pred_rank = elements[18]

            locus = elements[14]
            chr = locus.split(':')[0]
            start = int(locus.split(':')[1].split('-')[0])
            end = int(locus.split(':')[1].split('-')[1])
            strand = elements[5]
            seq = elements[15].strip('*')

            if start > end:
                start, end = end, start
            # DR gives start codon coordinade inclusive, but stop codon exclusive
            if strand == '+':
                end = end - 1       # exclude stop codon for gtf compatibility
            elif strand == '-':
                start = start + 1   # exclude stop codon for gtf compatibility

            classification, gene_id, frame_shift, incomplete, pseudogene = classify_orf(start, end, strand, chr, genes)

            if len(entries) < n: # Top n from DR
                entries.append((locus, strand, classification, frame_shift, gene_id, incomplete, pseudogene, seq, float(pred_rank)))
            else:
                entries.sort(key=lambda x: x[-1])
                if float(pred_rank) < entries[-1][-1]:
                    entries[-1] = (locus, strand, classification, frame_shift, gene_id, incomplete, pseudogene, seq, float(pred_rank))
    return [(locus, strand, classification, frame_shift, gene_id, incomplete, pseudogene, seq) for (locus, strand, classification, frame_shift, gene_id, incomplete, pseudogene, seq, pred_rank) in entries] 


def write_output(ORFS, out_file):
    with open(out_file, 'w') as fasta:
        for (locus, strand, classification, frame_shift, gene_id, incomplete, pseudogene, seq) in ORFS:
            length = len(seq)
            header = f'>{locus.replace(':', '_')}_{strand}_{length}_{classification}'
            if classification != "intergenic":
                header += f'_{frame_shift}_{gene_id}'
            if pseudogene:
                header += f'_pseudogene'
            if incomplete:
                header += f'_incomplete'
            fasta.write(f'{header}\n{seq}\n')


def classify_ORFs(mode, orfs_path, gtf_path, n_ORFs, out_dir):
    genes = parse_gtf(gtf_path)
    if mode == "REP":
        # Takes REPARATION .fasta
        ORFS = parse_reparation_fasta(orfs_path, genes)
    elif mode == "DR":
        # Takes DeepRibo .csv
        ORFS = parse_deepribo_csv(orfs_path, genes, n=n_ORFs)
    elif mode == "CUSTOM":
        # Takes standardized .fasta
        # >CHROMOSOME:START-END_STRAND
        ORFS = parse_fasta(orfs_path, genes)
    classified_path = os.path.join(out_dir, 'classified_orfs')
    if not os.path.exists(classified_path):
        os.mkdir(classified_path)
    classified_file = os.path.join(classified_path, 'classified_ORFs.fasta')
    write_output(ORFS, classified_file)
    return classified_file