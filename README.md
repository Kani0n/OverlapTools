# OverlapTools (OvlpTls)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

OverlapTools is a comprehensive Python pipeline for annotating and analyzing ORFs in genomic sequences. It provides a complete workflow for classifying ORFs by their overlap with known genes, finding homologous sequences, aligning codons, and analyzing coding potential.

## Overview

OverlapTools processes ORF sequences through a multi-step analysis pipeline:

1. **Classify ORFs** - Categorize ORFs based on their overlap with annotated genes
   - *(optional)* **Score ORFs** - Gene-likeness scores (Prodigal, Balrog, PGAP, Bakta) for every classified ORF
2. **Find Homologs** - Search for homologous sequences in BLAST database and annotate/filter them by taxonomy
3. **Align Codons** - Create high-quality codon alignments from homologous sequences
4. **Run RNAcode** - Detect coding regions in alignments using RNAcode
5. **Run OLGenie** - Analyze frameshift signals in overlapping genes using OLGenie

`run_all` runs all steps in one call.

---

## Installation

OverlapTools runs on Linux (x86-64) and is used directly from the cloned repository:

```bash
git clone https://github.com/Kani0n/OverlapTools.git
cd OverlapTools
```

All dependencies (Python, pyrodigal, PyTorch, pyhmmer, Bakta, BLAST+, TaxonKit, seqkit, MUSCLE, PAL2NAL, RNAcode, Perl and R incl. the R packages needed by the OLGenie bootstrap script) are listed in `env.yml`:

```bash
conda env create -f env.yml     # or: mamba env create -f env.yml
conda activate OverlapTools
python3 OvlpTls.py --help
```

OverlapTools also needs data that is not part of the repository; see **[SETUP.md](SETUP.md)**:

- `score_ORFs`: models and databases (Balrog models, PGAP HMM library ~3 GB, Bakta light database ~1.5 GB download), downloaded automatically on first use into `OverlapTools/db/` (or `--db-dir`) or beforehand
- `find_homologs`: a taxonomy dump for TaxonKit (e.g. the GTDB taxdump, given with `--taxdump`) and a BLAST nucleotide database whose taxids belong to that taxonomy, built with `build_blastdb` from the GTDB representative genomes

---

## Global Options

- `--out-dir`, `-o` - Output directory (default: ./out)
- `--log-level` - Logging level: `DEBUG`, `INFO`, `WARNING` or `ERROR` (default: `INFO`). Log messages are written to stderr; with `DEBUG` the full traceback is shown on errors.

---

## Analyses

### 1. `classify_ORFs`

**Purpose:** Classify ORFs by their overlap relationship with known genes from a GTF annotation file.

**Input Types:**
- **REP** - Reparation
- **DR** - DeepRibo
- **CUSTOM** - Custom FASTA (>Chr_start-end)

**How it works:**
- Parses GTF gene annotation file to extract CDS ranges
- Compares each ORF start/end positions against CDS ranges
- Classifies relationships as:
  - `within` - ORF is completely contained within a CDS
  - `overlapping` - ORF partially overlaps with a CDS
  - `intergenic` - ORF is between CDS regions
- Detects frameshift relationships (=, ss12, ss13, sas11, sas12, sas13)
- Identifies pseudogenes and incomplete coding sequences

**Required Arguments:**
- `--input` - FASTA or csv (only in DeepRibo case) file containing ORF sequences
- `--gtf-path` - Path to GTF genome annotation file
- `--mode` - Input Type (REP, DR, or CUSTOM)
- `--out-dir` - Output directory (default: ./out)

**Optional Arguments:**
- `--n-orfs` - Number of ORFs to process (default: 100)

**Output:**
- `classified_orfs/classified_ORFs.fasta` - FASTA file with classified ORF annotations in header line
---

### `score_ORFs`

**Purpose:** Give every ORF one or more "gene-likeness" scores, including ORFs that gene finders never consider themselves (alternative start codons, short ORFs, ORFs nested in other genes).

**Scoring methods** (`--score`, several can be combined):

| Method | Score | Needs |
|---|---|---|
| `prodigal` | Prodigal gene score (coding + start score, natural-log likelihood ratio). Prodigal's node scoring re-implemented for any start/stop pair with the parameters of pyrodigal trained on the genome; reproduces pyrodigal's scores exactly for ATG/GTG/TTG starts. | `--genome` |
| `balrog` | Balrog ORF score: `(0.975·gene_prob + 0.254·TIS_prob + start_weight − 0.979) · length_nt`, with Balrog's pretrained protein gene model and translation-initiation-site model (BalrogCPP v0.5). Positive = more gene-like than Balrog's threshold. | `--genome`; models are downloaded on first use |
| `pgap` | Best bitscore of a significant hit (E ≤ 0.01) of the ORF protein against PGAP's protein family HMM library (NCBIfams, TIGRFAMs, Pfam), 0 without a hit. | HMM library (~3 GB), downloaded on first use |
| `bakta` | Best alignment bitscore of the evidence `bakta_proteins` bases the annotation on (UniRef50 protein clusters, expert systems such as AMRFinderPlus/BlastRules), 0 if Bakta calls the protein hypothetical. | Bakta light DB (~1.5 GB download), downloaded on first use |

Rare start codons get the weight of the nearest normal start codon (ATA/ATC/ATT → ATG, CTG → GTG) in `prodigal` and `balrog`, and are translated as Met for `pgap` and `bakta`.

**Required Arguments:**
- `--input` - Classified ORFs FASTA (output from classify_ORFs; ORF coordinates are parsed from the headers)
- `--score` - One or more of `prodigal`, `balrog`, `pgap`, `bakta`
- `--genome` - Genome FASTA the ORFs come from (for `prodigal` and `balrog`)
- `--out-dir` - Output directory (default: ./out)

**Optional Arguments:**
- `--db-dir` - Directory for models and databases, including the default BLAST database (default: `OverlapTools/db`)
- `--n-threads` - Threads for `pgap` and `bakta` (default: 5)

**Output:**
- `scores/ORF_scores.tsv` - ORF header and one column per selected score
---

### 2. `find_homologs`

**Purpose:** Search for homologous sequences to ORFs using tBLASTn protein search and put the hits into taxonomic context.

**How it works:**
- Without `--level`: searches the ORF amino acid sequences against the whole BLAST nucleotide database with tBLASTn (E-value cutoff 1e-5)
- With `--level`, the taxonomy restriction is part of the search itself:
  1. Reference taxid: the one given with `--taxid`; otherwise a preliminary search of the whole database (`-max_target_seqs 5`, temporary `tBLASTn_reference.tsv`) and the taxid of its overall top hit (lowest E-value, highest bitscore)
  2. TaxonKit determines the reference's taxon at `--level` and lists all taxids below it in the taxonomy dump; those present in the BLAST database are written to a temporary `search_taxids.txt`
  3. tBLASTn searches only these sequences (`-taxidlist search_taxids.txt -no_taxid_expansion`, threads split by query with `-mt_mode 1` for subsets below 280 Mbp, otherwise split over the database with at most one thread per subject sequence)
- Annotates every hit with its taxonomic context relative to the reference taxid (`--taxid` or the top hit): the most specific shared rank (`species`, `genus`, `family`, `order`, `class`, `phylum`, `kingdom`), `none` if no rank is shared, or `unclassified` if the hit's taxid could not be resolved
- Requires a pre-built BLAST nucleotide database with taxids (BLAST v5, `makeblastdb -parse_seqids -taxid_map`) whose taxids belong to the TaxonKit taxonomy dump

**Notes on the restricted search:**
- E-values are computed for the size of the searched subset, so the same alignment gets a smaller E-value than in a whole-database search (bitscores are identical) and the 1e-5 cutoff is more permissive
- Because only the taxon is searched, its hits are no longer crowded out by the 500-subject limit per query
- Giving `--taxid` avoids the preliminary whole-database search, which takes as long as an unrestricted search

**Required Arguments:**
- `--input` - Classified ORFs FASTA (aa-sequences) (output from classify_ORFs)
- `--taxdump` - TaxonKit taxonomy dump directory (`nodes.dmp`, `names.dmp`, ...) matching the taxids of the BLAST database
- `--out-dir` - Output directory (default: ./out)

**Optional Arguments:**
- `--blastdb` - BLAST database (path prefix; default: `db/blast/GTDBtax_r232` in `--db-dir`)
- `--n-threads` - Number of threads to be used by tBLASTn (default: 5)
- `--level` - Restrict the search to the reference's taxon at this level (`kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`; default: whole database)
- `--taxid` - Reference taxid for the search restriction and the taxonomic context (default: inferred from the top hit). Must exist in the taxonomy dump used by TaxonKit.

**Output:**
- `homologs/tBLASTn.tsv` - tBLASTn hits (tab-separated, outfmt 6) with the columns
  `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids` and a 14th column `tax_context`

Intermediate files (raw tBLASTn output `tBLASTn_raw.tsv`, the taxid list `search_taxids.txt` and the preliminary search `tBLASTn_reference.tsv`) are removed after the run.
---

### 3. `align_codons`

**Purpose:** Create high-quality codon alignments from BLAST homolog results for comparative genomics analysis.

**Pipeline Steps:**
1. **Filter BLAST results** - Keeps top N homologs per ORF (based on E-value)
2. **Extract sequences** - Uses blastdbcmd to retrieve homologous sequences
3. **Clean sequences** - Removes duplicates and sequences with internal stop codons
4. **Align proteins** - Uses MUSCLE to align amino acid sequences
5. **Generate codon alignments** - Uses PAL2NAL to convert protein alignment to codon alignment
6. **Format output** - Removes special characters from headers for downstream tools

**Required Arguments:**
- `--input` - BLAST results file (output from find_homologs)
- `--out-dir` - Output directory (default: ./out)

**Optional Arguments:**
- `--blastdb` - BLAST database (path prefix; default: `db/blast/GTDBtax_r232` in `--db-dir`)
- `--n-homologs` - Maximum homologs to keep per ORF (default: 100)

**Dependencies:**
- MUSCLE (multiple sequence alignment)
- PAL2NAL (codon alignment from protein alignment)
- seqkit (sequence processing)

**Output:**
- `align/ID.aln` files - Codon alignments in CLUSTAL format
- `align/ID.info` files - Metadata including classification and frameshift information

---

### 4. `run_RNAcode`

**Purpose:** Analyze codon alignments to detect coding potential and protein-coding signals.

**How it works:**
- Runs RNAcode on each codon alignment file
- Detects highly conserved codon usage patterns indicative of protein-coding sequences
- Can process single files or batch process entire directories

**Required Arguments:**
- `--input` - Codon alignment .aln file or directory containing .aln files
- `--out-dir` - Output directory (default: ./out)

**Dependencies:**
- RNAcode (conservation and codon usage analysis tool)

**Output:**
- `RNAcode/ID.txt` files - RNAcode analysis results with coding potential scores
---

### 5. `run_OLGenie`

**Purpose:** Detect frameshift signals in overlapping genes to identify non-standard genetic code usage.

**How it works:**
- Reads frameshift information from .info metadata files
- Skips ORFs in same frame as reference CDS (frameshift: "=")
- Intergenic ORFs are run in all frame modes (ss12, ss13, sas11, sas12, sas13)
- Converts codon alignments to FASTA format
- Runs OLGenie with appropriate frameshift parameters
- Detects overlapping genes that may use alternative translation frames

**Frameshift Modes:**
- `ss12` - Same strand, frame shift 1→2
- `ss13` - Same strand, frame shift 1→3
- `sas11` - Opposite Strand, frame shift 1→1
- `sas12` - Opposite Strand, frame shift 1→2
- `sas13` - Opposite Strand, frame shift 1→3
- `=` - Same strand, No frameshift (skipped)

**Required Arguments:**
- `--input` - Codon alignment file or directory containing .aln files
- `--out-dir` - Output directory (default: ./out)

**Dependencies:**
- OLGenie (overlapping gene analysis tool, [chasewnelson/OLGenie](https://github.com/chasewnelson/OLGenie), Chase W. Nelson). The scripts are shipped with OverlapTools in `src/OLGenie/`:
  - `OLGenie.pl` - OLGenie itself (Perl, only core modules)
  - `OLGenie_bootstrap.R` - bootstrapping of the OLGenie results (R packages `readr`, `stringr`, `boot`; installed via `env.yml`)

**Output:**
- `OLGenie/ID.tsv` files - OLGenie results with overlapping gene annotations
---

### `run_all`

**Purpose:** Run the whole pipeline in one call: `classify_ORFs` → `score_ORFs` (only if `--score` is given) → `find_homologs` → `align_codons` → `run_RNAcode` → `run_OLGenie`, each step on the output of the previous one, all in `--out-dir`.

**Required Arguments:**
- `--input` - ORF predictions (as for classify_ORFs)
- `--mode` - Input type (`REP`, `DR` or `CUSTOM`)
- `--gtf-path` - Genome annotation GTF
- `--taxdump` - TaxonKit taxonomy dump matching the BLAST database
- `--genome` - Genome FASTA (only with `--score prodigal` and/or `balrog`)

**Optional Arguments:** all options of the single steps: `--n-orfs`, `--score`, `--db-dir`, `--blastdb`, `--level`, `--taxid`, `--n-threads`, `--n-homologs`, `--out-dir`, `--log-level`

**Output:** the outputs of all steps (see Output Directory Structure)

---

### `build_blastdb`

**Purpose:** Build the BLAST nucleotide database with taxids used by `find_homologs` and `align_codons` (setup, see [SETUP.md](SETUP.md#4-blast-database---blastdb)).

**How it works:** gives every sequence of each genome FASTA the taxid of its assembly (GCA_/GCF_ accession in the file name, looked up in the taxonomy dump's `taxid.map`), combines all genomes and runs `makeblastdb -parse_seqids -taxid_map`.

**Required Arguments:**
- `--input` - Directory with genome FASTA files (`.fna`/`.fna.gz`, searched recursively), e.g. the unpacked GTDB representative genomes
- `--taxdump` - Taxonomy dump with `taxid.map` (assembly accession → taxid)

**Optional Arguments:**
- `--blastdb` - Path prefix of the database to build (default: `db/blast/GTDBtax_r232` in `--db-dir`)
- `--keep-intermediate` - Keep the combined FASTA (`<blastdb>.combined.fna`) and the sequence → taxid map (`<blastdb>.taxid.map`)

---

## Usage Examples

### Run the whole pipeline
```bash
python3 OvlpTls.py run_all \
  --input predictions.csv --mode DR --gtf-path genome.gtf --n-orfs 500 \
  --genome genome.fna --score prodigal balrog pgap bakta \
  --taxdump /path/to/gtdb-taxdump \
  --level family --taxid 533388576 \
  --n-threads 20 --out-dir results
```

### Basic workflow - Classify ORFs
```bash
python3 OvlpTls.py classify_ORFs \
  --input orfs.fasta \
  --gtf-path genome.gtf \
  --mode CUSTOM \
  --n-orfs 50 \
  --out-dir results
```

### Score ORFs
```bash
python3 OvlpTls.py score_ORFs \
  --input results/classified_orfs/classified_ORFs.fasta \
  --genome genome.fna \
  --score prodigal balrog pgap bakta \
  --n-threads 20 \
  --out-dir results
```

### Find homologs
```bash
python3 OvlpTls.py find_homologs \
  --input results/classified_orfs/classified_ORFs.fasta \
  --taxdump /path/to/gtdb-taxdump \
  --n-threads 20 \
  --out-dir results
```

### Find homologs, keep only hits from the same genus as a given taxid
```bash
python3 OvlpTls.py find_homologs \
  --input results/classified_orfs/classified_ORFs.fasta \
  --taxdump /path/to/gtdb-taxdump \
  --level genus \
  --taxid 562 \
  --out-dir results \
  --log-level DEBUG
```

### Create codon alignments
```bash
python3 OvlpTls.py align_codons \
  --input results/homologs/tBLASTn.tsv \
  --n-homologs 50 \
  --out-dir results
```

### Run RNAcode analysis
```bash
python3 OvlpTls.py run_RNAcode \
  --input results/align \
  --out-dir results
```

### Run OLGenie analysis
```bash
python3 OvlpTls.py run_OLGenie \
  --input results/align \
  --out-dir results
```

---

## Pipeline Architecture

### Module Organization

- **OvlpTls.py** - Main entry point with argument parsing, logging setup and workflow orchestration
- **src/classify.py** - ORF classification logic against GTF annotations
- **src/scores/** - ORF scoring (score_ORFs): `scoring.py` runs the selected methods, `orfs.py` locates the ORFs on the genome, `prodigal.py`, `balrog.py`, `pgap.py`, `bakta.py` implement the scores
- **src/homologs.py** - BLAST database searching, homolog detection and taxonomic filtering (TaxonKit)
- **src/align.py** - Sequence alignment and codon-level processing
- **src/blastdb.py** - Building the BLAST database with taxids (build_blastdb)
- **src/tools.py** - interface for RNAcode and OLGenie analyses
- **src/OLGenie/** - bundled third-party OLGenie scripts (`OLGenie.pl`, `OLGenie_bootstrap.R`)

### Data Flow

```
Input ORFs (.FASTA)
    ↓
[classify_ORFs] → Classified ORFs (.FASTA)
    ↓
[find_homologs] → BLAST Results (.tsv, with taxonomic context)
    ↓
[align_codons] → Codon Alignments (.aln)
    ├→ [run_RNAcode] → Coding Potential Results
    └→ [run_OLGenie] → Frameshift Analysis Results
```

---

## Output Directory Structure

```
out/
├── classified_orfs/
│   └── classified_ORFs.fasta
├── scores/
│   └── ORF_scores.tsv
├── homologs/
│   └── tBLASTn.tsv
├── align/
│   ├── [ORF_ID].aln
│   └── [ORF_ID].info
├── RNAcode/
│   └── [ORF_ID].txt
└── OLGenie/
    └── [ORF_ID]/
        ├── OLGenie[_frame].tsv
        └── OLGenie[_frame]_bootstrap.tsv
```

---

## Requirements

All requirements can be installed with `env.yml` (see [Installation](#installation)).

### External Tools
- **BLAST+** - tBLASTn and blastdbcmd for sequence searching and retrieval
- **TaxonKit** - Taxonomy lookup for the taxonomic context/filtering in find_homologs
- **MUSCLE (v5)** - Multiple sequence alignment
- **PAL2NAL** - Codon alignment from protein alignment
- **seqkit** - Sequence processing utilities
- **RNAcode** - Coding region detection (optional)
- **Perl** - runs the bundled `src/OLGenie/OLGenie.pl` (optional)
- **R** with `readr`, `stringr`, `boot` - runs the bundled `src/OLGenie/OLGenie_bootstrap.R` (optional)

### Python Libraries
- Standard library (argparse, csv, logging, os, subprocess, ...)
- **pyrodigal**, **pytorch**, **pyhmmer**, **bakta** - ORF scoring in score_ORFs

---

## Error Handling

The tool includes comprehensive input validation:
- Checks file existence before processing
- Validates directory paths for BLAST databases
- Enforces required arguments per analysis type
- Provides clear error messages for missing inputs
- Failing external tools stop the run with an `ERROR` log message and exit code 1

---

## Notes

- All analyses assume aminoacid sequences as input
- BLAST database must be in nucleotide format
- GTF files should follow standard GTF/GFF3 format
- Output files are overwritten if they already exist
- Processing time scales quickly with --n-orfs, --n-homologs, --n-threads respectively

---

## License

OverlapTools is free software, licensed under the [GNU General Public License v3.0](LICENSE) (GPL-3.0-or-later).

## Third-party code and data

- **OLGenie** (`src/OLGenie/`): `OLGenie.pl` and `OLGenie_bootstrap.R` by Chase W. Nelson, [chasewnelson/OLGenie](https://github.com/chasewnelson/OLGenie), GPL-3.0. Included unchanged.
- **Prodigal / pyrodigal**: `src/scores/prodigal.py` re-implements Prodigal's node scoring following [pyrodigal](https://github.com/althonos/pyrodigal) (Martin Larralde) and [Prodigal](https://github.com/hyattpd/Prodigal) (Doug Hyatt), both GPL-3.0.
- **Balrog**: `src/scores/balrog.py` reproduces the ORF score of [BalrogCPP](https://github.com/salzberg-lab/BalrogCPP) (Markus J. Sommer & Steven L. Salzberg, MIT) using its pretrained models from [balrog_models](https://github.com/Markusjsommer/balrog_models) (MIT), which are downloaded on first use.
- **PGAP HMMs**: the `pgap` score downloads NCBI's protein family HMM library ([hmm_PGAP](https://ftp.ncbi.nlm.nih.gov/hmm/current/)) on first use.
- **Bakta**: the `bakta` score uses [Bakta](https://github.com/oschwengers/bakta) (GPL-3.0) and downloads its light database on first use.

Models and databases are stored in `db/` (or `--db-dir`) and are not part of the repository.

## Citation

If you use OverlapTools, please cite it via the [`CITATION.cff`](CITATION.cff) file (GitHub: "Cite this repository") and the tools it builds on: Prodigal, Balrog, Bakta, PGAP, BLAST+, TaxonKit, MUSCLE, PAL2NAL, RNAcode and OLGenie.
