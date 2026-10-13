# Setup

OverlapTools needs a conda environment plus some external data that is not part of the repository:

| Data | Needed for | Size | How |
|---|---|---|---|
| Scoring models and databases (`db/`) | `score_ORFs` | ~7 GB | downloaded automatically on first use, or beforehand (step 2) |
| GTDB taxonomy dump | `find_homologs` (`--taxdump`) | ~150 MB unpacked | download (step 3) |
| BLAST nucleotide database with taxids | `find_homologs`, `align_codons` (`--blastdb`) | ~160 GB for GTDB R232 | build with `build_blastdb` (step 4) |

## 1. Environment

```bash
git clone https://github.com/Kani0n/OverlapTools.git
cd OverlapTools
conda env create -f env.yml     # or: mamba env create -f env.yml
conda activate OverlapTools
python3 OvlpTls.py --help
```

All commands below are run in the `OverlapTools` folder with the environment activated.

## 2. Scoring models and databases (`score_ORFs`)

`score_ORFs` downloads what a scoring method needs the first time that method is used, into `db/` in the OverlapTools folder (or the directory given with `--db-dir`). Only the methods you select with `--score` are downloaded. `prodigal` needs no download.

To download everything beforehand (e.g. on a machine with internet access before running on a cluster node without), run:

```bash
# balrog: gene and TIS models (< 1 MB)
python3 -c "
import torch
torch.hub.set_dir('db/balrog')
for model in ('geneTCN', 'tisTCN'):
    torch.hub.load('Markusjsommer/balrog_models', model, trust_repo=True, verbose=False)
"

# pgap: NCBI protein family HMM library (~3 GB)
mkdir -p db/pgap
curl -L -o db/pgap/hmm_PGAP.LIB https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB

# bakta: Bakta light database (~1.5 GB download, ~4 GB unpacked)
bakta_db download --output db/bakta --type light
```

The resulting layout:

```
db/
├── balrog/
│   ├── Markusjsommer_balrog_models_main/
│   └── trusted_list
├── pgap/
│   └── hmm_PGAP.LIB
├── bakta/
│   └── db-light/
└── blast/                 # BLAST database, built in step 4
    └── GTDBtax_r232.*
```

To keep the databases elsewhere, download into that directory instead of `db/` and pass it with `--db-dir`.

## 3. Taxonomy dump (`find_homologs --taxdump`)

`find_homologs` resolves taxids with TaxonKit, so it needs the taxonomy dump whose taxids the BLAST database uses. For GTDB, download the release matching your database from [shenwei356/gtdb-taxdump](https://github.com/shenwei356/gtdb-taxdump/releases), e.g. GTDB R232:

```bash
curl -LO https://github.com/shenwei356/gtdb-taxdump/releases/download/v0.7.0/gtdb-taxdump-R232.tar.gz
tar -xzf gtdb-taxdump-R232.tar.gz      # -> gtdb-taxdump-R232/ (nodes.dmp, names.dmp, taxid.map, ...)
```

Pass the unpacked folder with `--taxdump gtdb-taxdump-R232`. `taxid.map` maps genome assemblies (GCA_/GCF_ accessions) to their taxids, which is also how to find the taxid of your own genome for `--taxid`:

```bash
grep GCF_039605695.1 gtdb-taxdump-R232/taxid.map
echo "Salmonella enterica" | taxonkit name2taxid --data-dir gtdb-taxdump-R232 --show-rank
```

## 4. BLAST database (`--blastdb`)

`find_homologs` and `align_codons` search a BLAST nucleotide database whose sequences carry taxids from the taxonomy dump of step 3. OverlapTools builds it from the GTDB representative genomes (one genome per GTDB species) with `build_blastdb`.

**Disk space** (GTDB R232): ~180 GB download, ~180 GB unpacked genomes, ~670 GB temporary combined FASTA during the build, ~160 GB final database. The download and the unpacked genomes can be deleted after the build.

### 4.1 Download the genomes

Use the GTDB release that matches the taxonomy dump of step 3:

```bash
curl -LO https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/gtdb_genomes_reps_r232.tar.gz
tar -xzf gtdb_genomes_reps_r232.tar.gz   # -> gtdb_genomes_reps_r232/database/GC[AF]/.../GCF_000006945.2_genomic.fna.gz
```

### 4.2 Build the database

```bash
python3 OvlpTls.py build_blastdb \
  --input gtdb_genomes_reps_r232 \
  --taxdump gtdb-taxdump-R232
```

The database is built at the default location `db/blast/GTDBtax_r232` (in `--db-dir`), which `find_homologs`, `align_codons` and `run_all` use when `--blastdb` is not given. To build it elsewhere, e.g. on a disk with more space, pass the path prefix with `--blastdb` and give the same `--blastdb` to the later steps.

What `build_blastdb` does:

1. Reads the assembly → taxid map of the taxonomy dump (`gtdb-taxdump-R232/taxid.map`, e.g. `GCF_000006945.2  1002421060`)
2. Finds all `.fna`/`.fna.gz` files below `--input`, takes the assembly accession (GCA_/GCF_) from each file name and gives every sequence of the file the assembly's taxid. The sequences have their own accessions (e.g. `NC_003197.2`), so this per-sequence map is what makeblastdb needs.
3. Writes all genomes into one combined FASTA (`<blastdb>.combined.fna`) and the map (`<blastdb>.taxid.map`). Files without accession or without taxid and duplicate sequence accessions are skipped with a warning.
4. Runs `makeblastdb -dbtype nucl -parse_seqids -taxid_map <blastdb>.taxid.map -blastdb_version 5`
5. Removes the combined FASTA and the map (keep them with `--keep-intermediate`)

Check the result:

```bash
blastdbcmd -db db/blast/GTDBtax_r232 -info               # sequences and bases
blastdbcmd -db db/blast/GTDBtax_r232 -tax_info | head    # number of taxids (one per genome)
```

### Using another database

Any BLAST nucleotide database works if it is in BLAST v5 format with taxids (`makeblastdb -parse_seqids -taxid_map ...`) and its taxids come from the taxonomy dump passed with `--taxdump`.

## 5. Check the setup

```bash
# scoring (prodigal needs only the genome)
python3 OvlpTls.py score_ORFs --input classified_ORFs.fasta --genome genome.fna --score prodigal balrog pgap bakta --out-dir test

# homolog search
python3 OvlpTls.py find_homologs --input classified_ORFs.fasta --taxdump gtdb-taxdump-R232 --level genus --taxid <taxid> --out-dir test
```
