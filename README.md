## nanoID: Recovery of Exact Sequence Variants from Long-Read Amplicon Data

### What is nanoID?
`nanoID` is a bioinformatics pipeline for recovering amplicon sequence variants (ASVs) from noisy long-read amplicon sequencing data. It supports Oxford Nanopore Technologies (ONT) R10.4.1 and Pacific Biosciences (PacBio) sequencing data and provides workflows for both per-sample ASV recovery (`nanoid condens`) and multi-sample ASV profiling (`nanoid profile`).

### How does nanoID work?

#### 🔬 Per-sample ASV recovery: `nanoID condens`

```
Reads
  │
  ├─ Optional: primer trimming, quality/length filtering, orientation correction
  ▼
Split reads into N disjoint subsets
  │
  ├─ Near-neighbor search (VSEARCH)
  ├─ Consensus generation (abPOA)
  └─ Consensus sequences ("conseqs")
  ▼
Cross-split shared-neighbor graph
  • Nodes: unique conseqs
  • Edges: shared contributing reads
  ▼
Graph-based denoising
(constrained abundance ascent, retain nodes without allowable ascent)
  ▼
Candidate ASVs
  ▼
Abundance estimation
(read matching, Levenshtein distance, + EM)
  ▼
Chimera filtering
(UCHIME3 + custom filter)
  ▼
Final ASVs and abundances
```

1. **Near-neighbor search**  
   For each read, `nanoID condens` identifies closely related reads using pairwise sequence identity. By default, 4 neighbors are used for consensus generation.

2. **Consensus generation**  
   Consensus sequences are generated from, potentially, overlapping read neighborhoods using adaptive banded partial-order alignment. Because neighborhoods can overlap, a read may contribute to multiple consensus sequences.

3. **Graph-based denoising**  
Consensus sequences are counted and represented as nodes in a shared-neighbor graph, with edges connecting consensus sequences that share contributing reads. Candidate ASVs are identified by constrained abundance ascent, in which less abundant conseqs are linked to sufficiently more abundant neighboring conseqs.

To reduce false-positive ASVs, reads are processed in multiple disjoint splits. Consensus sequences and shared-neighbor graphs are generated independently for each split and consolidated.

4. **Abundance estimation**  
   Reads are matched to candidate ASVs using Levenshtein distance, and their proportional assignments are estimated by expectation-maximization (EM).

5. **Chimera filtering**  
  Candidate ASVs are screened for chimeras using UCHIME and a more stringent prefix-suffix matching algorithm, yielding a final set of non-chimeric ASVs and their estimated abundances.


### 📊 Multi-sample integration: `nanoID profile`

```
ASVs from all samples
      +
Pre-denoising conseqs
  ▼
Global ASV catalogue
  ▼
Sample-wise ASV rescue
(minimum abundance threshold)
  ▼
Expanded ASV sets
  ▼
EM quantification
  ▼
ASV abundance matrix

```
`nanoID profile` improves consistency of ASV detection across samples by rescuing ASVs that are present in the global ASV set and supported by pre-denoising consensus sequences within a sample. The expanded sample-specific ASV sets are then re-quantified. 

Optionally, ASVs can be clustered into high-identity operational taxonomic units (OTUs), leveraging the high accuracy of ASV sequences rather than clustering noisy reads directly.

## 🚀 Installation

#### Requirements
- nanoID requires Python 3.10–3.13 for compatibility with pyabpoa
- nanoID requires the following external programs:
    - VSEARCH ([https://github.com/torognes/vsearch])
    - Cutadapt, optional ([https://github.com/marcelm/cutadapt])
    - Emu, optional ([https://github.com/treangenlab/emu])

#### Installation

```
conda create -n nanoid \
    -c conda-forge -c bioconda \
    "python>=3.10,<3.14" vsearch cutadapt
conda activate nanoid
wget https://github.com/dietertourlousse/nanoID/releases/download/v0.1/nanoid-0.1.tar.gz
pip install nanoid-0.1.tar.gz
```

##### Verify installation

After installation, confirm that `nanoid` is available:

```sh
nanoid -h
```

```
usage: nanoid [-h] [-v] <command> ...

nanoID: Long‑read amplicon denoising and profiling tool

options:
  -h, --help     show this help message and exit
  -v, --version  Show program version and exit.

Subcommands:
  <command>
    condens      Workflow for recovering ASVs.
    profile      Workflow for profiling based on recovered ASVs.
```

---

## ⚡ Quick start

🖥️ **Basic usage**

Run `nanoid condens` on a single sample:
```sh
nanoid condens -i input.fastq -o nanoid_condens
```

For downstream joint profiling across multiple samples, write all `nanoid` outputs to a single directory.

```sh
mkdir nanoid

## Using a shell loop
for fq in *.fastq; do
    nanoid condens -i "$fq" -o nanoid_condens
done

## Using GNU Parallel
parallel -j 4 \
    "nanoid condens -i {} -o nanoid_condens" \
    ::: *.fastq
```

Once all samples have been processed, run `nanoid profile` on the combined output directory.
```sh
nanoid profile -i nanoid_condens -o nanoid_profile
```

---

🖥️ **Detailed usage**

`nanoid condens -h`

```
General settings:
  --input_fastq INPUT_FASTQ, -i INPUT_FASTQ
                        Input FASTQ file. (default: None)
  --output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
                        Output directory. (default: None)
  --randseed RANDSEED, -s RANDSEED
                        Random seed used for fastq splitting / subsampling. (default: 21336)
  --threads THREADS, -t THREADS
                        Number of threads. (default: 12)
  --remove_intermediates
                        Add this flag to remove intermediate files. (default: False)

Cutadapt parameters:
  --cutadapt_f_primer CUTADAPT_F_PRIMER
                        Forward primer. (default: AGRGTTYGATYHTGGCTCAG)
  --cutadapt_r_primer_rc CUTADAPT_R_PRIMER_RC
                        Reverse primer (reverse-complemented). (default: AAGTCGTAACAAGGTARCCG)
  --cutadapt_f_primer_overlap CUTADAPT_F_PRIMER_OVERLAP
                        Min overlap for forward primer (None: use primer length). (default: None)
  --cutadapt_r_primer_overlap CUTADAPT_R_PRIMER_OVERLAP
                        Min overlap for reverse primer (None: use primer length). (default: None)
  --cutadapt_f_primer_errors CUTADAPT_F_PRIMER_ERRORS
                        Max mismatches for forward primer. (default: 2)
  --cutadapt_r_primer_errors CUTADAPT_R_PRIMER_ERRORS
                        Max mismatches for reverse primer. (default: 2)
  --cutadapt_mean_quality CUTADAPT_MEAN_QUALITY
                        Min mean read quality after trimming. (default: 16)
  --cutadapt_minimum_length CUTADAPT_MINIMUM_LENGTH
                        Min read length after primer trimming. (default: 1200)
  --cutadapt_maximum_length CUTADAPT_MAXIMUM_LENGTH
                        Max read length after primer trimming. (default: 1800)
  --cutadapt_maximum_expected_errors CUTADAPT_MAXIMUM_EXPECTED_ERRORS
                        Max expected errors after primer trimming. (default: 99)
  --skip_cutadapt       Add this flag to skip primer trimming with Cutadapt. (default: False)

Splitting parameters:
  --fastq_splits FASTQ_SPLITS
                        Number of FASTQ splits. (default: 2)
  --reads_per_split READS_PER_SPLIT
                        Number of reads per FASTQ split (0: even splits). (default: 0)

Condens parameters:
  --neighbors_n NEIGHBORS_N
                        Number of neighbors per consensus. (default: 4)
  --kappa KAPPA         Abundance dominance threshold for graph acsent. (default: 4.0)
  --min_conseq_size MIN_CONSEQ_SIZE
                        Min size (counts) conseq, per split. (default: 8)

Additional chimera/bimera detection parameters:
  --skip_additional_chimera_check
                        Add this flag to SKIP additional bimera/chimera filtering (CAUTION: aggressive default settings, may lead to false-positives).
                        (default: False)
  --abskew ABSKEW       Min abundance ratio parent to query. (default: 2.0)
  --min_parent_len_frac MIN_PARENT_LEN_FRAC
                        Min prefix/suffix lenght fraction of parent to match query. (default: 0.05)
  --min_parent_div MIN_PARENT_DIV
                        Min divergence between parents. (default: 0.03)
  --min_query_cov MIN_QUERY_COV
                        Min fraction of query covered by parent prefix + suffix. (default: 0.5)
  --allow_one_off ALLOW_ONE_OFF
                        Add this flag to allow one error between query and parent for prefix and suffix. (default: False)
```

`nanoid profile -h`

```
  -i INPUT_DIRECTORY, --input_directory INPUT_DIRECTORY
                        Input directory containing FASTQ files and FASTA files. (default: None)
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Output directory (default: None)
  --fastq_pattern FASTQ_PATTERN
                        Substring pattern for FASTQ files from input directory. (default: _split.all.fastq)
  --fasta_pattern FASTA_PATTERN
                        Substring pattern for FASTA files from input directory. (default: _nanoid_final.fasta)
  --min_size_uniques MIN_SIZE_UNIQUES
                        Min total size (counts) to retain ASV in global pool. (default: 2)
  --min_prev_uniques MIN_PREV_UNIQUES
                        Min prevalence to retain ASV in global pool. (default: 1)
  --recruit_min_size RECRUIT_MIN_SIZE
                        Min size (count, per split) for ASV to recruit/rescue to per-sample ASVs based on conseqs prior to denoising. (default: 2)
  --cluster_id CLUSTER_ID
                        Identity threshold for OTU clustering. (default: 0.99)
  --run_emu_quant       Add this flag to Emu quantification of OTUs. (default: False)
  --include_seqs_emu_db {all,representatives}
                        Sequences to include in Emu database. (default: representatives)
  -t THREADS, --threads THREADS
                        Number of threads. (default: 12)
  -f, --force           Overwrite output directory if it already exists. (default: False)
  -y, --yes             Assume yes for overwrite prompts. (default: False)
  -h, --help            Show this help message and exit.
```

---
##### Output

* Fasta file `{basename}_nanoid_final.fasta` of ASV sequences with USEARCH/VSEARCH-style size annotations.
* Log file `{basename}_nanoid_condens.log`

##### Output

* Fasta file `nanoid_asvs.fasta` of ASV sequences.
* Count table `nanoid_asvs_cts.tsv` ASV abundances.
* Log file `nanoid_profile.log`

---

## Third-party attributions

`nanoID` rely on several third‑party tools and we recommend citing the original publications of these tools.

* Rognes T, Flouri T, Nichols B, Quince C, Mahé F (2016).
  VSEARCH: a versatile open source tool for metagenomics.
  PeerJ 4:e2584.
  doi:[10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)

* Curry KD, Wang Q, Nute MG, Tyshaieva A, Reeves E, Soriano S, Wu Q, Graeber E, Finzer P, Mendling W, Savidge T, Villapol S, Dilthey A, Treangen TJ (2022).
  Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data.
  Nat Methods 19(7):845-853.
  doi:[10.1038/s41592-022-01520-4](https://doi.org/10.1038/s41592-022-01520-4)

* Marcel M (2011).
  Cutadapt removes adapter sequences from high-throughput sequencing reads.
  EMBnet Journal 17(1):10-12.
  doi:[10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)

---

## Citation

A manuscript describing nanoID is in preparation. In the meantime, please cite this repository when using nanoID.

---

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a complete list of changes across releases.

---

## Disclaimer

Portions of the codebase were developed with assistance from AI tools, specifically Claude AI (Sonnet 4.6) and Microsoft 365 Copilot (GPT‑5).

---

## License

Copyright (C) 2026 Dieter Tourlousse

This project is licensed under the GNU General Public License
version 3 or later (GPL-3.0-or-later).

See the LICENSE file for details.
