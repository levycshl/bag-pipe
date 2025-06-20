<p align="center">
  <img src="docs/logo_bag-pipe.png" height="256" alt="bag-pipe logo"/>
</p>

# Bag-Pipe <!-- omit in toc -->

End-to-end pipeline that converts BAG paired-end FASTQ files → gene-resolved
read tables suitable for repertoire analysis and QC.

* **Stream-based** — constant memory, even on >100 GB BAMs  
* **Modular** — each step is a CLI sub-command (`trim`, `extract`, `summarise`)
  so you can slot them into any workflow manager  
* **Snakemake workflow included** for single-command reproducibility  
* **Pure Python** core; heavy lifting (HISAT2, samtools, GNU sort) stays outside  
* **MIT/BSD licensed** — free for academic & commercial use

| | build | tests | license |
|---|---|---|---|
| **Status** | ![CI badge](https://img.shields.io/badge/CI-passing-brightgreen?logo=github) | ![pytest](https://img.shields.io/badge/tests-100%25-brightgreen) | ![MIT](https://img.shields.io/badge/license-MIT-blue) |

---

## Table of contents <!-- omit in toc -->
1. [Quick install](#quick-install)
2. [CLI overview](#cli-overview)
3. [Minimal example](#minimal-example)
4. [Full Snakemake pipeline](#full-snakemake-pipeline)
5. [Configuration](#configuration)
6. [Output files](#output-files)
7. [Development & testing](#development--testing)
8. [Citation](#citation)
9. [License](#license)

---

## Quick install

```bash
# 1 – Clone the public repo
git clone https://github.com/YourOrg/bag-pipe.git
cd bag-pipe

# 2 – Create the exact software stack declared in environment.yml
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate bag-pipe                  # env name is set inside the file

# 3 – Install *this* checkout into the env
#     └─ pip reads setup.py, builds a wheel, and
#        pulls read-break from GitHub automatically
pip install .

# 4 – Smoke-test the CLIs (should print help text)
bag_pipe trim --help
bag_pipe extract --help
bag_pipe summarise --help

````

Everything is also pinned in **`environment.yml`**;
`mamba env create -f environment.yml` reproduces the exact toolchain.

---

## CLI overview

| Command              | Purpose                                                                   |
| -------------------- | ------------------------------------------------------------------------- |
| `bag_pipe trim`      | Clip + barcode-tag FASTQ pairs via **read-break**.                        |
| `bag_pipe extract`   | Stream a (sorted) BAM/SAM → 13-column TSV of per-alignment features.      |
| `bag_pipe summarise` | Roll up varietal-tags, map to gene models, emit VT summaries per barcode. |

Run any sub-command with `--help` for full argument reference.

---

## Minimal example

```bash
# 1. Clip & tag raw reads
bag_pipe trim \
  --r1 fastq/R1.fq.gz \
  --r2 fastq/R2.fq.gz \
  --parse-config configs/bagpipe_2N.yaml \
  --sample-id demo \
  --outdir trimmed/

# 2. Align + sort (HISAT2 → BAM)
hisat2 -x /refs/hg19 -1 trimmed/demo.R1.fastq.gz -2 trimmed/demo.R2.fastq.gz \
  | samtools sort -o alignment/demo.sorted.bam -
samtools index alignment/demo.sorted.bam

# 3. Feature extraction for chr14
samtools view -h alignment/demo.sorted.bam chr14 \
  | bag_pipe extract /dev/stdin tagmap/demo.chr14.tsv.gz

# 4. VT summarisation
bag_pipe summarise \
  --chrom chr14 \
  --tsv-in  tagmap/demo.chr14.tsv.gz \
  --ref-flat refs/hg19.refFlat.txt \
  --tsv-out results/demo.chr14.vt.tsv.gz
```

---

## Full Snakemake pipeline

```bash
snakemake -j 8 --use-conda          # default DAG (all chromosomes, all samples)
```

* **`config/pipeline_params.yaml`** — paths to FASTQ, splice sites, indexes
* **`workflow/Snakefile`** — modular rules: `trim_and_tag`, `hisat2_align`,
  `sort_bam`, `extract_features_chr`, `sort_features_chr`, `summarise_chr`
* Final artefacts land in `results/{sample}.{chr}.vt.tsv.gz`.

**TODO** It would be nice if we had a `workflow/README.md` that offers cluster profiles and per-rule resource hints.

---

## Configuration

| Key in `pipeline_params.yaml` | Example                                         | Meaning                                                      |
| ----------------------------- | ----------------------------------------------- | ------------------------------------------------------------ |
| `hisat2_index`                | `/data/refs/hisat2/ucsc.hg19`                   | basename of HISAT2 index                                     |
| `splice_sites`                | `/data/refs/hisat2/hg19.refseq.splicesites.txt` | known junctions                                              |
| `chromosomes`                 | `['chr14', 'chr21', 'chr22']`                   | chromosomes to process                                       |
| `samples`                     | YAML list                                       | each item needs `id`, `fastq_r1`, `fastq_r2`, `parse_config` |

---

## Output files

| Path pattern                       | Description                                  |
| ---------------------------------- | -------------------------------------------- |
| `trimmed/{sample}.R1.fastq.gz`     | clipped + tagged reads                       |
| `alignment/{sample}.sorted.bam`    | coordinate-sorted BAM                        |
| `tagmap/{sample}.{chr}.sorted.tsv` | per-alignment feature table                  |
| `results/{sample}.{chr}.vt.tsv.gz` | VT-level summary (one line per surviving VT) |

**TODO** Column definitions will live in **`docs/output_schema.md`**.

---

## Development & testing

```bash
# lint + unit tests + tiny integration test
pip install -e .[dev]
ruff check .
pytest -q
```

CI runs the same checks plus a miniature Snakemake DAG on ≤50 kB toy data.

---

## Citation

> 

Please cite the companion tool **read-break** separately if you use the trimming module.

---

## License

Distributed under the **MIT License** (see `LICENSE`).
