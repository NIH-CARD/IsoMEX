# IsoMEX

IsoMEX is a stand-alone Python tool that converts IsoSeq/pigeon make-seurat output files into the [standardized Matrix Exchange (MEX)](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices) format, matching the [filtered feature-barcode matrix](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) produced by the 10x Genomics CellRanger pipeline. It fixes the duplicate gene entry issue found in the existing IsoSeq workflow (e.g. [isoseq.how](https://isoseq.how/)), ensuring consistent gene and transcript annotations across multiple samples. This allows for seamless integration with downstream single cell/single nuclei tools, including multi-sample UMAP clustering.

## Overview

The isomex.py script:
- Loads two tab-delimited files, [outputs from `pigeon make-seurat`](https://isoseq.how/classification/pigeon-output.html), with a common basename (Note that while they have the ".csv" extension, they are actually tab-delimited.):
  - `<basename>.info.csv`
  - `<basename>.annotated.info.csv`
- Merges these files on the `id` column.
- Optionally filters rows based on [SQANTI3 isoform categories](https://isoseq.how/classification/categories.html) provided via command-line arguments.
- Aggregates counts by grouping on the `gene` or `transcript` column and by cell barcode (`BC`).
- Sums the counts and converts the result into a MEX-format directory (containing `matrix.mtx`, `features.tsv`, and `barcodes.tsv`).
- Gzips the output files (producing files such as `matrix.mtx.gz`).


It produces two sets of output directories with the filtered feature-barcode matrix in the 10× Genomics cellranger MEX format (i.e. a directory containing `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv,gz`). One output is a gene-level summary and the other a transcript-level summary.

## Features

- **Fixes Duplicate Gene Entries:**  
  Ensures consistent gene and transcript annotations across samples, preventing issues caused by duplicate gene entries in the original Iso-Seq workflow.
  
- **10× Genomics Compatible:**  
  Generates output files (`matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`) that conform to CellRanger's MEX format, making outputs comparable to the standard 10x Genomics pipeline.

- **Filter by Isoform Classification Category:**  
  Allows filtering based on SQANTI3 isoform classification categories, such as `full-splice_match`, `novel_in_catalog`, or `antisense`. See the full list of categories: [Iso-Seq Classification](https://isoseq.how/classification/categories.html).

- **Consistent Gene and Transcript Feature Annotation:**
  Uses gene and transcript mapping files to ensure feature names and IDs are preserved while still allowing novel genes and transcripts to be included in the output.
  
- **Compressed Output:**  
  Outputs are gzipped (`.gz`) for efficient storage, and using MEX files instead of CSVs significantly reduces file size.
  
- **Cluster-Ready:**  
  Includes example SLURM array job scripts for running IsoMEX on HPC systems like Biowulf.

---

## **Installation**

### **Dependencies**
- Python 3.6+
- [pandas](https://pandas.pydata.org/)
- [gffutils](https://daler.github.io/gffutils/)
- Standard Python libraries: `argparse`, `os`, `gzip`, `shutil`

Install dependencies using:

```bash
pip install pandas gffutils
```

---

## **Usage**

The basic command-line usage is:

```bash
python isomex.py <basename> --gene_map <gene_map.txt> --transcript_map <transcript_map.txt> [--filter_category "cat1,cat2"] --output_dir <output_directory>
```

### **Example**
If you have:
- `sample1.info.csv` and `sample1.annotated.info.csv`
- A gene mapping file: `gene_map.txt`
- A transcript mapping file: `transcript_map.txt`

Run:

```bash
python isomex.py sample1 --gene_map gene_map.txt --transcript_map transcript_map.txt --output_dir filtered_feature_bc_matrix
```

This will generate two directories (`gene_filtered_feature_bc_matrix/` and `transcript_filtered_feature_bc_matrix/`) with MEX-formatted output:
- `matrix.mtx.gz`
- `features.tsv.gz`
- `barcodes.tsv.gz`

---

## **Generating Gene & Transcript Mapping Files**
IsoMEX requires **two separate mapping files** before running:
1. **Gene Map File** (`gene_map.txt`)  
   - A tab-delimited file containing `gene_id` and `gene_name`.
   
2. **Transcript Map File** (`transcript_map.txt`)  
   - A tab-delimited file containing `transcript_id` and `transcript_name`.

### **Generating Mapping Files from a GTF**
To create these files from an annotation GTF file, use the included utility script:

```bash
python utils/generate_map.py annotation.forPigeon.gtf gene_map.txt transcript_map.txt
```

This will extract:
- `gene_id` and `gene_name` → `gene_map.txt`
- `transcript_id` and `transcript_name` → `transcript_map.txt`

---

## **Running IsoMEX on an HPC Cluster (SLURM)**
IsoMEX includes a **SLURM job array script** for processing multiple samples efficiently.

### **Submitting the Job**
1. **Prepare a sample list** (`samples.txt`) with one sample per line:
   ```
   sample1
   sample2
   sample3
   ```

2. **Submit the SLURM job array**:
   ```bash
   sbatch slurm/isomex.sh
   ```
