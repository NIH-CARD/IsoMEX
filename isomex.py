#!/usr/bin/env python3

"""
This script loads two text files that share the same basename:
  - <basename>.info.csv
  - <basename>.annotated.info.csv

It then produces two sets of output MEX directories (one for gene-level and one for transcript-level summaries)
following the 10x Genomics cellranger MEX format (matrix.mtx, features.tsv, barcodes.tsv).

https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices

Can filter on SQANTI3 isoform classification categories:
https://isoseq.how/classification/categories.html
array(['full-splice_match', 'novel_in_catalog', 'incomplete-splice_match',
       'novel_not_in_catalog', 'genic', 'fusion', 'intergenic',
       'antisense', 'moreJunctions'], dtype=object)

Usage:
    python script.py <basename> [--isoform_category CAT1,CAT2] [--output_dir OUTPUT]

Example:
    python script.py sample1 --isoform_category "full-splice_match,novel_in_catalog" --output_dir results
"""

import argparse
import os
import pandas as pd
import gzip
import shutil

def load_data(base_path):
    """
    Load the two CSV files (assumed to be tab-delimited) and merge on the 'id' column.
    """
    info_file = base_path + ".info.csv"
    annotated_file = base_path + ".annotated.info.csv"
    
    # Read the txt files - note that they are actually tab delimited, even though they have the ".csv" extension
    info_df = pd.read_csv(info_file, sep='\t')
    annotated_df = pd.read_csv(annotated_file, sep='\t')
    
    # Merge on the "id" column; if both files have overlapping column names (other than id),
    # the ones from the annotated file will be suffixed.
    merged_df = pd.merge(info_df, annotated_df, on="id", suffixes=("", "_annotated"))
    return merged_df

def load_gene_map(gene_map_file):
    """
    Load a tab-delimited gene mapping file.
    Expected columns: gene_id, gene_name.
    Returns a dictionary mapping gene_name to (gene_id, gene_name).
    """
    df = pd.read_csv(gene_map_file, sep='\t')
    gene_map = dict(zip(df['gene_name'], zip(df['gene_id'], df['gene_name'])))
    return gene_map

def load_transcript_map(transcript_map_file):
    """
    Load a tab-delimited transcript mapping file.
    Expected columns: transcript_id, transcript_name.
    Returns a dictionary mapping transcript_name to (transcript_id, transcript_name).
    """
    df = pd.read_csv(transcript_map_file, sep='\t')
    transcript_map = dict(zip(df['transcript_name'], zip(df['transcript_id'], df['transcript_name'])))
    return transcript_map

def filter_data(df, isoform_categories):
    """
    Optionally filter the merged DataFrame by the 'category' column.
    """
    if isoform_categories:
        df = df[df['category'].isin(isoform_categories)]
    return df

def gzip_file(file_path):
    output_file = file_path + ".gz"
    with open(file_path, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(file_path)

def create_mex_matrices(df, group_by, output_prefix, gene_map=None, transcript_map=None):
    """
    Create a MEX-format output given:
      - df: the merged pandas DataFrame
      - group_by: the column to aggregate by ("gene" or "transcript")
      - output_prefix: output directory where files will be saved
      - gene_level_map: dictionary mapping gene names to (gene_id, gene_name); used if group_by=="gene".
      - transcript_level_map: dictionary mapping transcript names to (transcript_id, transcript_name); used if group_by=="transcript".
    
    The function automatically sets:
      - feature_type_label = "Gene Expression" if group_by=="gene", or
      - feature_type_label = "Transcript Expression" if group_by=="transcript".
    
    The function groups by the feature and barcode ('BC'), sums the 'count', and writes three files:
      matrix.mtx, features.tsv, and barcodes.tsv.
    """

    # Set feature_type_label and mapping dictionary automatically.
    if group_by == "gene":
        feature_type_label = "Gene Expression"
        mapping = gene_map
        outdir = "gene_" + str(output_prefix)
    elif group_by == "transcript":
        feature_type_label = "Transcript Expression"
        mapping = transcript_map
        outdir = "transcript_" + str(output_prefix)
    else:
        feature_type_label = group_by + " Expression"
        mapping = None
        outdir = str(group_by) + "_" + str(output_prefix)

    # Group by the feature (gene or transcript) and the cell barcode (BC)
    grouped = df.groupby([group_by, 'BC'])['count'].sum().reset_index()
    
    # Get sorted lists of unique features and barcodes for reproducible output order
    features = sorted(grouped[group_by].unique())
    barcodes = sorted(grouped['BC'].unique())
    
    # Create mappings (we will write the matrix in 1-indexed format per MEX convention)
    feature_to_index = {feat: i for i, feat in enumerate(features)}
    barcode_to_index = {bc: i for i, bc in enumerate(barcodes)}
    
    # Build coordinate lists for the matrix (using 1-indexing)
    rows, cols, data = [], [], []
    for _, row in grouped.iterrows():
        feat = row[group_by]
        bc = row['BC']
        cnt = row['count']
        rows.append(feature_to_index[feat] + 1)
        cols.append(barcode_to_index[bc] + 1)
        data.append(cnt)
    
    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)
    mm_file = os.path.join(outdir, "matrix.mtx")
    features_file = os.path.join(outdir, "features.tsv")
    barcodes_file = os.path.join(outdir, "barcodes.tsv")
    
    # Write matrix.mtx using pandas for the numeric portion and manual header insertion.
    matrix_df = pd.DataFrame({'row': rows, 'col': cols, 'data': data})
    with open(mm_file, 'w') as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write("%\n")
        f.write(f"{len(features)} {len(barcodes)} {len(data)}\n")
    matrix_df.to_csv(mm_file, sep=' ', mode='a', header=False, index=False)
    
    # Write features.tsv: each line has feature_id, feature_name, and feature type label.
    features_list = []
    for feat in features:
        if mapping and (feat in mapping):
            feat_id, feat_name = mapping[feat]
        else:
            feat_id, feat_name = feat, feat
        features_list.append([feat_id, feat_name, feature_type_label])
    features_df = pd.DataFrame(features_list)
    features_df.to_csv(features_file, sep='\t', header=False, index=False)
    
    # Write barcodes.tsv: one barcode per line.
    barcodes_modified = [str(bc) + "-1" for bc in barcodes]
    barcodes_df = pd.DataFrame(barcodes_modified)
    barcodes_df.to_csv(barcodes_file, sep='\t', header=False, index=False)
    
    # Gzip the three files
    gzip_file(mm_file)
    gzip_file(features_file)
    gzip_file(barcodes_file)

    print(f"{group_by} output written to directory: {outdir}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate CellRanger comparable gene and transcript level MEX matrices from pigeon output files.")
    parser.add_argument("base", help="Base path for input files (e.g., 'sample1' for sample1.info.csv and sample1.annotated.info.csv)")
    parser.add_argument("--gene_map", type=str, default=None,
                        help="Path to the gene mapping file (tab-delimited with columns: gene_id, gene_name)")
    parser.add_argument("--transcript_map", type=str, default=None,
                        help="Path to the transcript mapping file (tab-delimited with columns: transcript_id, transcript_name)")
    parser.add_argument("--filter_category", type=str, default=None,
                        help="Comma-separated list of categories to include (filters the 'category' column)")
    parser.add_argument("--output_dir", type=str, default="output",
                        help="Directory to store output MEX directories")
    args = parser.parse_args()
    
    # Process filter categories if provided.
    filter_categories = None
    if args.filter_category:
        filter_categories = [cat.strip() for cat in args.filter_category.split(",")]
    
    # Load and merge input files.
    df = load_data(args.base)
    df = filter_data(df, filter_categories)
    
    gene_level_map = None
    transcript_level_map = None
    if args.gene_map:
        gene_level_map = load_gene_map(args.gene_map)
    if args.transcript_map:
        transcript_level_map = load_transcript_map(args.transcript_map)
    
    # Create gene-level MEX matrix.
    create_mex_matrices(df, group_by="gene", output_prefix=args.output_dir,
                        gene_level_map=gene_level_map, transcript_level_map=transcript_level_map)
    
    # Create transcript-level MEX matrix.
    create_mex_matrices(df, group_by="transcript", output_prefix=args.output_dir,
                        gene_level_map=gene_level_map, transcript_level_map=transcript_level_map)
    
if __name__ == "__main__":
    main()