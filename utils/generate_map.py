#!/usr/bin/env python3

import argparse
import gffutils
import pandas as pd

def generate_mappings(gtf_file, gene_map_out, transcript_map_out):
    """
    Parses a GTF file and generates:
    - A gene mapping file (gene_id → gene_name)
    - A transcript mapping file (transcript_id → transcript_name)
    """

    # Create an in-memory database from the GTF file
    db = gffutils.create_db(gtf_file, dbfn=":memory:", force=True, keep_order=True,
                            merge_strategy="merge", sort_attribute_values=True)

    # Generate gene map dataframe
    gene_map = []
    for gene in db.features_of_type("gene"):
        gene_id = gene.attributes.get("gene_id", [""])[0]
        gene_name = gene.attributes.get("gene_name", [gene_id])[0]
        gene_map.append({"gene_id": gene_id, "gene_name": gene_name})
    gene_df = pd.DataFrame(gene_map)
    gene_df.to_csv("gene_map.txt", sep="\t", index=False)
    print(f"Gene map saved to {gene_map_out}")

    # Generate transcript map dataframe
    transcript_map = []
    for transcript in db.features_of_type("transcript"):
        transcript_id = transcript.attributes.get("transcript_id", [""])[0]
        transcript_name = transcript.attributes.get("transcript_name", [transcript_id])[0]
        transcript_map.append({"transcript_id": transcript_id, "transcript_name": transcript_name})
    transcript_df = pd.DataFrame(transcript_map)
    transcript_df.to_csv("transcript_map.txt", sep="\t", index=False)
    print(f"Transcript map saved to {transcript_map_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate gene and transcript mapping files from a GTF annotation.")
    parser.add_argument("gtf_file", help="Input GTF annotation file")
    parser.add_argument("gene_map_out", help="Output file for gene mappings")
    parser.add_argument("transcript_map_out", help="Output file for transcript mappings")
    
    args = parser.parse_args()
    generate_mappings(args.gtf_file, args.gene_map_out, args.transcript_map_out)