# Utilities for IsoMEX

This directory contains helper scripts for preprocessing files used in IsoMEX.

## Dependencies

- Python 3.6+
- [pandas](https://pandas.pydata.org/)
- [gffutils](https://daler.github.io/gffutils/)
- Standard Python libraries: `argparse`

Install the dependencies via pip (if not already installed):

```bash
pip install pandas gffutils
```

## Scripts

### **`generate_map.py`**  
Extracts **gene IDs & names** and **transcript IDs & names** from a **GTF/GFF annotation file**, generating two tab-separated mapping files.

**Usage:**
```bash
python generate_map.py annotation.forPigeon.gtf gene_map.txt transcript_map.txt
```

**Inputs:**
- `annotation.forPigeon.gtf` → Input GTF annotation file (can use the same one made with [pigeon prepare](https://isoseq.how/classification/workflow.html))

**Outputs:**
- `gene_map.txt` → Contains `gene_id` and `gene_name`
- `transcript_map.txt` → Contains `transcript_id` and `transcript_name`

---

## **How to Use These Scripts**
1. Ensure all dependencies are installed.
2. Run `generate_map.py` with your **GTF annotation file** to generate `gene_map.txt` and `transcript_map.txt`.
3. Use the output files as input for IsoMEX.

---

### **Future Utilities**
- More utility scripts may be added to handle **data processing, file formatting, and annotation refinement**.

---
