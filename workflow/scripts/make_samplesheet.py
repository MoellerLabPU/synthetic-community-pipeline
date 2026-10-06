#!/usr/bin/env python
"""
Build the sample sheet (config/samples.tsv): one row per sequenced sample,
with its two FASTQ files and its row from the "Metadata" tab.

Run from the top folder of the repository, with the conda env active:

    python workflow/scripts/make_samplesheet.py
"""

import sys
from pathlib import Path

import pandas as pd
import yaml

# Metadata columns to read: name in the spreadsheet -> name in the sample sheet.
# Only these columns are read from the tab.
METADATA_COLUMNS = {
    "Sample_ID": "sample_id",
    "Plate_ID": "plate_id",
    "Plate Type": "plate_type",
    "Well": "well",
    "Passage_Day": "passage_day",
    "Incubation Type": "incubation_type",
    "Passage_Dilution": "passage_dilution",
    "Community_Starting_Concentration": "starting_concentration",
    "Community_Concentration": "community_cfu_per_ml",
    "Community_OD600": "community_od600",
    "DNA.Extraction.Quant.nguL": "dna_ng_per_ul",
    "Library.Prep.Quant.nguL": "library_ng_per_ul",
}

config = yaml.safe_load(open("config/config.yaml"))

# --- 1. One row per R1 file ---------------------------------------------------
# File names look like 461_SCP_SI_01_A01_10488880_253GGLLT4_L1_R1.fastq.gz:
# a running number, then plate type (SCP_SI / SCP_NSI), plate number and well.
rows = []
for r1 in sorted(Path(config["reads_dir"]).glob("*_R1.fastq.gz")):
    r2 = r1.with_name(r1.name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
    if not r2.is_file():
        sys.exit(f"ERROR: no R2 file for {r1.name}")
    # Drop the running number, keep the next four fields: SCP, SI, 01, A01
    scp, incubation, plate, well = r1.name.split("_")[1:5]
    rows.append(
        {
            "sample": f"{scp}_{incubation}_{plate}_{well}",     # used in file names
            "sample_id": f"{scp}_{incubation}.{plate}.{well}",  # as in the spreadsheet
            "fq1": r1.name,
            "fq2": r2.name,
        }
    )
if not rows:
    sys.exit(f"ERROR: no *_R1.fastq.gz files found in {config['reads_dir']}")
reads = pd.DataFrame(rows)

# Every sample must appear exactly once.
duplicates = reads.loc[reads["sample"].duplicated(), "sample"]
if not duplicates.empty:
    sys.exit(f"ERROR: more than one R1 file for: {', '.join(duplicates)}")

# --- 2. Metadata ----------------------------------------------------------------
# usecols reads only the columns listed above.
# dtype=str keeps values exactly as typed (for example "01" and "1/10").
metadata = pd.read_excel(
    config["metadata_xlsx"],
    sheet_name=config["metadata_sheet"],
    usecols=list(METADATA_COLUMNS),
    dtype=str,
).rename(columns=METADATA_COLUMNS)
# The tab has empty rows below the table; drop them.
metadata = metadata.dropna(subset=["sample_id"])

# --- 3. Join and write ------------------------------------------------------------
# validate= makes pandas stop if a sample ID appears more than once on either side.
sheet = reads.merge(metadata, on="sample_id", how="left", validate="one_to_one")
unmatched = sheet.loc[sheet["plate_id"].isna(), "sample"]
if not unmatched.empty:
    sys.exit(f"ERROR: FASTQ samples not found in the metadata: {', '.join(unmatched)}")

sheet.sort_values("sample").to_csv(config["samples"], sep="\t", index=False, na_rep="NA")
print(f"Wrote {len(sheet)} samples to {config['samples']}")
print(sheet.groupby(["passage_day", "incubation_type"]).size().to_string())
