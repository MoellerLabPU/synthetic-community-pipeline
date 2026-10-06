#!/usr/bin/env python
"""
Build the sample sheet (config/samples.tsv) for the pipeline.

It does two things:
  1. Lists the FASTQ files in the reads folder and pairs R1 with R2.
  2. Looks up every sample in the "Metadata" tab of the experiment spreadsheet
     and copies over the columns that describe the experiment.

The sample sheet is built from the FASTQ files that actually exist. Wells that
are in the metadata but have no FASTQ files are reported and left out.

Usage (from the top folder of the repository, with the conda env active):

    python workflow/scripts/make_samplesheet.py --config config/config.yaml
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

# FASTQ file names look like this:
#   461_SCP_SI_01_A01_10488880_253GGLLT4_L1_R1.fastq.gz
#   ^   ^      ^  ^                         ^
#   |   |      |  well                      read (R1 or R2)
#   |   |      plate number
#   |   plate type (SCP_SI = shaking, SCP_NSI = non-shaking)
#   running number given by the sequencing facility
FASTQ_PATTERN = re.compile(
    r"^(?P<seq_number>\d+)_"
    r"(?P<plate_type>SCP_N?SI)_"
    r"(?P<plate_number>\d+)_"
    r"(?P<well>[A-H]\d{2})_"
    r".*_R(?P<read>[12])\.fastq\.gz$"
)

# Metadata columns to keep: name in the spreadsheet -> name in the sample sheet.
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


def list_fastq_pairs(reads_dir: Path) -> pd.DataFrame:
    """Return one row per sample with its R1 and R2 file names."""
    rows = {}
    for path in sorted(reads_dir.glob("*.fastq.gz")):
        match = FASTQ_PATTERN.match(path.name)
        if match is None:
            sys.exit(f"ERROR: FASTQ name does not fit the expected pattern: {path.name}")
        parts = match.groupdict()
        # "sample" is used in output file names, so it has no dots: SCP_SI_01_A01
        sample = f"{parts['plate_type']}_{parts['plate_number']}_{parts['well']}"
        # "sample_id" is written the way the spreadsheet writes it: SCP_SI.01.A01
        sample_id = f"{parts['plate_type']}.{parts['plate_number']}.{parts['well']}"
        row = rows.setdefault(sample, {"sample": sample, "sample_id": sample_id})
        column = f"fq{parts['read']}"
        if column in row:
            sys.exit(f"ERROR: more than one R{parts['read']} file for sample {sample}")
        row[column] = path.name

    pairs = pd.DataFrame(rows.values())
    if pairs.empty:
        sys.exit(f"ERROR: no *.fastq.gz files found in {reads_dir}")
    # Every sample must have both mates.
    unpaired = pairs[pairs[["fq1", "fq2"]].isna().any(axis=1)]
    if not unpaired.empty:
        sys.exit(f"ERROR: samples without both R1 and R2: {', '.join(unpaired['sample'])}")
    return pairs


def read_metadata(xlsx: Path, sheet: str) -> pd.DataFrame:
    """Read the metadata tab and keep the columns listed in METADATA_COLUMNS."""
    # dtype=str keeps values exactly as typed (for example "01" and "1/10").
    metadata = pd.read_excel(xlsx, sheet_name=sheet, dtype=str)
    missing = [col for col in METADATA_COLUMNS if col not in metadata.columns]
    if missing:
        sys.exit(f"ERROR: columns missing from tab '{sheet}': {', '.join(missing)}")
    metadata = metadata[list(METADATA_COLUMNS)].rename(columns=METADATA_COLUMNS)
    # The tab has empty rows below the table; drop them.
    metadata = metadata.dropna(subset=["sample_id"])
    if metadata["sample_id"].duplicated().any():
        sys.exit(f"ERROR: duplicated Sample_ID values in tab '{sheet}'")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml", help="pipeline config file")
    args = parser.parse_args()

    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    reads_dir = Path(config["reads_dir"])
    xlsx = Path(config["metadata_xlsx"])
    out_path = Path(config["samples"])

    # Check the inputs before doing anything.
    if not reads_dir.is_dir():
        sys.exit(f"ERROR: reads folder not found: {reads_dir}")
    if not xlsx.is_file():
        sys.exit(f"ERROR: metadata spreadsheet not found: {xlsx}")

    pairs = list_fastq_pairs(reads_dir)
    metadata = read_metadata(xlsx, config["metadata_sheet"])

    # Join FASTQ pairs to metadata. validate= makes pandas fail on duplicates.
    sheet = pairs.merge(metadata, on="sample_id", how="left", validate="one_to_one", indicator=True)
    not_in_metadata = sheet.loc[sheet["_merge"] == "left_only", "sample"]
    if not not_in_metadata.empty:
        sys.exit(f"ERROR: FASTQ samples not found in the metadata: {', '.join(not_in_metadata)}")
    sheet = sheet.drop(columns="_merge")

    # Report wells that have a sequencing concentration but no FASTQ files.
    sequenced = metadata[metadata["library_ng_per_ul"].notna() & (metadata["library_ng_per_ul"] != "NA")]
    no_fastq = sequenced[~sequenced["sample_id"].isin(sheet["sample_id"])]
    if not no_fastq.empty:
        per_plate = no_fastq["plate_id"].value_counts().sort_index()
        print(f"NOTE: {len(no_fastq)} wells have library data in the metadata but no FASTQ files:")
        for plate, count in per_plate.items():
            print(f"        {plate}: {count} wells")

    if sheet.empty:
        sys.exit("ERROR: the sample sheet is empty, nothing written")

    sheet = sheet.sort_values("sample")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sheet.to_csv(out_path, sep="\t", index=False, na_rep="NA")

    print(f"Wrote {len(sheet)} samples to {out_path}")
    print(sheet.groupby(["passage_day", "incubation_type"]).size().to_string())


if __name__ == "__main__":
    main()
