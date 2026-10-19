#!/usr/bin/env python
"""
Sum the per-contig output of `samtools coverage` per strain and compute, for
every sample and strain:

  mapped_reads              reads mapped to the strain's contigs
  mean_depth                mapped bases / genome length: mean depth over the
                            whole genome, positions without reads count as
                            zero (each contig weighted by its length)
  mean_depth_contig_average plain average of the contigs' mean depths, every
                            contig counting the same whatever its length
  breadth                   fraction of the genome covered by at least one read
  rel_abundance             mean_depth / sum of mean_depth over all strains
  rel_abundance_contig_average  the same, from mean_depth_contig_average

Metadata columns from the sample sheet are added to every row.

Usage: abundance.py --samples config/samples.tsv --out table.tsv --coverage <files...>

The coverage files are named <sample>.coverage.tsv; the sample name is read
from the file name.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument(
    "--samples", required=True, help="sample sheet (config/samples.tsv)"
)
parser.add_argument("--out", required=True, help="output table (TSV)")
parser.add_argument(
    "--coverage", nargs="+", required=True, help="samtools coverage files"
)
args = parser.parse_args()

rows = []
for path in args.coverage:
    sample = Path(path).name.split(".")[
        0
    ]  # SCP_SI_01_A01.coverage.tsv -> SCP_SI_01_A01

    cov = pd.read_csv(path, sep="\t")
    if cov.empty:
        sys.exit(f"ERROR: {path} has no rows")
    # Contig names are "<strain>|<contig>"; keep the strain part
    cov["strain"] = cov["#rname"].str.split("|").str[0]
    cov["length"] = cov["endpos"] - cov["startpos"] + 1
    # meandepth is per contig; multiply by length to get the mapped bases of
    # the contig, so that long contigs count more than short ones in the sum
    cov["mapped_bases"] = cov["meandepth"] * cov["length"]

    per_strain = cov.groupby("strain").agg(
        genome_length=("length", "sum"),
        mapped_reads=("numreads", "sum"),
        mapped_bases=("mapped_bases", "sum"),
        covered_bases=("covbases", "sum"),
        mean_depth_contig_average=("meandepth", "mean"),
    )
    per_strain["mean_depth"] = per_strain["mapped_bases"] / per_strain["genome_length"]
    per_strain["breadth"] = per_strain["covered_bases"] / per_strain["genome_length"]
    for depth, rel in [
        ("mean_depth", "rel_abundance"),
        ("mean_depth_contig_average", "rel_abundance_contig_average"),
    ]:
        total = per_strain[depth].sum()
        per_strain[rel] = per_strain[depth] / total if total > 0 else 0.0
    per_strain = per_strain.reset_index()
    per_strain.insert(0, "sample", sample)
    rows.append(per_strain)

table = pd.concat(rows, ignore_index=True)

# Add the metadata columns from the sample sheet
samples = pd.read_csv(args.samples, sep="\t", dtype=str)
metadata_columns = [
    "sample",
    "passage_day",
    "incubation_type",
    "passage_dilution",
    "starting_concentration",
    "plate_id",
    "well",
]
# indicator=True adds a "_merge" column saying whether each row found its
# sample in the sheet ("both") or not ("left_only")
table = table.merge(
    samples[metadata_columns],
    on="sample",
    how="left",
    validate="many_to_one",
    indicator=True,
)
if (table["_merge"] == "left_only").any():
    missing = table.loc[table["_merge"] == "left_only", "sample"].unique()
    sys.exit(f"ERROR: samples not in the sample sheet: {', '.join(missing)}")

columns = metadata_columns + [
    "strain",
    "genome_length",
    "mapped_reads",
    "mean_depth",
    "mean_depth_contig_average",
    "breadth",
    "rel_abundance",
    "rel_abundance_contig_average",
]
table = table[columns].sort_values(["sample", "strain"])
table.to_csv(args.out, sep="\t", index=False, float_format="%.6g")
print(
    f"Wrote {len(table)} rows ({table['sample'].nunique()} samples x {table['strain'].nunique()} strains) to {args.out}"
)
