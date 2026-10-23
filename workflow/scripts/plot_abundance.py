#!/usr/bin/env python
"""
Plot the relative abundance table from abundance.py with plotnine (ggplot2 for
Python). Four files are written to the output directory, two per abundance
measure (rel_abundance from the length-weighted depth, and
rel_abundance_contig_average from the plain average of contig depths):

composition_wells[_contig_average].png
    one stacked bar per well; rows = starting concentration x dilution,
    columns = incubation type x passage day
composition_means[_contig_average].png
    per strain, the mean over the wells of a day as a line, with one point
    per well; one panel per condition

Usage: plot_abundance.py --table abundance.tsv --outdir results/abundance/<preset>
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_text,
    facet_grid,
    facet_wrap,
    geom_col,
    geom_point,
    ggplot,
    labs,
    scale_colour_manual,
    scale_fill_manual,
    stat_summary,
    theme,
    theme_bw,
)

# Fixed colour per strain (in sorted strain order), the same in every plot
COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument(
    "--table", required=True, help="abundance table written by abundance.py (TSV)"
)
parser.add_argument("--outdir", required=True, help="directory for the two PNG files")
args = parser.parse_args()
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# dtype=str keeps "10000000" and "1/10" as written instead of as numbers
table = pd.read_csv(
    args.table, sep="\t", dtype={"starting_concentration": str, "passage_dilution": str}
)
strains = sorted(table["strain"].unique())
colour = dict(zip(strains, COLOURS))

# Panel labels. The stacked-bar figure uses rows = starting concentration x
# dilution and columns = incubation type x day (12 x 6 panels); the means
# figure uses one panel per full condition (24 panels, 4 per row).
table["concentration_dilution"] = (
    table["starting_concentration"] + " CFU/mL\ndilution " + table["passage_dilution"]
)
table["incubation_day"] = (
    table["incubation_type"] + "\nday " + table["passage_day"].astype(str)
)
table["condition"] = (
    table["incubation_type"]
    + "\n"
    + table["starting_concentration"]
    + " CFU/mL\n"
    + "dilution "
    + table["passage_dilution"]
)
n_rows = table["concentration_dilution"].nunique()
n_conditions = table["condition"].nunique()
# Replicate number (1, 2, ...) of each well within its condition and day, so
# that every panel has the same x positions whatever the wells are called
table["replicate"] = (
    table.groupby(["condition", "passage_day"])["sample"]
    .rank(method="dense")
    .astype(int)
)

# The same two figures for each abundance measure
MEASURES = {
    "rel_abundance": ("", "length-weighted depth"),
    "rel_abundance_contig_average": (
        "_contig_average",
        "plain average of contig depths",
    ),
}

for measure, (suffix, description) in MEASURES.items():
    # --- one stacked bar per well: rows = conditions, columns = days -----------
    wells = (
        ggplot(table, aes(x="factor(replicate)", y=measure, fill="strain"))
        + geom_col()  # bars stack by default
        + facet_grid("concentration_dilution ~ incubation_day")
        + scale_fill_manual(values=colour)
        + labs(
            x="replicate well",
            y="relative abundance",
            fill="strain",
            title=f"Relative abundance from {description}",
        )
        + theme_bw()
        + theme(strip_text_y=element_text(angle=0))  # condition labels horizontal
    )
    wells.save(
        outdir / f"composition_wells{suffix}.png",
        width=13,
        height=1.6 * n_rows + 1,
        dpi=150,
        verbose=False,
    )

    # --- mean per day as a line, one point per well ---------------------------
    means = (
        # factor(): only the days that exist appear on the x axis; group= joins
        # the points of one strain into a line
        ggplot(
            table,
            aes(x="factor(passage_day)", y=measure, colour="strain", group="strain"),
        )
        + stat_summary(
            fun_y=np.mean, geom="line", size=1
        )  # mean over the wells of a day
        + geom_point(alpha=0.7)  # one point per well
        + facet_wrap("condition", ncol=4)
        + scale_colour_manual(values=colour)
        + labs(
            x="passage day",
            y="relative abundance",
            colour="strain",
            title=f"Relative abundance from {description}",
        )
        + theme_bw()
    )
    means.save(
        outdir / f"composition_means{suffix}.png",
        width=12,
        height=3 * -(-n_conditions // 4) + 1,
        dpi=150,
        verbose=False,
    )
print(f"Wrote {2 * len(MEASURES)} figures to {outdir}")
