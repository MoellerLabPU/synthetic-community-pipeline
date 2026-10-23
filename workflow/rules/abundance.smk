# =============================================================================
# Relative abundance (pipeline step 8)
#
#   samtools_coverage  per-contig depth and breadth of one BAM
#   abundance_table    sums the contigs per strain and computes depth, breadth
#                      and relative abundance for all samples, joined to the
#                      metadata
#   abundance_plots    composition over time per condition
#
# Uses from the main Snakefile: config, RESULTS, MAP_SAMPLES, PRESET, SCRIPTS, SAMPLES_TSV
# =============================================================================


rule samtools_coverage:
    """Per-contig coverage of one BAM."""
    input:
        bam=f"{RESULTS}/mapping/{PRESET}/{{sample}}.bam",
        bai=f"{RESULTS}/mapping/{PRESET}/{{sample}}.bam.bai",
    output:
        f"{RESULTS}/coverage/{PRESET}/{{sample}}.coverage.tsv",
    threads: 1
    shell:
        """
        # One line per contig: numreads, covbases (bases with at least one
        # read), coverage (% of the contig covered), meandepth, ...
        samtools coverage {input.bam} > {output}
        """


rule abundance_table:
    """One table with depth, breadth and relative abundance per sample and strain."""
    input:
        coverage=expand(f"{RESULTS}/coverage/{PRESET}/{{sample}}.coverage.tsv", sample=MAP_SAMPLES),
        samples=SAMPLES_TSV,
    output:
        f"{RESULTS}/abundance/{PRESET}/abundance.tsv",
    threads: 1
    shell:
        """
        python {SCRIPTS}/abundance.py \
            --samples {input.samples} \
            --out {output} \
            --coverage {input.coverage}
        """


rule abundance_plots:
    """Composition over passage day, one panel per condition."""
    input:
        f"{RESULTS}/abundance/{PRESET}/abundance.tsv",
    output:
        # the script names the files; the rule only chooses the folder
        f"{RESULTS}/abundance/{PRESET}/composition_wells.png",
        f"{RESULTS}/abundance/{PRESET}/composition_means.png",
        f"{RESULTS}/abundance/{PRESET}/composition_wells_contig_average.png",
        f"{RESULTS}/abundance/{PRESET}/composition_means_contig_average.png",
    params:
        outdir=f"{RESULTS}/abundance/{PRESET}",
    threads: 1
    shell:
        """
        python {SCRIPTS}/plot_abundance.py \
            --table {input} \
            --outdir {params.outdir}
        """
