# =============================================================================
# Read mapping (pipeline steps 6-7)
#
#   bowtie2_map   map the trimmed reads of one sample to the combined reference,
#                 sort and index the BAM
#
# Uses from the main Snakefile: config, RESULTS, MAP_SAMPLES
# =============================================================================

# Is --very-sensitive switched on in the config? Compared as text so that
# "true"/"false" given on the command line (--config) work like the YAML booleans.
VERY_SENSITIVE = str(config["bowtie2"]["very_sensitive"]).lower() == "true"

# Name of the Bowtie2 preset, also used as a folder name so that runs with the
# two presets can sit side by side: results/mapping/sensitive/ and
# results/mapping/very_sensitive/
PRESET = "very_sensitive" if VERY_SENSITIVE else "sensitive"


rule bowtie2_map:
    """Map one sample's trimmed reads to the combined reference."""
    input:
        r1=f"{RESULTS}/fastp/trimmed_reads/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/fastp/trimmed_reads/{{sample}}_R2.fastq.gz",
        index=expand(f"{RESULTS}/reference/combined.{{n}}.bt2", n=["1", "2", "3", "4", "rev.1", "rev.2"]),
    output:
        bam=f"{RESULTS}/mapping/{PRESET}/{{sample}}.bam",
        bai=f"{RESULTS}/mapping/{PRESET}/{{sample}}.bam.bai",
    params:
        index=f"{RESULTS}/reference/combined",
        # "--very-sensitive", or nothing for Bowtie2's default (--sensitive)
        preset="--very-sensitive" if VERY_SENSITIVE else "",
    threads: config["resources"]["threads"]
    shell:
        """
        # bowtie2
        #   -x          index prefix
        #   -1, -2      trimmed R1 and R2
        #   --threads   CPU cores
        #   --no-unal   leave unmapped reads out of the BAM (smaller file)
        #   Bowtie2 prints its alignment summary (overall alignment rate) to
        #   stderr, which ends up in the job's SLURM log under logs/slurm/.
        # samtools sort
        #   -@          extra CPU cores
        #   -o          sorted BAM
        #   -           read the alignments from bowtie2 through the pipe
        bowtie2 \
            -x {params.index} \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            --no-unal \
            {params.preset} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam} \
            -

        samtools index -@ {threads} {output.bam}
        """
