# =============================================================================
# Read quality control (pipeline steps 1-3)
#
#   link_raw_reads   give the raw FASTQ files short names (symlinks, no copy)
#   fastqc_raw       FastQC on the raw reads
#   multiqc_raw      one report for all raw FastQC results          (step 1)
#   fastp            adapter and quality trimming                   (step 2)
#   fastqc_trimmed   FastQC on the trimmed reads
#   multiqc_trimmed  one report for trimmed FastQC + fastp results  (step 3)
#
# Uses from the main Snakefile: config, samples, SAMPLES, READS_DIR, RESULTS
# =============================================================================


# Tiny jobs that are run directly instead of being submitted to SLURM.
localrules:
    link_raw_reads,


rule link_raw_reads:
    """
    Link the raw FASTQ files under short names: <sample>_R1.fastq.gz.
    FastQC and MultiQC name their results after the input file, so this gives
    clean sample names in every report.
    """
    input:
        r1=lambda wc: READS_DIR / samples.at[wc.sample, "fq1"],
        r2=lambda wc: READS_DIR / samples.at[wc.sample, "fq2"],
    output:
        r1=f"{RESULTS}/reads/raw/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/reads/raw/{{sample}}_R2.fastq.gz",
    shell:
        """
        # -s symbolic link, -r path relative to the link, -f replace an old link
        ln -srf {input.r1} {output.r1}
        ln -srf {input.r2} {output.r2}
        """


rule fastqc_raw:
    """FastQC on the raw reads of one sample (R1 and R2 together)."""
    input:
        r1=f"{RESULTS}/reads/raw/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/reads/raw/{{sample}}_R2.fastq.gz",
    output:
        # FastQC names its outputs <input file name without .fastq.gz>_fastqc.*
        html=expand(f"{RESULTS}/qc/raw/fastqc/{{{{sample}}}}_{{read}}_fastqc.html", read=["R1", "R2"]),
        zip=expand(f"{RESULTS}/qc/raw/fastqc/{{{{sample}}}}_{{read}}_fastqc.zip", read=["R1", "R2"]),
    params:
        outdir=f"{RESULTS}/qc/raw/fastqc",
    threads: config["resources"]["fastqc"]["threads"]
    resources:
        mem_mb=config["resources"]["fastqc"]["mem_mb"],
        time=config["resources"]["fastqc"]["time"],
    shell:
        """
        # --outdir   folder for the .html and .zip results
        # --threads  CPU threads
        # --quiet    only print errors
        fastqc \
            --outdir {params.outdir} \
            --threads {threads} \
            --quiet \
            {input.r1} {input.r2}
        """


rule multiqc_raw:
    """Step 1: one MultiQC report for the raw reads of all samples."""
    input:
        expand(f"{RESULTS}/qc/raw/fastqc/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
    output:
        html=f"{RESULTS}/qc/raw/multiqc_report.html",
        data=directory(f"{RESULTS}/qc/raw/multiqc_report_data"),
    params:
        outdir=f"{RESULTS}/qc/raw",
        fastqc_dir=f"{RESULTS}/qc/raw/fastqc",
    threads: config["resources"]["multiqc"]["threads"]
    resources:
        mem_mb=config["resources"]["multiqc"]["mem_mb"],
        time=config["resources"]["multiqc"]["time"],
    shell:
        """
        # --outdir    where the report is written
        # --filename  name of the report
        # --title     title shown at the top
        # --force     overwrite an older report
        multiqc \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --title "SCP raw reads" \
            --force \
            {params.fastqc_dir}
        """


rule fastp:
    """Step 2: trim adapters and low-quality ends from one sample."""
    input:
        r1=f"{RESULTS}/reads/raw/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/reads/raw/{{sample}}_R2.fastq.gz",
    output:
        r1=f"{RESULTS}/reads/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/reads/trimmed/{{sample}}_R2.fastq.gz",
        json=f"{RESULTS}/qc/fastp/{{sample}}.fastp.json",
        html=f"{RESULTS}/qc/fastp/{{sample}}.fastp.html",
    params:
        min_length=config["fastp"]["min_length"],
    threads: config["resources"]["fastp"]["threads"]
    resources:
        mem_mb=config["resources"]["fastp"]["mem_mb"],
        time=config["resources"]["fastp"]["time"],
    shell:
        """
        # --in1, --in2             raw R1 and R2
        # --out1, --out2           trimmed R1 and R2
        # --detect_adapter_for_pe  find adapters from the read-pair overlap
        # --length_required        drop reads shorter than this after trimming
        # --json                   report for MultiQC
        # --html                   report to open in a browser
        # --report_title           title of the html report
        # --thread                 CPU cores
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --detect_adapter_for_pe \
            --length_required {params.min_length} \
            --json {output.json} \
            --html {output.html} \
            --report_title {wildcards.sample} \
            --thread {threads}
        """


rule fastqc_trimmed:
    """FastQC on the trimmed reads of one sample (R1 and R2 together)."""
    input:
        r1=f"{RESULTS}/reads/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/reads/trimmed/{{sample}}_R2.fastq.gz",
    output:
        html=expand(f"{RESULTS}/qc/trimmed/fastqc/{{{{sample}}}}_{{read}}_fastqc.html", read=["R1", "R2"]),
        zip=expand(f"{RESULTS}/qc/trimmed/fastqc/{{{{sample}}}}_{{read}}_fastqc.zip", read=["R1", "R2"]),
    params:
        outdir=f"{RESULTS}/qc/trimmed/fastqc",
    threads: config["resources"]["fastqc"]["threads"]
    resources:
        mem_mb=config["resources"]["fastqc"]["mem_mb"],
        time=config["resources"]["fastqc"]["time"],
    shell:
        """
        # --outdir   folder for the .html and .zip results
        # --threads  CPU threads
        # --quiet    only print errors
        fastqc \
            --outdir {params.outdir} \
            --threads {threads} \
            --quiet \
            {input.r1} {input.r2}
        """


rule multiqc_trimmed:
    """Step 3: one MultiQC report for the trimmed reads, with the fastp results."""
    input:
        expand(f"{RESULTS}/qc/trimmed/fastqc/{{sample}}_{{read}}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
        expand(f"{RESULTS}/qc/fastp/{{sample}}.fastp.json", sample=SAMPLES),
    output:
        html=f"{RESULTS}/qc/trimmed/multiqc_report.html",
        data=directory(f"{RESULTS}/qc/trimmed/multiqc_report_data"),
    params:
        outdir=f"{RESULTS}/qc/trimmed",
        fastqc_dir=f"{RESULTS}/qc/trimmed/fastqc",
        fastp_dir=f"{RESULTS}/qc/fastp",
    threads: config["resources"]["multiqc"]["threads"]
    resources:
        mem_mb=config["resources"]["multiqc"]["mem_mb"],
        time=config["resources"]["multiqc"]["time"],
    shell:
        """
        # --outdir    where the report is written
        # --filename  name of the report
        # --title     title shown at the top
        # --force     overwrite an older report
        multiqc \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --title "SCP trimmed reads" \
            --force \
            {params.fastqc_dir} {params.fastp_dir}
        """
