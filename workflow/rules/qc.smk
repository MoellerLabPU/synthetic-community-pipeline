# =============================================================================
# Read quality control (pipeline steps 1-3)
#
#   fastqc_raw       FastQC on the raw reads
#   multiqc_raw      one report for all raw FastQC results          (step 1)
#   fastp            adapter and quality trimming                   (step 2)
#   fastqc_trimmed   FastQC on the trimmed reads
#   multiqc_trimmed  one report for trimmed FastQC + fastp results  (step 3)
#
# Uses from the main Snakefile: config, samples, SAMPLES, READS_DIR, RESULTS
# =============================================================================


rule fastqc_raw:
    """FastQC on the raw reads of one sample (R1 and R2 together)."""
    input:
        r1=lambda wc: READS_DIR / samples.at[wc.sample, "fq1"],
        r2=lambda wc: READS_DIR / samples.at[wc.sample, "fq2"],
    output:
        # One folder per sample. FastQC names the files inside after the
        # FASTQ files: <FASTQ name>_fastqc.html and <FASTQ name>_fastqc.zip
        directory(f"{RESULTS}/qc/raw/fastqc/{{sample}}"),
    threads: config["resources"]["fastqc"]["threads"]
    resources:
        mem_mb=config["resources"]["fastqc"]["mem_mb"],
        time=config["resources"]["fastqc"]["time"],
    shell:
        """
        # FastQC does not create its output folder
        mkdir -p {output}

        # --outdir   folder for the .html and .zip results
        # --threads  CPU threads
        fastqc \
            --outdir {output} \
            --threads {threads} \
            {input.r1} {input.r2}
        """


rule multiqc_raw:
    """Step 1: one MultiQC report for the raw reads of all samples."""
    input:
        expand(f"{RESULTS}/qc/raw/fastqc/{{sample}}", sample=SAMPLES),
    output:
        # MultiQC's default names
        html=f"{RESULTS}/qc/raw/multiqc_report.html",
        data=directory(f"{RESULTS}/qc/raw/multiqc_data"),
    params:
        outdir=f"{RESULTS}/qc/raw",
        fastqc_dir=f"{RESULTS}/qc/raw/fastqc",
    threads: config["resources"]["multiqc"]["threads"]
    resources:
        mem_mb=config["resources"]["multiqc"]["mem_mb"],
        time=config["resources"]["multiqc"]["time"],
    shell:
        """
        # --outdir  where the report is written
        multiqc \
            --outdir {params.outdir} \
            {params.fastqc_dir}
        """


rule fastp:
    """Step 2: trim adapters and low-quality ends from one sample."""
    input:
        r1=lambda wc: READS_DIR / samples.at[wc.sample, "fq1"],
        r2=lambda wc: READS_DIR / samples.at[wc.sample, "fq2"],
    output:
        r1=f"{RESULTS}/reads/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/reads/trimmed/{{sample}}_R2.fastq.gz",
        json=f"{RESULTS}/qc/fastp/{{sample}}.fastp.json",
        html=f"{RESULTS}/qc/fastp/{{sample}}.fastp.html",
    threads: config["resources"]["fastp"]["threads"]
    resources:
        mem_mb=config["resources"]["fastp"]["mem_mb"],
        time=config["resources"]["fastp"]["time"],
    shell:
        """
        # --in1, --in2             raw R1 and R2
        # --out1, --out2           trimmed R1 and R2
        # --detect_adapter_for_pe  find adapters from the read-pair overlap
        # --json                   report for MultiQC
        # --html                   report to open in a browser
        # --thread                 CPU cores
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --detect_adapter_for_pe \
            --json {output.json} \
            --html {output.html} \
            --thread {threads}
        """


rule fastqc_trimmed:
    """FastQC on the trimmed reads of one sample (R1 and R2 together)."""
    input:
        r1=f"{RESULTS}/reads/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/reads/trimmed/{{sample}}_R2.fastq.gz",
    output:
        # One folder per sample, same layout as for the raw reads
        directory(f"{RESULTS}/qc/trimmed/fastqc/{{sample}}"),
    threads: config["resources"]["fastqc"]["threads"]
    resources:
        mem_mb=config["resources"]["fastqc"]["mem_mb"],
        time=config["resources"]["fastqc"]["time"],
    shell:
        """
        # FastQC does not create its output folder
        mkdir -p {output}

        # --outdir   folder for the .html and .zip results
        # --threads  CPU threads
        fastqc \
            --outdir {output} \
            --threads {threads} \
            {input.r1} {input.r2}
        """


rule multiqc_trimmed:
    """Step 3: one MultiQC report for the trimmed reads, with the fastp results."""
    input:
        expand(f"{RESULTS}/qc/trimmed/fastqc/{{sample}}", sample=SAMPLES),
        expand(f"{RESULTS}/qc/fastp/{{sample}}.fastp.json", sample=SAMPLES),
    output:
        # MultiQC's default names
        html=f"{RESULTS}/qc/trimmed/multiqc_report.html",
        data=directory(f"{RESULTS}/qc/trimmed/multiqc_data"),
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
        # --outdir  where the report is written
        multiqc \
            --outdir {params.outdir} \
            {params.fastqc_dir} {params.fastp_dir}
        """
