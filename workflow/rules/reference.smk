# =============================================================================
# Reference genomes (pipeline steps 4-5)
#
#   combine_genomes  one FASTA with all 5 genomes, contigs prefixed by strain  (step 4)
#   bowtie2_index    Bowtie2 index of the combined FASTA                        (step 4)
#   ani_skani        all-vs-all ANI between the 5 genomes with skani            (step 5)
#   fastani_genome_list  list of genome files for FastANI
#   ani_fastani      the same with FastANI                                      (step 5)
#
# Uses from the main Snakefile: config, RESULTS
# =============================================================================


rule combine_genomes:
    """
    Step 4a: put all genomes into one FASTA. Every contig name gets the strain ID
    in front, separated by "|", for example "PSL0003|012_PSL0003_..._00001", so
    that coverage per contig can later be summed per strain.
    """
    input:
        # genome FASTA files, in the order they are listed in the config
        list(config["genomes"].values()),
    output:
        fasta=f"{RESULTS}/reference/combined.fasta",
        fai=f"{RESULTS}/reference/combined.fasta.fai",
    params:
        # "PSL0003=/path/to/genome.fna.gz PSL0034=..." for the loop below
        genomes=" ".join(f"{name}={path}" for name, path in config["genomes"].items()),
    shell:
        """
        # Start from an empty file, then append one genome at a time
        : > {output.fasta}
        for entry in {params.genomes}; do
            name=${{entry%%=*}}    # text before "=": the strain ID
            path=${{entry#*=}}     # text after "=": the FASTA file
            # zcat -f reads both gzipped and plain FASTA.
            # On header lines (starting with ">") keep only the first word and
            # put the strain ID in front; sequence lines are copied as they are.
            zcat -f "$path" \
                | awk -v name="$name" '/^>/ {{print ">" name "|" substr($1, 2); next}} {{print}}' \
                >> {output.fasta}
        done

        # Index: lists every contig with its length (used for the per-genome sizes)
        samtools faidx {output.fasta}
        """


rule bowtie2_index:
    """Step 4b: build the Bowtie2 index of the combined FASTA."""
    input:
        f"{RESULTS}/reference/combined.fasta",
    output:
        # Bowtie2 writes six files with this prefix
        expand(f"{RESULTS}/reference/combined.{{n}}.bt2", n=["1", "2", "3", "4", "rev.1", "rev.2"]),
    params:
        prefix=f"{RESULTS}/reference/combined",
    threads: config["resources"]["threads"]
    shell:
        """
        # --threads  CPU cores
        bowtie2-build \
            --threads {threads} \
            {input} {params.prefix}
        """


rule ani_skani:
    """
    Step 5: average nucleotide identity (ANI) between every pair of genomes,
    with skani. Shows how similar the two Leuconostoc genomes are, which
    decides how hard it is to tell their reads apart.
    """
    input:
        list(config["genomes"].values()),
    output:
        ani=f"{RESULTS}/ani/skani.tsv",
        # skani also writes the aligned fraction of each pair next to it
        af=f"{RESULTS}/ani/skani.tsv.af",
    threads: config["resources"]["threads"]
    shell:
        """
        # triangle       compare every genome with every other
        # --full-matrix  write the ANI values as a genome x genome table
        # -t             CPU cores
        # -o             output file
        skani triangle \
            --full-matrix \
            -t {threads} \
            -o {output.ani} \
            {input}
        """


rule fastani_genome_list:
    """FastANI takes the genomes as a list file, one path per line."""
    input:
        list(config["genomes"].values()),
    output:
        f"{RESULTS}/ani/fastani_genome_list.txt",
    localrule: True
    shell:
        """
        printf "%s\n" {input} > {output}
        """


rule ani_fastani:
    """Step 5, same comparison with FastANI, for a second opinion."""
    input:
        f"{RESULTS}/ani/fastani_genome_list.txt",
    output:
        pairs=f"{RESULTS}/ani/fastani.tsv",
        matrix=f"{RESULTS}/ani/fastani.tsv.matrix",
    threads: config["resources"]["threads"]
    shell:
        """
        # --ql, --rl  query and reference lists: the same list, so all-vs-all
        # --matrix    also write the values as a lower-triangular table (.matrix)
        # -t          CPU cores
        # -o          output file: one line per pair (query, reference, ANI,
        #             matching fragments, total fragments)
        fastANI \
            --ql {input} \
            --rl {input} \
            --matrix \
            -t {threads} \
            -o {output.pairs}
        """
