# Synthetic community pilot: QC, mapping and relative abundance

A Snakemake pipeline that takes the shotgun reads from the Synthetic Community
Pilot (SCP) and reports how much of each of the five strains is in every sample.

| Strain  | Species (GTDB)              |
|---------|-----------------------------|
| PSL0003 | *Escherichia coli*          |
| PSL0034 | *Leuconostoc lactis*        |
| PSL0019 | *Leuconostoc lactis_A*      |
| PSL0037 | *Enterococcus_D gallinarum* |
| PSL0033 | *Phocaeicola sartorii*      |

## What the pipeline does

| Step | What happens                                                        | Status      |
|------|---------------------------------------------------------------------|-------------|
| 1    | FastQC + MultiQC on the raw reads                                   | written     |
| 2    | fastp trimming                                                      | written     |
| 3    | FastQC + MultiQC on the trimmed reads                               | written     |
| 4    | Combine the 5 genomes into one reference and build a Bowtie2 index  | written     |
| 5    | All-vs-all ANI between the 5 genomes (skani and FastANI)            | written     |
| 6    | Bowtie2 pilot on 12 samples: 2 presets x 2 MAPQ filters             | written     |
| 7    | Bowtie2 on all samples, sorted and indexed BAMs                     | written     |
| 8    | Coverage per genome, relative abundance, tables and plots           | written     |

The steps are being added one at a time. This table is updated as each one lands.

## Repository layout

```
config/config.yaml          all paths and settings (the only file you edit)
config/samples.tsv          one row per sample, built from the metadata spreadsheet
envs/environment.yml        the conda environment with every tool
profiles/slurm/config.yaml  tells Snakemake how to submit jobs to SLURM
workflow/Snakefile          loads the settings and the sample sheet
workflow/rules/qc.smk       steps 1-3: FastQC, fastp, MultiQC
workflow/rules/reference.smk  steps 4-5: combined reference, Bowtie2 index, ANI
workflow/rules/mapping.smk    steps 6-7: Bowtie2 mapping
workflow/rules/abundance.smk  step 8: coverage, relative abundance, plots
workflow/scripts/           small Python scripts called by the pipeline
results/                    everything the pipeline produces
logs/slurm/                 one log per job, written by SLURM
```

## 1. One-time setup

Do this once on the lab server (`cbsumoeller02`).

### 1.1 Get the code

```bash
cd /workdir/$USER
git clone https://github.com/MoellerLabPU/synthetic-community-pipeline.git
cd synthetic-community-pipeline
```

### 1.2 Install conda (skip if `conda --version` already works)

```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh     # accept the defaults, answer "yes" to init
```

Close the terminal and log in again so that `conda` is found.

### 1.3 Create the environment

```bash
conda env create -f envs/environment.yml
```

This installs Snakemake, the SLURM plugin, FastQC, MultiQC, fastp, Bowtie2,
samtools, skani, FastANI and the Python libraries into an environment called `synthcom`. It takes a
few minutes.

If you created the environment before and `envs/environment.yml` has changed
since, update it instead. This is needed after the ANI and abundance steps
were added: they brought in FastANI and plotnine (the plotting library), and
the mapping and plotting steps fail with "command not found" or
"No module named plotnine" until the environment is updated:

```bash
conda env update -f envs/environment.yml --prune
```

The server also has a system-wide Snakemake, but it cannot submit SLURM jobs,
so always use the one from this environment.

## 2. Every time you log in

```bash
cd /workdir/$USER/synthetic-community-pipeline
conda activate synthcom
snakemake --version        # should print 8.6 or higher
```

## 3. Configure

Open `config/config.yaml` and check the paths:

- the folder with the FASTQ files
- the metadata spreadsheet and the name of its tab
- the five genome FASTA files
- the output folder

Two settings control the mapping (steps 6-7):

- `bowtie2: very_sensitive`: `false` runs Bowtie2 with its default preset,
  `true` adds `--very-sensitive`. The BAM files go to
  `results/mapping/sensitive/` or `results/mapping/very_sensitive/`, so you
  can run both and compare them.
- `map_samples`: `all`, or a list of sample names from `config/samples.tsv`
  to map only those (for a pilot). The QC steps always run on every sample.

The `resources` block at the end sets the cores, memory and time limit used
by every job.

The sample sheet, `config/samples.tsv`, is built for you the first time you
run `snakemake`: it pairs the FASTQ files and looks up each sample in the
metadata tab. If the reads or the metadata change, delete it and it is built
again on the next run:

```bash
rm config/samples.tsv
```

## 4. Run

Snakemake has to keep running until the last job is done, which can take hours.
Start it inside `tmux` or `screen`, so the run survives if your connection
drops. Use whichever you prefer; both are installed on the server.

With `tmux`:

```bash
tmux new -s synthcom       # start a session called synthcom
conda activate synthcom
# detach (leave it running): Ctrl+B, then D
# come back later:           tmux attach -t synthcom
```

Or with `screen`:

```bash
screen -S synthcom         # start a session called synthcom
conda activate synthcom
# detach (leave it running): Ctrl+A, then D
# come back later:           screen -r synthcom
```

SLURM runs the individual jobs; `tmux` or `screen` only keeps the Snakemake
window that submits them alive.

**Always do a dry run first.** It lists the jobs that would run without starting any of them:

```bash
snakemake --profile profiles/slurm -n
```

If the list looks right, run it:

```bash
snakemake --profile profiles/slurm
```

Snakemake submits each job to SLURM for you and waits for them. You do not
write any `sbatch` scripts.

To run only part of the pipeline, name the output you want. Examples: the
raw read QC report only, or the ANI tables only:

```bash
snakemake --profile profiles/slurm results/qc/raw/multiqc_report.html
snakemake --profile profiles/slurm results/ani/skani.tsv results/ani/fastani.tsv
```

The mapping steps use the trimmed reads from step 2, so they pick up a
finished QC run in the same `results` folder without redoing it.

## 5. Check on a run

```bash
squeue -u $USER                    # your jobs that are queued or running
ls logs/slurm/                     # one folder per pipeline step
less logs/slurm/<step>/<file>.out  # the log of one job
```

If a job fails, look at its log in `logs/slurm/<step>/`. Fix the cause and
run the same `snakemake` command again: finished jobs are not repeated.

To stop a run, press Ctrl+C once in the Snakemake window. It cancels the jobs
it submitted.

## 6. Results

| File or folder                             | What it is                                         |
|--------------------------------------------|----------------------------------------------------|
| `results/qc/raw/multiqc_report.html`       | QC of the raw reads, all samples in one report     |
| `results/qc/trimmed/multiqc_report.html`   | QC of the trimmed reads, all samples in one report |
| `results/fastp/`                           | one fastp report per sample                        |
| `results/fastp/trimmed_reads/`             | trimmed reads, used for mapping                    |
| `results/reference/combined.fasta`       | the 5 genomes in one file, contigs prefixed with the strain ID |
| `results/reference/combined.*.bt2`       | Bowtie2 index of that file                         |
| `results/ani/skani.tsv`                   | ANI matrix between the genomes (skani)             |
| `results/ani/fastani.tsv`, `.matrix`      | the same with FastANI, as pairs and as a matrix    |
| `results/mapping/<preset>/`               | one sorted, indexed BAM per mapped sample          |
| `logs/slurm/bowtie2_map/`                 | Bowtie2's alignment summary for each sample (overall alignment rate) |
| `results/coverage/<preset>/`              | `samtools coverage` output per sample, one line per contig |
| `results/abundance/<preset>/abundance.tsv` | per sample and strain: mapped reads, mean depth, breadth, relative abundance, with the metadata |
| `results/abundance/<preset>/composition_wells.png` | one stacked bar per well; rows = conditions, columns = passage days |
| `results/abundance/<preset>/composition_means.png` | mean per day as a line with one point per well, one panel per condition |
| `results/abundance/<preset>/composition_*_contig_average.png` | the same two figures from rel_abundance_contig_average |

In the ANI tables, genomes are named by their FASTA file path. skani writes a
second file, `skani.tsv.af`, with the fraction of each genome that aligned to
the other.

### How coverage and abundance are calculated

`samtools coverage` reports, for every contig, the number of mapped reads, the
number of bases covered by at least one read (`covbases`) and the mean depth
(`meandepth`). The genomes are in several contigs each, so `abundance.py` sums
the contigs of each strain. For every sample and strain:

- **mapped bases** = sum over the strain's contigs of `meandepth` x contig
  length. Multiplying by the length first means a long contig counts for more
  than a short one, which is what makes the sum correct.
- **mean_depth** = mapped bases / genome length. The mean depth over the whole
  genome; positions without any read count as zero, and each contig counts in
  proportion to its length. This is the same number as `bedtools genomecov`
  gives (as used by Autometa) and as AlleleFlux's
  `average_coverage_including_zeros` and `length_weighted_coverage`.
- **mean_depth_contig_average** = plain average of the contigs' `meandepth`
  values, so a 500 bp contig counts as much as a 400 kb one. Given for
  comparison.
- **breadth** = covered bases / genome length. The fraction of the genome
  with at least one read.
- **rel_abundance** = mean_depth / sum of mean_depth over the five strains.
- **rel_abundance_contig_average** = the same from mean_depth_contig_average.

The two depths answer different questions. mean_depth asks "how many read
bases landed on this genome, per base of genome?", which is what abundance
means. mean_depth_contig_average asks "what is the typical depth of a contig
of this genome?", which depends on how the assembly happens to be split up.
Example with a genome of two contigs, one of 1,000,000 bp at 10x depth and one
of 1,000 bp at 100x: mean_depth = (10 x 1,000,000 + 100 x 1,000) / 1,001,000 =
10.1x, while mean_depth_contig_average = (10 + 100) / 2 = 55x. In this data the
short contigs (often repeats or rRNA fragments) carry more depth than the long
ones, so mean_depth_contig_average is 1.3x higher for E. coli, 1.8x for
Enterococcus and 2-50x for the near-absent strains. Because the factor differs
between strains it changes the composition: in the pilot, Enterococcus moves
from about 75% to about 81% and E. coli from 25% to 18% with the same reads.
Both are in the table; the plots are made for both.

A strain with a low rel_abundance and a low breadth (a few percent of the
genome) is at the detection limit rather than present.

## Changing how jobs are submitted

`profiles/slurm/config.yaml` holds the `sbatch` command and the limits:

- `jobs`: how many jobs may be queued or running at once
- `default-resources`: memory and time for steps that do not set their own

Memory, time and threads for individual steps are set in `config/config.yaml`.
