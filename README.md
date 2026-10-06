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
| 4    | Combine the 5 genomes into one reference and build a Bowtie2 index  | not written |
| 5    | All-vs-all ANI between the 5 genomes (skani)                        | not written |
| 6    | Bowtie2 pilot on 12 samples: 2 presets x 2 MAPQ filters             | not written |
| 7    | Bowtie2 on all samples, sorted and indexed BAMs                     | not written |
| 8    | Coverage per genome, relative abundance, tables and plots           | not written |

The steps are being added one at a time. This table is updated as each one lands.

## Repository layout

```
config/config.yaml          all paths and settings (the only file you edit)
config/samples.tsv          one row per sample, built from the metadata spreadsheet
envs/environment.yml        the conda environment with every tool
profiles/slurm/config.yaml  tells Snakemake how to submit jobs to SLURM
workflow/Snakefile          loads the settings and the sample sheet
workflow/rules/qc.smk       steps 1-3: FastQC, fastp, MultiQC
workflow/scripts/           small Python scripts called by the pipeline
results/                    everything the pipeline produces
logs/                       one log per job
```

## 1. One-time setup

Do this once on the lab server (`cbsumoeller02`).

### 1.1 Get the code

```bash
cd /workdir/$USER
git clone https://github.com/Sidduppal/synthetic-community-pipeline.git
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
samtools and skani into an environment called `synthcom`. It takes a few minutes.

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

Then build the sample sheet. This pairs the FASTQ files and looks up each
sample in the metadata tab:

```bash
python workflow/scripts/make_samplesheet.py --config config/config.yaml
```

It writes `config/samples.tsv` and prints how many samples it found per day and
incubation type. Run it again only if the reads or the metadata change.

## 4. Run

Start a `screen` session first, so the run survives if your connection drops:

```bash
screen -S synthcom         # detach with Ctrl+A then D, come back with: screen -r synthcom
conda activate synthcom
```

**Always do a dry run first.** It lists the jobs that would run and touches nothing:

```bash
snakemake --profile profiles/slurm -n
```

If the list looks right, run it:

```bash
snakemake --profile profiles/slurm
```

Snakemake submits each job to SLURM for you and waits for them. You do not
write any `sbatch` scripts.

To run only part of the pipeline, name the output you want. Example, the raw
read QC report only:

```bash
snakemake --profile profiles/slurm results/qc/raw/multiqc_report.html
```

## 5. Check on a run

```bash
squeue -u $USER                    # your jobs that are queued or running
ls logs/slurm/                     # one folder per pipeline step
less logs/slurm/<step>/<file>.out  # the log of one job
```

If a job fails, Snakemake prints the name of its log file. Fix the cause and
run the same `snakemake` command again: finished jobs are not repeated.

To stop a run, press Ctrl+C once in the Snakemake window. It cancels the jobs
it submitted.

## 6. Results

| File or folder                             | What it is                                         |
|--------------------------------------------|----------------------------------------------------|
| `results/qc/raw/multiqc_report.html`       | QC of the raw reads, all samples in one report     |
| `results/qc/trimmed/multiqc_report.html`   | QC after trimming, with the fastp trimming numbers |
| `results/qc/fastp/`                        | one fastp report per sample                        |
| `results/reads/trimmed/`                   | trimmed reads, used for mapping                    |

More rows are added as each step is written.

## Changing how jobs are submitted

`profiles/slurm/config.yaml` holds the `sbatch` command and the limits:

- `jobs`: how many jobs may be queued or running at once
- `default-resources`: memory and time for steps that do not set their own

Memory, time and threads for individual steps are set in `config/config.yaml`.
