# Stage 01: merge lanes, raw FastQC + MultiQC, overrepresented sequences.
#
#   snakemake -s workflow/01_raw_qc.smk --sdm conda -c 16
#
# INSPECT
#   results/qc/raw_multiqc/multiqc_report.html
#   results/qc/adapters/<sample>_overrepresented.tsv   (count, %, possible source)
# DECIDE
#   decisions/adapters/<sample>.fa  - sequences fastp should trim as adapters.
#   Pre-filled with the adapter-like hits only; add or remove sequences.
#   An empty file is fine (fastp then auto-detects adapters).
# Also check fastp_args in config/config.yaml before stage 02.

import glob
import os
import shutil

import pandas as pd

configfile: "config/config.yaml"

# ---- samples ---------------------------------------------------------------
# config/samples.tsv has three columns: sample  ID  species
samples_df = pd.read_csv(config["samples_tsv"], sep=r"\s+", engine="python", dtype=str)
SAMPLES = samples_df["sample"].tolist()

READS_DIR = config["reads_dir"]          # data/illumina
EXPECTED_LANES = config.get("expected_lanes")   # 2


def lane_files(sample, read):
    """All lane files for one sample and read (R1 or R2), sorted L001, L002, ...
    Stops with an error if files are missing or a lane is absent."""
    pattern = f"{READS_DIR}/{sample}_S*_L*_{read}_001.fastq.gz"
    files = sorted(glob.glob(pattern))
    if not files:
        raise ValueError(f"No files match {pattern}")
    if EXPECTED_LANES and len(files) != int(EXPECTED_LANES):
        raise ValueError(f"{sample} {read}: found {len(files)} lane(s), expected {EXPECTED_LANES}: {files}")
    return files


wildcard_constraints:
    sample="|".join(SAMPLES),
    read="R1|R2",


rule all:
    input:
        "results/qc/raw_multiqc/multiqc_report.html",
        expand("results/qc/adapters/{sample}_overrepresented.tsv", sample=SAMPLES),
    default_target: True


rule merge_lanes:
    """Concatenate all lanes into one R1 and one R2 per sample."""
    input:
        r1=lambda wc: lane_files(wc.sample, "R1"),
        r2=lambda wc: lane_files(wc.sample, "R2"),
    output:
        r1="data/merged/{sample}_R1.fastq.gz",
        r2="data/merged/{sample}_R2.fastq.gz",
    log:
        "logs/merge_lanes/{sample}.log",
    shell:
        "(cat {input.r1} > {output.r1} && cat {input.r2} > {output.r2}) 2> {log}"


rule fastqc_raw:
    input:
        "data/merged/{sample}_{read}.fastq.gz",
    output:
        zip="results/qc/raw_fastqc/{sample}_{read}_fastqc.zip",
        html="results/qc/raw_fastqc/{sample}_{read}_fastqc.html",
    log:
        "logs/fastqc_raw/{sample}_{read}.log",
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        "fastqc -t {threads} --outdir results/qc/raw_fastqc {input} &> {log}"


rule multiqc_raw:
    input:
        expand("results/qc/raw_fastqc/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=["R1", "R2"]),
    output:
        "results/qc/raw_multiqc/multiqc_report.html",
    log:
        "logs/multiqc/raw.log",
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc -f -o results/qc/raw_multiqc results/qc/raw_fastqc &> {log}"


rule overrepresented_sequences:
    """FastQC overrepresented sequences for R1+R2: a review table (count, %,
    possible source) plus a FASTA of the ones that look like adapters."""
    input:
        r1="results/qc/raw_fastqc/{sample}_R1_fastqc.zip",
        r2="results/qc/raw_fastqc/{sample}_R2_fastqc.zip",
    output:
        tsv="results/qc/adapters/{sample}_overrepresented.tsv",
        all_fa="results/qc/adapters/{sample}_overrepresented.fa",
        adapter_fa="results/qc/adapters/{sample}_adapter_like.fa",
    log:
        "logs/overrepresented/{sample}.log",
    params:
        regex=config["adapter_source_regex"],
    conda:
        "envs/utils.yaml"
    script:
        "scripts/overrepresented.py"


# When everything has finished: copy the adapter-like sequences into
# decisions/adapters/<sample>.fa for you to review. Existing files are never
# overwritten, so your edits are safe if you rerun this stage.
onsuccess:
    os.makedirs("decisions/adapters", exist_ok=True)
    for s in SAMPLES:
        src = f"results/qc/adapters/{s}_adapter_like.fa"
        dst = f"decisions/adapters/{s}.fa"
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.copyfile(src, dst)
            print(f"Created {dst}")
