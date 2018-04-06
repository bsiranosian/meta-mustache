configfile: "config.yaml"
from os.path import join
WC_fastqs = glob_wildcards(join(config['fastq_dir'], '{SRR}_1.fastq'))
SAMPLES = set(WC_fastqs.SRR)
super_deduper = config['super_deduper_exec']

rule all:
    input:
        expand('trimmed/{SRR}_nodup_PE1_val_1.fq',SRR=SAMPLES)
    run:
        print("Finished through alignment step!")


rule deduplicate_reads:
    input:
        r1="srr/{SRR}_1.fastq",
        r2="srr/{SRR}_2.fastq"
    output:
        "deduplicated/{SRR}_nodup_PE1.fastq",
        "deduplicated/{SRR}_nodup_PE2.fastq"
    shell:
        "{super_deduper} -1 {input.r1} -2 {input.r2} -p deduplicated/{wildcards.SRR}"

rule trim_galore:
    input: 
        "deduplicated/{SRR}_nodup_PE1.fastq",
        "deduplicated/{SRR}_nodup_PE2.fastq"
    output:
        "trimmed/{SRR}_nodup_PE1_val_1.fq",
        "trimmed/{SRR}_nodup_PE2_val_2.fq"
    shell:
        "module load trim_galore"
        "module load cutadapt"
        "trim_galore --fastqc --paired {input} -o trimmed"


rule mustache_align:
    input:
        "trimmed/{SRR}_nodup_PE1_val_1.fq",
        "trimmed/{SRR}_nodup_PE2_val_2.fq",
        config['ref_genome']
    output:
        "aligned/{SRR}.bam"
    shell:
        "mustache align paired {input} {output}"

rule mustache_find:
    input:
        "aligned/{SRR}.bam"
    output:
        "find_results/{SRR}/{SRR}.stats.tsv"
    shell:
        "mustache find paired {input} find_results/{wildcards.SRR}"