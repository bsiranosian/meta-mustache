configfile: "config.yaml"

rule deduplicate_reads:
    input:
        r1="srr/{SRR}_1.fastq",
        r2="srr/{SRR}_2.fastq"
    output:
        "deduplicated/{SRR}_nodup_PE1.fastq",
        "deduplicated/{SRR}_nodup_PE2.fastq"
    shell:
        "~/software/Super-Deduper/super_deduper -1 {input.r1} -2 {input.r2} -p deduplicated/{wildcards.SRR}"

rule trim_galore:
    input: 
        "deduplicated/{SRR}_nodup_PE1.fastq",
        "deduplicated/{SRR}_nodup_PE2.fastq"
    output:
        "trimmed/{SRR}_nodup_PE1_val_1.fq",
        "trimmed/{SRR}_nodup_PE2_val_2.fq"
    shell:
        "trim_galore --fastqc --paired {input} -o trimmed"


rule mustache_align:
    input:
        "trimmed/{SRR}_nodup_PE1_val_1.fq",
        "trimmed/{SRR}_nodup_PE2_val_2.fq"
        "/home/bsiranos/bs_scratch/genomes/GCF_000005845.2_ASM584v2_genomic.fna"
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