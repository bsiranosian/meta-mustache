configfile: "config.yaml"
from os.path import join
WC_fastqs = glob_wildcards(join(config['fastq_dir'], '{SRR}_1.fastq'))
SAMPLES = set(WC_fastqs.SRR)
# SAMPLES = ['SRR2239592']
super_deduper = config['super_deduper_exec']
blastdb_location = config['blastdb_location']

rule all:
    input:
        # expand('trimmed/{SRR}_nodup_PE1_val_1.fq',SRR=SAMPLES)
        # expand('aligned/{SRR}.bam', SRR=SAMPLES)
        expand("find_results/{SRR}/{SRR}.insertions_seqs_blast_results.tsv", SRR=SAMPLES)
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
        "module load trim_galore && module load cutadapt && trim_galore --fastqc --paired {input} -o trimmed"

rule metaqc:
    input:
        expand('trimmed/{SRR}_nodup_PE1_val_1.fq',SRR=SAMPLES)
    output: 
        "trimmed/multiqc_report.html"
    shell: 
        "module load multiqc && multiqc trimmed "

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
        "aligned/{SRR}.bam",
        config['ref_genome']
    output:
        "find_results/{SRR}/{SRR}.insertions_seqs.tsv"
    shell:
        "mustache find paired {input} {wildcards.SRR} --outdir find_results/{wildcards.SRR}"

rule prepare_blast:
    input:
        "find_results/{SRR}/{SRR}.insertions_seqs.tsv"
    output:
        "find_results/{SRR}/{SRR}.insertions_seqs_blast.fa"
    shell:
        """
        paste -d '\n' <(sed -e 's/^/>/' <(grep "M" {input} | cut -f 1,2,3,4,5,6,7,8 | tr "\t" "_"))  <(grep "M" {input} | cut -f 9) > {output}
        """
rule run_blast:
    input:
        "find_results/{SRR}/{SRR}.insertions_seqs_blast.fa"
    output:
        "find_results/{SRR}/{SRR}.insertions_seqs_blast_results.tsv"
    shell:
        """
        module load blast
        export BLASTDB="/labs/asbhatt/gsfs0/bsiranos/blast_db/nr/"
        export BLASTDB="/labs/asbhatt/gsfs0/bsiranos/blast_db/tsa_nr:$BLASTDB"
        export BLASTDB="/labs/asbhatt/gsfs0/bsiranos/blast_db/env_nr:$BLASTDB"
        echo -e 'qseqid\tsseqid\tstaxids\tsscinames\tscomnames\tsblastnames\tsalltitles\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' > {output}
        blastx -query {input} -db "nr tsa_nr env_nr" -outfmt '6 qseqid sseqid staxids sscinames scomnames sblastnames salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore' >> {output}
        """