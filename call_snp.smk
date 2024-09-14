# Define configuration
SAMPLES = ["sample1", "sample2"]  # 根据实际样本进行扩展
REF_GENOME = "data/reference/sorghum_reference.fa"
KNOWN_SNPS = "data/reference/known_snps.vcf"
KNOWN_INDELS = "data/reference/known_indels.vcf"
THREADS = 56  # 设置为 56 个线程，适配你的节点

rule all:
    input:
        expand("output/vcf/{sample}.g.vcf", sample=SAMPLES),
        "output/vcf/final_output.vcf"

# 参考基因组索引 (使用 bwa 和 samtools)
rule index_reference:
    input:
        REF_GENOME
    output:
        "data/reference/sorghum_reference.fa.bwt"
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa-mem2 index {input}"

rule create_ref_dict_and_fai:
    input:
        REF_GENOME
    output:
        "data/reference/sorghum_reference.dict",
        "data/reference/sorghum_reference.fa.fai"
    conda:
        "envs/gatk4.yaml"
    shell:
        """
        gatk CreateSequenceDictionary -R {input} -O {output[0]}
        samtools faidx {input}
        """

# Mapping reads 到参考基因组 (使用 bwa 和 samtools)
rule map_reads:
    input:
        ref="data/reference/sorghum_reference.fa",
        r1="data/samples/{sample}_R1.fastq.gz",
        r2="data/samples/{sample}_R2.fastq.gz"
    output:
        "output/bam/{sample}.bam"
    threads: THREADS
    conda:
        "envs/mapping.yaml"
    shell:
        """
        bwa-mem2 mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -Sb - > {output}
        """

rule sort_bam:
    input:
        "output/bam/{sample}.bam"
    output:
        "output/bam/{sample}_sorted.bam"
    threads: THREADS
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule mark_duplicates:
    input:
        "output/bam/{sample}_sorted.bam"
    output:
        bam="output/bam/{sample}_dedup.bam",
        metrics="output/bam/{sample}_metrics.txt"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}"

rule base_recalibrator:
    input:
        bam="output/bam/{sample}_dedup.bam",
        ref=REF_GENOME,
        known_snps=KNOWN_SNPS,
        known_indels=KNOWN_INDELS
    output:
        "output/bam/{sample}_recal_data.table"
    conda:
        "envs/gatk4.yaml"
    shell:
        """
        gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_snps} --known-sites {input.known_indels} -O {output}
        """

rule apply_bqsr:
    input:
        bam="output/bam/{sample}_dedup.bam",
        recal_data="output/bam/{sample}_recal_data.table",
        ref=REF_GENOME
    output:
        "output/bam/{sample}_recal.bam"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {input.recal_data} -O {output}"

rule haplotype_caller:
    input:
        bam="output/bam/{sample}_recal.bam",
        ref=REF_GENOME
    output:
        "output/vcf/{sample}.g.vcf"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output} -ERC GVCF"

rule combine_gvcfs:
    input:
        expand("output/vcf/{sample}.g.vcf", sample=SAMPLES),
        ref=REF_GENOME
    output:
        "output/vcf/combined.g.vcf"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk CombineGVCFs -R {input.ref} " + " ".join(f"-V {gvcf}" for gvcf in input[:-1]) + " -O {output}"

rule genotype_gvcfs:
    input:
        gvcf="output/vcf/combined.g.vcf",
        ref=REF_GENOME
    output:
        "output/vcf/final_output.vcf"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"
