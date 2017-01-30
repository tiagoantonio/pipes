rule ubam_fastq2bam:
    input:
        # lambda wildcards: config["samples"][wildcards.sample]
        fastq="{sample}.R1.fastq.gz",
        fastq2="{sample}.R2.fastq.gz"
    params:
        samplename="{sample}",
        libname="hd_nati",
        seqcenter="CEGH RUN_DATE=2016-10-20T00:00:00-0400"
    output: 
        temp("mapped_reads/{sample}_fastqtosam.bam")
    log:
        "logs/ubam_fastq2bam/{sample}.log"
    benchmark:
        "benchmarks/{sample}.ubam_fastq2bam.benchmark.txt"
    shell:
        "~/tools/jdk1.8.0_111/bin/java -Xmx24G -XX:ParallelGCThreads=20 -jar ~/tools/picard.jar FastqToSam FASTQ={input.fastq} "
        "FASTQ2={input.fastq2} OUTPUT={output} READ_GROUP_NAME=H0164.2 SAMPLE_NAME={params.samplename} LIBRARY_NAME={params.libname} "
        "PLATFORM_UNIT=H0164ALXX140820.2 PLATFORM=illumina SEQUENCING_CENTER={params.seqcenter}"

rule ubam_revertsam:
    input:
        revertsam="mapped_reads/{sample}_fastqtosam.bam",
    output: 
        temp("mapped_reads/{sample}_fastqtosam.revertsam.bam")
    log:
        "logs/ubam_revertsam/{sample}.log"
    benchmark:
        "benchmarks/{sample}.ubam_revertsam.benchmark.txt"
    shell:
        "~/tools/jdk1.8.0_111/bin/java -Xmx24G -XX:ParallelGCThreads=20 -jar ~/tools/picard.jar RevertSam I={input.revertsam} "
        "O={output} SANITIZE=true MAX_DISCARD_FRACTION=0.005 ATTRIBUTE_TO_CLEAR=XT ATTRIBUTE_TO_CLEAR=XN "
        "ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=true "
        "REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true"

rule ubam_markadapters:
    input:
        markadapters="mapped_reads/{sample}_fastqtosam.revertsam.bam"
    output: 
        temp("mapped_reads/{sample}_fastqtosam.revertsam.markadapters.bam")
    params:
        ubammetrics="{sample}_fastqtosam.revertsam.markadapters_metrics.txt"
    log:
        "logs/ubam_markadapters/{sample}.log"
    benchmark:
        "benchmarks/{sample}.ubam_markadapters.benchmark.txt"
    shell:
        "~/tools/jdk1.8.0_111/bin/java -Xmx24G -XX:ParallelGCThreads=20 -jar ~/tools/picard.jar MarkIlluminaAdapters "
        "I={input.markadapters} "
        "O={output} "
        "M={params.ubammetrics}"

rule mapping_bam2fastq:
    input:
        bammarkadapters="mapped_reads/{sample}_fastqtosam.revertsam.markadapters.bam"
    output: 
        temp("mapped_reads/{sample}_samtofastq_interleaved.fq")
    log:
        "logs/mapping_bam2fastq/{sample}.log"
    benchmark:
        "benchmarks/{sample}.mapping_bam2fastq.benchmark.txt"
    shell:
        "~/tools/jdk1.8.0_111/bin/java -Xmx24G -XX:ParallelGCThreads=20 -jar ~/tools/picard.jar SamToFastq I={input.bammarkadapters} FASTQ={output} "
        "CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true"

rule mapping_bwa:
    input:
        samfastq="mapped_reads/{sample}_samtofastq_interleaved.fq"
    output: 
        temp("mapped_reads/{sample}_bwaaligned.sam")
    log:
        "logs/mapping_bwa/{sample}.log"
    benchmark:
        "benchmarks/{sample}.mapping_bwa.txt"
    shell:
        "bwa mem -M -t 20 -p ~/reference/hg19-gatk/ucsc.hg19.fasta {input.samfastq} > {output}"

rule mapping_mergebam:
    input:
        bwamapped="mapped_reads/{sample}_bwaaligned.sam",
        unmapbam="mapped_reads/{sample}_fastqtosam.revertsam.markadapters.bam"
    output: 
        temp("mapped_reads/{sample}_merged.bam")
    log:
        "logs/mapping_mergebam/{sample}.log"
    benchmark:
        "benchmarks/{sample}.mapping_mergebam.txt"
    shell:
        "~/tools/jdk1.8.0_111/bin/java -Xmx24G -XX:ParallelGCThreads=20 -jar ~/tools/picard.jar MergeBamAlignment R=~/reference/hg19-gatk/ucsc.hg19.fasta "
        "UNMAPPED_BAM={input.unmapbam} ALIGNED_BAM={input.bwamapped} O={output} CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false "
        "CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant "
        "ATTRIBUTES_TO_RETAIN=XS"

rule mapping_mdup:
    input:
        mergedbam="mapped_reads/{sample}_merged.bam"
    output: 
        temp("mapped_reads/{sample}_mergebamalignment_markduplicates.bam")
    params:
        mdupmetrics="{sample}_mdupmetrics.txt"
    log:
        "logs/mapping_mdup/{sample}.log"
    benchmark:
        "benchmarks/{sample}.mapping_mdup.txt"
    shell:
        "~/tools/jdk1.8.0_111/bin/java -Xmx24G -XX:ParallelGCThreads=20 -jar ~/tools/picard.jar MarkDuplicates INPUT={input.mergedbam} OUTPUT={output} "
        "METRICS_FILE={params.mdupmetrics} OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true"

rule realign_knowintervals:
    input:
        mapdup="mapped_reads/{sample}_mergebamalignment_markduplicates.bam"
    output: 
        temp("mapped_reads/{sample}_merge_mduplicates_realigned.bam")
    log:
        "logs/realign_knowintervals/{sample}.log"
    benchmark:
        "benchmarks/{sample}.realign_knowintervals.txt"
    shell:
        "java -Xmx24G -Djava.io.tmpdir=tmpdir -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -I {input.mapdup} -R ~/reference/hg19-gatk/ucsc.hg19.fasta -T IndelRealigner -targetIntervals output.intervals -o {output} -known ~/reference/hg19-gatk/1000G_phase1.indels.hg19.sites.vcf --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4"

rule bqsr_firststep:
    input:
        realigned="mapped_reads/{sample}_merge_mduplicates_realigned.bam"
    output: 
        temp("mapped_reads/{sample}_recal_data_table")
    log:
        "logs/bqsr_firststep/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bqsr_firststep.txt"
    shell:
        "java -Xmx24G -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T BaseRecalibrator -R ~/reference/hg19-gatk/ucsc.hg19.fasta -I {input.realigned} -knownSites ~/reference/hg19-gatk/dbsnp_138.hg19.vcf -knownSites ~/reference/hg19-gatk/1000G_phase1.indels.hg19.sites.vcf  -o {output}",

rule bqsr_apply:
    input:
        realigned="mapped_reads/{sample}_merge_mduplicates_realigned.bam",
        recal_data="mapped_reads/{sample}_recal_data_table"
    output: 
        "mapped_reads/{sample}_merged_mdup_realigned_recal.bam"
    log:
        "logs/bqsr_apply/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bqsr_apply.txt"
    shell:
        "java -Xmx24G -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T PrintReads -R  ~/reference/hg19-gatk/ucsc.hg19.fasta -I {input.realigned} -BQSR {input.recal_data} -o {output}"

rule cleanup_bai:
    shell:
        "rm mapped_reads/"
        "{sample}_mergebamalignment_markduplicates.bai "
        "{sample}_merged.bai "
        "{sample}_merge_mduplicates_realigned.bai"

rule collect_enrich_metrics:
    input:
        recalbam="mapped_reads/{sample}_merged_mdup_realigned_recal.bam"
    output:
        "mapped_reads/{sample}_enrichment_stats.txt"
    log:
        "logs/collect_enrich_metrics/{sample}.log"
    benchmark:
        "benchmarks/{sample}.collect_enrich_metrics.txt"
    shell:
        "~/tools/jdk1.8.0_111/bin/java -Xmx24G -XX:ParallelGCThreads=20 -jar ~/tools/picard.jar CollectHsMetrics I={input} O={output} R=~/reference/hg19-gatk/ucsc.hg19.fasta BAIT_INTERVALS=~/reference/hg19-gatk/bed/sureselectV6/S07604514_list.interval_list TARGET_INTERVALS=~/reference/hg19-gatk/bed/sureselectV6/S07604514_list.interval_list"

rule snpcalling_haplotypecaller:
    input:
        recalbam="mapped_reads/{sample}_merged_mdup_realigned_recal.bam"
    output:
        "snp_calling/{sample}_raw_variants.vcf"
    log:
        "logs/snpcalling_haplotypecaller/{sample}.log"
    benchmark:
        "benchmarks/{sample}.snpcalling_haplotypecaller.txt"
    shell:
        "java -Xmx24G -XX:ParallelGCThreads=20 -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T HaplotypeCaller "
        "-R ~/reference/hg19-gatk/ucsc.hg19.fasta -I {input.recalbam} --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o {output}"

rule snpcalling_collectmetrics:
    input: 
        rawsnps="snp_calling/{sample}_raw_variants.vcf"
    output:
        "snp_calling/{sample}_snpmetrics"
    log:
        "logs/snpcalling_collectmetrics/{sample}.log"
    benchmark:
        "benchmarks/{sample}.snpcalling_collectmetrics.txt"
    shell:
        "~/tools/jdk1.8.0_111/bin/java -Xmx24G -XX:ParallelGCThreads=20 -jar ~/tools/picard.jar CollectVariantCallingMetrics INPUT={input.rawsnps} "
        "OUTPUT={output} DBSNP=~/reference/hg19-gatk/dbsnp_138.hg19.vcf TARGET_INTERVALS=/home/tantonio/reference/hg19-gatk/bed/sureselectV6/S07604514_list.interval_list "
        "SEQUENCE_DICTIONARY=~/reference/hg19-gatk/ucsc.hg19.dict"

rule filtering_gatkonlysnps:
    input:
        rawsnps="snp_calling/{sample}_raw_variants.vcf"
    output:
        "snp_filtered/{sample}_onlysnps.vcf"
    log:
        "logs/filtering_gatkonlysnps/{sample}.log"
    benchmark:
        "benchmarks/{sample}.filtering_gatkonlysnps.txt"
    shell:
        "java -Xmx24G -XX:ParallelGCThreads=20 -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T SelectVariants -R ~/reference/hg19-gatk/ucsc.hg19.fasta -V {input.rawsnps} -selectType SNP -o {output}"

rule filtering_gatkfilter:
    input:
        onlysnps="snp_filtered/{sample}_onlysnps.vcf"
    output:
        "snp_filtered/{sample}_passed.vcf"
    log:
        "logs/filtering_gatkfilter/{sample}.log"
    benchmark:
        "benchmarks/{sample}.filtering_gatkfilter.txt"
    shell:
        "java -Xmx24G -XX:ParallelGCThreads=20 -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T VariantFiltration -R ~/reference/hg19-gatk/ucsc.hg19.fasta -V {input.onlysnps} --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" -L ~/reference/hg19-gatk/bed/sureselectV6/S07604514_list.interval_list --filterName \"my_snp_filter\" -o {output}"

rule filtering_vcftools:
    input:
        passed="snp_filtered/{sample}_passed.vcf"
    output:
        "snp_filtered/{sample}_filteredgz.vcf.gz"
    log:
        "logs/filtering_vcftools/{sample}.log"
    benchmark:
        "benchmarks/{sample}.filtering_vcftools.txt"
    shell:
        "~/tools/vcftools_0.1.13/bin/vcftools --vcf {input.passed} --remove-filtered-all --max-alleles 2 --recode --stdout | ~/tools/tabix-0.2.6/bgzip -c > {output}"

rule filtering_tabix:
    input:
        filteredgz="snp_filtered/{sample}_filteredgz.vcf.gz"
    output:
        "snp_filtered/{sample}_filteredgz.vcf.gz.tbi"
    log:
        "logs/filtering_tabix/{sample}.log"
    benchmark:
        "benchmarks/{sample}.filtering_tabix.txt"
    shell:
        "~/tools/tabix-0.2.6/tabix {input.filteredgz}"

rule select_allindels:
    input:
        rawsnps="snp_calling/{sample}_raw_variants.vcf"
    output:
        "indels/{sample}_rawindels.vcf"
    log:
        "logs/select_allindels/{sample}.log"
    benchmark:
        "benchmarks/{sample}.select_allindels.txt"
    shell:
        "java -Xmx24G -XX:ParallelGCThreads=20 -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T SelectVariants -R ~/reference/hg19-gatk/ucsc.hg19.fasta -V {input.rawsnps} -selectType INDEL -o {output}"

rule filter_allindels:
    input:
        selectallindels="indels/{sample}_rawindels.vcf"
    output:
        "indels/{sample}_filteredindels.vcf"
    log:
        "logs/filter_allindels/{sample}.log"
    benchmark:
        "benchmarks/{sample}.filter_allindels.txt"
    shell:
        "java -Xmx24G -XX:ParallelGCThreads=20 -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T VariantFiltration -R ~/reference/hg19-gatk/ucsc.hg19.fasta -V {input.selectallindels} --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" -L ~/reference/hg19-gatk/bed/sureselectV6/S07604514_list.interval_list --filterName \"my_snp_filter\" -o {output}"

rule filtering_vcftools_allindels:
    input:
        passedindelsall="indels/{sample}_filteredindels.vcf"
    output:
        "indels/{sample}_passed_allindels.vcf"
    log:
        "logs/filtering_vcftools_allindels/{sample}.log"
    benchmark:
        "benchmarks/{sample}.filtering_vcftools_allindels.txt"
    shell:
        "~/tools/vcftools_0.1.13/bin/vcftools --vcf {input.passedindelsall} --remove-filtered-all --max-alleles 2 --recode --out {output}"

rule hist_allindels:
    input:
        filteredallindels="indels/{sample}_passed_allindels.vcf"
    output:
        "indels/{sample}_hist_allindels"
    log:
        "logs/hist_allindels/{sample}.log"
    benchmark:
        "benchmarks/{sample}.hist_allindels.txt"
    shell:
    	"~/tools/vcftools_0.1.13/bin/vcftools --vcf {input.filteredallindels} --hist-indel-len --out {output}"

rule select_novelindels:
    input:
        rawsnps="snp_calling/{sample}_raw_variants.vcf"
    output:
        "indels/{sample}_rawindels_novel.vcf"
    log:
        "logs/select_allindels_novel/{sample}.log"
    benchmark:
        "benchmarks/{sample}.select_allindels_novel.txt"
    shell:
        "java -Xmx24G -XX:ParallelGCThreads=20 -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T SelectVariants -R ~/reference/hg19-gatk/ucsc.hg19.fasta -V {input.rawsnps} -selectType INDEL --discordance ~/reference/hg19-gatk/dbsnp_138.hg19.vcf -o {output}"

rule filter_novelindels:
    input:
        selectnovelindels="indels/{sample}_rawindels_novel.vcf"
    output:
        "indels/{sample}_filteredindels_novel.vcf"
    log:
        "logs/filter_novelindels/{sample}.log"
    benchmark:
        "benchmarks/{sample}.filter_novelindels.txt"
    shell:
        "java -Xmx24G -XX:ParallelGCThreads=20 -jar /usr/local/genome/GenomeAnalysisTK-2.6-5-gba531bd/GenomeAnalysisTK.jar -T VariantFiltration -R ~/reference/hg19-gatk/ucsc.hg19.fasta -V {input.selectnovelindels} --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" -L ~/reference/hg19-gatk/bed/sureselectV6/S07604514_list.interval_list --filterName \"my_snp_filter\" -o {output}"

rule filtering_vcftools_novelindels:
    input:
        passedindelsnovel="indels/{sample}_filteredindels_novel.vcf"
    output:
        "indels/{sample}_passed_novelindels.vcf"
    log:
        "logs/filtering_vcftools_novelindels/{sample}.log"
    benchmark:
        "benchmarks/{sample}.filtering_vcftools_novelindels.txt"
    shell:
        "~/tools/vcftools_0.1.13/bin/vcftools --vcf {input.passedindelsnovel} --remove-filtered-all --max-alleles 2 --recode --out {output}"

rule hist_novelindels:
    input:
        filterednovelindels="indels/{sample}_passed_novelindels.vcf.recode.vcf"
    output:
        "indels/{sample}_hist_novelindels"
    log:
        "logs/hist_novelindels/{sample}.log"
    benchmark:
        "benchmarks/{sample}.hist_novelindels.txt"
    shell:
    	"~/tools/vcftools_0.1.13/bin/vcftools --vcf {input.filterednovelindels} --hist-indel-len --out {output}"