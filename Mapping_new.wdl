# 版本定义
version 1.0

task mapping{
    input{
        File DNB_read1_fn
        File DNB_read2_fn
        File cyclone_fastq_fn
        File ref_dir
        String ref_name
        Int threads
        String sample
        String out_dir
    }
    String Cyc_sam_file = "~{out_dir}_~{sample}_Cyclone.sam"
    String Cyc_bam_file = "~{out_dir}_~{sample}_Cyclone.sort.bam"
    String DNB_sam_file = "~{out_dir}_~{sample}_DNB.sam"
    String DNB_bam_file = "~{out_dir}_~{sample}_DNB.sort.bam"
    Int minimap_threads = threads/2
    Int bwa_threads = threads/2
    String ref_fn = "~{ref_dir}/~{ref_name}"
    command<<<
        mkdir -p "~{out_dir}_bwa"  && \
        start_time=$(date +%s)  && \
        bwa-mem2 mem -t ~{bwa_threads} -m 30G -R "@RG\tID:foo_lane\tPL:BGI\tLB:library\tSM:~{sample}_pcr_free" ~{ref_fn} ~{DNB_read1_fn} ~{DNB_read2_fn} -o ~{DNB_sam_file}  && \
        end_time1=$(date +%s)  && \
        conda run -n SNV samtools view -@ ~{bwa_threads} -bS ~{DNB_sam_file} | samtools sort -@ ~{bwa_threads} -o ~{DNB_bam_file} && \
        samtools index -@ ~{bwa_threads} ~{DNB_bam_file}  && \
        end_time=$(date +%s)  && \
        Runtime1=$((end_time1 - start_time))  && \
        Runtime2=$((end_time - end_time1))   && \
        echo "bwa runtime: ${Runtime1}; sort bam runtime: ${Runtime2}" > "~{out_dir}_bwa/bwa.log"  & bwa=$!

        # mkdir -p "~{out_dir}_minimap2" && \
        # start_time=$(date +%s) && \
        # /app/mm2-fast/minimap2 -ax map-ont --eqx --secondary=no -t ~{minimap_threads} ~{ref_fn} ~{cyclone_fastq_fn} -o ~{Cyc_sam_file} && \
        # end_time1=$(date +%s) && \
        # conda run -n SNV samtools view -@ ~{minimap_threads} -bS ~{Cyc_sam_file} | samtools sort -@ ~{minimap_threads} -o ~{Cyc_bam_file} && \
        # samtools index -@ ~{minimap_threads} ~{Cyc_bam_file} && \
        # end_time=$(date +%s) && \
        # Runtime1=$((end_time1 - start_time)) && \
        # Runtime2=$((end_time - end_time1)) && \
        # echo "minimap2 runtime: ${Runtime1}; sort bam runtime: ${Runtime2}" > "~{out_dir}_minimap2/minimap2.log" & minimap2=$!
        
        wait $bwa
        # wait $minimap2
    >>>

    output{
        File Cyc_bam_fn = "~{out_dir}_~{sample}_Cyclone.sort.bam"
        File Cyc_bam_idx_fn = "~{out_dir}_~{sample}_Cyclone.sort.bam.bai"
        File Cyc_log = "~{out_dir}_minimap2/minimap2.log"
        File DNB_bam_fn = "~{out_dir}_~{sample}_DNB.sort.bam"
        File DNB_bam_idx_fn = "~{out_dir}_~{sample}_DNB.sort.bam.bai"
        File DNB_log = "~{out_dir}_bwa/bwa.log"
    }
}

workflow map_pipline{
    input{
        File cyclone_fastq_fn
        File DNB_read1_fn
        File DNB_read2_fn
        File ref_dir
        String ref_name
        Int threads
        String sample
        String out_dir
    }

    call mapping{
        input:
            DNB_read1_fn = DNB_read1_fn,
            DNB_read2_fn = DNB_read2_fn,
            cyclone_fastq_fn = cyclone_fastq_fn,
            ref_dir = ref_dir,
            ref_name = ref_name,
            threads = threads,
            sample = sample,
            out_dir = out_dir
    }

    output{
        File Cyc_bam_fn = mapping.Cyc_bam_fn
        File Cyc_bam_idx_fn = mapping.Cyc_bam_idx_fn
        File Cyc_log = mapping.Cyc_log
        File DNB_bam_fn = mapping.DNB_bam_fn
        File DNB_bam_idx_fn = mapping.DNB_bam_idx_fn
        File DNB_log = mapping.DNB_log
    }
}