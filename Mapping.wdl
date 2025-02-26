# 版本定义
version 1.0

task minimap2{
    input{
        File cyclone_fastq_fn
        String ref_dir
        String ref_prefix
        Int threads
        String sample
        String out_dir
    }
    String sam_file = "~{out_dir}/~{sample}_Cyclone.sam"
    String bam_file = "~{out_dir}/~{sample}_Cyclone.sort.bam"
    String ref_fn = "~{ref_dir}/~{ref_prefix}"
    command<<<
        start_time=$(date +%s)
        echo "Starting task at ${start_time}" > "~{out_dir}/minimap2.log"
        /app/mm2-fast/minimap2 -ax map-ont --eqx --secondary=no -t ~{threads} ~{ref_fn} ~{cyclone_fastq_fn} -o ~{sam_file} 
        end_time1=$(date +%s)
        echo "Ending minimap2 task at ${end_time1}" >> "~{out_dir}/minimap2.log"
        bash -c "samtools view -@ ~{threads} -bS ~{sam_file} | samtools sort -@ ~{threads} -o ~{bam_file}"
        samtools index -@ ~{threads} ~{bam_file}
        end_time=$(date +%s)
        echo "Ending bam sort at ${end_time}" >> "~{out_dir}/minimap2.log"
        Runtime1=$((end_time1 - start_time)) 
        Runtime2=$((end_time - end_time1))
        echo "minimap2 runtime: ${Runtime1}; sort bam runtime: ${Runtime2}" >> "~{out_dir}/minimap2.log"
    >>>

    output{
        File bam_fn = "~{out_dir}/~{sample}_Cyclone.sort.bam"
        File bam_idx_fn = "~{out_dir}/~{sample}_Cyclone.sort.bam.bai"
        File Cyc_log = "~{out_dir}/minimap2.log"
    }
}

task bwa_mem2{
    input{
        File DNB_read1_fn
        File DNB_read2_fn
        String ref_dir
        String ref_prefix
        Int threads
        String sample
        String out_dir
    }

    String sam_file = "~{out_dir}/~{sample}_DNB.sam"
    String bam_file = "~{out_dir}/~{sample}_DNB.sort.bam"
    String ref_fn = "~{ref_dir}/~{ref_prefix}"
    command<<<
        start_time=$(date +%s)
        echo "Starting task at ${start_time}" > "~{out_dir}/bwa.log"
        bwa-mem2 mem -t ~{threads} -m 30G -R "@RG\tID:foo_lane\tPL:BGI\tLB:library\tSM:~{sample}_pcr_free" ~{ref_fn} ~{DNB_read1_fn} ~{DNB_read2_fn} -o ~{sam_file}
        end_time1=$(date +%s)
        echo "Ending bwa at ${end_time1}" >> "~{out_dir}/bwa.log"
        bash -c "samtools view -@ ~{threads} -bS ~{sam_file} | samtools sort -@ ~{threads} -o ~{bam_file}"
        samtools index -@ ~{threads} ~{bam_file}
        end_time=$(date +%s)
        echo "Ending sort at ${end_time}" >> "~{out_dir}/bwa.log"
        Runtime1=$((end_time1 - start_time)) 
        Runtime2=$((end_time - end_time1))
        echo "bwa runtime: ${Runtime1}; sort bam runtime: ${Runtime2}" >> "~{out_dir}/bwa.log"
    >>>
    output{
        File bam_fn = "~{out_dir}/~{sample}_DNB.sort.bam"
        File bam_idx_fn = "~{out_dir}/~{sample}_DNB.sort.bam.bai"
        File DNB_log = "~{out_dir}/bwa.log"
    }
}

workflow Mapping{
    input{
        File cyclone_fastq_fn
        File DNB_read1_fn
        File DNB_read2_fn
        String ref_dir
        String ref_prefix
        Int threads
        String sample
        String out_dir
    }
    call minimap2{
        input:
            cyclone_fastq_fn = cyclone_fastq_fn,
            ref_dir = ref_dir,
            ref_prefix = ref_prefix,
            threads = threads/2,
            sample = sample,
            out_dir = out_dir
    }
    
    call bwa_mem2{
        input:
            DNB_read1_fn = DNB_read1_fn,
            DNB_read2_fn = DNB_read2_fn,
            ref_dir = ref_dir,
            ref_prefix = ref_prefix,
            threads = threads/2,
            sample = sample,
            out_dir = out_dir
    }

    output{
        File Cyc_bam_fn = minimap2.bam_fn
        File Cyc_bam_idx_fn = minimap2.bam_idx_fn
        File DNB_bam_fn = bwa_mem2.bam_fn
        File DNB_bam_idx_fn = bwa_mem2.bam_idx_fn
        File Cyc_log = minimap2.Cyc_log
        File DNB_log = bwa_mem2.DNB_log
    }
}