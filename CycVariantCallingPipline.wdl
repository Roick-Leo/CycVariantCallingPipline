# 版本定义
version 1.0

task minimap2{
    input{
        File cyclone_fastq_fn
        File ref_fn
        Int threads
        String sample
        String out_dir
    }
    String bam_file = "~{out_dir}_~{sample}_Cyclone.bam"
    String sort_bam_file = "~{out_dir}_~{sample}_Cyclone.sort.bam"
    command<<<
        start_time=$(date +%s)
        echo "Starting task at ${start_time}" > "~{out_dir}/minimap2.log"  # 记录开始时间
        conda run -n SNV /app/mm2-fast/minimap2 -ax map-ont --eqx --secondary=no -t ~{threads} ~{ref_fn} ~{cyclone_fastq_fn} \
        | samtools view -@ ~{threads} -bS -o ~{bam_file}
        end_time1=$(date +%s)
        echo "Ending minimap2 task at ${end_time1}" >> "~{out_dir}/minimap2.log"
        conda run -n SNV samtools sort -@ ~{threads} ~{bam_file} -o ~{sort_bam_file}
        conda run -n SNV samtools index -@ ~{threads} ~{sort_bam_file}
        end_time=$(date +%s)  # 记录结束时间
        echo "Ending bam sort at ${end_time}" >> "~{out_dir}/minimap2.log"
        Runtime1=$((end_time1 - start_time)) 
        Runtime2=$((end_time - end_time1))  # 计算运行时长
        echo "minimap2 runtime: ${Runtime1}; sort bam runtime: ${Runtime2}" >> "~{out_dir}/minimap2.log"
    >>>

    output{
        File bam_fn = "~{out_dir}_~{sample}_Cyclone.sort.bam"
        File bam_idx_fn = "~{out_dir}_~{sample}_Cyclone.sort.bam.bai"
    }
}

task bwa_mem2{
    input{
        File DNB_read1_fn
        File DNB_read2_fn
        File ref_fn
        Int threads
        String sample
        String out_dir
    }
    String bam_file = "~{out_dir}_~{sample}_DNB.bam"
    String sort_bam_file = "~{out_dir}_~{sample}_DNB.sort.bam"
    command<<<
        start_time=$(date +%s)  # 记录开始时间
        echo "Starting task at ${start_time}" > "~{out_dir}/bwa.log"  # 记录开始时间
        conda run -n SNV bwa-mem2 mem -t ~{threads} -m 40G -R "@RG\tID:foo_lane\tPL:BGI\tLB:library\tSM:~{sample}_pcr_free" ~{ref_fn} ~{DNB_read1_fn} ~{DNB_read2_fn} | samtools view -@ ~{threads} -bS -o ~{bam_file}
        end_time1=$(date +%s)
        echo "Ending bwa at ${end_time1}" >> "~{out_dir}/bwa.log"
        conda run -n SNV samtools sort -@ ~{threads} ~{bam_file} -o ~{sort_bam_file}
        conda run -n SNV samtools index -@ ~{threads} ~{sort_bam_file}
        end_time=$(date +%s)  # 记录结束时间
        echo "Ending sort at ${end_time}" >> "~{out_dir}/bwa.log"
        Runtime1=$((end_time1 - start_time)) 
        Runtime2=$((end_time - end_time1))  # 计算运行时长
        echo "bwa runtime: ${Runtime1}; sort bam runtime: ${Runtime2}" >> "~{out_dir}/bwa.log"
    >>>
    output{
        File bam_fn = "~{out_dir}_~{sample}_DNB.sort.bam"
        File bam_idx_fn = "~{out_dir}_~{sample}_DNB.sort.bam.bai"
    }
}

task Cyc_sniffles{
    input{
        File cyclone_bam_fn
        File ref_fn
        String ref_version
    }
    command<<<
        conda run -n SV python /app/CycVariantCallingPipline/CycSV.py \
        
}

task GetCandidateBed{
    input{
        File input_fn
        File 
    }
    command<<<
        conda run -n SV python /app/CycVariantCallingPipline/SV/Generate_target_bed.py \
        -i 
}

task manta{
    input{
        File dnb_bam_fn
        File target_bed
        File ref_fn
        Int threads
        String sample
        String out_dir
    }
    command<<<
        TargetBam = "~{out_dir}_~{sample}.bam"
        samtools view -@ ~{threads} -Shb ~{dnb_bam_fn} -L ~{target_bed} -o ~{TargetBam}
        TargetSortBam = "~{out_dir}_~{sample}.sort.bam"
        samtools sort -@ ~{threads} ~{TargetBam} -o ~{TargetSortBam}
        samtools index -@ ~{threads} ~{TargetSortBam}
        mkdir -p ~{out_dir}
        conda run -n Manta python2 /app/manta-1.6.0.centos6_x86_64/bin/configManta.py \
        --bam ~{TargetSortBam} \
        --referenceFasta ~{ref_fn} \
        --runDir ~{out_dir} && python2 "~{out_dir}/runWorkflow.py" -m local -j ~{threads}
    >>>
    output{
        File out_vcf = ~{out_dir}/results/variants/diploidSV.vcf.gz
    }
}


task freesia{
    input{
        File dnb_bam_fn
        File cyclone_bam_fn
        File ref_fn
        Int threads
        String sample
        String out_dir
    }
    String env_dir = "/opt/conda/envs/Freesia/bin"
    String model_path = "/tools/models"
    String longreads_region = "/tools/data/Longreads.bed"
    command<<<
        samtools index -@ ~{threads} ~{dnb_bam_fn}
        samtools index -@ ~{threads} ~{cyclone_bam_fn} 
        samtools faidx ~{ref_fn} 
        /tools/Freesia.bin \
        --env_dir ~{env_dir} \
        --dnb_bam_fn ~{dnb_bam_fn} \
        --cyclone_bam_fn ~{cyclone_bam_fn} \
        --ref_fn ~{ref_fn} \
        --model_path ~{model_path} \
        --longreads_region ~{longreads_region} \
        --threads ~{threads} \
        --sample ~{sample} \
        --out_dir ~{out_dir}
    >>>

    output{
        File out_path = "~{out_dir}"
    }

    runtime{
        docker_url:"stereonote_hpc_external/helei1_d817bfeb3a2c444ea7dd384fa89f3613_private:latest"
        req_cpu:"~{threads}"
        req_memory:"90Gi"
    }
}

workflow Freesia{
    input{
        File dnb_bam_fn
        File cyclone_bam_fn
        File ref_fn
        Int threads
        String sample
        String out_dir
    }
    
    call freesia{
        input:
            dnb_bam_fn = dnb_bam_fn,
            cyclone_bam_fn = cyclone_bam_fn,
            ref_fn = ref_fn,
            threads = threads,
            sample = sample,
            out_dir = out_dir
    }
    
    output{
        File out_path = freesia.out_path
    }
}