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
        String Cyc_dir = out_dir
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
        String DNB_dir = out_dir
    }
}

task variant_calling{
    input{
        String dnb_bam_dir
        String dnb_bam_name
        String cyclone_bam_dir
        String cyclone_bam_name
        String ref_dir
        String ref_prefix
        String ref_version
        Int threads
        String prefix
        String out_dir
    }
    String ref_fn = "~{ref_dir}/~{ref_prefix}"
    String env_dir = "/app/conda/envs/SNV/bin"
    String model_path = "/app/freesia_v1/SNV/models"
    String dnb_bam_fn = "~{dnb_bam_dir}/~{dnb_bam_name}"
    String cyclone_bam_fn = "~{cyclone_bam_dir}/~{cyclone_bam_name}"

    Boolean valid_version = (ref_version == "38" || ref_version == "37")
    String longreads_region = if (ref_version == "38") then "/app/freesia_v1/SNV/data/GRCh38_no_alt.longreads.bed" else "/app/freesia_v1/SNV/data/hs37d5.longreads.bed"
    String tandem_repeats = if (ref_version == "38") then "/app/freesia_v1/SV/data/GRCh38_no_alt.trf.bed" else "/app/freesia_v1/SV/data/hs37d5.trf.bed"

    command<<<
        if [[ ${valid_version} == "false" ]]; then
            echo "ERROR: Unsupported ref_version. Please use '37' or '38'." >&2
            exit 1
        fi

        start_time=$(date +%s)
        echo "Starting SNV calling at ${start_time}" > "~{out_dir}/variant_calling.log"
        conda run -n SNV /app/freesia_v1/SNV/Freesia.bin \
            --env_dir ~{env_dir} \
            --dnb_bam_fn ~{dnb_bam_fn} \
            --cyclone_bam_fn ~{cyclone_bam_fn} \
            --ref_fn ~{ref_fn} \
            --model_path ~{model_path} \
            --longreads_region ~{longreads_region} \
            --threads ~{threads/2} \
            --sample "~{prefix}_Freesia_SNV" \
            --out_dir ~{out_dir}
        end_time1=$(date +%s)
        echo "Ending SNV calling at ${end_time1}" >> "~{out_dir}/variant_calling.log"
        
        export LD_LIBRARY_PATH=/app/conda/envs/SV/lib:$LD_LIBRARY_PATH
        echo "Starting SV calling at ${end_time1}" > "~{out_dir}/variant_calling.log"
        conda run -n SV /app/freesia_v1/SV/src/Get_candidate.bin \
            -i ~{cyclone_bam_fn} \
            -o ~{out_dir} \
            -p "~{prefix}_Freesia_Candidate_SV" \
            -r ~{ref_fn} \
            --tandem-repeats ~{tandem_repeats} \
            -t ~{threads}
        
        conda run -n PY2 /app/freesia_v1/SV/src/Get_intergration.bin \
            --input_bam ~{dnb_bam_fn} \
            --out_dir ~{out_dir} \
            --prefix ~{prefix} \
            --ref_file ~{ref_fn} \
            --threads ~{threads} \
            --candidate_vcf "~{out_dir}/~{prefix}_Freesia_Candidate_SV_Candidate.vcf" \
            --extension 1500 \
            -hq "~{out_dir}/~{prefix}_Freesia_Candidate_SV_HighQ.vcf"
        end_time=$(date +%s)
        echo "Ending SV calling at ${end_time}" >> "~{out_dir}/variant_calling.log"
    >>>
    output{
        File snv_vcf_file = "~{out_dir}/~{prefix}_Freesia_SNV.vcf.gz"
        File snv_vcf_idx = "~{out_dir}/~{prefix}_Freesia_SNV.vcf.gz.tbi"
        File sv_vcf_file = "~{out_dir}/~{prefix}_Freesia_SV.vcf.gz"
        File sv_vcf_idx = "~{out_dir}/~{prefix}_Freesia_SV.vcf.gz.tbi"
        File Variant_Calling_log = "~{out_dir}/variant_calling.log"
    }
}


workflow freesia{
    input{
        File cyclone_fastq_fn
        File DNB_read1_fn
        File DNB_read2_fn
        String ref_dir
        String ref_prefix
        String ref_version
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

    call variant_calling{
        input:
            dnb_bam_dir = bwa_mem2.DNB_dir,
            dnb_bam_name = "~{sample}_DNB.sort.bam",
            cyclone_bam_dir = minimap2.Cyc_dir,
            cyclone_bam_name = "~{sample}_Cyclone.sort.bam",
            ref_dir = ref_dir,
            ref_prefix = ref_prefix,
            ref_version = ref_version,
            threads = threads,
            prefix = sample,
            out_dir = out_dir
    }

    output{
        File Cyc_bam_fn = minimap2.bam_fn
        File Cyc_bam_idx_fn = minimap2.bam_idx_fn
        File DNB_bam_fn = bwa_mem2.bam_fn
        File DNB_bam_idx_fn = bwa_mem2.bam_idx_fn
        File Cyc_log = minimap2.Cyc_log
        File DNB_log = bwa_mem2.DNB_log
        File snv_vcf_file = variant_calling.snv_vcf_file
        File snv_vcf_idx = variant_calling.snv_vcf_idx
        File sv_vcf_file = variant_calling.sv_vcf_file
        File sv_vcf_idx = variant_calling.sv_vcf_idx
        File Variant_Calling_log = variant_calling.Variant_Calling_log
    }
}