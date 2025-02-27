# 版本定义
version 1.0

task SNV_calling{
    input{
        File dnb_bam_fn
        File cyclone_bam_fn
        String ref_dir
        String ref_prefix
        String ref_version
        Int threads
        String prefix
        String out_dir
    }
    String ref_fn = "~{ref_dir}/~{ref_prefix}"
    String env_dir = "/app/conda/envs/SNV/bin"
    String model_path = "/app/CycVariantCallingPipline/SNV/models"
    command<<<
        if (ref_version == "38") {
            longreads_region = "/app/CycVariantCallingPipline/SNV/data/GRCh38_no_alt.longreads.bed"
            tandem_repeats = "/app/CycVariantCallingPipline/SV/data/GRCh38_no_alt.trf.bed"
        } else if (use_longreads_region == "37") {
            longreads_region = "/app/CycVariantCallingPipline/SNV/data/hs37d5.longreads.bed"
            tandem_repeats = "/app/CycVariantCallingPipline/SV/data/hs37d5.trf.bed"
        } else {
            fail("The workflow is stopped because of the wrong ref version (chose frome '37' and '38').")
        }

        # /app/CycVariantCallingPipline/SNV/Freesia.bin \
        #     --env_dir ~{env_dir} \
        #     --dnb_bam_fn ~{dnb_bam_fn} \
        #     --cyclone_bam_fn ~{cyclone_bam_fn} \
        #     --ref_fn ~{ref_fn} \
        #     --model_path ~{model_path} \
        #     --longreads_region ${longreads_region} \
        #     --threads ~{threads/2} \
        #     --sample "~{prefix}_Freesia_SNV" \
        #     --out_dir ~{out_dir}
        
        conda run -n SV python /app/CycVariantCallingPipline/SV/src/Get_candidate_SV.py \
            -i ~{cyclone_bam_fn} \
            -o ~{out_dir} \
            -p "~{prefix}_Freesia_Candidate_SV" \
            -r ~{ref_fn} \
            --tandem-repeats ${tandem_repeats} \
            -t ~{threads}
        
        conda run -n Manta python /app/CycVariantCallingPipline/SV/src/Get_intergration/Get_intergration.py \
            --input_bam ~{dnb_bam_fn} \
            --out_dir ~{out_dir} \
            --prefix ~{prefix} \
            --ref_file ~{ref_fn} \
            --threads ~{threads} \
            --candidate_vcf "~{out_dir}/~{prefix}_Freesia_Candidate_SV_Candidate.vcf" \
            --extension 1500 \
            -hq "~{out_dir}/~{prefix}_Freesia_Candidate_SV_HighQ.vcf"
    >>>
    output{
        # File snv_vcf_file = "~{out_dir}/~{prefix}_Freesia_SNV.vcf.gz"
        # File snv_vcf_idx = "~{out_dir}/~{prefix}_Freesia_SNV.vcf.gz.tbi"
        File sv_vcf_file = "~{out_dir}/~{prefix}_Freesia_SV.vcf.gz"
        File sv_vcf_idx = "~{out_dir}/~{prefix}_Freesia_SV.vcf.gz.tbi"
    }
}

workflow Freesia{
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
    String dnb_bam_fn = "~{dnb_bam_dir}/{dnb_bam_name}"
    String cyclone_bam_fn = "~{cyclone_bam_dir}/{cyclone_bam_name}"
    
    call SNV_calling{
        input:
            dnb_bam_fn = dnb_bam_fn
            cyclone_bam_fn = cyclone_bam_fn
            ref_dir = ref_dir
            ref_prefix = ref_prefix
            ref_version = ref_version
            threads = threads
            prefix = prefix
            out_dir = out_dir
    }
    
    output{
        File sv_vcf_file = SNV_calling.sv_vcf_file
        File sv_vcf_file_idx = SNV_calling.sv_vcf_idx
    }
}