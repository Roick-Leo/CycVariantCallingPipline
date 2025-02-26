# 版本定义
version 1.0

task CycSV{
    input{
        File cyclone_bam_fn
        File ref_fn
        Int threads
        String out_dir
        String ref_version
    }
    
    command<<<
        if
        conda run -n SV python /app/CycVariantCallingPipline/CycSV.py \
        -i ~{cyclone_bam_fn} -o ~{out_dir} -r ~{ref_fn} --tandem-repeats /app/
    >>>

    output{
        File HighQ_vcf = "~{out_dir}/HighQ.vcf"
        File candidata_vcf = "~{out_dir}/Candidate.vcf"
    }
}
