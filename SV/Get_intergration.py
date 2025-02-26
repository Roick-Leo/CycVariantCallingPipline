import argparse
import os
import gzip
import re

def load_fai(fai_file):
    chrome_len = {}
    with open(fai_file, 'r') as f:
        for line in f:
            chrome, pos = line.strip().split()[0:2]
            chrome_len[chrome] = int(pos)
    return chrome_len

def process_positions(input_file, output_file, chrome_len, extension=1500):
    target_chr = [str(i) for i in range(1,23)] + ["X", "Y"] + ["chr" + str(i) for i in range(1,23)] + ["chrX", "chrY"]
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for lines in infile:
            if lines.startswith('#'):
                continue
            svtype = lines.split("SVTYPE=")[1].split(";")[0]
            if svtype == "DEL":
                sv_len = abs(int(re.split(r'[;\t]', lines.split("SVLEN=")[1])[0]))
            else:
                sv_len = 0
            line = lines.strip().split('\t')
            chr_name = line[0]
            pos = line[1]
            pos = int(pos)
            chr_len = chrome_len[chr_name]
            start = max(0, pos - extension)
            end = min(chr_len, pos + sv_len + extension)
            if chr_name in target_chr:
                outfile.write(f"{chr_name}\t{start}\t{end}\n")
                
def intergration(cyc_vcf_HighQ_file, cyc_vcf_file, dnb_vcf_file, out_vcf_file, extension = 1500):
    result = []
    with open(cyc_vcf_HighQ_file, "r") as highQ, open(cyc_vcf_file, "r") as Cyc_vcf, gzip.open(dnb_vcf_file, "rt") as dnb_vcf, open(out_vcf_file, "w") as out_vcf:
        cyc_dict = {}
        vcf_header = []
        header_info = []
        for line in Cyc_vcf:
            if line.startswith("##"):
                if "command" in line:
                    continue
                elif line.startswith("##FILTER"):
                    header_info.append(line.split("FILTER=<ID=")[1].split(",")[0])
                    vcf_header.append(line)
                elif line.startswith("##INFO"):
                    header_info.append(line.split("INFO=<ID=")[1].split(",")[0])
                    vcf_header.append(line)
                elif line.startswith("##FORMAT"):
                    header_info.append(line.split("FORMAT=<ID=")[1].split(",")[0])
                    vcf_header.append(line)
                else:
                    out_vcf.write(line)
                continue
            if line.startswith("#"):
                continue
            cyc_chrom, cyc_pos_start, cyc_pos_end, cyc_svtype, cyc_svlen = trans_line_feature(line, extension)
            if cyc_chrom in cyc_dict.keys():
                cyc_dict[cyc_chrom].append((cyc_chrom, cyc_pos_start, cyc_pos_end, cyc_svtype, cyc_svlen))
            else:
                cyc_dict[cyc_chrom] = [(cyc_chrom, cyc_pos_start, cyc_pos_end, cyc_svtype, cyc_svlen)]
        for dnb_line in dnb_vcf:
            if dnb_line.startswith("#"):
                if dnb_line.startswith("##"):
                    header_if = False
                    if dnb_line.startswith("##FILTER"):
                        header_if = dnb_line.split("FILTER=<ID=")[1].split(",")[0]
                    elif dnb_line.startswith("##INFO"):
                        header_if = dnb_line.split("INFO=<ID=")[1].split(",")[0]
                    elif dnb_line.startswith("##FORMAT"):
                        header_if = dnb_line.split("FORMAT=<ID=")[1].split(",")[0]
                    else:
                        pass
                    if header_if and header_if not in header_info:
                        vcf_header.append(dnb_line)
                continue
            dnb_chrom, dnb_pos_start, dnb_pos_end, dnb_svtype, dnb_svlen = trans_line_feature(dnb_line, extension)
            if dnb_chrom not in cyc_dict.keys():
                continue
            merged = False
            start_idx = 0
            for idx, sv_record in enumerate(cyc_dict[dnb_chrom]):
                d = dnb_pos_start - sv_record[1]
                if d > 1500:
                    start_idx = idx
                    continue
                if d < -1500:
                    break
                if dnb_pos_end >= sv_record[1] and dnb_pos_end <= sv_record[2] and dnb_svtype == sv_record[3]:
                    merged = True
                elif dnb_pos_start >= sv_record[1] and dnb_pos_start <= sv_record[2] and dnb_svtype == sv_record[3]:
                    merged = True
                else:
                    pass
                if merged:
                    break
            if start_idx > 0:
                cyc_dict[dnb_chrom] = cyc_dict[dnb_chrom][start_idx::]
            if merged and dnb_svtype != "BND":
                result.append(dnb_line.replace("Manta", "Freesia."))
        for line in sorted(vcf_header):
            out_vcf.write(line)
        for highq_line in highQ:
            if highq_line.startswith("##"):
                continue
            if highq_line.startswith("#"):
                out_vcf.write(highq_line)
                continue
            out_vcf.write(highq_line.replace("Sniffles2", "Freesia"))
        for line in result:
            out_vcf.write(line)
            
def trans_line_feature(line, extension): 
    ll = line.strip().split("\t")
    chrom = ll[0]
    pos = int(ll[1])
    svtype = line.split("SVTYPE=")[1].split(";")[0]
    if svtype == "DUP":
        svtype = "INS"
    if "SVLEN" in line:
        svlen = abs(int(re.split(r'[;\t]', line.split("SVLEN=")[1])[0]))
    else:
        svlen = 0
    pos_start = pos - max(100,min(0.3*svlen, extension))
    if svtype == "DEL":
        pos_end = pos + svlen + max(100, min(0.3*svlen, extension))
    else:
        pos_end = pos +max(100, min(0.3*svlen, extension))
    return chrom, pos_start, pos_end, svtype, svlen

def main(threads, bam_file, ref_file, HighQ_vcf_file, candidate_vcf_file, out_dir, prefix, extension=1500):
    out_target_bed = os.path.join(out_dir, prefix+"_target.bed")
    target_bam = os.path.join(out_dir, prefix+"_target.sort.bam")
    fai_file = ref_file + ".fai"
    chrome_len = load_fai(fai_file)
    manta_dir = os.path.join(out_dir, "manta")
    process_positions(candidate_vcf_file, out_target_bed, chrome_len, 1500)
    os.system(f"samtools view -@ {threads} -Shb {bam_file} -L {out_target_bed} | samtools sort -@ {threads} -o {target_bam} ; samtools index -@ {threads} {target_bam}")
    os.system(f"configManta.py --bam {target_bam} --referenceFasta {ref_file} --runDir {manta_dir} && {manta_dir}/runWorkflow.py -m local -j {threads}")
    dnb_vcf_file = f"{manta_dir}/results/variants/diploidSV.vcf.gz"
    out_vcf_file = os.path.join(out_dir, prefix+"_Freesia_SV.vcf")
    intergration(HighQ_vcf_file, candidate_vcf_file, dnb_vcf_file, out_vcf_file, extension)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='...')
    parser.add_argument('-i','--input_bam',dest='input_bam', type=str, help='input file')
    parser.add_argument('-o','--out_dir',dest='out_dir', type=str, help='outfile',default='out')
    parser.add_argument('-p','--prefix',dest='prefix', type=str, help='outfile',default='Demo')
    parser.add_argument('-r','--ref_file',dest='ref_file', type=str, help='ref file')
    parser.add_argument('-t','--threads', default= 15, type=int, help='threads')
    parser.add_argument('-c','--candidate_vcf', type=str, help='candidate_vcf')
    parser.add_argument('-hq','--HighQ_vcf', type=str, help='HighQ_vcf')
    parser.add_argument('-e','--extension',dest='extension', type=int, help='extension',default=1500)
    args = parser.parse_args()

    input_file = args.input_bam
    out_dir = args.out_dir
    prefix = args.prefix
    ref_file = args.ref_file
    threads = args.threads
    candidate_vcf_file = args.candidate_vcf
    HighQ_vcf_file = args.HighQ_vcf
    extension = args.extension

    main(threads, input_file, ref_file, HighQ_vcf_file, candidate_vcf_file, out_dir, prefix, extension)