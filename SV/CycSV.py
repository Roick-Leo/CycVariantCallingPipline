#!/zfssz2/ST_BIOINTEL/P17H10200N0325/01.User/linwei4/soft/miniforge3/envs/sniffles2/bin/python
# -*- encoding: utf-8 -*-

"""
Created on Mon Jan 06 18:17:25 2025

@File    :   get_supp2.py

@Author : linwei
"""
import pysam
from dataclasses import dataclass
import time
import sys
import statistics

import argparse
from argparse import Namespace
import os
import datetime

import itertools

import math

from typing import List, Tuple, Dict, Any
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from collections import defaultdict
############ 0.0 config ############
def get_config():
    VERSION="Sniffles2"
    BUILD="2.0.7"
    SNF_VERSION="S2_rc4"
    configName=Namespace(
        non_germline = False,
        phase = False,
        span=10000,
        hapcpu=4,
        NGSmapq_min=25,
        NGSalen_min=10,
        abnormal=3,
        abnormal_ratio=0.01,
        #disstart=10000,
        #dislen=500,
        #signal=3,
        #depth=0,
        mode="call_sample",
        minsupport = "auto",
        minsupport_auto_mult = None,
        minsvlen = 35,
        minsvlen_screen_ratio = 0.9,
        mapq = 25,
        no_qc = False,
        qc_stdev = True,
        qc_stdev_abs_max = 500,
        qc_coverage = 1,
        long_ins_length = 2500,
        long_del_length = 50000,
        long_del_coverage = 0.66,
        long_dup_length = 50000,
        long_dup_coverage = 1.33,
        max_splits_kb = 0.1,
        max_splits_base = 3,
        min_alignment_length = 1000,
        phase_conflict_threshold = 0.1,
        detect_large_ins = True,
        large_ins_threshold = 5000,
        cluster_binsize = 100,
        cluster_r = 2.5,
        cluster_repeat_h = 1.5,
        cluster_repeat_h_max = 1000,
        cluster_merge_pos = 150,
        cluster_merge_len = 0.33,
        cluster_merge_bnd = 1500,
        genotype_ploidy = 2,
        genotype_error = 0.05,
        sample_id = None,
        genotype_vcf = None,
        combine_high_confidence = 0.0,
        combine_low_confidence = 0.2,
        combine_low_confidence_abs = 2,
        combine_null_min_coverage = 5,
        combine_match = 250,
        combine_match_max = 1000,
        combine_separate_intra = False,
        combine_output_filtered = False,
        combine_exhaustive = False,
        combine_relabel_rare = False,
        combine_with_missing = False,
        output_rnames = False,
        no_consensus = False,
        no_sort = False,
        no_progress = False,
        quiet = False,
        max_del_seq_len = 50000,
        symbolic = False,
        allow_overwrite = True,
        dev_cache = False,
        dev_cache_dir = None,
        dev_debug_svtyping = False,
        dev_keep_lowqual_splits = False,
        dev_call_region = None,
        dev_dump_clusters = False,
        dev_merge_inline = False,
        dev_seq_cache_maxlen = 50000,
        consensus_max_reads = 20,
        consensus_max_reads_bin = 10,
        combine_consensus = False,
        dev_dump_coverage = False,
        dev_no_resplit = False,
        dev_no_resplit_repeat = False,
        dev_skip_snf_validation = False,
        low_memory = False,
        repeat = False,
        qc_nm = False,
        qc_nm_max = 0.2,
        coverage_updown_bins = 5,
        coverage_shift_bins = 3,
        coverage_shift_bins_min_aln_length = 1000,
        cluster_binsize_combine_mult = 5,
        cluster_resplit_binsize = 20,
        qc_strand = False,)
    config=configName
    if config.quiet:
        sys.stdout=open(os.devnull,"w")

    config.start_date=datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")

    config.sort=not config.no_sort
    #if config.low_memory:
    #    config.task_count_multiplier=64
    #else:
    #    config.task_count_multiplier=1
    config.task_count_multiplier=0

    config.version=VERSION
    config.build=BUILD
    config.snf_format_version=SNF_VERSION
    config.command=" ".join(sys.argv)

    if config.dev_call_region != None:
        region_contig,region_startend=config.dev_call_region.replace(",","").split(":")
        start,end=region_startend.split("-")
        config.dev_call_region=dict(contig=region_contig,start=int(start),end=int(end))

    #"--minsvlen" parameter is for final output filtering
    #for intermediate steps, a lower threshold is used to account for sequencing, mapping imprecision
    config.minsvlen_screen=int(config.minsvlen_screen_ratio*config.minsvlen)
    #config.minsupport_screen=max(1,int(0.333*config.minsupport*(config.cluster_binsize/100.0)))

    if config.minsupport!="auto":
        config.minsupport=int(config.minsupport)

    #--minsupport auto defaults
    config.minsupport_auto_base=1.5
    config.minsupport_auto_regional_coverage_weight=0.75

    if config.minsupport_auto_mult==None:
        if config.non_germline:
            config.minsupport_auto_mult=0.025
        else:
            config.minsupport_auto_mult=0.1

    if config.non_germline:
        config.qc_nm=True

    config.coverage_binsize=config.cluster_binsize
    config.coverage_binsize_combine=config.cluster_binsize*config.cluster_binsize_combine_mult


    #INS Consensus parameters
    #config.consensus_max_reads=20
    #config.consensus_max_reads_bin=10
    config.consensus_min_reads=4
    config.consensus_kmer_len=6
    config.consensus_kmer_skip_base=3
    config.consensus_kmer_skip_seqlen_mult=1.0/500.0
    config.consensus_low_threshold=0.0 #0.15

    #Large INS
    config.long_ins_rescale_base=1.66
    config.long_ins_rescale_mult=0.33

    #BND
    config.bnd_cluster_length=1000
    config.bnd_cluster_resplit=0

    #Genotyping
    config.genotype_format="GT:GQ:DR:DV"
    config.genotype_none=(".",".",0,0,0,None)
    config.genotype_null=(0,0,0,0,0,None)
    config.genotype_min_z_score=5
    if config.genotype_ploidy!=2:
        util.fatal_error("Currently only --genotype-ploidy 2 is supported")

    #SNF
    config.snf_block_size=10**5
    config.snf_combine_keep_open=True #Keep file handles open during .snf combining (might be an issue if the number of .snf files to merge is very large)

    #Combine
    config.combine_exhaustive=False
    config.combine_relabel_rare=False
    config.combine_overlap_abs=2500
    config.combine_min_size=100

    #Misc
    config.precise=25 #Max. sum of pos and length stdev for SVs to be labelled PRECISE
    config.tandem_repeat_region_pad=500
    config.id_prefix="Sniffles2."
    config.phase_identifiers=["1","2"]

    config.dev_profile=False

    config.workdir=os.getcwd()

    return config

############ 0.1 util ############
class util:
    @staticmethod
    def stdev(nums):
        nums = list(nums)
        return statistics.stdev(nums) if len(nums) > 1 else 0

    @staticmethod
    def variance(nums):
        nums = list(nums)
        return statistics.variance(nums) if len(nums) > 1 else 0
    
    @staticmethod
    def pvariance(nums):
        nums = list(nums)
        return statistics.pvariance(nums) if len(nums) > 1 else 0

    @staticmethod
    def median(nums):
        return int(statistics.median(nums))

    @staticmethod
    def median_or_mode(nums):
        nums = list(nums)
        top = util.most_common(nums)
        if len(top) > 1 and (top[0][0] - top[1][0] < 2):
            return util.median_noavg(nums)
        else:
            return util.median_modes(nums)

    @staticmethod
    def median_noavg(nums):
        nums = sorted(list(nums))
        mid = int(len(nums) / 2)
        return nums[mid]

    @staticmethod
    def median_modes(nums):
        max_count = 0
        counts = {}
        for n in nums:
            if n not in counts:
                counts[n] = 1
            else:
                counts[n] += 1
            max_count = max(max_count, counts[n])
        return util.median_noavg([k for k, n in counts.items() if max_count - n < 3])

    @staticmethod
    def mean(nums):
        nums = list(nums)
        return sum(nums) / len(nums)

    @staticmethod
    def mean_or_none(nums):
        nums = list(nums)
        if len(nums) == 0:
            return None
        else:
            return sum(nums) / len(nums)

    @staticmethod
    def mean_or_none_round(nums):
        r = util.mean_or_none(nums)
        if r is None:
            return r
        else:
            return round(r)

    @staticmethod
    def trim(nums, pct=25):
        nums = sorted(list(nums))
        trim_n = int(len(nums) / float(100.0) * pct)
        if trim_n > 0:
            return nums[trim_n:-trim_n]
        else:
            return nums

    @staticmethod
    def most_common(nums):
        counts = {}
        for n in nums:
            if n not in counts:
                counts[n] = 1
            else:
                counts[n] += 1
        return sorted(((counts[n], n) for n in counts), reverse=True)

    @staticmethod
    def most_common_top(nums):
        result = util.most_common(nums)
        return sorted(item for count, item in result if count == result[0][0])[0]

    @staticmethod
    def error(msg):
        sys.stderr.write("Sniffles2 Error: " + msg + "\n")
        sys.stderr.flush()

    @staticmethod
    def fatal_error(msg):
        util.error(msg + " (Fatal error, exiting.)")
        exit(1)


    @staticmethod
    def load_tandem_repeats(filename, padding):
        contigs_tr = {}
        unsorted = False
        with open(filename, "r") as handle:
            for line in handle:
                parts = line.split("\t")
                if len(parts) >= 3:
                    contig, start, end = parts[:3]
                    start = int(start)
                    end = int(end)
                    if contig not in contigs_tr:
                        contigs_tr[contig] = []
                    if len(contigs_tr[contig]) > 0:
                        last_start, last_end = contigs_tr[contig][-1]
                        if start < last_start:
                            unsorted = True
                    contigs_tr[contig].append((max(0, int(start) - padding), int(end) + padding))

        if unsorted:
            #print("Info: The tandem repeat annotations file was not sorted. Sorting it in-memory after loading... (please sort the .bed file once before to save time when running multiple samples)")
            sort_start = time.time()
            for contig in contigs_tr:
                contigs_tr[contig].sort()
            #print(f"Info: Optional sorting of input tandem repeat annotations took {time.time() - sort_start:.2f}s.")

        return contigs_tr

    @staticmethod
    def center(nums):
        return util.median_modes(nums)

TYPES=["INS","DEL","DUP","INV","BND"]


################ 1.get_SV_signal ################
############ 1.0 lead ############
@dataclass
class Lead:
    read_id: int=None
    read_qname: str=None
    contig: str=None
    ref_start: int=None
    ref_end: int=None
    qry_start: int=None
    qry_end: int=None
    strand: str=None
    mapq: int=None
    nm: float=None
    source: str=None
    svtype: str=None
    svlen: int=None
    seq: str=None
    svtypes_starts_lens: list=None

######## 1.0.1 ########    
OPTAB={pysam.CMATCH:     (1,1,0),
       pysam.CEQUAL:     (1,1,0),
       pysam.CDIFF:      (1,1,0),
       pysam.CINS:       (1,0,1),
       pysam.CDEL:       (0,1,1),
       pysam.CREF_SKIP:  (0,1,0),
       pysam.CSOFT_CLIP: (1,0,1),
       pysam.CHARD_CLIP: (0,0,0),
       pysam.CPAD:       (0,0,0)}
#      pysam.CBACK:      (0,0,0)}
OPLIST=[(0,0,0) for i in range(max(int(k) for k in OPTAB.keys())+1)]
for k,v in OPTAB.items():
    OPLIST[int(k)]=v

######## 1.0.2  ########
def CIGAR_analyze(cigar):
    buf=""
    readspan=0
    refspan=0
    clip_start=None
    clip=0
    for c in cigar:
        if c.isnumeric():
            buf+=c
        else:
            oplen=int(buf)
            h=False
            if c in "MIX=":
                readspan+=oplen
                h=True
            if c in "MDX=N":
                refspan+=oplen
                h=True
            if not h:
                if c in "SH":
                    if clip_start==None and readspan+refspan>0:
                        clip_start=clip
                    clip+=oplen
                else:
                    raise f"Unknown CIGAR operation: '{c}'"
            buf=""
    if clip_start==None:
        clip_start=clip
    return clip_start, clip-clip_start, refspan, readspan

############ 1.1 from CIGAR ############
def read_iterindels(read_id,read,contig,config,use_clips,read_nm):
    minsvlen=config.minsvlen_screen     #0.9*35
    longinslen=config.long_ins_length/2.0
    seq_cache_maxlen=config.dev_seq_cache_maxlen
    qname=read.query_name
    mapq=read.mapping_quality
    strand="-" if read.is_reverse else "+"
    CINS=pysam.CINS
    CDEL=pysam.CDEL
    CSOFT_CLIP=pysam.CSOFT_CLIP

    pos_read=0
    pos_ref=read.reference_start
    for op,oplength in read.cigartuples:
        add_read,add_ref,event=OPLIST[op]
        if event and oplength >= minsvlen:  #31
            if op==CINS:
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref,
                           pos_ref,
                           pos_read,
                           pos_read+oplength,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "INS",
                           oplength,
                           seq=read.query_sequence[pos_read:pos_read+oplength] if oplength <= seq_cache_maxlen else None)
            elif op==CDEL:
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref+oplength,
                           pos_ref,
                           pos_read,
                           pos_read,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "DEL",
                           -oplength)
            elif use_clips and op==CSOFT_CLIP and oplength >= longinslen:       #为了防止有些reads的首尾刚好在INS区段中
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref,
                           pos_ref,
                           pos_read,
                           pos_read+oplength,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "INS",
                           None,
                           seq=None)
        pos_read+=add_read*oplength
        pos_ref+=add_ref*oplength

############ 1.2 from split reads ############
######## 1.2.0  Handle the list of SA, the part that really determines the SV type ########
@dataclass
class SVCallBNDInfo:
    mate_contig: str
    mate_ref_start: int
    is_first: bool
    is_reverse:bool
def classify_splits(read,leads,config,main_contig):
    minsvlen_screen=config.minsvlen_screen
    maxsvlen_other=minsvlen_screen*5

    leads.sort(key=lambda ld: ld.qry_start)
    last=leads[0]
    last.svtypes_starts_lens=[]

    if last.qry_start >= config.long_ins_length*0.5:
        last.svtypes_starts_lens.append(("INS",last.ref_start,None))

    for i in range(1,len(leads)):
        curr=leads[i]
        curr.svtypes_starts_lens=[]
        qry_dist_abs=abs(curr.qry_start - last.qry_end)

        if curr.contig == last.contig:
            rev=(curr.strand == "-")
            fwd=not rev
            if curr.strand == last.strand:
                #
                #INS, DEL, DUP
                #
                if (fwd and (curr.qry_start - last.qry_end) >= minsvlen_screen
                    and (curr.ref_start - last.ref_end) < maxsvlen_other
                    and (curr.qry_start - last.qry_end) - (curr.ref_start - last.ref_end) >= minsvlen_screen):
                    #INS, FWD
                    svstart=curr.ref_start
                    svlen=(curr.qry_start - last.qry_end)
                    if svlen <= config.dev_seq_cache_maxlen:
                        curr.seq=read.query_sequence[last.qry_end:curr.qry_start]
                    else:
                        curr.seq=None
                    curr.svtypes_starts_lens.append(("INS",svstart,svlen))

                elif (rev and (curr.qry_start - last.qry_end) >= minsvlen_screen
                      and (last.ref_start - curr.ref_end) < maxsvlen_other
                      and (curr.qry_start - last.qry_end) - (last.ref_start - curr.ref_end) >= minsvlen_screen):
                    #INS, REV
                    svstart=last.ref_start
                    svlen=(curr.qry_start - last.qry_end)
                    if svlen <= config.dev_seq_cache_maxlen:
                        curr.seq=read.query_sequence[last.qry_end:curr.qry_start]
                    else:
                        curr.seq=None
                    curr.svtypes_starts_lens.append(("INS",svstart,svlen))

                elif (fwd and (curr.ref_start - last.ref_end) >= minsvlen_screen
                      and qry_dist_abs < maxsvlen_other
                      and (curr.ref_start - last.ref_end) - (curr.qry_start - last.qry_end) >= minsvlen_screen):
                        #DEL, FWD
                        svstart=curr.ref_start
                        svlen=(curr.ref_start - last.ref_end)
                        curr.svtypes_starts_lens.append(("DEL",svstart,-svlen))

                elif (rev and (last.ref_start - curr.ref_end) >= minsvlen_screen
                      and qry_dist_abs < maxsvlen_other
                      and (last.ref_start - curr.ref_end) - (curr.qry_start - last.qry_end) >= minsvlen_screen):
                        #DEL, REV
                        svstart=last.ref_start
                        svlen=(last.ref_start - curr.ref_end)
                        curr.svtypes_starts_lens.append(("DEL",svstart,-svlen))

                elif fwd and curr.ref_start <= last.ref_end:
                    if qry_dist_abs < maxsvlen_other:
                        #DUP, FWD
                        svstart=curr.ref_start
                        svlen=(last.ref_end - curr.ref_start)
                        if svlen >= minsvlen_screen:
                            curr.svtypes_starts_lens.append(("DUP",svstart,svlen))

                elif rev and last.ref_start <= curr.ref_end:
                    if qry_dist_abs < maxsvlen_other:
                        #DUP, REV
                        svstart=last.ref_start
                        svlen=(curr.ref_end - last.ref_start)
                        if svlen >= minsvlen_screen:
                            curr.svtypes_starts_lens.append(("DUP",svstart,svlen))

            elif qry_dist_abs < maxsvlen_other:
                #
                #INV
                #
                if fwd and curr.ref_start <= last.ref_start:
                    #CASE B
                    svstart=curr.ref_start
                    svlen=last.ref_start-curr.ref_start
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV",svstart,svlen))

                elif fwd and curr.ref_start > last.ref_start:
                    #CASE C
                    svstart=last.ref_start
                    svlen=curr.ref_start-last.ref_start
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV",svstart,svlen))

                elif rev and curr.ref_end >= last.ref_end:
                    #CASE A
                    svstart=last.ref_end
                    svlen=curr.ref_end-last.ref_end
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV",svstart,svlen))

                elif rev and curr.ref_end < last.ref_end:
                    #CASE D
                    svstart=curr.ref_end
                    svlen=last.ref_end-curr.ref_end
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV",svstart,svlen))

        elif qry_dist_abs < maxsvlen_other:
            #
            #BND
            #
            if curr.contig == main_contig:
                a,b=curr,last
            else:
                a,b=last,curr

            if a.contig == main_contig:
                is_first=a.qry_start < b.qry_start
                if is_first:
                    if a.strand=="+":
                        svstart=a.ref_end
                    else:
                        svstart=a.ref_start
                else:
                    if a.strand=="+":
                        svstart=a.ref_start
                    else:
                        svstart=a.ref_end
                a.svtypes_starts_lens.append(("BND",
                                            svstart,
                                            SVCallBNDInfo(b.contig,
                                                          b.ref_start,
                                                          is_first,
                                                          a.strand!=b.strand)))
        last=curr
######## 1.2.1 from supp_reads ########
def read_itersplits_bnd(read_id,read,contig,config,read_nm):
    assert(read.is_supplementary)
    #SA:refname,pos,strand,CIGAR,MAPQ,NM
    all_leads=[]
    supps=[part.split(",") for part in read.get_tag("SA").split(";") if len(part)>0]

    if len(supps) > config.max_splits_base + config.max_splits_kb*(read.query_length/1000.0):
        return

    if read.is_reverse:
        qry_start=read.query_length-read.query_alignment_end
    else:
        qry_start=read.query_alignment_start

    curr_lead=Lead(read_id,
                   read.query_name,
                   contig,
                   read.reference_start,
                   read.reference_start+read.reference_length,
                   qry_start,
                   qry_start+read.query_alignment_length,
                   "-" if read.is_reverse else "+",
                   read.mapping_quality,
                   read_nm,
                   "SPLIT_SUP",
                   "?")
    all_leads.append(curr_lead)

    prim_refname,prim_pos,prim_strand,prim_cigar,prim_mapq,prim_nm=supps[0]
    if prim_refname == contig:
        #Primary alignment is on this chromosome, no need to parse the supplementary
        return

    minpos_curr_chr=min(itertools.chain([read.reference_start],(int(pos) for refname,pos,strand,cigar,mapq,nm in supps if refname==contig)))
    if minpos_curr_chr < read.reference_start:
        #Only process splits once per chr (there may be multiple supplementary alignments on the same chr)
        return

    for refname,pos,strand,cigar,mapq,nm in supps:
        mapq=int(mapq)
        nm=int(nm)
        #if not config.dev_keep_lowqual_splits and mapq < config.mapq:
        #    continue

        is_rev=(strand=="-")

        try:
            readstart_fwd,readstart_rev,refspan,readspan=CIGAR_analyze(cigar)
        except Exception as e:
            util.error(f"Malformed CIGAR '{cigar}' with pos {pos} of read '{read.query_name}' ({e}). Skipping.")
            return

        pos_zero=int(pos)-1
        split_qry_start=readstart_rev if is_rev else readstart_fwd

        all_leads.append(Lead(read_id,
                              read.query_name,
                              refname,
                              pos_zero,
                              pos_zero + refspan,
                              split_qry_start,
                              split_qry_start+readspan,
                              strand,
                              mapq,
                              nm/float(readspan+1),
                              "SPLIT_SUP",
                              "?"))

    classify_splits(read,all_leads,config,contig)

    for lead in all_leads:
        for svtype, svstart, arg in lead.svtypes_starts_lens:
            if svtype=="BND":
                bnd = Lead(read_id=lead.read_id,
                           read_qname=lead.read_qname,
                           contig=lead.contig,
                           ref_start=svstart,
                           ref_end=svstart,
                           qry_start=lead.qry_start,
                           qry_end=lead.qry_end,
                           strand=lead.strand,
                           mapq=lead.mapq,
                           nm=lead.nm,
                           source=lead.source,
                           svtype=svtype,
                           svlen=config.bnd_cluster_length,
                           seq=None)
                bnd.bnd_info=arg
                yield bnd

######## 1.2.2 from has_sa but not supp_reads ########
def read_itersplits(read_id,read,contig,config,read_nm):
    #SA:refname,pos,strand,CIGAR,MAPQ,NM
    all_leads=[]
    supps=[part.split(",") for part in read.get_tag("SA").split(";") if len(part)>0]

    if len(supps) > config.max_splits_base + config.max_splits_kb*(read.query_length/1000.0):   #3+0.1*(query_len/1000),也就是10k以内为3，每多 10k supp 数目可以 +1
        return

    if read.is_reverse:
        qry_start=read.query_length-read.query_alignment_end
    else:
        qry_start=read.query_alignment_start

    curr_lead=Lead(read_id,
                   read.query_name,
                   contig,
                   read.reference_start,
                   read.reference_start+read.reference_length,
                   qry_start,
                   qry_start+read.query_alignment_length,
                   "-" if read.is_reverse else "+",
                   read.mapping_quality,
                   read_nm,
                   "SPLIT_PRIM",
                   "?")
    all_leads.append(curr_lead)

    for refname,pos,strand,cigar,mapq,nm in supps:
        mapq=int(mapq)
        nm=int(nm)
        #if not config.dev_keep_lowqual_splits and mapq < config.mapq:
        #    continue

        is_rev=(strand=="-")

        try:
            readstart_fwd,readstart_rev,refspan,readspan=CIGAR_analyze(cigar)
        except Exception as e:
            util.error(f"Malformed CIGAR '{cigar}' with pos {pos} of read '{read.query_name}' ({e}). Skipping.")
            return

        pos_zero=int(pos)-1
        split_qry_start=readstart_rev if is_rev else readstart_fwd

        all_leads.append(Lead(read_id,
                              read.query_name,
                              refname,
                              pos_zero,
                              pos_zero + refspan,
                              split_qry_start,
                              split_qry_start+readspan,
                              strand,
                              mapq,
                              nm/float(readspan+1), # -1/(read_span+1)
                              "SPLIT_SUP",
                              "?"))

        #QC on: 08Sep21; O.K.
        #cigarl=CIGAR_tolist(cigar)
        #assert(CIGAR_listrefspan(cigarl)==refspan)
        #assert(CIGAR_listreadspan(cigarl)==readspan)
        #assert(CIGAR_listreadstart_fwd(cigarl)==readstart_fwd)
        #assert(CIGAR_listreadstart_rev(cigarl)==readstart_rev)
        #End QC

    classify_splits(read,all_leads,config,contig)

    for lead_i, lead in enumerate(all_leads):
        for svtype, svstart, arg in lead.svtypes_starts_lens:
            min_mapq=min(lead.mapq,all_leads[max(0,lead_i-1)].mapq)
            if not config.dev_keep_lowqual_splits and min_mapq < config.mapq:
                continue

            if svtype=="BND":
                bnd = Lead(read_id=lead.read_id,
                           read_qname=lead.read_qname,
                           contig=lead.contig,
                           ref_start=svstart,
                           ref_end=svstart,
                           qry_start=lead.qry_start,
                           qry_end=lead.qry_end,
                           strand=lead.strand,
                           mapq=lead.mapq,
                           nm=lead.nm,
                           source=lead.source,
                           svtype=svtype,
                           svlen=config.bnd_cluster_length,
                           seq=None)
                bnd.bnd_info=arg
                yield bnd

            elif svtype!="NOSV":
                svlen=arg

                yield Lead(read_id=lead.read_id,
                           read_qname=lead.read_qname,
                           contig=lead.contig,
                           ref_start=svstart,
                           ref_end=svstart+svlen if svlen!=None and svtype!="INS" else svstart,
                           qry_start=lead.qry_start,
                           qry_end=lead.qry_end,
                           strand=lead.strand,
                           mapq=lead.mapq,
                           nm=lead.nm,
                           source=lead.source,
                           svtype=svtype,
                           svlen=svlen,
                           seq=lead.seq if svtype=="INS" else None)

############ 1. Call the above function to perform preliminary identification of SVs ---- ensure full detection. ############
class LeadProvider:
    def __init__(self,config,read_id_offset):
        self.config=config

        self.leadtab={}
        self.leadcounts={}

        for svtype in TYPES:
            self.leadtab[svtype]={}
            self.leadcounts[svtype]=0

        self.covrtab_fwd={}
        self.covrtab_rev={}
        self.covrtab_min_bin=None
        #self.covrtab_read_start={}
        #self.covrtab_read_end={}

        self.read_id=read_id_offset
        self.read_count=0

        self.contig=None
        self.start=None
        self.end=None
    # get leadtab,leadlab={svtype1:{pos1:{lead1,lead2},pos2:{lead3,lead4}},}
    def record_lead(self,ld,pos_leadtab):
        leadtab_svtype=self.leadtab[ld.svtype]
        if pos_leadtab in leadtab_svtype:
            leadtab_svtype[pos_leadtab].append(ld)
            lead_count=len(leadtab_svtype[pos_leadtab])
            if lead_count > self.config.consensus_max_reads_bin:
                ld.seq=None
        else:
            leadtab_svtype[pos_leadtab]=[ld]
            lead_count=1
        self.leadcounts[ld.svtype]+=1
    ####  遍历bam, 调用前述函数根据 cigar、split 信息进行初步 sv 判定，逐个输出 lead  ####
    '''
    # 传入 bam 以及区间信息, 读取 bam, 遍历 read
    # 对每条 reads, 使用函数 read_iterindels, 根据 read 中 CIGAR 值初步判断 read 上的 SV, 一个 SV 传出一个 lead
    # 如果 has_sa: is_supp 则使用函数 read_itersplits_bnd, 根据在 read 上相邻的两个 supp 比对的位置信息判断SV类型, 输出 BND, 以 lead 形式
    #              else: 使用函数 read_itersplits, 输出 SV
    # 调用该函数得到 多个 lead,  可以用 for 遍历
    # 注意：用该函数记得先调用 build_lead, 初始化 start/end 等信息
    '''
    def build_leadtab(self,contig,start,end,bam):
        if self.config.dev_cache:
            loaded_externals=self.dev_load_leadtab(contig,start,end)
            if loaded_externals!=False:
                return loaded_externals

        assert(self.contig==None)
        assert(self.start==None)
        assert(self.end==None)
        self.contig=contig
        self.start=start
        self.end=end
        self.covrtab_min_bin=int(self.start/self.config.coverage_binsize)*self.config.coverage_binsize

        externals=[]
        ld_binsize=self.config.cluster_binsize

        for ld in self.iter_region(bam,contig,start,end):
            ld_contig,ld_ref_start=ld.contig,ld.ref_start
            #TODO: Handle leads overlapping region ends (start/end)
            if contig==ld_contig and ld_ref_start >= start and ld_ref_start < end:
                pos_leadtab=int(ld_ref_start/ld_binsize)*ld_binsize
                self.record_lead(ld,pos_leadtab)
            else:
                externals.append(ld)

        if self.config.dev_cache:
            self.dev_store_leadtab(contig,start,end,externals)

        return externals

    def iter_region(self,bam,contig,start=None,end=None):
        leads_all=[]
        binsize=self.config.cluster_binsize
        coverage_binsize=self.config.coverage_binsize
        coverage_shift_bins=self.config.coverage_shift_bins
        coverage_shift_min_aln_len=self.config.coverage_shift_bins_min_aln_length
        long_ins_threshold=self.config.long_ins_length*0.5
        qc_nm=self.config.qc_nm
        phase=self.config.phase
        advanced_tags=qc_nm or phase    #False
        mapq_min=self.config.mapq
        alen_min=self.config.min_alignment_length   #1000

        for read in bam.fetch(contig,start,end,until_eof=False):
            #if self.read_count % 1000000 == 0:
            #    gc.collect()
            if read.reference_start < start or read.reference_start >= end: #整条reads起始点在区间，没有限制end
                continue

            self.read_id+=1
            self.read_count+=1

            alen=read.query_alignment_length
            if read.mapping_quality < mapq_min or read.is_secondary or alen < alen_min:
                continue

            has_sa=read.has_tag("SA")
            use_clips=self.config.detect_large_ins and not read.is_supplementary and not has_sa

            nm=-1
            curr_read_id=self.read_id
            if advanced_tags:
                if qc_nm: # nm 的差异是2.0.7是nm/reads_length,2.3.3是 mismatch/reads_length
                    if read.has_tag("NM"):
                        nm=read.get_tag("NM")/float(read.query_alignment_length+1)

                if phase:
                    curr_read_id=(self.read_id,str(read.get_tag("HP")) if read.has_tag("HP") else "NULL",str(read.get_tag("PS")) if read.has_tag("PS") else "NULL")

            #Extract small indels
            for lead in read_iterindels(curr_read_id,read,contig,self.config,use_clips,read_nm=nm): #nm=-1,默认advanced_tags=False
                yield lead

            #Extract read splits
            if has_sa:
                if read.is_supplementary:
                    for lead in read_itersplits_bnd(curr_read_id,read,contig,self.config,read_nm=nm):
                        yield lead
                else:
                    for lead in read_itersplits(curr_read_id,read,contig,self.config,read_nm=nm):
                        yield lead

            #Record in coverage table
            read_end=read.reference_start+read.reference_length
            assert(read_end==read.reference_end)
            #assert(read_end>=read.reference_start)
            if read.is_reverse:
                target_tab=self.covrtab_rev
            else:
                target_tab=self.covrtab_fwd
            covr_start_bin=(int(read.reference_start/coverage_binsize)+coverage_shift_bins*(alen>=coverage_shift_min_aln_len))*coverage_binsize
            covr_end_bin=(int(read_end/coverage_binsize)-coverage_shift_bins*(alen>=coverage_shift_min_aln_len))*coverage_binsize

            if covr_end_bin > covr_start_bin:
                self.covrtab_min_bin=min(self.covrtab_min_bin,covr_start_bin)
                target_tab[covr_start_bin]=target_tab[covr_start_bin]+1 if covr_start_bin in target_tab else 1

                if read_end <= self.end:
                    target_tab[covr_end_bin]=target_tab[covr_end_bin]-1 if covr_end_bin in target_tab else -1


################ 2. merge cluster all SV signal. ################
############ 2.0 cluster ############
@dataclass
class Cluster:
    id: str
    svtype: str
    contig: str
    start: int
    end: int
    seed: int
    leads: list
    repeat: bool
    leads_long: list
    varstart: float = None
    varend: float = None
    varsvlen: float = None
    regionstart: float = None
    regionend: float = None


    def compute_metrics(self,max_n=100):
        self.span=self.end-self.start

        n=min(len(self.leads),max_n)
        if n==0:
            self.mean_svlen=0
            self.stdev_start=0
            return

        step=int(len(self.leads)/n)
        if n>1:
            self.mean_svlen=sum(self.leads[i].svlen for i in range(0,len(self.leads),step))/float(n)
            self.stdev_start=statistics.stdev(self.leads[i].ref_start for i in range(0,len(self.leads),step))
        else:
            self.mean_svlen=self.leads[0].svlen
            self.stdev_start=0
    
@dataclass
class SVCallPostprocessingInfo:
    cluster: list

@dataclass
class SVCall:
    contig: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: int
    filter: str
    info: dict

    svtype: str
    svlen: int
    end: int
    genotypes: dict

    precise: bool
    support: int
    rnames: list

    qc: bool
    nm: float
    postprocess: SVCallPostprocessingInfo

    varstart : float = None 
    varend : float = None
    varsvlen : float = None
    regionstart: float = None
    regionend: float = None
    
    fwd: int = None
    rev: int = None

    coverage_upstream: int = None
    coverage_downstream: int = None
    coverage_start: int = None
    coverage_center: int = None
    coverage_end: int = None

    

    def set_info(self,k,v):
        self.info[k]=v

    def get_info(self,k):
        return self.info[k] if k in self.info else None

    def has_info(self,k):
        return k in self.info

    def finalize(self):
        self.postprocess=None

######## 2.0.1 calculate start and end of svlen_bin ########
def calculate_bounds(svtype,ref_start_mode,svlen_mode):
    if svtype=="INS":
        svstart=ref_start_mode
        svend=ref_start_mode
    elif svtype=="DEL":
        svstart=ref_start_mode+svlen_mode
        svend=ref_start_mode
    else:
        svstart=ref_start_mode
        svend=svstart+abs(svlen_mode)
    return svstart,svend

class CoverageCalculator:
    # 把 svcall，filed 加入到 目标bin 对应的 value 中
    @classmethod
    def add_request(cls, svcall, field, pos, requests_for_coverage, config):
        bin = int(pos / config.coverage_binsize) * config.coverage_binsize
        if bin not in requests_for_coverage:
            requests_for_coverage[bin] = []
        requests_for_coverage[bin].append((svcall, field))
    #构造 requests_for_coverage
    @classmethod
    def coverage_build_requests(cls, calls, lead_provider, config):      # calls 是  resplit 处理后输出的 new_cluster 的list,合并完SVlen的结果
        requests_for_coverage = {}
        for svcall in calls:
            start = svcall.pos      #svstart, 众数中位数
            if svcall.svtype == "INS":
                end = start + 1
            else:
                end = svcall.pos + abs(svcall.svlen)        #coverage_binsize=100
            #requests_for_coverage中是 {bin1:(svcall,"coverage_start"),bin2:(svcall,"coverage_center")}
            cls.add_request(svcall, "coverage_start", start - config.coverage_binsize, requests_for_coverage, config)   #向上扩100bp
            cls.add_request(svcall, "coverage_center", int((start + end) / 2), requests_for_coverage, config)
            cls.add_request(svcall, "coverage_end", end + config.coverage_binsize, requests_for_coverage, config)   #coverage_updown_bins=5
            cls.add_request(svcall, "coverage_upstream", start - config.coverage_binsize * config.coverage_updown_bins, requests_for_coverage,  config)
            cls.add_request(svcall, "coverage_downstream", end + config.coverage_binsize * config.coverage_updown_bins, requests_for_coverage,  config)
        return requests_for_coverage
    # 计算平均覆盖度
    @classmethod
    def coverage_fulfill(cls,requests_for_coverage, calls, lead_provider, config):
        if len(requests_for_coverage) == 0:
            return -1, -1

        start_bin = lead_provider.covrtab_min_bin   #start_bin 向内 300bp,整个染色体开始有比对的位置
        end_bin = int(lead_provider.end / config.coverage_binsize) * config.coverage_binsize
        coverage_fwd = 0
        coverage_rev = 0
        coverage_fwd_total = 0
        coverage_rev_total = 0
        n = 0

        for bin in range(start_bin, end_bin + config.coverage_binsize, config.coverage_binsize):    #对 bin 进行遍历
            n += 1

            if bin in lead_provider.covrtab_fwd:
                coverage_fwd += lead_provider.covrtab_fwd[bin]  #计算单个 bin 总的 reads 覆盖度

            if bin in lead_provider.covrtab_rev:
                coverage_rev += lead_provider.covrtab_rev[bin]  #计算单个 bin 总的 reads 覆盖度

            if bin in requests_for_coverage:        #
                coverage_total_curr = coverage_fwd + coverage_rev
                for svcall, field in requests_for_coverage[bin]:
                    setattr(svcall, field, coverage_total_curr)      #coverage 赋值给svcall的filed，比如coverage_start

            coverage_fwd_total += coverage_fwd      # 每个bin的 coverange 累加
            coverage_rev_total += coverage_rev      

        average_coverage_fwd = coverage_fwd_total / float(n) if n > 0 else 0
        average_coverage_rev = coverage_rev_total / float(n) if n > 0 else 0
        return average_coverage_fwd, average_coverage_rev

    @classmethod
    def coverage(cls,calls, lead_provider, config):
        requests_for_coverage = cls.coverage_build_requests(calls, lead_provider, config)
        return cls.coverage_fulfill(requests_for_coverage, calls, lead_provider, config)

############ 2.2 merge sv of same reads  ############
def merge_inner(cluster,threshold):
    read_seq={}
    for ld in cluster.leads:
        if not ld.read_qname in read_seq:
            read_seq[ld.read_qname]=[]
        read_seq[ld.read_qname].append(ld)

    cluster.leads=[]
    for qname in read_seq:
        read_seq[qname].sort(key=lambda k: k.ref_start)
        to_merge=read_seq[qname][0]

        curr_lead=to_merge

        last_ref_end=to_merge.ref_end
        last_qry_end=to_merge.qry_end
        last_ref_start=to_merge.ref_start
        last_qry_start=to_merge.qry_start

        for to_merge in read_seq[qname][1:]:
            merge=(threshold==-1) or ((abs(to_merge.ref_start-last_ref_end) < threshold or abs(to_merge.ref_start-last_ref_start) < threshold) and (abs(to_merge.qry_start-last_qry_end) < threshold or abs(to_merge.qry_start-last_qry_start) < threshold))
            if merge:
                curr_lead.svlen+=to_merge.svlen
                if to_merge.seq==None or curr_lead.seq==None:
                    curr_lead.seq=None
                else:
                    curr_lead.seq+=to_merge.seq
            else:
                cluster.leads.append(curr_lead)
                curr_lead=to_merge
            last_ref_end=to_merge.ref_end
            last_qry_end=to_merge.qry_end
            last_ref_start=to_merge.ref_start
            last_qry_start=to_merge.qry_start

        cluster.leads.append(curr_lead)
    return cluster

############ 2.3 partition to 20 bp bin  ############
def resplit(cluster,prop,binsize,merge_threshold_min,merge_threshold_frac):
    bins_leads={}
    for lead in cluster.leads:
        bin=int(abs(prop(lead))/binsize)*binsize
        if not bin in bins_leads:
            bins_leads[bin]=[lead]
        else:
            bins_leads[bin].append(lead)

    new_clusters=list(sorted(bins_leads.keys()))
    i=1
    while len(new_clusters) > 1 and i < len(new_clusters):
        last_cluster=new_clusters[i-1]
        curr_cluster=new_clusters[i]
        merge_threshold=max(merge_threshold_min,min(curr_cluster,last_cluster)*merge_threshold_frac)
        merge=abs(curr_cluster - last_cluster) <= merge_threshold # svlen差异<=max(35,min(svlen1,svlen2)*0.33)
        if merge:
            bins_leads[new_clusters[i]].extend(bins_leads[new_clusters[i-1]])
            new_clusters.pop(i-1)
            i=max(0,i-2)
        else:
            i+=1

    for cluster_index in new_clusters:
        new_cluster=Cluster(id=cluster.id+f".{cluster_index}",
                            svtype=cluster.svtype,
                            contig=cluster.contig,
                            start=cluster.start,
                            end=cluster.end,
                            seed=cluster.seed,
                            leads=bins_leads[cluster_index],
                            repeat=cluster.repeat,
                            leads_long=cluster.leads_long,
                            varstart=cluster.varstart,
                            varend=cluster.varend,
                            varsvlen=cluster.varsvlen,
                            regionstart=cluster.regionstart,
                            regionend=cluster.regionend,
                        )
        yield new_cluster

############ 2.3.b   ############
def resplit_bnd(cluster,merge_threshold):
    contigs_leads={}
    #Split by mate contig and mate ref start
    for lead in cluster.leads:
        if not lead.bnd_info.mate_contig in contigs_leads:
            contigs_leads[lead.bnd_info.mate_contig]={}

        bin=int(lead.bnd_info.mate_ref_start/merge_threshold)*merge_threshold if merge_threshold>0 else 0
        if not bin in contigs_leads[lead.bnd_info.mate_contig]:
            contigs_leads[lead.bnd_info.mate_contig][bin]=[lead]
        else:
            contigs_leads[lead.bnd_info.mate_contig][bin].append(lead)

    for contig in contigs_leads:
        bins=sorted(contigs_leads[contig])
        curr_leads=[]+contigs_leads[contig][bins[0]]
        last_bin=bins[0]
        for bin in bins[1:]:
            if bin - last_bin <= merge_threshold:
                curr_leads.extend(contigs_leads[contig][bin])
            else:
                if len(curr_leads):
                    new_cluster=Cluster(id=cluster.id+f".CHR2.{contig}.POS2.{bin}",
                                        svtype=cluster.svtype,
                                        contig=cluster.contig,
                                        start=cluster.start,
                                        end=cluster.end,
                                        seed=cluster.seed,
                                        leads=[k for k in curr_leads],
                                        repeat=cluster.repeat,
                                        leads_long=None,
                                        regionstart=cluster.regionstart,
                                        regionend=cluster.regionend,
                                        varstart=cluster.varstart,
                                        varend=cluster.varend,
                                        varsvlen=cluster.varsvlen,
                                        )
                    yield new_cluster
                curr_leads=[]+contigs_leads[contig][bin]
            last_bin=bin
        if len(curr_leads):
            new_cluster=Cluster(id=cluster.id+f".CHR2.{contig}.POS2.{bin}",
                                svtype=cluster.svtype,
                                contig=cluster.contig,
                                start=cluster.start,
                                end=cluster.end,
                                seed=cluster.seed,
                                leads=[k for k in curr_leads],
                                repeat=cluster.repeat,
                                leads_long=None,
                                regionstart=cluster.regionstart,
                                regionend=cluster.regionend,
                                varstart=cluster.varstart,
                                varend=cluster.varend,
                                varsvlen=cluster.varsvlen,
                                )
            curr_leads=[]
            yield new_cluster

def resolve_bnd(svcall,cluster,config):
    mate_contig=util.most_common_top([lead.bnd_info.mate_contig for lead in cluster.leads])
    selected=[lead for lead in cluster.leads if lead.bnd_info.mate_contig==mate_contig]
    mate_ref_start=util.center([lead.bnd_info.mate_ref_start for lead in selected])
    is_first=util.most_common_top([lead.bnd_info.is_first for lead in selected])
    is_reverse=util.most_common_top([lead.bnd_info.is_reverse for lead in selected])
    svcall.alt=(("N" if is_first else "") +
                ("]" if is_reverse else "[" ) +
                f"{mate_contig}:{mate_ref_start}" +
                ("]" if is_reverse else "[" ) +
                ("N" if not is_first else ""))
    svcall.support=len(set(k.read_qname for k in selected))
    cluster.leads=selected
    svcall.bnd_info=SVCallBNDInfo(mate_contig=mate_contig, mate_ref_start=mate_ref_start, is_first=is_first, is_reverse=is_reverse)
    svcall.set_info("CHR2",mate_contig)


############ 2. get_clusted -three_step. ############
def run_script(contig, start, end):
    os.system("/path/to/hap.sh {0} {1} {2}".format(contig, start, end))
def hapasm(regions, max_workers):
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for region in regions:
            contig, start, end = region[:3]
            futures.append(executor.submit(run_script, contig, start, end))
        for future in futures:
            future.result() 

def resolve_pre(svtype,leadtab_provider,config,tr):
    leadtab=leadtab_provider.leadtab[svtype]
    seeds=sorted(leadtab_provider.leadtab[svtype])

    if len(seeds)==0:
        return [],[],[]
    #print(f'len(seeds):{len(seeds)}')
    #Initialize tandem repeat region handling
    if tr != None:
        tr_index=0
        if len(tr)==0:
            tr=None
        else:
            tr_start,tr_end=tr[tr_index]

    clusters = []
    regions = []
    lowclusters = []
    for seed_index,seed in enumerate(seeds):
        if config.dev_call_region != None:
            if seed<config.dev_call_region["start"] or seed>config.dev_call_region["end"]:
                continue

        within_tr=False
        if tr!=None and tr_index < len(tr):
            while tr_end < seed and tr_index+1 < len(tr):
                tr_index+=1
                tr_start,tr_end=tr[tr_index]
            if seed > tr_start and seed < tr_end:
                within_tr=True

        if svtype=="INS":
            leads=[lead for lead in leadtab[seed] if lead.svlen!=None]
            leads_long=[lead for lead in leadtab[seed] if lead.svlen==None]
        else:
            leads=leadtab[seed]
            leads_long=None

        cluster=Cluster(id=f"CL.{svtype}.{leadtab_provider.contig}.{leadtab_provider.start}.{seed_index}",
                        svtype=svtype,
                        contig=leadtab_provider.contig,
                        start=seed,
                        end=seed+config.cluster_binsize,
                        seed=seed,
                        leads=leads,
                        repeat=within_tr or config.repeat,
                        leads_long=leads_long,
                        varstart=util.variance(min(lead.ref_start,lead.ref_end) for lead in leads), #util.variance(lead.ref_start for lead in leads)
                        varend=util.variance(max(lead.ref_start,lead.ref_end) for lead in leads), #util.variance(lead.ref_end for lead in leads)
                        varsvlen=util.variance(lead.svlen for lead in leads),
                        regionstart=min([lead.ref_start for lead in leads if len(leads)!=0]+[lead.ref_end for lead in leads if len(leads)!=0],default=None),
                        regionend=max([lead.ref_start for lead in leads if len(leads)!=0]+[lead.ref_end for lead in leads if len(leads)!=0],default=None),
                        )
        cluster.compute_metrics()
        clusters.append(cluster)
    return clusters,regions,lowclusters

    '''#regions = []
        #clusterss = []
        if len(cluster.leads)==0:
            continue
        if cluster.regionstart and cluster.regionend:
            region_start = cluster.regionstart - config.span
            region_end = cluster.regionend + config.span
            region_svlen = util.mean([lead.svlen for lead in leads]) if len([lead.svlen for lead in leads]) > 0 else 0
            cluster.compute_metrics()
            if (cluster.varstart > config.disstart and cluster.varsvlen > config.dislen):
                lowclusters.append(cluster)
                regions.append((cluster.contig,region_start,region_end,cluster.svtype,region_svlen,"disorder"))
            elif (len(set(k.read_qname for k in cluster.leads)) < config.signal):
                lowclusters.append(cluster)
                regions.append((cluster.contig,region_start,region_end,cluster.svtype,region_svlen,"lowsignal"))
            else:
                #print(f'signal: {len(set(k.read_qname for k in cluster.leads))}') 
                clusters.append(cluster)
        else: 
            #if times < 3:
            #    print(cluster)
            #print('No regionstart or regionend !!!')
            #times+=1
            cluster.compute_metrics()
            clusters.append(cluster)
    #print(f'After filtered disorder :{clusterss}')
    #print(svtype)
    #print(leadtab_provider.contig)
    #print(len(clusters), len(regions), len(lowclusters)) 
    if len(clusters) == 0 and len(regions) == 0 and len(lowclusters) == 0:
        return [],[],[]
    return clusters,regions,lowclusters
    '''

def resolve(clusters,svtype,leadtab_provider,config):
    #Merge clusters
    cluster_count_initial=len(clusters)
    #print(f'Merge svstart, cluster_count_initial:{cluster_count_initial}')
    i=0
    while i < len(clusters) - 1:
        curr_cluster=clusters[i]
        next_cluster=clusters[i+1]
        
        inner_dist=(next_cluster.start - curr_cluster.end)
        outer_dist=(next_cluster.end - curr_cluster.start) #<2.5*(min起始位点标准差)
        merge = inner_dist <= min(curr_cluster.stdev_start, next_cluster.stdev_start) * config.cluster_r  # <min(1000,平均长度*1.5)
        merge = merge or ( (config.repeat or curr_cluster.repeat or next_cluster.repeat) and outer_dist <= min(config.cluster_repeat_h_max, (abs(curr_cluster.mean_svlen)+abs(next_cluster.mean_svlen)) * config.cluster_repeat_h) )
        merge = merge or (svtype=="BND" and inner_dist <= config.cluster_merge_bnd)

        if merge:
            clusters.pop(i+1)
            curr_cluster.leads+=next_cluster.leads
            if svtype=="INS":
                curr_cluster.leads_long+=next_cluster.leads_long
            curr_cluster.end=next_cluster.end
            curr_cluster.repeat=curr_cluster.repeat or next_cluster.repeat
            curr_cluster.compute_metrics()
            i=max(0,i-2)
        i+=1
        
    if config.dev_dump_clusters:
        filename=f"{config.input}.clusters.{svtype}.{leadtab_provider.contig}.{leadtab_provider.start}.{leadtab_provider.end}.bed"
        #print(f"Dumping clusters to {filename}")
        with open(filename,"w") as h:
            for c in clusters:
                info=f"ID={c.id}, #LEADS={len(c.leads)}; "
                for ld in c.leads:
                    info+=f"(ref_start={ld.ref_start},svlen={ld.svlen},source={ld.source}); "
                h.write(f"{c.contig}\t{c.start}\t{c.end}\t\"{info}\"\n")
    #filter disorder and signal
    for cluster in clusters:
        if len(cluster.leads)==0:
            continue

        if svtype == "BND":
            if config.dev_no_resplit:
                yield cluster
            else:
                for new_cluster in resplit_bnd(cluster,merge_threshold=config.bnd_cluster_resplit):
                    yield new_cluster
        else:
            if svtype=="INS" or svtype=="DEL":
                if cluster.repeat:
                    merge_inner_threshold=-1
                else:
                    merge_inner_threshold=config.cluster_merge_pos

                merge_inner(cluster,merge_inner_threshold)	#合并同一qname下的相近SV

            if not config.dev_no_resplit_repeat and not config.dev_no_resplit:
                for new_cluster in resplit(cluster,			# 合并svlen
                                           prop=lambda lead: lead.svlen,
                                           binsize=config.cluster_resplit_binsize,
                                           merge_threshold_min=config.minsvlen,
                                           merge_threshold_frac=config.cluster_merge_len):
                    yield new_cluster
            else:
                yield cluster
#merge svlen_bin，cluster的最后一步，计算svstart
def call_from(cluster,config,keep_qc_fails,task):
    leads=cluster.leads

    svtype=cluster.svtype
    stdev_pos,stdev_len=None,None
    qc=True

    svlen=util.center(v.svlen for v in leads)

    if abs(svlen) < config.minsvlen_screen:
        return

    #Count inline events only once per read, but split events as individual alignments, as in coverage calculation
    #inline_qnames=set(k.read_qname for k in leads if k.source=="INLINE")
    #support=len(inline_qnames)+sum(1 for k in leads if k.source!="INLINE")
    if svtype=="INS" and svlen>=config.long_ins_length:
        support_long_set=set(lead.read_qname for lead in cluster.leads_long)
        support=len(set(k.read_qname for k in leads) | support_long_set)
        support_long=len(support_long_set)
    else:
        support=len(set(k.read_qname for k in leads))   # 有信号的lead的数目
        support_long=0
    ref_start=util.center(v.ref_start for v in leads)
    stdev_pos=util.stdev(util.trim((v.ref_start for v in leads)))
    if svtype!="BND":
        stdev_len=util.stdev(util.trim((v.svlen for v in leads)))
        precise=(stdev_pos+stdev_len < config.precise)
    else:
        precise=(stdev_pos < config.precise)
    svstart,svend=calculate_bounds(svtype,ref_start,svlen)
    qual=int(util.mean(v.mapq for v in leads))

    support_fwd=sum(lead.strand == "+" for lead in leads)
    support_rev=len(leads) - support_fwd

    filter="PASS"
    if config.qc_strand and (support_fwd==0 or support_rev==0):
        filter="STRAND"
        qc=False

    if config.qc_nm:
        nm_mean=util.mean(v.nm for v in leads)
        if nm_mean > config.qc_nm_max:
            filter="NM"
            qc=False
    else:
        nm_mean=-1

    if not keep_qc_fails and not qc:
        return

    svpi=SVCallPostprocessingInfo(cluster=cluster)

    if config.output_rnames or config.snf!=None:
        rnames=list(set(k.read_qname for k in leads))
    else:
        rnames=None

    svcall=SVCall(contig=cluster.contig,
                  pos=svstart,
                  id=f"{svtype}.{task.sv_id:X}S{task.id:X}",
                  ref="N",
                  alt=f"<{svtype}>",
                  qual=qual,
                  filter=filter,
                  info=dict(),
                  svtype=svtype,
                  svlen=svlen,
                  end=svend,
                  genotypes=dict(),
                  precise=precise,
                  support=support,
                  rnames=rnames,
                  postprocess=svpi,
                  qc=qc,
                  nm=nm_mean,
                  fwd=support_fwd,
                  rev=support_rev,
                  varstart=cluster.varstart,
                  varend=cluster.varend,
                  varsvlen=cluster.varsvlen,
                  regionstart=cluster.regionstart,
                  regionend=cluster.regionend,
                  )

    if svtype=="BND":
        resolve_bnd(svcall,cluster,config)
    elif svtype=="INS":
        svcall.set_info("SUPPORT_LONG",support_long)

    if stdev_pos!=None:
        svcall.set_info("STDEV_POS",stdev_pos)
    if stdev_len!=None:
        svcall.set_info("STDEV_LEN",stdev_len)

    #svcall.set_info("CLUSTER_ID",cluster.id)

    task.sv_id+=1

    yield svcall



################ 3. post_processing. ################
############ 3.0 prepare  ############
# rescale_support
def rescale_support(svcall,config):
    if svcall.svtype!="INS" or svcall.svlen < config.long_ins_length:
        return svcall.support
    else:
        base=svcall.support + svcall.get_info("SUPPORT_LONG")
        scale_factor=config.long_ins_rescale_mult*(float(svcall.svlen)/config.long_ins_length)
        return round(base*(config.long_ins_rescale_base+scale_factor))
# phase
def phase_sv(svcall,config):
    reads_phases={lead.read_id[0]: (lead.read_id[1],lead.read_id[2]) for lead in svcall.postprocess.cluster.leads}
    hp_list=util.most_common(hp for hp,ps in reads_phases.values())
    ps_list=util.most_common(ps for hp,ps in reads_phases.values())

    hp_support,hp=hp_list[0]
    ps_support,ps=ps_list[0]
    if hp==None:
        hp="NULL"
    if ps==None:
        ps="NULL"

    other_hp_support=sum(other_supp for other_supp, other_hp in hp_list if other_hp != hp and other_hp != "NULL")
    other_ps_support=sum(other_supp for other_supp, other_ps in ps_list if other_ps != ps and other_ps != "NULL")

    hp_filter="FAIL"
    if hp != "NULL" and hp_support > 0 and float(other_hp_support)/(hp_support+other_hp_support) < config.phase_conflict_threshold:
        hp_filter="PASS"

    ps_filter="FAIL"
    if ps != "NULL" and ps_support > 0 and float(other_ps_support)/(ps_support+other_ps_support) < config.phase_conflict_threshold:
        ps_filter="PASS"

    svcall.set_info("PHASE",f"{hp},{ps},{hp_support},{ps_support},{hp_filter},{ps_filter}")
    return (config.phase_identifiers.index(hp) if hp in config.phase_identifiers else None if hp_filter=="PASS" else None)

# genotype
def genotype_sv(svcall,config,phase):
    normalization_target=250
    hom_ref_p=config.genotype_error
    het_p=(1.0/config.genotype_ploidy) # - config.genotype_error
    hom_var_p=1.0 - config.genotype_error
    coverage=0
    def binomial_probability(k,n,p):
        try:
            #Binomial coef cancels out for likelihood ratios
            #return binomial_coef(n,k) * (p**k) * ((1.0-p)**(n-k))
            return (p**k) * ((1.0-p)**(n-k))
        except OverflowError:
            return 1.0

    def likelihood_ratio(q1,q2):
        if q1/q2>0:
            try:
                return math.log(q1/q2,10)
            except ValueError:
                return 0
        else:
            return 0
    #Count inline events only once per read, but split events as individual alignments, as in coverage calculation
    leads=svcall.postprocess.cluster.leads
    support=rescale_support(svcall,config)

    if svcall.svtype=="INS":
        coverage_list=[svcall.coverage_center]
    else:
        if svcall.svtype=="DUP":
            if False and svcall.coverage_start!=None and svcall.coverage_end!=None:
                if svcall.coverage_start > svcall.coverage_end:
                    coverage_list=[svcall.coverage_end]
                else:
                    coverage_list=[svcall.coverage_start]
            else:
                coverage_list=[svcall.coverage_start,svcall.coverage_end]
            coverage+=round(support*0.75)
        elif svcall.svtype=="INV":
            coverage_list=[svcall.coverage_upstream,svcall.coverage_downstream]
            coverage+=round(support*0.5)
        else:
            coverage_list=[svcall.coverage_start,svcall.coverage_center,svcall.coverage_end]

    coverage_list=[c for c in coverage_list if c!=None and c!=0]
    if len(coverage_list)==0:
        return
    coverage+=round(sum(coverage_list)/len(coverage_list))

    if support > coverage:
        coverage=support

    af=support / float(coverage)

    genotype_p=[((0,0),hom_ref_p),
                ((0,1),het_p),
                ((1,1),hom_var_p)]

    max_lead=max(support,coverage)
    if max_lead>normalization_target:
        norm=normalization_target/float(max_lead)
        normalized_support=round(support*norm)
        normalized_coverage=round(coverage*norm)
    else:
        normalized_support=support
        normalized_coverage=coverage

    genotype_likelihoods=[]
    for gt, p in genotype_p:
        q=binomial_probability(normalized_support,normalized_coverage,p)
        genotype_likelihoods.append((gt,q))
    genotype_likelihoods.sort(key=lambda k: k[1], reverse=True)

    sum_likelihoods=sum(q for gt,q in genotype_likelihoods)
    normalized_likelihoods=[ (gt,(q/sum_likelihoods)) for gt,q in genotype_likelihoods]

    gt1,q1=normalized_likelihoods[0]
    gt2,q2=normalized_likelihoods[1]
    qz=[q for gt,q in normalized_likelihoods if gt==(0,0)][0]
    genotype_z_score = min(60,int((-10) * likelihood_ratio(qz,q1)))
    genotype_quality = min(60,int((-10) * likelihood_ratio(q2,q1)))

    is_long_ins=(svcall.svtype=="INS" and svcall.svlen >= config.long_ins_length)
    if genotype_z_score < config.genotype_min_z_score and not config.non_germline and not is_long_ins:
        if svcall.filter=="PASS":
            svcall.filter="GT"

    if is_long_ins and gt1==(0,0):
        a,b=".","."
    else:
        a,b=gt1
    svcall.genotypes[0]=(a,b,genotype_quality,coverage-support,support,phase)
    svcall.set_info("AF",af)

# get consensus seq
def novel_from_reads(best_lead,other_leads,klen,skip,skip_repetitive,debug=False):
    consensus_min=2
    maxshift=klen
    minspan=0.2
    minalns=0.25
    minident=0.5
    minident_abs=5
    minbestdiff=3
    def iter_kmers(seq,klen,skip):
        for i in range(0,len(seq)-klen,skip):
            yield (i,seq[i:i+klen])
    alignments=[]
    anchors={}
    taboo=set()
    for i, kmer in iter_kmers(best_lead.seq,klen=klen,skip=skip_repetitive):
        if kmer in taboo:
            continue
        if kmer in anchors:
            del anchors[kmer]
            taboo.add(kmer)
            continue
        anchors[kmer]=i

    for leadi,lead in enumerate(other_leads):
        last_i=None
        last_j=None
        conseq=""
        span=0
        for j, kmer in iter_kmers(lead.seq,klen=klen,skip=skip):
            if not kmer in anchors:
                continue
            i=anchors[kmer]
            if abs(i-j) > maxshift:
                continue
            if last_i != None and i <= last_i:
                continue

            if last_i == None:
                if j>0:
                    conseq="-"*i
            else:
                fwd_i=i-last_i
                fwd_j=j-last_j
                if len(conseq)+fwd_j > len(best_lead.seq):
                    fwd_j=len(best_lead.seq)-len(conseq)

                if fwd_i == fwd_j and fwd_j > 0:
                    span+=(j-last_j)
                    m=0
                    for l in range(1,(j-last_j)+1):
                        if lead.seq[last_j+l] == best_lead.seq[last_i+l]:
                            m+=1
                    ident=m/float((j-last_j))
                    if ident >= minident:
                        conseq+=lead.seq[last_j:j][:fwd_j]
                    else:
                        conseq+="-"*(fwd_j)
                else:
                    conseq+="-"*(fwd_j)
            last_i=i
            last_j=j

        if len(conseq) < len(best_lead.seq):
            conseq+="-"*(len(best_lead.seq)-len(conseq))

        conseq_new=[]
        h=0
        while h < len(best_lead.seq):
            if conseq[h]=="-":
                conseq_new.append("-")
                h+=1
            else:
                buffer=[]
                ident=0
                while h < len(best_lead.seq) and conseq[h]!="-":
                    ident+=(best_lead.seq[h]==conseq[h])
                    buffer.append(conseq[h])
                    h+=1
                if ident/float(len(buffer)) > minident and ident>minident_abs:
                    conseq_new.append("".join(buffer))
                else:
                    conseq_new.append("-"*len(buffer))
        conseq="".join(conseq_new)

        if span/float(len(best_lead.seq)) > minspan:
            alignments.append(conseq)

    maxal=1
    for i in range(len(best_lead.seq)):
        maxal=max(maxal,len([best_lead.seq[i]]+[a[i] for a in alignments if not a[i] in "^_"]))
    maxal=float(maxal)

    flattened=""
    for i in range(len(best_lead.seq)):
        al=[a[i] for a in alignments if not a[i]=="-"]
        if len(al) < consensus_min or len(al)/maxal < minalns:
            flattened+=best_lead.seq[i]
        else:
            top=util.most_common([best_lead.seq[i]]+al)
            if len(top)>1 and top[0][0]-top[1][0] >= minbestdiff:
                flattened+=top[0][1]
            else:
                flattened+=best_lead.seq[i]

            #if top[0][0]/float(len(al)+1) < 0.75:
            #    flattened+=best_lead.seq[i]
            #else:
            #    flattened+=top[0][1]

    ##print("FLT",flattened)
    #if debug:
    #    #print("B",best_lead.seq)
    #    for o in other_leads:
    #        #print("O",o.seq)
    #    #print("F",flattened)
    #    #print("=====")
    return flattened

############ 3.1 filter stdev_svstart/end  ############
# 3.1 对 候选sv 进行 stdev_svstart、stdev_svlen的过滤，去掉波动太大的; svlen、coverage的过滤，去掉<50，长DEL中间coverage高，长DUP中间coverage低，INS的coverage<=1的
def qc_sv(svcall,config):
    if config.qc_stdev:
        stdev_pos=svcall.get_info("STDEV_POS")
        if stdev_pos > config.qc_stdev_abs_max:
            svcall.filter="STDEV_POS"
            return False
        if svcall.svtype!="BND" and stdev_pos / abs(svcall.svlen) > 2.0:
            svcall.filter="STDEV_POS"
            return False

        stdev_len = svcall.get_info("STDEV_LEN")
        if stdev_len != None:
            if svcall.svtype != "BND" and stdev_len / abs(svcall.svlen) > 1.0:
                svcall.filter="STDEV_LEN"
                return False
            if stdev_len > config.qc_stdev_abs_max:
                svcall.filter="STDEV_LEN"
                return False

    if abs(svcall.svlen) < config.minsvlen:
        svcall.filter="SVLEN_MIN"
        return False

    #if (svcall.coverage_upstream != None and svcall.coverage_upstream < config.qc_coverage) or (svcall.coverage_downstream != None and svcall.coverage_downstream < config.qc_coverage):
    if svcall.svtype != "DEL" and svcall.svtype != "INS" and (svcall.coverage_center != None and svcall.coverage_center < config.qc_coverage):
        svcall.filter="COV_MIN"
        return False

    if svcall.svtype == "DEL" and config.long_del_length != -1 and abs(svcall.svlen) >= config.long_del_length and not config.non_germline:
        if svcall.coverage_center != None and svcall.coverage_upstream != None and svcall.coverage_downstream != None and svcall.coverage_center > (svcall.coverage_upstream+svcall.coverage_downstream)/2.0 * config.long_del_coverage:
            svcall.filter="COV_CHANGE"
            return False
    elif svcall.svtype=="INS" and ( (svcall.coverage_upstream != None and svcall.coverage_upstream < config.qc_coverage) or (svcall.coverage_downstream != None and svcall.coverage_downstream < config.qc_coverage)):
        svcall.filter="COV_CHANGE"
        return False
    elif svcall.svtype == "DUP" and config.long_dup_length != -1 and abs(svcall.svlen) >= config.long_dup_length and not config.non_germline:
        if svcall.coverage_center != None and svcall.coverage_upstream != None and svcall.coverage_downstream != None and svcall.coverage_center < (svcall.coverage_upstream+svcall.coverage_downstream)/2.0 * config.long_dup_coverage:
            svcall.filter="COV_CHANGE"
            return False

    return True

############ 3.2 filter supp  ############
# 3.2 对 候选sv 进行 supp >= minsupp 的过滤
def qc_sv_support(svcall, coverage_global, config):
    #print(f'config.minsupport: {config.minsupport}')
    support = rescale_support(svcall, config)
    coverage_list = [svcall.coverage_upstream, svcall.coverage_downstream]
    coverage_list = [c for c in coverage_list if c is not None and c != 0]

    if not coverage_list:
        coverage_list = [svcall.coverage_start, svcall.coverage_center, svcall.coverage_end]
        coverage_list = [c for c in coverage_list if c is not None and c != 0]

    if not coverage_list:
        coverage_regional = coverage_global
    else:
        coverage_regional = round(sum(coverage_list) / len(coverage_list))
        if coverage_regional == 0:
            coverage_regional = coverage_global

    if config.minsupport == "auto":
        coverage = (coverage_regional * config.minsupport_auto_regional_coverage_weight + 
                    coverage_global * (1.0 - config.minsupport_auto_regional_coverage_weight))
        min_support = round(config.minsupport_auto_base + config.minsupport_auto_mult * coverage)
        minsupp = round(2+0.2*coverage) #int(config.minsupp)  ###max(int(config.minsupp),round(config.minsupport_auto_base + config.minsupport_auto_mult * coverage))
        if support < min_support:
            return False
        elif support < minsupp:
            svcall.filter = "SUPPORT_MIN"
            #print('candidate !!!')
            return False
        
#        if coverage_regional < config.depth:
#            svcall.filter = "LOW_DEPTH"
#            return False
#        elif support < min_supports:
#            svcall.filter = "SUPPORT_MIN"
#            return False
        
    else:
        if support < int(config.minsupport):
            svcall.filter = "SUPPORT_MIN_noauto"
            return False
    return True

############ 3.3 get GT  ############
#调用 genotype_sv 计算基因型，对INS 生成最佳 consensus seq
def annotate_sv(svcall,config):
    if config.phase:
        phase=phase_sv(svcall,config)
    else:
        phase=None

    genotype_sv(svcall,config,phase)

    if svcall.svtype=="INS" and not config.symbolic:
        merged_leads=[l for l in svcall.postprocess.cluster.leads if l.seq!=None]

        if len(merged_leads):
            best_lead=merged_leads[0]
            best_index=0
            best_diff=abs(len(best_lead.seq) - svcall.svlen) + abs(best_lead.ref_start - svcall.pos)*1.5
            for i,ld in enumerate(merged_leads):
                if i==0:
                    continue
                curr_diff=abs(len(ld.seq) - svcall.svlen) + abs(ld.ref_start - svcall.pos)*1.5
                if curr_diff < best_diff:
                    best_lead=ld
                    best_index=i
                    best_diff=curr_diff

            merged_leads.pop(best_index)
            #merged_leads_new=list()

            #for lead in merged_leads:
            #    curr_lendiff=abs(len(ld.seq) - len(best_lead.seq)) + abs(ld.ref_start - best_lead.ref_start)
            #    #if curr_lendiff < 25:
            #    merged_leads_new.append(lead)
            #merged_leads=merged_leads_new


            if len(merged_leads) >= config.consensus_min_reads and not config.no_consensus:
                klen=config.consensus_kmer_len
                skip=config.consensus_kmer_skip_base+int(len(best_lead.seq)*config.consensus_kmer_skip_seqlen_mult)
                skip_repetitive=skip

                svcall.alt=novel_from_reads(best_lead,merged_leads,klen=klen,skip=skip,skip_repetitive=skip_repetitive)
            else:
                svcall.alt=best_lead.seq

############ 3.4 filter GT/nm  ############
def qc_sv_post_annotate(svcall,config):
    if (len(svcall.genotypes)==0 or (svcall.genotypes[0][0]!="." and svcall.genotypes[0][0]+svcall.genotypes[0][1]<2)) and (svcall.coverage_center != None and svcall.coverage_center < config.qc_coverage):
        svcall.filter="COV_MIN"
        return False

    return True

############ * NGS_DEL  * ############

class GetValid:
    def __init__(self, bamfile, mapq_min, alen_min):
        self.bamfile = bamfile
        self.mapq_min = mapq_min
        self.alen_min = alen_min
        
    def get_abnormal_value(self,contig,start,end,raw_svtype,abnormal,abnormal_ratio):
        tlen = defaultdict(list)
        bam = pysam.AlignmentFile(self.bamfile, 'rb', require_index=True)
        for read in bam.fetch(contig,start,end,until_eof=False):
            alen = read.query_alignment_length
            if read.mapping_quality < self.mapq_min or read.is_secondary or alen < self.alen_min:
                continue
            if read.is_paired: 
                length = abs(read.template_length)
                if length not in tlen[read.query_name]:
                    tlen[read.query_name].append(length)
        bam.close()
        all_values = [value for sublist in tlen.values() for value in list(sublist)]
        if not all_values:
            return None
        mean_value = util.mean(all_values)
        std_dev =  util.stdev(all_values)
        threshold_upper = mean_value + 3 * std_dev
        threshold_lower = mean_value - 3 * std_dev
        mid_outliers = {k: v for k, v in tlen.items() if any(x > threshold_upper or x < threshold_lower for x in v)}
        outliers = defaultdict(list)
        for readid,temlen in mid_outliers.items():
            for lens in temlen:
                svlen = lens - mean_value 
                svtype = "DEL" #if svlen > 0 else "INS"
                if raw_svtype != svtype :
                    #print('svtype error')
                    continue
                else:
                    if abs(svlen) < 300:
                        #print('svtype or svlen Error')
                        #print(f"svlen:{svlen}")
                        continue
                    else:
                        outliers[readid].append((svlen,svtype))
        if len(outliers) < abnormal and len(outliers)/len(tlen) < abnormal_ratio:
            print('low frequency Error')
            outliers = None
        return outliers

    def split_cigar(self,cigar):
        numl = []
        num = 0
        cigarstrl = []
        for i in cigar:
            if i.isdigit():
                num = num*10+int(i)
            else:
                numl.append(num)
                num = 0
                cigarstrl.append(i)
        return numl,cigarstrl

    def get_soft_clip_ratio(self,contig,start,end):
        bam=pysam.AlignmentFile(self.bamfile, 'rb', require_index=True)
        allread = 0
        sclipread = 0
        for read in bam.fetch(contig,start,end,until_eof=False):
            alen = read.query_alignment_length
            if read.mapping_quality < self.mapq_min or read.is_secondary or alen < self.alen_min:
                continue
            allread += 1
            if read.cigarstring:
                numl,cigarl = self.split_cigar(read.cigarstring)
            else:
                continue
            if cigarl[0] == "S" or cigarl[-1] == 'S':
                sclipread += 1
        bam.close()
        sclipratio = (sclipread / allread) if allread else 0
        if sclipratio > 0.5 :
            print('soft_clip!!!!')
        return sclipratio

    def verify_abnormal_value(self, candidates, abnormal, abnormal_ratio):
        reserve_candidate=[]
        for candidate in candidates:
            svtype = candidate.svtype
            if svtype != "DEL":
                print("not DEL")
                continue
            elif abs(candidate.svlen) < 300:
                continue
            contig = candidate.contig
            region_start = candidate.pos - 100 - abs(candidate.svlen)
            region_end = candidate.pos + 100 + abs(candidate.svlen)
            
            #if self.get_soft_clip_ratio(contig, region_start, region_end) > 0.5:
            if self.get_abnormal_value(contig, region_start, region_end, svtype, abnormal, abnormal_ratio):
                reserve_candidate.append(candidate)
        return reserve_candidate
        
    def verify_softclip(self, candidates):
        reserve_candidate=[]
        for candidate in candidates:
            svtype = candidate.svtype
            contig = candidate.contig
            region_start = candidate.pos - 100 - abs(candidate.svlen)
            region_end = candidate.pos + 100 + abs(candidate.svlen)
            
            if self.get_soft_clip_ratio(contig, region_start, region_end) > 0.5:
            #if self.get_abnormal_value(contig, region_start, region_end, svtype, abnormal, abnormal_ratio):
                reserve_candidate.append(candidate)
                
        return reserve_candidate

############ 4 output vcf format ############
class VCF:
    def __init__(self,config,handle):
        self.config=config
        self.handle=handle
        self.call_count=0
        self.info_order=["SVTYPE","SVLEN","END","SUPPORT","RNAMES","COVERAGE","STRAND"]
        if config.qc_nm:
            self.info_order.append("NM")

        self.default_genotype=config.genotype_none

        if config.mode=="combine":
            self.genotype_format=config.genotype_format+":ID"
            self.default_genotype+=tuple(["NULL"])
        else:
            self.genotype_format=config.genotype_format

        self.reference_handle=None
        if config.mode == "call_sample":
            if config.sample_id is None:
                # config.sample_id,_=os.path.splitext(os.path.basename(config.input))
                config.sample_ids_vcf = [(0, "SAMPLE")]
            else:
                config.sample_ids_vcf = [(0, config.sample_id)]


    def format_info(self, k, v):
        if isinstance(v, float):
            return f"{k}={v:.3f}"
        elif isinstance(v, list):
            return f"{k}={','.join(v)}"
        else:
            return f"{k}={v}"

    def format_genotype(self, gt):
        if len(gt) == 6:
            a, b, qual, dr, dv, ps = gt
            if ps is not None and (a, b) == (0, 1):
                gt_sep = "|"
                if ps == 1:
                    a, b = b, a
            else:
                gt_sep = "/"
            return f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}"
        else:
            a, b, qual, dr, dv, ps, id = gt
            if ps is not None and (a, b) == (0, 1):
                gt_sep = "|"
                if ps == 1:
                    a, b = b, a
            else:
                gt_sep = "/"
            return f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}:{id}"

    def open_reference(self):
        if self.config.reference == None:
            return
        if not os.path.exists(self.config.reference+".fai") and not os.path.exists(self.config.reference+".gzi"):
            print(f"Info: Fasta index for {self.config.reference} not found. Generating with pysam.faidx (this may take a while)")
            pysam.faidx(self.config.reference)
        self.reference_handle=pysam.FastaFile(self.config.reference)

    def write_header(self,contigs_lengths):
        self.write_header_line("fileformat=VCFv4.2")
        self.write_header_line(f"source={self.config.version}_{self.config.build}")
        self.write_header_line('command="'+self.config.command+'"')
        self.write_header_line('fileDate="'+self.config.start_date+'"')
        for contig,contig_len in contigs_lengths:
            self.write_header_line(f"contig=<ID={contig},length={contig_len}>")

        self.write_header_line('ALT=<ID=INS,Description="Insertion">')
        self.write_header_line('ALT=<ID=DEL,Description="Deletion">')
        self.write_header_line('ALT=<ID=DUP,Description="Duplication">')
        self.write_header_line('ALT=<ID=INV,Description="Inversion">')
        self.write_header_line('ALT=<ID=BND,Description="Breakend; Translocation">')

        self.write_header_line('FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        self.write_header_line('FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
        self.write_header_line('FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">')
        self.write_header_line('FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">')
        self.write_header_line('FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">')

        self.write_header_line('FILTER=<ID=PASS,Description="All filters passed">')
        self.write_header_line('FILTER=<ID=GT,Description="Genotype filter">')
        self.write_header_line('FILTER=<ID=SUPPORT_MIN,Description="Minimum read support filter">')
        self.write_header_line('FILTER=<ID=STDEV_POS,Description="SV Breakpoint standard deviation filter">')
        self.write_header_line('FILTER=<ID=STDEV_LEN,Description="SV length standard deviation filter">')
        self.write_header_line('FILTER=<ID=COV_MIN,Description="Minimum coverage filter">')
        self.write_header_line('FILTER=<ID=COV_CHANGE,Description="Coverage change filter">')
        self.write_header_line('FILTER=<ID=STRAND,Description="Strand support filter">')
        self.write_header_line('FILTER=<ID=SVLEN_MIN,Description="SV length filter">')
        self.write_header_line('FILTER=<ID=NM,Description="Alignment noise level filter">')
        self.write_header_line('INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">')
        self.write_header_line('INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">')
        self.write_header_line('INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">')
        self.write_header_line('INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">')
        self.write_header_line('INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">')
        self.write_header_line('INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">')
        self.write_header_line('INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description="Number of reads supporting an INS/DEL SV (non-split events only)">')
        self.write_header_line('INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Number of soft-clipped reads putatively supporting the long insertion SV">')
        self.write_header_line('INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">')
        self.write_header_line('INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">')
        self.write_header_line('INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">')
        self.write_header_line('INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Coverages near upstream, start, center, end, downstream of structural variation">')
        self.write_header_line('INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">')
        self.write_header_line('INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">')
        self.write_header_line('INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="List of read support for all samples">')
        self.write_header_line('INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description="Number of reads that support the generated insertion (INS) consensus sequence">')
        self.write_header_line('INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads (if enabled with --output-rnames)">')
        self.write_header_line('INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">')
        self.write_header_line('INFO=<ID=NM,Number=.,Type=Float,Description="Mean number of query alignment length adjusted mismatches of supporting reads">')
        self.write_header_line('INFO=<ID=PHASE,Number=.,Type=String,Description="Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPORT,HAPLOTYPE_FILTER,PHASESET_FILTER">')

        samples_header="\t".join(sample_id for internal_id,sample_id in self.config.sample_ids_vcf)
        self.write_raw(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples_header}")

    def write_raw(self,text,endl="\n"):
        if self.config.vcf_output_bgz:
            self.handle.write(text.encode())
            self.handle.write(endl.encode())
        else:
            self.handle.write(text)
            self.handle.write(endl)

    def write_header_line(self,text):
        self.write_raw("##"+text)

    def write_call(self,call):
        #pysam coordinates are 0-based, VCF 1-based - therefore +1 is added to call.pos
        end=call.end + 1
        pos=call.pos + 1

        #Determine genotypes columns
        ac=0 #Allele count
        supvec=[]
        sample_genotypes=[]
        for internal_id, sample_id in self.config.sample_ids_vcf:
            if internal_id in call.genotypes and call.genotypes[internal_id]!=None:
                gt_curr=call.genotypes[internal_id]
                sample_genotypes.append(self.format_genotype(gt_curr))
                if gt_curr[0]!="." and gt_curr[4]>0: #Not non-genotype and has supporting reads
                    ac+=sum(call.genotypes[internal_id][:2])
                    supp="1"
                else:
                    supp="0"
            else:
                sample_genotypes.append(self.format_genotype(self.default_genotype))
                supp="0"
            supvec.append(supp)

        if len(self.config.sample_ids_vcf) > 1:
            call.set_info("AC",ac)
            call.set_info("SUPP_VEC","".join(supvec))
            if ac==0:
                call.filter="GT"

        #Output core SV attributes
        infos={"SVTYPE": call.svtype,
               "SVLEN": call.svlen,
               "END": end,
               "SUPPORT": call.support,
               "RNAMES": call.rnames if self.config.output_rnames else None,
               "COVERAGE": f"{call.coverage_upstream},{call.coverage_start},{call.coverage_center},{call.coverage_end},{call.coverage_downstream}",
               "STRAND": ("+" if call.fwd>0 else "") + ("-" if call.rev>0 else ""),
               "NM": call.nm}

        if call.svtype=="BND":
            infos["SVLEN"]=None
            infos["END"]=None

        infos_ordered=["PRECISE" if call.precise else "IMPRECISE"]
        infos_ordered.extend(self.format_info(k,infos[k]) for k in self.info_order if infos[k]!=None)
        info_str=";".join(infos_ordered)

        #Output call specific additional information
        for k in sorted(call.info):
            if call.info[k]==None:
                continue
            info_str+=";" + self.format_info(k,call.info[k])

        #if call.id==None:
        #    call.id=f"Sniffles2.{call.svtype}.{self.call_count+1:06}"

        #Resolve DEL sequence
        if not self.config.symbolic and call.svtype=="DEL" and self.reference_handle != None and abs(call.svlen) <= self.config.max_del_seq_len:
            try:
                call.ref=self.reference_handle.fetch(call.contig,call.pos,call.pos-call.svlen)
                call.alt="N"
            except KeyError:
                call.ref="N"
                call.alt=f"<{call.svtype}>"
            except ValueError:
                call.ref="N"
                call.alt=f"<{call.svtype}>"

        if self.config.symbolic:
            call.ref="N"
            call.alt=f"<{call.svtype}>"

        call.qual=max(0,min(60,call.qual))

        self.write_raw("\t".join(str(v) for v in [call.contig,pos,self.config.id_prefix+call.id,call.ref,call.alt,call.qual,call.filter,info_str,self.genotype_format]+sample_genotypes))
        self.call_count+=1

    def read_svs_iter(self):
        self.header_str=""
        line_index=0
        for line in self.handle:
            try:
                if isinstance(line,bytes):
                    line=line.decode("utf-8")
                line_index+=1
                line_strip=line.strip()
                if line_strip=="" or line_strip[0]=="#":
                    if line_strip[0]=="#":
                        self.header_str+=line_strip+"\n"
                    continue
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO=line.split("\t")[:8]
                info_dict={}
                for info_item in INFO.split(";"):
                    if "=" in info_item:
                        key,value=info_item.split("=")
                    else:
                        key,value=info_item,True
                    info_dict[key]=value
                call=SVCall(contig=CHROM,
                               pos=int(POS)-1,
                               id=line_index,
                               ref=REF,
                               alt=ALT,
                               qual=QUAL,
                               filter=FILTER,
                               info=info_dict,
                               svtype=None,
                               svlen=None,
                               end=None,
                               rnames=None,
                               qc=True,
                               postprocess=None,
                               genotypes=None,
                               precise=None,
                               support=0,
                               fwd=0,
                               rev=0,
                               nm=-1,
                               )
                if len(call.alt)>len(call.ref):
                    call.svtype="INS"
                    call.svlen=len(call.alt)
                    call.end=call.pos
                else:
                    call.svtype="DEL"
                    call.svlen=-len(call.ref)
                    call.end=call.pos+call.svlen

                if "SVTYPE" in info_dict:
                    call.svtype=info_dict["SVTYPE"]
                    if call.svtype=="TRA":
                        call.svtype="BND"

                if "SVLEN" in info_dict:
                    call.svlen=int(info_dict["SVLEN"])
                if "END" in info_dict:
                    call.end=int(info_dict["END"])

                if call.svtype=="BND":
                    bnd_parts=call.alt.replace("]","[").split("[")
                    if len(bnd_parts)>2:
                        mate_contig,mate_ref_start=bnd_parts[1].split(":")
                        call.bnd_info=SVCallBNDInfo(mate_contig=mate_contig, mate_ref_start=int(mate_ref_start), is_first=(call.alt[0]=="N"), is_reverse=("]" in call.alt))
                    else:
                        raise ValueError("BND ALT not formatted according to VCF 4.2 specifications")

                call.raw_vcf_line=line_strip
                call.raw_vcf_line_index=line_index
                yield call
            except Exception as e:
                util.fatal_error(f"Error parsing input VCF: Line {line_index}: {e}")

    def rewrite_genotype(self,svcall):
        parts_no_gt=svcall.raw_vcf_line.split("\t")[:8]
        gt_format=self.config.genotype_format
        if svcall.genotype_match_sv != None:
            if len(svcall.genotype_match_sv.genotypes)>0:
                gt=svcall.genotype_match_sv.genotypes[0]
            else:
                gt=svcall.genotypes[0]
        else:
            gt=svcall.genotypes[0] #0/0 or ./. depending on whether input SV had coverage in sample
        #parts=parts_no_gt + [gt_format,vcf.format_genotype(gt)]
        #gt_vcf=svcall.raw_vcf_line.split("\t")[9].split(":")[0]
        #parts= parts_no_gt + [gt_vcf] + [gt_format,vcf.format_genotype(gt)]
        parts=parts_no_gt + [gt_format,self.format_genotype(gt)]
        #parts[7]="NA"
        #parts[3]=f"REF_{len(parts[3])}"
        #parts[4]=f"ALT_{len(parts[4])}"
        self.write_raw("\t".join(parts))

    def rewrite_header_genotype(self,orig_header):
        header_lines=orig_header.split("\n")
        header_lines.insert(1,'##genotypeFileDate="'+self.config.start_date+'"')
        header_lines.insert(1,'##genotypeCommand="'+self.config.command+'"')
        header_lines.insert(1,f"##genotypeSource={self.config.version}_{self.config.build}")
        self.write_raw("\n".join(header_lines),endl="")

    def close(self):
        self.handle.close()

def get_VCF(config,outpath,svcalls,contig_lengths):
    config.vcf_output_bgz=False
    if os.path.exists(outpath) and not config.allow_overwrite:
        raise ValueError(f"Output file '{outpath}' already exists! Use --allow-overwrite to ignore this check and overwrite.")
    else:
        if not Path(outpath).parent.exists():
            Path(outpath).parent.mkdir(parents=True)
        if config.vcf_output_bgz:
            if not config.sort:
                raise ValueError(".gz (bgzip) output is only supported with sorting enabled")
            vcf_handle=pysam.BGZFile(outpath,"w")
        else:
            vcf_handle=open(outpath,"w")    
        vcf_out=VCF(config,vcf_handle)    
    if config.mode=="call_sample" or config.mode=="combine":
        if config.reference!=None:
            print(f"Opening for reading: {config.reference}")
        vcf_out.open_reference()    
    print(f"Opening for writing: {outpath} ")
    if config.mode!="genotype_vcf" and outpath != None:
        vcf_out.write_header(contig_lengths)
    if svcalls:
        for svcall in svcalls:
            vcf_out.write_call(svcall)



################ 整合应用  ################
@dataclass
class Task:
    id: int
    sv_id: int
    contig: str
    start: int
    end: int
    assigned_process_id: int=None
    lead_provider: object=None
    bam: object=None
    tandem_repeats: list=None
    genotype_svs: list=None

    def build_leadtab(self,config):             #get_SV_signal
        assert(self.lead_provider==None)

        self.bam=pysam.AlignmentFile(config.input, 'rb', require_index=True)
        total_mapped=self.bam.mapped
        if total_mapped==0:
            #Total mapped returns 0 for CRAM files
            config.task_read_id_offset_mult=10**9
        else:
            #BAM file
            config.task_read_id_offset_mult=10**math.ceil(math.log(total_mapped)+1)
        self.lead_provider=LeadProvider(config,self.id*config.task_read_id_offset_mult)
        externals=self.lead_provider.build_leadtab(self.contig,self.start,self.end,self.bam)
        #print(f'DEL num :{self.lead_provider.leadtab}')
        return externals,self.lead_provider.read_count

    def call_candidates(self,keep_qc_fails,config):
        candidates=[]
        lowcandidates=[]
        for svtype in TYPES:
            clusters,regions,lowclusters=resolve_pre(svtype,self.lead_provider,config,self.tandem_repeats)
            """if len(result) == 3:
                clusters,regions,lowclusters=result
            else:
                print('error in resolve_pre')
                print(svtype)
                print(self.lead_provider.leadtab[svtype])
                continue"""
                
            
            '''
            if regions:														# !!!限速步骤!!!! 
                hapasm(regions, config.hapcpu)
            '''
            for svcluster in resolve(clusters,svtype,self.lead_provider,config):     #输入的是某种SVtype下所有的
                for svcall in call_from(svcluster,config,keep_qc_fails,self):
                    candidates.append(svcall)
            for fcluster in resolve(lowclusters,svtype,self.lead_provider,config):
                for fsvcall in call_from(fcluster,config,keep_qc_fails,self):
                    lowcandidates.append(fsvcall)
                    
        self.coverage_average_fwd,self.coverage_average_rev=CoverageCalculator.coverage(candidates,self.lead_provider,config)
        self.coverage_average_total=self.coverage_average_fwd+self.coverage_average_rev
        #过滤
        return candidates,regions,lowcandidates

    def finalize_candidates(self,candidates,keep_qc_fails,config):
        passed=[]
        tregions=[]
        supp2=[]
        for svcall in candidates:
            svcall.qc=svcall.qc and qc_sv(svcall,config)
            if not keep_qc_fails and not svcall.qc:
                tregions.append((svcall.contig,svcall.regionstart-config.span if svcall.regionstart else 0,svcall.regionend+config.span if svcall.regionend else 0,svcall.svtype,svcall.svlen,"qc_sv"))
                continue
            svcall.qc=svcall.qc and qc_sv_support(svcall,self.coverage_average_total,config)
            if not keep_qc_fails and not svcall.qc:
                tregions.append((svcall.contig,svcall.regionstart-config.span if svcall.regionstart else 0,svcall.regionend+config.span if svcall.regionend else 0,svcall.svtype,svcall.svlen,"sv_support"))
                
                if svcall.filter == 'SUPPORT_MIN' :
                    #print('get supp < 3 ')
                    annotate_sv(svcall,config)
                    supp2.append(svcall)
                continue

            annotate_sv(svcall,config)
            
            svcall.qc=svcall.qc and qc_sv_post_annotate(svcall,config)
            if not keep_qc_fails and not svcall.qc:
                continue

            #svcall.finalize() #Remove internal information (not written to output) before sending to mainthread for VCF writing
            passed.append(svcall)
        '''
        if tregions:     # !!!限速步骤!!!! 
            hapasm(tregions, config.hapcpu)
        '''
        return passed,tregions,supp2


def get_svcall(config, contig, start, end,tandem_repeats):
    if config.snf is not None or config.no_qc:          #default no_qc=False
        qc = False
    else:
        qc = True
    # # 调用 build_leadtab 方法来设置 self.contig, self.start 和 self.end
    call = Task(id=1, sv_id=0, contig=contig, start=start, end=end, tandem_repeats=tandem_repeats)
    externals, read_count = call.build_leadtab(config)
    svcandidates,regions,lowcandidates = call.call_candidates(qc, config)
    #print(f'svcandidates:{len(svcandidates)}')
    svcalls,tregions,lowsvcalls = call.finalize_candidates(svcandidates, not qc, config)
    #lowsvcalls,lowtregions = call.finalize_candidates(lowcandidates, not qc, config)
    allregion = regions + tregions #+ lowtregions 
    #### 在这里加上验证步骤，把有二代支持的再放且不参与后面的过滤。不管其他的候选，只管第一次过滤出来的。
    #validSV = GetValid(config.NGSbam, config.NGSmapq_min, config.NGSalen_min)
    #valid_candidates = validSV.verify_abnormal_value(lowsvcalls, config.abnormal, config.abnormal_ratio)
    #valid_candidates = validSV.verify_softclip(lowsvcalls)
    
    if not config.no_qc:    #no_qc=False
        svcalls = [s for s in svcalls if s.qc]
    if config.sort:
        svcalls = sorted(svcalls, key=lambda svcall: svcall.pos)    #根据pos排序
    return svcalls,allregion,lowsvcalls

def process_svcall(outfile,svcalls):
    with open(outfile,'w') as f:
        for svcall in svcalls:
            #print('\t'.join(map(str, [svcall.contig, svcall.pos, svcall.svtype, svcall.svlen,svcall.support, svcall.genotypes])))
            f.write('\t'.join(map(str, [svcall.contig, svcall.pos + 1,svcall.svtype, svcall.svlen, svcall.varstart, svcall.varend, svcall.varsvlen,svcall.regionstart,svcall.regionend])) + '\n')

def getallregions(outfile,allregions):
    with open(outfile,'w') as f:
        for region in allregions:
            #print(region[0])
            f.write('\t'.join(map(str, region)) + '\n')


def process_contig(args: Tuple[str, str, Dict[str, Any], Any]) -> List[SVCall] :
    contig_str, input_file, tandem_repeats, config = args
    print(contig_str)
    svcalls = []
    allregions = []
    valid_candidates = []
    bam = pysam.AlignmentFile(input_file, 'rb', require_index=True)
    startpos = 0
    contig_length = bam.get_reference_length(contig_str)
    endpos = contig_length - 1
    svcall,allregion,valid_candidate = get_svcall(config, contig_str, startpos, endpos, tandem_repeats)
    svcalls.extend(svcall)
    allregions.extend(allregion)
    valid_candidates.extend(valid_candidate)
    bam.close()
    return svcalls,allregions,valid_candidates

def process_bam_file(input_file: str, output_file: str, config, contig_tandem_repeats: Dict[str, Any],cpu):
    bam = pysam.AlignmentFile(input_file, 'rb', require_index=True)
    contigs = [str(contig.contig) for contig in bam.get_index_statistics()]
    contig_lengths = [(str(contig.contig),bam.get_reference_length(str(contig.contig))) for contig in bam.get_index_statistics()]
    args = [(contig, input_file, contig_tandem_repeats.get(contig, None), config) for contig in contigs]
    bam.close()
    analysis_start_time=time.time()
    pool = Pool(processes=cpu)
    results = pool.map(process_contig, args)
    pool.close()
    pool.join()
    svcalls = []
    allregions = []
    valid_candidates = []
    finalSVcalls = []
    for svcall,allregion,valid_candidate in results:
        svcalls.extend(svcall)
        #print(len(svcalls))
        allregions.extend(allregion)
        valid_candidates.extend(valid_candidate)
    print(f"Took {time.time()-analysis_start_time:.2f}s.")
    processtime=time.time()
    process_svcall(output_file, svcalls)
    print(f'valid_candidates:{len(valid_candidates)}')
    #getallregions(config.outregion, allregions) 
    #finalSVcalls = svcalls.extend(valid_candidates)
    get_VCF(config,config.filtervcf,valid_candidates,contig_lengths)
    #get_VCF(config,config.vcf,svcalls,contig_lengths)
    get_VCF(config,config.vcf,svcalls,contig_lengths)
    if config.truvari:
        os.system("/data/VarTeam/linwei/02_scripts/truvari.sh sniffles207_signal{0}_filtered {1}".format(config.signal,config.outpath))
    print(f"Took {time.time()-processtime:.2f}s to process svcall.")


if __name__ == '__main__':
    import sys
    import time
    import argparse
    
    parses=argparse.ArgumentParser(description='get_SV')
    parses.add_argument("-i","--input", type=str, help="For single-sample calling: A coordinate-sorted and indexed .bam/.cram (BAM/CRAM format) file containing aligned reads. - OR - For multi-sample calling: Multiple .snf files (generated before by running Sniffles2 for individual samples with --snf)")
    parses.add_argument("-N","--NGSbam", type=str, help="NGS bam")
    parses.add_argument("-v","--vcf", metavar="OUT.vcf", type=str, help="VCF output filename to write the called and refined SVs to. If the given filename ends with .gz, the VCF file will be automatically bgzipped and a .tbi index built for it.")
    parses.add_argument("-f","--filtervcf", metavar="OUT.vcf", type=str, help="Filt supp3 VCF output filename to write the called and refined SVs to. If the given filename ends with .gz, the VCF file will be automatically bgzipped and a .tbi index built for it.")
    parses.add_argument("--snf", metavar="OUT.snf", type=str, help="Sniffles2 file (.snf) output filename to store candidates for later multi-sample calling", required=False)
    parses.add_argument("-r","--reference", metavar="reference.fasta", type=str, help="(Optional) Reference sequence the reads were aligned against. To enable output of deletion SV sequences, this parameter must be set.")
    parses.add_argument("--tandem-repeats", metavar="IN.bed", type=str, help="(Optional) Input .bed file containing tandem repeat annotations for the reference genome.", default=None)
    parses.add_argument('-l','--location',dest='location',action='store',help='chr_strat_end_types')
    parses.add_argument('-s','--disstart',dest='disstart',type=int,action='store',help='disorder_start',default=1000000000000)
    parses.add_argument('-e','--dislen',dest='dislen',type=int,action='store',help='disorder_svlen',default=100000000000)
    parses.add_argument('-signal','--signal',dest='signal',type=int,action='store',help='signal nums',default=0)
    parses.add_argument('-d','--depth',dest='depth',type=int,action='store',help='real(no) depth',default=0)
    parses.add_argument('-o',"--outpath",  type=str, help="outpath")
    parses.add_argument('-t',"--threads",  type=int, help="threads", default=8)
    parses.add_argument("--minsupport", metavar="auto", type=str, help="Minimum number of supporting reads for a SV to be reported (default: automatically choose based on coverage)", default="auto")
    parses.add_argument("--minsupp", metavar="auto", type=int, help="Minimum number of supporting reads for a SV to be reported (default: automatically choose based on coverage)", default=0)
    
    args=parses.parse_args()

    config = get_config()
    config.input = args.input
    config.NGSbam = args.NGSbam
    config.snf = args.snf
    config.reference = args.reference
    config.disstart = args.disstart
    config.dislen = args.dislen
    config.signal = args.signal
    config.depth = args.depth
    config.truvari = args.truvari
    config.vcf = args.outpath + '/HighQ.vcf'
    config.filtervcf = args.outpath + '/Candidate.vcf'
    location=args.location
    config.minsupport = args.minsupport
    config.minsupp = args.minsupp
    
    config.outpath = args.outpath 
    config.outregion = args.outpath+'/region_filtered_signal{0}_to_hapasm.info'.format(config.signal)

    tandem_repeats_file = args.tandem_repeats
    if tandem_repeats_file is None:
        contig_tandem_repeats = {}
    else:
        contig_tandem_repeats = util.load_tandem_repeats(tandem_repeats_file,500)
    cpu = args.threads
    outfile = args.outpath+'/variance_signal{0}_filteredSV.info'.format(config.signal)
    print(config.input)
    process_bam_file(config.input, outfile, config, contig_tandem_repeats, cpu)
