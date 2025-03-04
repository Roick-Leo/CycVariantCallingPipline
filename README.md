# freesia_v1
![License](https://img.shields.io/badge/license-GPLv3-blue.svg)<br>
Contact: Lei He<br>
Email: helei1@genomic.cn

## Introduction
freesia_v1 is a germline small variant and structural variant calling program for hybrid data (LRS & SRS).

## Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Features
- 🌟 SNV calling
- ⚡ SV calling

## Installation
### Option 1. Docker pre-built image:
A pre-built docker image is available [here](crpi-y7bygzzrulh9scmj.cn-hangzhou.personal.cr.aliyuncs.com/roick/freesia_v1:v1.0.0.0). With it you can run freesia_v1 using a single command.<br>
<br>
Caution: Absolute path is needed for both INPUT_FILE and OUTPUT_DIR.

```bash
docker run \
    -v /path/to/SRS_R1.fq.gz:/path/to/SRS_R1.fq.gz \
    -v /path/to/SRS_R2.fq.gz:/path/to/SRS_R2.fq.gz \
    -v /path/to/LRS.fq.gz:/path/to/LRS.fq.gz \
    -v /path/to/ref_dir:/path/to/ref_dir \
    -v /path/to/freesia_v1.json:/app/freesia_v1.json \
    -v /path/to/freesia_v1_new.wdl:/app/freesia_v1.wdl \
    roick/freesia_v1:v1.0.0.0 \
    java -jar /app/cromwell.jar run /app/freesia_v1.wdl -i /app/freesia_v1.json
```

### Option 2. Docker Dockerfile
```bash
# clone freesia_v1
git clone https://github.com/Roick-Leo/freesia_v1.git
cd freesia_v1

# build a docker image named roick/freesia_v1:v1.0.0.0
docker build -f ./Dockerfile -t roick/freesia_v1:v1.0.0.0 .
```

## Usage
### Step1. update the input json
```json
{
    "freesia.cyclone_fastq_fn": "/path/to/LRS.fq.gz",
    "freesia.DNB_read1_fn": "/path/to/SRS_R1.fq.gz",
    "freesia.DNB_read2_fn": "/path/to/SRS_R2.fq.gz",
    "freesia.ref_dir": "/path/to/ref_dir",
    "freesia.ref_prefix": "hs37d5.fa",
    "freesia.ref_version": "37",
    "freesia.threads": 30,
    "freesia.sample": "Sample_ID",
    "freesia.out_dir": "/path/to/output_dir"
}
```
### Step2. run 
```bash
docker run \
    -v /path/to/SRS_R1.fq.gz:/path/to/SRS_R1.fq.gz \
    -v /path/to/SRS_R2.fq.gz:/path/to/SRS_R2.fq.gz \
    -v /path/to/LRS.fq.gz:/path/to/LRS.fq.gz \
    -v /path/to/ref_dir:/path/to/ref_dir \
    -v /path/to/freesia_v1.json:/app/freesia_v1.json \
    -v /path/to/freesia_v1_new.wdl:/app/freesia_v1.wdl \
    roick/freesia_v1:v1.0.0.0 \
    java -jar /app/cromwell.jar run /app/freesia_v1.wdl -i /app/freesia_v1.json
```
## Contributing
- He Lei (helei1@genomics.cn)
- Xin Huang (huangxin7@genomics.cn)
- Wei Lin (linwei4@genomics.cn)

## License
Freesia2 is distributed under the GPLv3 license. Source code will we released shortly.