FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN sed -i s@/archive.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list
RUN sed -i s@/security.ubuntu.com/@/mirrors.aliyun.com/@g /etc/apt/sources.list
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get clean -y
RUN apt-get update && apt-get install -y \
    wget bzip2 ca-certificates git gnupg \
    lsb-release \
    openjdk-11-jre \
    sudo \
    locales \
    unzip curl \
    build-essential time \
    zlib1g-dev libncurses5-dev libboost-all-dev libbz2-dev libcurl4-openssl-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN ln -fs /usr/share/zoneinfo/Asia/Shanghai /etc/localtime && dpkg-reconfigure --frontend noninteractive tzdata

WORKDIR /app
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
RUN curl -L https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 | tar -jxvf -
ENV PATH="/app/minimap2-2.28_x64-linux/:$PATH"
ENV PATH="/app/manta-1.6.0.centos6_x86_64/bin:$PATH"

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p /app/conda
RUN rm miniconda.sh
ENV PATH=/app/conda/bin:$PATH
RUN echo ". /app/conda/etc/profile.d/conda.sh" >> ~/.bashrc

RUN wget https://github.com/broadinstitute/cromwell/releases/download/88/cromwell-88.jar -O /app/cromwell.jar
ENV PATH="/app/cromwell.jar:$PATH"

RUN echo "channels:" > /app/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/" >> /app/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/" >> /app/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2/" >> /app/conda/.condarc && \
    echo "  - defaults" >> /app/conda/.condarc && \
    echo "show_channel_urls: true" >> /app/conda/.condarc

# 将Conda添加到系统路径
ENV PATH="/app/conda/bin:$PATH"

RUN git clone --recursive https://github.com/bwa-mem2/mm2-fast.git mm2-fast && \
    cd mm2-fast && \
    make clean && make
RUN git clone --recursive https://github.com/bwa-mem2/bwa-mem2 && \
    cd bwa-mem2 && \
    make clean && make
ENV PATH="/app/bwa-mem2:$PATH"
RUN curl -L https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 | tar -jxvf - && \
    cd samtools-1.21 && \
    ./configure && make && make install
RUN curl -L https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 | tar -jxvf - && \
    cd bcftools-1.21 && \
    ./configure && make && make install
ENV PATH="/app/samtools-1.21:$PATH"
ENV PATH="/app/bcftools-1.21:$PATH"

RUN conda create -n SNV
RUN conda create -n SV
RUN conda create -n PY2

SHELL ["conda", "run", "-n", "SNV", "/bin/bash", "-c"]
RUN conda install -c bioconda -c conda-forge clair3=1.0.4 python=3.9.0 -y
RUN conda install -c conda-forge boost -y
WORKDIR /app/conda/envs/SNV/bin/preprocess/realign
RUN g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp
RUN g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp -I /app/conda/envs/SNV/include -L /app/conda/envs/SNV/lib

SHELL ["conda", "run", "-n", "PY2", "/bin/bash", "-c"]
RUN conda install -c bioconda -c conda-forge python=2.7 hap.py -y

SHELL ["conda", "run", "-n", "SV", "/bin/bash", "-c"]
RUN conda install -c bioconda -c conda-forge python=3.10.15 sniffles=2.0.7 -y
SHELL ["/bin/bash", "-c"]

RUN git clone https://github.com/Roick-Leo/CycVariantCallingPipline.git /app
ENV PATH="/app/freesia_v1:$PATH"
WORKDIR /workflows
EXPOSE 8000
RUN apt-get purge -y openjdk-11-jre && \
    apt-get autoremove -y && apt-get autoclean -y
RUN apt-get update && apt-get install -y software-properties-common && add-apt-repository ppa:openjdk-r/ppa && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && apt-get update && \
    apt-get install -y openjdk-17-jdk
CMD ["bash"]