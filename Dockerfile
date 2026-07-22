FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg2-1 AS build

# set bash as the shell of choice
SHELL ["/bin/bash", "-c"]

# install build-essential for compiling
# install Rust & Cargo to facilitate install of dupblaster
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        bzip2 \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        wget \
        zlib1g-dev && \
    curl https://sh.rustup.rs -sSf | bash -s -- -y && \
    echo 'source $HOME/.cargo/env' >> $HOME/.bashrc && \
    source $HOME/.bashrc && \
    cargo install dupblaster

ENV SAMTOOLS_VERSION=1.23.1

	# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make && \
    make DESTDIR=/samtools_install install && \
    make -C htslib-${SAMTOOLS_VERSION} DESTDIR=/samtools_install install

FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg2-1

# copy in dupblaster
COPY --from=build /root/.cargo/bin/dupblaster /usr/local/bin/dupblaster

# copy in samtools
COPY --from=build /samtools_install/usr/local/bin/* /usr/local/bin/

# update, install packages required for samtools, clear lists
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libncurses5-dev \
        libssl1.1 \
        zlib1g && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

ENV PYTHONDONTWRITEBYTECODE=1
ENV VERSION=0.5.3

WORKDIR /align_genotype

COPY src src/
COPY LICENSE pyproject.toml README.md ./
RUN pip install .
