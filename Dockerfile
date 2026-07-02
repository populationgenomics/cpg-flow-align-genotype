FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg2-1

ENV PYTHONDONTWRITEBYTECODE=1
ENV VERSION=0.4.19

WORKDIR /align_genotype

COPY src src/
COPY LICENSE pyproject.toml README.md ./

# Install samtools and pip packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential samtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    pip install .

# set bash as the shell of choice
SHELL ["/bin/bash", "-c"]

# install Rust & Cargo to facilitate streamed running of dupblaster
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y && \
    echo 'source $HOME/.cargo/env' >> $HOME/.bashrc && \
    source $HOME/.bashrc \
    cargo install dupblaster
