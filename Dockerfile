FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg2-1 AS build

# set bash as the shell of choice
SHELL ["/bin/bash", "-c"]

# install build-essential for compiling
# install Rust & Cargo to facilitate install of dupblaster
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential && \
    curl https://sh.rustup.rs -sSf | bash -s -- -y && \
    echo 'source $HOME/.cargo/env' >> $HOME/.bashrc && \
    source $HOME/.bashrc && \
    cargo install dupblaster

FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg2-1

COPY --from=build /root/.cargo/bin/dupblaster /usr/local/bin/dupblaster

ENV PYTHONDONTWRITEBYTECODE=1
ENV VERSION=0.4.19

# Install samtools and pip packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends samtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /align_genotype

COPY src src/
COPY LICENSE pyproject.toml README.md ./
RUN pip install .
