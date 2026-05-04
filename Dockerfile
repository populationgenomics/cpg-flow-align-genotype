FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.138.cpg1-1

ENV PYTHONDONTWRITEBYTECODE=1
ENV VERSION=0.4.9

WORKDIR /align_genotype

COPY src src/
COPY LICENSE pyproject.toml README.md ./

# Install samtools and pip packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends samtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir . && \
    pip install --no-cache-dir pysam
