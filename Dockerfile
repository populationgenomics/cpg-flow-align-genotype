FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.134.cpg2-1
FROM australia-southeast1-docker.pkg.dev/cpg-common/images-dev/cpg_flow:5455de219db190a9a7aa59920b486d06685efeee

ENV PYTHONDONTWRITEBYTECODE=1

WORKDIR /align_genotype

COPY src src/
COPY LICENSE pyproject.toml README.md ./

# pip install but don't retain the cache files
RUN pip install --no-cache-dir .
