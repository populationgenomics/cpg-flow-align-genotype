import click
from cpg_flow.resources import STANDARD
from cpg_utils import Path, to_path
from cpg_utils.config import get_gcp_project, image_path
from cpg_utils.hail_batch import Batch, command, get_batch
from google.cloud import storage
from loguru import logger

TEN_GB = 10 * 1024 ** 3

def find_files_to_index(
    path: Path,
) -> dict[str, list[tuple[str, str]]]:
    """Finds bam and cram files and queues them for indexing if they are not already indexed."""
    client = storage.Client()
    bucket_name = path.as_uri().split('/')[2]
    prefix = '/'.join(path.as_uri().split('/')[3:])
    bucket = client.bucket(bucket_name, user_project=get_gcp_project())

    # Find all files with the given extensions
    files_to_index = {'samtools': [], 'bcftools': []}
    blobs = bucket.list_blobs(prefix=prefix)
    blob_names = [f'gs://{bucket_name}/{blob.name}' for blob in blobs]
    for blob_name in blob_names:
        # Check if the index file exists
        if blob_name.endswith('.bam'):
            index_blob_name = f'{blob_name}.bai'
        elif blob_name.endswith('.cram'):
            index_blob_name = f'{blob_name}.crai'
        elif blob_name.endswith('.vcf.gz'):
            index_blob_name = f'{blob_name}.tbi'
        else:
            continue
        if index_blob_name in blob_names:
            continue
        if blob_name.endswith(('.bam', '.cram')):
            files_to_index['samtools'].append((blob_name, index_blob_name))
        elif blob_name.endswith('.vcf.gz'):
            files_to_index['bcftools'].append((blob_name, index_blob_name))

    return files_to_index


def index_with_samtools(
    b: Batch,
    input_files: list[tuple[str, str]],
    disk_size: str | None,
):
    """
    Index bam or cram files using samtools.
    """
    for file_to_index_path, index_file_path in input_files:
        input_file = b.read_input(file_to_index_path)
        j = b.new_bash_job(f'samtools index {file_to_index_path}')
        j.image(image_path('samtools'))
        # Set resource requirements
        storage_gb = int(disk_size) if disk_size else 20 if to_path(file_to_index_path).stat().st_size < TEN_GB else 100
        nthreads = 8
        res = STANDARD.set_resources(
            j,
            ncpu=nthreads,
            storage_gb=storage_gb,
        )
        cmd = f'samtools index -@ {res.get_nthreads() - 1} {input_file} -o {j.output}'
        j.command(command(cmd, monitor_space=True))
        b.write_output(j.output, index_file_path)


def index_with_bcftools(
    b: Batch,
    input_files: list[tuple[str, str]],
    disk_size: str | None,
):
    """
    Index VCF.GZ files using bcftools.
    """
    for file_to_index_path, index_file_path in input_files:
        input_file = b.read_input(file_to_index_path)
        j = b.new_bash_job(f'bcftools index {file_to_index_path}')
        j.image(image_path('bcftools'))
        # Set resource requirements
        storage_gb = int(disk_size) if disk_size else 20 if to_path(file_to_index_path).stat().st_size < TEN_GB else 100
        nthreads = 8
        res = STANDARD.set_resources(
            j,
            ncpu=nthreads,
            storage_gb=storage_gb,
        )
        cmd = f'bcftools index -@ {res.get_nthreads() - 1} {input_file} -o {j.output}'
        j.command(command(cmd, monitor_space=True))
        b.write_output(j.output, index_file_path)


@click.command()
@click.option(
    '--input-path',
    help='Path to the input files',
)
@click.option('--disk-size', help='Specify the disk size in gb')
def main(input_path: str, disk_size: str | None):
    # Find all bam, cram, and vcf.gz files that need to be indexed
    files_to_index = find_files_to_index(to_path(input_path))
    # Index the files
    logger.info(f'Found {len(files_to_index["samtools"])} BAM/CRAM files to index')
    logger.info(f'Found {len(files_to_index["bcftools"])} VCF.GZ files to index')
    index_with_samtools(get_batch(), files_to_index['samtools'], disk_size)
    index_with_bcftools(get_batch(), files_to_index['bcftools'], disk_size)
    get_batch().run()


if __name__ == '__main__':
    main()
