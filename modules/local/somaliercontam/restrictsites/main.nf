process SOMALIERCONTAM_RESTRICTSITES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/df/df323a03567e752b8e7a9954f7692f9074f67dbd8127d2cdf563937eb591e1e0/data'
        : 'community.wave.seqera.io/library/python_numpy_pandas_matplotlib-base:3781e1b411a92e2b'}"

    input:
    tuple val(meta), path(sites)
    tuple val(meta2), path(bed)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('somalier_contam'), eval('somalier_contam.py --version | sed "s/.* //"'), emit: versions_somalier_contam, topic: versions
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.panel"
    """
    somalier_contam.py restrict-sites \\
        --sites ${sites} \\
        --bed ${bed} \\
        --output ${prefix}.vcf \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.panel"
    """
    touch ${prefix}.vcf
    """
}
