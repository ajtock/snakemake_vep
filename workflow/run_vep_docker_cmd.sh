#!/bin/sh

source ~/.zshrc

CWD=$(pwd)
echo ${CWD}

docker run --rm -v ${CWD}/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep \
    vep \
    --cache \
    --offline \
    --input_file input/homo_sapiens_GRCh38.vcf \
    --output_file output/test_ouput_GRCh38.vcf \
    --format vcf \
    --vcf \
    --custom custom/clinvar_20221217.vcf.gz,ClinVar,vcf,exact,0,ID,ALLELEID,CLNDISDB,CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,CLNVC,CLNVCSO,GENEINFO,MC,ORIGIN,RS
