#!/bin/sh

# Usage:
# vcf_vep_anno.sh vep_data/input/deanno/{vcfgz_prefix}.consensus.variant-calls.deanno.vcf vep_data/output/{vcfgz_prefix}.consensus.variant-calls.anno_vep.vcf

source ~/.zshrc

INPUT=$1
OUTPUT=$2
INPUT_VEP=$( echo ${INPUT} | awk '{ sub("vep_data/", ""); print $0 }' )
OUTPUT_VEP=$( echo ${OUTPUT} | awk '{ sub("vep_data/", ""); print $0 }' )
SHCWD=$( pwd )

docker run --rm -v ${SHCWD}/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep \
    vep \
    --cache \
    --offline \
    --input_file ${INPUT_VEP} \
    --output_file ${OUTPUT_VEP} \
    --format vcf \
    --vcf \
    --custom custom/clinvar_20221217.vcf.gz,ClinVar,vcf,exact,0,ID,ALLELEID,CLNDISDB,CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,CLNVC,CLNVCSO,GENEINFO,MC,ORIGIN,RS


## For Snakefile:
# AWK_CMD = r"""{ sub("vep_data/", ""); print $0 }"""
# shell: """echo {input} | awk {AWK_CMD:q}"""
# input_vep=re.sub("vep_data/", "", "vep_data/input/deanno/{vcfgz_prefix}.consensus.variant-calls.deanno.vcf")

