#!/bin/sh

# Works in Snakemake shell directive:
#awk "{{ FS="\\t"; OFS="\\t" }} /^#+/ {{ print \$0 }}"
# Does not work in Snakemake shell directive (gsub pattern replacement seems to be the problem):
#awk "{{ FS="\\t"; OFS="\\t" }} /^[^#]/ {{ gsub("chr","",\$1); print \$0 }}" {input} > {output.main}
# Hence need shell script vcf_deanno.sh to pass to shell directive of rule vcf_deanno

# Usage:
# vcf_deanno.sh vep_data/input/{vcfgz_prefix}.consensus.variant-calls.vcf.gz

VCFGZ=$1
VCF=$( echo ${VCFGZ} | awk '{ sub(".gz", ""); print $0 }' )
VCF_HEADER=${VCF}_header
VCF_MAIN=${VCF}_main

gunzip --keep --stdout ${VCFGZ} > ${VCF} && \
awk '{ FS="\t"; OFS="\t" } /^#+/ { print $0 }' ${VCF} > ${VCF_HEADER} && \
awk '{ FS="\t"; OFS="\t" } /^[^#]/ { gsub("chr", "", $1); gsub(";ANN=[^;]*;", ";", $8); gsub("TopEffect[^;]*;", "" ,$8); print $0 }' ${VCF} \
    > ${VCF_MAIN} && \
cat ${VCF_HEADER} ${VCF_MAIN}
