#!/usr/bin/env python3

# Author: Andy Tock
# Date: 03/01/2023

# Usage:
# conda activate 6run_analysis
# ./combine_reshape_vcf.py -d /Users/ajtock/dnanexus -s "run*_EU" -m sample_manifest/ori02_ts_sample_manifest_2022_09_26_EW.csv -c clinvar/clinvar_result
# conda deactivate

import os
import glob
import argparse
import re
import pandas as pd
import gc


def create_parser():
    """
    Capture input command-line arguments
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--inDir", type=str, default="/Users/ajtock/dnanexus/snakemake_vep/workflow/vep_data/output",
                        help="Path to the directory containing VEP-annotated VCF files. Default: Users/ajtock/dnanexus/snakemake_vep/workflow/vep_data/output")
    parser.add_argument("-m", "--sampleManifest", type=str, default="/Users/ajtock/dnanexus/snakemake_vep/workflow/resources/manifest/ori02_ts_sample_manifest_2022_09_26_EW.csv",
                        help="CSV of sample manifest. Default: /Users/ajtock/dnanexus/snakemake_vep/workflow/resources/manifest/ori02_ts_sample_manifest_2022_09_26_EW.csv")
    parser.add_argument("--mode", default="client")
    parser.add_argument("--port", default=50301)
    parser.add_argument("--host", default="127.0.0.1")
    return parser

parser = create_parser().parse_args()
print(parser)


gene_list = ['ABL1', 'AKT1', 'ALK', 'APC', 'ATM', 'BRAF', 'CDH1', 'CDKN2A',
             'CSF1R', 'CTNNB1', 'EGFR', 'ERBB2', 'ERBB4', 'EZH2', 'FBXW7',
             'FGFR1', 'FGFR2', 'FGFR3', 'FLT3', 'GNA11', 'GNAQ', 'GNAS', 'HNF1A',
             'IDH1', 'IDH2', 'JAK2', 'JAK3', 'KDR', 'KIT', 'KRAS', 'MET',
             'MLH1', 'MPL', 'NOTCH1', 'NPM1', 'NRAS', 'PDGFRA', 'PIK3CA', 'PTEN',
             'PTPN11', 'RB1', 'RET', 'SMAD4', 'SMARCB1', 'SMO', 'SRC', 'STK11',
             'TP53', 'VHL'] #HRAS (No single gene variants found on ClinVar)
print("Number of genes in gene_list: %s" % len(gene_list))


def cat_vcf_files(indir, info_cols):
    """
    Load and concatenate MUT files from specified directories.
    :param indir: Directory containing VEP-annotated VCF files
    :param info_cols: dictionary of field:dtype for INFO fields to
    store as distinct columns.
    :return: DataFrame of concatenated MUT files
    """
    #indir=parser.inDir
    #info_cols={
    #           "END":"int",
    #           "HIAF":"float",
    #           "HICNT":"int",
    #           "HICOV":"int",
    #           "MQ":"float",
    #           "MSI":"float",
    #           "MSILEN":"float",
    #           "NM":"float",
    #           "PMEAN":"float",
    #           "QUAL":"float",
    #           "TYPE":"str",
    #           "SAMPLE":"str",
    #           "CSQ":"str"
    #          }
    # End args
    file_list = glob.glob(indir + "/*calls.anno_vep.vcf", recursive=True)
    header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT".split()
    vcf = pd.DataFrame()
    for file_path in file_list:
        vcf_x = pd.read_csv(file_path, sep="\t", comment="#", names=header)
        # Concatenate VCFs
        vcf = pd.concat(objs=[vcf, vcf_x],
                        axis=0,
                        ignore_index=True)

    # Convert INFO into a dictionary of subfields
    vcf["INFO"] = vcf["INFO"].str.split(";") \
        .apply(lambda x: dict([y.split("=") for y in x]))
    # Add given INFO subfields as columns
    if info_cols is not None:
        for field, dtype in info_cols.items():
            vcf[field] = vcf["INFO"].apply(lambda x: x.get(field, None))
            vcf[field] = vcf[field].astype(dtype)

    # Convert CSQ into a list of subfields
    vcf["CSQ_ANN_LIST"] = vcf["CSQ"].str.split(",")
    vcf["CSQ_LIST"] = vcf["CSQ_ANN_LIST"] \
        .apply(lambda x: [y.split("|") for y in x])
    vcf["Consequence"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[1] for y in x]))).keys())))
    vcf["IMPACT"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[2] for y in x]))).keys())))
    vcf["SYMBOL"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[3] for y in x]))).keys())))
    vcf["Gene"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[4] for y in x]))).keys())))
    vcf["Feature_type"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[5] for y in x]))).keys())))
    vcf["Feature"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[6] for y in x]))).keys())))
    vcf["BIOTYPE"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[7] for y in x]))).keys())))

    # Convert GT into list and add columns containing list elements
    vcf["GT"] = vcf["GT"].str.split(":")
    vcf["GENOTYPE"] = [x[0] for x in vcf["GT"]]
    vcf["REF_DEPTH"] = [int(x[1].split(",")[0]) for x in vcf["GT"]]
    vcf["ALT_DEPTH"] = [int(x[1].split(",")[1]) for x in vcf["GT"]]
    vcf["VARIANT_DEPTH"] = [int(x[5]) for x in vcf["GT"]]
    vcf["TOTAL_DEPTH"] = [int(x[3]) for x in vcf["GT"]]
    vcf["TOTAL_DEPTH_LESS_NO_CALLS"] = vcf["ALT_DEPTH"] + vcf["REF_DEPTH"]
    vcf["NO_CALLS"] = vcf["TOTAL_DEPTH"] - vcf["TOTAL_DEPTH_LESS_NO_CALLS"]
    vcf["vaf"] = vcf["ALT_DEPTH"] / (vcf["ALT_DEPTH"] + vcf["REF_DEPTH"])

    # Remove unnecessary columns
    vcf.drop(columns=["CSQ_ANN_LIST","CSQ_LIST"], inplace=True)

    # Write to TSV
    vcf.to_csv("results/concat_sample_VCFs.tsv",
               na_rep="NaN", sep="\t", header=True, index=False)

    return vcf


# Reproduce Canonical SPDI
# seq_id_map = {
#     'ABL1': 'NC_000009.12', 'AKT1': 'NC_000014.9', 'ALK': 'NC_000002.12',
#     'APC': 'NC_000005.10', 'ATM': 'NC_000011.10', 'BRAF': 'NC_000007.14',
#     'CDH1': 'NC_000016.10', 'CDKN2A': 'NC_000009.12', 'CSF1R': 'NC_000005.10',
#     'CTNNB1': 'NC_000003.12', 'EGFR': 'NC_000007.14', 'ERBB2': 'NC_000017.11',
#     'ERBB4': 'NC_000002.12', 'EZH2': 'NC_000007.14', 'FBXW7': 'NC_000004.12',
#     'FGFR1': 'NC_000008.11', 'FGFR2': 'NC_000010.11', 'FGFR3': 'NC_000004.12',
#     'FLT3': 'NC_000013.11', 'GNA11': 'NC_000019.10', 'GNAQ': 'NC_000009.12',
#     'GNAS': 'NC_000020.11', 'HNF1A': 'NC_000012.12', #'HRAS': '',
#     'IDH1': 'NC_000002.12', 'IDH2': 'NC_000015.10', 'JAK2': 'NC_000009.12',
#     'JAK3': 'NC_000019.10', 'KDR': 'NC_000004.12', 'KIT': 'NC_000004.12',
#     'KRAS': 'NC_000012.12', 'MET': 'NC_000007.14', 'MLH1': 'NC_000003.12',
#     'MPL': 'NC_000001.11', 'NOTCH1': 'NC_000009.12', 'NPM1': 'NC_000005.10',
#     'NRAS': 'NC_000001.11', 'PDGFRA': 'NC_000004.12', 'PIK3CA': 'NC_000003.12',
#     'PTEN': 'NC_000010.11', 'PTPN11': 'NC_000012.12', 'RB1': 'NC_000013.11',
#     'RET': 'NC_000010.11', 'SMAD4': 'NC_000018.10', 'SMARCB1': 'NC_000022.11',
#     'SMO': 'NC_000007.14', 'SRC': 'NC_000020.11', 'STK11': 'NC_000019.10',
#     'TP53': 'NC_000017.11', 'VHL': 'NC_000003.12'}
# sequence_id = seq_id_map[gene]
#     gene_data['canonical_id'] = (sequence_id + ':'
#                                  + gene_data.start.astype(str) + ':'
#                                  + gene_data.ref + ':' + gene_data.alt)
# https://www.ncbi.nlm.nih.gov/variation/notation/
# Control?: EGI-02-006-189


def merge_mut_mani(vcf_DF, mani_csv):
    """
    Extend concatenated VCF DataFrame with additional information from
    sample manifest CSV, including run, "Subject ID", "cat" (categorical response variable;
    "CRC+", "Polyps+" or "Control"), and "Mean On-Target Duplex Depth" to be used for
    filtering (>= 1000).
    Output filtered VAF values as matrices.
    :param vcf_DF: Concatenated MUT DataFrame
    :param mani_csv: path to sample manifest CSV file
    """
    #vcf_DF = cat_DF
    #mani_csv = parser.sampleManifest
    # End args
    mani_DF = pd.read_csv(mani_csv)
    merge_DF = pd.merge(left=vcf_DF, right=mani_DF,
                        how="inner", left_on="SAMPLE", right_on="Library Name")
    merge_DF["locus"] = merge_DF["CHROM"].astype(str) + "_" + \
                        merge_DF["POS"].astype(str) + "_" + \
                        merge_DF["END"].astype(str)
    merge_DF_dp = merge_DF.loc[
        (merge_DF["Mean On-Target Duplex Depth"] >= 1000) & \
        (merge_DF["ALT_DEPTH"] >= 3) & \
        (merge_DF["FILTER"] == "PASS")
    ]
    #
    file_list = glob.glob(cvdir + "/*.txt", recursive=True)
    cv_DF_concat = pd.DataFrame()
    for file_path in file_list:
        cv_DF = pd.read_csv(file_path, sep="\t")
        cv_DF_concat = pd.concat(objs=[cv_DF_concat, cv_DF],
                                 axis=0,
                                 ignore_index=True)
    #
    variant = [re.sub(".+:", "", x) for x in list(cv_DF_concat["Name"])]
    variant_dna = [re.sub(" .+", "", x) for x in variant]
    variant_dna2 = [re.sub(" \\(p\\..+\\)", "", x) for x in variant]
    featureID = [re.sub("\\(.+", "", x) for x in list(cv_DF_concat["Name"])]
    featureID = [re.sub(":.+", "", x) for x in featureID]
    clinical_significance = [re.sub("\\(Last reviewed:.+", "", x) \
                             for x in list(cv_DF_concat["Clinical significance (Last reviewed)"])]
    cv_DF_concat["variant_dna"] = variant_dna
    cv_DF_concat["variant_dna2"] = variant_dna2
    cv_DF_concat["featureID"] = featureID
    cv_DF_concat["clinical_significance"] = clinical_significance
    cv_DF_concat_pathogenic = cv_DF_concat[cv_DF_concat["clinical_significance"].str.contains("pathogenic", case=False)]
    cv_DF_concat_pathogenic_not_conflicting = cv_DF_concat_pathogenic[
        -cv_DF_concat_pathogenic["clinical_significance"].str.contains("Conflicting", case=False)
    ]
    #
    merge_DF_dp_patho = merge_DF_dp[merge_DF_dp["TopEffectGeneName"].isin(list(cv_DF_concat_pathogenic_not_conflicting["Gene(s)"]))]
    # not_LOW
    merge_DF_dp_patho_CRC_not_LOW = merge_DF_dp_patho.loc[
        (merge_DF_dp_patho["cat"] == "CRC+") & \
        (merge_DF_dp_patho["TopEffectImpact"] != "LOW")
    ]
    merge_DF_dp_patho_Control = merge_DF_dp_patho.loc[
        (merge_DF_dp_patho["cat"] == "Control")
    ]
    merge_DF_dp_patho_Polyps = merge_DF_dp_patho.loc[
        (merge_DF_dp_patho["cat"] == "Polyps+")
    ]
    #
    vaf_DF_dp_patho_CRC_not_LOW = merge_DF_dp_patho_CRC_not_LOW.pivot_table(index=["locus", "contig", "start", "end", "TopEffectGeneName", "TopEffectFeatureId"],
                                                                            columns=["Sample ID", "Subject ID", "cat", "run"],
                                                                            values="vaf",
                                                                            fill_value="NaN")
    vaf_DF_dp_patho_Control = merge_DF_dp_patho_Control.pivot_table(index=["locus", "contig", "start", "end", "TopEffectGeneName", "TopEffectFeatureId"],
                                                                    columns=["Sample ID", "Subject ID", "cat", "run"],
                                                                    values="vaf",
                                                                    fill_value="NaN")
    vaf_DF_dp_patho_Polyps = merge_DF_dp_patho_Polyps.pivot_table(index=["locus", "contig", "start", "end", "TopEffectGeneName", "TopEffectFeatureId"],
                                                                  columns=["Sample ID", "Subject ID", "cat", "run"],
                                                                  values="vaf",
                                                                  fill_value="NaN")
    #
    vaf_DF_dp_patho_CRC_not_LOW_merge_Control = pd.merge(left=vaf_DF_dp_patho_CRC_not_LOW,
                                                         right=vaf_DF_dp_patho_Control,
                                                         how="left",
                                                         on=["locus", "contig", "start", "end", "TopEffectGeneName", "TopEffectFeatureId"])
    vaf_DF_dp_patho_CRC_not_LOW_merge_Control_Polyps = pd.merge(left=vaf_DF_dp_patho_CRC_not_LOW_merge_Control,
                                                                right=vaf_DF_dp_patho_Polyps,
                                                                how="left",
                                                                on=["locus", "contig", "start", "end", "TopEffectGeneName", "TopEffectFeatureId"])

    vaf_DF_dp_patho_CRC_not_LOW_merge_Control.to_csv("vaf_DF_dp_patho_CRC_not_LOW_merge_Control.tsv",
                                                     na_rep="NaN", sep="\t", header=True, index=True)
    vaf_DF_dp_patho_CRC_not_LOW_merge_Control_Polyps.to_csv("vaf_DF_dp_patho_CRC_not_LOW_merge_Control_Polyps.tsv",
                                                            na_rep="NaN", sep="\t", header=True, index=True)
    #
    # HIGH
    merge_DF_dp_patho_CRC_HIGH = merge_DF_dp_patho.loc[
        (merge_DF_dp_patho["cat"] == "CRC+") & \
        (merge_DF_dp_patho["TopEffectImpact"] == "HIGH")
    ]
    #
    vaf_DF_dp_patho_CRC_HIGH = merge_DF_dp_patho_CRC_HIGH.pivot_table(index=["locus", "contig", "start", "end", "TopEffectGeneName", "TopEffectFeatureId"],
                                                                      columns=["Sample ID", "Subject ID", "cat", "run"],
                                                                      values="vaf",
                                                                      fill_value="NaN")
    #
    vaf_DF_dp_patho_CRC_HIGH_merge_Control = pd.merge(left=vaf_DF_dp_patho_CRC_HIGH,
                                                      right=vaf_DF_dp_patho_Control,
                                                      how="left",
                                                      on=["locus", "contig", "start", "end", "TopEffectGeneName", "TopEffectFeatureId"])
    vaf_DF_dp_patho_CRC_HIGH_merge_Control_Polyps = pd.merge(left=vaf_DF_dp_patho_CRC_HIGH_merge_Control,
                                                             right=vaf_DF_dp_patho_Polyps,
                                                             how="left",
                                                             on=["locus", "contig", "start", "end", "TopEffectGeneName", "TopEffectFeatureId"])
    #
    vaf_DF_dp_patho_CRC_HIGH_merge_Control.to_csv("vaf_DF_dp_patho_CRC_HIGH_merge_Control.tsv",
                                                  na_rep="NaN", sep="\t", header=True, index=True)
    vaf_DF_dp_patho_CRC_HIGH_merge_Control_Polyps.to_csv("vaf_DF_dp_patho_CRC_HIGH_merge_Control_Polyps.tsv",
                                                         na_rep="NaN", sep="\t", header=True, index=True)



def main():
    # Concatenate VCF files
    cat_DF = cat_vcf_files(indir=parser.inDir,
                           info_cols={
                                      "END":"int",
                                      "HIAF":"float",
                                      "HICNT":"int",
                                      "HICOV":"int",
                                      "MQ":"float",
                                      "MSI":"float",
                                      "MSILEN":"float",
                                      "NM":"float",
                                      "PMEAN":"float",
                                      "QUAL":"float",
                                      "TYPE":"str",
                                      "SAMPLE":"str",
                                      "CSQ":"str"
                                     })

    #
    # Merge MUT and sample manifest tables, and write filtered MUT-derived
    # tables of VAFs as matrices
    merge_mut_mani(mut_DF=cat_DF,
                   mani_csv=parser.inDir + "/" + parser.sampleManifest,
                   cvdir=os.path.join(parser.inDir, parser.cvDir))


if __name__ == "__main__":
    main()
