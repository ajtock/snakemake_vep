#!/usr/bin/env python3

# Author: Andy Tock
# Date: 03/01/2023

# Usage:
# conda activate 6run_analysis
# ./combine_reshape_vcf.py -d /Users/ajtock/dnanexus/snakemake_vep/workflow/vep_data/output -m /Users/ajtock/dnanexus/snakemake_vep/workflow/resources/manifest/ori02_ts_sample_manifest_2022_09_26_EW.csv
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

    # Convert CSQ into a list in which each element is one annotation containing multiple subfields
    vcf["CSQ_ANN_LIST"] = vcf["CSQ"].str.split(",")
    # Convert each CSQ_ANN_LIST element into a list of subfields
    vcf["CSQ_LIST"] = vcf["CSQ_ANN_LIST"] \
        .apply(lambda x: [y.split("|") for y in x])
    # Extract given subfields from CSQ_LIST as separate columns
    vcf["consequence"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[1] for y in x]))).keys())))
    vcf["impact"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[2] for y in x]))).keys())))
    vcf["gene_symbol"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[3] for y in x]))).keys())))
    vcf["gene_ID"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[4] for y in x]))).keys())))
    vcf["feature_type"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[5] for y in x]))).keys())))
    vcf["feature_ID"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[6] for y in x]))).keys())))
    vcf["biotype"] = vcf["CSQ_LIST"] \
        .apply(lambda x: "__".join(list(dict.fromkeys(list(filter(None, [y[7] for y in x]))).keys())))

    # Convert GT into list and add columns containing list elements
    vcf["GT"] = vcf["GT"].str.split(":")
    vcf["genotype"] = [x[0] for x in vcf["GT"]]
    vcf["genotype2"] = ["hom" if x == "1/1" else "het" if x in ["0/1", "1/0"] else "NaN" for x in vcf["genotype"]]
    vcf["ref_depth"] = [int(x[1].split(",")[0]) for x in vcf["GT"]]
    vcf["alt_depth"] = [int(x[1].split(",")[1]) for x in vcf["GT"]]
    vcf["variant_depth"] = [int(x[5]) for x in vcf["GT"]]
    vcf["total_depth"] = [int(x[3]) for x in vcf["GT"]]
    vcf["total_depth_less_no_calls"] = vcf["alt_depth"] + vcf["ref_depth"]
    vcf["no_calls"] = vcf["total_depth"] - vcf["total_depth_less_no_calls"]
    vcf["vaf"] = vcf["alt_depth"] / (vcf["alt_depth"] + vcf["ref_depth"])

    # Remove unnecessary columns
    vcf.drop(columns=["CSQ_ANN_LIST","CSQ_LIST"], inplace=True)

    # Write to TSV
    vcf.to_csv("results/concat_sample_VCFs.tsv",
               na_rep="NaN", sep="\t", header=True, index=False)

    return vcf


def merge_vcf_mani(vcf_DF, mani_csv):
    """
    Extend concatenated VCF DataFrame with additional information from
    sample manifest CSV, including "Subject ID", "run_no", "cat" (categorical response variable;
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
                        merge_DF["END"].astype(str) + "_" + \
                        merge_DF["gene_symbol"].astype(str) + "_" + \
                        merge_DF["ALT"].astype(str)
    merge_DF_dp = merge_DF.loc[
        (merge_DF["Mean On-Target Duplex Depth"] >= 1000)
    ]
    merge_DF_dp_CRC_Control = merge_DF_dp.loc[
        (merge_DF_dp["cat"].isin(["CRC+", "Control"]))
    ]

    #merge_DF_dp = merge_DF.loc[
    #    (merge_DF["Mean On-Target Duplex Depth"] >= 1000) & \
    #    (merge_DF["alt_depth"] >= 3) & \
    #    (merge_DF["gene_symbol"].str.contains("|".join(gene_list)))
    ##    (merge_DF["vaf"] < 0.25) & \
    ##    (~merge_DF["FILTER"].str.contains("|".join(["NM8.0","Q10","q22.5"])))
    #]

    # Get "pathogenic" variants based on VEP-derived
    # ClinVar clinical significance annotations
    merge_DF_dp_patho = merge_DF_dp.loc[
        (merge_DF_dp["alt_depth"] >= 3) & \
        (merge_DF_dp["gene_symbol"].str.contains("|".join(gene_list))) & \
        #(merge_DF_dp["vaf"] < 0.25) & \
        (~merge_DF_dp["FILTER"].str.contains("|".join(["NM8.0", "Q10", "q22.5"]))) & \
        (merge_DF_dp["CSQ"].str.contains("pathogenic", case=False)) & \
        (~merge_DF_dp["CSQ"].str.contains("Conflicting", case=False))
    ]
    merge_DF_dp_CRC_Control_patho = merge_DF_dp_CRC_Control.loc[
        (merge_DF_dp_CRC_Control["alt_depth"] >= 3) & \
        (merge_DF_dp_CRC_Control["gene_symbol"].str.contains("|".join(gene_list))) & \
        # (merge_DF_dp_CRC_Control["vaf"] < 0.25) & \
        (~merge_DF_dp_CRC_Control["FILTER"].str.contains("|".join(["NM8.0", "Q10", "q22.5"]))) & \
        (merge_DF_dp_CRC_Control["CSQ"].str.contains("pathogenic", case=False)) & \
        (~merge_DF_dp_CRC_Control["CSQ"].str.contains("Conflicting", case=False))
        ]

    # Reshape as matrices of VAF values
    vaf_DF_dp = merge_DF_dp.pivot_table(index=["locus"],
                                        columns=["Sample ID", "Subject ID", "cat", "run_no"],
                                        values="vaf",
                                        fill_value="NaN")
    vaf_DF_dp_patho = merge_DF_dp_patho.pivot_table(index=["locus"],
                                                    columns=["Sample ID", "Subject ID", "cat", "run_no"],
                                                    values="vaf",
                                                    fill_value="NaN")
    vaf_DF_dp_CRC_Control = merge_DF_dp_CRC_Control.pivot_table(index=["locus"],
                                                                columns=["Sample ID", "Subject ID",
                                                                         "cat", "run_no"],
                                                                values="vaf",
                                                                fill_value="NaN")
    vaf_DF_dp_CRC_Control_patho = merge_DF_dp_CRC_Control_patho.pivot_table(index=["locus"],
                                                                            columns=["Sample ID", "Subject ID",
                                                                                     "cat", "run_no"],
                                                                            values="vaf",
                                                                            fill_value="NaN")

    # Merge VAF matrices and remove duplicate columns
    vaf_DF_dp_diff_cols = vaf_DF_dp.columns.difference(vaf_DF_dp_patho.columns)
    vaf_DF_dp_patho_merge_vaf_DF_dp = pd.merge(left=vaf_DF_dp_patho,
                                               right=vaf_DF_dp[vaf_DF_dp_diff_cols],
                                               how="left",
                                               on=["locus"])
    vaf_DF_dp_CRC_Control_diff_cols = vaf_DF_dp_CRC_Control.columns.difference(vaf_DF_dp_CRC_Control_patho.columns)
    vaf_DF_dp_CRC_Control_patho_merge_vaf_DF_dp_CRC_Control = pd.merge(left=vaf_DF_dp_CRC_Control_patho,
                                                                       right=vaf_DF_dp_CRC_Control[vaf_DF_dp_CRC_Control_diff_cols],
                                                                       how="left",
                                                                       on=["locus"])
    vaf_DF_dp_patho.to_csv("results/vaf_DF_dp_patho.tsv",
                           na_rep="NaN", sep="\t", header=True, index=True)
    vaf_DF_dp_patho_merge_vaf_DF_dp.to_csv("results/vaf_DF_dp_patho_merge_vaf_DF_dp.tsv",
                                           na_rep="NaN", sep="\t", header=True, index=True)
    vaf_DF_dp_CRC_Control_patho.to_csv("results/vaf_DF_dp_CRC_Control_patho.tsv",
                                       na_rep="NaN", sep="\t", header=True, index=True)
    vaf_DF_dp_CRC_Control_patho_merge_vaf_DF_dp_CRC_Control.to_csv("results/vaf_DF_dp_CRC_Control_patho_merge_vaf_DF_dp_CRC_Control.tsv",
                                                                   na_rep="NaN", sep="\t", header=True, index=True)


def main():
    # Concatenate VCF files and convert into DataFrame with
    # additional columns corresponding to subfields of
    # INFO, CSQ (VEP-generated), and GT fields
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


    # Merge MUT and sample manifest tables, and write filtered MUT-derived
    # tables of VAFs as matrices
    merge_vcf_mani(vcf_DF=cat_DF,
                   mani_csv=parser.sampleManifest)


if __name__ == "__main__":
    main()
