import os
import sys
import argparse
import pysam
import yaml
import json
import numpy as np
import pandas as pd
from collections import defaultdict, namedtuple

# from Ardal import *


def dataframe_to_dict(df):
    ## convert dataframe to a dictionary
    result_dict = {}

    for allele_id, allele_data in df.items():
        allele_dict = allele_data.to_dict()
        result_dict[allele_id] = allele_dict

    return result_dict


def merge_dicts(dicts):
    merged_dict = {}

    for current_dict in dicts:
        for allele_id, allele_data in current_dict.items():
            if allele_id not in merged_dict:
                merged_dict[allele_id] = {}

            for sample_id, value in allele_data.items():
                if sample_id not in merged_dict[allele_id]:
                    merged_dict[allele_id][sample_id] = 0
                merged_dict[allele_id][sample_id] += value

    return merged_dict


def dict_to_dataframe(input_dict):
    df = pd.DataFrame.from_dict(input_dict, orient='index')
    df = df.transpose()
    df = df.fillna(0)
    return df


def mergeAMs(allele_matrix, database_matrix):
    
    dataframes = [allele_matrix, database_matrix]
    
    dicts = []
    ## Iterate through each DataFrame
    for df in dataframes:
        d = dataframe_to_dict(df)
        dicts.append(d)

    merged_dicts = merge_dicts(dicts)

    # print(dict_to_dataframe(merged_dicts))
    return dict_to_dataframe(merged_dicts)


def read_json(jfile):
    with open(jfile, 'r') as fin:
        return json.load(fin)
    

def countSNPs(ard_mat, sample_ids):

    counts_dict = defaultdict(dict)
    allale_mapping = ard_mat.toDict()

    for id in sample_ids:
        unique_alleles = ard_mat.unique([id])
        total_alleles = allale_mapping[id]
        counts_dict[id]["total_snps"] = len(total_alleles)
        counts_dict[id]["unique_snps"] = list(unique_alleles)

    with open("allele_stats.json", 'w') as fout:
        json.dump(counts_dict, fout, indent=4)


def isolate_and_save(allele_matrix, output_prefix):
    """Isolates binary matrix, converts to uint8, and saves array and headers."""

    try:
        # Convert directly to uint8 NumPy array
        matrix_array = allele_matrix.values.astype(np.uint8)

        # Save NumPy array
        np.save(f"{output_prefix}_matrix.npy", matrix_array)

        # Store headers
        headers = {
            "guids": allele_matrix.index.tolist(),
            "alleles": allele_matrix.columns.tolist()
        }

        with open(f"{output_prefix}_headers.json", "w") as f:
            json.dump(headers, f, indent=4)

        print(f"Matrix saved as '{output_prefix}_matrix.npy', headers as '{output_prefix}_headers.json'")

    except Exception as e:
        print(f"An error occurred: {e}")


def runmerged(args):
    allele_dict = defaultdict(dict)

    vcf_path = args.merged_vcf
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)


    for record in vcf.fetch():

        for sample_id in samples:
            try:
                call = record.samples[sample_id]
            except KeyError:
                continue

            if call['GT'] == (0, 0):
                allele_dict[f"{record.chrom}.{record.pos}.{record.ref}"][sample_id] = 0 

            ## Handles multi allelic sites
            for allele_index in call['GT']:
                if allele_index != 0 and allele_index != None and len(record.alleles[allele_index]) == 1 and record.qual >= args.qual:
                    allele_id = f"{record.chrom}.{record.pos}.{record.alleles[allele_index]}"

                    if allele_id not in allele_dict:
                        for s in samples:
                            allele_dict[allele_id][s] = 0 

                    allele_dict[allele_id][sample_id] = 1

    allele_matrix = pd.DataFrame.from_dict(allele_dict).fillna(0)

    if args.database:
        try:
            database_matrix = pd.read_csv(args.database, index_col=0)
            allele_matrix = mergeAMs(allele_matrix, database_matrix)
        except Exception as e:
            print(f"Failed to merge database {args.database}: {e}")

    isolate_and_save(allele_matrix, args.prefix)

    allele_matrix.to_csv(f"{args.prefix}.csv")
    # allele_matrix.to_parquet(f"{args.prefix}.parquet", compression=None)
    # df = pd.read_csv("./core_allele_matrix.csv")
    # df.to_parquet("./core_allele_matrix.parquet", compression=None)


def runVCF(args):

    allele_dict = defaultdict(dict)

    samples = []

    for vcf_path in args.vcfs:
        sample_id = vcf_path.split("/")[-1].split(".")[0]

        if os.path.exists(f"{sample_id}.bed"):
            continue

        samples.append(sample_id)
        vcf = pysam.VariantFile(vcf_path)

        print(f"Stating {sample_id}...")

        # Create a list to store BED file rows for this VCF
        bed_rows = []

        for record in vcf.fetch():
            try:
                call = record.samples[sample_id]
            except KeyError:
                continue

            # if call['GT'] == (0, 0):
            #     allele_dict[f"{record.chrom}.{record.pos}.{record.ref}"][sample_id] = 0

            # Handles multi allelic sites
            for allele_index in call['GT']:
                if allele_index != 0 and allele_index is not None and len(record.alleles[allele_index]) == 1 and record.qual >= args.qual:
                    alt_nucleotide = record.alleles[allele_index]
                    allele_id = f"{record.chrom}.{record.pos}.{alt_nucleotide}"

                    # if allele_id not in allele_dict:
                    #     for s in samples:
                    #         allele_dict[allele_id][s] = 0

                    # allele_dict[allele_id][sample_id] = 1

                    # Add a row to the BED data for this SNP
                    bed_row = [
                        record.chrom,
                        record.pos - 1,  # BED format is 0-based, VCF is 1-based
                        record.pos,
                        allele_id,
                    ]
                    bed_rows.append(bed_row)

        # Write the BED file for this sample
        bed_df = pd.DataFrame(bed_rows, columns=["chrom", "start", "end", "description"])
        bed_df.to_csv(f"{sample_id}.bed", sep="\t", index=False, header=False)

    # allele_matrix = pd.DataFrame.from_dict(allele_dict).fillna(0)

    # if args.database:
    #     try:
    #         database_matrix = pd.read_csv(args.database, index_col=0)
    #         allele_matrix = mergeAMs(allele_matrix, database_matrix)
    #     except Exception as e:
    #         print(f"Failed to merge database {args.database}: {e}")

    # isolate_and_save(allele_matrix, args.prefix)

    # allele_matrix.to_csv(f"{args.prefix}.csv")
    # allele_matrix.to_parquet(f"{args.prefix}.parquet", compression=None)
    # df = pd.read_csv("./core_allele_matrix.csv")
    # df.to_parquet("./core_allele_matrix.parquet", compression=None)

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('--merged_vcf', default=None, action='store', help='Path to a merged VCF file.')
    parser.add_argument('--vcfs', default=None, nargs="*", action='store', help='Path to a directory containing VCF files.')
    parser.add_argument('--qual', action='store', default=30, type=int, help='Quality score threshold.')
    parser.add_argument('--database', action='store', default=False, help='Allele matrix database.')
    parser.add_argument('--prefix', action='store', default="allele_matrix", help='Prefix for output file.')
    args = parser.parse_args(argv)

    if args.merged_vcf != None and args.vcfs != None:
        print("Can't use --vcf and --merged_vcf in the same run.")
        sys.exit(1)
    elif args.merged_vcf == None and args.vcfs == None:
        print("--vcf or --merged_vcf required.")
        sys.exit(1)

    return args


if __name__=="__main__":
    args = parse_args(sys.argv)

    if args.merged_vcf != None:
        runmerged(args)
    elif args.vcfs != None:
        runVCF(args)
