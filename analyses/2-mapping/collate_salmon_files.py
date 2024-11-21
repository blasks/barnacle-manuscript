#!/usr/bin/env python
# AUTHOR: Stephen Blaskowski
# CREATE DATE: 18 February 2021


import argparse
import glob
import json
import numpy as np 
import os
import pandas as pd


def handle_arguments():
    '''
    returns argument parser
    '''
    description = '''
        This program collates quant.sf files that are the output from 
        Salmon read mapping. The TPM and raw read values from multiple
        mappings are collected into a single gzipped csv file.

        Example usage: ./collate_salmon_files.py input_dir output_dir
        '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input_dir', type=str, help='Directory containing multiple Salmon output directories')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('-c', '--column', default='Both', type=str, help='"Both", "TPM" or "NumReads"')
    parser.add_argument('-n', '--norm', type=str, help='Filepath with normalization factors. Must contain "SampleID", and "NormFactor" columns corresponding to samples represented in input directory')
    return parser


def build_norm_map(norm_factors_filepath):
    norm_df = pd.read_csv(norm_factors_filepath)
    norm_dict = dict(zip(norm_df['SampleID'], norm_df['NormFactor']))
    return norm_dict


def import_quant_files(path_list, sample_list, quant_col, norm_factors=None):
    if quant_col == 'Both':
        columns = ['Name', 'Length', 'TPM', 'NumReads']
    else:
        columns = ['Name', 'Length', quant_col]
    df = pd.DataFrame
    for i, filepath in enumerate(path_list):
        print('{} of {}: Reading data from {}'.format(i, len(path_list), filepath))
        new_df = pd.read_csv(filepath, sep='\t', usecols=columns)
        # add in SampleId column
        new_df['SampleID'] = sample_list[i]
        # normalize read counts (if applicable)
        if norm_factors is not None:
            new_df['NormedReads'] = new_df['NumReads'] * norm_factors[sample_list[i]]
        if i == 0:
            df = new_df
        else:
            df = pd.concat([df, new_df], ignore_index=True)
    return df


def main():
    # parse arguments
    parser = handle_arguments()
    args = parser.parse_args()
    # check output path
    output_path = '{}/collated_salmon_data.csv.gz'.format(args.output_dir)
    if os.path.isfile(output_path):
        raise FileExistsError('A file by the name of {} already exists.'.format(output_path))
    print('Collating {} values from quant.sf files found under {}\n'.format(args.column, args.input_dir))
    # define quant filepaths
    path_quant_files = '{}/*/quant.sf'.format(args.input_dir)
    # collect sample names and quant files
    path_list_quant_files = []
    sample_list = []
    for filepath in glob.glob(path_quant_files):
        path_list_quant_files.append(filepath)
        sample_list.append(os.path.basename(os.path.dirname(filepath)))
    # generate norm factor list (if necessary)
    if args.norm is not None:
        norm_dict = build_norm_map(args.norm)
        # iterate through sample list and parse each dataset
        quant_df = import_quant_files(path_list_quant_files, sample_list, args.column, norm_dict)
    else:
        quant_df = import_quant_files(path_list_quant_files, sample_list, args.column)
    # extract genome name and gene id from mapping name
    quant_df['GenomeName'], quant_df['GeneID'] = list(zip(*[name.rsplit("_", 1) for name in  quant_df['Name']]))
    # rename Salmon data columns
    quant_df.rename(columns={'Name':'MappingName', 'Length':'GeneLength'}, inplace=True)
    # save dataframe to output directory
    quant_df.to_csv(output_path, index=False, compression='gzip')
    print('Values successfully collated.')

    # parse metadata info
    output_path = '{}/salmon_metadata.csv'.format(args.output_dir)
    # check output path
    if os.path.isfile(output_path):
        raise FileExistsError('A file by the name of {} already exists.'.format(output_path))
    print('Collating metadata from meta_info.json files found under {}\n'.format(args.input_dir))
    path_metadata = '{}/*/aux_info/meta_info.json'.format(args.input_dir)
    json_list = []
    for filepath in glob.glob(path_metadata):
        # parse json
        with open(filepath) as file:
            json_file = json.load(file)
        # append SampleID
        json_file['SampleID'] = os.path.basename(os.path.dirname(os.path.dirname(filepath)))
        # compile list
        json_list.append(json_file)
    # make DataFrame
    metadata_df = pd.DataFrame(json_list)
    # save dataframe to output directory
    metadata_df.to_csv(output_path, index=False)
    print('Metadata successfully collated.')


if __name__ == "__main__":
    main()

