#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 2021

@author: Palash Sashittal
"""

import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from collections import Counter
import datetime as dt
from venn import venn, generate_petal_labels, draw_venn, generate_colors

def is_mutation_present(seq, target_seq):
    if seq.find(target_seq) != -1:
        return True
    else:
        return False

def probe_fasta(fasta_file, primer_name_list, primer_dict):

    seq_dict = {}
    entry = False
    with open(fasta_file,'r') as inp:
        for line in inp:
            if line.startswith('>'):
                if entry:
                    mutation_attendance = [is_mutation_present(seq, primer_dict[primer_name]) for primer_name in primer_name_list]
                    seq_dict[seq_name] = mutation_attendance
                    
                seq_name = line.lstrip('>').rstrip('\n')
                seq = ''
                entry = True
            else:
                seq += line.rstrip('\n')
                    
    mutation_attendance = [is_mutation_present(seq, primer_dict[primer_name]) for primer_name in primer_name_list]
    seq_dict[seq_name] = mutation_attendance

    return seq_dict

def get_days_from_start(str_date, start_date = dt.date(2020, 12, 1)):
    if len(str.split(str_date, '-')) == 3:
        curr_year, curr_month, curr_day = map(int, str.split(str_date, '-'))
        curr_date = dt.date(curr_year, curr_month, curr_day)
        return (curr_date - start_date).days
    else:
        return np.nan

def get_condensed_lineage(lineage, resolution = 3):
    if lineage != lineage:
        return lineage
    split_lineage = lineage.split('.')
    if len(split_lineage) < resolution:
        return lineage
    else:
        return '.'.join(split_lineage[:resolution])

def get_probe_dataframe(seq_dict, target_name_list):
    target_data = []
    for seq_name, val in seq_dict.items():
        virus_name = seq_name.split('|')[0]
        target_data.append([virus_name] + val)

    return pd.DataFrame(target_data, columns = ['Virus name'] + target_name_list).set_index(['Virus name'])

def plot_temporal_variation(df_target, df_probe_meta, figure_prefix, startdate):
    for probe_lineage in df_target['lineage'].dropna().unique():
        probe_target_name_list = list(df_target[df_target['lineage'] == probe_lineage]['id'])

        resolution = len(probe_lineage.split('.'))
        df_probe_meta['condensed lineage'] = df_probe_meta['Pango lineage'].apply(lambda x: get_condensed_lineage(x, resolution))
        for target_name in probe_target_name_list:
            selected_data = df_probe_meta[df_probe_meta[target_name]]
            selected_data[probe_lineage] = selected_data['condensed lineage'].apply(lambda x: 'True positive' if x == probe_lineage else 'False positive')

            unselected_data = df_probe_meta[~df_probe_meta[target_name]]
            unselected_data[probe_lineage] = unselected_data['condensed lineage'].apply(lambda x: x == probe_lineage)
            unselected_data = unselected_data[unselected_data[probe_lineage]]
            unselected_data[probe_lineage] = 'False negative'

            selected_data = pd.concat([selected_data, unselected_data])

            fig, ax = plt.subplots(1,1,figsize=(8,4))
            sns.histplot(data=selected_data, x='time', hue=probe_lineage, multiple='stack', ax=ax, bins = 50, palette={'True positive':sns.color_palette()[0], 'False positive':sns.color_palette()[1], 'False negative': sns.color_palette()[2]})
            plt.xticks(fontsize=20);
            plt.ylabel(f'# of {target_name} samples', fontsize = 20)
            plt.xlabel(f'days since {startdate}', fontsize = 20)
            plt.yticks(fontsize=20);

            plt.savefig(f'{figure_prefix}temporal_variation_{probe_lineage}_{target_name}.pdf', bbox_inches='tight')

def plot_specificity_sensitivity(df_target, df_target_analysis, df_probe_meta, figure_prefix, startdate):
    time_bins = np.arange(0, max(df_probe_meta['time']), 20)
    
    for probe_lineage in df_target['lineage'].dropna().unique():
        probe_target_name_list = list(df_target[df_target['lineage'] == probe_lineage]['id'])
        
        resolution = len(probe_lineage.split('.'))
        df_probe_meta['condensed lineage'] = df_probe_meta['Pango lineage'].apply(lambda x: get_condensed_lineage(x, resolution))
        for target_name in probe_target_name_list:
            selected_data = df_probe_meta[df_probe_meta[target_name]]
            selected_data[probe_lineage] = selected_data['condensed lineage'].apply(lambda x: 'True positive' if x == probe_lineage else 'False positive')

            nlineage_counts, _ = np.histogram(df_probe_meta[df_probe_meta['condensed lineage'] == probe_lineage]['time'], bins = time_bins)
            true_positive_counts, _ = np.histogram(selected_data[selected_data[probe_lineage] == 'True positive']['time'], bins = time_bins)
            false_positive_counts, _ = np.histogram(selected_data[selected_data[probe_lineage] == 'False positive']['time'], bins = time_bins)
            negative_counts, _ = np.histogram(df_probe_meta[df_probe_meta['condensed lineage'] != probe_lineage]['time'], bins = time_bins)

            sensitivity = true_positive_counts / nlineage_counts if nlineage_counts !=0 else np.nan
            specificity = (negative_counts - false_positive_counts) / negative_counts if negative_counts !=0 else np.nan


            fig, ax = plt.subplots(1,1,figsize=(8,4))
            df_plot = pd.DataFrame({'time':list((time_bins[:-1] + time_bins[1:])/2)*2, 'measure': list(sensitivity) + list(specificity), 'kind': ['sensitivity']*len(sensitivity) + ['specificity']*len(specificity)})
            sns.lineplot(data=df_plot, x='time', y='measure', hue='kind')
            plt.xticks(fontsize=20);
            plt.ylabel(f'sensitivity/specificity', fontsize = 20)
            plt.xlabel(f'days since {startdate}', fontsize = 20)
            plt.yticks(fontsize=20);
            plt.ylim([0,1])

            ax2 = plt.twinx()
            df_total = pd.DataFrame({'time':list((time_bins[:-1] + time_bins[1:])/2), 'measure': list(nlineage_counts)})
            ax2 = sns.lineplot(data=df_total, ax=ax2, x='time', y='measure', color='black')
            ax2.set_ylabel(f'# of samples of {probe_lineage}', fontsize = 20)
            plt.yticks(fontsize=20);

            plt.savefig(f'{figure_prefix}sensitivity_specificity_{probe_lineage}_{target_name}.pdf', bbox_inches = 'tight')


def analyze_probe_lineage(df_target_meta, probe_lineage, target_name_list, probe_target_name_list):
    analyze_probe_lineage_data = []
    ntargets = len(probe_target_name_list)

    resolution = len(probe_lineage.split('.'))

    df_target_meta['condensed lineage'] = df_target_meta['Pango lineage'].apply(lambda x: get_condensed_lineage(x, resolution))
    nlineage = len(df_target_meta[df_target_meta['condensed lineage'] == probe_lineage])
    nnegative = len(df_target_meta) - nlineage

    case_id = 0
    for case in itertools.product([True, False], repeat=ntargets):
        df_select = df_target_meta.copy()
        for target_id in range(ntargets):
            df_select = df_select[df_select[probe_target_name_list[target_id]] == case[target_id]]

        tp = len(df_select[df_select['condensed lineage'] == probe_lineage])
        fp = len(df_select[df_select['condensed lineage'] != probe_lineage])
        fn = nlineage - tp
        tn = nnegative - fp

        sensitivity = tp / nlineage if nlineage != 0 else np.nan
        specificity = tn / nnegative if nnegative !=0 else np.nan

        case_entry = [probe_lineage]
        for target in target_name_list:
            if target in probe_target_name_list:
                case_entry.append(case[probe_target_name_list.index(target)])
            else:
                case_entry.append(np.nan)

        analyze_probe_lineage_data.append(case_entry + [case_id, nlineage, nnegative, tp, fp, fn, tn, sensitivity, specificity])

        case_id += 1

    df_target_meta = df_target_meta.drop(columns = ['condensed lineage'])

    return analyze_probe_lineage_data

def merge_probe_metadata(df_probe_data, df_meta):
    df_probe_data = df_probe_data.loc[~df_probe_data.index.duplicated(keep='first')]
    df_meta = df_meta.loc[~df_meta.index.duplicated(keep='first')]

    return pd.concat([df_meta, df_probe_data], axis=1)

def get_analysis_results(df_target, df_probe_meta, target_name_list):
    analysis_data = []

    for probe_lineage in df_target['lineage'].dropna().unique():
        probe_target_name_list = list(df_target[df_target['lineage'] == probe_lineage]['id'])
        analysis_data += analyze_probe_lineage(df_probe_meta, probe_lineage, target_name_list, probe_target_name_list)

    return pd.DataFrame(analysis_data, columns = ['lineage'] + target_name_list + ['case_id', 'P', 'N', 'TP', 'FP', 'FN', 'TN', 'sensitivity', 'specificity'])

def get_petal_labels(probe_lineage, df_target_analysis, probe_target_name_list):
    petal_labels = {}

    df_labels = df_target_analysis[df_target_analysis['lineage'] == probe_lineage][probe_target_name_list + ['TP', 'FP']]
    assert(len(df_labels) == 2**len(probe_target_name_list))
    for _, row in df_labels.iterrows():
        case_key = ''
        for target_name in probe_target_name_list:
            if row[target_name] == True:
                case_key += '1'
            else:
                case_key += '0'

        petal_labels[case_key + '0'] = row['FP']
        petal_labels[case_key + '1'] = row['TP']

    del petal_labels['0000']

    return petal_labels

def plot_venn_diagrams(df_target, df_target_analysis, figure_prefix):
    for probe_lineage in df_target['lineage'].dropna().unique():
        probe_target_name_list = list(df_target[df_target['lineage'] == probe_lineage]['id'])
        petal_labels = get_petal_labels(probe_lineage, df_target_analysis, probe_target_name_list)

        fig, ax = plt.subplots(1, 1, figsize=(8,8))
        draw_venn(petal_labels = petal_labels, dataset_labels = probe_target_name_list + [probe_lineage], colors = generate_colors(n_colors = len(probe_target_name_list) + 1),
                  hint_hidden = False, fontsize = 14, legend_loc = 'best', ax = ax, figsize = (10, 10))

        plt.savefig(f'{figure_prefix}_venn_{probe_lineage}.pdf', bbox_inches = 'tight')


def string_to_date(date_string):
    curr_year, curr_month, curr_day = map(int, str.split(date_string, '-'))
    return dt.date(curr_year, curr_month, curr_day)

def impose_time_constraints(df_probe_meta, startdate, enddate):

    if enddate:
        end_date = string_to_date(enddate)
        df_probe_meta['time'] = df_probe_meta['Collection date'].apply(lambda x: get_days_from_start(x, start_date = end_date))
        df_probe_meta = df_probe_meta[df_probe_meta['time'] <= 0]

    start_date = string_to_date(startdate)
    df_probe_meta['time'] = df_probe_meta['Collection date'].apply(lambda x: get_days_from_start(x, start_date = start_date))
    df_probe_meta = df_probe_meta[df_probe_meta['time'] >= 0]

    return df_probe_meta

def main(args):

    df_target = pd.read_csv(args.targets)
    target_name_list = list(df_target['id'])
    if args.probe:
        target_dict = dict(zip(df_target['id'], df_target['sequence']))

        seq_dict = probe_fasta(args.fasta, target_name_list, target_dict)
        df_probe_data = get_probe_dataframe(seq_dict, target_name_list)

        if args.outprobe:
            df_probe_data.to_csv(args.outprobe)

        if args.outprobemeta:
            if not args.meta:
                print('cannot output probemeta without meta as input')
            else:
                df_meta = pd.read_csv(args.meta, sep='\t', usecols=['Virus name', 'Collection date', 'Pango lineage']).set_index('Virus name')
                df_probe_meta = merge_probe_metadata(df_probe_data, df_meta)
                df_probe_meta.to_csv(args.outprobemeta)


    if args.analyze:
        if not args.probe:
            df_probe_data = pd.read_csv(args.inprobe).set_index('Virus name')

        df_meta = pd.read_csv(args.meta, sep='\t', usecols=['Virus name', 'Collection date', 'Pango lineage']).set_index('Virus name')

        df_probe_meta = merge_probe_metadata(df_probe_data, df_meta)

        df_probe_meta = impose_time_constraints(df_probe_meta, args.startdate, args.enddate)
        #df_probe_meta = pd.concat([df_meta, df_probe_data], axis=1)
        #df_target_meta.to_csv('temp_file.tsv', sep='\t')

        if args.outprobemeta:
            df_probe_meta.to_csv(args.outprobemeta)

        df_target_analysis = get_analysis_results(df_target, df_probe_meta, target_name_list)
        df_target_analysis.to_csv(args.outanalysis, index=False)

        if args.outhtml:
            html_buf = df_target_analysis.to_html()
            html_file = open(args.outhtml, 'w')
            html_file.write(html_buf)
            html_file.close()

    if args.plot:
        if not args.analyze:
            df_target_analysis = pd.read_csv(args.inanalysis)
            df_probe_meta = pd.read_csv(args.inprobemeta)
            df_probe_meta = impose_time_constraints(df_probe_meta, args.startdate, args.enddate)

        # venn diagram
        plot_venn_diagrams(df_target, df_target_analysis, args.prefix)

        # temporal variation plot
        #plot_temporal_variation(df_target, df_probe_meta, args.prefix, args.startdate)
        #plot_specificity_sensitivity(df_target, df_target_analysis, df_probe_meta, args.prefix, args.startdate)


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--probe', help='probe target sequences', action='store_true', default=False)
    parser.add_argument('--meta', type=str, help='metadata of the SARS-CoV-2 genomes')
    parser.add_argument('--fasta', type=str, help='fasta file with SARS-CoV-2 genomes')
    parser.add_argument('--targets', type=str, help='csv file with target sequences and their probe lineages')
    parser.add_argument('--outprobe', type=str, help='output file name for updated metadata file')
    parser.add_argument('--condmeta', type=str, help='output file name for condensed form of metadata file')

    parser.add_argument('--analyze', help='analyze target sequences', action='store_true', default=False)
    parser.add_argument('--inprobe', type=str, help='input metadata file of probed target sequences')
    parser.add_argument('--outanalysis', type=str, help='output csv file with analysis results')
    parser.add_argument('--outprobemeta', type=str, help='output csv file with probing results along with metadata')
    parser.add_argument('--outhtml', type=str, help='output html file with anlaysis results')
    parser.add_argument('--country', type=str, help='focus on a specific country for analysis')

    parser.add_argument('--startdate', type=str, help='date (YYYY-MM-DD) after which the analysis should be performed [2020-12-01]', default='2020-12-01')
    parser.add_argument('--enddate', type=str, help='date (YYYY-MM-DD) until which the analysis should be performed [last entry in the metadata file]')
    parser.add_argument('--plot', help='plot the analysis results', action='store_true', default=False)
    parser.add_argument('--inanalysis', type=str, help='csv file with analysis results')
    parser.add_argument('--inprobemeta', type=str, help='csv file with probing results along with metadata (output of analyze mode)')
    parser.add_argument('--prefix', type=str, help='prefix for output files of plotting')

    # specificity and sensitivity plots for every combination
    # temporal variation plot for every combination of target sequences

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    if args.probe or args.analyze or args.plot:
        main(args)
    else:
        print('INVALID INPUT')
        print('must perform one of the following')
        print('1. probe')
        print('2. analyze')
        print('3. plot')
