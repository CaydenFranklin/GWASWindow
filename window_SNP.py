import pandas as pd
import numpy as np
import sys, getopt

# TODO: Currently hardcoded, will add as arguments
#dataset_1 = "/Volumes/BossDaddy/GWASAnalysis/FINAL_FIXED_uk_sig_into_clumped_bbj_replication.txt"
#dataset_2 = "/Volumes/BossDaddy/GWASAnalysis/Meta_Analysis/plink_african_formatted_cutoff.meta"
#final_file_name = "/Volumes/BossDaddy/GWASAnalysis/Meta_Analysis/uk_BBJ_AFR_windowed_50kb_replication_ordered_unique.txt"


def check_snps(uk_bbj_sig_snps):
    significant_snp_df = pd.read_csv(uk_bbj_sig_snps, sep='\t', skiprows=1, names=["CHR", "BP", "OR", "SE", "P.x", "A1",
                                                                                   "A2", "SNP", "P.y"])
    return significant_snp_df


def chunk_afr(chromosomes, data_iterator, chunks):
    chr_to_data_dict = {}
    # Each chunk is in dataframe format
    i = 1
    for chunk in data_iterator:
        data_chunk = chunk.copy(deep=True)
        for chromosome in chromosomes:
            filtered_chunk = data_chunk.loc[data_chunk['CHR'] == chromosome]
            # prevents addition of empty dataframes
            if filtered_chunk.shape[0] != 0:
                # creates empty list if this is the first entry from this chromosome read in
                if chromosome not in chr_to_data_dict.keys():
                    chr_to_data_dict[chromosome] = []
                chr_to_data_dict[chromosome].append(filtered_chunk)
        print("Chunk " + str(i) + " read in of " + str(round(chunks, 0)))
        i += 1
    # concatenates all dataframes belonging to one chromosome into one
    for chrm in chr_to_data_dict:
        chr_to_data_dict[chrm] = pd.concat(chr_to_data_dict[chrm])
    return chr_to_data_dict


def find_sig_snps(window_size, sig_snp_df, sample2_df):
    row_list = []
    for bp in sig_snp_df.BP:
        sample2_df_subset = sample2_df[np.abs(sample2_df.BP - bp) <= window_size].copy()
        sample2_df_subset['BP_dataset1'] = bp
        sample2_df_subset = pd.merge(sample2_df_subset, sig_snp_df, left_on=sample2_df_subset['BP_dataset1'],
                                     right_on=sig_snp_df['BP'],
                                     suffixes=['_d2', '_d1'])
        row_list.append(sample2_df_subset)
    subset_df = pd.concat(row_list)
    filtered_df = subset_df.loc[subset_df.groupby('BP_dataset1').P.idxmin()]
    return filtered_df


def main(argv):
    # both input files
    dataset_1 = ''
    dataset_2 = ''
    #file for results to be written to
    outputfile = ''
    # maximum base pair distance between significant SNPs on each dataset
    window_size = 0
    try:
        # A and B used instead of 1 and 2 because shortargs only allows single char args
        opts, args = getopt.getopt(argv, "ha:b:o:w:")
    except getopt.GetoptError:
        print('window_SNP.py -a <dataset_a> -b <dataset_b> -o <outputfile> -w <window size in BP>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('window_SNP.py -a <dataset_a> -b <dataset_b> -o <outputfile> -w <window size in BP>')
            sys.exit()
        elif opt == "-a":
            dataset_1 = arg
        elif opt == "-b":
            dataset_2 = arg
        elif opt == "-o":
            outputfile = arg
        elif opt == "-w":
            window_size = int(arg)
    print('Dataset A: ', dataset_1)
    print('Dataset B: ', dataset_2)
    print('Output File : ', outputfile)
    print('Window Size : ', window_size)

    chunk_size = 100000

    sig_snp_df = check_snps(dataset_1)
    chromosomes = sig_snp_df['CHR'].unique()
    headers = ['CHR', 'BP', 'SNP', 'A1', 'A2', 'N', 'P', 'P(R)', 'BETA', 'BETA(R)', 'Q', 'I']
    data_iterator = pd.read_csv(dataset_2, sep=' ', chunksize=chunk_size, names=headers, skiprows=1)
    num_lines = sum(1 for row in open(dataset_2, 'r'))
    chr_sig_replicated_snps = {}
    afr_chr_to_data = chunk_afr(chromosomes, data_iterator, num_lines / chunk_size)
    output_df = pd.DataFrame()
    for chrm in chromosomes:
        sample2_chr_set = afr_chr_to_data[chrm]
        sample1_chr_set = sig_snp_df[sig_snp_df['CHR'] == chrm]
        chr_sig_replicated_snps[chrm] = find_sig_snps(window_size, sample1_chr_set, sample2_chr_set)
        output_df = pd.concat([output_df, chr_sig_replicated_snps[chrm]])
    # changes order to move bp_uk to 2nd column and renames column names
    output_df = output_df.drop(columns=['BP_dataset1', 'SNP_d2', 'CHR_d1'])
    #renaming columns for consistancy
    output_df = output_df.rename(columns={'key_0': 'BP_dataset1', 'BP_d2': 'BP_dataset2', 'CHR_d2': 'CHR',
                                          'BETA': 'BETA_d2', 'BETA(R)': 'BETA(R)_d2', 'OR': 'OR_d1',
                                          'P.x': 'P.x_d1', 'P.y': 'P.y_d1', 'P': 'P_d2', 'P(R)': 'P(R)_d2'})
    # sorts data and removes duplicate rows
    output_df = output_df.sort_values(by=['CHR', 'BP_dataset2']).drop_duplicates()
    output_df['BP_distance'] = np.abs(output_df['BP_dataset1'] - output_df['BP_dataset2'])
    # Beta is ln(OR)
    output_df['BETA_d1'] = np.log(output_df['OR_d1'])
    headers = ['CHR', 'BP_dataset1', 'BP_dataset2', 'BP_distance', 'A1_d1', 'A2_d1', 'A1_d2', 'A2_d2', 'N', 'P_d2', 'P(R)_d2',
               'P.x_d1', 'P.y_d1', 'BETA_d1', 'BETA_d2']
    output_df = output_df[headers]
    f = open(outputfile, 'w')
    f.write("dataset_1 " + dataset_1.split('/')[-1] + '\n')
    f.write("dataset_2 " + dataset_2.split('/')[-1] + '\n')
    f.write("window_size " + str(window_size) + '\n')
    f.write(output_df.to_string(index=False))
    f.close()

if __name__ == '__main__':
    main(sys.argv[1:])