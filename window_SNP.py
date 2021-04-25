import pandas as pd
import numpy as np
import sys, getopt


dataset_1 = "/Volumes/BossDaddy/GWASAnalysis/pruned_uk_bb_SNPs.txt"
dataset_2 = "/Volumes/BossDaddy/GWASAnalysis/Meta_Analysis/plink_african_formatted_cutoff.meta"
final_file_name = "/Volumes/BossDaddy/GWASAnalysis/WindowOutput/uk_AFR_windowed_10kb.txt"


def check_snps(uk_bbj_sig_snps):
    significant_snp_df = pd.read_csv(uk_bbj_sig_snps, sep='\t', index_col=False, skiprows=1, names=["CHR", "BP", "SNP",
                                                                                                    "BETA", "SE", "P",
                                                                                                    "A1", "A2"])
    significant_snp_df.CHR = significant_snp_df.CHR.replace('X', 24)
    significant_snp_df.CHR = significant_snp_df.CHR.astype(int)
    significant_snp_df = significant_snp_df.dropna()
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
    count = 0
    length = len(sig_snp_df.BP)
    print("Number of D1 SNPs to be compared on next chromosome " + str(length))
    for index, row in sig_snp_df.iterrows():
        sample2_df_window_subset = sample2_df[np.logical_and(np.abs(sample2_df.BP - row.BP) <= window_size,
                                                             sample2_df.P <= 5 * (10 ** -6))].copy()
        sample2_df_window_subset['BP_dataset1'] = row.BP
        row_list.append(sample2_df_window_subset)
    subset_df = pd.concat(row_list)
    merged_df = pd.merge(subset_df, sig_snp_df,
                         left_on=subset_df['BP_dataset1'],
                         right_on=sig_snp_df['BP'],
                         suffixes=['_d2', '_d1'])
    filtered_df = merged_df.loc[merged_df.groupby('BP_dataset1').P_d2.idxmin()]
    return filtered_df


def main(argv):
    # both input files
    dataset_1 = ''
    dataset_2 = ''
    # file for results to be written to
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
        print("Finished merge for chromosome " + str(chrm))
    # changes order to move bp_uk to 2nd column and renames column names
    output_df = output_df.drop(columns=['BP_dataset1', 'CHR_d1'])
    # renaming columns for consistancy
    output_df = output_df.rename(columns={'key_0': 'BP_dataset1', 'BP_d2': 'BP_dataset2', 'CHR_d2': 'CHR',
                                          'P(R)': 'P(R)_d2', 'SNP': 'SNP_d1'})
    # sorts data and removes duplicate rows
    output_df = output_df.sort_values(by=['CHR', 'BP_dataset2']).drop_duplicates()
    output_df['BP_distance'] = output_df['BP_dataset1'] - output_df['BP_dataset2']
    headers = ['CHR', 'SNP_d1', 'BP_dataset1', 'BP_dataset2', 'BP_distance', 'A1_d1', 'A2_d1', 'A1_d2', 'A2_d2', 'N',
               'P_d1', 'P_d2', 'P(R)_d2', 'BETA_d1', 'BETA_d2']
    output_df = output_df[headers]
    f = open(outputfile, 'w')
    f.write("dataset_1 " + dataset_1.split('/')[-1] + '\n')
    f.write("dataset_2 " + dataset_2.split('/')[-1] + '\n')
    f.write("window_size " + str(window_size) + '\n')
    f.write(output_df.to_csv(index=False, sep='\t'))
    f.close()


if __name__ == '__main__':
    main(sys.argv[1:])
