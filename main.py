import pandas as pd
import numpy as np

# TODO: Currently hardcoded, will add as arguments
dataset_1 = "/Volumes/BossDaddy/GWASAnalysis/ALL_FIXED_uk_sig_into_clumped_bbj_replication.txt"
dataset_2 = "/Volumes/BossDaddy/GWASAnalysis/Meta_Analysis/plink_african_formatted_cutoff.meta"
final_file_name = "/Volumes/BossDaddy/GWASAnalysis/Meta_Analysis/uk_BBJ_AFR_windowed_replication_ordered_unique.txt"


def check_snps(uk_bbj_sig_snps):
    significant_snp_df = pd.read_csv(uk_bbj_sig_snps, sep='\t', skiprows=1, names=["CHR", "BP", "P.x", "SNP", "P.y"])
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
        row_list.append(sample2_df_subset)
    subset_df = pd.concat(row_list)
    filtered_df = subset_df.loc[subset_df.groupby('BP_dataset1').P.idxmin()]
    return filtered_df


if __name__ == '__main__':
    chunk_size = 100000
    # maximum base pair distance between significant SNPs on each dataset
    window_size = 50000
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
    headers.insert(1, 'BP_dataset1')
    output_df = output_df[headers]
    output_df = output_df.rename(columns={"BP": "BP_dataset2"})

    # sorts data and removes duplicate rows
    output_df = output_df.sort_values(by=['CHR', 'BP_dataset2']).drop_duplicates()
    f = open(final_file_name, 'w')
    f.write("dataset_1 " + dataset_1.split('/')[-1] + '\n')
    f.write("dataset_2 " + dataset_2.split('/')[-1] + '\n')
    f.write(output_df.to_string(index=False))
    f.close()
