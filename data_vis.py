import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics

dataset_1 = "/Volumes/BossDaddy/GWASAnalysis/WindowOutput/f_uk_AFR_windowed_0kb.txt"
dataset_2 = "/Volumes/BossDaddy/GWASAnalysis/WindowOutput/f_uk_AFR_windowed_10kb.txt"
dataset_3 = "/Volumes/BossDaddy/GWASAnalysis/WindowOutput/f_uk_AFR_windowed_50kb.txt"
dataset_4 = "/Volumes/BossDaddy/GWASAnalysis/WindowOutput/f_uk_AFR_windowed_100kb.txt"
dataset_5 = "/Volumes/BossDaddy/GWASAnalysis/WindowOutput/f_uk_AFR_windowed_200kb.txt"



def plot_scatter(d, ax, i, name):
    df = pd.read_csv(d, sep='\t', index_col=False, skiprows=4, names=["CHR", "SNP_d1", "BP_dataset1",
                                                                      "BP_dataset2",
                                                                      "BP_distance", "A1_d1", "A2_d1",
                                                                      "A1_d2", "A2_d2", "N", "P_d1", "P_d2",
                                                                      "P(R)_d2", "BETA_d1", "BETA_d2"])
    df['BP_distance'] = df['BP_distance'].astype('float')
    # To plot outliers, comment the below out
    # if i == 3:
    #     df = df.loc[np.logical_and(np.abs(df['BETA_d1']) <= .05, np.abs(df['BETA_d2']) < .5)]
    # elif i >= 4:
    #     df = df.loc[np.logical_and(np.abs(df['BETA_d1']) <= .08, np.abs(df['BETA_d2']) < .5)]
    scatter = ax.scatter(np.abs(df['BETA_d1']), np.abs(df['BETA_d2']),
                         c=np.abs(df['BP_distance']))
    m, b = np.polyfit(np.abs(df['BETA_d1']), np.abs(df['BETA_d2']), 1)
    df['predicted_b2'] = np.abs(df['BETA_d1']) * m + b

    r2 = round(sklearn.metrics.r2_score(np.abs(df['BETA_d2']), df['predicted_b2']), 2)
    df['outlier_stat'] = np.abs(df['BETA_d2']) - np.abs(df['predicted_b2'])

    ax.plot(np.abs(df['BETA_d1']), df['predicted_b2'], color='red')
    ax.set_xlabel("abs(BETA) UK Biobank")
    ax.set_ylabel("abs(BETA) AFR Meta")
    ax.set_title('BETA vs BETA (' + str(name) + ')')
    if i == 0:
        for index, row in df.iterrows():
            k = 1
            if abs(row['BETA_d1']) > .055:
                k = -7
            ax.annotate(row['SNP_d1'], (abs(row['BETA_d1']), abs(row['BETA_d2'])), textcoords="offset points",
                        xytext=(k * 10, 0))
        ax.legend(*scatter.legend_elements(num=7), loc="center right", title="BP Distance")
    else:
        df = df.sort_values(by=['outlier_stat'], ascending=False)
        df.reset_index()
        top_x = np.ceil(.05 * df.shape[0])
        n = 0
        for index, row in df.head(int(top_x)).iterrows():
            n += 1
            j = (n % 2 * 2) - 1
            ax.annotate(row['SNP_d1'], (abs(row['BETA_d1']), abs(row['BETA_d2'])), textcoords="offset points",
                        xytext=(5, j * 5 * np.log10(row['outlier_stat'])))
        ax.legend(*scatter.legend_elements(num=7), loc="upper right", title="BP Distance")
    ax.annotate('$R^2$' + ': ' + str(r2), xy=(0.8, 0.3), xycoords='axes fraction')

def plot_hist(d, ax, name):
    df = pd.read_csv(d, sep='\t', index_col=False, skiprows=4, names=["CHR", "SNP_d1", "BP_dataset1",
                                                                      "BP_dataset2",
                                                                      "BP_distance", "A1_d1", "A2_d1",
                                                                      "A1_d2", "A2_d2", "N", "P_d1", "P_d2",
                                                                      "P(R)_d2", "BETA_d1", "BETA_d2"])
    ax.hist(np.abs(df['BP_distance']), density=1, edgecolor='black')
    ax.set_xlabel("Distance (BP)")
    ax.set_ylabel("Occurrences")
    ax.set_title('Distribution of BP Distances' + ' (' + str(name) + ')')

if __name__ == '__main__':
    d_arr = [dataset_1, dataset_2, dataset_3, dataset_4, dataset_5]
    n_arr = ["Exact Match", "10kb", "50kb", "100kb", "200kb"]
    fig1, f1_axes = plt.subplots(2, 3, figsize=(15, 10))
    f1_axes[1][2].set_visible(False)
    f1_axes[1][0].set_position([0.24, 0.125, 0.228, 0.343])
    f1_axes[1][1].set_position([0.55, 0.125, 0.228, 0.343])

    for i in range(len(d_arr)):
        d = d_arr[i]
        name = n_arr[i]
        ax = f1_axes.flatten()[i]
        plot_scatter(d, ax, i, name)
    fig1.savefig('/Volumes/BossDaddy/GWASAnalysis/WindowOutput/' + 'BETA vs BETA (All + Outliers)')
    fig1.clf()

    fig2, f2_axes = plt.subplots(2, 2, figsize=(15, 10))
    for z in range(len(d_arr)-1):
        d = d_arr[z+1]
        name = n_arr[z+1]
        ax = f2_axes.flatten()[z]
        plot_hist(d, ax, name)
    fig2.suptitle("Distribution of BP Distances")
    fig2.savefig('/Volumes/BossDaddy/GWASAnalysis/WindowOutput/' + 'Histogram (All)')
    fig2.clf()


