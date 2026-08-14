import matplotlib.pyplot as plt
import seaborn as sns # pip include seaborn
import numpy as np

raw_data = { # 10227 entries
  15: 779,
  16: 5710,
  17: 0,
  18: 3738
}

def main():
    data = []
    for num_additions, count in raw_data.items():
        data.extend([num_additions] * count)

    # create histogram
    plt.figure(figsize=(9, 6)) # Optional: Adjust figure size
    #ax = sns.histplot(data, bins=sorted(raw_data.keys()), kde=True, color='red', edgecolor='black', facecolor='skyblue', discrete=True)
    ax = sns.histplot(data, bins=sorted(raw_data.keys()), kde=False, color='red', edgecolor='black', facecolor='skyblue', discrete=True)

    # add labels and title
    ax.set_xlabel("Num Additions After Reduction")
    ax.set_ylabel("Num FMM Algorithms")
    #ax.set_title("Addition Count Distribution")
    ax.set_xticks(sorted(raw_data.keys())) #  x-axis ticks to match addition count

    # save plot to png file
    plt.savefig('fmm_addition_distribution_222_20260808_192215.png')

    # display plot
    #plt.show()

if __name__ == "__main__":
    main()
