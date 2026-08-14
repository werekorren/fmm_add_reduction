import matplotlib.pyplot as plt
import seaborn as sns # pip include seaborn
import numpy as np

raw_data = { # 4464 entries
  26: 4,
  27: 98,
  28: 216,
  29: 267,
  30: 523,
  31: 817,
  32: 1088,
  33: 589,
  34: 613,
  35: 98,
  36: 95,
  37: 19,
  38: 25,
  39: 9,
  40: 3
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
    plt.savefig('fmm_addition_distribution_224_20260809_130652.png')

    # display plot
    #plt.show()

if __name__ == "__main__":
    main()
