import matplotlib.pyplot as plt
import seaborn as sns # pip include seaborn
import numpy as np

raw_data = { # 1631 entries
  40: 1,
  41: 5,
  42: 15,
  43: 9,
  44: 41,
  45: 85,
  46: 100,
  47: 259,
  48: 205,
  49: 290,
  50: 268,
  51: 176,
  52: 79,
  53: 24,
  54: 29,
  55: 24,
  56: 12,
  57: 5,
  58: 4
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
    plt.savefig('fmm_addition_distribution_227_20260809_152550.png')

    # display plot
    #plt.show()

if __name__ == "__main__":
    main()
