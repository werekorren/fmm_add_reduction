import matplotlib.pyplot as plt
import seaborn as sns # pip include seaborn
import numpy as np

raw_data = { # 256 entries
  55: 5,
  56: 14,
  57: 10,
  58: 45,
  59: 29,
  60: 23,
  61: 30,
  62: 32,
  63: 43,
  64: 11,
  65: 5,
  66: 8,
  67: 1
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
    plt.savefig('fmm_addition_distribution_234_20260809_122929.png')

    # display plot
    #plt.show()

if __name__ == "__main__":
    main()
