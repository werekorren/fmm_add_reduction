import matplotlib.pyplot as plt
import seaborn as sns # pip include seaborn
import numpy as np

raw_data = { # 8378 entries
  28: 2,
  29: 107,
  30: 362,
  31: 378,
  32: 791,
  33: 1251,
  34: 1855,
  35: 1088,
  36: 1134,
  37: 399,
  38: 351,
  39: 263,
  40: 165,
  41: 98,
  42: 66,
  43: 40,
  44: 16,
  45: 9,
  46: 2,
  47: 0,
  48: 0,
  49: 1
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
    plt.savefig('fmm_addition_distribution_225_20260809_140650.png')

    # display plot
    #plt.show()

if __name__ == "__main__":
    main()
