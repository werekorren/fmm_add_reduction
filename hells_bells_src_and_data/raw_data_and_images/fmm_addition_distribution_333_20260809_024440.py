import matplotlib.pyplot as plt
import seaborn as sns # pip include seaborn
import numpy as np

raw_data = { # 15316 entries
  59: 1,
  60: 18,
  61: 149,
  62: 654,
  63: 1021,
  64: 1156,
  65: 2135,
  66: 2068,
  67: 2072,
  68: 2081,
  69: 1254,
  70: 923,
  71: 845,
  72: 418,
  73: 258,
  74: 132,
  75: 66,
  76: 36,
  77: 16,
  78: 9,
  79: 1,
  80: 2,
  81: 1
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
    plt.savefig('fmm_addition_distribution_333_20260809_024440.png')

    # display plot
    #plt.show()

if __name__ == "__main__":
    main()
