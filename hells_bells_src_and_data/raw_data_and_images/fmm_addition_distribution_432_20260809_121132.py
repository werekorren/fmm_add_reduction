import matplotlib.pyplot as plt
import seaborn as sns # pip include seaborn
import numpy as np

raw_data = { # 342 entries
  54: 3,
  55: 3,
  56: 28,
  57: 12,
  58: 64,
  59: 32,
  60: 25,
  61: 48,
  62: 36,
  63: 44,
  64: 17,
  65: 11,
  66: 10,
  67: 5,
  68: 2,
  69: 0,
  70: 0,
  71: 1,
  72: 0,
  73: 0,
  74: 1
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
    plt.savefig('fmm_addition_distribution_432_20260809_121132.png')

    # display plot
    #plt.show()

if __name__ == "__main__":
    main()
