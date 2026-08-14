import matplotlib.pyplot as plt
import seaborn as sns # pip include seaborn
import numpy as np

raw_data = { # 1499 entries
  51: 1,
  52: 5,
  53: 4,
  54: 23,
  55: 18,
  56: 68,
  57: 98,
  58: 135,
  59: 180,
  60: 218,
  61: 285,
  62: 230,
  63: 176,
  64: 52,
  65: 2,
  66: 2,
  67: 1,
  68: 1
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
    plt.savefig('fmm_addition_distribution_228_20260809_153226.png')

    # display plot
    #plt.show()

if __name__ == "__main__":
    main()
