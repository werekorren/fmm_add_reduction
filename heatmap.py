import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter
from matplotlib.ticker import MaxNLocator

# Use output file name if provided command-line
if len(sys.argv) > 1:
    heatmap_file_name = sys.argv[1]
else:
    heatmap_file_name = "heatmap.png"

# Also use input file name if provided command-line
if len(sys.argv) > 2:
    input_data_file_name = sys.argv[2]
else:
    input_data_file_name = "additions_before_and_after_reduction.txt"

# Load x y pairs from text file
data = np.loadtxt(input_data_file_name)

# Validate inputs
if data.ndim != 2 or data.shape[1] != 2:
    raise ValueError("Input file must contain exactly two columns: x y.")
x = data[:, 0]
y = data[:, 1]

# Validate x and y
if x.ndim != 1 or y.ndim != 1:
    raise ValueError("Both x and y must be 1D vectors.")
if len(x) != len(y):
    raise ValueError("x and y must have equal length.")
if not np.all(np.equal(np.mod(x, 1), 0)) or not np.all(np.equal(np.mod(y, 1), 0)):
    raise ValueError("Assumes x and y contain integer values.")

x = x.astype(int)
y = y.astype(int)

# Integer coordinate ranges
x_min, x_max = x.min(), x.max()
y_min, y_max = y.min(), y.max()

# Exact count matrix: counts[i, j] = number of times (x=i+x_min, y=j+y_min) occurs
counts = np.zeros((x_max - x_min + 1, y_max - y_min + 1), dtype=float)
for xi, yi in zip(x, y):
    counts[xi - x_min, yi - y_min] += 1

# Optional smoothing for a more solid / less pixelated appearance
sigma = 0.5
counts_smooth = gaussian_filter(counts, sigma=sigma)
# no smoothing
#counts_smooth = counts
#interpolation = "nearest"

#
vmax=counts.max()

# white → light blue → blue → light green → green → yellow → orange → red → purple → black
colors = [
    (1.0, 1.0, 1.0),   # white (low)
    (0.1, 0.1, 0.4),   # blue
    (0.7, 0.7, 1.0),   # light blue
    (0.7, 1.0, 0.7),   # light green
    (0.1, 0.4, 0.1),   # green
    (1.0, 1.0, 0.0),   # yellow
    (1.0, 0.6, 0.0),   # orange
    (1.0, 0.0, 0.0),   # red
    (0.5, 0.0, 0.5),   # purple
    (0.0, 0.0, 0.0),   # black (high)
    (0.0, 0.0, 0.0)    # black (high)
]
cmap = LinearSegmentedColormap.from_list("custom_heatmap", colors, N=256)

# Plot
fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(
    counts_smooth.T,
    origin="lower",
    extent=[x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5],
    aspect="auto",
    cmap=cmap,
    interpolation="bicubic",
    vmin=0,
    vmax=counts.max(),   # colorbar top corresponds to true max count
)

# Axes limits
#ax.set_xlim(x_min - 0.5, x_max + 0.5)
#ax.set_ylim(y_min - 0.5, y_max + 0.5)
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# Nice ticks including min and max
x_ticks = MaxNLocator(nbins=6, integer=True).tick_values(x_min, x_max)
x_ticks = np.unique(np.concatenate(([x_min], x_ticks, [x_max])))
x_ticks = x_ticks[(x_ticks >= x_min) & (x_ticks <= x_max)]
ax.set_xticks(x_ticks)

y_ticks = MaxNLocator(nbins=6, integer=True).tick_values(y_min, y_max)
y_ticks = np.unique(np.concatenate(([y_min], y_ticks, [y_max])))
y_ticks = y_ticks[(y_ticks >= y_min) & (y_ticks <= y_max)]
ax.set_yticks(y_ticks)

ax.set_xlabel("Num Additions Before Reduction")
ax.set_ylabel("Num Additions After Reduction")
#ax.set_title("Addition Reduction Heatmap")

# Colorbar labeled in counts
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Num FMM Schemes")
cbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()
plt.savefig(heatmap_file_name, dpi=1200)
#plt.show()