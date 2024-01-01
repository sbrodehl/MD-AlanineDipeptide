import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Function to load histograms and calculate mean and std
def load_histograms_and_calculate_stats(directory):
    histograms = []
    temps = []
    concs = []

    for filename in os.listdir(directory):
        if filename.endswith('.npy'):
            # Correctly parse the filename to extract temp and conc
            # Assuming filename format: "histogram_tempXXX_concYYY.npy"
            parts = filename.split('_')
            temp = int(parts[1][4:])
            conc = int(parts[2][5:7])
            data = np.load(os.path.join(directory, filename), allow_pickle=True).item()
            histograms.append(data['histogram'])
            temps.append(temp)
            concs.append(conc)

    histograms = np.array(histograms)
    temps = np.array(temps)
    concs = np.array(concs)

    mean_histogram = np.mean(histograms, axis=0)
    std_histogram = np.std(histograms, axis=0)

    return mean_histogram, std_histogram, temps, concs, histograms

# Function to plot histogram (mean or std)
def plot_histogram(hist, title, colorbar=True):
    plt.figure()
    plt.imshow(hist.T, extent=[-180, 180, -180, 180], origin='lower', norm=LogNorm())
    if colorbar:
        plt.colorbar()
    plt.xlabel('$\phi$')
    plt.ylabel('$\psi$')
    plt.title(title)
    plt.show()

# Function to plot 2D scatter plot for a specific bin
# Function to find bins with top 10 biggest standard deviations
def find_top_std_bins(std_hist, top_n=10):
    # Flatten the std_hist array and get the indices of the top N values
    flat_indices = np.argsort(std_hist.ravel())[-top_n:]
    # Convert flat indices to 2D indices
    return np.unravel_index(flat_indices, std_hist.shape)

# Modified function to plot 2D scatter plot for specific bins
def plot_2d_bin_function(histograms, temps, concs, std_hist, top_n=30):
    top_bins = find_top_std_bins(std_hist, top_n)
    for bin_x, bin_y in zip(*top_bins):
        bin_values = histograms[:, bin_x, bin_y]

        plt.figure()
        plt.scatter(concs, temps, c=bin_values, cmap='viridis')
        plt.colorbar(label='Bin Value')
        plt.xlabel('Concentration')
        plt.ylabel('Temperature')
        plt.title(f'Bin Value at ({bin_x}, {bin_y}) for Each Temp and Conc')
        plt.show()


def plot_3d_bin_function(histograms, temps, concs, std_hist, top_n=30):
    top_bins = find_top_std_bins(std_hist, top_n)
    for bin_x, bin_y in zip(*top_bins):
        bin_values = histograms[:, bin_x, bin_y]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(concs, temps, bin_values, c=bin_values, cmap='viridis', depthshade=False)
        ax.set_xlabel('Concentration')
        ax.set_ylabel('Temperature')
        ax.set_zlabel('Bin Value')
        ax.set_title(f'3D View - Bin Value at ({bin_x}, {bin_y}) for Each Temp and Conc')
        plt.show()


def plot_max_std_bins_in_blocks(histograms, temps, concs, std_hist, block_size=45):
    # Calculate the number of blocks
    num_blocks = 360 // block_size

    for i in range(num_blocks):
        for j in range(num_blocks):
            # Define the block boundaries
            x_start, x_end = i * block_size, (i + 1) * block_size
            y_start, y_end = j * block_size, (j + 1) * block_size

            # Extract the block from the std_hist
            block = std_hist[x_start:x_end, y_start:y_end]

            # Find the position of the maximum std within this block
            max_pos = np.unravel_index(np.argmax(block), block.shape)

            # Convert local block position to global position
            global_max_pos = (max_pos[0] + x_start, max_pos[1] + y_start)

            # Extract bin values for this position
            bin_values = histograms[:, global_max_pos[0], global_max_pos[1]]

            # 2D scatter plot for the bin
            plt.figure()
            plt.scatter(concs, temps, c=bin_values, cmap='viridis')
            plt.colorbar(label='Bin Value')
            plt.xlabel('Concentration')
            plt.ylabel('Temperature')
            plt.title(f'Bin Value at ({global_max_pos[0]}, {global_max_pos[1]}) for Each Temp and Conc')
            plt.show()



directory = 'hist'
mean_hist, std_hist, temps, concs, histograms = load_histograms_and_calculate_stats(directory)


plot_histogram(mean_hist, 'Mean Histogram')
plot_histogram(std_hist, 'Standard Deviation Histogram')

plot_max_std_bins_in_blocks(histograms, temps, concs, std_hist)


