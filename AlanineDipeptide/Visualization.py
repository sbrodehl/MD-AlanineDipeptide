import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree
from tqdm import tqdm


def plot_frame_transitions(df, frames, title, lag_time=2, phi_offset=60, psi_offset=-90):
    phi_s, phi_e = -180 + phi_offset, 180 + phi_offset
    psi_s, psi_e = -180 + psi_offset, 180 + psi_offset

    # "improved" x/y axis ticks
    plt.xticks(np.arange(-180, 180, 30), np.arange(phi_s, phi_e, 30))
    plt.yticks(np.arange(-180, 180, 30), np.arange(psi_s, psi_e, 30))

    # axis labels (right order?)
    plt.xlabel('$\phi$')
    plt.ylabel('$\psi$')

    # 1:1 aspect ratio
    plt.axes().set_aspect('equal')
    # remove grid lines
    plt.axes().grid(False)

    if len(title):
        plt.title(title)

    sd = df.as_matrix()
    plt.scatter(sd[0], sd[1], color='red')

    for key in tqdm(frames.keys(), desc="Plotting transition paths"):
        prev = df.iloc[key]
        c = np.random.rand(3,)
        for step in range(lag_time):
            curr = df.iloc[key + step]
            plt.plot([prev[0], curr[0]], [prev[1], curr[1]], color=c, linewidth=1.0)
            prev = curr

    return plt


def plot_groups(df, c_centers, phi_offset=60, psi_offset=-90):
    phi_s, phi_e = -180 + phi_offset, 180 + phi_offset
    psi_s, psi_e = -180 + psi_offset, 180 + psi_offset

    tree = cKDTree(c_centers)
    points = np.c_[df['phi'], df['psi']]  # fancy zipping
    queries = tree.query(points)[1]
    plt.figure()
    plt.scatter(df['phi'], df['psi'], c=queries, s=5, linewidths=0.25)

    # "improved" x/y axis ticks
    plt.xlim([-180, 180])
    plt.ylim([-180, 180])
    # "improved" x/y axis ticks
    plt.xticks(np.arange(-180, 180, 30), np.arange(phi_s, phi_e, 30))
    plt.yticks(np.arange(-180, 180, 30), np.arange(psi_s, psi_e, 30))

    # axis labels (right order?)
    plt.ylabel('$\psi$')
    plt.xlabel('$\phi$')

    # 1:1 aspect ratio
    plt.axes().set_aspect('equal')
    # remove grid lines
    plt.axes().grid(False)

    return plt


def create_ramachandran_plot(df, title, colorbar=True, phi_offset=60, psi_offset=-90):
    phi_s, phi_e = -180 + phi_offset, 180 + phi_offset
    psi_s, psi_e = -180 + psi_offset, 180 + psi_offset

    plt.hist2d(df[0], df[1],
               range=[[-180, 180], [-180, 180]],  # not really necessary
               bins=360,
               # cmap='viridis',  # cf. http://matplotlib.org/examples/color/colormaps_reference.html
               norm=LogNorm())
    if colorbar:
        plt.colorbar()

    # "improved" x/y axis ticks
    plt.xticks(np.arange(-180, 180, 30), np.arange(phi_s, phi_e, 30))
    plt.yticks(np.arange(-180, 180, 30), np.arange(psi_s, psi_e, 30))

    # axis labels (right order?)
    plt.xlabel('$\phi$')
    plt.ylabel('$\psi$')

    # 1:1 aspect ratio
    plt.axes().set_aspect('equal')
    # remove grid lines
    plt.axes().grid(False)

    if len(title):
        plt.title(title)

    return plt


def plot_voronoi_ridges(ridges, linewidth=1):
    for segment in ridges:
        plt.plot([segment[0][0], segment[1][0]], [segment[0][1], segment[1][1]], linewidth=linewidth, color='black')


def plot_sampled_frames(frames, linewidth=1, n=-1):
    sample_frames = pd.DataFrame(frames)
    sd = sample_frames.as_matrix()
    plt.scatter(sd[0][:n], sd[1][:n], color='red', linewidth=linewidth)
