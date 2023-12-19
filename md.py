import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import glob
import random
import MDAnalysis

from scipy.spatial import Voronoi

from AlanineDipeptide.Computations import *
from AlanineDipeptide.Utilities import *
from AlanineDipeptide.Visualization import *


if __name__ == '__main__':

    DEBUG = True
    samples = 20

    # first, run a long simulation
    #run_macro_simulation(debug=DEBUG)

    df = pd.DataFrame(read_xvg('data/md_long_nojump_rama.xvg'))
    df = apply_offset(df)

    if df.shape[0] < 1000:
        print("Not enough steps saved for analysis.")
        print("Please configure .mdp for GROMACS simulation.")
        exit(1)

    # draw point density
    create_ramachandran_plot(df, "")
    print("Saving figure 'Ramachandran-Plot'")
    plt.savefig('plots/Ramachandran-Plot', dpi=300)
    plt.clf()

    ms_centers = clustering(df)
    vor = Voronoi(ms_centers)

    # save cluster centers to disk
    with open('data/mean_shift_clusters', 'wb') as outfile:
        pickle.dump(ms_centers, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    # draw point density
    create_ramachandran_plot(df, "")
    # draw cluster centers
    plt.scatter(ms_centers.T[0], ms_centers.T[1], linewidths=6, color='black')
    print("Saving figure 'Ramachandran-Meanshift-Plot'")
    plt.savefig('plots/Ramachandran-Meanshift-Plot', dpi=300)
    plt.clf()

    create_ramachandran_plot(df, "")
    # draw cluster centers
    plt.scatter(ms_centers.T[0], ms_centers.T[1], linewidths=6, color='black')

    ridge_lines = create_voronoi_ridges(vor)

    # draw voronoi border
    plot_voronoi_ridges(ridge_lines)
    print("Saving figure 'Ramachandran-Voronoi-Plot'")
    plt.savefig('plots/Ramachandran-Voronoi-Plot', dpi=300)
    plt.clf()

    # load sampled points from disk or sample new points
    sample_files = glob.glob("data/" + str(samples) + "_transition_states_*")
    if len(sample_files):
        with open(sample_files[0], 'rb') as outfile:
            sampled_frames = pickle.load(outfile)
    else:
        # choose only interesting ridges to sample from
        # sampled_frames = sample_boundaries(df, ridge_lines[3:6:1], dist=15, samples=samples)
        # sample from all borders
        sampled_frames = sample_boundaries(df, ridge_lines, dist=15, samples=samples)

    # draw point density
    create_ramachandran_plot(df, "")
    # draw voronoi border
    plot_voronoi_ridges(ridge_lines)
    # draw samples
    plot_sampled_frames(sampled_frames)

    print("Saving figure 'Ramachandran-Transitions-Plot'")
    plt.savefig('plots/Ramachandran-Transitions-Plot', dpi=300)
    plt.clf()

    plot_voronoi_ridges(ridge_lines)
    plot_sampled_frames(sampled_frames, linewidth=1)
    plot_frame_transitions(df, sampled_frames, "")

    print("Saving figure 'Ramachandran-Transitions-Paths-Plot'")
    plt.savefig('plots/Ramachandran-Transitions-Paths-Plot', dpi=300)
    plt.clf()

    paths_files = glob.glob("data/" + str(samples) + "_micro_sims_*")
    if len(paths_files):
        with open(paths_files[0], 'rb') as outfile:
            f_all = pickle.load(outfile)
    else:
        # run small simulations
        frame_keys = sampled_frames.keys()
        # Load simulation results
        u = MDAnalysis.Universe('data/md_long.gro', 'data/md_long_corr.xtc')

        # save frame ids as pdb files
        save_frames_as_pdb(sampled_frames, u)
        # run small simulations and save convert to numpy
        f_all = run_micro_simulations(frame_keys, n=20, debug=DEBUG)

    # plot interesting voronoi region
    plot_voronoi_ridges(ridge_lines)

    # plot all micro transition paths
    traj = samples
    keys = list(f_all.keys())
    random.shuffle(keys)
    for idx in tqdm(keys[:traj], desc="Plotting micro transition paths"):
        c = np.random.rand(3,)
        for t in f_all[idx]:
            npt = pd.DataFrame(np.array(t))
            npt = apply_offset(npt)
            # plot big transition
            prev = df.iloc[idx]
            curr = df.iloc[idx + 1]
            plt.plot([prev[0], curr[0]], [prev[1], curr[1]], color='black', linewidth=2)
            plt.plot([prev[0], curr[0]], [prev[1], curr[1]], color=c, linewidth=1)
            # plot micro transition
            plt.scatter(npt[0][0], npt[1][0], color='black', linewidth=2)
            plt.scatter(npt[0][0], npt[1][0], color=c, linewidth=1)
            plt.plot(npt[0], npt[1], color=c, linewidth=1)
            # plot difference
            plt.plot([npt[0][npt.shape[0]-1], curr[0]], [npt[1][npt.shape[0]-1], curr[1]], linestyle='--', color=c)

    # compute transition matrices and save them
    t_mats = {}
    for i in tqdm(range(1, 4), desc="Computing transition matrices"):
        trans_mat, time_spent = transition_matrix(df, ms_centers, lag_time=i)
        t_mats[i] = {'transition': trans_mat, 'time': time_spent}

    # save transition matrices to disk
    with open('data/transition_matrices', 'wb') as outfile:
        pickle.dump(t_mats, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    print("Saving figure 'Micro-Transitions-Paths-Plot'")
    plt.savefig('plots/Micro-Transitions-Paths-Plot', dpi=300)
    plt.clf()
