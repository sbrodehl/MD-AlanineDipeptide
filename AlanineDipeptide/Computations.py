import shlex
import subprocess
import datetime
import os
import pickle
import time
from sympy.geometry import *

import numpy as np
from scipy.spatial import cKDTree
from sklearn.cluster import MeanShift, estimate_bandwidth
from tqdm import tqdm

from AlanineDipeptide.Utilities import read_xvg


def clustering(df):
    # MeanShift clustering
    # The following bandwidth can be automatically detected using
    bandwidth = estimate_bandwidth(df.as_matrix(), quantile=0.2, n_samples=500)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True).fit(df.as_matrix())

    return ms.cluster_centers_


def sample_boundaries(df, ridges, dist=15, samples=10):
    # now sample frames and check if those points are near the voronoi ridges
    delta_dist = dist  # in degree
    frames = {}
    times = []
    with tqdm(total=samples, desc="sampling frames") as pbar:
        # sample transition 'positions'
        while len(frames) < samples:
            start = time.time()  # measure the time spend searching
            # choose random frame number
            rnd = np.random.randint(0, df.shape[0] - 1)
            frame = df.iloc[rnd]
            sample_point = Point(frame[0], frame[1])
            # check if the angles are in the defined range of the ridges
            for seg in ridges:
                rnd_offset = np.random.randint(-0.2 * delta_dist, 0.2 * delta_dist)
                ridge = Segment(seg[0], seg[1])
                if ridge.distance(sample_point) <= delta_dist + rnd_offset:
                    frames[rnd] = frame  # apped frame to our list
                    times.append(time.time() - start)  # save the time spend
                    pbar.update(1)
                    break

    with open("../data/" + str(samples) + "_transition_states_" + datetime.datetime.now().strftime("%Y%m%d%H%M%S"), 'wb') as outfile:
        pickle.dump(frames, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    return frames


def transition_matrix(df, ms_centers, lag_time=1):
    if lag_time == 0:
        print("lag_time is 0!")
        return

    tree = cKDTree(ms_centers)
    points = np.c_[df[0], df[1]]  # fancy zipping
    queries = tree.query(points, n_jobs=2)[1]
    trans_mat = np.zeros(shape=(len(ms_centers), len(ms_centers)))
    for idx in range(df.shape[0] - lag_time):
        old_state = queries[idx]
        new_state = queries[idx + lag_time]
        trans_mat[old_state, new_state] += 1

    # count occurences
    occurences = np.zeros(shape=(len(ms_centers)))
    for label in queries:
        if label == -1:
            continue
        occurences[label] += 1

    # normalization
    # T(A,B)=T(A,B)/N(B)
    for row in range(0, len(ms_centers)):
        for col in range(0, len(ms_centers)):
            t_a_b = trans_mat[row, col]
            t_b = occurences[col]
            trans_mat[row, col] = t_a_b / t_b

    time_spent = occurences / occurences.sum()

    return np.array(trans_mat).T, time_spent


def create_voronoi_ridges(vor):
    line_segments = []
    for simplex in vor.ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            line_segments.append([(x, y) for x, y in vor.vertices[simplex]])

    ptp_bound = vor.points.ptp(axis=0)

    center = vor.points.mean(axis=0)
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.any(simplex < 0):
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex

            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[i] + direction * ptp_bound.max()

            line_segments.append([(vor.vertices[i, 0], vor.vertices[i, 1]),
                                  (far_point[0], far_point[1])])

    return line_segments


def run_micro_simulations(frames, n=10, debug=False):
    sims = {}
    pbar = tqdm(total=len(frames)*n, desc="Simulation microstate transitions")
    for idx in frames:
        runs = []
        for run in range(n):
            # copy topol file to simulations folder
            cmd = 'cp ../data/md_long.top ../simulation/short/topol.top'
            proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = proc.communicate()
            if proc.returncode:
                raise Exception(err)

            # copy gro file to simulations folder
            gro = '../data/frames/' + 'f' + str(idx) + '.gro'
            cmd = 'cp ' + gro + ' ../simulation/short/md.gro'
            proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = proc.communicate()
            if proc.returncode:
                raise Exception(err)

            # run the simulation
            cmd = 'bash --login ../simulation/short/short_mpi.sh'
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            if proc.returncode:
                raise Exception(err)

            # clean up afterwards
            if not debug:
                cmd = 'bash --login ../simulation/short/cleanup.sh'
                proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                if proc.returncode:
                    raise Exception(err)

            # copy rama file to data folder
            rama = '../data/xvg/' + 'f' + str(idx) + '.rama'
            cmd = 'cp ../simulation/short/rama.xvg ' + rama
            proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out, err) = proc.communicate()
            if proc.returncode:
                raise Exception(err)

            # read in the xvg and save the data
            runs.append(read_xvg(rama))
            pbar.update(1)

        # create big array with all the data
        sims[idx] = runs

    pbar.close()

    with open('../data/' + str(len(frames)) + '_micro_sims_' + datetime.datetime.now().strftime("%Y%m%d%H%M%S"), 'wb') as outfile:
        pickle.dump(sims, outfile)

    return sims


def run_macro_simulation(debug=False):
    # run the simulation
    cmd = 'bash --login simulation/long/task_mpi.sh'
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    print("STDOUT:", proc.stdout)
    print("STDERR:", proc.stderr)
    print('task_mpi.sh done')
    (out, err) = proc.communicate()
    if proc.returncode:
       raise Exception(err)
    print('new command')
    # copy gro file of the simulation
    import os

    # Get the current working directory
    current_directory = os.getcwd()

    # Print the current working directory
    print("Current Working Directory:", current_directory)
    cmd = 'cp simulation/long/md.gro data/md_long.gro'
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    if proc.returncode:
        raise Exception(err)

    # copy topol file of the simulation
    cmd = 'cp simulation/long/topol.top data/md_long.top'
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    if proc.returncode:
        raise Exception(err)

    # copy xtc trajectory file of the simulation
    cmd = 'cp simulation/long/md_corr.xtc data/md_long_corr.xtc'
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    if proc.returncode:
        raise Exception(err)

    # save rama file to data folder
    cmd = 'cp simulation/long/rama.xvg data/md_long_nojump_rama.xvg'
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    if proc.returncode:
        raise Exception(err)

    # clean up afterwards
    if not debug:
        cmd = 'bash --login simulation/long/cleanup.sh'
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        if proc.returncode:
            raise Exception(err)
