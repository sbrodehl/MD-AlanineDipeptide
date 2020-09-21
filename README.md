# MD Simulations of the Alanine Dipeptide
Molecular dynamics simulation and analysis of conformational states and energy landscape of the alanine dipeptide.

![][logo]

We calculate Ramachandran plots for the alanine dipeptide dihedral angle distributions to find stable conformations using GROMACS.
Transition states are then further analysed to refine the energy landscape. 

### Prerequisites

- [Python3](https://www.python.org)
- [GROMACS](http://www.gromacs.org) 4.6.7\*\*
- [MDAnalysis](http://www.mdanalysis.org/)
- [Matplotlib](http://matplotlib.org)
- [NumPy](http://www.numpy.org)
- [SciPy](https://scipy.org)
- [SymPy](http://www.sympy.org)
- [scikit-learn](http://scikit-learn.org/)
- [pandas](http://pandas.pydata.org)
- [tqdm](https://github.com/tqdm/tqdm)

\*\*\) newer version can be used, but scripts need to be changed

To reproduce this work run the script [md.py](md.py), which runs all needed simulations and creates plots and saves intermediate output.


## Walkthrough

1. First, we run a macro simulation with GROMACS defined in [AlanineDipeptide/Computations.py#L177 (Python)](AlanineDipeptide/Computations.py#L177) which executes [simulation/long/task_mpi.sh (GROMACS)](simulation/long/task_mpi.sh).
    We can look at the dihedral angles ϕ, Ψ at every step of the simulation and thus compute their distributions.

    ![][ramachandran-plot]

    One can directly see, that there are some "preferred" ϕ-Ψ pairs, which are the conformational states of the alanine dipeptide.
    These are the angles, which are geometrically preferable, given the molecule structure.

2. One can compute the amount of conformational states and their positions using a clustering algorithm (such as [MeanShift](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html)) - black dots - and the regions between those states using [Voronoi](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Voronoi.html) ridges - black lines.
    
    Of course, this is only an approximation.
    A longer simulation will produce slightly different statistics and therefore other (more) conformal states and transition regions.

    ![][ramachandran-cluster-voronoi]

3. An interesting question one can now ask is:
    > What happens to the geometry of the molecule inbetween conformational states?

    To get an answer we sample states - red dots - from the _transition regions_ and run lots of small simulations using these states as starting configurations.

    ![][ramachandran-transitions]
    
    The colored lines show the _trajectories_ taken.
    Most of the states are drawn to their _closest_ conformational state, but sometimes they cross the _border_ and take on another conformational state.
    
    With this (generated sample) data one can create a [Transition Matrix](https://en.wikipedia.org/wiki/Transition_rate_matrix) T_ij which in our case is an 5x5 matrix:
    
    |         |   S1  |   S2  |   S3  |   S4  |   S5  |
    | ------- | ----- | ----- | ----- | ----- | ----- |
    |  **S1** | 0.116 | 0.067 | 0.064 | 0.750 | 0.004 |
    |  **S2** | 0.021 | 0.638 | 0.337 | 0.004 | 0.001 |
    |  **S3** | 0.019 | 0.314 | 0.472 | 0.195 | 0.000 |
    |  **S4** | 0.061 | 0.001 | 0.054 | 0.882 | 0.001 |
    |  **S5** | 0.059 | 0.026 | 0.007 | 0.280 | 0.628 |
    
    which describes the transition from state i to state j.
    This approximation can then be used in the [master equation](https://en.wikipedia.org/wiki/Master_equation) (which describes the time evolution of a system).

## References

Here is a (probably incomplete) list with related material, which was used in this project. 

- Cornell, Wendy D., et al. **"A second generation force field for the simulation of proteins, nucleic acids, and organic molecules."** Journal of the American Chemical Society 117.19 (1995): 5179-5197.
- Hornak, Viktor, et al. **"Comparison of multiple Amber force fields and development of improved protein backbone parameters."** Proteins: Structure, Function, and Bioinformatics 65.3 (2006): 712-725.
- Lindorff‐Larsen, Kresten, et al. **"Improved side‐chain torsion potentials for the Amber ff99SB protein force field."** Proteins: Structure, Function, and Bioinformatics 78.8 (2010): 1950-1958.
- Buchete, Nicolae-Viorel, and Gerhard Hummer. **"Coarse master equations for peptide folding dynamics."** The Journal of Physical Chemistry B 112.19 (2008): 6057-6069.
- Beauchamp, Kyle A., et al. **"MSMBuilder2: modeling conformational dynamics on the picosecond to millisecond scale."** Journal of chemical theory and computation 7.10 (2011): 3412-3419.
- Sriraman, Saravanapriyan, Ioannis G. Kevrekidis, and Gerhard Hummer. **"Coarse master equation from Bayesian analysis of replica molecular dynamics simulations."** The Journal of Physical Chemistry B 109.14 (2005): 6479-6484.
- Seeber, Michele, et al. **"Wordom: a user‐friendly program for the analysis of molecular structures, trajectories, and free energy surfaces."** Journal of computational chemistry 32.6 (2011): 1183-1194.
- Van der Spoel, D., E. Lindahl, and B. Hess. **"GROMACS User Manual version 4.6. 7."** (2014).
- Humphrey, William, Andrew Dalke, and Klaus Schulten. **"VMD: visual molecular dynamics."** Journal of molecular graphics 14.1 (1996): 33-38.
- Michaud‐Agrawal, Naveen, et al. **"MDAnalysis: a toolkit for the analysis of molecular dynamics simulations."** Journal of computational chemistry 32.10 (2011): 2319-2327.


[logo]: https://github.com/sbrodehl/MD-AlanineDipeptide/raw/master/plots/Alanine_Dipeptide_3D_Logo.png "Alanine Dipeptide"
[ramachandran-plot]: plots/Alanine_Dipeptide_Ramachandran_Plot.png "Ramachandran plot"
[ramachandran-cluster-voronoi]: plots/Alanine_Dipeptide_Cluster_Voronoi.png "Ramachandran plot with cluster centers and voronoi regions"
[ramachandran-transitions]: plots/Alanine_Dipeptide_State_Transitions.png "Ramachandran plot with cluster centers, voronoi regions and transition states and their trajectories"
