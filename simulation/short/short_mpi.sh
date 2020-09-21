#!/usr/bin/env bash

# cd to the folder of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${DIR}

# dont keep gmx backups
export GMX_MAXBACKUP=-1

echo "--- final simulation ----"
grompp_mpi -o md.tpr -f md.mdp -c md.gro || exit 1
mpirun -np 8 mdrun_mpi -v -deffnm md || exit 1

echo "--- write phi and psi angles to rama.xvg ---"
g_rama_mpi -f md.trr -s md.tpr || exit 1

echo "--- convert trajectory (otherwise, atom is split on borders of box) ---"
echo "--- also save to xtc instead of trr (smaller file) ---"
trjconv_mpi -pbc nojump -f md.trr -o md_corr.xtc || exit 1
# rm md.trr || exit 1

echo "SUCCESS, results in:"
echo "md_corr.xtc: corrected trajectory"
echo "rama.xvg: plot of dihedral angles"

