#!/usr/bin/env bash

# cd to the folder of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${DIR}

# dont keep gmx backups
export GMX_MAXBACKUP=-1

echo "--- converting pdb file ---"
pdb2gmx_mpi -f ala2.pdb -water tip3p -ff amber03 || exit 1 #also try amber99

echo "--- setting box size ---"
editconf_mpi -o box.gro -f conf.gro -bt cubic -d 1.2 || exit 1

echo "--- creating water molecules ---"
genbox_mpi -o sol.gro -cp box.gro -cs spc216.gro -p topol.top || exit 1

echo "--- adding NaCl in physiologial concentration ---"
grompp_mpi -o iongen.tpr -c sol.gro -f em.mdp || exit 1
echo "13\n" | genion_mpi -o ionized.gro -s iongen.tpr -p topol.top -conc 0.15 || exit 1

echo "--- running energy minimization ----"
grompp_mpi -o em.tpr -f em.mdp -c ionized.gro || exit 1
mpirun -np 8 mdrun_mpi -deffnm em || exit 1

echo "--- init temperature ----"
grompp_mpi -o nvt.tpr -f nvt.mdp -c em.gro || exit 1
mpirun -np 8 mdrun_mpi -deffnm nvt || exit 1

echo "--- init pressure ----"
grompp_mpi -o npt.tpr -f npt.mdp -c nvt.gro || exit 1
mpirun -np 8 mdrun_mpi -deffnm npt || exit 1

echo "--- final simulation ----"
grompp_mpi -o md.tpr -f md.mdp -c npt.gro || exit 1
mpirun -np 8 mdrun_mpi -v -deffnm md || exit 1

echo "--- write phi and psi angles to rama.xvg ---"
g_rama_mpi -f md.trr -s md.tpr || exit 1

echo "--- convert trajectory (otherwise, atom is split on borders of box) ---"
echo "--- also save to xtc instead of trr (smaller file) ---"
trjconv_mpi -pbc nojump -f md.trr -o md_corr.xtc || exit 1
# rm md.trr || exit 1

echo "SUCCESS, results in:"
echo "md.gro: final configuration"
echo "md_corr.xtc: corrected trajectory"
echo "md.edr: energy for trajectory"
echo "md.cpt: data needed to continue simulation"
echo "rama.xvg: plot of dihedral angles"

