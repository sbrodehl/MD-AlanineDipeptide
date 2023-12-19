#!/usr/bin/env bash

# cd to the folder of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "${DIR}"
export PATH=$PATH:/usr/local/gromacs/bin

# don't keep gmx backups
export GMX_MAXBACKUP=-1

echo "--- converting pdb file ---"
#gmx pdb2gmx -f ala2.pdb -water tip3p -ff amber03 || exit 1 #also try amber99
gmx pdb2gmx -v -f ala2.pdb -water tip3p -ff amber03 || exit 1


echo "--- setting box size ---"
gmx editconf -o box.gro -f conf.gro -bt cubic -d 1.2 || exit 1

echo "--- creating water molecules ---"
gmx solvate -o sol.gro -cp box.gro -cs spc216.gro -p topol.top || exit 1

echo "--- adding NaCl in physiological concentration ---"
gmx grompp -o iongen.tpr -c sol.gro -f em.mdp || exit 1
echo "13\n" | gmx genion -o ionized.gro -s iongen.tpr -p topol.top -conc 0.15 || exit 1

echo "--- running energy minimization ----"
gmx grompp -o em.tpr -f em.mdp -c ionized.gro || exit 1
#mpirun -np 1 gmx mdrun -deffnm em || exit 1
gmx mdrun -nt 1 -deffnm em || exit 1





echo "--- init temperature ----"
#gmx grompp -o nvt.tpr -f nvt.mdp -c em.gro || exit 1
gmx grompp -o nvt.tpr -f nvt.mdp -c em.gro -p topol.top -r em.gro || exit 1


#mpirun -np 1 gmx mdrun -deffnm nvt || exit 1
gmx mdrun -nt 1 -deffnm nvt || exit 1

echo "--- init pressure ----"
gmx grompp -o npt.tpr -f npt.mdp -c nvt.gro -p topol.top -r nvt.gro || exit 1

#mpirun -np 1 gmx mdrun -deffnm npt || exit 1
gmx mdrun -nt 1 -deffnm npt || exit 1

echo "--- final simulation ----"
gmx grompp -o md.tpr -f md.mdp -c npt.gro || exit 1
# mpirun -np 1 gmx mdrun -v -deffnm md || exit 1
gmx mdrun -nt 1 -v -deffnm md

echo "--- write phi and psi angles to rama.xvg ---"
gmx rama -f md.trr -s md.tpr || exit 1

echo "--- convert trajectory (otherwise, atom is split on borders of box) ---"
echo "0\n" | gmx trjconv -pbc nojump -f md.trr -o md_corr.xtc || exit 1

echo "SUCCESS, results in:"
echo "md.gro: final configuration"
echo "md_corr.xtc: corrected trajectory"
echo "md.edr: energy for trajectory"
echo "md.cpt: data needed to continue simulation"
echo "rama.xvg: plot of dihedral angles"
