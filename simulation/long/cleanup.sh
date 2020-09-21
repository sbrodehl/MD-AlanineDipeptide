#!/usr/bin/env bash

# cd to the folder of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${DIR}

rm *.cpt *.top \#*.top.*\# *.gro \#*.gro.*\# *.itp \#*.itp.*\# *.tpr \#*.tpr.*\# mdout.mdp \#mdout.mdp.*\# *.edr \#*.edr.*\# *.trr \#*.trr.*\# *.log \#*.log.*\# > /dev/null 2>&1

exit 0
