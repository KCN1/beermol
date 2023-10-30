## BeerMol

A Python script displaying molecular geometry from .xyz file or .log file of Gaussian or Orca in real time in a HyperChem style. It can be useful for geometry optimization jobs.

Can be run together with a calculation using a simple script, for example:

#!/bin/bash  
export GIN="$@"  
export GOUT="${GIN%.*}.log"  
touch $GOUT  
exec python3 ~/beermol/main.py $GOUT &  
exec g09 "$@"