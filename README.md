## BeerMol

Displays the molecular geometry from an .xyz file or a .log file of Gaussian, Orca or Priroda in real time 
in a HyperChem style. It can be useful for geometry optimization jobs. GAMESS / Firefly support coming soon.

Can be run alongside with the calculation using a simple script, for example:

#!/bin/sh  
export GIN="$@"  
export GOUT="${GIN%.*}.log"  
touch $GOUT  
exec python3 ~/beermol/main.py $GOUT &  
exec g09 "$@"

BeerMol requires VTK (Visualisation Toolkit). To install, type in terminal:  

pip install vtk