## BeerMol

Displays molecular geometry from an .xyz file or a .log file of Gaussian, Orca or Priroda in real time 
in HyperChem style. It can be useful for geometry optimization jobs. GAMESS / Firefly support coming soon.

Can be run alongside with a calculation using a simple script, for example:

#!/bin/sh  
export GIN="$@"  
export GOUT="${GIN%.*}.log"  
touch $GOUT  
exec python3 ~/beermol/main.py $GOUT &  
exec g09 "$@"

BeerMol requires VTK, SciPy and PyQT. To install, type in terminal:  

pip install -r requirements.txt