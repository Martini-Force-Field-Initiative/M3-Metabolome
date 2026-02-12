## Topology files 

This directory contains topology files in the 
Gromacs included topology style for each metabolite 
in the Martini 3 Metabolome.

The [molecules](molecules) folder contains individual
molecular topologies, listed by resname. If a nucleotide 
topology is used, the 
[`martini_v3.0.0_nucleotide_ffbonded_v1.itp`](bonded/martini_v3.0.0_nucleotide_ffbonded_v1.itp)
contained in the [bonded](bonded) subdirectory must be 
included in the top level topology file **beforehand**. 
For example, to use the ATP topology:

```commandline
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_nucleotide_ffbonded_v1.itp"
#include "ATP.itp"
...
```