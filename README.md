# $\pi^0$ -> $\gamma$ $\gamma$
Quick and dirty analysis script for producing isotropic $\pi^0$ -> $\gamma$ $\gamma$ decay using $\pi^0$ from NEUT and boosting to lab frame

## Instructions
Requires ROOT and NEUT to run. Edit the Makefile variables to point to your locations.

Make simply by typing `make`.

Do `./pi0gammagamma NAME_OF_ROOT_FILE.root`, where `NAME_OF_ROOT_FILE.root` is the location of a NEUT output file.

This will produce `NAME_OF_ROOT_FILE.root_output.root` which contains a `TTree` with the physics variables.
