# FerroDispCalc

Calculator of polarization displacements and polarization of ferroelectrics in molecular dynamics simulations. \
For better performence, the code is written in C++.

## Code Structure

```
.
├── example   // example of how to use code, in todo list
├── ferrodispcalc  // main folder
│   ├── basic.hpp   // This file contains basic functions to read lammps dump files.
│   ├── build_neighbor_list.py  // Build neighbor list for ABO3 perovskite
│   ├── get_averaged_structure.cpp  // Calculate the average structure from LAMMPS dump file
│   ├── get_polarization.cpp  // Calculate the polarization of each perovskite unit call using LAMMPS dump file
│   ├── get_polarization_displacement.cpp  // Calculate the polarization displacement from LAMMPS dump file
│   ├── PolarPolter.py
│   └── type_map.py  // store the common type map
├── install.sh  // installation script
├── LICENSE
├── README.md
└── test  // files related to testing
    ├── BTO
    │   ├── 1.lammpstrj
    │   ├── Main.py
    │   └── parameter.py
    ├── BTO_large_cell
    │   └── head.lammpstrj
    ├── PIN-PMN-PT
    │   ├── 1.lammpstrj
    │   ├── Main.py
    │   └── parameter.py
    └── PSZO
        ├── 1.lammpstrj
        ├── Main.py
        └── parameter.py
```

## Installation

1. Prerequisites:
    * C++ compiler
    * Eigen (For matrix operation)

2. Edit the following variables in <kbd>install.sh</kbd>:
```Bash
####### general seeting #######
EIGEN=/home/liulab/eigen-3.4.0  # path to your eigen header
CXX=gcc                         # your c++ compiler
CXXFLAGS="-O3" # or "-O3 -march=native" for example
INSTALL_DIR=$(pwd)
###############################
```

3. Run the <kbd>install.sh</kbd>:
```Bash
chmod +x install.sh
./install.sh make    # compile the source code
./install.sh clean   # clean the compiled code
./install.sh test    # test the executable
./install.sh install # installation
```
## Example

## Todo list
0. testing case
1. interface with Python3
2. more easy to use
3. doc with spinx