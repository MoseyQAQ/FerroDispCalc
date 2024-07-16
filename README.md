# FerroDispCalc

Calculator of polarization displacements and polarization of ferroelectrics in molecular dynamics simulations. 

This code is designed for: 
* build neighbor list for ferrorelectrics (ABO3 perovskite)
* calculate the polarization and displacement from a LAMMPS dump file. 

For better performence, the main part of code is written in C++.

## Code Structure

```
.
├── example  //example of how to use code, in todo list
├── ferrodispcalc // python library
│   ├── build_neighbor_list.py  //Build neighbor list for ABO3 perovskite
│   ├── PolarPolter.py
│   └── type_map.py  //store the common type map
├── install.sh //installation script
├── LICENSE
├── README.md
├── src  //c++ source code
│   ├── basic.hpp //This file contains basic functions to read lammps dump files.
│   ├── get_averaged_structure.cpp //Calculate the average structure from LAMMPS dump file
│   ├── get_polarization.cpp //Calculate the polarization of each perovskite unit call using LAMMPS dump file
│   └── get_polarization_displacement.cpp //Calculate the polarization displacement from LAMMPS dump file
└── test //files related to testing
    ├── BaTiO3
    ├── compare.py
    ├── com_sys.py
    ├── PIN-PMN-PT
    ├── PIN-PMN-PT_multi_domain
    ├── PSZO
    ├── README.md
    └── run_test.sh 
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
1. interface with Python3
2. more easy to use
    * in get_a: support output of POSCAR format. support flexible slice of frames
3. doc with spinx