#!/bin/bash

# Test get_a
$1/get_a traj.lammpstrj test.xsf type_map_file 0 10 1 > /dev/null

# compare the output with the reference
python3 -B ../compare.py ref.xsf test.xsf stru

# clean up
rm -f test.xsf