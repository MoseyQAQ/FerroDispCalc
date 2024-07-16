#!/bin/bash

##################### Test for get_d #####################
# build_neighbor
python3 -B build_neighbor.py

# test get_d
$1/get_d traj.lammpstrj B.dat test_B.disp  0 1 1 > /dev/null
$1/get_d traj.lammpstrj A.dat test_A.disp  0 1 1 > /dev/null
cat test_A.disp test_B.disp > test.disp

# compare the output with the reference
python3 -B ../compare.py ref.disp test.disp disp

# clean up
rm -f test.disp test_A.disp test_B.disp A.dat 
##################### End Test get_d #####################

##################### Test for get_p #####################
# test get_p
$1/get_p traj.lammpstrj BA.dat B.dat test.polar 0 1 1 > /dev/null

# compare the output with the reference
python3 -B ../compare.py ref.polar test.polar polar

# clean up
rm -f test.polar BA.dat B.dat
##################### End Test get_p #####################
