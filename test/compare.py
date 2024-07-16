import sys
from com_sys import compare_stru, compare_disp, compare_polar

ref_file = sys.argv[1]
test_file = sys.argv[2]
compare_type = sys.argv[3]

if compare_type == 'stru':
    compare_stru(ref_file, test_file)
elif compare_type == 'disp':
    compare_disp(ref_file, test_file)
elif compare_type == 'polar':
    compare_polar(ref_file, test_file)
else:
    raise ValueError('Unknown compare type: %s' % compare_type)