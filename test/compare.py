import sys
from com_sys import compare_stru, compare_disp, compare_polar
from pathlib import Path

ref_file = sys.argv[1]
test_file = sys.argv[2]
compare_type = sys.argv[3]

if compare_type == 'stru':
    status, distance = compare_stru(ref_file, test_file)
elif compare_type == 'disp':
    status, distance = compare_disp(ref_file, test_file)
elif compare_type == 'polar':
    status, distance = compare_polar(ref_file, test_file)
else:
    raise ValueError('Unknown compare type: %s' % compare_type)

print(f"{str(Path(test_file).absolute()):<80} {'Pass' if status else 'Fail':10}  Max Error: {distance:.4f}")