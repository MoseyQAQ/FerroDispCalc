g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) basic.cpp get_averaged_structure.cpp binding.cpp -o fdc$(python3-config --extension-suffix)