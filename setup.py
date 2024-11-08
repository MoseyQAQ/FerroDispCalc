from setuptools import setup, find_namespace_packages, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
from glob import glob
import os
import sys

file_list = [
    'src/ferrodispcalc/basic.cpp',
    'src/ferrodispcalc/get_averaged_structure.cpp',
    'src/ferrodispcalc/get_displacement.cpp',
    'src/ferrodispcalc/get_polarization.cpp',
    'src/ferrodispcalc/binding.cpp',

]

def get_include_path():
    python_executable_path = sys.executable
    env_base = os.path.dirname(os.path.dirname(python_executable_path))
    include_path = os.path.join(env_base, 'include', 'eigen3')
    if os.path.exists(include_path):
        return include_path
    else:
        return None
    
eigen_include_path = [get_include_path()] if get_include_path() is not None else []

fdc_ext_models = [
    Pybind11Extension(
        "ferrodispcalc.core",
        file_list,
        include_dirs=eigen_include_path,
    )
]

setup(
    ext_modules=fdc_ext_models,
    cmdclass={'build_ext': build_ext},
)