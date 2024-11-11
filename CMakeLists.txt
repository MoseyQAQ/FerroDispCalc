cmake_minimum_required(VERSION 3.17)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

# Find Python and pybind11
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Find Eigen
find_package(Eigen3 REQUIRED NO_MODULE)

# Define the sources and headers
set(SOURCES
    src/basic.cpp
    src/get_averaged_structure.cpp
    src/get_displacement.cpp
    src/get_polarization.cpp
    src/binding.cpp
)

set(HEADERS
    src/basic.hpp
    src/get_averaged_structure.hpp
    src/get_displacement.hpp
    src/get_polarization.hpp
)

# Add module
python_add_library(_core MODULE ${SOURCES} WITH_SOABI)
target_link_libraries(_core PRIVATE pybind11::headers Eigen3::Eigen)

# Install the module
install(TARGETS _core DESTINATION ${SKBUILD_PROJECT_NAME})