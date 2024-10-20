cmake_minimum_required(VERSION 3.14)
project(
    ferrodispcalc
    VERSION 0.1
    LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Python and pybind11
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Include directories for Eigen and project src
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)

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
pybind11_add_module(ferrodispcalc MODULE ${SOURCES} ${HEADERS})

# Link libraries if necessary (Eigen is header-only, usually no need to link)
# target_link_libraries(ferrodispcalc PRIVATE Eigen3::Eigen)

# Install the module
install(TARGETS ferrodispcalc LIBRARY DESTINATION ferrodispcalc)
