####### general seeting #######
EIGEN=/home/liulab/eigen-3.4.0
CXX=g++
CXXFLAGS="-O3" # or "-O3 -march=native" for example
INSTALL_DIR=$(pwd) # default to the current directory
###############################

####### compile #######
function make() {
    echo "Compiling ..."
    mkdir -p build/bin
    $CXX $CXXFLAGS -I $EIGEN -o build/bin/get_averaged_structure src/get_averaged_structure.cpp --std=c++11
    $CXX $CXXFLAGS -I $EIGEN -o build/bin/get_polarization_displacement src/get_polarization_displacement.cpp --std=c++11
    $CXX $CXXFLAGS -I $EIGEN -o build/bin/get_polarization src/get_polarization.cpp --std=c++11

    cd build/bin
    ln -s get_averaged_structure get_a
    ln -s get_polarization_displacement get_d
    ln -s get_polarization get_p
    cd ../..

    cp -r ferrodispcalc build/
}

####### clean #######
function clean() {
    echo "Cleaning ..."
    rm -rf build
}

####### install #######
function install() {
    echo "Installing the package into: $INSTALL_DIR"
    cp -r build/* $INSTALL_DIR
    echo "export PATH=$INSTALL_DIR/bin:\$PATH" >> ~/.bashrc
    echo "export PYTHONPATH=$INSTALL_DIR:\$PYTHONPATH" >> ~/.bashrc
    echo "re-open the terminal or run 'source ~/.bashrc' to make the changes effective"
}

####### testing #######
# test the compiled program

####### main #######
if [ $# -eq 0 ]; then
    echo "Usage: ./install.sh [make|clean|install|test]"
elif [ $# -eq 1 ]; then
    if [ $1 == "make" ]; then
        make
    elif [ $1 == "install" ]; then
        install
    elif [ $1 == "test" ]; then
        #test, to be implemented
        echo "Testing ..."
    elif [ $1 == "clean" ]; then
        clean
    else
        echo "Usage: ./install.sh [make|clean|install|test]"
    fi
else
    echo "Usage: ./install.sh [make|clean|install|test]"
fi