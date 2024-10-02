#ifndef BASIC_HPP
#define BASIC_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <iterator>
#include <unordered_map>

/* --------------------------------- data strcture --------------------------------- */
struct Frame {
    Eigen::Matrix3d cell;
    Eigen::MatrixXd coords;
};

struct Traj {
    std::vector<Frame> frames;
    std::vector<int> atom_types;
};

/* --------------------------------- lmp-dump IO --------------------------------- */

/* skip n lines in the file */
void skip_lines(std::ifstream &file, int n) {
    std::string line;
    for (int i = 0; i < n; ++i) {
        std::getline(file, line);
    }
}

/* read cell matrix from lammps dump file */
Eigen::Matrix3d read_cell(std::ifstream &file) {
    Eigen::Matrix3d cell;
    std::string line;
    std::vector<double> bounds;

    // Read three lines for x, y, z dimensions
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        std::vector<double> lineData((std::istream_iterator<double>(iss)), std::istream_iterator<double>());
        bounds.insert(bounds.end(), lineData.begin(), lineData.end());
    }

    // Extract bounds and tilt factors from the read data
    double xlo_bound = bounds[0], xhi_bound = bounds[1], xy = bounds[2];
    double ylo_bound = bounds[3], yhi_bound = bounds[4], xz = bounds[5];
    double zlo_bound = bounds[6], zhi_bound = bounds[7], yz = bounds[8];

    // Calculate lo and hi values adjusted for tilt (triclinic corrections)
    double xlo = xlo_bound - std::min({0.0, xy, xz, xy + xz});
    double xhi = xhi_bound - std::max({0.0, xy, xz, xy + xz});
    double ylo = ylo_bound - std::min(0.0, yz);
    double yhi = yhi_bound - std::max(0.0, yz);

    // Set the matrix values
    cell(0, 0) = xhi - xlo;
    cell(0, 1) = 0;
    cell(0, 2) = 0;
    cell(1, 0) = xy;
    cell(1, 1) = yhi - ylo;
    cell(1, 2) = 0;
    cell(2, 0) = xz;
    cell(2, 1) = yz;
    cell(2, 2) = zhi_bound - zlo_bound;

    return cell;
}

/* Read coordinates from lammps dump file */
Eigen::MatrixXd read_coords(std::ifstream &file, int n_atoms) {
    Eigen::MatrixXd coords(n_atoms, 3); // 3 columns for x, y, z
    std::string line;
    for (int i = 0; i < n_atoms; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        std::vector<double> lineData((std::istream_iterator<double>(iss)), std::istream_iterator<double>());
        coords(i, 0) = lineData[2];
        coords(i, 1) = lineData[3];
        coords(i, 2) = lineData[4];
    }
    return coords;
}

/* get number of frames */
int get_nframes(std::string filename) {
    std::ifstream file(filename);
    std::string line;
    int nframes = 0;
    while (std::getline(file, line)) {
        if (line.find("ITEM: TIMESTEP") != std::string::npos) {
            nframes++;
        }
    }
    return nframes;
}

/* get position of each frame */
std::vector<std::streampos> get_frame_positions(std::string filename) {
    // initialize file stream
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    // define variables
    std::string line;
    std::vector<std::streampos> frame_positions;

    // write positions of each frame to a vector
    // BE CAREFUL: THE "ITEM: TIMESTEP" LINE HAS BEEN SKIPPED
    while (std::getline(file, line)) {
        if (line.find("ITEM: TIMESTEP") != std::string::npos) {
            frame_positions.push_back(file.tellg());
        }
    }

    return frame_positions;
}

/* get number of atoms */
int get_natoms(std::string filename) {
    std::ifstream file(filename);
    std::string line;
    skip_lines(file, 3);
    std::getline(file, line);
    int n_atoms = std::stoi(line);

    return n_atoms;
}

/* read single frame */
Frame read_single_frame(std::ifstream &file, int n_atoms) {
    Frame frame;
    
    // Since we are already at the beginning of the frame, 
    // we only need to skip 4 lines to read the cell and coordinates
    skip_lines(file, 4);
    frame.cell = read_cell(file);
    skip_lines(file, 1);
    frame.coords = read_coords(file, n_atoms);
    
    return frame;
}

/* read atom types in intger */
std::vector<int> read_atom_types(std::string filename, int n_atoms) {
    // The atom types here is integer, not string
    std::ifstream file(filename);
    std::string line;
    std::vector<int> atom_types(n_atoms);
    skip_lines(file, 9);

    for (int i = 0; i < n_atoms; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        std::vector<int> lineData((std::istream_iterator<int>(iss)), std::istream_iterator<int>());
        atom_types[i] = lineData[1];
    }

    return atom_types;
}

/* read selected frames */
Traj read_selected_frames(std::string filename, int n_atoms, std::vector<std::streampos> frame_positions, std::vector<int> select) {

    // initialize file stream and vector of frames
    std::ifstream file(filename);
    std::vector<Frame> frames;
    Traj traj;

    // read selected frames
    for (int i = 0; i < select.size(); i++) {
        file.seekg(frame_positions[select[i]]);
        Frame frame = read_single_frame(file, n_atoms);
        frames.push_back(frame);
    }
    traj.frames = frames;

    // read atom types
    traj.atom_types = read_atom_types(filename, n_atoms);

    return traj;
}

#endif 