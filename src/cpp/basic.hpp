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
void skip_lines(std::ifstream &file, int n);

/* read cell matrix from lammps dump file */
Eigen::Matrix3d read_cell(std::ifstream &file);

/* Read coordinates from lammps dump file */
Eigen::MatrixXd read_coords(std::ifstream &file, int n_atoms);

/* get number of frames */
int get_nframes(std::string filename);

/* get position of each frame */
std::vector<std::streampos> get_frame_positions(std::string filename);

/* get number of atoms */
int get_natoms(std::string filename);

/* read single frame */
Frame read_single_frame(std::ifstream &file, int n_atoms);

/* read atom types in intger */
std::vector<int> read_atom_types(std::string filename, int n_atoms);

/* read selected frames */
Traj read_selected_frames(std::string filename, int n_atoms, std::vector<std::streampos> frame_positions, std::vector<int> select);

/* calculate the neighbor coord after applying the PBC */
Eigen::RowVector3d apply_pbc(Eigen::RowVector3d neighbor, Eigen::RowVector3d center, Eigen::Matrix3d cell);

#endif 