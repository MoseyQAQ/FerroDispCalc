/*
get_polarization_displacement.cpp:
    This program reads a LAMMPS dump file and a neighbor list file, and calculates the polarization displacement of cations.
    The neighbor list file contains the indices of the neighbors of each cation, starting from 0. It can be generated using 
    "build_neighbor_list.py"
    Eigen library is required for matrix operations.

Compile:
    g++ get_polarization_displacement.cpp -O3 -o get_polarization_displacement -I /path/to/eigen

Usage:
    ./get_polarization_displacement traj_file nl_file output_file start_frame end_frame step
    OPTIONS:
        traj_file: LAMMPS dump file
        nl_file: neighbor list file, it can be generated using "build_neighbor_list.py"
        output_file: output file, each line contains the original coordinates and displacements of a cation.
        start_frame: the first frame to be read, starting from 0
        end_frame: the last frame to be read
        step: read every step frames

Author: Denan Li
Last modified: 2024-07-12
Email: lidenan@westlake.edu.cn

Todo list:
    1. more testing cases
    2. support read in .xsf file
    3. support output an averaged polarization displacement.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <Eigen/Dense>
#include <iomanip>
#include "basic.hpp"

/* define struct for polarization data */
struct polarization_data
{
    Eigen::MatrixXd original_coords; // original coordinates of cations
    Eigen::MatrixXd displacements; // displacements of cations
};

// define functions
std::vector<std::vector<int>> parse_nl_file(std::string nl_file);
polarization_data get_polarization_displacement_in_one_frame(Frame frame, std::vector<std::vector<int>> neighbor_list);

/* Main function */
int main(int argc, char** argv) {
    // setup input parameters
    std::string traj_file = argv[1];
    std::string nl_file = argv[2];
    std::string output_file = argv[3];
    int start_frame = std::stoi(argv[4]);
    int end_frame = std::stoi(argv[5]);
    int step = std::stoi(argv[6]);

    // check whether the output file exists, if so, exit
    std::ifstream check_file(output_file);
    if (check_file.good() == true) {
        std::cerr << "File already exists: " << output_file << std::endl;
        exit(1);
    }

    // read input file
    int natoms = get_natoms(traj_file);
    std::vector<std::streampos> frame_pos = get_frame_positions(traj_file);
    std::vector<std::vector<int>> neighbor_list = parse_nl_file(nl_file);
    
    // print information
    std::cout << "Number of atoms: " << natoms << std::endl;
    std::cout << "Number of frames: " << frame_pos.size() << std::endl;
    
    // read selected frames
    std::vector<Frame> frames = read_selected_frames(traj_file, natoms, frame_pos, start_frame, end_frame, step);
    std::cout << "Frames read: " << frames.size() << std::endl;

    // get polarization displacement data in all frames
    std::vector<polarization_data> data(frames.size());
    for (int i = 0; i < frames.size(); i++) {
        data[i] = get_polarization_displacement_in_one_frame(frames[i], neighbor_list);
    }

    // write output file
    std::ofstream file(output_file);
    file << std::fixed << std::setprecision(16);
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].original_coords.rows(); j++) {
            file << data[i].original_coords.row(j) << " " << data[i].displacements.row(j) << std::endl;
        }
    }

    return 0;
}

/* parse the neighbor list file */
std::vector<std::vector<int>> parse_nl_file(std::string nl_file) {
    // initialize the neighbor list
    std::ifstream file(nl_file);
    std::string line;
    std::vector<std::vector<int>> neighbor_list;

    // read the neighbor list
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<int> lineData((std::istream_iterator<int>(iss)), std::istream_iterator<int>());
        neighbor_list.push_back(lineData);
    }

    return neighbor_list;
}

/* calculate the polarization displacement in one frame */
polarization_data get_polarization_displacement_in_one_frame(Frame frame, std::vector<std::vector<int>> neighbor_list) {
    // initialize the original coordinates and displacements
    int natoms = neighbor_list.size();
    Eigen::MatrixXd original_coords = Eigen::MatrixXd::Zero(natoms, 3);
    Eigen::MatrixXd displacements = Eigen::MatrixXd::Zero(natoms, 3);

    // loop over all center atoms
    for (int i = 0; i < neighbor_list.size(); i++) {
        Eigen::RowVector3d center = frame.coords.row(neighbor_list[i][0]);
        Eigen::RowVector3d neighbor = Eigen::RowVector3d::Zero();

        // loop over all neighbors of the center atom
        for (int j = 1; j < neighbor_list[i].size(); j++) {
            // check pbc
            Eigen::RowVector3d neighbor_coord = frame.coords.row(neighbor_list[i][j]);
            Eigen::RowVector3d diff = neighbor_coord - center;
            Eigen::RowVector3d diff_frac = diff * frame.cell.inverse();
            Eigen::RowVector3d neighbor_coord_frac = neighbor_coord * frame.cell.inverse();

            for (int k = 0; k < 3; k++) {
                if (diff_frac(k) > 0.5) {
                    neighbor_coord_frac(k) -= 1;
                } else if (diff_frac(k) < -0.5) {
                    neighbor_coord_frac(k) += 1;
                }
            }

            neighbor_coord = neighbor_coord_frac * frame.cell;
            neighbor += neighbor_coord;
        }

        // calculate the average position of the neighbors
        neighbor /= neighbor_list[i].size() - 1;
        Eigen::RowVector3d displacement = center - neighbor;
        original_coords.row(i) = center;
        displacements.row(i) = displacement;
    }

    polarization_data data;
    data.original_coords = original_coords;
    data.displacements = displacements;
    return data;
}