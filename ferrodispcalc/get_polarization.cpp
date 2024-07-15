/*
get_polarization.cpp:
    This program read a LAMMPS dump file and two neighbor list files, and calculate the polarization of each perovskite unit call.
    The polarization of each unit cell, P, is calculated by the formula (could be found at: PHYSICAL REVIEWB95,094102 (2017)):
        P = (bec * B + bec * 0.5X +  bec * 0.125A)
    where B, X, and A are the coordinates of B site cations, X site anions, and A site cations, respectively.
    The bec of each atom type is read from a bec file.

Compile:
    g++ get_polarization.cpp -o get_polarization -I /path/to/eigen

Usage:
    ./get_polarization input_file a_nl_file x_nl_file output_file start_frame end_frame step
    OPTIONS:
        input_file: LAMMPS dump file
        a_nl_file: neighbor list file for A site cations
        x_nl_file: neighbor list file for X site anions
        output_file: output file, each line contains the polarization of a unit cell.
        start_frame: the first frame to be read, starting from 0
        end_frame: the last frame to be read
        step: read every step frames

Author: Denan Li
Last modified: 2024-07-14
Email: lidenan@westlake.edu.cn

Todo list:
    1. more testing cases
    2. support read in .xsf file
    3. support output an averaged polarization.
    4. flexible IO
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <Eigen/Dense>
#include <iomanip>
#include "basic.hpp"

std::vector<double> get_bec(std::string bec_file);
Eigen::MatrixXd get_polarization_in_one_frame(Frame frame, 
                                              std::vector<std::vector<int>> a_neighbor_list, 
                                              std::vector<std::vector<int>> x_neighbor_list,
                                              std::vector<double> bec_of_atoms,
                                              int total_atoms_num);

int main(int argc, char** argv) {
    std::string input_file = argv[1];
    std::string a_nl_file = argv[2];
    std::string x_nl_file = argv[3];
    std::string output_file = argv[4];
    int start_frame = std::stoi(argv[5]);
    int end_frame = std::stoi(argv[6]);
    int step = std::stoi(argv[7]);
    std::string type_map_file = "type_map_file";
    std::string bec_file = "bec_file";

    int natoms = get_natoms(input_file);
    std::vector<std::streampos> frame_pos = get_frame_positions(input_file);
    std::vector<std::vector<int>> a_neighbor_list = parse_neighbor_list_file(a_nl_file);
    std::vector<std::vector<int>> x_neighbor_list = parse_neighbor_list_file(x_nl_file);
    
    // build bec
    std::vector<int> atom_types = read_atom_types(input_file, natoms);
    std::vector<std::string> type_map = get_type_map(type_map_file);
    std::vector<double> bec = get_bec(bec_file);
    std::vector<double> bec_of_atoms;
    for (int i = 0; i < natoms; i++) {
        bec_of_atoms.push_back(bec[atom_types[i]-1]);
    }

    // read selected frames
    std::vector<Frame> frames = read_selected_frames(input_file, natoms, frame_pos, start_frame, end_frame, step);

    // get polarization data in all frames
    std::vector<Eigen::MatrixXd> polarization(frames.size());
    for (int i = 0; i < frames.size(); i++) {
        polarization[i] = get_polarization_in_one_frame(frames[i], a_neighbor_list, x_neighbor_list, bec_of_atoms, natoms);
    }

    // write output file
    std::ofstream file(output_file);
    file << std::fixed << std::setprecision(16);
    for (int i = 0; i < polarization.size(); i++) {
        for (int j = 0; j < polarization[i].rows(); j++) {
            file << polarization[i].row(j) << std::endl;
        }
    
    }

    return 0;
}

Eigen::MatrixXd get_polarization_in_one_frame(Frame frame, 
                                              std::vector<std::vector<int>> a_neighbor_list, 
                                              std::vector<std::vector<int>> x_neighbor_list,
                                              std::vector<double> bec_of_atoms,
                                              int total_atoms_num) {
    // check if the number of atoms in a_neighbor_list and x_neighbor_list are the same
    if (a_neighbor_list.size() != x_neighbor_list.size()) {
        std::cerr << "The number of atoms in a_neighbor_list and x_neighbor_list are different" << std::endl;
        exit(1);
    }

    int natoms = a_neighbor_list.size();
    
    Eigen::MatrixXd polarization(natoms, 3);

    // loop over all B site cations
    for (int i = 0; i < a_neighbor_list.size(); i++) {
        Eigen::RowVector3d term_a = Eigen::RowVector3d::Zero();
        Eigen::RowVector3d term_b = Eigen::RowVector3d::Zero();
        Eigen::RowVector3d term_x = Eigen::RowVector3d::Zero();
        Eigen::RowVector3d center = frame.coords.row(a_neighbor_list[i][0]);

        // loop over all A site cations
        for (int j = 1; j < a_neighbor_list[i].size(); j++) {

            // check pbc
            Eigen::RowVector3d neighbor_coord = frame.coords.row(a_neighbor_list[i][j]);
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
            term_a += bec_of_atoms[a_neighbor_list[i][j]] * neighbor_coord;
        }

        // loop over all X site anions
        for (int j = 1; j < x_neighbor_list[i].size(); j++) {
            Eigen::RowVector3d neighbor_coord = frame.coords.row(x_neighbor_list[i][j]);
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
            term_x += bec_of_atoms[x_neighbor_list[i][j]] * neighbor_coord;
        }

        // calculate the polarization
        term_b = bec_of_atoms[a_neighbor_list[i][0]] * center;
        polarization.row(i) = term_b + term_x * 0.5 + term_a * 0.125;
    }

    double volume = std::abs(frame.cell.determinant());
    double conversion_factor = 1.602176E-19 * 1.0E-10 * 1.0E30;
    polarization = 0.2 * polarization * total_atoms_num * conversion_factor / volume; // 0.2 here convert total_atoms_num to number of unit cells.

    return polarization;
}

/* parse the bec file */
std::vector<double> get_bec(std::string bec_file) {
    // initialize file stream and variables
    std::ifstream file(bec_file);
    std::string line;
    std::vector<double> bec;

    // parse the file
    std::getline(file, line);
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
        bec.push_back(std::stod(token));
    }

    return bec;
}