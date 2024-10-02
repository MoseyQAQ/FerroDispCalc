#ifndef GET_AVERAGED_STRUCTURE_HPP
#define GET_AVERAGED_STRUCTURE_HPP

#include <Eigen/Dense>
#include <vector>
#include "basic.hpp"

Eigen::Matrix3d get_avg_cell(std::vector<Frame> frames);
Eigen::MatrixXd get_avg_coords(std::vector<Frame> frames, int natoms);
std::tuple<Eigen::Matrix3d, Eigen::MatrixXd, std::vector<int>> get_averaged_structure(std::string input_file,
                                                                                       std::string type_map_file,
                                                                                       std::string output_file,
                                                                                       std::vector<int> frames_to_read);

#endif