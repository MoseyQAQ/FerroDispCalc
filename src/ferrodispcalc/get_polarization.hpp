#ifndef GET_POLARIZATION_HPP
#define GET_POLARIZATION_HPP

#include <Eigen/Dense>
#include <vector>
#include "basic.hpp"

Eigen::MatrixXd get_polarization_in_one_frame(Frame frame, 
                                              std::vector<std::vector<int>> ba_neighbor_list, 
                                              std::vector<std::vector<int>> bx_neighbor_list,
                                              std::vector<double> atomic_bec);
std::vector<Eigen::MatrixXd> get_polarization(std::string input_file, 
                                              std::vector<std::vector<int>> ba_neighbor_list, 
                                              std::vector<std::vector<int>> bx_neighbor_list,
                                              std::vector<double> atomic_bec,
                                              std::vector<int> frames_to_read);
#endif