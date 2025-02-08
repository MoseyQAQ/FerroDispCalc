#ifndef GET_DISPLACEMENT_HPP
#define GET_DISPLACEMENT_HPP

#include <Eigen/Dense>
#include "basic.hpp"
#include <vector>
Eigen::MatrixXd get_displacement_in_one_frame(Frame frame, std::vector<std::vector<int>> neighbor_list);
std::vector<Eigen::MatrixXd> get_displacement(std::string input_file, 
                                              std::vector<std::vector<int>> neighbor_list, 
                                              std::vector<int> frames_to_read);

#endif 