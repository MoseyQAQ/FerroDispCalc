#ifndef GET_DISPLACEMENT_HPP
#define GET_DISPLACEMENT_HPP

#include <Eigen/Dense>
#include "basic.hpp"

/* define struct for polarization data */
struct displacement_data
{
    Eigen::MatrixXd original_coords; // original coordinates of cations
    Eigen::MatrixXd displacements; // displacements of cations
};

displacement_data get_displacement_in_one_frame(Frame frame, std::vector<std::vector<int>> neighbor_list);
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> get_displacement(std::string input_file, std::vector<std::vector<int>> neighbor_list, std::vector<int> frames_to_read);

#endif 