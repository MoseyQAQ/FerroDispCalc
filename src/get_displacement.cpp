#include "get_displacement.hpp"
#include "basic.hpp"

Eigen::MatrixXd get_displacement_in_one_frame(Frame frame, std::vector<std::vector<int>> neighbor_list)
{ 
     // initialize the original coordinates and displacements
    int natoms = neighbor_list.size();
    Eigen::MatrixXd displacements = Eigen::MatrixXd::Zero(natoms, 3);

    // loop over all center atoms
    for (int i = 0; i < neighbor_list.size(); i++) {
        Eigen::RowVector3d center = frame.coords.row(neighbor_list[i][0]);
        Eigen::RowVector3d neighbor = Eigen::RowVector3d::Zero();

        // loop over all neighbors of the center atom
        for (int j = 1; j < neighbor_list[i].size(); j++) {
            Eigen::RowVector3d neighbor_coord = apply_pbc(frame.coords.row(neighbor_list[i][j]), center, frame.cell);
            neighbor += neighbor_coord;
        }

        // calculate the average position of the neighbors
        neighbor /= neighbor_list[i].size() - 1;
        Eigen::RowVector3d displacement = center - neighbor;
        displacements.row(i) = displacement;
    }

    return displacements;
}

std::vector<Eigen::MatrixXd> get_displacement(std::string input_file, std::vector<std::vector<int>> neighbor_list, std::vector<int> frames_to_read)
{   
    int natoms = get_natoms(input_file);
    std::vector<std::streampos> frame_pos = get_frame_positions(input_file);
    Traj traj = read_selected_frames(input_file, natoms, frame_pos, frames_to_read);

    std::vector<Eigen::MatrixXd> data(traj.frames.size());
    for (int i = 0; i < traj.frames.size(); i++) {
        data[i] = get_displacement_in_one_frame(traj.frames[i], neighbor_list);
    }

    return data;
}