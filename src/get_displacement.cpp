#include "get_displacement.hpp"
#include "basic.hpp"

displacement_data get_displacement_in_one_frame(Frame frame, std::vector<std::vector<int>> neighbor_list)
{
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
            Eigen::RowVector3d neighbor_coord = apply_pbc(frame.coords.row(neighbor_list[i][j]), center, frame.cell);
            neighbor += neighbor_coord;
        }

        // calculate the average position of the neighbors
        neighbor /= neighbor_list[i].size() - 1;
        Eigen::RowVector3d displacement = center - neighbor;
        original_coords.row(i) = center;
        displacements.row(i) = displacement;
    }

    displacement_data data;
    data.original_coords = original_coords;
    data.displacements = displacements;
    return data;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> get_displacement(std::string input_file, std::vector<std::vector<int>> neighbor_list, std::vector<int> frames_to_read)
{   
    /*
    We should change the func:
    1. now, after calculating disp, we need to store the original coords and displacements in a struct
    2. to return, we nned to store the original coords and displacements in a matrix
    3. this is a bit tricky, because we need to store the original coords and displacements in a matrix for all frames
    */
    int natoms = get_natoms(input_file);
    std::vector<std::streampos> frame_pos = get_frame_positions(input_file);
    Traj traj = read_selected_frames(input_file, natoms, frame_pos, frames_to_read);

    std::vector<displacement_data> data(traj.frames.size());
    for (int i = 0; i < traj.frames.size(); i++) {
        data[i] = get_displacement_in_one_frame(traj.frames[i], neighbor_list);
    }

    Eigen::MatrixXd original_coords = Eigen::MatrixXd::Zero(traj.frames.size() * natoms, 3);
    Eigen::MatrixXd displacements = Eigen::MatrixXd::Zero(traj.frames.size() * natoms, 3);
    for (int i = 0; i < traj.frames.size(); i++) {
        original_coords.block(i * natoms, 0, natoms, 3) = data[i].original_coords;
        displacements.block(i * natoms, 0, natoms, 3) = data[i].displacements;
    }

    return std::make_tuple(original_coords, displacements);
}