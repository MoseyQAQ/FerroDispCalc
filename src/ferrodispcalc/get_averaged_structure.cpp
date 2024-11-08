#include "get_averaged_structure.hpp"
#include "basic.hpp"

Eigen::Matrix3d get_avg_cell(std::vector<Frame> frames) {
    Eigen::Matrix3d avg_cell = Eigen::Matrix3d::Zero();
    for (int i = 0; i < frames.size(); ++i) {
        avg_cell += frames[i].cell;
    }
    avg_cell /= frames.size();
    return avg_cell;
}

Eigen::MatrixXd get_avg_coords(std::vector<Frame> frames, int natoms) {
    // initialize the averaged coordinates
    Eigen::MatrixXd avg_coords = frames[0].coords;
    Eigen::MatrixXd coord_init = frames[0].coords;

    for (int i = 1; i < frames.size(); ++i) {

        // check the periodic boundary conditions
        Eigen::MatrixXd coord_current = frames[i].coords;
        Eigen::MatrixXd coord_current_frac = coord_current * frames[i].cell.inverse();
        Eigen::MatrixXd coord_diff = coord_current - coord_init;
        Eigen::MatrixXd coord_diff_frac = coord_diff * frames[i].cell.inverse();

        for (int j = 0; j < natoms; ++j) {
            for (int k = 0; k < 3; ++k) {
                if (coord_diff_frac(j, k) > 0.5) {
                    coord_current_frac(j, k) -= 1;
                } else if (coord_diff_frac(j, k) < -0.5) {
                    coord_current_frac(j, k) += 1;
                }
            }
        }

        coord_current = coord_current_frac * frames[i].cell;
        avg_coords += coord_current;
    }

    avg_coords /= frames.size();
    return avg_coords;
}

std::tuple<Eigen::Matrix3d, Eigen::MatrixXd, std::vector<std::string>> get_averaged_structure(std::string input_file,
                                                                                       std::vector<std::string> type_map,
                                                                                       std::vector<int> frames_to_read) {
    // read input file
    int natoms = get_natoms(input_file);
    std::vector<int> atom_type = read_atom_types(input_file, natoms);
    std::vector<std::streampos> frame_pos = get_frame_positions(input_file);

    // map atom_type int to type_map
    std::vector<std::string> symbols;
    for (int i = 0; i < atom_type.size(); ++i) {
        symbols.push_back(type_map[atom_type[i]-1]);
    }

    // read selected frames
    Traj traj = read_selected_frames(input_file, natoms, frame_pos, frames_to_read);
    Eigen::Matrix3d avg_cell = get_avg_cell(traj.frames);
    Eigen::MatrixXd avg_coords = get_avg_coords(traj.frames, natoms);

    return std::make_tuple(avg_cell, avg_coords, symbols);
}