#include "basic.hpp"
#include "get_polarization.hpp"

Eigen::MatrixXd get_polarization_in_one_frame(Frame frame, 
                                              std::vector<std::vector<int>> ba_nl, 
                                              std::vector<std::vector<int>> bx_nl,
                                              std::vector<double> atomic_bec)
{
    int natoms = ba_nl.size();
    Eigen::MatrixXd polarization(natoms, 3);

    // loop over all B site cations
    for (int i = 0; i < ba_nl.size(); i++) {
        Eigen::RowVector3d term_a = Eigen::RowVector3d::Zero();
        Eigen::RowVector3d term_b = Eigen::RowVector3d::Zero();
        Eigen::RowVector3d term_x = Eigen::RowVector3d::Zero();
        Eigen::RowVector3d center = frame.coords.row(ba_nl[i][0]);

        // loop over all A site cations
        for (int j = 1; j < ba_nl[i].size(); j++) {

            // check pbc
            Eigen::RowVector3d neighbor_coord = apply_pbc(frame.coords.row(ba_nl[i][j]), center, frame.cell);
            term_a += atomic_bec[ba_nl[i][j]] * neighbor_coord;
        }

        // loop over all X site anions
        for (int j = 1; j < bx_nl[i].size(); j++) {
            Eigen::RowVector3d neighbor_coord = apply_pbc(frame.coords.row(bx_nl[i][j]), center, frame.cell);
            term_x += atomic_bec[bx_nl[i][j]] * neighbor_coord;
        }

        // calculate the polarization
        term_b = atomic_bec[ba_nl[i][0]] * center;
        polarization.row(i) = term_b + term_x * 0.5 + term_a * 0.125;
    }

    double volume = std::abs(frame.cell.determinant());
    double conversion_factor = 1.602176E-19 * 1.0E-10 * 1.0E30; // convert to C/m^2
    polarization = ba_nl.size() * polarization * conversion_factor / volume; // here, we assume the volume is the same for all unit cells.

    return polarization;
}

std::vector<Eigen::MatrixXd> get_polarization(std::string input_file, 
                                              std::vector<std::vector<int>> ba_neighbor_list, 
                                              std::vector<std::vector<int>> bx_neighbor_list,
                                              std::vector<double> atomic_bec,
                                              std::vector<int> frames_to_read)
{
    int natoms = get_natoms(input_file);
    std::vector<std::streampos> frame_pos = get_frame_positions(input_file);
    Traj traj = read_selected_frames(input_file, natoms, frame_pos, frames_to_read);

    std::vector<Eigen::MatrixXd> data(traj.frames.size());
    for (int i = 0; i < traj.frames.size(); i++) {
        data[i] = get_polarization_in_one_frame(traj.frames[i], ba_neighbor_list, bx_neighbor_list, atomic_bec);
    }

    return data;
}