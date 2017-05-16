#ifndef DATA_HPP
#define DATA_HPP

#include <fstream>
#include <string>
#include <memory>
#include <armadillo>
#include "params.hpp"
#include "trap.hpp"
#include "ion.hpp"

namespace ionmd {


/**
 * Class for managing simulation data output.
 */
class DataWriter
{
private:
    /// Path to data files.
    std::string path;

    std::ofstream traj_file;

public:
    /**
     * Initialize data output.
     * @param params
     * @param trap
     * @param ions
     * @param overwrite Overwrite existing data.
     */
    DataWriter(const SimParams &params, const Trap &trap,
               const std::vector<Ion> &ions, bool overwrite=false);

    ~DataWriter();

    /**
     * Write trajectory updates to disk.
     */
    void write_trajectories(const arma::mat &data);
};

}  // namespace ionmd

#endif
