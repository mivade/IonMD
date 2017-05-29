#ifndef DATA_HPP
#define DATA_HPP

#include <fstream>
#include <string>
#include <memory>
#include <vector>
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

    /// Trajectory data file
    std::ofstream traj_file;

    std::vector<std::ofstream> ion_files;

    const size_t buffer_size;
    unsigned int buffer_pos;

public:
    /// Buffer for storing ion data
    arma::mat buffer;

    /**
     * Initialize data output.
     * @param params
     * @param trap
     * @param ions
     * @param overwrite Overwrite existing data.
     */
    DataWriter(params_ptr params, trap_ptr trap,
               const std::vector<Ion> &ions, bool overwrite=false);

    ~DataWriter();

    inline void write(const unsigned int &ion, const arma::vec &data);

    /**
     * Insert trajectory data.
     */
    inline void update_buffer(const unsigned int &ion, const arma::vec &data)
    {
        for (unsigned int j = 0; j < 3; j++) {
            buffer(buffer_pos, 3*ion + j) = data[j];
        }

        // if (++buffer_pos == buffer_size)
        // {
        //     buffer.save(traj_file, arma::raw_binary);
        //     buffer_pos = 0;
        // }
    }
};

}  // namespace ionmd

#endif
