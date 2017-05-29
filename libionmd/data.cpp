#include <iostream>
#include <sstream>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <ionmd/data.hpp>

using namespace ionmd;
namespace fs = boost::filesystem;


DataWriter::DataWriter(params_ptr params, trap_ptr trap,
                       const std::vector<Ion> &ions, bool overwrite)
    : path(params->path), buffer_size(params->buffer_size), buffer_pos(0)
{
    // Create output directory
    if (fs::exists(path))
    {
        if (!overwrite) {
            auto msg = "The path " + path + " already exists";
            throw (std::runtime_error(msg));
        }

        if (!fs::is_directory(path)) {
            auto msg = "Path " + path + " is not a directory";
            throw(std::runtime_error(msg));
        }
    }
    else {
        fs::create_directory(path);
    }
    this->path = path;

    // Output sim params
    fs::path p_path = path;
    p_path /= "params.json";
    std::ofstream pout(p_path.c_str());
    pout << params->to_json() << std::endl;
    pout.close();

    // Output trap params
    fs::path trap_path = path;
    trap_path /= "trap.json";
    std::ofstream tout(trap_path.c_str());
    tout << trap->to_json() << std::endl;
    tout.close();

    // Output initial ion positions
    fs::path ions_path = path;
    ions_path /= "ions-init.csv";
    std::ofstream ions_out(ions_path.c_str());
    ions_out << "m,Z,position,velocity,acceleration\n";
    for (const auto &ion: ions) {
        ions_out << ion.m << "," << ion.Z << "\n";
    }
    ions_out.close();

    // Create stream for writing trajectory data
    fs::path traj_filename = path;
    // traj_filename /= "trajectories.bin";
    // traj_file.open(traj_filename.c_str(), std::ios::out | std::ios::binary);

    // Allocate ion data buffer
    buffer.zeros(params->num_steps, 3 * params->buffer_size);
}


DataWriter::~DataWriter()
{
    // traj_file.close();
}


inline void DataWriter::write(const unsigned int &ion, const arma::vec &data)
{
    ion_files[ion] << data;
}
