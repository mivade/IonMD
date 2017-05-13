#include <iostream>
#include <ionmd/data.hpp>

using namespace ionmd;


DataWriter::DataWriter(std::string filename, bool overwrite)
{
    auto mode = overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL;
    file = H5Fcreate(filename.c_str(), mode, H5P_DEFAULT, H5P_DEFAULT);

    // Group for storing ion data; for now, store initial and final data. This
    // lets us reconstruct a simulation from where we stopped since it contains
    // velocities as well as positions.
    create_group("/ions");

    //handles["trajectories"] = H5D
}


DataWriter::~DataWriter()
{
    if (file > 0) {
        H5Fclose(file);
    }
}


auto DataWriter::create_group(std::string name) -> hid_t
{
    auto handle = H5Gcreate(file, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (handle > 0) {
        groups[name] = handle;
    }
    return handle;
}


auto DataWriter::get_group(std::string key) -> hid_t
{
    auto handle = groups.find(key);
    if (handle != groups.end()) {
        return handle->second;
    }
    else {
        return -1;
    }
}


auto DataWriter::get_dataset(std::string key) -> hid_t
{
    auto handle = datasets.find(key);
    if (handle != datasets.end()) {
        return handle->second;
    }
    else {
        return -1;
    }
}


auto DataWriter::store_ions(std::string dset, const std::vector<Ion> &ions) -> void
{
    auto handle = get_dataset(dset);
    if (handle < 0) {

    }
}
