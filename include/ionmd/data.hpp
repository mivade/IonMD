#ifndef DATA_HPP
#define DATA_HPP

#include <string>
#include <vector>
#include <map>
#include "hdf5.h"
#include "ion.hpp"

namespace ionmd {


/**
 * Class for managing simulation data output.
 */
class DataWriter
{
private:
    hid_t file;
    std::map<std::string, hid_t> groups;
    std::map<std::string, hid_t> datasets;

    hid_t create_group(std::string name);

    hid_t get_group(std::string key);
    hid_t get_dataset(std::string key);

public:
    /**
     * Opens and initializes an HDF5 file for output.
     * @param filename Full filename for the data file.
     * @param overwrite Overwrite existing data.
     */
    DataWriter(std::string filename, bool overwrite=false);

    /** Closes file. */
    ~DataWriter();

    /**
     * Write ion settings to disk.
     * @param dset Name of HDF5 dataset
     * @param ions
     */
    void store_ions(std::string dset, const std::vector<Ion> &ions);
};

}  // namespace ionmd

#endif
