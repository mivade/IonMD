#define CATCH_CONFIG_MAIN

#include <memory>
#include <boost/filesystem.hpp>
#include <ionmd/data.hpp>
#include "catch.hpp"

using namespace ionmd;
namespace fs = boost::filesystem;


TEST_CASE("data can be written", "[data]")
{
    auto filename = "test.h5";
    auto writer = DataWriter(filename, true);

    SECTION("file exists") {
        REQUIRE(fs::exists(filename));
    }

    SECTION("ions group created")
    {

    }

    // fs::remove(filename);
}
