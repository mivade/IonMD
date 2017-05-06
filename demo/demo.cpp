#include <iostream>
#include <vector>

#include <ionmd/simulation.hpp>
#include <ionmd/params.hpp>
#include <ionmd/trap.hpp>
#include <ionmd/constants.hpp>

using std::cout;
using std::endl;

using namespace ionmd;


int main(int argc, char *argv[])
{
    cout << "IonMD demo\n" << "==========" << endl;

    auto params = SimParams();
    auto trap = Trap();

    params.coulomb_enabled = true;
    params.secular_enabled = true;

    trap.U_ec = 5;

    cout << params.to_string() << endl
         << trap.to_string() << endl;

    Simulation sim(params, trap);

    auto m = double(constants::amu * 40);
    sim.add_ion(m, 1, {0, 0, -20e-6});
    sim.add_ion(m, 1, {0, 0, 20e-6});

    sim.run();

    cout << "Complete!" << endl;

    return 0;
}
