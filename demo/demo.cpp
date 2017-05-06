#include <iostream>
#include <vector>
#include <memory>

#include <ionmd/simulation.hpp>
#include <ionmd/params.hpp>
#include <ionmd/trap.hpp>
#include <ionmd/laser.hpp>
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
    params.doppler_enabled = true;

    trap.U_ec = 5;

    cout << params.to_string() << endl
         << trap.to_string() << endl;

    Simulation sim(params, trap);

    double beta = 2e-22;
    double F0 = 1.3e-19;
    auto laser = Laser(beta, F0, {0, 0, 1});
    lasers_ptr lasers;
    lasers.push_back(std::make_shared<Laser>(laser));

    auto m = double(constants::amu * 40);
    std::vector<Ion> ions;

    // FIXME: do this less stupidly
    auto params_p = std::make_shared<SimParams>(params);
    auto trap_p = std::make_shared<Trap>(trap);
    ions.push_back(Ion(params_p, trap_p, lasers, m, 1, {0, 0, -100e-6}));
    ions.push_back(Ion(params_p, trap_p, lasers, m, 1, {0, 0, 100e-6}));

    sim.set_ions(ions);
    sim.run();

    cout << "Complete!" << endl;

    return 0;
}
