#include <iostream>
#include <vector>

#include <ionmd/simulation.hpp>
#include <ionmd/params.hpp>
#include <ionmd/trap.hpp>

using std::cout;
using std::endl;

using namespace ionmd;


int main(int argc, char *argv[]) {
    cout << "IonMD demo\n" << "==========" << endl;

    cout << "Using default simulation and trap parameters..." << endl;
    auto p = std::make_shared<SimParams>();
    auto trap = std::make_shared<Trap>();

    p->coulomb_enabled = false;
    p->verbosity = 2;

    Simulation sim(*p.get(), *trap.get());

    sim.add_ion(40, 1, {0, 0, -20});
    sim.add_ion(40, 1, {0, 0, 20});

    sim.run();

    return 0;
}
