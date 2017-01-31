#include <iostream>
#include <vector>

#include <armadillo>

#include <ionmd/ionmd.hpp>
#include <ionmd/params.hpp>
#include <ionmd/trap.hpp>
#include <ionmd/ion.hpp>

using std::cout;
using std::endl;
using arma::vec;

using namespace ionmd;

typedef std::vector<Ion> ions_t;


ions_t make_ions(params_ptr p, trap_ptr trap) {
    ions_t ions;
    ions.push_back(Ion(p, trap, 40., 1., vec{{0., 0., -20e-6}}));
    ions.push_back(Ion(p, trap, 40., 1., vec{{0., 0., 20e-6}}));
    return ions;
}


int main(int argc, char *argv[]) {
    cout << "IonMD demo\n" << "==========" << endl;

    cout << "Using default simulation and trap parameters..." << endl;
    auto p = std::make_shared<SimParams>();
    auto trap = std::make_shared<Trap>();

    p->verbosity = 2;

    auto ions = make_ions(p, trap);
    Simulation sim(*p.get(), *trap.get(), ions);

    cout << "Running simulation..." << endl;
    sim.run();
    cout << "Simulation complete." << endl;

    return 0;
}
