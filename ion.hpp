/*
  This file is part of IonMD.

  IonMD is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  IonMD is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU General Public License
  along with IonMD.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ION_HPP
#define ION_HPP

#include <armadillo>
#include "params.hpp"

using arma::vec;
using arma::dot;
using arma::mat;

struct Ion {
    vec x, v, a;
    double m, Z;  // mass and charge
    int index;
    bool doppler_coolable;  // can it be Doppler cooled?
    Params *p;

    Ion(Params *p, double m, double Z, vec x0);

    void update(double t, mat forces);

    // Pretty printing
    void print_position();

private:
    vec doppler_force();
    vec coulomb_force(mat forces);
    vec secular_force();
    vec micromotion_force(double t);
    vec stochastic_force();
};

#endif
