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

#ifndef MINIMIZE_HPP
#define MINIMIZE_HPP

//#include <vector>
#include "ionmd.hpp"

//using std::vector;

// Struct for passing data to minimization function
typedef struct MinData {
    int i, fcalls;
    Params *p;
    Ion **ions;
} MinData;

// Old-style minimization by treating all ions as laser cooled
void minimize(Ion **ions, Params *p);

// Minimization main controlling function for nlopt
void minimize(double *x0, Ion **ions, Params *p);

// Writes initial positions to a file
void writeInitPos(Ion **ions, Params *p);
    
// Minimization function for nlopt
//double minfunc(const vector<double> &x, vector<double> &grad, void *_data);

// Potential energy functions called by minfunc
double UTrap(int i, Ion **ions, Params *p);
double ULaser(int i, Ion **ions, Params *p);
double UCoulomb(int i, Ion **ions, Params *p);

#endif
