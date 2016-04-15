#ifndef MAP // Include guard
#define MAP

#include <iostream>
#include <cmath>

#include "optimize.hpp"

/** Create mesh for flat plate problem
 *
 * @tparam V Array-like container, 2-Dimensional
 * 
 * @param Nt  Total vertical mesh points
 * @param Nbl Boundary vertical mesh points
 * @param Mt  Total horizontal mesh points
 * @param v   Array to which we will write the grid points
 * 
 * @pre Nbl < Nt
 * @pre size(v) == Nt*Mt
 * @pre for all i, s.t. 0<=i<=Nt*Mt, have size(v[i]) == 2
 * 
 * @post v Contains the x,y coordinates of the grid points
 */
template <typename V>
void make_flat_plate(int Nt, int Nbl, int Mt, V& v) {
  /* --- SETUP --- */
  // int Nt  = 34; // Total vertical mesh points
  // int Nbl = 22; // Boundary layer vertical points

  float H = 0.5;        // Domain height
  float L = 0.1;        // Plate length
  float y_fm = 5.4e-4;  // 
  float r_min = 0.2 / (Nbl-2); // Minimum spacing ratio
  float dy_min = r_min*y_fm; // Minimum spacing

  // Space warping operator
  auto y_i = [&y_fm,&Nbl](int i, float k) {
    return y_fm*(exp(k*(i-2)/(Nbl-2))-1)/(exp(k)-1);
  };

  /* --- FIND PARAMETER VALUES --- */
  // Objective function for k_fm
  auto f_fm = [&r_min,&Nbl](float k) {
    return r_min-(exp(k/(Nbl-2))-1)/(exp(k)-1);
  };
  // Find k_fm
  float k_fm = secant(2.0,3.0,f_fm);
  float dy_cm = y_fm*( (exp(k_fm*(Nbl-1)/(Nbl-2))-1)/(exp(k_fm)-1) - 1);

  // Objective function for k_cm
  auto f_cm = [&dy_cm,&H,&y_fm,&Nt,&Nbl](float k) {
    return dy_cm - (H-y_fm)*(exp(k/(Nt-Nbl))-1)/(exp(k)-1);
  };
  // Find k_cm
  float k_cm = secant(7.0,8.0,f_cm);
  

}

#endif // MAP
