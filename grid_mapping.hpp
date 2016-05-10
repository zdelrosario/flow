#ifndef MAP // Include guard
#define MAP

#include <iostream> // std::cout
#include <cmath>    // pow()
#include <string>   // std::string

#include "optimize.hpp"

/** Create mesh for flat plate problem
 *
 * @tparam V Array-like container, 2-Dimensional
 * 
 * @param Nt  Total vertical mesh points
 * @param Nbl Boundary vertical mesh points
 * @param Mt  Total horizontal mesh points
 * @param v   Array to which we will write the grid points
 * @param pad Number of cells to pad forward and aft of plate
 * 
 * @pre Nbl < Nt
 * @pre size(v) == Nt*Mt
 * @pre for all i, s.t. 0<=i<=Nt*Mt, have size(v[i]) == 2
 * 
 * @post v Contains the x,y coordinates of the grid points
 */
template <typename V>
void make_flat_plate(int Nt, int Nbl, int Mt, V& v, short pad) {
  /* --- SETUP --- */
  // int Nt  = 34; // Total vertical mesh points
  // int Nbl = 22; // Boundary layer vertical points

  float H = 0.1;        // Domain height
  float L = 0.1;        // Plate length
  float y_fm = 5.4e-4;  // Height of fine mesh
  float r_min = 0.2 / (Nbl-2); // Minimum spacing ratio
  float alp = 2;        // Expansion coefficient

  float dx = L / (Mt-2*pad);

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
  
  /* --- BUILD GRID --- */
  // Fine mesh vertical points
  auto y_i = [&y_fm,&Nbl,&k_fm](int i) {
    return y_fm*(exp(k_fm*(i-2)/(Nbl-2))-1)/(exp(k_fm)-1);
  };
  // Coarse mesh vertical points
  auto y_j = [&](int j) {
    return (H-y_fm)*(exp(k_cm*(j-Nbl)/(Nt-Nbl))-1)/(exp(k_cm)-1)+y_fm;
  };

  // Vertical point function
  auto y_mesh = [&](int i) {
    if (i<=Nbl) {
      return y_i(i);
    }
    else {
      return y_j(i);
    }
  };
  // Horizontal point function
  auto x_mesh = [&](int j) {
    if (j==0)
      return float(0);
    else if (j<=pad) {
      float temp = 0;
      for (int i=0; i<j; ++i) 
        temp += pow(alp,pad-i);
      return float(temp * dx);
    }
    else if (j>=Mt-pad-1) {
      float temp = 0;
      for (int i=0; i<pad; ++i) 
        temp += pow(alp,pad-i);
      for (int i=0; i<(j-Mt+pad+1); ++i)
        temp += pow(alp,j-Mt+pad);
      return float((temp + j-pad)*dx);
    }
    else {
      float temp = 0;
      for (int i=0; i<pad; ++i) 
        temp += pow(alp,pad-i);
      return float((temp + j-pad)*dx);
    }
  };

  // Store grid values
  for (int i=0; i<Nt; ++i) {    // vertical index
    for (int j=0; j<Mt; ++j) {  // horizontal index
      // std::cout << v[i*Mt+j][0] << "->";
      // Handle x values
      v[i*Mt+j][0] = x_mesh(j);
      // Handle y values
      v[i*Mt+j][1] = y_mesh(Nt-i-1);
      // std::cout << v[i*Mt+j][0] << std::endl;
    }
  }

}

#endif // MAP
