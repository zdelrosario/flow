#ifndef MAP // Include guard
#define MAP

#include <iostream>
#include <cmath>
#include <string>

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

  float H = 0.1;        // Domain height
  float L = 0.1;        // Plate length
  float y_fm = 5.4e-4;  // Height of fine mesh
  float r_min = 0.2 / (Nbl-2); // Minimum spacing ratio
  float dy_min = r_min*y_fm; // Minimum spacing

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

  // Piecewise vertical points
  auto y_mesh = [&](int i) {
    if (i<=Nbl) {
      return y_i(i);
    }
    else {
      return y_j(i);
    }
  };
  // Store grid values
  for (int i=0; i<Nt; ++i) {
    for (int j=0; j<Mt; ++j) {
      // std::cout << v[i*Mt+j][0] << "->";
      v[i*Mt+j][0] = L/Mt * j;  // X-value
      v[i*Mt+j][1] = y_mesh(i); // Y-value
      // std::cout << v[i*Mt+j][0] << std::endl;
    }
  }

}

/** Grid file writeout
 * @brief Writes a set of x,y points to a 
 *        formatted data file
 *
 * @tparam V Array-like container, 2-Dimensional
 * 
 * @param outputfile String which defines output filename
 * @param v Vector which defines gridpoints
 * 
 * @post A file with name 'outputfile' is written to the
 *       local directory with the gridpoints
 */
template <typename V>
void writeout(std::string outputfile, V& v) {
  std::ofstream f_out(outputfile.c_str());
  // f_out.precision(5);
  // Write out elements
  for (auto it=v.begin(); it!=v.end(); ++it) {
    f_out << (*it)[0] << "," << (*it)[1] << std::endl;
  }
}

#endif // MAP
