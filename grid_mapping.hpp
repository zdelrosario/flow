#ifndef MAP // Include guard
#define MAP

#include <iostream>
#include <cmath>

#include "optimize.hpp"

void make_flat_plate() {
  /* --- SETUP --- */
  int Nt  = 34; // Total vertical mesh points
  int Nbl = 22; // Boundary layer vertical points

  float eps = 1e-9;

  float L = 0.1;        // Plate length
  float y_fm = 5.4e-4;  // 
  float r_min = 0.2 / (Nbl-2); // Minimum spacing ratio
  float dy_min = r_min*y_fm; // Minimum spacing

  // 
  auto y_i = [&y_fm,&Nbl](int i, float k) {
              return y_fm*(exp(k*(i-2)/(Nbl-2))-1)/(exp(k)-1);
            };
  // 
  auto f = [&r_min,&Nbl](float k) {
              return r_min-(exp(k/(Nbl-2))-1)/(exp(k)-1);
            };
  // 
  auto df= [&Nbl](float k) {
              return -(exp(k/(Nbl-2))/(Nbl-2) 
                       -exp(k)*(exp(k/(Nbl-2))-1)/pow(exp(k)-1,2)) / 
                      (exp(k)-1);
            };
  /* --- FIND PARAMETER VALUES --- */
  // Find K
  float k = secant(2.0,3.0,f);
  // DEBUG
  std::cout << "k=" << k << std::endl;
  std::cout << "f(k)=" << f(k) << std::endl;

}

#endif // MAP
