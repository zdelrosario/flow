#include "structured_grid.hpp"
#include "grid_mapping.hpp"

#include <iostream>
#include <array>
#include <string>
#include <fstream>

using coord = std::array<float,2>; // Grid coordinates
using value = std::array<int,2>;   // State vector values

typedef StructuredGrid<coord,value>::size_type size_type;

// DEBUG -- Print an array
template <typename V>
void print_array(size_type n, size_type m, V v) {
  // row
  for (size_type i=0; i<n; ++i) {
    // col
    for (size_type j=0; j<m-1; ++j) {
      std::cout << v[i*m+j] << ",";
    }
    std::cout << v[i*m+n] << std::endl;
  }
  return;
}

int main() {
  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  int Nt = 34; // Total vertical points
  int Nbl= 22; // Number of boundary layer points
  int Mt = 34; // Total horizontal points

  std::vector<coord> X(Nt*Mt);  // Generate grid points for
  make_flat_plate(Nt,Nbl,Mt,X); // boundary layer simulation

  /* --- TEST GRID HANDLING --- */
  value u_inf;  u_inf[0] = 1; u_inf[1] = 1;   // Uniform flow condition
  value u_wall; u_wall[0] = 0; u_wall[1] = 0; // Wall (dirichlet) condition
  std::vector<value> v( (Nt-2)*(Mt-2), u_inf );
  // Set dirichlet condition on bottom
  for (int j = 0; j<Mt-2; ++j) {
    v[(Nt-3)*(Mt-2)+j] = u_wall;
  }

  StructuredGrid<coord,value> grid(Nt,Mt,v,X);
  StructuredGrid<coord,value>::Access val = grid.access();

  value t = val(2,2);
  std::cout << "(" << t[0] << "," << t[1] << ")" << std::endl;

  // DEBUG -- write values to file
  grid.write_grid("grid.dat");      // grid points
  grid.write_values("values.dat");  // cell values

  // grid.printv();

  return 0;
}