#include "structured_grid.hpp"
#include "grid_mapping.hpp"
#include "eflux.hpp"

#include <iostream>
#include <array>
#include <string>
#include <fstream>

using coord = std::array<float,2>; // Grid coordinates
using value = std::array<float,4>; // State vector values

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
  /* --- SOLVER PARAMETERS --- */
  // float     h = 1e-3;       // time step
  // size_type iter_max = 1e3; // max iterations
  // size_type n = 0;          // current iterations

  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  int Nt = 34; // Total vertical cells
  int Nbl= 22; // Number of boundary layer cells
  int Mt = 34; // Total horizontal cells

  std::vector<coord> X((Nt-1)*(Mt-1));  // Generate grid points for
  make_flat_plate((Nt-1),Nbl,(Mt-1),X); // boundary layer simulation

  /* --- SET UP GRID --- */
  value u_inf = {1,1,0,1};    // Uniform flow condition
  value u_wall = {1,0,0,1};   // Wall (dirichlet) condition

  std::vector<value> V( (Nt-2)*(Mt-2), u_inf );
  // Set dirichlet condition on bottom
  for (int j = 0; j<Mt-2; ++j) {
    V[(Nt-3)*(Mt-2)+j] = u_wall;
  }
  StructuredGrid<coord,value> grid(Nt,Mt,V,X);

  /* --- RESERVE SPACE FOR RK4 --- */


  /* --- RUN SOLVER --- */
  // DEBUG -- single step


  /* --- FILE OUTPUT --- */
  grid.write_grid("grid.dat");      // grid points
  grid.write_values("values.dat");  // cell values

  return 0;
}