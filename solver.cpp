#include "structured_grid.hpp"
#include "grid_mapping.hpp"
#include "eflux.hpp"

#include <iostream>
#include <array>
#include <string>
#include <fstream>

using scalar = float;             // Specify solver precision

using coord = std::array<scalar,2>; // Grid coordinates
using value = std::array<scalar,4>; // State vector values
using flag  = std::array<bool,4>; // State vector values

typedef StructuredGrid<scalar,coord,value,flag>::size_type size_type;

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
  // Physical parameters
  scalar u_inf   = 1;
  scalar rho_inf = 1;
  scalar v_inf   = 0;
  scalar e_inf   = 1;
  // Time integration parameters
  // scalar     h = 1e-3;       // time step
  // size_type iter_max = 1e3; // max iterations
  // size_type n = 0;          // current iterations
  // Discretization parameters
  int Nt = 36; // Total vertical cells
  int Nbl= 23; // Number of boundary layer cells
  int Mt = 38; // Total horizontal cells

  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  std::vector<coord> X((Nt-1)*(Mt-1));  // Generate grid points for
  make_flat_plate((Nt-1),Nbl,(Mt-1),X); // boundary layer simulation

  /* --- SET UP GRID --- */
  // Flow conditions
  value U_inf = {rho_inf,rho_inf*u_inf,rho_inf*v_inf,rho_inf*e_inf}; // Inlet
  // value U_wall = {rho_inf,0,0,rho_inf*e_inf};   // Wall state
  value U_wall = {2*rho_inf,0,0,2*rho_inf*e_inf};   // Wall state
  // Boundary conditions
  flag B_wall = {1,0,0,1}; // Dirichlet in momentum, neumann in density and energy
  flag B_dir  = {0,0,0,0}; // Full dirichlet condition
  flag B_neu  = {1,1,1,1}; // Full neumann condition
  // Reserve space for cell values
  std::vector<value> V( Nt*Mt, U_inf );
  // Set dirichlet value on bottom
  for (int j = 1; j<Mt; ++j) {
    V[(Nt-1)*Mt+j] = U_wall;
  }

  // Specify boundary condition flag vectors
  std::vector<flag> left_b(Nt,B_dir);   // Dirichlet inlet
  std::vector<flag> right_b(Nt,B_neu);  // Neumann outlet
  right_b[Nt-1] = B_wall;               // plate extends to right end of domain
  std::vector<flag> top_b(Mt-2,B_neu);  // Neumann top
  std::vector<flag> bot_b(Mt-2,B_wall); // Wall bottom

  // Define grid
  StructuredGrid<scalar,coord,value,flag> grid(Nt,Mt,V,X,
                                        left_b,right_b,top_b,bot_b);

  // DEBUG -- Check bc handling
  auto val = grid.access();
  val(35,37).print(); std::cout<<std::endl;

  /* --- RESERVE SPACE FOR RK4 --- */
  value zeros = {0,0,0,0};
  std::vector<value> v0( (Nt-2)*(Mt-2) );

  /* --- RUN SOLVER --- */
  std::fill(v0.begin(),v0.end(),zeros);
  // DEBUG -- single step
  // eflux(grid.cell_begin(),grid.cell_end(),v0);

  /* --- FILE OUTPUT --- */
  grid.write_grid("solution.grid.dat");      // grid points
  grid.write_values("solution.val.dat");  // cell values

  return 0;
}