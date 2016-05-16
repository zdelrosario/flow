#include "structured_grid.hpp"
#include "grid_mapping.hpp"
#include "eflux.hpp"
#include "nsflux.hpp"

#include <iostream>   // std::cout, std::endl
#include <array>      // std::array
#include <valarray>   // std::valarray, std::begin, std::end
#include <algorithm>  // std::transform

using scalar = float;             // Specify solver precision

using coord = std::array<scalar,2>; // Grid coordinates
using value = std::valarray<scalar>;// State vector
using flag  = std::array<unsigned,4>;   // Boundary condition flag

typedef StructuredGrid<scalar,coord,value,flag> GridType;
typedef GridType::size_type size_type;

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
  scalar u_inf   = 68.93;
  scalar rho_inf = 1.1462;
  scalar v_inf   = 0;
  scalar e_inf   = 298537;
  // Time integration parameters
  scalar     h = 1e-6;       // time step
  // size_type iter_max = 1e3; // max iterations
  // size_type n = 0;          // current iterations
  // Discretization parameters
  int Nt = 36; // Total vertical cells
  int Nbl= 23; // Number of boundary layer cells
  int Mt = 40; // Total horizontal cells

  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  std::vector<coord> X((Nt-1)*(Mt-1));    // Generate grid points for
  make_flat_plate((Nt-1),Nbl,(Mt-1),X,4); // boundary layer simulation

  /* --- SET UP GRID --- */
  // Flow conditions
  value U_inf = {rho_inf,rho_inf*u_inf,rho_inf*v_inf,rho_inf*e_inf}; // Inlet
  value U_wall = {rho_inf,0,0,rho_inf*e_inf};   // Wall state
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
  GridType grid(Nt,Mt,V,X,left_b,right_b,top_b,bot_b);
  auto val = grid.access();

  // DEBUG -- Check bc handling
  // val(35,37).print(); std::cout<<std::endl;

  // DEBUG -- Check access of RK stages
  // val(1,1,4).print(); std::cout<<std::endl;

  /* --- RUN SOLVER --- */
  grid.fill_stages({0,0,0,0}); // Zero out the RK stages

  // DEBUG -- single iteration of RK
  size_type i=34,j=1; // Bottom left
  // size_type i=1,j=1; // Top left
  value y = {0,0,0,0};   value f = {0,0,0,0};
  value k1 = {0,0,0,0};  value k2 = {0,0,0,0};
  value k3 = {0,0,0,0};  value k4 = {0,0,0,0};
  // Stage 1
  eflux(grid.cell_begin(),grid.cell_end(),grid.cell_begin(0));
  nsflux(grid.cell_begin(),grid.cell_end(),grid.cell_begin(0));
  auto flux_it = grid.cell_begin(0);
  auto  out_it = grid.cell_begin(4);
  for (auto it=grid.cell_begin(); it!=grid.cell_end(); ++it) {
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + scalar(0.5)*h*f;
    // Iterate
    ++flux_it;
    ++out_it;
  }
std::cout << "Result of stage 1" << std::endl;
val(i,j,4).print(); std::cout<<std::endl;
  // Stage 2
  eflux(grid.cell_begin(4),grid.cell_end(4),grid.cell_begin(1));
  nsflux(grid.cell_begin(4),grid.cell_end(4),grid.cell_begin(1));
  flux_it = grid.cell_begin(1);
   out_it = grid.cell_begin(4);
  for (auto it=grid.cell_begin(); it!=grid.cell_end(); ++it) {
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + scalar(0.5)*h*f;
    // Iterate
    ++flux_it;
    ++out_it;
  }
std::cout << "Result of stage 2" << std::endl;
val(i,j,4).print(); std::cout<<std::endl;
  // Stage 3
  eflux(grid.cell_begin(4),grid.cell_end(4),grid.cell_begin(2));
  nsflux(grid.cell_begin(4),grid.cell_end(4),grid.cell_begin(2));
  flux_it = grid.cell_begin(2);
   out_it = grid.cell_begin(4);
  for (auto it=grid.cell_begin(); it!=grid.cell_end(); ++it) {
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + h*f;
    // Iterate
    ++flux_it;
    ++out_it;
  }
std::cout << "Result of stage 3" << std::endl;
val(i,j,4).print(); std::cout<<std::endl;
  // Stage 4
  eflux(grid.cell_begin(4),grid.cell_end(4),grid.cell_begin(3));
  nsflux(grid.cell_begin(4),grid.cell_end(4),grid.cell_begin(3));
  auto k1_it = grid.cell_begin(0);
  auto k2_it = grid.cell_begin(1);
  auto k3_it = grid.cell_begin(2);
  auto k4_it = grid.cell_begin(3);
  out_it = grid.cell_begin(4); // DEBUG -- should write to (-1)
  for (auto it=grid.cell_begin(); it!=grid.cell_end(); ++it) {
    // Apply RK stage
    y  = (*it).value();
    k1 = (*k1_it).value();    k2 = (*k2_it).value();
    k3 = (*k3_it).value();    k4 = (*k4_it).value();
    (*out_it) = y + h/scalar(6.0)*(k1+scalar(2)*k2+scalar(2)*k3+k4);
    // Iterate
    ++flux_it;
    ++out_it;
  }

  // DEBUG -- Check the results of Euler flux
  std::cout << "end of RK step" << std::endl;
  std::cout<<"y(t=0)="; val(i,j).print(); std::cout<<std::endl;   // y(t=0)
  std::cout<<"k1=";     val(i,j,0).print(); std::cout<<std::endl; // k1
  std::cout<<"k2=";     val(i,j,1).print(); std::cout<<std::endl; // k2
  std::cout<<"k3=";     val(i,j,2).print(); std::cout<<std::endl; // k3
  std::cout<<"k4=";     val(i,j,3).print(); std::cout<<std::endl; // k4
  std::cout<<"y(t=h)="; val(i,j,4).print(); std::cout<<std::endl; // y(t=h)

  /* --- FILE OUTPUT --- */
  grid.write_grid("solution.grid.dat");      // grid points
  grid.write_values("solution.val.dat");  // cell values

  return 0;
}