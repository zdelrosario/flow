#include "structured_grid.hpp"
#include "grid_mapping.hpp"
#include "eflux.hpp"
#include "nsflux.hpp"
#include "rk4_viscous.hpp"

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
// DEBUG -- Print a vector
template <typename V>
void print_vector(V v) {
  std::cout<<"(";
  for (auto it=v.begin(); it+1!=v.end(); ++it) {
    std::cout<<*it<<",";
  }
  std::cout<<v[v.size()-1]<<")"<<std::endl;
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
  scalar     h = 1e-8;       // time step
  // size_type iter_max = 1e3; // max iterations
  // size_type n = 0;          // current iterations
  // Discretization parameters
  int Nt = 36; // Total vertical cells
  int Nbl= 23; // Number of boundary layer cells
  int Mt = 40; // Total horizontal cells
  int buf = 4; // Freestream buffer cells

  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  std::vector<coord> X((Nt-1)*(Mt-1));    // Generate grid points for
  make_flat_plate((Nt-1),Nbl,(Mt-1),X,buf); // boundary layer simulation

  /* --- SET UP GRID --- */
  // Flow conditions
  value U_inf = {rho_inf,rho_inf*u_inf,rho_inf*v_inf,rho_inf*e_inf}; // Inlet
  value U_wall = {rho_inf,0,0,rho_inf*e_inf};   // Wall state
  // Boundary conditions
  flag B_wall = {1,2,2,1}; // Mirror momentum, neumann in density and energy
  flag B_in   = {3,3,3,3}; // Inlet condition
  flag B_out  = {4,4,4,4}; // Outlet condition
  flag B_mir  = {1,1,2,1}; // Vertical mirror
  // Reserve space for cell values
  std::vector<value> V( Nt*Mt, U_inf );
  // Set dirichlet value on bottom
  // for (int j = 1; j<Mt; ++j) { // DEPRECIATED
  //   V[(Nt-1)*Mt+j] = U_inf;
  // }

  // Specify boundary condition flag vectors
  std::vector<flag> left_b(Nt,B_in);    // Subsonic inlet
  std::vector<flag> right_b(Nt,B_out);  // Subsonic outlet
  std::vector<flag> top_b(Mt-2,B_out);  // Subsonic outlet
  std::vector<flag> bot_b(Mt-2,B_wall); // Wall bottom
  for (int i=0; i<buf; ++i) {      // Freestream mirror
    bot_b[i] = B_mir;
    bot_b[Mt-3-i] = B_mir;
  }
  right_b[Nt-1] = B_mir;

  // Define grid
  GridType grid(Nt,Mt,V,X,left_b,right_b,top_b,bot_b);
  auto val = grid.access();

  // DEBUG -- Check bc flags
  // for (auto it=bot_b.begin(); it!=bot_b.end(); ++it) {
  //   print_vector(*it);
  // }

  // DEBUG -- Check bc handling
  // std::cout<<"(1,1):"<<std::endl;
  // val(1,1).print(); std::cout<<std::endl;
  // std::cout<<"(1,0):"<<std::endl;
  // val(1,0).print(); std::cout<<std::endl;
  // std::cout<<"(0,1):"<<std::endl;
  // val(0,1).print(); std::cout<<std::endl;

  // std::cout<<"(34,1):"<<std::endl;
  // val(34,1).print(); std::cout<<std::endl;
  // std::cout<<"(34,0):"<<std::endl;
  // val(34,0).print(); std::cout<<std::endl;
  // std::cout<<"(35,1):"<<std::endl;
  // val(35,1).print(); std::cout<<std::endl;
  // std::cout<<"(35,2):"<<std::endl;
  // val(35,2).print(); std::cout<<std::endl;
  // std::cout<<"(35,3):"<<std::endl;
  // val(35,3).print(); std::cout<<std::endl;
  // std::cout<<"(35,4):"<<std::endl;
  // val(35,4).print(); std::cout<<std::endl;
  // std::cout<<"(35,5):"<<std::endl;
  // val(35,5).print(); std::cout<<std::endl;
  // std::cout<<"(35,6):"<<std::endl;
  // val(35,6).print(); std::cout<<std::endl;

  // DEBUG -- Check access of RK stages
  // val(1,1,4).print(); std::cout<<std::endl;

  /* --- RUN SOLVER --- */
  grid.fill_stages({0,0,0,0}); // Zero out the RK stages

  // DEBUG -- single iteration of RK

    // DEBUG -- Check the results of Euler flux
  size_type i=34,j=1; // Bottom left
  // size_type i=1,j=1; // Top left
  // size_type i=2,j=2; // Interior
  // size_type i=34,j=1+buf; // Plate front

  std::cout << "i="<<i<<",j="<<j<<std::endl;
  std::cout << "end of RK step" << std::endl;
  std::cout<<"y(t=0)="; val(i,j).print(); std::cout<<std::endl;   // y(t=0)
  rk4(h,val);
  std::cout<<"k1=";     val(i,j,0).print(); std::cout<<std::endl; // k1
  std::cout<<"k2=";     val(i,j,1).print(); std::cout<<std::endl; // k2
  std::cout<<"k3=";     val(i,j,2).print(); std::cout<<std::endl; // k3
  std::cout<<"k4=";     val(i,j,3).print(); std::cout<<std::endl; // k4
  std::cout<<"y(t=h)="; val(i,j).print(); std::cout<<std::endl; // y(t=h)

  /* --- FILE OUTPUT --- */
  grid.write_grid("solution.grid.dat");   // grid points
  grid.write_values("solution.val.dat");  // cell values

  return 0;
}