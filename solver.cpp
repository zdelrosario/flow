#include "structured_grid.hpp"  // Grid/cell handling
#include "grid_mapping.hpp"     // Grid generation
#include "eflux.hpp"            // Euler fluxes
#include "nsflux.hpp"           // Viscous fluxes
#include "rk4_viscous.hpp"      // RK4 time integration
#include "get_time.hpp"         // Code timing

#include <iostream>   // std::cout, std::endl
#include <array>      // std::array
#include <valarray>   // std::valarray, std::begin, std::end
#include <algorithm>  // std::transform

using scalar = double;             // Specify solver precision

using coord = std::array<scalar,2>;     // Grid coordinates
using value = std::valarray<scalar>;    // State vector
using flag  = std::array<unsigned,4>;   // Boundary condition flag

typedef StructuredGrid<scalar,coord,value,flag> GridType;
typedef GridType::size_type size_type;

int main() {
  /* --- SOLVER PARAMETERS --- */
  // Physical parameters
  scalar u_inf   = 68.93;   // Horizontal velocity
  scalar rho_inf = 1.1462;  // Density
  scalar v_inf   = 0;       // Vertical velocity
  scalar e_inf   = 298537;  // Internal energy
  // Time integration parameters
  // scalar h = 1e-8;          // fixed timestep
  size_type iter_max = 2e4; // max iterations
  size_type n = 0;          // current iterations
  size_type stride = 1e2;   // iteration stride for console printback
  scalar res_min = 5e-3;    // residual convergence tolerance
  // Discretization parameters
  int Nt = 36; // Total vertical cells
  int Nbl= 23; // Number of boundary layer cells
  int Mt = 40; // Total horizontal cells
  int buf = 4; // Freestream buffer cells

  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  std::vector<coord> X((Nt-1)*(Mt-1));      // Generate grid points for
  make_flat_plate((Nt-1),Nbl,(Mt-1),X,buf); // boundary layer simulation

  /* --- SET UP GRID --- */
  // Flow conditions
  value U_inf = {rho_inf,rho_inf*u_inf,rho_inf*v_inf,rho_inf*e_inf}; // Inlet
  // Boundary conditions
  flag B_wall = {{1,2,2,1}}; // Mirror momentum, neumann in density and energy
  flag B_in   = {{3,3,3,3}}; // Inlet condition
  flag B_out  = {{4,4,4,4}}; // Outlet condition
  flag B_mir  = {{1,1,2,1}}; // Vertical mirror
  // Reserve space for cell values
  std::vector<value> V( Nt*Mt, U_inf );

  // Specify boundary condition flag vectors
  std::vector<flag> left_b(Nt,B_in);    // Subsonic inlet
  std::vector<flag> right_b(Nt,B_out);  // Subsonic outlet
  std::vector<flag> top_b(Mt-2,B_out);  // Subsonic outlet
  std::vector<flag> bot_b(Mt-2,B_wall); // Wall bottom
  // Freestream mirror
  for (int i=0; i<buf; ++i) {
    bot_b[i] = B_mir;         // leading mirror
    // bot_b[Mt-3-i] = B_mir;    // trailing mirror
  }
  right_b[Nt-1] = B_wall;

  // Define grid
  GridType grid(Nt,Mt,V,X,left_b,right_b,top_b,bot_b);
  auto val = grid.access();

  /* --- RUN SOLVER --- */
  scalar res = 1e10;
  uint64 T_0 = GetTimeMs64();
  uint64 T;
  double dT;
  while ((n<iter_max) and (res>res_min)) {
    // Zero out the RK stages
    grid.fill_stages({0,0,0,0});
    // Take RK4 time step
    res = rk4_local(val);
    // res = euler_local(val);
    // Compute current time
    T = GetTimeMs64();
    dT = double(T-T_0)/1e3/60.; // minutes
    // Print back residual, execution time
    if (n % stride == 0) {
      std::cout << "n=" << n << ", res=" << res;
      std::cout << ", dT=" << dT << "min" << std::endl;
    }
    // Iterate the counter
    ++n;
  }

  // Print final iteration
  std::cout << "n=" << n << ", res=" << res;
  std::cout << ", dT=" << dT << "min" << std::endl;

  // Print status messages
  if (res <= res_min) {
    std::cout << "Convergence tolerance reached!" << std::endl;
  }
  else
    std::cout << "Convergence not reached..." << std::endl;
  if (n >= iter_max) {
    std::cout << "Iteration limit reached." << std::endl;
  }
  else
    std::cout << "Iteration limit not reached" << std::endl;
  
  /* --- FILE OUTPUT --- */
  grid.write_grid("solution.grid.dat");   // grid points
  grid.write_values("solution.val.dat");  // cell values

  return 0;
}