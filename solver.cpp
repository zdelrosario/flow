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
#include <cmath>      // sin, cos

using scalar = double;             // Specify solver precision

using coord = std::array<scalar,2>;     // Grid coordinates
using value = std::valarray<scalar>;    // State vector
using flag  = std::array<unsigned,4>;   // Boundary condition flag

typedef StructuredGrid<scalar,coord,value,flag> GridType;
typedef GridType::size_type size_type;

// 
// SUPERSONIC BUMP
// 
void bump() {
  /* --- SOLVER PARAMETERS --- */
  scalar rho_inf = 1.1462;  // Density
  scalar u_inf   = 7.72e2;   // Horizontal velocity
  scalar p_e     = 1.2199e5; // Exit pressure
  scalar e_inf   = 5.6406e5;  // Internal energy
  scalar t       = 0.12;   // Thickness
  // Solver parameters
  bool restart_flag = 0;    // Load the restart file?
  size_type iter_max = 1e6; // max iterations
  size_type n = 0;          // current iterations
  size_type stride = 1e3;   // iteration stride for console printback
  size_type restart = 1e4;  // iteration stride for restart file
  short  res_type = -1;     // Density
  scalar res_min = 1e-3;    // residual convergence tolerance
  // Discretization parameters
  int Nt = 40; // Total vertical cells
  int Nbl= 20; // Number of boundary layer cells
  int Mt = 40; // Total horizontal cells
  int buf = 5; // Freestream buffer cells

  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  std::vector<coord> X((Nt-1)*(Mt-1));      // Generate grid points for
  // make_flat_plate((Nt-1),Nbl,(Mt-1),X,buf); // boundary layer simulation
  make_bump(Nt-1,Nbl,Mt-1,X,buf,t); // Bump

  /* --- SET UP GRID --- */
  // Flow conditions
  value U_inf = {rho_inf,rho_inf*u_inf,0,rho_inf*e_inf}; // Inlet
  // Boundary conditions
  flag B_wall = {{1,2,2,1}}; // Mirror momentum, neumann in density and energy
  flag B_in   = {{3,3,3,3}}; // Inlet condition
  flag B_pres = {{4,4,4,4}}; // Outlet condition
  // flag B_pres = {{5,5,5,5}}; // Pressure set
  flag B_mir  = {{1,1,2,1}}; // Vertical mirror
  // Reserve space for cell values
  std::vector<value> V( Nt*Mt, U_inf );

  // Specify boundary condition flag vectors
  std::vector<flag> left_b(Nt,B_in);    // Subsonic inlet
  // std::vector<flag> right_b(Nt,B_out);  // Subsonic outlet
  std::vector<flag> right_b(Nt,B_pres);  // Pressure set
  // std::vector<flag> top_b(Mt-2,B_out);  // Subsonic outlet
  std::vector<flag> top_b(Mt-2,B_pres);  // Pressure set
  std::vector<flag> bot_b(Mt-2,B_wall); // Wall bottom
  // Freestream mirror
  for (int i=0; i<buf; ++i) {
    bot_b[i] = B_mir;         // leading mirror
    bot_b[Mt-3-i] = B_mir;    // trailing mirror
  }
  right_b[Nt-1] = B_mir;

  /* --- LOAD RESTART IF APPLICABLE --- */
  if (restart_flag) {
    std::cout << "Reading restart file..." << std::endl;
    std::vector<value> V_temp;
    readin_val("restart.val.dat", V_temp);
    // Overwrite interior values
    for (size_type i=0; i<size_type(Nt-2); ++i){
      for (size_type j=0; j<size_type(Mt-2); ++j) {
        V[(i+1)*Mt+(j+1)] = V_temp[i*(Mt-2)+j];
      }
    }
  }
  else {
    std::cout << "Uniform initial conditions assumed..." << std::endl;
  }

  // Define grid
  GridType grid(Nt,Mt,V,X,left_b,right_b,top_b,bot_b,p_e);
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
    res = rk4_local(val,res_type);
    // res = euler_local(val);
    // Compute current time
    T = GetTimeMs64();
    dT = double(T-T_0)/1e3/60.; // minutes
    // Print back residual, execution time
    if (n % stride == 0) {
      std::cout << "n=" << n << ", res=" << res;
      std::cout << ", dT=" << dT << "min" << std::endl;
    }
    // Write out restart file
    if ((n % restart == 0) and (n>0)) {
      grid.write_grid("restart.grid.dat");   // grid points
      grid.write_values("restart.val.dat");  // cell values
    }
    // Iterate the counter
    ++n;
  }
  // Print final iteration
  std::cout << "n=" << n << ", res=" << res;
  std::cout << ", dT=" << dT << "min" << std::endl;
  // Print convergence message
  if (res <= res_min) {
    std::cout << "Convergence tolerance reached!" << std::endl;
  }
  else
    std::cout << "Convergence not reached..." << std::endl;
  // Print iteration message
  if (n >= iter_max) {
    std::cout << "Iteration limit reached." << std::endl;
  }
  else
    std::cout << "Iteration limit not reached" << std::endl;
  /* --- FILE OUTPUT --- */
  grid.write_grid("solution.grid.dat");   // grid points
  grid.write_values("solution.val.dat");  // cell values
  return;
}

// 
// FLAT PLATE BOUNDARY LAYER
// 
void flat_plate() {
  /* --- SOLVER PARAMETERS --- */
  // Physical parameters
  scalar u_inf   = 68.93;   // Horizontal velocity
  scalar rho_inf = 1.1462;  // Density
  scalar v_inf   = 0;       // Vertical velocity
  scalar e_inf   = 298537;  // Internal energy
  scalar p_e     = 9.725e4; // Exit pressure
  // Solver parameters
  bool restart_flag = 0;    // Load the restart file?
  size_type iter_max = 1e4; // max iterations
  size_type n = 0;          // current iterations
  size_type stride = 1e3;   // iteration stride for console printback
  size_type restart = 1e4;  // iteration stride for restart file
  short  res_type = -1;     // Max over all fluxes
  scalar res_min = 1e-2;    // residual convergence tolerance
  // Discretization parameters
  int Nt = 40; // Total vertical cells
  int Nbl= 30; // Number of boundary layer cells
  int Mt = 40; // Total horizontal cells
  int buf = 5; // Freestream buffer cells

  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  std::vector<coord> X((Nt-1)*(Mt-1));      // Generate grid points for
  make_flat_plate((Nt-1),Nbl,(Mt-1),X,buf); // boundary layer simulation

  /* --- SET UP GRID --- */
  // Flow conditions
  value U_inf = {rho_inf,rho_inf*u_inf,rho_inf*v_inf,rho_inf*e_inf}; // Inlet
  // Boundary conditions
  flag B_wall = {{1,2,2,1}}; // Mirror momentum, neumann in density and energy
  flag B_in   = {{3,3,3,3}}; // Inlet condition
  flag B_pres = {{4,4,4,4}}; // Outlet condition
  // flag B_pres = {{5,5,5,5}}; // Pressure set
  flag B_mir  = {{1,1,2,1}}; // Vertical mirror
  // Reserve space for cell values
  std::vector<value> V( Nt*Mt, U_inf );

  // Specify boundary condition flag vectors
  std::vector<flag> left_b(Nt,B_in);    // Subsonic inlet
  // std::vector<flag> right_b(Nt,B_out);  // Subsonic outlet
  std::vector<flag> right_b(Nt,B_pres);  // Pressure set
  // std::vector<flag> top_b(Mt-2,B_out);  // Subsonic outlet
  std::vector<flag> top_b(Mt-2,B_pres);  // Pressure set
  std::vector<flag> bot_b(Mt-2,B_wall); // Wall bottom
  // Freestream mirror
  for (int i=0; i<buf; ++i) {
    bot_b[i] = B_mir;         // leading mirror
    // bot_b[Mt-3-i] = B_mir;    // trailing mirror
  }
  right_b[Nt-1] = B_mir;

  /* --- LOAD RESTART IF APPLICABLE --- */
  if (restart_flag) {
    std::cout << "Reading restart file..." << std::endl;
    std::vector<value> V_temp;
    readin_val("restart.val.dat", V_temp);
    // Overwrite interior values
    for (size_type i=0; i<size_type(Nt-2); ++i){
      for (size_type j=0; j<size_type(Mt-2); ++j) {
        V[(i+1)*Mt+(j+1)] = V_temp[i*(Mt-2)+j];
      }
    }
  }
  else {
    std::cout << "Uniform initial conditions assumed..." << std::endl;
  }

  // Define grid
  GridType grid(Nt,Mt,V,X,left_b,right_b,top_b,bot_b,p_e);
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
    res = rk4_local(val,res_type);
    // res = euler_local(val);
    // Compute current time
    T = GetTimeMs64();
    dT = double(T-T_0)/1e3/60.; // minutes
    // Print back residual, execution time
    if (n % stride == 0) {
      std::cout << "n=" << n << ", res=" << res;
      std::cout << ", dT=" << dT << "min" << std::endl;
    }
    // Write out restart file
    if ((n % restart == 0) and (n>0)) {
      grid.write_grid("restart.grid.dat");   // grid points
      grid.write_values("restart.val.dat");  // cell values
    }
    // Iterate the counter
    ++n;
  }
  // Print final iteration
  std::cout << "n=" << n << ", res=" << res;
  std::cout << ", dT=" << dT << "min" << std::endl;
  // Print convergence message
  if (res <= res_min) {
    std::cout << "Convergence tolerance reached!" << std::endl;
  }
  else
    std::cout << "Convergence not reached..." << std::endl;
  // Print iteration message
  if (n >= iter_max) {
    std::cout << "Iteration limit reached." << std::endl;
  }
  else
    std::cout << "Iteration limit not reached" << std::endl;
  /* --- FILE OUTPUT --- */
  grid.write_grid("solution.grid.dat");   // grid points
  grid.write_values("solution.val.dat");  // cell values
  return;
}

// 
// OBLIQUE SHOCK
// 
void oblique_shock() {
  /* --- SOLVER PARAMETERS --- */
  // Physical parameters
  scalar theta   = 15.0 * 3.14/180.0; // Turn angle (radians)
  scalar c_sound = 347.19;  // Speed of sound
  scalar mach    = 2;       // Desired mach number
  scalar u_inf   = mach*c_sound*cos(-theta);// Horizontal velocity
  scalar v_inf   = mach*c_sound*sin(-theta);// Vertical velocity
  scalar rho_inf = 1.1462;  // Density
  scalar e_inf   = 453193;  // Internal energy
  scalar p_e     = 9.725e4; // Exit pressure
  // Solver parameters
  size_type iter_max = 1e6; // max iterations
  size_type n = 0;          // current iterations
  size_type stride  = 1e3;  // iteration stride for console printback
  size_type restart = 1e4;  // iteration stride for restart file
  short  res_type = -1;     // Max over all fluxes
  scalar res_min = 1e-8;    // residual convergence tolerance
  // Discretization parameters
  int Nt = 50; // Total vertical cells
  int Mt = 50; // Total horizontal cells

  /* --- FLAT PLATE BOUNDARY LAYER GRID --- */
  std::vector<coord> X((Nt-1)*(Mt-1));      // Generate grid points for
  make_rectangle((Nt-1), (Mt-1), X); // Simple rectangle domain

  /* --- SET UP GRID --- */
  // Flow conditions
  value U_inf = {rho_inf,rho_inf*u_inf,rho_inf*v_inf,rho_inf*e_inf}; // Inlet
  // Boundary conditions
  flag B_in   = {{0,0,0,0}}; // Supersonic Inlet condition
  flag B_out  = {{1,1,1,1}}; // Supersonic Outlet condition
  flag B_wall = {{1,1,2,1}}; // Slip wall
  // Reserve space for cell values
  std::vector<value> V( Nt*Mt, U_inf );

  // Specify boundary condition flag vectors
  std::vector<flag> left_b(Nt,B_in);    // Supersonic inlet
  std::vector<flag> top_b(Mt-2,B_in);   // Supersonic inlet
  std::vector<flag> right_b(Nt,B_out);  // Supersonic outlet
  std::vector<flag> bot_b(Mt-2,B_wall); // Wall bottom
  right_b[Nt-1] = B_wall; // 

  // Define grid
  GridType grid(Nt,Mt,V,X,left_b,right_b,top_b,bot_b,p_e);
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
    res = rk4_local(val,res_type);
    // res = euler_local(val);
    // Compute current time
    T = GetTimeMs64();
    dT = double(T-T_0)/1e3/60.; // minutes
    // Print back residual, execution time
    if (n % stride == 0) {
      std::cout << "n=" << n << ", res=" << res;
      std::cout << ", dT=" << dT << "min" << std::endl;
    }
    // Write out restart file
    if (n % restart == 0) {
      grid.write_grid("restart.grid.dat");   // grid points
      grid.write_values("restart.val.dat");  // cell values
    }
    // Iterate the counter
    ++n;
  }
  // Print final iteration
  std::cout << "n=" << n << ", res=" << res;
  std::cout << ", dT=" << dT << "min" << std::endl;
  // Print convergence message
  if (res <= res_min) {
    std::cout << "Convergence tolerance reached!" << std::endl;
  }
  else
    std::cout << "Convergence not reached..." << std::endl;
  // Print iteration message
  if (n >= iter_max) {
    std::cout << "Iteration limit reached." << std::endl;
  }
  else
    std::cout << "Iteration limit not reached" << std::endl;
  /* --- FILE OUTPUT --- */
  grid.write_grid("solution.grid.dat");   // grid points
  grid.write_values("solution.val.dat");  // cell values
  return;
}

int main() {
  bump();
  // flat_plate();
  // oblique_shock();
  return 0;
}