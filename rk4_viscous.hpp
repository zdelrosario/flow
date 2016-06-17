#include <vector>
#include <cmath>
#include <valarray>

#include "eflux.hpp"
#include "nsflux.hpp"
#include "gas_dynamics.hpp"

// Jameson timestep estimate
/** Local time step estimate
 *
 * @param c Cell
 * 
 * @return h Local time step approximately
 *           corresponding to CFL=1
 */
template <typename C>
typename C::CellScalar timestep(C c) {
  typename C::CellScalar a,b;
  // TODO -- project horizontal and vertical velocities
  // for curved cells; currently assuming rectangular
  // elements
  a = (std::abs(uf(c.value()))+cf(c.value())) * c.dx();
  b = (std::abs(vf(c.value()))+cf(c.value())) * c.dy();
  return c.dx()*c.dy() / (a+b);
}

/** Take an RK4 timestep with eflux and nsflux
 *  Uses a local CFL estimate for timestep
 * 
 * @param access Grid access handle
 * @param res_type Residual type; -1=all entries, 0=density
 * 
 * @post Grid values updated based on fluxes
 *       cell_begin(-1) contains y(t+h)
 *       cell_begin(4)  contains y(t)
 * @return res Residual after computation
 */
template <typename A>
typename A::GridScalar rk4_local(A access, short res_type) {
  using value  = typename A::GridValue;
  using scalar = typename A::GridScalar;
  // short res_type = 0; // -1=all entries, 0=density

  std::vector<scalar> H;
  scalar res = -1;
  scalar h_min = 1, h_max = -1;
  value y = {0,0,0,0};   value f = {0,0,0,0};
  value k1 = {0,0,0,0};  value k2 = {0,0,0,0};
  value k3 = {0,0,0,0};  value k4 = {0,0,0,0};
  
  // Stage 1
  // Write k1 to cell_begin(0)
  eflux(access.cell_begin(),access.cell_end(),access.cell_begin(0));
  nsflux(access.cell_begin(),access.cell_end(),access.cell_begin(0));
  
  // Stage 2
  auto flux_it = access.cell_begin(0);
  auto  out_it = access.cell_begin(4);
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Local time step approximation
    H.push_back(timestep(*it));
    h_min = std::min(H.back(),h_min); // DEBUG
    h_max = std::max(H.back(),h_max);
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + scalar(0.5)*(H.back())*f;
    // Iterate
    ++flux_it;
    ++out_it;
  }
  // Write k2 to cell_begin(1)
  eflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(1));
  nsflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(1));

  // Stage 3
  flux_it = access.cell_begin(1);
   out_it = access.cell_begin(4);
  auto h_it = H.begin();
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + scalar(0.5)*(*h_it)*f;
    // Iterate
    ++flux_it;
    ++out_it;
    ++h_it;
  }
  // Write k3 to cell_begin(2)
  eflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(2));
  nsflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(2));
  
  // Stage 4
  flux_it = access.cell_begin(2);
   out_it = access.cell_begin(4);
     h_it = H.begin();
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + (*h_it)*f;
    // Iterate
    ++flux_it;
    ++out_it;
    ++h_it;
  }
  // Write k4 to cell_begin(3)
  eflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(3));
  nsflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(3));

  // Add all stages
  auto k1_it = access.cell_begin(0);
  auto k2_it = access.cell_begin(1);
  auto k3_it = access.cell_begin(2);
  auto k4_it = access.cell_begin(3);
  auto old_it = access.cell_begin(4);
   out_it = access.cell_begin();
     h_it = H.begin();
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Apply RK stage
    y  = (*it).value();
    k1 = (*k1_it).value();    k2 = (*k2_it).value();
    k3 = (*k3_it).value();    k4 = (*k4_it).value();
    f = (*h_it)/scalar(6.0)*(k1+scalar(2)*k2+scalar(2)*k3+k4);

     // update residual
    if (res_type == -1)
      res = std::max( (std::abs(f)).max(), res );
    else if (res_type == 0)
      res = std::max( std::abs(f[0]), res );
    else
      assert(false);

    (*old_it) = y; // store old value
    (*out_it) = y + f;
    // Iterate
    ++flux_it;
    ++out_it;
    ++old_it;
    ++k1_it;
    ++k2_it;
    ++k3_it;
    ++k4_it;
    ++h_it;
  }

  // DEBUG -- print min and max timesteps
  // std::cout << "h_min = " << h_min << std::endl;
  // std::cout << "h_max = " << h_max << std::endl;
  // return residual
  return res;
}

template <typename A>
typename A::GridScalar euler_local(A access, short res_type) {
  using value  = typename A::GridValue;
  using scalar = typename A::GridScalar;

  std::vector<scalar> H;
  scalar res = -1;
  scalar h_min = 1, h_max = -1;
  value y = {0,0,0,0};   value f = {0,0,0,0};
  value k1 = {0,0,0,0};  value k2 = {0,0,0,0};
  value k3 = {0,0,0,0};  value k4 = {0,0,0,0};
  
  // Stage 1
  // Write k1 to cell_begin(0)
  eflux(access.cell_begin(),access.cell_end(),access.cell_begin(0));
  nsflux(access.cell_begin(),access.cell_end(),access.cell_begin(0));
  
  // Stage 2
  auto flux_it = access.cell_begin(0);
  auto  out_it = access.cell_begin();
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Local time step approximation
    H.push_back(timestep(*it));
    h_min = std::min(H.back(),h_min); // DEBUG
    h_max = std::max(H.back(),h_max);
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + (H.back())*f;
    
    // update residual
    if (res_type == -1)
      res = std::max( (std::abs(f)).max(), res );
    else if (res_type == 0)
      res = std::max( std::abs(f[0]), res );
    else
      assert(false);
    
    // Iterate
    ++flux_it;
    ++out_it;
  }

  return res;
}