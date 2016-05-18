#include "eflux.hpp"
#include "nsflux.hpp"

template <typename S, typename A>
void rk4(S h, A access) {
  using value  = typename A::GridValue;
  using scalar = typename A::GridScalar;
  value y = {0,0,0,0};   value f = {0,0,0,0};
  value k1 = {0,0,0,0};  value k2 = {0,0,0,0};
  value k3 = {0,0,0,0};  value k4 = {0,0,0,0};
  
  // Stage 1
  eflux(access.cell_begin(),access.cell_end(),access.cell_begin(0));
  nsflux(access.cell_begin(),access.cell_end(),access.cell_begin(0));
  auto flux_it = access.cell_begin(0);
  auto  out_it = access.cell_begin(4);
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + scalar(0.5)*h*f;
    // Iterate
    ++flux_it;
    ++out_it;
  }

  // Stage 2
  eflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(1));
  nsflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(1));
  flux_it = access.cell_begin(1);
   out_it = access.cell_begin(4);
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + scalar(0.5)*h*f;
    // Iterate
    ++flux_it;
    ++out_it;
  }

  // Stage 3
  eflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(2));
  nsflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(2));
  flux_it = access.cell_begin(2);
   out_it = access.cell_begin(4);
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Apply RK stage
    y = (*it).value();
    f = (*flux_it).value();
    (*out_it) = y + h*f;
    // Iterate
    ++flux_it;
    ++out_it;
  }

  // Stage 4
  eflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(3));
  nsflux(access.cell_begin(4),access.cell_end(4),access.cell_begin(3));

  // Add all stages
  auto k1_it = access.cell_begin(0);
  auto k2_it = access.cell_begin(1);
  auto k3_it = access.cell_begin(2);
  auto k4_it = access.cell_begin(3);
   out_it = access.cell_begin();
  for (auto it=access.cell_begin(); it!=access.cell_end(); ++it) {
    // Apply RK stage
    y  = (*it).value();
    k1 = (*k1_it).value();    k2 = (*k2_it).value();
    k3 = (*k3_it).value();    k4 = (*k4_it).value();
    (*out_it) = y + h/scalar(6.0)*(k1+scalar(2)*k2+scalar(2)*k3+k4);
    // Iterate
    ++flux_it;
    ++out_it;
    ++k1_it;
    ++k2_it;
    ++k3_it;
    ++k4_it;
  }
}