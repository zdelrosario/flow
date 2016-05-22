#ifndef NSFLUX // Include guard
#define NSFLUX
/** Viscous Flux Calculator
 */

#include <vector>     // data handling
#include <algorithm>  // std::copy(), std::transform()
#include <iostream>   // debug
#include <cmath>      // pow()
#include <valarray>   // std::valarray, std::begin, std::end
#include <functional>   // std::function

#include "gas_dynamics.hpp" // gas dynamics equations

using size_type = unsigned;

// 
// DUAL GRID DERIVATIVES
// 

// DEBUG -- estimated cell dimensions, need to use
//          grid mapping

/* --- HORIZONTAL DERIVATIVES --- */
/* --- du_dx --- */
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dx_r(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25/c.dx())*(
        // Upper right
          f(c.value(1,0))  - f(c.value(0,0)) + 
          f(c.value(1,-1)) - f(c.value(0,-1)) +
        // Lower right
          f(c.value(1,0)) - f(c.value(0,0)) + 
          f(c.value(1,1)) - f(c.value(0,1))
        );
}
template <typename Cell, class Fcn>
typename Cell::CellScalar df_dx_l(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25/c.dx())*(
        // Upper left
          f(c.value(0,0))  - f(c.value(-1,0)) + 
          f(c.value(0,-1)) - f(c.value(-1,-1)) + 
        // Lower left
          f(c.value(0,0)) - f(c.value(-1,0)) + 
          f(c.value(0,1)) - f(c.value(-1,1))
        );
}
template <typename Cell, class Fcn>
typename Cell::CellScalar df_dx_t(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25/c.dx())*(
        // Upper right
          f(c.value(1,-1))  - f(c.value(-1,-1)) + 
          f(c.value(1,0)) - f(c.value(-1,0))
        );
}
template <typename Cell, class Fcn>
typename Cell::CellScalar df_dx_b(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25/c.dx())*(
        // Upper left
          f(c.value(1,1)) - f(c.value(-1,1)) + 
          f(c.value(1,0)) - f(c.value(-1,0))
        );
}

/* --- VERTICAL DERIVATIVES --- */
/* --- du_dy --- */
template <typename Cell, class Fcn>
typename Cell::CellScalar df_dy_t(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25/c.dy())*(
        // Upper right
          f(c.value(1,-1)) - f(c.value(1,0)) + 
          f(c.value(0,-1)) - f(c.value(0,0)) +
        // Upper left
          f(c.value(0,-1))  - f(c.value(0,0)) + 
          f(c.value(-1,-1)) - f(c.value(-1,0))
        );
}
template <typename Cell, class Fcn>
typename Cell::CellScalar df_dy_b(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25/c.dy())*(
        // Bottom right
          f(c.value(1,0)) - f(c.value(1,1)) + 
          f(c.value(0,0)) - f(c.value(0,1)) +
        // Bottom left
          f(c.value(0,0))  - f(c.value(0,1)) + 
          f(c.value(-1,0)) - f(c.value(-1,1))
        );
}
template <typename Cell, class Fcn>
typename Cell::CellScalar df_dy_l(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25/c.dy())*(
          f(c.value(-1,-1)) - f(c.value(-1,1)) + 
          f(c.value(0,-1)) - f(c.value(0,1))
        );
}
template <typename Cell, class Fcn>
typename Cell::CellScalar df_dy_r(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25/c.dy())*(
          f(c.value(1,-1)) - f(c.value(1,1)) + 
          f(c.value(0,-1)) - f(c.value(0,1))
        );
}

// 
// VISCOUS STRESS TERMS
// 
template <typename Cell>
typename Cell::CellScalar tau_xx_r(Cell c) {
  using scalar = typename Cell::CellScalar;
  auto u = [](typename Cell::CellValue val) { return uf(val); };
  auto v = [](typename Cell::CellValue val) { return vf(val); };
  return scalar(muf(c.value())) * (
          scalar(0.5)*df_dx_r(c,u) - 
          scalar(1.5)*df_dy_r(c,v) );
}
template <typename Cell>
typename Cell::CellScalar tau_xx_l(Cell c) {
  using scalar = typename Cell::CellScalar;
  auto u = [](typename Cell::CellValue val) { return uf(val); };
  auto v = [](typename Cell::CellValue val) { return vf(val); };
  return scalar(muf(c.value())) * (
          scalar(0.5)*df_dx_l(c,u) - 
          scalar(1.5)*df_dy_l(c,v) );
}

template <typename Cell>
typename Cell::CellScalar tau_yy_t(Cell c) {
  using scalar = typename Cell::CellScalar;
  auto u = [](typename Cell::CellValue val) { return uf(val); };
  auto v = [](typename Cell::CellValue val) { return vf(val); };
  return scalar(muf(c.value())) * (
          scalar(0.5)*df_dy_t(c,v) - 
          scalar(1.5)*df_dx_t(c,u));
}
template <typename Cell>
typename Cell::CellScalar tau_yy_b(Cell c) {
  using scalar = typename Cell::CellScalar;
  auto u = [](typename Cell::CellValue val) { return uf(val); };
  auto v = [](typename Cell::CellValue val) { return vf(val); };
  return scalar(muf(c.value())) * (
          scalar(0.5)*df_dy_b(c,v) - 
          scalar(1.5)*df_dx_b(c,u));
}

template <typename Cell>
typename Cell::CellScalar tau_xy_r(Cell c) {
  using scalar = typename Cell::CellScalar;
  auto u = [](typename Cell::CellValue val) { return uf(val); };
  auto v = [](typename Cell::CellValue val) { return vf(val); };
  return scalar(muf(c.value())) * (df_dx_r(c,v) - df_dy_r(c,u));
}
template <typename Cell>
typename Cell::CellScalar tau_xy_l(Cell c) {
  using scalar = typename Cell::CellScalar;
  auto u = [](typename Cell::CellValue val) { return uf(val); };
  auto v = [](typename Cell::CellValue val) { return vf(val); };
  return scalar(muf(c.value())) * (df_dx_l(c,v) - df_dy_l(c,u));
}
template <typename Cell>
typename Cell::CellScalar tau_xy_t(Cell c) {
  using scalar = typename Cell::CellScalar;
  auto u = [](typename Cell::CellValue val) { return uf(val); };
  auto v = [](typename Cell::CellValue val) { return vf(val); };
  return scalar(muf(c.value())) * (df_dx_t(c,v) - df_dy_t(c,u));
}
template <typename Cell>
typename Cell::CellScalar tau_xy_b(Cell c) {
  using scalar = typename Cell::CellScalar;
  auto u = [](typename Cell::CellValue val) { return uf(val); };
  auto v = [](typename Cell::CellValue val) { return vf(val); };
  return scalar(muf(c.value())) * (df_dx_b(c,v) - df_dy_b(c,u));
}

/** Viscous Horizontal Flux
 */
template <typename Cell>
typename Cell::CellValue fv(Cell c) {
  typename Cell::CellValue res(c.value().size()); // reserve some space
  typename Cell::CellScalar ks = typename Cell::CellScalar(K);
  // resolve typenames with lambda function
  auto fcn = [](typename Cell::CellValue var) { return Tf(var); };
  res[0] = 0;
  res[1] = tau_xx_r(c) - tau_xx_l(c);
  res[2] = tau_xy_r(c) - tau_xy_l(c);
  // I'm using the cell velocity, rather than an average here -- is this OK?
  res[3] = tau_xx_r(c)*uf(c.value()) - tau_xx_l(c)*uf(c.value()) + 
           tau_xy_r(c)*vf(c.value()) - tau_xy_l(c)*vf(c.value()) - 
           ks*df_dx_r(c,fcn) + ks*df_dx_l(c,fcn);

  return res;
}

/** Viscous Vertical Flux
 */
template <typename Cell>
typename Cell::CellValue gv(Cell c) {
  typename Cell::CellValue res(c.value().size()); // reserve some space
  typename Cell::CellScalar ks = typename Cell::CellScalar(K);
  // resolve typenames with lambda function
  auto fcn = [](typename Cell::CellValue var) { return Tf(var); };
  res[0] = 0;
  res[1] = tau_xy_t(c) - tau_xy_b(c);
  res[2] = tau_yy_t(c) - tau_yy_b(c);
  // I'm using the cell velocity, rather than an average here -- is this OK?
  res[3] = tau_yy_t(c)*vf(c.value()) - tau_yy_b(c)*vf(c.value()) + 
           tau_xy_t(c)*uf(c.value()) - tau_xy_b(c)*uf(c.value()) - 
           ks*df_dy_t(c,fcn) + ks*df_dy_b(c,fcn);

  return res;
}

// 
// VISCOUS FLUX
// 
/** Computes the 2D viscous fluxes on a structured grid
 *  @brief Computes the viscous fluxes on a 
 *         structured grid via a dual-grid scheme
 * @tparam CellIter Cell iterator
 * @tparam value    State vector
 * 
 * @param cell_begin Beginning cell iterator
 * @param cell_end   Ending cell iterator
 * @param w          Container of flux vectors to which
 *                    we will write the euler fluxes
 * 
 * @pre std::distance(cell_begin,cell_end) == w.size()
 * 
 * @post w contains the euler fluxes for each cell
 */
// template <typename CellIter, typename Value>
template <typename CellIter>
void nsflux(CellIter cell_begin, CellIter cell_end, CellIter stage_begin) {
  /* --- SETUP --- */
  using Value = typename CellIter::Value;
  using scalar = typename Value::value_type;
  size_type cs = (*cell_begin).value().size();
  // Intermediate vectors
  Value Fx(cs), Fy(cs);

  /* --- ITERATE OVER CELLS --- */
  for ( ; cell_begin!=cell_end; ++cell_begin) {
    // Dereference cell
    auto c = *cell_begin;   // Read
    auto z = *stage_begin;  // Write
    // DEBUG -- Print cell index
// std::cout << "cell index=" << c.idx() << " (" << c.iy() << "," << c.jx() << ")";
// std::cout << std::endl;
    /* --- COMPUTE FLUXES --- */
    Fx = fv(c);
    Fy = gv(c);
    // Scale the fluxes
    Fx= scalar(1/c.dx())*Fx;
    Fy= scalar(1/c.dy())*Fy;
    // Add the result to the writeout vector
    z = scalar(-1)*(Fx+Fy)+z.value();
    // Iterate stage_begin to keep up with cell_begin
    ++stage_begin;
  }
  return;
}

#endif // NSFLUX
