#ifndef EFLUX // Include guard
#define EFLUX
/** Euler Flux Calculator
 */

#include <vector>     // data handling
#include <algorithm>  // std::abs(), std::max()
#include <iostream>   // debug
#include <cmath>      // pow()
#include <valarray>   // std::valarray, std::begin, std::end

#include "gas_dynamics.hpp" // gas dynamics equations

using size_type = unsigned;

double k2 = 1;
double k4 = 1.0/32.0; // Jameson recommends 1/32 for transonic flows
double c4 = 2;

/** Horizontal Wave Speed
 */
template <typename Value>
double wave_x(const Value& w) {
  return std::abs(w[1]/w[0])+cf(w);
}

/** Vertical Wave Speed
 */
template <typename Value>
double wave_y(const Value& w) {
  return std::abs(w[2]/w[0])+cf(w);
}

/** Horizontal flux
 * @brief Calculates the physical horizontal flux based on a state vector
 * 
 * @param w   Input state vector
 */
template <typename Value>
Value f(const Value& w) {
  Value out; out.resize(w.size());
  // Calculate pressure
  double P = pf(w);
  // Calculate flux elements
  out[0] = w[1];
  out[1] = pow(w[1],2)/w[0]+P;
  out[2] = w[1]*w[2]/w[0];
  out[3] = w[1]*(w[3]+P)/w[0];
  return out;
}

/** Vertical flux
 * @brief Calculates the physical vertical flux based on a state vector
 * 
 * @param w   Input state vector
 */
template <typename Value>
Value g(const Value& w) {
  Value out; out.resize(w.size());
  // Calculate pressure
  double P = pf(w);
  // Calculate flux elements
  out[0] = w[2];
  out[1] = w[1]*w[2]/w[0];
  out[2] = pow(w[2],2)/w[0] + P;
  out[3] = w[2]*(w[3]+P)/w[0];
  return out;
}

// 
// CELL HELPER FUNCTIONS
// 
/** Horizontal State Difference
 */
template <typename Cell>
typename Cell::CellValue dw_dx(Cell c) {
  using value = typename Cell::CellValue;
  return value(c.value(0,1)) - value(c.value());
}
/** Vertical State Difference
 */
template <typename Cell>
typename Cell::CellValue dw_dy(Cell c) {
  using value = typename Cell::CellValue;
  return value(c.value(-1,0)) - value(c.value());
}

/** Horizontal Speed Max
 */
template <typename Cell>
typename Cell::CellScalar rx(Cell c) {
  using value  = typename Cell::CellValue;
  return std::max( 
                wave_x(value(c.value(0,1))),
                wave_x(value(c.value())) );
}
/** Vertical Speed Max
 */
template <typename Cell>
typename Cell::CellScalar ry(Cell c) {
  using value  = typename Cell::CellValue;
  return std::max( 
                wave_y(value(c.value(-1,0))),
                wave_y(value(c.value())) );
}

/** Horizontal Pressure Sensor
 */
template <typename Cell>
typename Cell::CellScalar sx(Cell c) {
  return std::abs( (pf(c.value(0,1))-2*pf(c.value())+pf(c.value(0,-1))) / 
                   (pf(c.value(0,1))+2*pf(c.value())+pf(c.value(0,-1))) );
}
/** Vertical Pressure Sensor
 */
template <typename Cell>
typename Cell::CellScalar sy(Cell c) {
  return std::abs( (pf(c.value(1,0))-2*pf(c.value())+pf(c.value(-1,0))) / 
                   (pf(c.value(1,0))+2*pf(c.value())+pf(c.value(-1,0))) );
}

// Boundary helper
template <typename Cell>
typename Cell::CellScalar bcx(Cell c) {
  return (1.0-abs(c.bx()));
}
template <typename Cell>
typename Cell::CellScalar bcy(Cell c) {
  return (1.0-abs(c.by()));
}

// Dissipative Coefficients
template <typename Cell>
typename Cell::CellScalar eps2_x(Cell c) {
  // Standard JST
  // return k2 * sx(c) * rx(c);
  // Disable on boundary
  return k2 * sx(c) * rx(c) * (1-abs(c.bx()));
  // Disable
  // (void) c; return 0.;
}
template <typename Cell>
typename Cell::CellScalar eps2_y(Cell c) {
  // Standard JST
  // return k2 * sy(c) * ry(c);
  // Disable on boundary
  return k2 * sy(c) * ry(c) * (1-abs(c.by()));
  // Disable
  // (void) c; return 0.;
}
template <typename Cell>
typename Cell::CellScalar eps4_x(Cell c) {
  // Standard JST
  // return std::max(0.,k4*rx(c)-c4*eps2_x(c));

  // disable on boundary
  return k4*rx(c) * 
         (1-abs(c.bx())); 

  // disable on or near boundary
  // return k4*rx(c) * 
  //        (1-abs(c.neighbor(0,+1).bx())) *
  //        (1-abs(c.bx())) * 
  //        (1-abs(c.neighbor(0,-1).bx())); 

  // Disable
  (void) c; return 0.;
}
template <typename Cell>
typename Cell::CellScalar eps4_y(Cell c) {
  // Standard JST
  // return std::max(0.,k4*ry(c)-c4*eps2_y(c));

  // disable on boundary
  return k4*ry(c) * 
         (1-abs(c.by()));

  // disable on or near boundary
  // return k4*ry(c) * 
  //        (1-abs(c.neighbor(+1,0).by())) *
  //        (1-abs(c.by())) * 
  //        (1-abs(c.neighbor(-1,0).by()));

  // Disable
  // (void) c; return 0.;
}

// 4th Order Difference -- artificial dissipation
template <typename Cell>
typename Cell::CellValue d4x(Cell c) {
  // Standard JST
  return dw_dx(c.neighbor(0,1))
     -2.*dw_dx(c)
        +dw_dx(c.neighbor(0,-1));
  // Disable boundary terms
  // return dw_dx(c.neighbor(0,+1))*bcx(c.neighbor(0,+1))
  //    -2.*dw_dx(c)               *bcx(c)
  //       +dw_dx(c.neighbor(0,-1))*bcx(c.neighbor(0,-1));
}
template <typename Cell>
typename Cell::CellValue d4y(Cell c) {
  // Standard JST
  return dw_dy(c.neighbor(1,0))
     -2.*dw_dy(c)
        +dw_dy(c.neighbor(-1,0));
  // Disable boundary terms
  // return dw_dy(c.neighbor(+1,0))*bcy(c.neighbor(+1,0))
  //    -2.*dw_dy(c)               *bcy(c)
  //       +dw_dy(c.neighbor(-1,0))*bcy(c.neighbor(-1,0));
}

/** Jameson Horizontal Flux
 */
template <typename Cell>
typename Cell::CellValue fj(Cell c) {
  using scalar = typename Cell::CellScalar;
  using value  = typename Cell::CellValue;
  // Central flux + O(dx^2) dissipation + O(dx^4) dissipation
  return scalar(0.5)*( f(value(c.value(0,1))) + f(value(c.value())) ) 
    - eps2_x(c) * dw_dx(c)
    + eps4_x(c) * d4x(c);
}

/** Jameson Vertical Flux
 */
template <typename Cell>
typename Cell::CellValue gj(Cell c) {
  using scalar = typename Cell::CellScalar;
  using value  = typename Cell::CellValue;
  // Central flux + O(dy^2) dissipation + O(dy^4) dissipation
  return scalar(0.5)*( g(value(c.value(-1,0))) + g(value(c.value())) ) 
    - eps2_y(c) * dw_dy(c)
    + eps4_y(c) * d4y(c);
}

// 
// EULER FLUX
// 
/** Computes the 2D euler fluxes on a structured grid
 *  @brief Implements the JST scheme on a 
 *         2D structured grid
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
void eflux(CellIter cell_begin, CellIter cell_end, CellIter stage_begin) {
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
    Fx = fj(c) + scalar(-1)*fj(c.neighbor(0,-1));
    Fy = gj(c) + scalar(-1)*gj(c.neighbor(1,0));
    // Scale the fluxes by spatial discretization
    Fx= scalar(1/c.dx())*Fx;
    Fy= scalar(1/c.dy())*Fy;
    // Add the result to the writeout vector
    z = scalar(-1)*(Fx+Fy)+z.value();
    // Iterate stage_begin to keep up with cell_begin
    ++stage_begin;
  }
  return;
}

#endif // EFLUX
