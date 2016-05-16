#ifndef EFLUX // Include guard
#define EFLUX
/** Euler Flux Calculator
 */

#include <vector>     // data handling
#include <algorithm>  // std::copy(), std::transform()
#include <iostream>   // debug
#include <cmath>      // pow()
#include <valarray>   // std::valarray, std::begin, std::end

#include "gas_dynamics.hpp" // gas dynamics equations

using size_type = unsigned;

float eps = 0.25;  // artificial dissipation constant

/** Horizontal Wave Speed
 */
template <typename Value>
float wave_x(const Value& w) {
  return std::abs(w[1]/w[0])+cf(w);
}

/** Vertical Wave Speed
 */
template <typename Value>
float wave_y(const Value& w) {
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
  float P = pf(w);
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
  float P = pf(w);
  // Calculate flux elements
  out[0] = w[1];
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
  return typename Cell::CellValue(c.value(1,0)) - typename Cell::CellValue(c.value());
}
/** Vertical State Difference
 */
template <typename Cell>
typename Cell::CellValue dw_dy(Cell c) {
  return typename Cell::CellValue(c.value(0,-1)) - typename Cell::CellValue(c.value());
}

/** Jameson Horizontal Flux
 */
template <typename Cell>
typename Cell::CellValue fj(Cell c) {
  return typename Cell::CellScalar(0.5)*(
            f(typename Cell::CellValue(c.value(1,0))) 
          + f(typename Cell::CellValue(c.value()))      // physical flux
          ) + typename Cell::CellScalar(-eps/2)*(       // average wave speed
            wave_x(typename Cell::CellValue(c.value(1,0)))
          + wave_x(typename Cell::CellValue(c.value()))
          ) * dw_dx(c);   // state vector difference
}

/** Jameson Vertical Flux
 */
template <typename Cell>
typename Cell::CellValue gj(Cell c) {
  return typename Cell::CellScalar(0.5)*(
            g(typename Cell::CellValue(c.value(0,-1))) 
          + g(typename Cell::CellValue(c.value()))      // physical flux
          ) + typename Cell::CellScalar(-eps/2)*(       // average wave speed
            wave_y(typename Cell::CellValue(c.value(0,-1)))
          + wave_y(typename Cell::CellValue(c.value()))
          ) * dw_dy(c);   // state vector difference
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
    Fx = fj(c) + scalar(-1)*fj(c.neighbor(-1,0));
    Fy = gj(c) + scalar(-1)*gj(c.neighbor(0,1));
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
