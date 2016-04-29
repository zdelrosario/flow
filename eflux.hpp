#ifndef EFLUX // Include guard
#define EFLUX
/** Euler Flux Calculator
 */

#include <vector>     // data handling
#include <algorithm>  // std::copy(), std::transform()
#include <iostream>   // debug
#include <cmath>      // pow()
#include <valarray>   // std::valarray, std::begin, std::end

using size_type = unsigned;

float gam = 1.4; // adiabatic coefficient; assume calorically perfect gas
float K2 = 1;
float K4 = 1/32;
float C4 = 2;

// 
// VECTOR OPERATIONS
// 
/* Vector Addition v + w */
template <typename Value>
Value add(const Value& v, const Value& w) {
  Value u; u.resize(v.size());
  std::transform(std::begin(v),std::end(v),
                 std::begin(w),std::begin(u),
                 [](float a, float b){return a+b;});
  return u;
}
/* Vector Addition v - w */
template <typename Value>
Value sub(const Value& v, const Value& w) {
  Value u; u.resize(v.size());
  std::transform(std::begin(v),std::end(v),
                 std::begin(w),std::begin(u),
                 [](float a, float b){return a-b;});
  return u;
}
/* Scalar-Vector Multiplication c*v */
template <typename Scalar, typename Value>
Value mul(Scalar c, const Value& w) {
  Value u; u.resize(w.size());
  std::transform(std::begin(w),std::end(w),std::begin(u),
                 [&c](Scalar a){return c*a;});
  return u;
}

// 
// VECTOR HELPER FUNCTIONS
// 
/** Pressure
 * @brief Calculates the pressure based on a state vector
 */ 
template <typename Value>
float pf(const Value& w) {
  return (gam-1)*(w[3]-(pow(w[1],2)+pow(w[2],2))/(2*w[0]));
}

/** Speed of sound
 * @brief Calculates the local speed of sound based on a state vector
 */
template <typename Value>
float cf(const Value& w) {
  return sqrt( gam*pf(w)/w[0] );
}

/** Horizontal Wave Speed
 */
template <typename Value>
float wave_x(const Value& w) {
  return std::abs(w[1]/w[0])+cf(w);
}

/** Horizontal Wave Speed
 */
template <typename Value>
float wave_y(const Value& w) {
  return std::abs(w[2]/w[0])+cf(w);
}

/** Horizontal flux
 * @brief Calculates the horizontal flux based on a state vector
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
 * @brief Calculates the vertical flux based on a state vector
 * 
 * @param w   Input state vector
 * @param out Output flux vector
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
/** Horizontal Shock Sensor
 */
template <typename Cell>
float sense_x(Cell c) {
  return std::abs( (pf(c.value(1,0))+pf(c.value(-1,0))-2*pf(c.value())) / 
              (pf(c.value(1,0))+pf(c.value(-1,0))+2*pf(c.value())) );
}
/** Vertical Shock Sensor
 */
template <typename Cell>
float sense_y(Cell c) {
  return std::abs( (pf(c.value(0,1))+pf(c.value(0,-1))-2*pf(c.value())) / 
              (pf(c.value(0,1))+pf(c.value(0,-1))+2*pf(c.value())) );
}
/** Horizontal State Difference
 */
template <typename Cell>
typename Cell::CellValue dw_dx(Cell c) {
  return sub(typename Cell::CellValue(c.value(1,0)),
             typename Cell::CellValue(c.value()));
}
/** Vertical State Difference
 */
template <typename Cell>
typename Cell::CellValue dw_dy(Cell c) {
  return sub(typename Cell::CellValue(c.value(0,1)),
             typename Cell::CellValue(c.value()));
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
template <typename CellIter, typename Value>
void eflux(CellIter cell_begin, CellIter cell_end, std::vector<Value>& W) {
  /* --- SETUP --- */
  size_type cs = (*cell_begin).value().size();
  // Intermediate vectors
  Value h(cs),d(cs),Fx(cs), k(cs),e(cs),Fy(cs), temp(cs);
  float e2,e4,d2,d4,Sx,Sy,Rx,Ry;

  /* --- ITERATE OVER CELLS --- */
  for ( ; cell_begin!=cell_end; ++cell_begin) {
    // Dereference cell
    auto c = *cell_begin;
    // DEBUG -- Print cell index
// std::cout << "cell index=" << c.idx() << " (" << c.iy() << "," << c.jx() << ")";
// std::cout << std::endl;
    /* --- COMPUTE +1/2 BOUNDARIES --- */
    // Compute coefficients
    Sx = std::max(sense_x(c),sense_x(c.neighbor(1,0)));
    Sy = std::max(sense_y(c),sense_y(c.neighbor(0,1)));
    Rx = std::max(wave_x(c.value()),wave_x(c.value(1,0)));
    Ry = std::max(wave_y(c.value()),wave_y(c.value(0,1)));
    e2 = K2*Sx*Rx;
    e4 = std::max(float(0.),K4*Rx-C4*e2);
    d2 = K2*Sy*Ry;
    d4 = std::max(float(0.),K4*Ry-C4*d2);
    // Compute vector differences
    temp = sub(add(dw_dx(c.neighbor(1,0)),dw_dx(c.neighbor(-1,0))),mul(2,dw_dx(c)));
    d = sub(mul(e2,dw_dx(c)),mul(e4,temp));
    temp = sub(add(dw_dy(c.neighbor(0,1)),dw_dy(c.neighbor(0,-1))),mul(2,dw_dy(c)));
    e = sub(mul(d2,dw_dy(c)),mul(d4,temp));
    // Compute fluxes
    h = sub(mul( 0.5, add( f(Value(c.value(1,0))), f(Value(c.value())) ) ),d);
    k = sub(mul( 0.5, add( g(Value(c.value(0,1))), g(Value(c.value())) ) ),e);
    Fx= h;
    Fy= k;
    /* --- COMPUTE -1/2 BOUNDARIES --- */
    // Compute coefficients
    Sx = std::max(sense_x(c),sense_x(c.neighbor(-1,0)));
    Sy = std::max(sense_y(c),sense_y(c.neighbor(0,-1)));
    Rx = std::max(wave_x(c.value()),wave_x(c.value(-1,0)));
    Ry = std::max(wave_y(c.value()),wave_y(c.value(0,-1)));
    e2 = K2*Sx*Rx;
    e4 = std::max(float(0.),K4*Rx-C4*e2);
    d2 = K2*Sy*Ry;
    d4 = std::max(float(0.),K4*Ry-C4*d2);
    // Compute vector differences
    temp = sub(add(dw_dx(c),dw_dx(c.neighbor(-2,0))),mul(2,dw_dx(c.neighbor(-1,0))));
    d = sub(mul(e2,dw_dx(c.neighbor(-1,0))),mul(e4,temp));
    temp = sub(add(dw_dy(c),dw_dy(c.neighbor(0,-2))),mul(2,dw_dy(c.neighbor(0,-1))));
    e = sub(mul(d2,dw_dy(c.neighbor(0,-1))),mul(d4,temp));
    // Compute fluxes
    h = sub(mul( 0.5, add( f(Value(c.value(-1,0))), f(Value(c.value())) ) ),d);
    k = sub(mul( 0.5, add( g(Value(c.value(0,-1))), g(Value(c.value())) ) ),e);
    Fx= sub(Fx,h);
    Fy= sub(Fy,k);
    // Scale the fluxes
    Fx= mul(1/c.dx(),Fx);
    Fy= mul(1/c.dy(),Fy);
    // Add the result to the writeout vector
    W[c.idx()] = add(add(Fx,Fy),W[c.idx()]);
  }
  return;
}

#endif // EFLUX
