#ifndef NSFLUX // Include guard
#define NSFLUX
/** Viscous Flux Calculator
 */

#include <vector>     // data handling
#include <algorithm>  // std::copy(), std::transform()
#include <iostream>   // debug
#include <cmath>      // pow()
#include <valarray>   // std::valarray, std::begin, std::end

using size_type = unsigned;

// float gam = 1.4; // adiabatic coefficient; assume calorically perfect gas
float C1 = 1.458e-6; // Sutherland's law C_1
float S  = 110.4;    // Sutherland's law S
float Cp = 1.005e3;  // Isobaric specific heat capacity of air @ 300K
float Pr = 0.713;    // Prandtl number of air @ 300K
float R  = 287.058;  // Ideal gas constant of air
float K  = 0.0271;   // Thermal conductivity of air @ 300 K

// 
// STATE HELPER FUNCTORS
// 
/* Pressure */ 
struct pf {
  template <typename Value>
  float operator()(const Value& w) {
    return (gam-1)*(w[3]-(pow(w[1],2)+pow(w[2],2))/(2*w[0]));
  }
};

/* Temperature */
struct Tf {
  template <typename Value>
  float operator()(const Value& w) {
    return pf(w)/w[0]/R;
  }
};

/* Horizontal velocity */
struct uf {
  template <typename Value>
  float operator()(const Value& w) {
    return w[1]/w[0];
  }
};

/* Vertical velocity */
struct vf {
  template <typename Value>
  float operator()(const Value& w) {
    return w[2]/w[0];
  }
};

/* Speed of sound */
struct cf {
  template <typename Value>
  float operator()(const Value& w) {
    return sqrt( gam*pf(w)/w[0] );
  }
};

/* Viscosity */
struct muf {
  template <typename Value>
  float operator()(const Value& w) {
    Tf T;
    return C1*pow(T(w),1.5) / (T(w)+S);
  }
};

// 
// DUAL GRID DERIVATIVES
// 
/* --- HORIZONTAL DERIVATIVES --- */
/* --- du_dx --- */
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dx_r(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25)*(
        // Upper right
          f(c.value(1,0))  - f(c.value()) + 
          f(c.value(1,-1)) - f(c.value(0,-1)) +
        // Lower right
          f(c.value(1,0)) - f(c.value()) + 
          f(c.value(1,1)) - f(c.value(0,1))
        );
}
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dx_l(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25)*(
        // Upper left
          f(c.value())     - f(c.value(-1,0)) + 
          f(c.value(0,-1)) - f(c.value(-1,-1)) + 
        // Lower left
          f(c.value())    - f(c.value(-1,0)) + 
          f(c.value(0,1)) - f(c.value(-1,1))
        );
}
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dx_t(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25)*(
        // Upper right
          f(c.value(1,-1))  - f(c.value(-1,-1)) + 
          f(c.value(1,0)) - f(c.value(-1,0))
        );
}
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dx_b(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25)*(
        // Upper left
          f(c.value(1,1)) - f(c.value(-1,1)) + 
          f(c.value(1,0)) - f(c.value(-1,0))
        );
}

/* --- VERTICAL DERIVATIVES --- */
/* --- du_dy --- */
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dy_t(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25)*(
        // Upper right
          f(c.value(1,-1)) - f(c.value(1,0)) + 
          f(c.value(0,-1)) - f(c.value(0,0)) +
        // Upper left
          f(c.value(0,-1))  - f(c.value(0,0)) + 
          f(c.value(-1,-1)) - f(c.value(-1,0))
        );
}
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dy_b(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25)*(
        // Bottom right
          f(c.value(1,0)) - f(c.value(1,1)) + 
          f(c.value(0,0)) - f(c.value(0,1)) +
        // Bottom left
          f(c.value(0,0))  - f(c.value(0,1)) + 
          f(c.value(-1,0)) - f(c.value(-1,1))
        );
}
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dy_l(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25)*(
          f(c.value(-1,-1)) - f(c.value(-1,1)) + 
          f(c.value(0,-1)) - f(c.value(0,1))
        );
}
template <typename Cell, typename Fcn>
typename Cell::CellScalar df_dy_r(Cell c, Fcn f) {
  return typename Cell::CellScalar(0.25)*(
          f(c.value(1,-1)) - f(c.value(1,1)) + 
          f(c.value(0,-1)) - f(c.value(0,1))
        );
}

// 
// VISCOUS STRESS TERMS
// 
template <typename Cell>
typename Cell::CellScalar tau_xx_r(Cell c) {
  muf mu;
  return typename Cell::CellScalar(mu(c.value())) * (
          typename Cell::CellScalar(0.5)*df_dx_r(c,uf()) - 
          typename Cell::CellScalar(1.5)*df_dy_r(c,vf()));
}
template <typename Cell>
typename Cell::CellScalar tau_xx_l(Cell c) {
  muf mu;
  return typename Cell::CellScalar(mu(c.value())) * (
          typename Cell::CellScalar(0.5)*df_dx_l(c,uf()) - 
          typename Cell::CellScalar(1.5)*df_dy_l(c,vf()));
}

template <typename Cell>
typename Cell::CellScalar tau_yy_t(Cell c) {
  muf mu;
  return typename Cell::CellScalar(mu(c.value())) * (
          typename Cell::CellScalar(0.5)*df_dy_t(c,vf()) - 
          typename Cell::CellScalar(1.5)*df_dx_t(c,uf()));
}
template <typename Cell>
typename Cell::CellScalar tau_yy_b(Cell c) {
  muf mu;
  return typename Cell::CellScalar(mu(c.value())) * (
          typename Cell::CellScalar(0.5)*df_dy_b(c,vf()) - 
          typename Cell::CellScalar(1.5)*df_dx_b(c,uf()));
}

template <typename Cell>
typename Cell::CellScalar tau_xy_r(Cell c) {
  muf mu;
  return typename Cell::CellScalar(mu(c.value())) * (
          typename Cell::CellScalar(0.5)*df_dx_r(c,vf()) - 
          typename Cell::CellScalar(0.5)*df_dy_r(c,uf()));
}
template <typename Cell>
typename Cell::CellScalar tau_xy_l(Cell c) {
  muf mu;
  return typename Cell::CellScalar(mu(c.value())) * (
          typename Cell::CellScalar(0.5)*df_dx_l(c,vf()) - 
          typename Cell::CellScalar(0.5)*df_dy_l(c,uf()));
}
template <typename Cell>
typename Cell::CellScalar tau_xy_t(Cell c) {
  muf mu;
  return typename Cell::CellScalar(mu(c.value())) * (
          typename Cell::CellScalar(0.5)*df_dx_t(c,vf()) - 
          typename Cell::CellScalar(0.5)*df_dy_t(c,uf()));
}
template <typename Cell>
typename Cell::CellScalar tau_xy_b(Cell c) {
  muf mu;
  return typename Cell::CellScalar(mu(c.value())) * (
          typename Cell::CellScalar(0.5)*df_dx_b(c,vf()) - 
          typename Cell::CellScalar(0.5)*df_dy_b(c,uf()));
}

/** Horizontal flux
 * @brief Calculates the physical horizontal flux based on a state vector
 * 
 * @param w   Input state vector
 */
// template <typename Value>
// Value f(const Value& w) {
//   Value out; out.resize(w.size());
//   // Calculate pressure
//   float P = pf(w);
//   // Calculate flux elements
//   out[0] = w[1];
//   out[1] = pow(w[1],2)/w[0]+P;
//   out[2] = w[1]*w[2]/w[0];
//   out[3] = w[1]*(w[3]+P)/w[0];
//   return out;
// }

/** Vertical flux
 * @brief Calculates the physical vertical flux based on a state vector
 * 
 * @param w   Input state vector
 */


// 
// CELL HELPER FUNCTIONS
// 
/** Horizontal State Difference
 */
// template <typename Cell>
// typename Cell::CellValue dw_dx(Cell c) {
//   return typename Cell::CellValue(c.value(1,0)) - typename Cell::CellValue(c.value());
// }
/** Vertical State Difference
 */
// template <typename Cell>
// typename Cell::CellValue dw_dy(Cell c) {
//   return typename Cell::CellValue(c.value(0,1)) - typename Cell::CellValue(c.value());
// }

/** Viscous Horizontal Flux
 */
template <typename Cell>
typename Cell::CellValue fv(Cell c) {
  typename Cell::CellValue res(c.value().size()); // reserve some space
  uf u; vf v;
  typename Cell::CellScalar ks = typename Cell::CellScalar(K);

  res[0] = 0;
  res[1] = tau_xx_r(c) - tau_xx_l(c);
  res[2] = tau_xy_r(c) - tau_xy_l(c);
  // I'm using the cell velocity, rather than an average here -- is this OK?
  res[3] = tau_xx_r(c)*u(c.value()) - tau_xx_l(c)*u(c.value()) + 
           tau_xy_r(c)*v(c.value()) - tau_xy_l(c)*v(c.value()) - 
           ks*df_dx_r(c,Tf()) + ks*df_dx_l(c,Tf());

  return res;
}

/** Viscous Vertical Flux
 */
template <typename Cell>
typename Cell::CellValue gv(Cell c) {
  typename Cell::CellValue res(c.value().size()); // reserve some space
  uf u; vf v;
  typename Cell::CellScalar ks = typename Cell::CellScalar(K);

  res[0] = 0;
  res[1] = tau_xy_t(c) - tau_xy_b(c);
  res[2] = tau_yy_t(c) - tau_yy_b(c);
  // I'm using the cell velocity, rather than an average here -- is this OK?
  res[3] = tau_yy_t(c)*v(c.value()) - tau_yy_b(c)*v(c.value()) + 
           tau_xy_t(c)*u(c.value()) - tau_xy_b(c)*u(c.value()) - 
           ks*df_dy_t(c,Tf()) + ks*df_dy_b(c,Tf());

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
