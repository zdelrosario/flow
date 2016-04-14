#ifndef NEWTON // Include guard
#define NEWTON

#include <iostream>

/** Newton-Rhapson rootfinding method
 * @brief Finds the root of a function via
 *        successive local linear approximation
 *
 * @tparam V  Input & output type for F and dF
 * @tparam F  Function from V->V
 * @tparam dF Function from V->V
 *
 * @param max Maximum iteraion count (default=1000)
 * @param eps Convergence tolerance, as measured
 *            by absolute function value (default=1e-6)
 * 
 * @pre x0 within the region of convergence for
 *      Newton-Rhapson on f
 * @pre df = d/dx(f)
 * 
 * @return res Root of provided f
 */
template <typename V, typename F, typename dF>
V newton(V x0, F f, dF df, double eps=1e-6, unsigned long max=1000) {
  // Setup
  V fk = f(x0); V fpk= df(x0); V delta;
  unsigned long it = 0;
  // Continue while not converged and 
  // under iteration limit
  while ((std::abs(fk)>eps)&&(it<max)) {
    // Take newton step
    delta = -fk/fpk;
    x0    = x0 + delta;
    // Update function values
    fk  = f(x0);
    fpk = df(x0);
    // Iterate the counter
    ++it;
  }
  return x0;
}

/** Secant method
 */
template <typename V, typename F>
V secant(V x0, V x1, F f, double eps=1e-6, unsigned long max=1000) {
  // Setup
  V f0 = f(x0); V f1 = f(x1); V delta;
  unsigned long it = 0;
  // Continue while not converged and 
  // under iteration limit
  while ((std::abs(f1)>eps)&&(it<max)) {
    // Take secant step
    delta = -f1 * (x1-x0)/(f1-f0);
    x0 = x1;
    x1 = x1 + delta;
    // Update function values
    f0 = f1;
    f1 = f(x1);
    // Iterate the counter
    ++it;
  }
  std::cout << "it=" << it << std::endl;
  return x1;
}

#endif // NEWTON

/** TEST CODE
#include <iostream>
#include <cmath>

#include "newton.hpp"

int main() {
  float x0 = 1.0;
  float xs = newton( x0, 
                     [](float x) {return cos(x);}, 
                     [](float x) {return -sin(x);} );

  std::cout << xs << std::endl;
}
 */
