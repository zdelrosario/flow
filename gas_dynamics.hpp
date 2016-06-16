#ifndef GAS // Include guard
#define GAS

double gam = 1.4; // adiabatic coefficient; assume calorically perfect gas
double C1 = 1.458e-6; // Sutherland's law C_1
double S  = 110.4;    // Sutherland's law S
double Cp = 1.005e3;  // Isobaric specific heat capacity of air @ 300K
double Pr = 0.713;    // Prandtl number of air @ 300K
double R  = 287.058;  // Ideal gas constant of air
double K  = 0.0271;   // Thermal conductivity of air @ 300 K

// 
// STATE HELPER FUNCTIONS
// 
/* Pressure */ 
template <typename Value>
double pf(const Value& w) {
  return (gam-1)*(w[3]-(pow(w[1],2)+pow(w[2],2))/(2*w[0]));
}

/* Temperature */
template <typename Value>
double Tf(const Value& w) {
  return pf(w)/w[0]/R;
}

/* Horizontal velocity */
template <typename Value>
double uf(const Value& w) {
  return w[1]/w[0];
}

/* Vertical velocity */
template <typename Value>
double vf(const Value& w) {
  return w[2]/w[0];
}

/* Speed of sound */
template <typename Value>
double cf(const Value& w) {
  return sqrt( gam*pf(w)/w[0] );
}

/* Viscosity */
template <typename Value>
double muf(const Value& w) {
  // Sutherland's Law
  // return C1*pow(Tf(w),1.5) / (Tf(w)+S);
  // Fixed at 20 C
  // return 1.82e-5;
  // Inviscid
  (void) w; return 0.;
}

#endif // GAS
