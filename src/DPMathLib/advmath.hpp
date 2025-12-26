/**
 * This is the header containing the declaration of the function
 * of advanced mathematical tools not present in the C++ standard libararies
 *
 *  - Integration algorithms
 *  - Special functions definitions
 */
#pragma once

#ifndef __cplusplus
  #define __cplusplus
#endif

#include "darkpack.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <future>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <iostream>


namespace advmath
{

  using real_t = dpk::real_t;
  using complex_floating_t = std::complex<real_t>;


  inline real_t getdiscrepancy(const real_t res1, const real_t res2)
  {
    real_t denom = res1 + res2;
    return (denom != 0.) ? std::abs((res1 - res2) / denom) : 0.;
  }

  real_t integrate_gauss(real_t xmin, real_t xmax, std::function<real_t(real_t)> f, int order = 5);

  real_t integrate_gauss_comparative(real_t xmin, real_t xmax, std::function<real_t(real_t)> f,
                                     real_t maxdiscrepancy = -1.0e-3, real_t* discr = nullptr);

  real_t integrate_trap(real_t xmin, real_t xmax, std::function<real_t(real_t)> f, real_t maxdiscrepancy = -1.0e-3,
                        real_t* discr = nullptr);

  real_t find_zero_secant(real_t xmin, real_t xmax, std::function<real_t(real_t)> f, real_t tol);

  complex_floating_t digamma(complex_floating_t z);

  real_t I0(real_t x);
  real_t I1(real_t x);
  real_t K0(real_t x);
  real_t K1(real_t x);
  real_t K2(real_t x);
  real_t K3(real_t x);
  real_t K4(real_t x);
  real_t Lbessel(real_t x);
  real_t Mbessel(real_t x);
  real_t Nbessel(real_t x);
  real_t expcor(real_t x);
  real_t K0exp(real_t x, real_t z);
  real_t K1exp(real_t x, real_t z);
  real_t K2exp(real_t x, real_t z);
  real_t K0exp_numdenom(real_t num_x, real_t num_z, real_t denom);
  real_t K1exp_numdenom(real_t num_x, real_t num_z, real_t denom);
  real_t K2exp_numdenom(real_t num_x, real_t num_z, real_t denom);
  real_t besselj0(real_t x);
  real_t besselj1(real_t x);



} // End of namespace advmath
