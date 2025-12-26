#include "advmath.hpp"
template <typename T> static inline T SQUARE(const T x) { return x * x; }

namespace advmath
{
  real_t integrate_gauss(real_t xmin, real_t xmax, std::function<real_t(real_t)> f, int order)
  {
    constexpr const real_t abs5[5] = {-0.9061798459, -0.5384693101, 0., 0.5384693101, 0.9061798459};
    constexpr const real_t wei5[5] = {0.2369268851, 0.4786286705, 0.5688888888, 0.4786286705, 0.2369268851};

    constexpr const real_t abs6[6] = {-0.9324695142, -0.6612093865, -0.2386191861,
                                      0.2386191861,  0.6612093865,  0.9324695142};
    constexpr const real_t wei6[6] = {0.1713244924, 0.3607615730, 0.4679139346,
                                      0.4679139346, 0.3607615730, 0.1713244924};

    constexpr const real_t abs7[7] = {-0.9491079123, -0.7415311856, -0.4058451514, 0.,
                                      0.4058451514,  0.7415311856,  0.9491079123};
    constexpr const real_t wei7[7] = {0.1294849662, 0.2797053915, 0.3818300505, 0.4179591837,
                                      0.3818300505, 0.2797053915, 0.1294849662};

    constexpr const real_t abs8[8] = {-0.9602898565, -0.7966664774, -0.5255324099, -0.1834346425,
                                      0.1834346425,  0.5255324099,  0.7966664774,  0.9602898565};
    constexpr const real_t wei8[8] = {0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834,
                                      0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363};

    constexpr const real_t wei9[9] = {0.3302393550012598, 0.1806481606948574, 0.1806481606948574,
                                      0.0812743883615744, 0.0812743883615744, 0.3123470770400029,
                                      0.3123470770400029, 0.2606106964029354, 0.2606106964029354};
    constexpr const real_t abs9[9] = {0.,
                                      -0.8360311073266358,
                                      0.8360311073266358,
                                      -0.9681602395076261,
                                      0.9681602395076261,
                                      -0.3242534234038089,
                                      0.3242534234038089,
                                      -0.6133714327005904,
                                      0.6133714327005904};

    const real_t *abs, *wei;

    switch (order)
    {
    case 5:
      abs = abs5;
      wei = wei5;
      break;

    case 6:
      abs = abs6;
      wei = wei6;
      break;

    case 7:
      abs = abs7;
      wei = wei7;
      break;

    case 8:
      abs = abs8;
      wei = wei8;
      break;

    case 9:
      abs = abs9;
      wei = wei9;
      break;

    default:
      abs = abs5;
      wei = wei5;
      order = 5;
      break;
    };

    std::future<real_t> fut[9];
    // Parameters for the linear transformation of the interval
    const real_t m = 0.5 * (xmax - xmin);
    const real_t q = 0.5 * (xmax + xmin);

    for (int i = 0; i < order; i++)
    {
      // Transforming the integration interval
      real_t abscissa = m * abs[i] + q;
      fut[i] = std::async(std::launch::async, f, abscissa);
    }

    real_t integral = 0.;

    for (int i = 0; i < order; i++)
      integral += wei[i] * fut[i].get();

    return integral * m;
  }

  real_t integrate_gauss_comparative(real_t xmin, real_t xmax, std::function<real_t(real_t)> f, real_t maxdiscrepancy,
                                     real_t* err)
  {
    real_t res1 = integrate_gauss(xmin, xmax, f, 5);
    real_t res2 = integrate_gauss(xmin, xmax, f, 6);
    real_t discrepancy = getdiscrepancy(res1, res2);

    int order = 7;
    while (discrepancy > maxdiscrepancy && order <= 9)
    {
      res1 = res2;
      res2 = integrate_gauss(xmin, xmax, f, order++);
      discrepancy = getdiscrepancy(res1, res2);
    };

    if (err != nullptr)
      *err = discrepancy;
    return res2;
  }

  real_t integrate_trap(real_t xmin, real_t xmax, std::function<real_t(real_t)> f, real_t maxdiscrepancy, real_t* discr)
  {
    const int MAXPOINTS = 257;
    int npoints = 3, prev_npoints = 0, i, j, newpoint;
    real_t res1 = 0., res2 = 0., deltax, abscissa, evaluation, discrepancy;
    real_t abscissas_sorted[MAXPOINTS];
#ifdef DEBUG2
    real_t abscissas_unsorted[MAXPOINTS], evaluations_unsorted[MAXPOINTS];
#endif
    std::future<real_t> evaluation_future[MAXPOINTS];

    // First evaluations
    deltax = (xmax - xmin) / (npoints - 1.); // First value for the stepsize
#ifdef DEBUG2
    printf("stepsize is %.4e\n", deltax);
#endif

    // #pragma omp for
    for (i = 0; i < npoints; i++)
    {
      abscissa = xmin + i * deltax;

#ifdef DEBUG
      printf("Adding value %d\n", i);
#endif

      abscissas_sorted[i] = abscissa;

#ifdef DEBUG2
      abscissas_unsorted[i] = abscissa;
#endif
#ifdef DEBUG
      printf("Creating asyncronous threads\n");
#endif

      evaluation_future[i] = std::async(std::launch::async, f, abscissa);
    }

#ifdef DEBUG
    printf("Created asyncronous threads\n");
#endif

    for (i = 0; i < npoints; i++)
    {

#ifdef DEBUG
      printf("Getting values from thread %d\n", i);
#endif

      evaluation = evaluation_future[i].get();

#ifdef DEBUG
      printf("Setting evaluation\n");
#endif

      evaluation = std::isnormal(evaluation) ? evaluation : 0.0;
      res1 += ((i == 0 || i == npoints - 1) ? 0.5 : 1.) * evaluation;

#ifdef DEBUG2
      evaluations_unsorted[i] = evaluation;
#endif
    }

    res2 = res1;    // res2 contains the sum of the values
    res1 *= deltax; // First estimate

    //     int prev_start = prev_npoints;
    prev_npoints = npoints;
    npoints = (npoints - 1) * 2 + 1;

    do // Next evaluations
    {
      if (npoints > MAXPOINTS)
      {
        if (discr != nullptr)
          *discr = discrepancy;
        return res1;
        break;
      }

#ifdef DEBUG
      //       printf("Entering cycle with npoints = %d, prev_npoints = %d, prev_start = %d\n", npoints,
      //       prev_npoints, prev_start);
      printf("Entering cycle with npoints = %d, prev_npoints = %d\n", npoints, prev_npoints);
      // Printing out the abscissas
      for (int k = 0; k < prev_npoints; k++)
        printf("abscissas_sorted[%d] = %.4e\n", k, abscissas_sorted[k]);
#endif

      j = 0;
      for (i = 0; i < prev_npoints - 1; i++)
      {
        newpoint = prev_npoints + j;

        // Note: deltax is the old stepsize
        if (getdiscrepancy(std::fabs((abscissas_sorted[i] - abscissas_sorted[i + 1])), deltax) > 1.0e-5)
        {
          printf("i = %d, error on the old stepsize: %.4e %4.e\n", i, abscissas_sorted[i + 1] - abscissas_sorted[i],
                 deltax);
#ifdef DEBUG
          for (k = 0; k < newpoint - 1; k++)
            printf("abscissas_sorted[%d] = %.4e\n", k, abscissas_sorted[k]);
#endif
#ifdef DEBUG2
          for (k = 0; k < newpoint - 1; k++)
            printf("abscissas_unsorted[%d] = %.4e\n", k, abscissas_unsorted[k]);
#endif
          exit(2);
        }

        abscissa = (abscissas_sorted[i] + abscissas_sorted[i + 1]) / 2.; // Getting the halfpoint

        if (newpoint >= npoints)
        {
          printf("Error on the array size!\n");
          exit(2);
        }
#ifdef DEBUG
        printf("Adding value %d 1\n", newpoint);
#endif
        abscissas_sorted[newpoint] = abscissa;
#ifdef DEBUG
        printf("Adding value %d 2\n", newpoint);
#endif
#ifdef DEBUG2
        abscissas_unsorted[newpoint] = abscissa;
#endif

        // Creating asyncronous threads
        evaluation_future[newpoint] = std::async(std::launch::async, f, abscissa);

        j++;
      }
      // Checking if all the evaluation are correct
      if (j + prev_npoints != npoints)
      {
        printf("Error on j + prev_npoints= %d\n", j + prev_npoints);
        exit(1);
      }

      deltax /= 2.; // Halving the stepsize
#ifdef DEBUG
      printf("deltax = %.4e, ", deltax);
#endif

      j = 0;
      for (i = 0; i < prev_npoints - 1; i++)
      {
        newpoint = prev_npoints + j;
        evaluation = evaluation_future[newpoint].get();
        evaluation = std::isnormal(evaluation) ? evaluation : 0.0;
        res2 += evaluation;
#ifdef DEBUG
        printf("Adding value %d 3\n", newpoint);
#endif
#ifdef DEBUG2
        evaluations_unsorted[newpoint] = evaluation;
#endif
        j++;
      }

#ifdef DEBUG
      printf("Getting res2, ");
#endif
      res2 *= deltax;

#ifdef DEBUG
      printf(" the discrepancy...");
#endif
      discrepancy = getdiscrepancy(res1, res2);
      if (discrepancy < maxdiscrepancy)
      {
        if (discr != nullptr)
          *discr = discrepancy;
        return res2;
        break;
      }

#ifdef DEBUG2
      printf("result is %.4e +- %.4e (%.4e). Bad comparison, resetting things\n", kinfac_onlysqrtS * res2,
             kinfac_onlysqrtS * discrepancy * res2, discrepancy);
#endif
      res1 = res2;
      res2 /= deltax;
#ifdef DEBUG
      printf("Sorting the array of abscissas\n");
#endif
      std::sort(&abscissas_sorted[0], &abscissas_sorted[npoints]);
#ifdef DEBUG
      for (k = 0; k < npoints; k++)
        printf("abscissas_sorted[%d] = %.4e\n", k, abscissas_sorted[k]);
  #ifdef DEBUG2
      for (k = 0; k < npoints; k++)
        printf("abscissas_unsorted[%d] = %.4e\n", k, abscissas_unsorted[k]);
  #endif
#endif
      //       prev_start=prev_npoints;
      prev_npoints = npoints;
      npoints = (npoints - 1) * 2 + 1;
    } while (true);
  }

  real_t find_zero_secant(real_t xmin, real_t xmax, std::function<real_t(real_t)> f, real_t tol)
  {
    if (xmax < xmin)
      std::swap(xmin, xmax);

    real_t fmin = f(xmin);
    real_t fmax = f(xmax);
    real_t xnew, fnew;
    const int MAXITER = 1000;
    int iter = 0;

    if (fmin * fmax > 0.)
    { // No zero crossing
      return NAN;
    }

    while (iter < MAXITER)
    {
      xnew = xmax - fmax * (xmax - xmin) / (fmax - fmin);
      fnew = f(xnew);

      if (std::abs(fnew) < tol)
        return xnew;

      if (fnew * fmax < 0.)
      {
        xmin = xmax;
        fmin = fmax;
      }

      xmax = xnew;
      fmax = fnew;

      iter++;
    }
  }

  complex_floating_t digamma(complex_floating_t z)
  {

    // Source
    // Paul Godfrey (2022). Special Functions math library
    // (https://www.mathworks.com/matlabcentral/fileexchange/978-special-functions-math-library), MATLAB Central
    // File Exchange.
    if (round(z.real()) == z.real() && round(z.imag()) == z.imag() && z.real() <= 0. && z.imag() == 0.)
      return INFINITY;

    const complex_floating_t zz = z;
#ifdef DEBUG
    std::cout << "Argument is " << zz << std::endl;
#endif
    const bool reflect = (z.real() < 0.5);
    if (reflect)
    {
#ifdef DEBUG
      std::cout << "Reflecting z\n";
#endif
      z = 1. - z;
    }

    const real_t g = 607. / 128.;

    const real_t c[15] = {0.99999999999999709182,    57.156235665862923517,     -59.597960355475491248,
                          14.136097974741747174,     -0.49191381609762019978,   .33994649984811888699e-4,
                          .46523628927048575665e-4,  -.98374475304879564677e-4, .15808870322491248884e-3,
                          -.21026444172410488319e-3, .21743961811521264320e-3,  -.16431810653676389022e-3,
                          .84418223983852743293e-4,  -.26190838401581408670e-4, .36899182659531622704e-5};

    complex_floating_t n = 0., d = 0.;

    // Running on all the array except the 1st element
    for (int k = 14; k >= 1; k--)
    {
      complex_floating_t dz = 1. / (z + static_cast<complex_floating_t>(k) - 1.);
      // I replaced k -2 -> k -1 because the size is decreased by 1, so it is k
      complex_floating_t dd = c[k] * dz;
#ifdef DEBUG
      std::cout << "k = " << k << " dz = " << dz << " dd = " << dd << std::endl;
#endif
      d = d + dd;
      n = n - dd * dz;
    }
    d = d + c[0]; // First element of the array

    const complex_floating_t gg = z + g - 0.5;

    complex_floating_t digamma = log(gg) + (n / d - g / gg) - reflect * M_PI / tan(M_PI * zz);

#ifdef DEBUG
    std::cout << "\n result is \n";
#endif


    return digamma;
  }


  /*-----------------------------------------------------------------------*/
  /* Modified Bessel functions of second type and of order 0,1,2 */
  /*-----------------------------------------------------------------------*/

  real_t I0(real_t x)
  /* calculates the Modified Bessel functions of first type and of order 0 of x */
  {
    real_t y;

    real_t absx = fabs(x);

    if (absx < 3.75)
    {
      y = x / 3.75;
      y = y * y;
      return 1. +
             y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
    }
    else
    {
      y = 3.75 / absx;
      return (exp(absx) / sqrt(absx)) *
             (0.39894228 +
              y * (0.1328592e-1 +
                   y * (0.225319e-2 +
                        y * (-0.157565e-2 +
                             y * (0.916281e-2 +
                                  y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
    }
  }

  real_t I1(real_t x)
  /* calculates the Modified Bessel functions of first type and of order 1 of x */
  {
    real_t y, tmp;
    real_t I1 = 0.;

    real_t absx = fabs(x);

    if (absx < 3.75)
    {
      y = x / 3.75;
      y = y * y;
      I1 = absx *
           (0.5 + y * (0.87890594 +
                       y * (0.51498869 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
    }
    else
    {
      y = 3.75 / absx;
      tmp = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
      tmp = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * tmp))));
      I1 *= (exp(absx) / sqrt(absx));
    }
    if (x < 0.)
      return -I1;
    else
      return I1;
  }

  real_t K0(real_t x)
  /* calculates the Modified Bessel functions of second type and of order 0 of x */
  {
#ifdef GSL_BESSEL_FUNCTIONS
    gsl_sf_result result;
    gsl_sf_bessel_K0_e(x, &result);
    return result.val;
#else
    real_t y;

    if (x <= 2.)
    {
      y = x * x / 4.;
      return (-log(x / 2.) * I0(x)) +
             (-0.57721566 +
              y * (0.42278420 +
                   y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5))))));
    }
    else
    {
      y = 2. / x;
      return (exp(-x) / sqrt(x)) *
             (1.25331414 +
              y * (-0.7832358e-1 +
                   y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3))))));
    }
#endif
  }

  real_t K1(real_t x)
  /* calculates the Modified Bessel functions of second type and of order 1 of x */
  {
#ifdef GSL_BESSEL_FUNCTIONS
    gsl_sf_result result;
    gsl_sf_bessel_K1_e(x, &result);
    return result.val;
#else
    real_t y;

    if (x <= 2.)
    {
      y = x * x / 4.;
      return (log(x / 2.) * I1(x)) +
             (1. / x) *
                 (1. + y * (0.15443144 +
                            y * (-0.67278579 +
                                 y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
    }
    else
    {
      y = 2. / x;
      return (exp(-x) / sqrt(x)) *
             (1.25331414 +
              y * (0.23498619 + y * (-0.3655620e-1 +
                                     y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
    }
#endif
  }

  real_t K2(real_t x)
  /* calculates the Modified Bessel functions of second type and of order 2 of x */
  {
#ifdef GSL_BESSEL_FUNCTIONS
    gsl_sf_result result;
    gsl_sf_bessel_Kn_e(2, x, &result);
    return result.val;
#else
    return K0(x) + 2. / x * K1(x);
#endif
  }

  real_t K3(real_t x)
  /* calculates the Modified Bessel functions of second type and of order 3 of x */ { return K1(x) + 4. / x * K2(x); }

  real_t K4(real_t x)
  /* calculates the Modified Bessel functions of second type and of order 4 of x */ { return K2(x) + 6. / x * K3(x); }

  real_t Lbessel(real_t x)
  /* calculates L(x)=K2(x)/x */ { return K2(x) / x; }

  real_t Mbessel(real_t x)
  /* calculates M(x)=(3*K3(x)+K1(x))/4x */ { return (0.75 * K3(x) + 0.25 * K1(x)) / x; }

  real_t Nbessel(real_t x)
  /* calculates N(x)=(K4(x)+K2(x))/2x */ { return (0.5 * K4(x) + 0.5 * K2(x)) / x; }

  real_t K0exp(real_t x, real_t z)
  /* calculates the extended Modified Bessel functions of second type and of order 0 of x and z */
  /* K0exp(x,z)=K0(x)*exp(z)*sqrt(z) */
  {
    real_t y;

    if (x <= 2.)
    {
      y = x * x / 4.;
      return ((-log(x / 2.) * I0(x)) +
              (-0.57721566 +
               y * (0.42278420 +
                    y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5))))))) *
             exp(z) * sqrt(z);
    }
    else
    {
      y = 2. / x;
      return (exp(z - x) / sqrt(x / z)) *
             (1.25331414 +
              y * (-0.7832358e-1 +
                   y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3))))));
    }
  }

  real_t K1exp(real_t x, real_t z)
  /* calculates the extended Modified Bessel functions of second type and of order 1 of x and z */
  /* K1exp(x,z)=K1(x)*exp(z)*sqrt(z) */
  {
    real_t y;

    if (x <= 2.)
    {
      y = x * x / 4.;
      return ((log(x / 2.) * I1(x)) +
              (1. / x) *
                  (1. + y * (0.15443144 +
                             y * (-0.67278579 +
                                  y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))))) *
             exp(z) * sqrt(z);
    }
    else
    {
      y = 2. / x;
      return (exp(z - x) / sqrt(x / z)) *
             (1.25331414 +
              y * (0.23498619 + y * (-0.3655620e-1 +
                                     y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
    }
  }

  real_t K2exp(real_t x, real_t z)
  /* calculates the extended Modified Bessel functions of second type and of order 2 of x and z */
  /* K2exp(x,z)=K2(x)*exp(z)*sqrt(z) */ { return K0exp(x, z) + 2. / x * K1exp(x, z); }

  real_t K0exp_numdenom(real_t num_x, real_t num_z, real_t denom)
  /* calculates the extended Modified Bessel functions of second type and of order 0 of x and z */
  /* K0exp(x,z)=K0(x)*exp(z)*sqrt(z) with x=num_x/denom, z =num_z/denom*/
  {
    real_t y;
    const real_t x = num_x / denom;
    const real_t z = num_z / denom;

    if (std::isnormal(x) && x <= 2.)
    {
      y = x * x / 4.;
      return ((-log(x / 2.) * I0(x)) +
              (-0.57721566 +
               y * (0.42278420 +
                    y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5))))))) *
             exp(z) * sqrt(z);
    }
    else
    {
      y = 2. / num_x * denom;
      return (std::exp((num_z - num_x) / denom) * sqrt(num_z / num_x)) *
             (1.25331414 +
              y * (-0.7832358e-1 +
                   y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3))))));
    }
  }

  real_t K1exp_numdenom(real_t num_x, real_t num_z, real_t denom)
  /* calculates the extended Modified Bessel functions of second type and of order 1 of x and z */
  /* K1exp(x,z)=K1(x)*exp(z)*sqrt(z) with x=num_x/denom, z =num_z/denom*/
  {
    real_t y;

    const real_t x = num_x / denom;
    const real_t z = num_z / denom;

    if (std::isnormal(x) && x <= 2.)
    {
      y = x * x / 4.;
      return ((log(x / 2.) * I1(x)) +
              (1. / x) *
                  (1. + y * (0.15443144 +
                             y * (-0.67278579 +
                                  y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))))) *
             exp(z) * sqrt(z);
    }
    else
    {
      y = 2. / num_x * denom;
      return (std::exp((num_z - num_x) / denom) * sqrt(num_z / num_x)) *
             (1.25331414 +
              y * (0.23498619 + y * (-0.3655620e-1 +
                                     y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
    }
  }

  real_t K2exp_numdenom(real_t num_x, real_t num_z, real_t denom)
  /* calculates the extended Modified Bessel functions of second type and of order 2 of x and z */
  /* K2exp(x,z)=K2(x)*exp(z)*sqrt(z) with x=num_x/denom, z =num_z/denom*/
  {
    return K0exp_numdenom(num_x, num_z, denom) + 2. / num_x * denom * K1exp_numdenom(num_x, num_z, denom);
  }

  real_t besselj0(real_t x)
  /* Bessel function J0 */
  {
    real_t P1 = 1.e0;
    real_t P2 = -.1098628627e-2;
    real_t P3 = .2734510407e-4;
    real_t P4 = -.2073370639e-5;
    real_t P5 = .2093887211e-6;

    real_t Q1 = -.1562499995e-1;
    real_t Q2 = .1430488765e-3;
    real_t Q3 = -.6911147651e-5;
    real_t Q4 = .7621095161e-6;
    real_t Q5 = -.934945152e-7;

    real_t R1 = 57568490574.e0;
    real_t R2 = -13362590354.e0;
    real_t R3 = 651619640.7e0;
    real_t R4 = -11214424.18e0;
    real_t R5 = 77392.33017e0;
    real_t R6 = -184.9052456e0;

    real_t S1 = 57568490411.e0;
    real_t S2 = 1029532985.e0;
    real_t S3 = 9494680.718e0;
    real_t S4 = 59272.64853e0;
    real_t S5 = 267.8532712e0;
    real_t S6 = 1.e0;
    real_t ax, xx, y, z, result;
    ax = (x < 0) ? -x : x;

    if (ax < 8.0)
    {
      y = x * x;
      result = (R1 + y * (R2 + y * (R3 + y * (R4 + y * (R5 + y * R6))))) /
               (S1 + y * (S2 + y * (S3 + y * (S4 + y * (S5 + y * S6)))));
      return result;
    }
    else
    {
      z = 8. / ax;
      y = z * z;
      xx = ax - .785398164;
      result = sqrt(.636619772 / ax) * (cos(xx) * (P1 + y * (P2 + y * (P3 + y * (P4 + y * P5)))) -
                                        z * sin(xx) * (Q1 + y * (Q2 + y * (Q3 + y * (Q4 + y * Q5)))));
      return result;
    }
  }
  real_t besselj1(real_t x)
  /* Bessel function J1 */
  {
    real_t P1 = 1.e0;
    real_t P2 = .183105e-2;
    real_t P3 = -.3516396496e-4;
    real_t P4 = .2457520174e-5;
    real_t P5 = -.240337019e-6;

    real_t Q1 = .04687499995e0;
    real_t Q2 = -.2002690873e-3;
    real_t Q3 = .8449199096e-5;
    real_t Q4 = -.88228987e-6;
    real_t Q5 = .105787412e-6;

    real_t R1 = 72362614232.e0;
    real_t R2 = -7895059235.e0;
    real_t R3 = 242396853.1e0;
    real_t R4 = -2972611.439e0;
    real_t R5 = 15704.48260e0;
    real_t R6 = -30.16036606e0;

    real_t S1 = 144725228442.e0;
    real_t S2 = 2300535178.e0;
    real_t S3 = 18583304.74e0;
    real_t S4 = 99447.43394e0;
    real_t S5 = 376.9991397e0;
    real_t S6 = 1.e0;
    real_t ax, xx, y, z, signe, resultat;
    ax = (x < 0) ? -x : x;
    signe = (x < 0) ? -1. : 1.;

    if (ax < 8.0)
    {
      y = x * x;
      resultat = x * (R1 + y * (R2 + y * (R3 + y * (R4 + y * (R5 + y * R6))))) /
                 (S1 + y * (S2 + y * (S3 + y * (S4 + y * (S5 + y * S6)))));
      return resultat;
    }
    else
    {
      z = 8. / ax;
      y = z * z;
      xx = ax - 2.356194491;
      resultat = sqrt(.636619772 / ax) * (cos(xx) * (P1 + y * (P2 + y * (P3 + y * (P4 + y * P5)))) -
                                          z * sin(xx) * (Q1 + y * (Q2 + y * (Q3 + y * (Q4 + y * Q5)))));
      return (signe * resultat);
    }
  }
} // end of namespace advmath
