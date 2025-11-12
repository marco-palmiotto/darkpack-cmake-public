#include "propagation.hpp"


template <typename T> static inline T SQUARE(const T x) { return x * x; }

template <typename T> static inline T SQRT(const T x) { return pow(x, 0.5); }

namespace __SPEC_LIB_NAME__
{

  Propagation_param_t::Propagation_param_t(const Param_t& input, enum propagation_models p_model,
                                           enum halo_profiles halo_model)
          : ind_param(input)
  {
    // Fetch the dark matter masses from the PPPC4DMID data
    ind_param.fetch_pppc4dmid_data();
    ind_param.read_fermi_data();

    // Fetch the cross section data for antiproton production
    ind_param.read_file(FILENAME_H_ON_H, 1, h_on_h_xsection);
    ind_param.read_file(FILENAME_HE_ON_H, 1, he_on_h_xsection);
    ind_param.read_file(FILENAME_HE_ON_HE, 1, he_on_he_xsection);
    ind_param.read_file(FILENAME_H_ON_HE, 1, h_on_he_xsection);

    switch (halo_model)
    {
    case (halo_profiles::NFW):
      halo_profile = &nfw;
      break;
    case (halo_profiles::BURKERT):
      halo_profile = &burkert;
      break;
    case (halo_profiles::EINASTO):
      halo_profile = &einasto;
      break;
    default:
      std::cerr << "Invalid halo model specified. Using NFW profile by default." << std::endl;
      halo_profile = &nfw;
    }

    switch (p_model)
    {
    case (propagation_models::MAX):
      equation_parameters = {0.0016, 0.46, 15.0, 1.0, 0.014};
      break;
    case (propagation_models::MED):
      equation_parameters = {0.0112, 0.7, 4.0, 12.0, 52.9};
      break;
    case (propagation_models::MIN):
      equation_parameters = {0.0016, 0.85, 1.0, 13.5, 22.4};
      break;
    default:
      std::cerr << "Invalid propagation model specified. Using 'med' parameters by default." << std::endl;
      equation_parameters = {0.0112, 0.7, 4.0, 12.0, 52.9};
    }
  };

  real_t Propagation_param_t::max_likelihood_aux_func(const std::vector<real_t>& logJ, const std::vector<real_t>& xtra,
                                                      const std::vector<std::vector<real_t>>& spec)
  {
    // Auxiliary function to be used in the maximization of the likelihood (needed redefinition to make the optimization
    // compatible with both fermi-LAT and AMS-02)
    real_t sigma_j = xtra[0];
    auto dsph = static_cast<int>(xtra[1]);

    real_t logJ0 = logJ[0];

    return ind_param.likelihood_one_dsph(dsph, logJ0);
  }

  real_t Propagation_param_t::max_likelihood_one_dsph(const int& dsph, const real_t& logJ)
  {
    // Calls max_likelihood_aux_func and returns the maximum log likelihood for one dSph by optimizing over logJ
    auto logJ_factors = ind_param.get_logJ();
    real_t logJ_obs = logJ_factors.at(dsph).at(0);
    real_t sigma_j = logJ_factors.at(dsph).at(1);

    std::vector<real_t> logJ_opt;
    logJ_opt.reserve(1);
    real_t lhmax;

    std::vector<real_t> xtra;
    xtra.reserve(2);
    xtra[0] = sigma_j;
    xtra[1] = dsph;

    std::vector<real_t> xlim_min = {1.0e-3};
    std::vector<real_t> xlim_max = {logJ_obs + 10. * sigma_j};

    std::vector<std::vector<real_t>> spect;

    IDpowell(1, &Propagation_param_t::max_likelihood_aux_func, xtra, spect, xlim_min, xlim_max, &lhmax, logJ_opt,
             1.0e-3, "max");

    return lhmax;
  }

  int Propagation_param_t::IDpowell(int n, input_function_t func, const std::vector<real_t>& xtra,
                                    const std::vector<std::vector<real_t>>& spect, const std::vector<real_t>& xlim_min,
                                    const std::vector<real_t>& xlim_max, real_t* fmin, std::vector<real_t>& xmin,
                                    real_t ftol, const std::string& option)
  /* finds the minimum (option="min") or maximum (option="max") with a precions ftol of the function func(real_t x[],
   * const std::array<real_t>& xtra) having n free parameters x[] with a lower bound xlim_min and an upper bound
   * xlim_max. The minimum (maximum) value of the function is saved in fmin and the corresponding values of the free
   * parameters in xmin. xtra[] is an array containing the fixed parameters of function func. */
  {
    int nmax = 5;

    srand((unsigned int)time(NULL));
    std::vector<real_t> x0;
    x0.reserve(nmax);
    int i;

    real_t factor = 0.;
    if (option == "min")
      factor = 1.;
    else if (option == "max")
      factor = -1.;

    for (i = 0; i < n; i++)
      x0[i] = (xlim_min[i] + xlim_max[i]) / 2.;
    int test = IDpowellaux(n, func, xtra, spect, x0, xlim_min, xlim_max, fmin, xmin, ftol, factor);
    if (option != "max")
      *fmin = -*fmin;

    return test;
  }

  /*----------------------------------------------------------------*/

  int Propagation_param_t::IDpowellaux(int n, input_function_t func, const std::vector<real_t>& xtra,
                                       const std::vector<std::vector<real_t>>& spect, const std::vector<real_t>& x0,
                                       const std::vector<real_t>& xlim_min, const std::vector<real_t>& xlim_max,
                                       real_t* fmin, std::vector<real_t>& xmin, real_t ftol, real_t factor)
  /* finds the minimum (factor=1) or maximum (factor=-1) using a Powell method with a precions ftol of the function
   * func(real_t x[], const std::array<real_t>& xtra) having n free parameters x[] with a lower bound xlim_min and an
   * upper bound xlim_max. The minimum value of the function is saved in fmin and the corresponding values of the free
   * parameters in xmin. xtra[] is an array containing the fixed parameters of function func. The computation starts
   * with the values x0[] of the free parameters x[]. */
  {
    real_t TINY = 1.0e-70;
    int ITMAX = 500;
    int nmax = 5;

    int i, j, iter;
    std::vector<std::vector<real_t>> vect; // set of vectors
    vect.reserve(nmax);

    /* initializes set of vectors */
    for (i = 0; i < n; i++)
    {
      vect[i].reserve(nmax);
      for (j = 0; j < n; j++)
      {
        if (i == j)
          vect[i][j] = 1;
        else
          vect[i][j] = 0;
      }
    }

    for (j = 0; j < n; j++)
      xmin[j] = x0[j]; // initializes minimum position

    std::vector<real_t> vecti; // current vector in the iteration
    vecti.reserve(nmax);
    real_t fmintemp;
    int imax;
    real_t fdeltamax;
    *fmin = factor * (this->*func)(x0, xtra, spect);
    real_t fiter;
    std::vector<real_t> xiter;
    xiter.reserve(nmax);
    int testbrent;

    std::vector<real_t> meanvect, xextrapol;
    meanvect.reserve(nmax);
    xextrapol.reserve(nmax);

    real_t fextrapol;
    real_t t;

    for (iter = 0; iter <= ITMAX; iter++)
    {

      fiter = *fmin;
      for (j = 0; j < n; j++)
        xiter[j] = xmin[j];
      fdeltamax = 0;
      for (i = 0; i < n; i++)
      {
        for (j = 0; j < n; j++)
          vecti[j] = vect[j][i];
        fmintemp = *fmin;
        testbrent = IDbrentmethod(n, func, xtra, spect, xmin, vecti, xlim_min, xlim_max, fmin, factor);
        if (testbrent == 0)
        {
#ifdef DEBUG
          std::cout << "Brent method unsuccessful\n";
#endif
          return 0;
        }

        if ((fmintemp - *fmin) > fdeltamax)
        {
          imax = i;
          fdeltamax = fmintemp - *fmin;
        }
      }


      if (2. * (fiter - (*fmin)) <= ftol * (fabs(fiter) + fabs(*fmin)) + TINY)
        return 0;

      if (iter == ITMAX)
      {
#ifdef DEBUG
        std::cout << "Powell exceeding maximum iterations.\n";
#endif
        return 0;
      }

      for (j = 0; j < n; j++)
      {
        // check if it is worth to keep the mean deplacement vector
        meanvect[j] = xmin[j] - xiter[j];
        xextrapol[j] = xmin[j] + meanvect[j];
      }

      fextrapol = factor * (this->*func)(xextrapol, xtra, spect);
      if (fextrapol < fiter)
      {
        t = 2. * (fiter - 2. * (*fmin) + fextrapol) * pow(fiter - (*fmin) - fdeltamax, 2.) -
            fdeltamax * pow(fiter - fextrapol, 2.);

        if (t < 0.)
        {
          IDbrentmethod(n, func, xtra, spect, xmin, meanvect, xlim_min, xlim_max, fmin, factor);
          vect[j][imax] = vect[j][n];
          vect[j][n] = meanvect[j];
        }
      }
    }

    return 1;
  }


  /*--------------------------------------------------------------------*/


  int Propagation_param_t::IDbrentmethod(int n, input_function_t func, const std::vector<real_t>& xtra,
                                         const std::vector<std::vector<real_t>>& spect, std::vector<real_t>& xmin,
                                         const std::vector<real_t>& xi, const std::vector<real_t>& xlim_min,
                                         const std::vector<real_t>& xlim_max, real_t* fmin, real_t factor)
  /* minimizes function func along vector xi */
  {
    real_t a0, b0, x0, xmin1d;
    std::vector<real_t> xitemp;
    xitemp.reserve(n);
    int i;

    IDbraket(n, func, xtra, spect, &a0, &x0, &b0, xmin, xi, xlim_min, xlim_max, factor);

    int test = IDbrentmethod1D(func, xtra, spect, a0, x0, b0, fmin, &xmin1d, factor, xmin, xi, n);
    for (i = 0; i < n; i++)
      xmin[i] = xmin[i] + xmin1d * xi[i];

    return test;
  }

  /*--------------------------------------------------------------------*/

  real_t Propagation_param_t::IDbrentmethod1Dfunc(input_function_t func, const std::vector<real_t>& xtra,
                                                  const std::vector<std::vector<real_t>>& spect, real_t x,
                                                  const std::vector<real_t>& xmini, const std::vector<real_t>& xi,
                                                  int n)
  /* function f(x)=func(xmini+x*xi, xtra) */
  {
    int i;
    std::vector<real_t> xitemp;
    xitemp.reserve(n);
    for (i = 0; i < n; i++)
      xitemp[i] = xmini[i] + x * xi[i];

    return (this->*func)(xitemp, xtra, spect);
  }

  int Propagation_param_t::IDbrentmethod1D(input_function_t func, const std::vector<real_t>& xtra,
                                           const std::vector<std::vector<real_t>>& spect, real_t a0, real_t x0,
                                           real_t b0, real_t* fmin, real_t* xmin, real_t factor,
                                           std::vector<real_t> xmini, std::vector<real_t> xi, int n)
  /* finds the minimum of DDbrentmethod1Dfunc for x between a0 and b0, starting with an initial value x=x0 */
  {
    int brent_iter_max = 1000;
    real_t ZEPS = 1.0e-70;
    real_t CGOLD = 0.1;
    real_t tol = 1.0e-3;

    real_t tol1, tol2;
    real_t a, b;

    if (a0 < b0)
    {
      a = a0;
      b = b0;
    }
    else
    {
      a = b0;
      b = a0;
    }

    real_t x, w, v; // points used in the parabolic fit
    real_t u;
    real_t pa1, pa2; // parabola parameters pa2*X^2+pa1*X+pa0
    real_t pamin;    // parabola minimum
    real_t xm;
    real_t fx, fw, fv, fu;
    real_t e = 0, etmp, d;
    int iter = 0;

    x = w = v = x0;
    fx = fw = fv = factor * IDbrentmethod1Dfunc(func, xtra, spect, x0, xmini, xi, n);
    xm = 0.5 * (a + b);

    while (iter < brent_iter_max && fabs(x - xm) + 0.5 * (b - a) > 2 * tol1)
    {
      tol1 = tol * fabs(x) + ZEPS;
      tol2 = e * tol1;

      if (fabs(e) > tol1)
      {
        pa2 = ((fw - fx) / (w - x) - (fv - fx) / (v - x)) / (w - v); // parabolic fit
        pa1 = (fv - fx) / (v - x) - pa2 * (v + x);
        pamin = -0.5 * pa1 / pa2;
        etmp = e;
        e = d;
        if (fabs(x - pamin) > 0.5 * etmp || pamin < a || pamin > b)
        {
          // cannot use the parabolic fit
          if (x >= xm)
            e = a - x;
          else
            e = b - x;
          d = CGOLD * e;
        }
        else
        {
          d = pamin - x;
          u = pamin;
          if ((u - a) < tol2 || (b - u) < tol2)
          {
            // pamin too close to borders
            d = fabs(tol1);
            if (xm < x)
              d = -d;
          }
        }
      }
      else
      {
        if (x >= xm)
          e = a - x;
        else
          e = b - x;
        d = CGOLD * e;
      }

      if (fabs(d) >= tol1)
        u = x + d;
      else
      {
        if (d < 0)
          u = x - fabs(tol1);
        else
          u = x + fabs(tol1);
      }
      fu = factor * IDbrentmethod1Dfunc(func, xtra, spect, u, xmini, xi, n);

      if (fu <= fx)
      {
        if (u >= x)
          a = x;
        else
          b = x;
        v = w;
        w = x;
        x = u;
        fv = fw;
        fw = fx;
        fx = fu;
      }
      else
      {
        if (u < x)
          a = u;
        else
          b = u;
        if (fu <= fw || w == x)
        {
          v = w;
          w = u;
          fv = fw;
          fw = fu;
        }
        else if (fu <= fv || v == x || v == w)
        {
          v = u;
          fv = fu;
        }
      }

      xm = 0.5 * (a + b);

      iter++;
    }


    *xmin = x;
    *fmin = fx;

    if (iter == brent_iter_max)
    {
      std::cout << "WARNING: too many iterations in Brent's method, precision decreased\n";
      return 0;
    }
    else
      return 1;
  }

  /*--------------------------------------------------------------------*/

  void Propagation_param_t::IDbraket(int n, input_function_t func, const std::vector<real_t>& xtra,
                                     const std::vector<std::vector<real_t>>& spect, real_t* a0, real_t* x0, real_t* b0,
                                     std::vector<real_t> xinit, std::vector<real_t> xi,
                                     const std::vector<real_t>& xlim_min, const std::vector<real_t>& xlim_max,
                                     real_t factor)
  /* find values  of a0, b0 and x0 required in function IDbrentmethod1D */
  {
    int BRAKETMAX = 50;
    int nmax = 5;

    int i;
    real_t a0tmp, b0tmp;

    real_t xr, frand;
    std::vector<real_t> xrand, xa, xb;
    xrand.reserve(nmax);
    xa.reserve(nmax);
    xb.reserve(nmax);
    real_t fa, fb, fmin;

    int afirst = 1;
    int bfirst = 1;
    for (i = 0; i < n; i++)
    {
      if (xi[i] != 0)
      {
        a0tmp = (xlim_min[i] - xinit[i]) / xi[i];
        b0tmp = (xlim_max[i] - xinit[i]) / xi[i];

        if (a0tmp > *a0 || afirst)
        {
          *a0 = a0tmp;
          afirst = 0;
        }
        if (b0tmp < *b0 || bfirst)
        {
          *b0 = b0tmp;
          bfirst = 0;
        }
      }
    }

    if (*a0 > *b0)
    {
      *a0 = 0;
      *b0 = 0;
      *x0 = 0;
    }
    else
    {
      for (i = 0; i < n; i++)
      {
        xa[i] = xinit[i] + *a0 * xi[i];
        xb[i] = xinit[i] + *b0 * xi[i];
      }
      fa = factor * (this->*func)(xa, xtra, spect);
      fb = factor * (this->*func)(xb, xtra, spect);
      if (fa < fb)
      {
        *x0 = *a0;
        fmin = fa;
      }
      else
      {
        *x0 = *b0;
        fmin = fb;
      }
      int irand;

      int icount = 0;
      for (irand = 0; irand < BRAKETMAX; irand++)
      {
        xr = *a0 + irand * (*b0 - *a0) / BRAKETMAX;
        for (i = 0; i < n; i++)
          xrand[i] = xinit[i] + xr * xi[i];
        frand = factor * (this->*func)(xrand, xtra, spect);
        if (frand < fmin)
        {
          *x0 = xr;
          fmin = frand;
          icount++;
        }
        if (icount == 5)
          break;
      }
    }
  }

  // Propagation_param_t::~Propagation_param_t() = default;
} // namespace __SPEC_LIB_NAME__