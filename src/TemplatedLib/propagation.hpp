#pragma once

#ifndef PROPAGATION__SPEC_LIB_NAME__
  #define PROPAGATION__SPEC_LIB_NAME__

  #ifndef __cplusplus
    #define __cplusplus
  #endif

  #include "RunningSM.hpp"
  #include "advmath.hpp"
  #include "config.hpp"
  #include "correspondance.hpp"
  #include "indirect.hpp"
  #include "process.hpp"
  #include "process_1to2.hpp"
  #include <array>
  #include <fstream>
  #include <iostream>
  #include <sstream>
  #include <vector>

namespace __SPEC_LIB_NAME__
{

  class Propagation_param_t //< Class containing the parameters and methods for indirect detection
  {
    using input_function_t = real_t (Propagation_param_t::*)(const std::vector<real_t>&, const std::vector<real_t>&,
                                                             const std::vector<std::vector<real_t>>&);
    using density_profile_t = real_t (*)(real_t, real_t);

private:
    std::array<real_t, 5> equation_parameters; //!< Parameters for the propagation equation
    real_t r_earth;                            //!< Distance between the Earth and the galactic center
    real_t rho_chi_local;                      //!< Local density of dark matter
    density_profile_t halo_profile;            //!< Pointer to the dark matter halo profile function
    Indirectparam_t ind_param;                 //!< Indirect detection parameters, for storing the dark matter masses
    std::vector<std::vector<real_t>>
        h_on_h_xsection; //!< Table of differential cross sections for p + p -> pbar + X processes
    std::vector<std::vector<real_t>>
        he_on_h_xsection; //!< Table of differential cross sections for He + p -> pbar + X processes
    std::vector<std::vector<real_t>>
        he_on_he_xsection; //!< Table of differential cross sections for He + He -> pbar + X processes
    std::vector<std::vector<real_t>>
        h_on_he_xsection; //!< Table of differential cross sections for p + He -> pbar + X processes


public:
    enum halo_profiles : short int
    {
      NFW = 1,
      BURKERT = 2,
      EINASTO = 3
    };
    enum propagation_models : short int
    {
      MIN = 1,
      MED = 2,
      MAX = 3
    };
    /**
     * @brief Construct a new Propagation_param_t object
     *
     * @param input Param_t object containing numerical inputs.
     *
     * @param propagation_model String that defines the set of values for the equation parameters. It can be set to
     * "min", "med" or "max".
     *
     * @param halo_model Integer that defines the dark matter halo profile. It can be set to 1 (NFW), 2 (Burkert) or 3
     * (Einasto).
     *
     */
    Propagation_param_t(const Param_t& input, enum propagation_models p_model, enum halo_profiles halo_model);

    Propagation_param_t() = default;

    /**
     * @brief Destructor for the Propagation_param_t class.
     *
     */
    ~Propagation_param_t();

    static inline real_t nfw(real_t r, real_t r_s)
    {
      real_t rho_s = 0.32; // in GeV/cm^3
      return rho_s / ((r / r_s) * (1 + r / r_s) * (1 + r / r_s));
    };
    static inline real_t burkert(real_t r, real_t r_s)
    {
      real_t rho_s = 1.57; // in GeV/cm^3
      return rho_s / (1 + r / r_s) / (1 + pow(r / r_s, 2.0));
    };
    static inline real_t einasto(real_t r, real_t r_s)
    {
      real_t alpha = 0.22;
      real_t rho_s = 0.11; // in GeV/cm^3
      return rho_s * exp(-2.0 / alpha * (pow(r / r_s, alpha) - 1.0));
    };
    /**
     * @brief Obtains the maximum likelihood for a given dSph by varying the J-factor.
     *
     * @param logJ The log J-factor.
     *
     * @param sigma_j The uncertainty on the log J-factor.
     *
     */
    real_t max_likelihood_aux_func(const std::vector<real_t>& logJ, const std::vector<real_t>& xtra,
                                   const std::vector<std::vector<real_t>>& spec);

    real_t max_likelihood_one_dsph(const int& dsph, const real_t& logJ);

    //--------------Optimization Methods------------------//


    int IDpowell(int n, input_function_t func, const std::vector<real_t>& xtra,
                 const std::vector<std::vector<real_t>>& spect, const std::vector<real_t>& xlim_min,
                 const std::vector<real_t>& xlim_max, real_t* fmin, std::vector<real_t>& xmin, real_t ftol,
                 const std::string& option);

    int IDpowellaux(int n, input_function_t func, const std::vector<real_t>& xtra,
                    const std::vector<std::vector<real_t>>& spect, const std::vector<real_t>& x0,
                    const std::vector<real_t>& xlim_min, const std::vector<real_t>& xlim_max, real_t* fmin,
                    std::vector<real_t>& xmin, real_t ftol, real_t factor);


    int IDbrentmethod(int n, input_function_t func, const std::vector<real_t>& xtra,
                      const std::vector<std::vector<real_t>>& spect, std::vector<real_t>& xmin,
                      const std::vector<real_t>& xi, const std::vector<real_t>& xlim_min,
                      const std::vector<real_t>& xlim_max, real_t* fmin, real_t factor);

    real_t IDbrentmethod1Dfunc(input_function_t func, const std::vector<real_t>& xtra,
                               const std::vector<std::vector<real_t>>& spect, real_t x,
                               const std::vector<real_t>& xmini, const std::vector<real_t>& xi, int n);



    int IDbrentmethod1D(input_function_t func, const std::vector<real_t>& xtra,
                        const std::vector<std::vector<real_t>>& spect, real_t a0, real_t x0, real_t b0, real_t* fmin,
                        real_t* xmin, real_t factor, std::vector<real_t> xmini, std::vector<real_t> xi, int n);



    void IDbraket(int n, input_function_t func, const std::vector<real_t>& xtra,
                  const std::vector<std::vector<real_t>>& spect, real_t* a0, real_t* x0, real_t* b0,
                  std::vector<real_t> xinit, std::vector<real_t> xi, const std::vector<real_t>& xlim_min,
                  const std::vector<real_t>& xlim_max, real_t factor);

  }; // class Propagation_param_t
} // namespace __SPEC_LIB_NAME__

#endif // PROPAGATION__SPEC_LIB_NAME__