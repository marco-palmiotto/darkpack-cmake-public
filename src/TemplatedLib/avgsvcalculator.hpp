#pragma once

#ifndef __cplusplus
  #define __cplusplus
#endif

#include <list>

#include "advmath.hpp"
#include "correspondance.hpp"
#include "process.hpp"
#include <future>
#include <gsl/gsl_deriv.h>

// #define DEBUG

namespace __SPEC_LIB_NAME__
{
  constexpr const short unsigned int NlowTsv =
      10; //!< Maximum order of the Taylor expansion for \f$<\sigma v>\f$ at low T
  constexpr const std::array<std::array<real_t, NlowTsv + 1>, NlowTsv + 1> CoefflowTsv = {
      {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {-3, 3. / 2, 0, 0, 0, 0, 0, 0, 0, 0, 0},
       {6, -3, 15. / 8, 0, 0, 0, 0, 0, 0, 0, 0},
       {-(75. / 8), 75. / 16, -(15. / 16), 35. / 16, 0, 0, 0, 0, 0, 0, 0},
       {23445. / 2048, -(1485. / 256), -(1575. / 256), -(525. / 64), 315. / 128, 0, 0, 0, 0, 0, 0},
       {-(17505. / 2048), 19395. / 4096, 11925. / 512, 9975. / 512, -(4725. / 512), 693. / 256, 0, 0, 0, 0, 0},
       {-(222885. / 32768), 13095. / 8192, -(878175. / 16384), -(74025. / 2048), 89775. / 4096, -(10395. / 1024),
        3003. / 1024, 0, 0, 0, 0},
       {1661715. / 32768, -(1264815. / 65536), 3173175. / 32768, 1800225. / 32768, -(666225. / 16384), 197505. / 8192,
        -(45045. / 4096), 6435. / 2048, 0, 0, 0},
       {-(1379496825. / 8388608), 32645025. / 524288, -(76137975. / 524288), -(8594775. / 131072), 16202025. / 262144,
        -(1465695. / 32768), 855855. / 32768, -(96525. / 8192), 109395. / 32768, 0, 0},
       {-(13671950879025. / 4294967296), -(2855855475. / 16777216), 186553125. / 1048576, 45808875. / 1048576,
        -(77352975. / 1048576), 35644455. / 524288, -(6351345. / 131072), 1833975. / 65536, -(1640925. / 131072),
        230945. / 65536, 0},
       {-(38822473644075. / 8589934592), -(43047242435475. / 8589934592), -(10584016875. / 67108864),
        275065875. / 4194304, 412279875. / 8388608, -(170176545. / 2097152), 154459305. / 2097152,
        -(13610025. / 262144), 31177575. / 1048576, -(3464175. / 262144),
        969969. / 262144}}}; //!< Coefficients of the Taylor expansion for \f$<\sigma v>\f$ at low T

  constexpr const unsigned short int NlowTsv_split =
      4; //!< Maximum order for the Taylor expansion for \f$<\sigma v>\f$ at low T for small mass splitting
  constexpr const real_t Qnij[NlowTsv_split + 1][NlowTsv_split + 1][NlowTsv_split + 1] = {
      {{1., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}},
      {{0.75, 1., 0., 0., 0.}, {1.5, 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}},
      {{-0.09375, 0.75, 1., 0., 0.},
       {2.625, 1.5, 0., 0., 0.},
       {1.875, 0., 0., 0., 0.},
       {0., 0., 0., 0., 0.},
       {0., 0., 0., 0., 0.}},
      {{0.117188, -0.09375, 0.75, 1., 0.},
       {1.17188, 2.625, 1.5, 0., 0.},
       {6.09375, 1.875, 0., 0., 0.},
       {2.1875, 0., 0., 0., 0.},
       {0., 0., 0., 0., 0.}},
      {{0., 0.117188, -0.09375, 0.75, 1.},
       {0., 1.17188, 2.625, 1.5, 0.},
       {0., 6.09375, 1.875, 0., 0.},
       {0., 2.1875, 0., 0., 0.},
       {2.46094, 0., 0., 0., 0.}}}; //!< $Q_{n,i,j}$ coefficients for the Taylor expansion for \f$<\sigma v>\f$ at low T
                                    //!< for small mass splitting

  // Using the notation phitilde_n = beta_n + eta*lambda_n
  constexpr const std::array<real_t, NlowTsv_split + 1> sv_beta = {
      1., -3.75, 8.90625, -16.5234, 25.1147}; //!< $\beta_n$ coefficients for the Taylor expansion for \f$<\sigma v>\f$
                                              //!< at low T for small mass splitting
  constexpr const std::array<real_t, NlowTsv_split + 1> sv_lambda = {
      -2., 2.75, 5.15625, -41.3672, 145.085}; //!< $\lambda_n$ coefficients for the Taylor expansion for \f$<\sigma
                                              //!< v>\f$ at low T for small mass splitting

  class AvgSvCalculator
  {
protected:
    Param_t input;                                   //!< Input parameters
    RunningSM run;                                   //!< Runnig class
    std::shared_ptr<std::vector<Process2to2>> p_ptr; //!< Pointer to the list of processes
    csl::InitSanitizer<bool> Weffcuts;               //!< Flag to enable/disable the kinematical cuts on Weff

    std::vector<real_t> sqrtStab, //!< Vector of the centre-of-mass energies for each process
        pefftab,                  //!< Vector of the \f$p_{eff}\f$ for each process
        g2_Wefftab;               //!< Vector of the \f$g^2 W_{eff}\f$ for each process

    csl::InitSanitizer<real_t> T_lim_sigmav; //!< Threshold temperature for the calculation of \f$<\sigma v>\f$ with the
                                             //!< algorithm for high temperatures
    short int orderT;
    std::array<int, corr::SIZEPHYSICALBSM> bsm_particles_sorted;
    std::vector<real_t> ytab;
    std::vector<real_t> Weffderiv;
    std::vector<real_t> Weffderiv_err;
    std::vector<real_t> TaylorCoeffSVT; // Taylor coefficients of \f$<\sigma v>\f$ (T) at low T
    bool wastherenosplitting, werealldeltaxinull;


    bool runningenabled;

    void sort_bsm_spectrum();
    void init_vector(std::shared_ptr<std::vector<Process2to2>>&);

public:
    AvgSvCalculator(const Param_t&);
    AvgSvCalculator(const Param_t&, std::shared_ptr<std::vector<Process2to2>>);

    void changeInput(const Param_t& input_ext, bool hasMassSpectrumChanged = true);

    virtual ~AvgSvCalculator()
    {
#ifdef DEBUG
      std::cout << "Called ~AvgSvCalculator()\n";
#endif
    };

    AvgSvCalculator(const AvgSvCalculator&);
    AvgSvCalculator& operator=(const AvgSvCalculator&);

    bool operator==(const AvgSvCalculator&) const;
    bool operator!=(const AvgSvCalculator& other) const { return !(*this == other); };

    void print_procs(std::ostream& out = std::cout, const bool verbose = false) const;
    void print_g2_wefftable(std::ostream& out = std::cout) const;
    void print_coefficients(std::ostream& out = std::cout) const;
    virtual void print(std::ostream& out = std::cout) const
    {
      print_procs(out);
      print_g2_wefftable(out);
      print_coefficients(out);
    };

    /**
     * @brief Prints the processes kinematically allowed at the energy Ecm
     *
     * @param Ecm Centre-of-mass energy
     * @param out Stream where to print the processes
     */
    inline void print_contributing_processes(const real_t& Ecm, std::ostream& out = std::cout) const
    {
      for (const auto& proc : *p_ptr)
      {
        if (proc.isAllowed(input, Ecm))
        {
          out << proc.getMname() << '\n';
        }
      }
    };

    inline void runAtScale(const real_t& Q)
    {
      if (runningenabled)
      {
        run.HandleParamRunning(input, Q);
#ifdef DEBUG
        std::cout << "Running is enabled\n";
#endif
      }
#ifdef DEBUG
      else
        std::cout << "Running is disabled\n";
#endif
    };

    real_t get_g2_Weff(const real_t& sqrtS);
    real_t get_g2_dWeff_dcos(const real_t& sqrtS, const real_t& cosine);

    real_t getWeff(const real_t& sqrtS);
    real_t getdWeff_dcos(const real_t& sqrtS, const real_t& cosine);

    real_t getdWeff_dcos_fromW12(const real_t& sqrtS, const real_t& cosine);
    real_t getWeff_fromW12(const real_t& sqrtS);
    // Getters

    inline size_t getN() const { return p_ptr->size(); };
    inline Process2to2 getProcess(const size_t index) const { return (*p_ptr)[index]; };
    inline Param_t getInput() const { return input; };
    inline std::shared_ptr<std::vector<Process2to2>> getProcList() const { return p_ptr; }
    inline int getLSP() const { return input.getLightestBSMpart(); };
    inline real_t getMassLBSM() const { return input.getLightestBSMmass(); };
    inline real_t getMass(const int& part) const { return input.masses_vector[part]; };

    inline size_t getWefftabsize() const { return g2_Wefftab.size(); };

    inline void get_g2_WeffTabElement(const size_t n, real_t& peff, real_t& sqrts, real_t& g2_weff) const
    {
      peff = pefftab.at(n);
      sqrts = sqrtStab.at(n);
      g2_weff = g2_Wefftab.at(n);
    };

    inline void setOrderT(const short int n) { orderT = n; };

    inline real_t get_TaylorCoeffSVT(const int n) const { return TaylorCoeffSVT.at(n); }

    void setLimitingValues();

    void setTaylorCoeffSVx_nosplitting();
    void setTaylorCoeffSVx_splitting(const real_t deltaxi[], const int Mrho, const int Meta);

    // Setters
    void setWeffcuts(bool x);
    void setProcessList(std::shared_ptr<std::vector<Process2to2>>);
    inline void setT_lim_sigmav(const real_t T) { T_lim_sigmav = T; };

    inline void enablerunning(const bool x = true) { runningenabled = x; };

public:
    void tabulateValues(const real_t& sqrtSmax, const size_t& Nmax);
    real_t getAverageSigmav_coan_noopti(const real_t& T);

    // Implementing here the SuperIso algorithm for the MSSM
    // It works well at any value of T because it performs a "constant" estrapolation at low T
    real_t getAverageSigmav_coan_siso(const real_t& T);

    // The following is the model-independent function
    // Supposed to work well at low temperatures as well
    // It relies on the following methods
    real_t getAverageSigmav_coan(const real_t& T);
    real_t getAverageSigmav_coan_hightemp(const real_t& T);
    real_t getAverageSigmav_coan_lowtemp(const real_t& T);
    real_t getAverageSigmav_coan_lowtemp_degenerate_mass(const real_t& T);

public:
    // Barebone interface function that we define here once for all
    inline real_t getAverageSigmav(const real_t& T) { return getAverageSigmav_coan_hightemp(T); };

  }; // End of AvgSvCalculator

} // End of namespace __SPEC_LIB_NAME__

inline std::ostream& operator<<(std::ostream& out, const __SPEC_LIB_NAME__::AvgSvCalculator& var)
{
  var.print(out);
  return out;
}

#ifdef DEBUG
  #undef DEBUG
#endif
