#pragma once

#ifndef __cplusplus
  #define __cplusplus
#endif

#include "RunningSM.hpp"
#include "advmath.hpp"
#include "correspondance.hpp"
#include "gsl/gsl_sf_psi.h"

// #define DEBUG_COPY
// #define DEBUG
namespace __SPEC_LIB_NAME__
{
  struct Insertion
  {
    int field;
    bool part;

    Insertion(const int i) : field(i), part(true) {}

    Insertion(const int i, const bool b) : field(i), part(b) {}
  };

  class Process2to2
  {
    // Data of the class
private:
    csl::InitSanitizer<int> p[4];   //!< Array of particle fields in the process.
    csl::InitSanitizer<bool> ap[4]; //!< Array of boolean flags indicating if each particle is an antiparticle.

    csl::InitSanitizer<std::string> Key; //!< Key matching this process in the hash table.
    csl::InitSanitizer<bool> Exists;     //!< Flag indicating if the process exists in the library.

    csl::InitSanitizer<Cfptr_t> sumSquaredAmpl{
        "sumSquaredAmpl_Cfptr"};                //!< Function pointer to the sum of squared amplitudes.
    csl::InitSanitizer<short int> Sf34{"Sf34"}; //!< Symmetry factor for the final state.
    csl::InitSanitizer<short int> CPfac{"CPfac"};
    csl::InitSanitizer<short int> combinFac{"dof"};

    RunningSM* runptr;
    csl::InitSanitizer<bool> isRunDataExternal{"isRunDataExternal"}; // Basic
    csl::InitSanitizer<bool> isRunningExternal{"isRunningExternal"}; // Basic
    csl::InitSanitizer<bool> haveToFreerunptr{"haveToFreerunptr"};   // Basic

    std::array<csl::InitSanitizer<int>, 4> part_g{csl::InitSanitizer<int>("g_1"), csl::InitSanitizer<int>("g_2"),
                                                  csl::InitSanitizer<int>("g_3"), csl::InitSanitizer<int>("g_4")};

    // Private methods needed by the constructor or the setters
    std::string findKey();
    void handleSetup();

    // Private methods needed in computations
    bool EP13calculation(const Param_t& input, const real_t& sqrts, real_t& p1, real_t& p3, real_t& E1,
                         real_t& E3) const;
    bool setKinematics(Param_t& input, const real_t& sqrts, const real_t& ctheta, real_t& p1, real_t& p3,
                       real_t sij[5][5]) const;

    // Constructors
public:
    // Empty constructor, destructor
    Process2to2();
    ~Process2to2();

    /**
      @brief Constructor for Process2to2.
      @param particles Array of 4 integers representing particle fields.
      @param antiparticles Array of 4 booleans indicating if each particle is an antiparticle.
    */
    Process2to2(const std::array<int, 4>&, const std::array<bool, 4>&);

    Process2to2(const std::array<Insertion, 4>&);

    /**
      @brief Set the n-th particle in the process.
      @param n Particle index (1-4).
      @param ip Particle field.
      @param iap True if antiparticle, false otherwise.
      @return Short integer status code.
    */
    short int set(short int n, const int& ip, const bool iap);

    // Copy constructor and operators
    Process2to2(const Process2to2& other);
    bool operator==(const Process2to2&) const;
    bool operator!=(const Process2to2&) const;
    Process2to2& operator=(const Process2to2&);

    /** @brief Checks if the process exists in the hash table.
        @return True if the process exists, false otherwise.
    */
    inline bool checkExistance() const { return Exists; };

    inline int getField(const size_t i) const { return p[i]; };

    inline bool getMatter(const size_t i) const { return ap[i]; };

    std::string getKey() const;
    std::string getName() const;
    std::string getMname() const;
    inline real_t getMass(const short int i, const Param_t& input) const { return input.masses_vector[p[i - 1]]; };
    inline short int getSf34() const { return Sf34; };
    inline short int getDof() const { return combinFac; };

    inline void print(std::ostream& out = std::cout) const
    {
      out << getMname() << "\nSf_34 = " << getSf34() << " dof = " << getDof() << " Cpfac " << CPfac.get() << '\n'
          << part_g[0] << part_g[1] << part_g[2] << part_g[3] << std::endl;
    }

    inline Cfptr_t getSumSquaredAmpl_ptr() const { return sumSquaredAmpl; };

    // setters
    inline bool isComplete() const
    {
      for (int i = 0; i < 4; i++)
        if (p[i] < 0)
          return false;
      return true;
    };

    inline void setRunningData(RunningSM* runin)
    {
      runptr = runin;
      isRunDataExternal = true;
    };
    inline void setRunningExternal() { isRunningExternal = true; };
    inline void setRunningInternal() { isRunningExternal = false; };
    inline void setRunDataInternal(Param_t& input)
    {
      runptr = new RunningSM(input);
      isRunDataExternal = false;
    };


    // Methods to perform computations
    inline void handleRunning(Param_t& input, const real_t& Ecm)
    {
      if (runptr == nullptr)
      {
        runptr = new RunningSM(input);
        haveToFreerunptr = true;
      }
      runptr->HandleParamRunning(input, Ecm);
#ifdef DEBUG
      std::cout << "parameter running handled: input is \n";
      input.Print();
#endif
    };

    /**
     * @brief Returns the sum of the sum of the squared amplitudes for the instance of the class.
     * @param input numerical parameters
     * @param sqrts centre-of-mass energy
     * @param ctheta cosine of the angle
     * @return real_t \f$\sum |\mathcal{M}|^2\f$ for the instance of the class
     */
    real_t getSumSquaredAmpl(Param_t& input, const real_t& sqrts, const real_t& ctheta);

    /**
     * @brief Returns the sum of the sum of the squared amplitudes, avereaged with the combinatorial factor, for the
     * instance of the class.
     * @param input numerical parameters
     * @param sqrts centre-of-mass energy
     * @param ctheta Cosine of the angle between particle 1 and 3.
     * @return real_t \f$\sum |\bar \mathcal{M}|^2\f$  for the instance of the class
     */
    inline real_t getAvgSquaredAmpl(Param_t& input, const real_t& sqrts, const real_t& ctheta)
    {
      return getSumSquaredAmpl(input, sqrts, ctheta) / combinFac;
    };

    /**
      @brief  Computes the process-dependent contributions to dWeff/dcostheta
              I.e.:
                (8 PI) g^2_LBSM * sqrtS * peff *dWeff/dcostheta = p12*p34/Sf34*CPfac*(sum |M|^2)
              for the given inputs.
              Running is assumed to be handled externally
      @param input Param_t object containing numerical inputs.
      @param sqrts Centre of mass energy.
      @param ctheta Cosine of the angle between particle 1 and 3.
      @return real_t
     */
    real_t get_process_dependent_contrib_dweff_dcos(Param_t& input, const real_t& sqrts, const real_t& ctheta) const;

    /**
     * @brief  This function computes [ g^2_LBSM * dWeff/dcostheta] for the given inputs
     *
      @param input Param_t object containing numerical inputs.
      @param sqrts Centre of mass energy.
      @param ctheta Cosine of the angle between particle 1 and 3.
      @return Contribution to dWeff/d(cos(theta)).
    */
    real_t get_g2_dweff_dcos(Param_t& input, const real_t& sqrts, const real_t& ctheta);

    /**
      @brief Get the contribution to dWeff/d(cos(theta)) for the process.
             We divide by (g_DM)^2 the result of the method get_g2_dweff_dcos
      @param input Param_t object containing numerical inputs.
      @param sqrts Centre of mass energy.
      @param ctheta Cosine of the angle between particle 1 and 3.
      @return Contribution to dWeff/d(cos(theta)).
    */
    real_t getDiffWeffContrib(Param_t& input, const real_t& Ecm, const real_t& ctheta);

    real_t getDiffCrossSection(Param_t& input, const real_t& sqrts, const real_t& ctheta);
    real_t getTotalCrossSection(Param_t& input, const real_t& sqrts, real_t* discr = nullptr);

    /**
      @brief Function defined for testing purposes. Get the Diff W12 Contrib object.
            It returns
            g_1 g_2 p1 W_12->34,
            using the identity
            g_1 g_2 p1 W_12->34 = (g_1 g_2 p_1)*4* sqrts p1 * dsigma_12->34 / dcostheta
      @param input Param_t object containing numerical inputs.
      @param sqrts Centre of mass energy.
      @param ctheta Cosine of the angle between particle 1 and 3.
      @param wantWeffPrefac
      @return real_t
     */
    real_t getDiffW12Contrib(Param_t& input, const real_t& sqrts, const real_t& ctheta,
                             const bool wantWeffPrefac = false);
    real_t getTotalW12Contrib(Param_t& input, const real_t& sqrts, real_t* discr = nullptr);

    /**
     * @brief Computes the relative speed in the initial state
     *
     * @param input Param_t object containing numerical inputs.
     * @param sqrts Centre of mass energy.
     * @return real_t Value of the relative speed in the intial state
     */
    real_t getInitialSpeed(Param_t& input, const real_t& sqrts) const;

    /**
     * @brief Computes the relative speed in the final state
     *
     * @param input Param_t object containing numerical inputs.
     * @param sqrts Centre of mass energy.
     * @return real_t Value of the relative speed in the final state
     */
    real_t getFinalSpeed(Param_t& input, const real_t& sqrts) const;

    // Friend classes and functions
    friend real_t getContributionSingle(const Process2to2 p, Param_t input, const real_t Ecm, const real_t ctheta);
    friend real_t getContribution(const std::shared_ptr<std::vector<Process2to2>>& vp, Param_t input, const real_t Ecm,
                                  const real_t ctheta, const bool cut);
  };

  inline std::ostream& operator<<(std::ostream& out, const Process2to2& proc)
  {
    proc.print(out);
    return out;
  }

} // end of namespace __SPEC_LIB_NAME__
