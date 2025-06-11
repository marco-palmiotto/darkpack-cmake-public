#pragma once

#ifndef __cplusplus
  #define __cplusplus
#endif

#ifndef RUNNINGSM__SPEC_LIB_NAME__
  #define RUNNINGSM__SPEC_LIB_NAME__

  #include "config.hpp"
  #include "correspondance.hpp"
  #include "params_new.hpp"

namespace __SPEC_LIB_NAME__
{
  /**
    @struct RunningSM
    @brief Structure for handling the running of Standard Model parameters.
    This struvture contains methods and data members for calculating running coupling constants
    and particle masses in the Standard Model.
  */
  struct RunningSM
  {
public:
    /**
      @enum ParticlesList
      @brief Enumeration of Standard Model particles.
    */
    enum ParticlesList : short int
    {
      UP = 0,
      DOWN,
      STRANGE,
      CHARM,
      BEAUTY,
      TOP,
      EL,
      MU,
      TAU,
      NUE,
      NUMU,
      NUTAU,
      GLUON,
      W,
      Z,
      PHOTON,
      HIGGS
    };

private:
    std::array<std::string, 17> ParticlesNames = {"UP",    "DOWN", "STRANGE", "CHARM",  "BEAUTY", "TOP",
                                                  "EL",    "MU",   "TAU",     "NUE",    "NUMU",   "NUTAU",
                                                  "GLUON", "W",    "Z",       "PHOTON", "HIGGS"};

    // Parameters needed for the running
    csl::InitSanitizer<real_t> alphas_MZ{
        "alphas_MZ"}; //!< Parameter needed for running's computation: \f$\alpha_s(M_Z)\f$. LHA block SMINPUT, N=3
    csl::InitSanitizer<real_t> mass_Z{
        "mass_Z"}; //!< Parameter needed for running's computation: \f$M_Z\f$. LHA block SMINPUT, N=4

    csl::InitSanitizer<real_t> m_d2GeV{
        "m_d2GeV"}; //!< Parameter needed for running's computation: \f$m_d(2\mathrm{GeV})\f$. LHA block SMINPUT, 21
    csl::InitSanitizer<real_t> m_u2GeV{
        "m_u2GeV"}; //!< Parameter needed for running's computation: \f$m_u(2\mathrm{GeV})\f$. LHA block SMINPUT, 22
    csl::InitSanitizer<real_t> m_s2GeV{
        "m_s2GeV"}; //!< Parameter needed for running's computation: \f$m_s(2\mathrm{GeV})\f$. LHA block SMINPUT, 23

    csl::InitSanitizer<real_t> mc_mc{
        "mc_mc"}; //!< Parameter needed for running's computation: \f$m_c(m_c)\f$. LHA block SMINPUT, 24

    csl::InitSanitizer<real_t> mb_mb{
        "mb_mb"}; //!< Parameter needed for running's computation: \f$m_b(m_b)\f$. LHA block SMINPUT, 5

    csl::InitSanitizer<real_t> mass_top_pole{"mass_top_pole"};
    /**< @brief Parameter needed for running's computation: \f$m_{t,\mathrm{pole}}\f$. LHA block SMINPUT, 6
     */

    bool runlightquarks; //!< Boolean parameter to determine if the light quark masses will or not be ran
    bool runcharm;       //!< Boolean parameter to determine if the charm quark mass will or not be ran
    bool higgsloops;     //!< Boolean parameter to determine if the Higgs loop factor correction is enabled or not

    csl::InitSanitizer<real_t> mass_b_pole{"mass_b_pole"}; //!< Member to store the pole mass of the b-quark
    csl::InitSanitizer<real_t> mass_c_pole{"mass_c_pole"}; //!< Member to store the pole mass of the c-quark

    csl::InitSanitizer<real_t> QuarkMass[6];
    /**< @brief Input parameters: array containg the values of the input quark masses.
     *
     *   The top quark mass is assigned with GetMtopMtop().
     */

    csl::InitSanitizer<real_t> QuarkMassQ[6];
    /**< @brief Input parameters: array containg the energy scales corresponding to the values of the input quark masses
     */

    csl::InitSanitizer<real_t>
        LastQuarkMass[6]; //!< Member for optimization: array of the last computed values of the quark masses
    csl::InitSanitizer<real_t> LastQuarkMassQ[6];
    /**< @brief Member for optimization: array of the energyn scales of the last computed values of the quark masses
     */

    csl::InitSanitizer<real_t> alphas{"alphas"}; //!< Member for optimization: last computed value of \f$\alpha_s\f$
    csl::InitSanitizer<real_t> Q_alphas{
        "Q_alphas"}; //!< Member for optimization: energy scale of the last computed value of \f$\alpha_s\f$

    csl::InitSanitizer<real_t> Lambda3{"Lambda3"};
    csl::InitSanitizer<real_t> alphasMZ_Lambda3{"alphasMZ_Lambda3"};

    csl::InitSanitizer<real_t> Lambda4{"Lambda4"};
    csl::InitSanitizer<real_t> alphasMZ_Lambda4{"alphasMZ_Lambda4"};

    csl::InitSanitizer<real_t> Lambda5{"Lambda5"};
    csl::InitSanitizer<real_t> alphasMZ_Lambda5{"alphasMZ_Lambda5"};

    csl::InitSanitizer<real_t> Lambda6{"Lambda6"};
    csl::InitSanitizer<real_t> alphasMZ_Lambda6{"alphasMZ_Lambda6"};

public:
    /**
      @brief Calculates the strong coupling constant at a given energy scale.
      @param Q The energy scale.
      @param mtop_pole Optional top quark pole mass (default is -1.0, i.e. the corresponding class member value is
      used).
      @return The value of alpha_strong at the given scale.
    */
    real_t AlphaStrong(real_t Q, real_t mtop_pole = -1.0);

    /**
     * @brief Computes the quantity \f$m_t(m_t)\f$
     *
     * @return real_t \f$m_t(m_t)\f$
     */
    real_t GetMtopMtop();

    /**
     * @brief Set the names for all the elements of the arrays of csl::InitSanitizer<> among the members of the class
     */
    constexpr inline void setSanitNames()
    {
      constexpr const unsigned short int n_arrays = 4, n_quarks = 6;
      const char* listofprefixes[n_arrays] = {"Mass_", "Q_", "LastMass_", "LastQ_"};
      csl::InitSanitizer<real_t>* listofarrays[n_arrays] = {QuarkMass, QuarkMassQ, LastQuarkMass, LastQuarkMassQ};
      for (unsigned short int i = 0; i < n_arrays; i++)
      {
        for (unsigned short int j = 0; j < n_quarks; j++)
          listofarrays[i][j].setName(static_cast<std::string>(listofprefixes[i]) + ParticlesNames[j]);
      }
    };

    /**
     * @brief Initializes the RunningSM object with input parameters.
     * @param input Reference to a Param_t object containing initial parameter values.
     */
    void init(const Param_t& input);

    /**
      @brief Constructor that initializes the RunningSM object with input parameters.
      @param input Reference to a Param_t object containing initial parameter values.
    */
    RunningSM(const Param_t& input)
    {
      setSanitNames();
      init(input);
    };

    /**
      @brief Default constructor that initializes the RunningSM object with PDG values.
      @note It would be nice to constexpr this constructor
    */
    RunningSM();

    ~RunningSM()
    {
#ifdef DEBUG
      std::cout << "Called ~RunningSM\n";
#endif
    };

    /**
      @brief Equality comparison operator for deep comparison of RunningSM objects.
      @param other The RunningSM object to compare with.
      @return true if the objects are equal, false otherwise.
    */
    bool operator==(const RunningSM&) const;

    /**
      @brief Inequality comparison operator for RunningSM objects.
      @param other The RunningSM object to compare with.
      @return true if the objects are not equal, false otherwise.
    */
    bool operator!=(const RunningSM& other) const { return !(*this == other); };

    /**
     * @brief Prints all the data members
     * @param out the stream where to print (default std::cout)
     */
    void print(std::ostream& out = std::cout) const;

    /**
     * @brief This function allows to enable or disable the running of the light quarks (u d s) masses.
     *
     * If the running of the light quarks' masses is enabled, also the running of the charm is enabled.
     *
     * @param light true to enable, false to disable
     */
    inline void RunLightQuarks(const bool light = true)
    {
      runlightquarks = light;
      runcharm = runlightquarks ? true : runcharm;
    };

    /**
     * @brief This function allows to enable or disable the running of the charm quark mass.
     *
     * If the running of the charm mass is disabled, also the running of its lighter quarks is disabled.
     *
     * @param runvalue true to enable, false to disable
     */
    inline void RunCharmMass(const bool runvalue = true)
    {
      runcharm = runvalue;
      runlightquarks = runcharm ? runlightquarks : false;
    };

    /**
     * @brief Enables or disables the Higgs loop factor correction in the running of quark masses
     * @param x true to enable, false to disable
     * @see Formula (2.11) from the micromega 1.3 manual
     */
    inline void setHiggsLoopFactor(const bool x = true) { higgsloops = x; };

    /**
     * @brief Returns the value of the top pole mass, stored as a data member of the instance.
     *
     * @return real_t this->mass_top_pole
     */
    inline real_t GetTopPoleMass() const { return mass_top_pole; };

    /**
     * @brief Returns the value of \f$m_b(m_b)\f$ stored as a data member of the instance.
     *
     * @return real_t this->mb_mb
     */
    inline real_t GetMbMb() const { return mb_mb; };

    /**
     * @brief Computes the running quark mass at the energy Qfin, from a given running quark mass quark_mass at energy
     * Qinit
     *
     * @param quark_mass Input quark mass at the scale Qinit
     * @param Qinit Energy scale of the input mass
     * @param Qfin Energy scale of the output mass
     * @return real_t The quark mass at the scale Qfin
     */
    real_t RunQuarkMass(real_t quark_mass, const real_t& Qinit, const real_t& Qfin);

    /**
     * @brief Computes the pole mass of the c-quark at 2-loops
     *
     * @return real_t The value of the pole mass of the c-quark at 2-loops
     */
    real_t GetMcPole();

    /**
     * @brief Computes the pole mass of the b-quark at 2-loops
     *
     * @return real_t The value of the pole mass of the b-quark at 2-loops
     */
    real_t GetMbPole();


    /**
     * @brief Computes the pole mass of the c-quark at 1-loop
     *
     * @return real_t The value of the pole mass of the c-quark at 1-loop
     */
    real_t GetMcPole1Loop();

    /**
     * @brief Computes the pole mass of the b-quark at 1-loop
     *
     * @return real_t The value of the pole mass of the b-quark at 1-loop
     */
    real_t GetMbPole1Loop();

    /**
     * @brief Computes the pole mass of the c-quark at 3-loops
     *
     * @return real_t The value of the pole mass of the c-quark at 3-loops
     */
    real_t GetMcPole3Loops();

    /**
     * @brief computes the 1S b mass
     *
     * @return real_t the 1S b mass
     */
    real_t GetMb1S();

    /**
     * @brief Computes the \f$m_c(m_c)\f$ mass from the pole mass of the c-quark
     *
     * @param mcpole pole mass of the c-quark, as input
     * @param loop number of desired loops of the computation
     * @return real_t \f$m_c(m_c)\f$
     * @note the number of loops actually used can be 1, 2, or 3. The fallback value is 3.
     */
    real_t GetMcMcFromPole(const real_t& mcpole, const int loop);

    /**
     * @brief Computes the \f$m_b(m_b)\f$ mass from the pole mass of the b-quark
     *
     * @param mbpole pole mass of the b-quark, as input
     * @param loop number of desired loops of the computation
     * @return real_t \f$m_b(m_b)\f$
     * @note the number of loops actually used can be 1 or 2. The fallback value is 2.
     */
    real_t GetMbMbFromPole(const real_t& mbpole, const int loop);

    /**
      @brief Calculates the mass of a quark at a given energy scale.
      @param part The quark type from the ParticlesList enum.
      @param Qf The energy scale.
      @return The mass of the quark at the given scale.
    */
    real_t GetQuarkMass(enum ParticlesList part, const real_t& Qf);

    /**
      @brief Updates the input parameters with values at a given energy scale.
      @param input Reference to a Param_t object to be updated.
      @param Q The energy scale.
    */
    void HandleParamRunning(Param_t& input, const real_t& Q);
  };

  void TestRunningMasses(RunningSM* const run, Param_t& input);

} // End of namespace __SPEC_LIB_NAME__

/**
 * @brief Prints all the data members of the __SPEC_LIB_NAME__::RunningSM in input
 * @param out the stream where to print (default std::cout)
 * @param var the __SPEC_LIB_NAME__::RunningSM in input
 * @return the reference at the print stream
 */
inline std::ostream& operator<<(std::ostream& out, const __SPEC_LIB_NAME__::RunningSM& var)
{
  var.print(out);
  return out;
}

#endif // RUNNINGSM__SPEC_LIB_NAME__