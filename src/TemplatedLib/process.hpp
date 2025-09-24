#pragma once

#ifndef PROCESS__SPEC_LIB_NAME__
  #define PROCESS__SPEC_LIB_NAME__

  #ifndef __cplusplus
    #define __cplusplus
  #endif

  #include "RunningSM.hpp"
  #include "config.hpp"
  #include "correspondance.hpp"


// #define DEBUG_COPY
// #define DEBUG
namespace __SPEC_LIB_NAME__
{
  /**
   * @brief This structure is defined to have an alternative wa of defining a particle in a process.
   *
   */
  struct Insertion
  {
    int field; //!< The field of the particle
    bool part; //!< True if the particle is a particle, false if it is an antiparticle

    /**
     * @brief Construct a new Insertion object, supposing that field excitation is a particle, not an antiparticle.
     *
     * @param i The integer corresponding to the field of the particle
     */
    Insertion(const int i) : field(i), part(true) {}

    /**
     * @brief Construct a new Insertion object
     *
     * @param i The integer corresponding to the field of the particle
     * @param b True if the particle is a particle, false if it is an antiparticle
     */
    Insertion(const int i, const bool b) : field(i), part(b) {}
  };

  class Process2to2
  {
    // Data of the class
private:
    std::array<csl::InitSanitizer<int>, 4> p; //!< Array of particle fields in the process.
    std::array<csl::InitSanitizer<bool>, 4>
        ap; //!< Array of boolean flags indicating if each particle is an antiparticle.

    csl::InitSanitizer<std::string> Key; //!< Key matching this process in the hash table.
    csl::InitSanitizer<bool> Exists;     //!< Flag indicating if the process exists in the library.

    csl::InitSanitizer<Cfptr_t> sumSquaredAmpl{
        "sumSquaredAmpl_Cfptr"};                    //!< Function pointer to the sum of squared amplitudes.
    csl::InitSanitizer<short int> Sf34{"Sf34"};     //!< Symmetry factor for the final state.
    csl::InitSanitizer<short int> CPfac{"CPfac"};   //!< CP symmetry factor for the reaction
    csl::InitSanitizer<short int> combinFac{"dof"}; //!< Combinatorial factor for the reaction

    RunningSM*
        runptr; //!< Pointer to the RunningSM instace to be used for the running in the current instance of Process2to2
    csl::InitSanitizer<bool> isRunDataExternal{
        "isRunDataExternal"}; //!< Flag indicating if the RunningSM data are external
    csl::InitSanitizer<bool> isRunningExternal{
        "isRunningExternal"}; //!< Flag indicating if the RunningSM instance is owned by something other than the
                              //!< current instance of Process2to2
    csl::InitSanitizer<bool> haveToFreerunptr{
        "haveToFreerunptr"}; //!< Flag indicating if the RunningSM instance has to be freed

    std::array<csl::InitSanitizer<int>, 4> part_g{
        csl::InitSanitizer<int>("g_1"), csl::InitSanitizer<int>("g_2"), csl::InitSanitizer<int>("g_3"),
        csl::InitSanitizer<int>("g_4")}; //!< Array containing the degrees of freedom of each particle in the process.


    // Private methods needed by the constructor or the setters
    /**
     * @brief Find the key for the process in the hash table.
     *
     * @return The key for the process in the hash table.
     */
    std::string findKey();

    /**
     * @brief Performs the setup of the process.
     */
    void handleSetup();

    // Private methods needed in computations
    /**
     * @brief Computes the energy and the momentum (in absoluite value) of the particles 1 and 3 in the center of mass
     * frame.
     *
     * @param input The numerical parameters.
     * @param sqrts The centre-of-mass energy.
     * @param p1 The momentum of particle 1.
     * @param p3 The momentum of particle 3.
     * @param E1 The energy of particle 1.
     * @param E3 The energy of particle 3.
     * @return true If the process is kinematically allowed and the calculation was successful.
     * @return false Otherwise.
     */
    bool EP13calculation(const Param_t& input, const real_t& sqrts, real_t& p1, real_t& p3, real_t& E1,
                         real_t& E3) const;

    /**
     * @brief Set the kinematical variables for the process.
     *
     * @param input The numerical parameters.
     * @param sqrts The centre-of-mass energy.
     * @param ctheta The cosine of the angle between particles 1 and 3 in the centre-of-mass frame.
     * @param p1 The momentum of particle 1.
     * @param p3 The momentum of particle 3.
     * @param sigdfj Matrix containing the scalar products of the momenta of the particles in the process, in the for
     * s_ij = p_i * p_j.
     * @return true If the process is kinematically allowed and the calculation was successful.
     * @return false Otherwise.
     */
    bool setKinematics(Param_t& input, const real_t& sqrts, const real_t& ctheta, real_t& p1, real_t& p3,
                       real_t sij[5][5]) const;

    /**
     * @brief Set the kinematical variables for the process.
     *
     * @param input The numerical parameters.
     * @param s The value of the Mandelstam variable s.
     * @param t The value of the Mandelstam variable t.
     * @param p1 The momentum of particle 1.
     * @param p3 The momentum of particle 3.
     * @param sigdfj Matrix containing the scalar products of the momenta of the particles in the process, in the for
     * s_ij = p_i * p_j.
     * @return true If the process is kinematically allowed and the calculation was successful.
     * @return false Otherwise.
     */
    bool set_kinematics_from_t(Param_t& input, const real_t& s, const real_t& t, real_t& p1, real_t& p3,
                               real_t sij[5][5]) const;

    // Constructors
public:
    /**
     * @brief Empty constructor for the Process2to2 class.
     *
     */
    Process2to2();

    /**
     * @brief Destructor for the Process2to2 class.
     *
     */
    ~Process2to2();

    /**
     * @brief Construct a new Process 2to 2 object
     *
     * @param field The fields of the particles in the process.
     * @param isParticle True if the particle is a particle, false if it is an antiparticle.
     */
    Process2to2(const std::array<int, 4>& field, const std::array<bool, 4>& isParticle);

    /**
     * @brief Construct a new Process 2to 2 object from an array of Insertion objects.
     *
     * @param insertArray
     */
    Process2to2(const std::array<Insertion, 4>& insertArray);

    /**
      @brief Set the n-th particle in the process.
      @param n Particle index (1-4).
      @param ip Particle field.
      @param iap True if antiparticle, false otherwise.
      @return Short integer status code.
    */
    short int set(short int n, const int& ip, const bool iap);

    // Copy constructor and operators
    /**
     * @brief Copy-constructor for the Process2to2 class.
     *
     * @param other
     */
    Process2to2(const Process2to2& other);

    /**
     * @brief Comparison operator for the Process2to2 class.
     *
     * @return true if the processes have the same key
     * @return false otherwise
     */
    bool operator==(const Process2to2&) const;

    /**
     * @brief Comparison operator for the Process2to2 class.
     *
     * @return true if the processes have different keys
     * @return false otherwise
     */
    bool operator!=(const Process2to2&) const;

    /**
     * @brief Assignment operator for the Process2to2 class. Constructs a new Process2to2 object, with the same key
     *        of the argument given.
     *
     * @param other
     * @return Process2to2&
     */
    Process2to2& operator=(const Process2to2& other);

    /** @brief Checks if the process exists in the hash table.
        @return True if the process exists, false otherwise.
    */
    inline bool checkExistance() const { return Exists; };

    /**
     * @brief Get the Field object corresponding to the i-th particle in the process.
     *
     * @param i Number of the field (1-4).
     * @return The enumeration corresponding to the field of the i-th particle in the process.
     */
    inline int getField(const size_t i) const { return p[i]; };

    /**
     * @brief Returns true if the i-th particle is a particle, false if it is an antiparticle.
     *
     * @param i Number of the field (1-4).
     * @return true if the i-th particle is a particle
     * @return false otherwise
     */
    inline bool getMatter(const size_t i) const { return ap[i]; };

    /**
     * @brief Returns the key of the process in the hash table.
     *
     * @return The key of the process in the hash table.
     */
    std::string getKey() const;

    /**
     * @brief Returns the name of the process.
     *
     * @return The name of the process.
     */
    std::string getName() const;

    /**
     * @brief Same as getName().
     */
    std::string getMname() const;

    /**
     * @brief Returns the mass of the i-th particle in the process.
     *
     * @param i  Number of the field (1-4).
     * @param input Param_t object containing numerical inputs.
     * @return real_t The value of the mass of the i-th particle in the process.
     */
    inline real_t getMass(const short int i, const Param_t& input) const { return input.masses_vector[p[i - 1]]; };

    /**
     * @brief Get the symmetry factor for the final state.
     *
     * @return short int the symmetry factor for the final state.
     */
    inline short int getSf34() const { return Sf34; };

    /**
     * @brief Get the degrees of freedom for the process.
     *
     * @return short int the degrees of freedom for the process.
     */
    inline short int getDof() const { return combinFac; };

    /**
     * @brief Returns the minimum energy for producing the final state.
     *
     * @param input Param_t object containing numerical inputs.
     * @return real_t The sum of the masses of the final state particles.
     */
    inline real_t compute_final_state_treshold(const Param_t& input) const
    {
      return input.masses_vector[p[2]] + input.masses_vector[p[3]];
    };

    /**
     * @brief Prints the name of the process and its specifications to the given output stream.
     *
     * @param out Output stream to print to (default is std::cout).
     */
    inline void print(std::ostream& out = std::cout) const
    {
      out << getMname() << "\nSf_34 = " << getSf34() << " dof = " << getDof() << " Cpfac " << CPfac.get() << '\n'
          << part_g[0] << part_g[1] << part_g[2] << part_g[3] << std::endl;
    }

    /**
     * @brief Returns true if the process is kinematically allowed at a certain centre-of-mass energy, otherwise returns
     * false.
     *
     * @param input Param_t object containing numerical inputs.
     * @param sqrts centre-of-mass energy
     * @return true
     * @return false
     */
    inline bool isAllowedAtDefinedEcm(const Param_t& input, const real_t& sqrts) const
    {
      return ((sqrts < getMass(3, input) + getMass(4, input)) || (sqrts < getMass(1, input) + getMass(2, input)));
    };

    /**
     * @brief Returns true if the mass of the final state is larger or equal than the mass in the initial state.
     *
     * @param input Param_t object containing numerical inputs.
     * @return true
     * @return false
     */
    inline bool isAllowedAtZeroMomentum(const Param_t& input) const
    {
      return (getMass(3, input) + getMass(4, input)) <= (getMass(1, input) + getMass(2, input));
    };

    /**
     * @brief Get the function pointer to the sum of squared amplitudes.
     *
     * @return Cfptr_t the function pointer to the sum of squared amplitudes.
     */
    inline Cfptr_t getSumSquaredAmpl_ptr() const { return sumSquaredAmpl; };

    // setters

    /**
     * @brief Returns true if the process is complete, i.e. all the particles are defined.
     *
     */
    inline bool isComplete() const
    {
      for (int i = 0; i < 4; i++)
        if (p[i] < 0)
          return false;
      return true;
    };

    /**
     * @brief Enables external running, by taking the pointer of the RunningSM instance as input.
     *
     * @param runin the pointer of the external RunningSM instance to be used for running
     */
    inline void setRunningData(RunningSM* runin)
    {
      runptr = runin;
      isRunDataExternal = true;
    };

    /**
     * @brief Sets to true the flag indicating that the running is external.
     *
     */
    inline void setRunningExternal() { isRunningExternal = true; };

    /**
     * @brief Sets to false the flag indicating that the running is external.
     *
     */
    inline void setRunningInternal() { isRunningExternal = false; };

    /**
     * @brief Sets the running as internal, by creating a new instance of RunningSM.
     *
     * @param input The parameters needed to create the RunningSM instance.
     */
    inline void setRunDataInternal(Param_t& input)
    {
      runptr = new RunningSM(input);
      isRunDataExternal = false;
    };


    // Methods to perform computations

    /**
     * @brief Handles the running of the parameters.
     *
     * @param input Param_t object containing numerical inputs.
     * @param Ecm The centre-of-mass energy.
     */
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
     * @param input Param_t object containing numerical inputs.
     * @param sqrts centre-of-mass energy
     * @param ctheta cosine of the angle
     * @return real_t \f$\sum |\mathcal{M}|^2\f$ for the instance of the class
     */
    real_t getSumSquaredAmpl(Param_t& input, const real_t& sqrts, const real_t& ctheta);

    /**
     * @brief Returns the sum of the sum of the squared amplitudes, avereaged with the combinatorial factor, for the
     * instance of the class.
     * @param input Param_t object containing numerical inputs.
     * @param sqrts centre-of-mass energy
     * @param ctheta Cosine of the angle between particle 1 and 3.
     * @return real_t \f$\sum |\bar M|^2\f$  for the instance of the class
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
      @return real_t Contribution to dWeff/d(cos(theta)).
    */
    real_t get_g2_dweff_dcos(Param_t& input, const real_t& sqrts, const real_t& ctheta);

    /**
      @brief Get the contribution to dWeff/d(cos(theta)) for the process.
             We divide by (g_DM)^2 the result of the method get_g2_dweff_dcos
      @param input Param_t object containing numerical inputs.
      @param sqrts Centre of mass energy.
      @param ctheta Cosine of the angle between particle 1 and 3.
      @return real_t Contribution to dWeff/d(cos(theta)).
    */
    real_t getDiffWeffContrib(Param_t& input, const real_t& Ecm, const real_t& ctheta);

    /**
       @brief Computes the differential cross-section at the given centre-of-mass energy Ecm and \f$cos(\theta_{13})\f$.

       @param input Param_t object containing numerical inputs.
       @param sqrts Centre of mass energy.
       @param ctheta Cosine of the angle between particle 1 and 3.
       @return real_t Differential cross-section.
     */
    real_t getDiffCrossSection(Param_t& input, const real_t& sqrts, const real_t& ctheta);


    /**
     @brief Computes the differential cross-section at the given centre-of-mass energy Ecm and \f$t\f$.

     @param input Param_t object containing numerical inputs.
     @param sqrts Centre of mass energy.
     @param t Value of the Mandelstam variable t.
     @return real_t Differential cross-section.
   */
    real_t compute_dsigma_dt(Param_t& input, const real_t& sqrts, const real_t& t);

    /**
       @brief Computes the total cross-section at the given centre-of-mass energy Ecm.

       @param input Param_t object containing numerical inputs.
       @param sqrts Centre of mass energy.
       @param ctheta Cosine of the angle between particle 1 and 3.
       @return real_t Differential cross-section.
     */
    real_t getTotalCrossSection(Param_t& input, const real_t& sqrts, real_t* discr = nullptr);

    /**
      @brief Function defined for testing purposes. Get the Diff W12 Contrib object.
            It returns
            $$ g_1 g_2 p_1 W_{12\to 34} $$
            using the identity
            $$ g_1 g_2 p_1 W_{12\to 34} = 4(g_1 g_2 p_1) \sqrt{s} p_1 \frac{d \sigma_{12\to 34}}{d \cos(\theta)} $$

      @param input Param_t object containing numerical inputs.
      @param sqrts Centre of mass energy.
      @param ctheta Cosine of the angle between particle 1 and 3.
      @param wantWeffPrefac Optional boolean to indicate if the prefactor of the effective potential is wanted
      (default=false).
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
  }; // End of the class Process2to2

  /**
   * @brief Overloading of the << operator on an output stram reference, to print the Process2to2 object.
   *
   * @param out The reference to the output stream.
   * @param proc The process to be printed.
   * @return std::ostream& The modified output stream.
   */
  inline std::ostream& operator<<(std::ostream& out, const Process2to2& proc)
  {
    proc.print(out);
    return out;
  }

} // end of namespace __SPEC_LIB_NAME__

#endif // PROCESS__SPEC_LIB_NAME__