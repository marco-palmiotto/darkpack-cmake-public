#pragma once

#ifndef PROCESS1TO2__SPEC_LIB_NAME__
  #define PROCESS1TO2__SPEC_LIB_NAME__

  #ifndef __cplusplus
    #define __cplusplus
  #endif

  #include "RunningSM.hpp"
  #include "advmath.hpp"
  #include "config.hpp"
  #include "correspondance.hpp"
  #include "gsl/gsl_sf_psi.h"
  #include "process.hpp"

// #define DEBUG_COPY
// #define DEBUG
namespace __SPEC_LIB_NAME__
{
  /**
   * @brief Checks whether there are decaying processes for a given particle in corr::squaredampl_1to2
   */
  // bool decays(int particle);

  class Process_1to2
  {
    // Data of the class
private:
    csl::InitSanitizer<int> p[3];   //!< Array of particle fields in the process.
    csl::InitSanitizer<bool> ap[3]; //!< Array of boolean flags indicating if each particle is an antiparticle.

    csl::InitSanitizer<std::string> Key; //!< Key matching this process in the hash table.
    csl::InitSanitizer<bool> Exists;     //!< Flag indicating if the process exists in the library.

    csl::InitSanitizer<CXXfptr_t> sumSquaredAmpl{
        "sumSquaredAmpl_CXXfptr"};                  //!< Function pointer to the sum of squared amplitudes.
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

    std::array<csl::InitSanitizer<int>, 3> part_g{
        csl::InitSanitizer<int>("g_1"),
        csl::InitSanitizer<int>("g_2"),
        csl::InitSanitizer<int>("g_3"),
    }; //!< Array containing the degrees of freedom of each particle in the process.


    // Private methods needed by the constructor or the setters
    /**
     * @brief Find the key for the process in the hash table.
     *
     * @return The key for the process in the hash table.
     */
    std::string find_key();

    /**
     * @brief Performs the setup of the process.
     */
    void handle_setup();

public:
    /**
     * @brief Empty constructor for the Process_1to2 class.
     *
     */
    Process_1to2();

    /**
     * @brief Destructor for the Process_1to2 class.
     *
     */
    ~Process_1to2();

    /**
     * @brief Construct a new Process_1to2 object
     *
     * @param field The fields of the particles in the process.
     * @param isParticle True if the particle is a particle, false if it is an antiparticle.
     */
    Process_1to2(const std::array<int, 3>& field, const std::array<bool, 3>& is_particle);

    /**
     * @brief Construct a new Process_1to2 object from an array of Insertion objects.
     *
     * @param insertArray
     */
    Process_1to2(const std::array<Insertion, 3>& insert_array);

    /**
      @brief Set the n-th particle in the process.
      @param n Particle index (1-3).
      @param ip Particle field.
      @param iap True if antiparticle, false otherwise.
      @return Short integer status code.
    */
    short int set(short int n, const int& ip, const bool iap);

    // Copy constructor and operators
    /**
     * @brief Copy-constructor for the Process_1to2 class.
     *
     * @param other
     */
    Process_1to2(const Process_1to2& other);

    /**
     * @brief Comparison operator for the Process_1to2 class.
     *
     * @return true if the processes have the same key
     * @return false otherwise
     */
    bool operator==(const Process_1to2& other) const;

    /**
     * @brief Comparison operator for the Process_1to2 class.
     *
     * @return true if the processes have different keys
     * @return false otherwise
     */
    bool operator!=(const Process_1to2& other) const;

    /**
     * @brief Assignment operator for the Process_1to2 class. Constructs a new Process2to2 object, with the same key
     *        of the argument given.
     *
     * @param other
     * @return Process_1to2&
     */
    Process_1to2& operator=(const Process_1to2& other);

    /** @brief Checks if the process exists in the hash table.
        @return True if the process exists, false otherwise.
    */
    inline bool check_existance() const { return Exists; };

    /**
     * @brief Returns true if the process is complete, i.e. all the particles are defined.
     *
     */
    inline bool is_complete() const
    {
      for (int i = 0; i < 3; i++)
        if (p[i] < 0)
          return false;
      return true;
    };

    /**
     * @brief Get the Field object corresponding to the i-th particle in the process.
     *
     * @param i Number of the field (1-3).
     * @return The enumeration corresponding to the field of the i-th particle in the process.
     */
    inline int get_field(const size_t i) const { return p[i]; };

    /**
     * @brief Returns true if the i-th particle is a particle, false if it is an antiparticle.
     *
     * @param i Number of the field (1-3).
     * @return true if the i-th particle is a particle
     * @return false otherwise
     */
    inline bool get_matter(const size_t i) const { return ap[i]; };

    /**
     * @brief Returns the key of the process in the hash table.
     *
     * @return The key of the process in the hash table.
     */
    std::string get_key() const;

    /**
     * @brief Returns the name of the process.
     *
     * @return The name of the process.
     */
    std::string get_name() const;

    /**
     * @brief Same as get_name().
     */
    std::string get_mname() const;

    /**
     * @brief Returns the mass of the i-th particle in the process.
     *
     * @param i  Number of the field (1-3).
     * @param input Param_t object containing numerical inputs.
     * @return real_t The value of the mass of the i-th particle in the process.
     */
    inline real_t get_mass(const short int i, const Param_t& input) const { return input.masses_vector[p[i - 1]]; };

    /**
     * @brief Get the symmetry factor for the final state.
     *
     * @return short int the symmetry factor for the final state.
     */
    inline short int get_sf34() const { return Sf34; };

    /**
     * @brief Get the degrees of freedom for the process.
     *
     * @return short int the degrees of freedom for the process.
     */
    inline short int get_dof() const { return combinFac; };

    inline bool is_allowed_at_zero_momentum(const Param_t& input) const
    {
      return (get_mass(2, input) + get_mass(3, input)) <= get_mass(1, input);
    };

    bool set_kinematics(Param_t& input, real_t sij[5][5], real_t m_vec[3], real_t& E1, real_t& E2, real_t& p2) const;

    inline void handle_running(Param_t& input, const real_t& Ecm)
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
     * @brief Calculates the partial width associated with the process.
     *
     * @param input Param_t object containing numerical inputs.
     *
     * @return real_t \f$\Gamma\f$ for the instance of the class.
     */
    real_t get_partial_width(Param_t& input);

    /**
     * @brief Calculates the branching ratio associated to the process.
     *
     * @param input Param_t object containing numerical inputs.
     *
     * @return real_t Branching ratio for the instance of the class.
     */
    real_t get_branching_ratio(Param_t& input);
  }; // end of class Process_1to2

} // end of namespace __SPEC_LIB_NAME__

#endif // PROCESS1TO2__SPEC_LIB_NAME__