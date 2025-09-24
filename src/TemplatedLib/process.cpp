#include "process.hpp"
#include "advmath.hpp"
#include "gsl/gsl_sf_psi.h"

#include "copy_to_c_struct.hpp"

namespace __SPEC_LIB_NAME__
{
// Constant to convert GeV^-2 to pbarn (surfaces)
#ifndef GEVTOPBARN
  #define GEVTOPBARN (0.389379 * 1.0e+9)
#endif

// Constant to convert GeV^-2 to cm^3/s (sigma*v)
#ifndef GEVTOCM3S
  #define GEVTOCM3S 1.1673288750358198e-05
#endif

#define BASIC_INITIALISATION                                                                                           \
  for (int i = 0; i < 4; i++)                                                                                          \
  {                                                                                                                    \
    p[i] = -1;                                                                                                         \
    ap[i] = false;                                                                                                     \
  }                                                                                                                    \
  Exists = false;                                                                                                      \
  sumSquaredAmpl = nullptr;                                                                                            \
  runptr = nullptr;                                                                                                    \
  isRunDataExternal = false;                                                                                           \
  isRunningExternal = false;                                                                                           \
  haveToFreerunptr = false;

#define masses_p(X) (input.masses_vector[p[X]])

  template <typename T> inline void swap(T* x, T* y)
  {
    T temp = *x;
    *x = *y;
    *y = temp;
  }

  template <typename T> static inline T SQUARE(const T x) { return x * x; }

  Process2to2::Process2to2(void){BASIC_INITIALISATION}

  Process2to2::~Process2to2(void)
  {
#ifdef DEBUG_COPY
    std::cout << "Calling ~Process2to2\n";
#endif
    Key.get().clear();
    if (haveToFreerunptr)
      delete runptr;
  }

  std::string Process2to2::findKey()
  {
    std::string key;
    using namespace corr;
    if (!isComplete())
    {
#ifdef DEBUG
      std::cout << "Process is not complete\n";
#endif
      Exists = false;
      key.push_back(EMPTYCHAR);
      return key;
    }
    key.reserve(8);

    for (auto i = 0; i < 4; i++)
    {
      if (!ap[i])
        key.push_back(ANTICHAR);
      key.push_back(EMPTYCHAR + p[i]);
    }
    auto it = squaredampl.find(key);
    if (it != squaredampl.end())
    {
#ifdef DEBUG
      std::cout << "Found key\n";
#endif
      Key = key;
      return key;
    }
    // If the key is not found, permutation have to be used
    key.clear();
    std::array<int, 4> v = {1, 0, 2, 3};
    for (auto i : v)
    {
      if (!ap[i])
        key.push_back(ANTICHAR);
      key.push_back(EMPTYCHAR + p[i]);
    }
    it = squaredampl.find(key);
    if (it != squaredampl.end())
    {
#ifdef DEBUG
      std::cout << "Found key swapping p0 and p1\n";
      std::cout << "Names before swapping are " << getMname() << /* ", " << getSname() << */ std::endl;
#endif
      Key = key;
      swap(&p[0], &p[1]);
      swap(&ap[0], &ap[1]);
#ifdef DEBUG
      std::cout << "Names after swapping are " << getMname() << /* ", " << getSname() << */ std::endl;
      scanf("%*c");
#endif
      return key;
    }

    key.clear();
    v = {0, 1, 3, 2};
    for (auto i : v)
    {
      if (!ap[i])
        key.push_back(ANTICHAR);
      key.push_back(EMPTYCHAR + p[i]);
    }
    it = squaredampl.find(key);
    if (it != squaredampl.end())
    {
#ifdef DEBUG
      std::cout << "Found key swapping p2 and p3\n";
      std::cout << "Names before swapping are " << getMname() << /* ", " << getSname() << */ std::endl;
#endif
      Key = key;
      swap(&p[2], &p[3]);
      swap(&ap[2], &ap[3]);
#ifdef DEBUG
      std::cout << "Names after swapping are " << getMname() << /* ", " << getSname() << */ std::endl;
      scanf("%*c");
#endif
      return key;
    }

    key.clear();
    v = {1, 0, 2, 3};
    for (auto i : v)
    {
      if (!ap[i])
        key.push_back(ANTICHAR);
      key.push_back(EMPTYCHAR + p[i]);
    }
    it = squaredampl.find(key);
    if (it != squaredampl.end())
    {
#ifdef DEBUG
      std::cout << "Found key swapping p0 p1, p2 p3\n";
      std::cout << "Names before swapping are " << getMname() << /* ", " << getSname() << */ std::endl;
#endif
      Key = key;
      swap(&p[0], &p[1]);
      swap(&ap[0], &ap[1]);
      swap(&p[2], &p[3]);
      swap(&ap[2], &ap[3]);
#ifdef DEBUG
      std::cout << "Names after swapping are " << getMname() << /* ", " << getSname() << */ std::endl;
      scanf("%*c");
#endif
      return key;
    }

#ifdef DEBUG
    std::cout << "Key not found: the process does not exist\n";
#endif

    Exists = false;
    key.clear();
    key.push_back(EMPTYCHAR);
    return key;
  }

  void Process2to2::handleSetup()
  {
    if (isComplete())
    {
      Exists = true;
      Key = findKey();

      if (Exists)
      {
#ifdef DEBUG
        std::cout << "Setting Sf34\n";
#endif
        Sf34 = (p[2] == p[3] && ((ap[2] && ap[3]) || (!ap[2] && !ap[3]))) ? 2 : 1;
#ifdef DEBUG
        std::cout << "Assigning values...";
#endif
        auto it = corr::squaredampl.find(Key);
#ifdef DEBUG
        std::cout << "... squaredamp";
#endif
        sumSquaredAmpl = std::get<0>(it->second);
#ifdef DEBUG
        std::cout << "... combfac";
#endif
        combinFac = std::get<1>(it->second);
#ifdef DEBUG
        std::cout << "... cpfac";
#endif
        CPfac = std::get<2>(it->second);

#ifdef DEBUG
        std::cout << "... runptr, extdata, extrunning,";
#endif
        runptr = nullptr;
        isRunDataExternal = false;
        isRunningExternal = false;
#ifdef DEBUG
        std::cout << "... helicity dofs are \n";
#endif

        for (size_t i = 0; i < 4; i++)
        {
          part_g[i] = corr::part_hel_dof[p[i]];
#ifdef DEBUG
          std::cout << part_g[i] << ' ';
#endif
        }

#ifdef DEBUG
        std::cout << "\nExiting handleSetup of " << getName() << "\n";
#endif
        return;
      }
      else
      {
#ifdef DEBUG
        std::cout << "Process does not exist: resetting variables\n";
#endif
        Sf34 = 0;
        sumSquaredAmpl = nullptr;
        combinFac = 0;
        CPfac = 0;

        runptr = nullptr;
        isRunDataExternal = false;
        isRunningExternal = false;
        return;
      }
    }
#ifdef DEBUG
    std::cout << "Process is not complete\n";
#endif
  }

  Process2to2::Process2to2(const std::array<int, 4>& ip, const std::array<bool, 4>& iap)
  {
#ifdef DEBUG
    std::cout << "Constructor called\n";
#endif
    BASIC_INITIALISATION
#ifdef DEBUG
    std::cout << "Pushing back values\n";
#endif
    for (int i = 0; i < 4; i++)
    {
      p[i] = ip[i];
      ap[i] = iap[i];
    }
#ifdef DEBUG
    std::cout << "Data pushed back\n";
#endif
    handleSetup();
  }

  Process2to2::Process2to2(const std::array<Insertion, 4>& v)
  {
#ifdef DEBUG
    std::cout << "Constructor called\n";
#endif
    BASIC_INITIALISATION
#ifdef DEBUG
    std::cout << "Pushing back values\n";
#endif
    for (int i = 0; i < 4; i++)
    {
      p[i] = v[i].field;
      ap[i] = v[i].part;
    }
#ifdef DEBUG
    std::cout << "Pushing back values\n";
#endif
    handleSetup();
  }

  short int Process2to2::set(short int n, const int& ip, const bool iap)
  {
    if (n > 4 || n <= 0)
      return 0;
    n--;

#ifdef DEBUG
    std::cout << "Assigning function values\n";
#endif
    p[n] = ip;
    ap[n] = iap;
#ifdef DEBUG
    std::cout << "Calling handleSetup()\n";
#endif
    handleSetup();
#ifdef DEBUG
    std::cout << "Exiting set()\n";
#endif
    return n;
  }

  std::string Process2to2::getKey() const { return Key; }

  std::string Process2to2::getName() const
  {
    std::ostringstream sout;
    for (auto i = 0; i < 2; i++)
    {
      if (!ap[i])
        sout << "anti_";
      sout << corr::part_names[p[i]] << ' ';
    }

    sout << "-> ";

    for (auto i = 2; i < 4; i++)
    {
      if (!ap[i])
        sout << "anti_";
      sout << corr::part_names[p[i]] << ' ';
    }

    std::string name = sout.str();
    name.erase(name.end() - 1); // erasing the last space
    return name;
  }

  std::string Process2to2::getMname() const
  {
    std::ostringstream sout;
    for (auto i = 0; i < 2; i++)
    {
      if (!ap[i])
        sout << "anti_";
      sout << corr::part_names[p[i]] << '_';
    }
    sout << "to_";
    for (auto i = 2; i < 4; i++)
    {
      if (!ap[i])
        sout << "anti_";
      sout << corr::part_names[p[i]] << '_';
    }
    std::string name = sout.str();
    name.erase(name.end() - 1); // erasing the last '-'
    return name;
  }

  bool Process2to2::EP13calculation(const Param_t& input, const real_t& sqrts, real_t& p1, real_t& p3, real_t& E1,
                                    real_t& E3) const
  {
    // In this function we suppose that the running is already taken care of

    if ((sqrts < getMass(3, input) + getMass(4, input)) || (sqrts < getMass(1, input) + getMass(2, input)))
      return false;
    real_t s = SQUARE(sqrts);
    E1 = (s + SQUARE(getMass(1, input)) - SQUARE(getMass(2, input))) / 2.0 / sqrts;
    if (E1 < 0.0)
      return false;
    E3 = (s + SQUARE(getMass(3, input)) - SQUARE(getMass(4, input))) / 2.0 / sqrts;
    if (E3 < 0.0)
      return false;
    real_t temp = SQUARE(E1) - SQUARE(getMass(1, input));
    if (temp < 0.0)
      return false;
    p1 = std::sqrt(temp);
    temp = SQUARE(E3) - SQUARE(getMass(3, input));
    if (temp < 0.0)
      return false;
    p3 = std::sqrt(temp);

    return true;
  }

  bool Process2to2::setKinematics(Param_t& input, const real_t& sqrts, const real_t& ctheta, real_t& p1, real_t& p3,
                                  real_t sij[5][5]) const
  {
    // In this function we suppose that the running is already taken care of
    // The most convenient thing is compute t and express all the sij in
    // terms of s, t, and the masses

    real_t E1, E3;

    if (!EP13calculation(input, sqrts, p1, p3, E1, E3))
      return false;
    real_t s = SQUARE(sqrts);

    sij[1][3] = E1 * E3 - p1 * p3 * ctheta;
    real_t t = SQUARE(getMass(1, input)) - 2 * (sij[1][3]) + SQUARE(getMass(3, input));

    sij[1][2] = 0.5 * (s - SQUARE(getMass(1, input)) - SQUARE(getMass(2, input)));
    sij[1][4] = 0.5 * (s + t - SQUARE(getMass(2, input)) - SQUARE(getMass(3, input)));
    sij[2][3] = 0.5 * (s + t - SQUARE(getMass(1, input)) - SQUARE(getMass(4, input)));
    sij[2][4] = -0.5 * (t - SQUARE(getMass(2, input)) - SQUARE(getMass(4, input)));
    sij[3][4] = 0.5 * (s - SQUARE(getMass(3, input)) - SQUARE(getMass(4, input)));

    return true;
  }

  real_t Process2to2::getSumSquaredAmpl(Param_t& input, const real_t& sqrts, const real_t& ctheta)
  {
    // This methods returns the sum of the squared amplitude for the instance of the class.
    // Therefore, it takes as input the numerical parameters, the centre-of-mass energy, and the
    // cosine of the angle.

#ifdef DEBUG
    std::cout << "Entering getSumSquaredAmpl\n";
#endif
    if (!isRunningExternal)
    {
#ifdef DEBUG
      std::cout << "Running is internal\n";
#endif
      handleRunning(input, sqrts);
    }
#ifdef DEBUG
    else
    {
      std::cout << "Running is external\n";
    }
#endif
    real_t p1, p3, sij[5][5];
    if (!setKinematics(input, sqrts, ctheta, p1, p3, sij))
      return 0.;

    update_kinematics(input, sij);

    cparam_s c_input = copy_to_c_struct(input);

    const auto SumSquaredAmpl_temp = sumSquaredAmpl(&c_input);

    return (std::isnormal(SumSquaredAmpl_temp.real) && SumSquaredAmpl_temp.real > 0.) ? SumSquaredAmpl_temp.real : 0.;
  }

  real_t Process2to2::get_process_dependent_contrib_dweff_dcos(Param_t& input, const real_t& Ecm,
                                                               const real_t& ctheta) const
  {
    // The following function computes the process-dependent contributions to dWeff/dcostheta, i.e.:
    // (8 PI) g^2_LBSM * sqrtS * peff *dWeff/dcostheta = p12*p34/Sf34*CPfac*(sum |M|^2)
    // for the given inputs

    // Running is assumed to be handled externally

    real_t p1, p3, sij[5][5];
    if (!setKinematics(input, Ecm, ctheta, p1, p3, sij))
      return 0.;
    update_kinematics(input, sij);

    cparam_s c_input = copy_to_c_struct(input);
    real_t SumSquaredAmpl = sumSquaredAmpl(&c_input).real;

    SumSquaredAmpl = (std::isnormal(SumSquaredAmpl) && SumSquaredAmpl > 0.) ? SumSquaredAmpl : 0.0;

    const real_t p12 = (Ecm > getMass(1, input) + getMass(2, input))
                           ? std::sqrt((SQUARE(Ecm) - std::pow(getMass(1, input) + getMass(2, input), 2)) *
                                       (SQUARE(Ecm) - std::pow(getMass(1, input) - getMass(2, input), 2))) /
                                 2.0 / Ecm
                           : 0.0;
    const real_t p34 = (Ecm > getMass(3, input) + getMass(4, input))
                           ? std::sqrt((SQUARE(Ecm) - std::pow(getMass(3, input) + getMass(4, input), 2)) *
                                       (SQUARE(Ecm) - std::pow(getMass(3, input) - getMass(4, input), 2))) /
                                 2.0 / Ecm
                           : 0.0;

    const real_t prefac = p12 * p34 / Sf34 * CPfac;

    return (std::isnormal(prefac) && prefac > 0.) ? prefac * SumSquaredAmpl : 0.0;
  }


  real_t Process2to2::get_g2_dweff_dcos(Param_t& input, const real_t& Ecm, const real_t& ctheta)
  {
    // This function computes [ g^2_LBSM * dWeff/dcostheta] for the given inputs
    if (!isRunningExternal)
      handleRunning(input, Ecm);
    real_t p1, p3, sij[5][5];
    if (!setKinematics(input, Ecm, ctheta, p1, p3, sij))
      return 0.;
    update_kinematics(input, sij);

    const real_t peff = 0.5 * std::sqrt(SQUARE(Ecm) - 4 * SQUARE(input.getLightestBSMmass()));

    cparam_s c_input = copy_to_c_struct(input);

    real_t SumSquaredAmpl = sumSquaredAmpl(&c_input).real;

    SumSquaredAmpl = (std::isnormal(SumSquaredAmpl) && SumSquaredAmpl > 0.) ? SumSquaredAmpl : 0.0;

    const real_t p12 = (Ecm > getMass(1, input) + getMass(2, input))
                           ? std::sqrt((SQUARE(Ecm) - std::pow(getMass(1, input) + getMass(2, input), 2)) *
                                       (SQUARE(Ecm) - std::pow(getMass(1, input) - getMass(2, input), 2))) /
                                 2.0 / Ecm
                           : 0.0;
    const real_t p34 = (Ecm > getMass(3, input) + getMass(4, input))
                           ? std::sqrt((SQUARE(Ecm) - std::pow(getMass(3, input) + getMass(4, input), 2)) *
                                       (SQUARE(Ecm) - std::pow(getMass(3, input) - getMass(4, input), 2))) /
                                 2.0 / Ecm
                           : 0.0;

    const real_t prefac = p12 * p34 / (Sf34 * Ecm * peff) * CPfac;

    const real_t newprefac = 1. / 8. / M_PI;

    return (std::isnormal(prefac) && prefac > 0.) ? newprefac * prefac * SumSquaredAmpl : 0.0;
  }

  real_t Process2to2::getDiffWeffContrib(Param_t& input, const real_t& Ecm, const real_t& ctheta)
  {
    // Returns dWeff/dcostheta(Ecm, ctheta) for the given input
    // We divide by (g_DM)^2 the result of the method get_g2_dweff_dcos
    return get_g2_dweff_dcos(input, Ecm, ctheta) / SQUARE(input.getLightestBSMdof());
  }


  real_t Process2to2::getDiffCrossSection(Param_t& input, const real_t& Ecm, const real_t& ctheta)
  {
    // Returns the differential cross-section at the given centre-of-mass energy Ecm and cosine
    // of the angle ctheta

    if (!isRunningExternal)
      handleRunning(input, Ecm);
    real_t p1, p3, sij[5][5];
    if (!setKinematics(input, Ecm, ctheta, p1, p3, sij))
      return 0.;
    update_kinematics(input, sij);

    cparam_s c_input = copy_to_c_struct(input);
    real_t SumSquaredAmpl = sumSquaredAmpl(&c_input).real;

    SumSquaredAmpl = (std::isnormal(SumSquaredAmpl)) ? SumSquaredAmpl : 0.0;

    const real_t kinfact = 1.0 / combinFac / 64.0 / M_PI / Ecm / Ecm / p1;
    const real_t jacobian = 2. * p3;
    return (std::isnormal(kinfact)) ? kinfact * jacobian * SumSquaredAmpl : 0.0;
  }

  real_t Process2to2::getTotalCrossSection(Param_t& input, const real_t& Ecm, real_t* discr)
  {
    if (!isRunningExternal)
      handleRunning(input, Ecm);

    real_t p1, p3, sij[5][5];

    // Returning 0 if the process is kinematically forbidden
    if (!setKinematics(input, Ecm, 0., p1, p3, sij))
      return 0.;

    update_kinematics(input, sij);

    cparam_s c_input = copy_to_c_struct(input);

    // Defining a lambda function for the differential cross section
    // Everything has to be passed by value because we use multithreading
    auto lambda_opti = [=, this](real_t cosine) mutable
    {
      real_t p1_local, p3_local, sij_local[5][5];
      if (!setKinematics(input, Ecm, cosine, p1_local, p3_local, sij_local))
        return 0.;
      update_kinematics(input, sij_local);
      c_input = update_c_kinematics(input, c_input);
      const real_t kinfac_angle = p3_local / p1_local;
      const real_t SumSquaredAmpl = getSumSquaredAmpl_ptr()(&c_input).real;
      return (std::isnormal(SumSquaredAmpl) && SumSquaredAmpl > 0.) ? SumSquaredAmpl * kinfac_angle : 0.0;
    };

    real_t integral, discrepancy;
    const real_t maxdiscrepancy = 5.0e-3;

    integral = advmath::integrate_gauss_comparative(-1., 1., lambda_opti, maxdiscrepancy, &discrepancy);

    if (discrepancy > maxdiscrepancy)
      integral = advmath::integrate_trap(-1., 1., lambda_opti, maxdiscrepancy, &discrepancy);
    if (discr != nullptr)
      *discr = discrepancy;

    const real_t kinfac_onlysqrtS = 1.0 / combinFac / 32.0 / M_PI / Ecm / Ecm;

    return integral * kinfac_onlysqrtS;
  }

  real_t Process2to2::getDiffW12Contrib(Param_t& input, const real_t& Ecm, const real_t& ctheta,
                                        const bool wantWeffPrefac)
  {
    /*   Function for testing. It returns
     *   g_1 g_2 p1 W_12->34,
     *   using the identity
     *   g_1 g_2 p1 W_12->34 = (g_1 g_2 p_1)*4* sqrts p1 * dsigma_12->34 / dcostheta
     */

    if (!isRunningExternal)
      handleRunning(input, Ecm);
    real_t p1, p3, sij[5][5];
    if (!setKinematics(input, Ecm, ctheta, p1, p3, sij))
      return 0.;
    update_kinematics(input, sij);

    cparam_s c_input = copy_to_c_struct(input);
    real_t SumSquaredAmpl = sumSquaredAmpl(&c_input).real;

    SumSquaredAmpl = (std::isnormal(SumSquaredAmpl)) ? SumSquaredAmpl : 0.0;

    const real_t xsec_kinfact = 1.0 / combinFac / 64.0 / M_PI / Ecm / Ecm / p1;
    const real_t xsec_jacobian = 2. * p3;
    const real_t W_kinfac = 4. * Ecm * p1;
    const real_t Weff_prefac = (wantWeffPrefac ? p1 * part_g[0] * part_g[1] : 1.);
    return (std::isnormal(xsec_kinfact)) ? Weff_prefac * xsec_kinfact * xsec_jacobian * W_kinfac * SumSquaredAmpl : 0.0;
  }

  real_t Process2to2::getTotalW12Contrib(Param_t& input, const real_t& sqrts, real_t* discr)
  {
    real_t p1, p3, sij[5][5];
    // Returning 0 if the process is kinematically forbidden
    if (!setKinematics(input, sqrts, 0., p1, p3, sij))
      return 0.;
    update_kinematics(input, sij);
    return 4. * sqrts * p1 * getTotalCrossSection(input, sqrts, discr);
  }

  real_t Process2to2::getInitialSpeed(Param_t& input, const real_t& sqrts) const
  {
    real_t ctheta = 0., p1, p3, sij[5][5];
    setKinematics(input, sqrts, ctheta, p1, p3, sij);
    return std::sqrt(1. - SQUARE(input.masses_vector[p[0]] * input.masses_vector[p[1]]) / SQUARE(sij[1][2]));
  }

#define FINAL_SPEED std::sqrt(1. - SQUARE(input.masses_vector[p[2]] * input.masses_vector[p[3]]) / SQUARE(sij[3][4]))
  real_t Process2to2::getFinalSpeed(Param_t& input, const real_t& sqrts) const
  {
    real_t ctheta = 0., p1, p3, sij[5][5];
    setKinematics(input, sqrts, ctheta, p1, p3, sij);
    return FINAL_SPEED;
  }

  Process2to2::Process2to2(const Process2to2& other)
          : Key(other.Key), Exists(other.Exists), sumSquaredAmpl(other.sumSquaredAmpl), Sf34(other.Sf34),
            CPfac(other.CPfac), combinFac(other.combinFac), runptr(nullptr), isRunDataExternal(false),
            isRunningExternal(false), haveToFreerunptr(false)
  {
    // I need to deep copy the C-style vectors
    for (int i = 0; i < 4; i++)
    {
      p[i] = other.p[i];
      ap[i] = other.ap[i];
      part_g[i] = other.part_g[i];
    }
  }

  bool Process2to2::operator==(const Process2to2& other) const { return this->Key.get() == other.Key.get(); }

  bool Process2to2::operator!=(const Process2to2& other) const { return !(*this == other); }

  Process2to2& Process2to2::operator=(const Process2to2& other)
  {
#ifdef DEBUG_COPY
    std::cout << "Called Process2to2::operator=\n";
#endif
    if (this == &other)
      return *this;

    // I need to deep copy the C-style vectors
    for (int i = 0; i < 4; i++)
    {
      p[i] = other.p[i];
      ap[i] = other.ap[i];
      part_g[i] = other.part_g[i];
    }

    // Copying all the other parameters
    Key = other.Key;
    Exists = other.Exists;
    sumSquaredAmpl = other.sumSquaredAmpl;
    Sf34 = other.Sf34;
    CPfac = other.CPfac;
    combinFac = other.combinFac;
    runptr = nullptr;
    isRunDataExternal = false;
    isRunningExternal = false;
    haveToFreerunptr = false;

    return *this;
  }


#undef BASIC_INITIALISATION

} // end of namespace __SPEC_LIB_NAME__
