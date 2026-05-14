#include "process_1to2.hpp"
#include "copy_to_c_struct.hpp"

template <typename T> static inline T SQRT(const T x) { return pow(x, 0.5); }

namespace __SPEC_LIB_NAME__
{

#ifndef GEVTOPBARN
  #define GEVTOPBARN (0.389379 * 1.0e+9)
#endif

// Constant to convert GeV^-2 to pbarn (surfaces)
#ifndef GEVTOPBARN
  #define GEVTOPBARN (0.389379 * 1.0e+9)
#endif

// Constant to convert GeV^-2 to cm^3/s (sigma*v)
#ifndef GEVTOCM3S
  #define GEVTOCM3S 1.1673288750358198e-05
#endif

  // bool decays(int particle)
  //{
  // for (const auto& it : corr::squaredampl_1to2)
  //{
  // const auto& current_process = it.first;
  // if (current_process.size() < 3)
  // continue;
  // if (current_process[0] == EMPTYCHAR + particle || current_process[1] == EMPTYCHAR + particle)
  // return true;
  //}
  // return false;
  //}

#define BASIC_INITIALISATION                                                                                           \
  for (int i = 0; i < 3; i++)                                                                                          \
  {                                                                                                                    \
    p[i] = -1;                                                                                                         \
    ap[i] = false;                                                                                                     \
  }                                                                                                                    \
  Exists = false;                                                                                                      \
  partWidth = nullptr;                                                                                                 \
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

  Process_1to2::Process_1to2(void){BASIC_INITIALISATION}

  Process_1to2::~Process_1to2(void)
  {
#ifdef DEBUG_COPY
    std::cerr << "Calling ~Process2to2\n";
#endif
    Key.get().clear();
    if (haveToFreerunptr)
      delete runptr;
  }

  std::string Process_1to2::find_key()
  {
    std::string key;
    using namespace corr;
    if (!is_complete())
    {
#ifdef DEBUG
      std::cerr << "Process is not complete\n";
#endif
      Exists = false;
      key.push_back(EMPTYCHAR);
      return key;
    }
    key.reserve(6);

    for (auto i = 0; i < 3; i++)
    {
      if (!ap[i])
        key.push_back(ANTICHAR);
      key.push_back(EMPTYCHAR + p[i]);
    }
    auto it = squaredampl_1to2.find(key);
    if (it != squaredampl_1to2.end())
    {
#ifdef DEBUG
      std::cerr << "Found key\n";
#endif
      Key = key;
      return key;
    }
    // If the key is not found, permutation have to be used
    key.clear();
    std::array<int, 3> v = {1, 0, 2};
    for (auto i : v)
    {
      if (!ap[i])
        key.push_back(ANTICHAR);
      key.push_back(EMPTYCHAR + p[i]);
    }
    it = squaredampl_1to2.find(key);
    if (it != squaredampl_1to2.end())
    {
#ifdef DEBUG
      std::cerr << "Found key swapping p0 and p1\n";
      std::cerr << "Names before swapping are " << get_mname() << /* ", " << getSname() << */ std::endl;
#endif
      Key = key;
      swap(&p[0], &p[1]);
      swap(&ap[0], &ap[1]);
#ifdef DEBUG
      std::cerr << "Names after swapping are " << get_mname() << /* ", " << getSname() << */ std::endl;
      scanf("%*c");
#endif
      return key;
    }

    key.clear();
    v = {2, 1, 0};
    for (auto i : v)
    {
      if (!ap[i])
        key.push_back(ANTICHAR);
      key.push_back(EMPTYCHAR + p[i]);
    }
    it = squaredampl_1to2.find(key);
    if (it != squaredampl_1to2.end())
    {
#ifdef DEBUG
      std::cerr << "Found key swapping p0 and p2\n";
      std::cerr << "Names before swapping are " << get_mname() << /* ", " << getSname() << */ std::endl;
#endif
      Key = key;
      swap(&p[0], &p[2]);
      swap(&ap[0], &ap[2]);
#ifdef DEBUG
      std::cerr << "Names after swapping are " << get_mname() << /* ", " << getSname() << */ std::endl;
      scanf("%*c");
#endif
      return key;
    }


#ifdef DEBUG
    std::cerr << "Key not found: the process does not exist\n";
#endif

    Exists = false;
    key.clear();
    key.push_back(EMPTYCHAR);
    return key;
  }

  void Process_1to2::handle_setup()
  {
    if (is_complete())
    {
      Exists = true;
      Key = find_key();

      if (Exists)
      {
#ifdef DEBUG
        std::cerr << "Setting Sf34\n";
#endif
        Sf34 = (p[1] == p[2] && ((ap[1] && ap[2]) || (!ap[1] && !ap[2]))) ? 2 : 1;
#ifdef DEBUG
        std::cerr << "Assigning values...";
#endif
        auto it = corr::squaredampl_1to2.find(Key);
#ifdef DEBUG
        std::cerr << "... partial width";
#endif
        partWidth = std::get<0>(it->second);
#ifdef DEBUG
        std::cerr << "... runptr, extdata, extrunning,";
#endif
        runptr = nullptr;
        isRunDataExternal = false;
        isRunningExternal = false;
#ifdef DEBUG
        std::cerr << "... helicity dofs are \n";
#endif

        for (size_t i = 0; i < 3; i++)
        {
          part_g[i] = corr::part_hel_dof[p[i]];
#ifdef DEBUG
          std::cerr << part_g[i] << ' ';
#endif
        }

#ifdef DEBUG
        std::cerr << "\nExiting handleSetup of " << get_name() << "\n";
#endif
        return;
      }
      else
      {
#ifdef DEBUG
        std::cerr << "Process does not exist: resetting variables\n";
#endif
        Sf34 = 0;
        partWidth = nullptr;
        runptr = nullptr;
        isRunDataExternal = false;
        isRunningExternal = false;
        return;
      }
    }
#ifdef DEBUG
    std::cerr << "Process is not complete\n";
#endif
  }


  Process_1to2::Process_1to2(const std::array<int, 3>& field, const std::array<bool, 3>& is_particle)
  {
#ifdef DEBUG
    std::cerr << "Constructor called\n";
#endif
    BASIC_INITIALISATION
#ifdef DEBUG
    std::cerr << "Pushing back values\n";
#endif
    for (int i = 0; i < 3; i++)
    {
      p[i] = field[i];
      ap[i] = is_particle[i];
    }
#ifdef DEBUG
    std::cerr << "Data pushed back\n";
#endif
    handle_setup();
  }

  Process_1to2::Process_1to2(const std::array<Insertion, 3>& v)
  {
#ifdef DEBUG
    std::cerr << "Constructor called\n";
#endif
    BASIC_INITIALISATION
#ifdef DEBUG
    std::cerr << "Pushing back values\n";
#endif
    for (int i = 0; i < 3; i++)
    {
      p[i] = v[i].field;
      ap[i] = v[i].part;
    }
#ifdef DEBUG
    std::cerr << "Pushing back values\n";
#endif
    handle_setup();
  }

  short int Process_1to2::set(short int n, const int& ip, const bool iap)
  {
    if (n > 3 || n <= 0)
      return 0;
    n--;

#ifdef DEBUG
    std::cerr << "Assigning function values\n";
#endif
    p[n] = ip;
    ap[n] = iap;
#ifdef DEBUG
    std::cerr << "Calling handleSetup()\n";
#endif
    handle_setup();
#ifdef DEBUG
    std::cerr << "Exiting set()\n";
#endif
    return n;
  }

  std::string Process_1to2::get_key() const { return Key; }

  std::string Process_1to2::get_name() const
  {
    std::ostringstream sout;
    for (auto i = 0; i < 1; i++)
    {
      if (!ap[i])
        sout << "anti_";
      sout << corr::part_names[p[i]] << ' ';
    }

    sout << "-> ";

    for (auto i = 1; i < 3; i++)
    {
      if (!ap[i])
        sout << "anti_";
      sout << corr::part_names[p[i]] << ' ';
    }

    std::string name = sout.str();
    name.erase(name.end() - 1); // erasing the last space
    return name;
  }

  std::string Process_1to2::get_mname() const
  {
    std::ostringstream sout;
    for (auto i = 0; i < 1; i++)
    {
      if (!ap[i])
        sout << "anti_";
      sout << corr::part_names[p[i]] << '_';
    }
    sout << "to_";
    for (auto i = 1; i < 3; i++)
    {
      if (!ap[i])
        sout << "anti_";
      sout << corr::part_names[p[i]] << '_';
    }
    std::string name = sout.str();
    name.erase(name.end() - 1); // erasing the last '-'
    return name;
  }

  bool Process_1to2::set_kinematics(Param_t& input, real_t s_ij[5][5], real_t m_vec[3], real_t& E1, real_t& E2,
                                    real_t& p2) const
  {
    // The following function assumes that running is already taken care of
    for (size_t i = 1; i < 4; i++)
    {
      m_vec[i - 1] = get_mass(i, input);
    }

    E1 = (m_vec[0] * m_vec[0] + m_vec[1] * m_vec[1] - m_vec[2] * m_vec[2]) / (2 * m_vec[0]);
    E2 = (m_vec[0] * m_vec[0] + m_vec[2] * m_vec[2] - m_vec[1] * m_vec[1]) / (2 * m_vec[0]);
    p2 = E1 * E1 - m_vec[1] * m_vec[1];
    s_ij[1][2] = m_vec[0] * E1;
    s_ij[1][3] = m_vec[0] * E2;
    s_ij[2][3] = E1 * E2 + p2;
    update_kinematics(input, s_ij);

    return true;
  }


  void Process_1to2::three_body_approximation(const Param_t& input, real_t& part_width_temp, const real_t m_vec[3])
  {
    if (!((p[1] == corr::Z && p[2] == corr::Z) || (p[1] == corr::W && p[2] == corr::W)))
      return;

    real_t sin_theta = sin(input.theta_W);
    real_t x = pow(m_vec[1], 2.) / pow(m_vec[0], 2.);
    real_t g_fermi = input.Gfermi;

    part_width_temp *= 3 * g_fermi * pow(m_vec[1], 4.) * SQRT(2.) /
                       (2 * (pow(M_PI * m_vec[0], 2.)) * SQRT(fabs(1 - 4 * x)) * (1 - 4 * x + 12 * pow(x, 2.)));

    switch (p[1])
    {
    case (corr::Z):
      part_width_temp *= 2 * (7. / 12. - 10. / 9. * pow(sin_theta, 2.) + 40. / 9. * pow(sin_theta, 4.));
      break;
    default:
      break;
    }
    printf("Process %s has a partial width of %e\n", get_name().c_str(), part_width_temp);
    return;
  }

  real_t Process_1to2::compute_partial_width(Param_t& input)
  {
    real_t s_ij[5][5];
    real_t m_vec[3];
    real_t E1, E2, p2;

    for (size_t i = 1; i < 4; i++)
    {
      m_vec[i - 1] = get_mass(i, input);
    }

    if (!isRunningExternal)
      handle_running(input, m_vec[0]);
    if (!set_kinematics(input, s_ij, m_vec, E1, E2, p2))
      return 0.;

    real_t part_width_temp = partWidth(input).real();

    if (!is_allowed_at_zero_momentum(input))
    {
      if (p[1] == corr::W || p[1] == corr::Z)
        three_body_approximation(input, part_width_temp, m_vec);
      else
        return 0;
    }
    printf("Process %s has a partial width of %e\n", get_name().c_str(), part_width_temp);
    // std::cerr << "Partial width for process " << get_name() << " is " << part_width_temp << " GeV\n";
    return (std::isnormal(part_width_temp) && part_width_temp > 0.) ? part_width_temp : 0.;
  }
  /*
  real_t Process_1to2::get_sum_squared_ampl(Param_t& input)
  {
    real_t s_ij[5][5];
    real_t m_vec[3];
    real_t E1, E2, p2;
    real_t part_width = compute_partial_width(input);

    if (!set_kinematics(input, s_ij, m_vec, E1, E2, p2))
      return 0.;

    real_t sum_squared_amp = part_width * 8 * M_PI * pow(M, 2.) / SQRT(fabs(p2));
    return (std::isnormal(part_width_temp) && part_width_temp > 0.) ? sum_squared_amp : 0;
  }
  */
  real_t Process_1to2::compute_branching_ratio(Param_t& input)
  {
    real_t s_ij[5][5];
    real_t m_vec[3];
    real_t E1, E2, p2;

    real_t partial_width = compute_partial_width(input);
    //* pow(sin(input.theta_W), 2.);
    real_t total_width = input.widths_vector[get_field(1)];

    return (std::isnormal(total_width) && total_width > 0.) ? partial_width / total_width : 0.;
  }

  bool Process_1to2::operator==(const Process_1to2& other) const { return this->Key.get() == other.Key.get(); }

  bool Process_1to2::operator!=(const Process_1to2& other) const { return !(*this == other); }

  Process_1to2& Process_1to2::operator=(const Process_1to2& other)
  {
#ifdef DEBUG_COPY
    std::cerr << "Called Process_1to2::operator=\n";
#endif
    if (this == &other)
      return *this;

    // I need to deep copy the C-style vectors
    for (int i = 0; i < 3; i++)
    {
      p[i] = other.p[i];
      ap[i] = other.ap[i];
      part_g[i] = other.part_g[i];
    }

    // Copying all the other parameters
    Key = other.Key;
    Exists = other.Exists;
    Sf34 = other.Sf34;
    runptr = nullptr;
    isRunDataExternal = false;
    isRunningExternal = false;
    haveToFreerunptr = false;

    return *this;
  }
  // Copy constructor
  Process_1to2::Process_1to2(const Process_1to2& other)
          : Key(other.Key), Exists(other.Exists), partWidth(other.partWidth), Sf34(other.Sf34), runptr(nullptr),
            isRunDataExternal(false), isRunningExternal(false), haveToFreerunptr(false)
  {
    // I need to deep copy the C-style vectors
    for (int i = 0; i < 3; i++)
    {
      p[i] = other.p[i];
      ap[i] = other.ap[i];
      part_g[i] = other.part_g[i];
    }
  }



#undef BASIC_INITIALISATION
} // namespace __SPEC_LIB_NAME__
  // namespace __SPEC_LIB_NAME__
  // end of namespace __SPEC_LIB_NAME__
