#include "process_1to2.hpp"
#include "copy_to_c_struct.hpp"

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

  Process_1to2::Process_1to2(void){BASIC_INITIALISATION}

  Process_1to2::~Process_1to2(void)
  {
#ifdef DEBUG_COPY
    std::cout << "Calling ~Process2to2\n";
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
      std::cout << "Process is not complete\n";
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
      std::cout << "Found key\n";
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
      std::cout << "Found key swapping p0 and p1\n";
      std::cout << "Names before swapping are " << get_mname() << /* ", " << getSname() << */ std::endl;
#endif
      Key = key;
      swap(&p[0], &p[1]);
      swap(&ap[0], &ap[1]);
#ifdef DEBUG
      std::cout << "Names after swapping are " << get_mname() << /* ", " << getSname() << */ std::endl;
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
      std::cout << "Found key swapping p0 and p2\n";
      std::cout << "Names before swapping are " << get_mname() << /* ", " << getSname() << */ std::endl;
#endif
      Key = key;
      swap(&p[0], &p[2]);
      swap(&ap[0], &ap[2]);
#ifdef DEBUG
      std::cout << "Names after swapping are " << get_mname() << /* ", " << getSname() << */ std::endl;
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

  void Process_1to2::handle_setup()
  {
    if (is_complete())
    {
      Exists = true;
      Key = find_key();

      if (Exists)
      {
#ifdef DEBUG
        std::cout << "Setting Sf34\n";
#endif
        Sf34 = (p[1] == p[2] && ((ap[1] && ap[2]) || (!ap[1] && !ap[2]))) ? 2 : 1;
#ifdef DEBUG
        std::cout << "Assigning values...";
#endif
        auto it = corr::squaredampl_1to2.find(Key);
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
        std::cout << "\nExiting handleSetup of " << get_name() << "\n";
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


  Process_1to2::Process_1to2(const std::array<int, 3>& field, const std::array<bool, 3>& is_particle)
  {
#ifdef DEBUG
    std::cout << "Constructor called\n";
#endif
    BASIC_INITIALISATION
#ifdef DEBUG
    std::cout << "Pushing back values\n";
#endif
    for (int i = 0; i < 3; i++)
    {
      p[i] = field[i];
      ap[i] = is_particle[i];
    }
#ifdef DEBUG
    std::cout << "Data pushed back\n";
#endif
    handle_setup();
  }

  Process_1to2::Process_1to2(const std::array<Insertion, 3>& v)
  {
#ifdef DEBUG
    std::cout << "Constructor called\n";
#endif
    BASIC_INITIALISATION
#ifdef DEBUG
    std::cout << "Pushing back values\n";
#endif
    for (int i = 0; i < 3; i++)
    {
      p[i] = v[i].field;
      ap[i] = v[i].part;
    }
#ifdef DEBUG
    std::cout << "Pushing back values\n";
#endif
    handle_setup();
  }

  short int Process_1to2::set(short int n, const int& ip, const bool iap)
  {
    if (n > 3 || n <= 0)
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
    handle_setup();
#ifdef DEBUG
    std::cout << "Exiting set()\n";
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
    if (!Exists || !is_allowed_at_zero_momentum(input))
      return false;

    for (size_t i = 1; i < 4; i++)
    {
      m_vec[i] = get_mass(1, input);
    }

    E1 = (m_vec[0] * m_vec[0] + m_vec[1] * m_vec[1] - m_vec[2] * m_vec[2]) / (2 * m_vec[0]);
    E2 = (m_vec[0] * m_vec[0] + m_vec[2] * m_vec[2] - m_vec[1] * m_vec[1]) / (2 * m_vec[0]);
    p2 = E1 * E1 - m_vec[1] * m_vec[1];
    s_ij[2][3] = E1 * E2 + p2;
    update_kinematics(input, s_ij);

    return true;
  }

  real_t Process_1to2::get_partial_width(Param_t& input)
  {
    real_t s_ij[5][5];
    real_t m_vec[3];
    real_t E1, E2, p2;

    if (!set_kinematics(input, s_ij, m_vec, E1, E2, p2))
      return 0.;
    if (!isRunningExternal)
      handle_running(input, m_vec[0]);

    real_t sum_sq_ampl = sumSquaredAmpl(input).real();
    const auto partial_width_temp = sum_sq_ampl * p2 / m_vec[0] / m_vec[0] / 8.0 / M_PI;
    return (std::isnormal(partial_width_temp) && partial_width_temp > 0.) ? partial_width_temp : 0.;
  }

  real_t Process_1to2::get_branching_ratio(Param_t& input)
  {
    real_t s_ij[5][5];
    real_t m_vec[3];
    real_t E1, E2, p2;

    real_t partial_width = get_partial_width(input);
    real_t total_width = input.widths_vector[get_field(0)];

    return (std::isnormal(total_width) && total_width > 0.) ? partial_width / total_width : 0.;
  }

  bool Process_1to2::operator==(const Process_1to2& other) const { return this->Key.get() == other.Key.get(); }

  bool Process_1to2::operator!=(const Process_1to2& other) const { return !(*this == other); }

  Process_1to2& Process_1to2::operator=(const Process_1to2& other)
  {
#ifdef DEBUG_COPY
    std::cout << "Called Process_1to2::operator=\n";
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
  // Copy constructor
  Process_1to2::Process_1to2(const Process_1to2& other)
          : Key(other.Key), Exists(other.Exists), sumSquaredAmpl(other.sumSquaredAmpl), Sf34(other.Sf34),
            CPfac(other.CPfac), combinFac(other.combinFac), runptr(nullptr), isRunDataExternal(false),
            isRunningExternal(false), haveToFreerunptr(false)
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
