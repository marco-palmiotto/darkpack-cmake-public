#include "params_new.hpp"
#include "compare.hpp"
#include "correspondance.hpp"
#include "leshouchesfrommarty.hpp"

#define MIN(X, Y) ((X < Y) ? X : Y)

namespace __SPEC_LIB_NAME__
{
  Param_t::Param_t(const std::string filename)
  {
#ifdef LESHOUCESFROMMARTY_INCLUDED
    empty_init();
    (*this) = readmodule::ReadLHA(filename);
#else
    std::cerr << "Error in calling Param_t::Param_t(const std::string filename)! \n"
              << "Is leshouchesfrommarty.h empty?\n";
#endif
  }

  bool darkpackparam_t::compare_with_darkpackparam(const darkpackparam_t& other) const
  {
    return (this == &other ||
            ((Running_scale == other.Running_scale) && (alpha_str_Mz == other.alpha_str_Mz) &&
             (m_Z_pole == other.m_Z_pole) && (mbmb == other.mbmb) && (mtop_pole == other.mtop_pole) &&
             (alpha_em_MZ_inv == other.alpha_em_MZ_inv) && (Gfermi == other.Gfermi) &&
             (m_tau_pole == other.m_tau_pole) && (m_el_pole == other.m_el_pole) && (m_mu_pole == other.m_mu_pole) &&
             (m_d2GeV == other.m_d2GeV) && (m_u2GeV == other.m_u2GeV) && (m_s2GeV == other.m_s2GeV) &&
             (m_c_m_c == other.m_c_m_c) && (ckm_lambda == other.ckm_lambda) && (ckm_A == other.ckm_A) &&
             (ckm_rhobar == other.ckm_rhobar) && (ckm_etabar == other.ckm_etabar) && (part_lbsm == other.part_lbsm) &&
             (m_lbsm == other.m_lbsm) && (em_charge_lbsm == other.em_charge_lbsm) && (g_lbsm == other.g_lbsm)));
  }

  bool darkpackparam_t::operator==(const darkpackparam_t& other) const { return compare_with_darkpackparam(other); }

  unsigned int Param_t::change_values(const Param_t& other)
  {
    if (!compare_with_darkpackparam(other))
      return 1;

    if (part_lbsm != other.part_lbsm)
      return 2;

    this->param_t::operator=(other);
    return 0;
  }

  bool Param_t::operator==(const Param_t& other) const
  {
#ifdef DEBUG
    bool temp = compare_param_t(*this, other);
    if (!temp)
    {
      std::cout << "The 2 param_t are different\n";
    }
#endif
    return (compare_param_t(*this, other) && darkpackparam_t::operator==(other));
  }

  // Reads data from the file "namefile" in LHA format and copies the values in the
  // data structure
  void Param_t::ReadSMparams(const std::string namefile)
  {
    mty::lha::LHAFileData data = mty::lha::Reader::readFile(namefile);
#ifdef DEBUG
    std::cout << data << std::endl;
#endif

    // Standard Model Parameters in LHA conventuions
    // If the value is not provided in the file, the default pdg value is assigned
    alpha_em_MZ_inv = get_value(data, "SMINPUTS", 1, pdgValue::alpha_em_MZ_inv);
    Gfermi = get_value(data, "SMINPUTS", 2, pdgValue::Gfermi);
    alpha_str_Mz = get_value(data, "SMINPUTS", 3, pdgValue::alpha_str_Mz);
    m_Z_pole = get_value(data, "SMINPUTS", 4, pdgValue::m_Z_pole);
    mbmb = get_value(data, "SMINPUTS", 5, pdgValue::mb_mb);
    mtop_pole = get_value(data, "SMINPUTS", 6, pdgValue::m_top_pole);

    m_tau_pole = get_value(data, "SMINPUTS", 7, pdgValue::m_tau_pole);

    m_el_pole = get_value(data, "SMINPUTS", 11, pdgValue::m_el_pole);

    m_mu_pole = get_value(data, "SMINPUTS", 13, pdgValue::m_mu_pole);

    // quark masses at 2GeV
    m_d2GeV = get_value(data, "SMINPUTS", 21, pdgValue::m_d2GeV);
    m_u2GeV = get_value(data, "SMINPUTS", 22, pdgValue::m_u2GeV);
    m_s2GeV = get_value(data, "SMINPUTS", 23, pdgValue::m_s2GeV);

    m_c_m_c = get_value(data, "SMINPUTS", 24, pdgValue::m_c_m_c);

    // CKM parameters
    ckm_lambda = get_value(data, "VCKMIN", 1, pdgValue::ckm_lambda);
    ckm_A = get_value(data, "VCKMIN", 2, pdgValue::ckm_A);
    ckm_rhobar = get_value(data, "VCKMIN", 3, pdgValue::ckm_rhobar);
    ckm_etabar = get_value(data, "VCKMIN", 4, pdgValue::ckm_etabar);

#ifdef DEBUG
    std::cout << "SM parameters have been read, the darkpackparam_t struct is\n" << std::endl;
    Print();
#endif
  }


  void Param_t::setLightestBSMparticle()
  {
    using namespace corr;
    if (!aremassesassigned)
      AssignMassesVector();
#ifdef DEBUG
    std::cout << "Assigning the first value\n";
#endif
    size_t count = 0;
    part_lbsm = corr::Part_t(bsm_particles[count]);
    m_lbsm = std::abs(masses_vector[bsm_particles[count]].get());
    em_charge_lbsm = corr::part_charge.at(bsm_particles[count]);
    g_lbsm = corr::part_hel_dof.at(bsm_particles[count]);

    if (SIZEPHYSICALBSM == 1)
      return;
    for (count = 1; count < SIZEPHYSICALBSM; count++)
    {
#ifdef DEBUG
      std::cout << "Assigning value " << count << std::endl;
#endif
      bool isless = (std::abs(masses_vector[bsm_particles[count]].get()) < m_lbsm.get());
      if (isless)
      {
        part_lbsm = corr::Part_t(bsm_particles[count]);
        m_lbsm = std::abs(masses_vector[bsm_particles[count]].get());
        em_charge_lbsm = corr::part_charge.at(bsm_particles[count]);
        g_lbsm = corr::part_hel_dof.at(bsm_particles[count]);
      }
    }
  }

  void Param_t::AssignMassesVector()
  {
    using namespace corr;
    for (int i = 1; i <= TOTAL_PARTICLES; i++)
    {
      masses_vector[i] = std::abs(getMassFirst(i, *this));
    }
    aremassesassigned = true;
  };

  void Param_t::refresh()
  {
    AssignMassesVector();
    setLightestBSMparticle();
  }

  void darkpackparam_t::print(std::ostream& out) const
  {
    out << "Printing the darkpackparam_t struct\n";
    Running_scale.print(out);
    alpha_str_Mz.print(out);
    m_Z_pole.print(out);
    mbmb.print(out);
    mtop_pole.print(out);
    alpha_em_MZ_inv.print(out);
    Gfermi.print(out);
    m_tau_pole.print(out);
    m_el_pole.print(out);
    m_mu_pole.print(out);
    m_d2GeV.print(out);
    m_u2GeV.print(out);
    m_s2GeV.print(out);
    m_c_m_c.print(out);
    ckm_lambda.print(out);
    ckm_A.print(out);
    ckm_rhobar.print(out);
    ckm_etabar.print(out);

    if (aremassesassigned)
      for (size_t i = 1; i <= corr::TOTAL_PARTICLES; i++)
      {
        masses_vector[i].print(out);
      }
    else
      out << "The masses_vector is not initialised\n";

    if (arewidthsassigned)
      for (size_t i = 1; i <= corr::TOTAL_PARTICLES; i++)
      {
        widths_vector[i].print(out);
      }
    else
      out << "The widths_vector is not initialised\n";

    part_lbsm.hasValue() ? (out << "The lightest BSM particle is " << corr::part_names[part_lbsm] << '\n')
                         : (out << "The lightest BSM particle is not set" << '\n');

    m_lbsm.print(out);
    g_lbsm.print(out);
  }

  void Param_t::Print(std::ostream& out) const
  {
    darkpackparam_t::print(out);
    param_t::print(out);
  }
} // End of namespace __SPEC_LIB_NAME__
