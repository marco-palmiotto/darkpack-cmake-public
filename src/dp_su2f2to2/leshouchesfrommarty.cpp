#include "leshouchesfrommarty.hpp"

namespace su2f2to2::readmodule
{

  void AssignSMparams(Param_t& param, const std::string namefile)
  {

    mty::lha::LHAFileData data = mty::lha::Reader::readFile(namefile);

    double thetaW, Mw;
    computeThetaWandMw(1. / param.alpha_em_MZ_inv.get(), param.Gfermi.get(), param.m_Z_pole.get(), thetaW, Mw);

    param.theta_W = thetaW;
    param.e_em = std::sqrt(4 * M_PI / param.alpha_em_MZ_inv);
    param.m_Z = get_value(data, "MASS", 23, param.m_Z_pole.get());
    param.m_W = get_value(data, "MASS", 24, Mw);
    param.m_h = get_value(data, "MASS", 25, pdgValue::m_h);

#ifdef DEBUG
    std::cout << "Parameters theta_W, e_m, m_Z, m_W, m_h have been initialised:\n";
    param.Print();
#endif

    param.m_d = param.m_d2GeV.get();
    param.m_u = param.m_u2GeV.get();
    param.m_s = param.m_s2GeV.get();
    param.m_c = param.m_c_m_c.get();

#ifdef DEBUG
    std::cout << "Parameters m_u, m_d, m_s, m_c have been initialised:\n";
    param.Print();
#endif

    param.Running_scale = -2.;
    RunningSM run(param);

#ifdef DEBUG
    std::cout << "Assigning quark masses and g_s at the scale of m_Z\n";
#endif
    run.HandleParamRunning(param, param.m_Z);

#ifdef DEBUG
    std::cout << "Light quark masses have been initialised:\n";
    param.Print();
#endif


    param.m_e = param.m_el_pole.get();
    param.m_mu = param.m_mu_pole.get();
    param.m_tau = param.m_tau_pole.get();

    param.reg_prop = 0;

    // CKM matrix assignment
    real_t s12, s13, s23, c12, c13, c23;
    complex_t expid;
    //     const complex_t I{0.0,1.0};
    complex_t CKM[4][4];

    s12 = param.ckm_lambda;
    s23 = param.ckm_A * param.ckm_lambda * param.ckm_lambda;

    const complex_t rhopIeta = {param.ckm_rhobar.get(), param.ckm_etabar.get()};

    s13 = std::abs(param.ckm_A * std::pow(param.ckm_lambda, 3.) *
                   (rhopIeta)*std::sqrt(1. - std::pow(param.ckm_A * param.ckm_lambda * param.ckm_lambda, 2.)) /
                   std::sqrt(1. - param.ckm_lambda * param.ckm_lambda) /
                   (1. - std::pow(param.ckm_A * param.ckm_lambda * param.ckm_lambda, 2.) * (rhopIeta)));


    expid = ((complex_t)(param.ckm_A) * std::pow(param.ckm_lambda, 3.) *
             (rhopIeta)*std::sqrt(1. - std::pow(param.ckm_A * param.ckm_lambda * param.ckm_lambda, 2.)) /
             std::sqrt(1. - param.ckm_lambda * param.ckm_lambda) /
             (1. - std::pow(param.ckm_A * param.ckm_lambda * param.ckm_lambda, 2.) * (rhopIeta))) /
            s13;
#ifdef DEBUG
    std::cout << "exp( i*delta ) = " << expid << std::endl;
#endif

    c12 = std::sqrt(1. - s12 * s12);
    c13 = std::sqrt(1. - s13 * s13);
    c23 = std::sqrt(1. - s23 * s23);

    CKM[1][1] = c12 * c13;
    CKM[1][2] = s12 * c13;
    CKM[1][3] = s13 / expid;
    //     IMCKM[1][3]=(s13/expid).imag();

    CKM[2][1] = -s12 * c23 - c12 * s23 * s13 * expid;
    //     IMCKM[2][1]=(-s12*c23-c12*s23*s13*expid).imag();

    CKM[2][2] = c12 * c23 - s12 * s23 * s13 * expid;
    //     IMCKM[2][2]=(c12*c23-s12*s23*s13*expid).imag();

    CKM[2][3] = s23 * c13;

    CKM[3][1] = s12 * s23 - c12 * c23 * s13 * expid;
    //     IMCKM[3][1]=(s12*s23-c12*c23*s13*expid).imag();

    CKM[3][2] = -c12 * s23 - s12 * c23 * s13 * expid;
    //     IMCKM[3][2]=(-c12*s23-s12*c23*s13*expid).imag();

    CKM[3][3] = c23 * c13;

    param.V_ud = CKM[1][1].real();
    param.V_us = CKM[1][2].real();
    param.V_ub = CKM[1][3];
    param.V_cd = CKM[2][1];
    param.V_cs = CKM[2][2];
    param.V_cb = CKM[2][3].real();
    param.V_td = CKM[3][1];
    param.V_ts = CKM[3][2];
    param.V_tb = CKM[3][3].real();
  }

  void ReadBSMparams(Param_t& param, const std::string namefile)
  {
    mty::lha::LHAFileData data = mty::lha::Reader::readFile(namefile);

#ifdef DEBUG
    std::cout << data << std::endl;
#endif
    // BSM parameters
    /* Put here the code to read the custom BSM blocks */
    param.theta_12_L = get_value(data, "SU2LEPTONIC", 1);
    param.theta_12_R = get_value(data, "SU2LEPTONIC", 2);
    param.theta_13_L = get_value(data, "SU2LEPTONIC", 3);
    param.theta_13_R = get_value(data, "SU2LEPTONIC", 4);
    param.theta_23_L = get_value(data, "SU2LEPTONIC", 5);
    param.theta_23_R = get_value(data, "SU2LEPTONIC", 6);
    param.g_f = get_value(data, "SU2LEPTONIC", 7);
    param.m_V = get_value(data, "SU2LEPTONIC", 8);
    param.m_phi = get_value(data, "SU2LEPTONIC", 9);
    param.m_chi_dm_2 = param.m_chi_dm_1 = param.m_cdm = get_value(data, "SU2LEPTONIC", 10);
  }

  void AssignBSMparams(Param_t& param)
  {
    param.pi = M_PI;

    // ----  Assigning widths ----
    param.Gamma_V = width_V(param).real();
    param.widths_vector[corr::V] = param.Gamma_V;
    param.Gamma_V_3 = width_V_3(param).real();
    param.widths_vector[corr::V_3] = param.Gamma_V_3;

    // You can replace the ones below with the ones from
    // a spectrum calculator
    param.widths_vector[corr::h] = width_h(param).real();
    param.widths_vector[corr::t] = width_t(param).real();
    param.widths_vector[corr::W] = width_W(param).real();
    param.widths_vector[corr::Z] = width_Z(param).real();

    for (auto elem : param.widths_vector)
      if (!elem.hasValue())
        elem = 0.;

    param.arewidthsassigned = true;
  }


  Param_t ReadLHA(const std::string name)
  {
    Param_t param;
#ifdef DEBUG
    std::cout << "Resetting widths\n";
#endif
    for (auto i = 1; i <= corr::TOTAL_PARTICLES; i++)
      param.widths_vector[i] = 0.;
#ifdef DEBUG
    std::cout << "Reading SM parameters\n";
#endif
    param.ReadSMparams(name);
#ifdef DEBUG
    std::cout << "Assigning SM parameters\n";
#endif
    AssignSMparams(param, name);
#ifdef DEBUG
    std::cout << "Reading BSM parameters\n";
#endif
    ReadBSMparams(param, name);
#ifdef DEBUG
    std::cout << "Assigning BSM parameters\n";
#endif
    AssignBSMparams(param);
#ifdef DEBUG
    std::cout << "Assigned BSM parameters\n";
    param.Print();
    std::cout << "Assigning mass vector\n";
#endif
    param.AssignMassesVector();
#ifdef DEBUG
    param.Print();
    for (auto i = 1; i <= corr::TOTAL_PARTICLES; i++)
      std::cout << static_cast<std::string>("m[") + corr::part_names[i] + "] = " << param.masses_vector[i] << std::endl;
    std::cout << "Getting lightest bsm particle mass\n";
#endif
    param.setLightestBSMparticle();
#ifdef DEBUG
    param.Print();
    std::cout << "Exiting read function\n";
#endif
    return param;
  }

} // End of namespace su2f2to2::readmodule
