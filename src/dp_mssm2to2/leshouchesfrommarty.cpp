/*  This file is not strictly necessary for the MSSM:
 *  the inputs can in fact be read thanks to the functions
 *  in leshouces.h, imported from SuperIso.
 *
 *  However, we provide this code, that shows how to read
 *  inputs directly from a lha file
 *
 *
 *  This file, as all the model-specific files, is intended
 *  to be modified by the user in case of need, for this reason
 *  we are leaving some commented code.
 *
 *  The commented code in AssignSMparams is how the input
 *  should be handled, according to the LHA conventions
 *  and our definition of the additional parameters
 *  that define the difference between a Param_t structure
 *  and a param_t structure.
 *
 *  The uncommented code treats the parameters the same way
 *  as leshouces.h, reproducing exactly the results of
 *  its functions.
 *
 */

#include "leshouchesfrommarty.hpp"

// #define DEBUG
namespace mssm2to2::readmodule
{
  void AssignSMparams(Param_t& param, const std::string& namefile)
  {
    mty::lha::LHAFileData data = mty::lha::Reader::readFile(namefile);
#ifdef DEBUG
    std::cout << "Data read\n";
#endif
    //     double thetaW, m_W;
    //     computeThetaWandMw(1./param.alpha_em_MZ_inv.get(), param.Gfermi.get(), param.m_Z_pole.get(), thetaW,
    //     m_W);

    //     param.e_em = std::sqrt(4*M_PI/param.alpha_em_MZ_inv);
    param.e_em = std::sqrt(4 * M_PI / pdg2016Value::alpha_em_inv);

    //     param.m_Z = get_value(data, "MASS", 23, param.m_Z_pole.get());
    param.m_Z = get_value(data, "SMINPUTS", 4, get_value(data, "MASS", 23, pdg2016Value::m_Z));

    //     param.m_W = get_value(data, "MASS", 24, m_W);
    param.m_W = get_value(data, "MASS", 24, pdg2016Value::m_W);

    param.m_h = get_value(data, "MASS", 25);
#ifdef WIDTH_Z
    param.Gamma_Z = get_value(data, "ALLDECAYS", 23, pdg2016Value::Gamma_Z);
#endif
#ifdef WIDTH_W
    param.Gamma_W = get_value(data, "ALLDECAYS", 24, pdg2016Value::Gamma_W);
#endif

    //     param.theta_W = thetaW;
    param.theta_W = std::acos(param.m_W / param.m_Z);

    //     // From now on, no more default parameters are called

#ifdef DEBUG
    std::cout << "Parameters theta_W, e_m, m_Z, m_W, m_h have been initialised:\n";
    param.Print();
#endif

    //     param.m_d = param.m_d2GeV.get();
    param.m_d = get_value(data, "MASS", 1, get_value(data, "SMINPUTS", 21, pdg2016Value::m_d2GeV));
    param.m_d2GeV = param.m_d;
    //     param.m_u = param.m_u2GeV.get();
    param.m_u = get_value(data, "MASS", 2, get_value(data, "SMINPUTS", 22, pdg2016Value::m_u2GeV));
    param.m_u2GeV = param.m_u;
    //     param.m_s = param.m_s2GeV.get();
    param.m_s = get_value(data, "MASS", 3, get_value(data, "SMINPUTS", 23, pdg2016Value::m_s2GeV));
    param.m_s2GeV = param.m_s;
    //     param.m_c = param.m_c_m_c.get();
    param.m_c = get_value(data, "MASS", 2, get_value(data, "SMINPUTS", 22, pdg2016Value::m_c_m_c));
    param.m_c_m_c = param.m_c;
#ifdef DEBUG
    std::cout << "Parameters m_u, m_d, m_s, m_c have been initialised\nAssigning g_s, m_t and m_b";
#endif

    //     param.Running_scale=-2.;
    //     RunningSM run(param);
    //
    //     // Assigning quark masses and g_s at the scale of m_Z
    //     run.HandleParamRunning(param, param.m_Z);

    param.g_s = std::sqrt(4 * M_PI * get_value(data, "SMINPUTS", 3, pdg2016Value::alpha_str_Mz));
    param.mtop_pole = get_value(data, "SMINPUTS", 6, get_value(data, "MASS", 6, pdg2016Value::m_top_pole));
    param.m_b = get_value(data, "SMINPUTS", 5, pdg2016Value::mb_mb);
    param.mbmb = param.m_b;
    // CKM matrix assignment
    real_t s12, s13, s23, c12, c13, c23;
    complex_t expid;
    // const complex_t I{0.0, 1.0};
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

    CKM[2][1] = -s12 * c23 - c12 * s23 * s13 * expid;

    CKM[2][2] = c12 * c23 - s12 * s23 * s13 * expid;

    CKM[2][3] = s23 * c13;

    CKM[3][1] = s12 * s23 - c12 * c23 * s13 * expid;

    CKM[3][2] = -c12 * s23 - s12 * c23 * s13 * expid;

    CKM[3][3] = c23 * c13;

    param.V_ud = CKM[1][1].real();
    param.V_us = CKM[1][2].real();
    param.V_ub_mod = s13;
    param.delta_wolf = std::log(expid).imag();
    param.V_cd = CKM[2][1];
    param.V_cs = CKM[2][2];
    param.V_cb = CKM[2][3].real();
    param.V_td = CKM[3][1];
    param.V_ts = CKM[3][2];
    param.V_tb = CKM[3][3].real();
#ifdef DEBUG
    std::cout << "Defining a running structure to compute mtop(mtop)\n";
#endif
    param.Running_scale = -2.;
    param.A_t = 0.; // dummy values for these variables so the running can be performed
    param.A_b = 0.;
    param.A_tau = 0.;
    param.m_t = 0.;
    RunningSM run(param);
    param.m_t = run.GetMtopMtop();

#ifdef DEBUG
    std::cout << "Light quark masses have been initialised:\n";
    param.Print();
#endif

#ifdef DEBUG
    std::cout << "Assigning leptons' masses\n";
#endif
    //     param.m_e = param.m_el_pole.get();
    param.m_e = get_value(data, "MASS", 11, get_value(data, "SMINPUTS", 11, pdg2016Value::m_el_pole));
    //     param.m_mu = param.m_mu_pole.get();
    param.m_mu = get_value(data, "MASS", 13, get_value(data, "SMINPUTS", 13, pdg2016Value::m_mu_pole));
    //     param.m_tau = param.m_tau_pole.get();
    param.m_tau = get_value(data, "MASS", 15, get_value(data, "SMINPUTS", 7, pdg2016Value::m_tau_pole));

    param.reg_prop = 0;
  }

  void ReadBSMparams(Param_t& param, const std::string& namefile)
  {
    // BSM parameters (the MSSM-specific parameters in this case)

    auto data = mty::lha::Reader::readFile(namefile);

    //     const double m0 = get_value(data,"MINPAR", 1);//
    //     const double m1 = get_value(data,"MINPAR", 2);//
    const double tanbeta = get_value(data, "MINPAR", 3);
    param.beta = std::atan(tanbeta);

    //     const double signmu = get_value(data,"MINPAR", 4);//
    //     const double A0 = get_value(data,"MINPAR", 5);//

    param.m_H0 = get_value(data, "MASS", 35);
    param.m_A0 = get_value(data, "MASS", 36);
    param.m_Hp = get_value(data, "MASS", 37);
    param.m_sd_L = get_value(data, "MASS", 1000001);
    param.m_sd_R = get_value(data, "MASS", 2000001);
    param.m_su_L = get_value(data, "MASS", 1000002);
    param.m_su_R = get_value(data, "MASS", 2000002);
    param.m_ss_L = get_value(data, "MASS", 1000003);
    param.m_ss_R = get_value(data, "MASS", 2000003);
    param.m_sc_L = get_value(data, "MASS", 1000004);
    param.m_sc_R = get_value(data, "MASS", 2000004);
    param.m_sb_1 = get_value(data, "MASS", 1000005);
    param.m_sb_2 = get_value(data, "MASS", 2000005);
    param.m_st_1 = get_value(data, "MASS", 1000006);
    param.m_st_2 = get_value(data, "MASS", 2000006);
    param.m_se_L = get_value(data, "MASS", 1000011);
    param.m_se_R = get_value(data, "MASS", 2000011);
    param.m_snu_e = get_value(data, "MASS", 1000012);
    param.m_smu_L = get_value(data, "MASS", 1000013);
    param.m_smu_R = get_value(data, "MASS", 2000013);
    param.m_snu_mu = get_value(data, "MASS", 1000014);
    param.m_stau_1 = get_value(data, "MASS", 1000015);
    param.m_stau_2 = get_value(data, "MASS", 2000015);
    param.m_snu_tau = get_value(data, "MASS", 1000016);
    param.m_sG = get_value(data, "MASS", 1000021);
    param.m_N_1 = get_value(data, "MASS", 1000022);
    param.m_N_2 = get_value(data, "MASS", 1000023);
    param.m_N_3 = get_value(data, "MASS", 1000025);
    param.m_N_4 = get_value(data, "MASS", 1000035);
    param.m_C_1 = get_value(data, "MASS", 1000024);
    param.m_C_2 = get_value(data, "MASS", 1000037);

    // Get block for a matrix
    std::vector<long double> temp = data.getValues("NMIX");
    if (temp.empty())
    {
      std::cout << "Error" << std::endl;
      return;
    }
    size_t N = static_cast<size_t>(std::round(std::sqrt(temp.size())));
    if (N != 4)
    {
      std::cout << "Error while reading the NMIX block" << std::endl;
      return;
    }
    csl::InitSanitizer<real_t> Nmat[4][4];
    for (size_t i = 0; i != 4; ++i)
    {
      for (size_t j = 0; j != 4; ++j)
      {
        Nmat[i][j] = static_cast<real_t>(temp[N * i + j]);
      }
    }
    temp.clear();

    temp = data.getValues("UMIX");
    if (temp.empty())
    {
      std::cout << "Error" << std::endl;
      return;
    }
    N = static_cast<size_t>(std::round(std::sqrt(temp.size())));
    if (N != 2)
    {
      std::cout << "Error while reading the NMIX block" << std::endl;
      return;
    }
    csl::InitSanitizer<real_t> U[2][2];
    for (size_t i = 0; i != 2; ++i)
    {
      for (size_t j = 0; j != 2; ++j)
      {
        U[i][j] = static_cast<real_t>(temp[N * i + j]);
      }
    }
    temp.clear();

    temp = data.getValues("VMIX");
    if (temp.empty())
    {
      std::cout << "Error" << std::endl;
      return;
    }
    N = static_cast<size_t>(std::round(std::sqrt(temp.size())));
    if (N != 2)
    {
      std::cout << "Error while reading the NMIX block" << std::endl;
      return;
    }
    csl::InitSanitizer<real_t> V[2][2];
    for (size_t i = 0; i != 2; ++i)
    {
      for (size_t j = 0; j != 2; ++j)
      {
        V[i][j] = static_cast<real_t>(temp[N * i + j]);
      }
    }
    temp.clear();

    temp = data.getValues("STOPMIX");
    if (temp.empty())
    {
      std::cout << "Error" << std::endl;
      return;
    }
    N = static_cast<size_t>(std::round(std::sqrt(temp.size())));
    if (N != 2)
    {
      std::cout << "Error while reading the NMIX block" << std::endl;
      return;
    }
    csl::InitSanitizer<real_t> StopMix[2][2];
    for (size_t i = 0; i != 2; ++i)
    {
      for (size_t j = 0; j != 2; ++j)
      {
        StopMix[i][j] = static_cast<real_t>(temp[N * i + j]);
      }
    }
    temp.clear();

    temp = data.getValues("SBOTMIX");
    if (temp.empty())
    {
      std::cout << "Error" << std::endl;
      return;
    }
    N = static_cast<size_t>(std::round(std::sqrt(temp.size())));
    if (N != 2)
    {
      std::cout << "Error while reading the NMIX block" << std::endl;
      return;
    }
    csl::InitSanitizer<real_t> SbottomMix[2][2];
    for (size_t i = 0; i != 2; ++i)
    {
      for (size_t j = 0; j != 2; ++j)
      {
        SbottomMix[i][j] = static_cast<real_t>(temp[N * i + j]);
      }
    }
    temp.clear();

    temp = data.getValues("STAUMIX");
    if (temp.empty())
    {
      std::cout << "Error" << std::endl;
      return;
    }
    N = static_cast<size_t>(std::round(std::sqrt(temp.size())));
    if (N != 2)
    {
      std::cout << "Error while reading the NMIX block" << std::endl;
      return;
    }
    csl::InitSanitizer<real_t> StauMix[2][2];
    for (size_t i = 0; i != 2; ++i)
    {
      for (size_t j = 0; j != 2; ++j)
      {
        StauMix[i][j] = static_cast<real_t>(temp[N * i + j]);
      }
    }
    temp.clear();

    // Assigning Matrices now
    param.N_B1 = -static_cast<complex_t>(Nmat[0][0]);
    param.N_B2 = -static_cast<complex_t>(Nmat[1][0]);
    param.N_B3 = -static_cast<complex_t>(Nmat[2][0]);
    param.N_B4 = -static_cast<complex_t>(Nmat[3][0]);

    param.N_W1 = -static_cast<complex_t>(Nmat[0][1]);
    param.N_W2 = -static_cast<complex_t>(Nmat[1][1]);
    param.N_W3 = -static_cast<complex_t>(Nmat[2][1]);
    param.N_W4 = -static_cast<complex_t>(Nmat[3][1]);

    param.N_d1 = -static_cast<complex_t>(Nmat[0][2]);
    param.N_d2 = -static_cast<complex_t>(Nmat[1][2]);
    param.N_d3 = -static_cast<complex_t>(Nmat[2][2]);
    param.N_d4 = -static_cast<complex_t>(Nmat[3][2]);

    param.N_u1 = -static_cast<complex_t>(Nmat[0][3]);
    param.N_u2 = -static_cast<complex_t>(Nmat[1][3]);
    param.N_u3 = -static_cast<complex_t>(Nmat[2][3]);
    param.N_u4 = -static_cast<complex_t>(Nmat[3][3]);

    param.U_Wm1 = static_cast<complex_t>(U[0][0]);
    param.U_d1 = static_cast<complex_t>(U[0][1]);
    param.U_Wm2 = static_cast<complex_t>(U[1][0]);
    param.U_d2 = static_cast<complex_t>(U[1][1]);

    param.V_Wp1 = static_cast<complex_t>(V[0][0]);
    param.V_u1 = static_cast<complex_t>(V[0][1]);
    param.V_Wp2 = static_cast<complex_t>(V[1][0]);
    param.V_u2 = static_cast<complex_t>(V[1][1]);

    param.U_sb_00 = static_cast<complex_t>(SbottomMix[0][0]);
    param.U_sb_01 = static_cast<complex_t>(SbottomMix[0][1]);
    param.U_sb_10 = static_cast<complex_t>(SbottomMix[1][0]);
    param.U_sb_11 = static_cast<complex_t>(SbottomMix[1][1]);

    param.U_st_00 = static_cast<complex_t>(StopMix[0][0]);
    param.U_st_01 = static_cast<complex_t>(StopMix[0][1]);
    param.U_st_10 = static_cast<complex_t>(StopMix[1][0]);
    param.U_st_11 = static_cast<complex_t>(StopMix[1][1]);

    param.U_stau_00 = -static_cast<complex_t>(StauMix[0][0]);
    param.U_stau_01 = -static_cast<complex_t>(StauMix[0][1]);
    param.U_stau_10 = -static_cast<complex_t>(StauMix[1][0]);
    param.U_stau_11 = -static_cast<complex_t>(StauMix[1][1]);

    param.alpha = static_cast<real_t>(data.getValue("ALPHA", 0).value());

    param.mu_h = get_value(data, "HMIX", 1);
    //     const double tanbetaQ = get_value(data, "HMIX", 2);
    //
    //     const double gprimeQ = get_value(data, "GAUGE", 1);
    //     const double gQ = get_value(data, "GAUGE", 2);
    //
    //     const double Au = get_value_mat(data, "AU", 1);
    //     const double Ac = get_value_mat(data, "AU", 2);
    param.A_t = get_value_mat(data, "AU", 3);

    //     const double Ad = get_value_mat(data, "AD", 1);
    //     const double As = get_value_mat(data, "AD", 2);
    param.A_b = get_value_mat(data, "AD", 3);

    //     const double Ae = get_value_mat(data, "AE", 1);
    //     const double Amu= get_value_mat(data, "AE", 2);
    param.A_tau = get_value_mat(data, "AE", 3);

    param.M_1 = get_value(data, "MSOFT", 1);
    param.M_2 = get_value(data, "MSOFT", 2);
    //     const double M3 = get_value(data, "MSOFT", 3);
    param.M_eL = get_value(data, "MSOFT", 31);
    //     const double MmuL = get_value(data, "MSOFT", 32);
    param.M_tauL = get_value(data, "MSOFT", 33);
    param.M_eR = get_value(data, "MSOFT", 34);
    //     const double MmuR = get_value(data, "MSOFT", 35);
    param.M_tauR = get_value(data, "MSOFT", 36);
    param.M_q1L = get_value(data, "MSOFT", 41);
    //     const double MqL2 = get_value(data, "MSOFT", 42);
    param.M_q3L = get_value(data, "MSOFT", 43);
    param.M_quR = get_value(data, "MSOFT", 44);
    //     const double McR = get_value(data, "MSOFT", 45);
    param.M_qtR = get_value(data, "MSOFT", 46);
    param.M_qdR = get_value(data, "MSOFT", 47);
    //     const double MsR = get_value(data, "MSOFT", 48);
    param.M_qbR = get_value(data, "MSOFT", 49);
  }

  void ReadWidths(Param_t& param, const std::string& namefile)
  {
    auto data = mty::lha::Reader::readFile(namefile);

#ifdef DEBUG
    std::cout << "Entering reading widths function\n";
#endif
#ifdef WIDTH_h
    param.Gamma_h = get_value(data, "ALLDECAYS", 25, width_h(param).real());
#endif
#ifdef WIDTH_TOP
    param.Gamma_top = get_value(data, "ALLDECAYS", 6, width_t(param).real());
#endif
#ifdef WIDTH_H0
    param.Gamma_H0 = get_value(data, "ALLDECAYS", 35, width_H0(param).real());
#endif
#ifdef WIDTH_A0
    param.Gamma_A0 = get_value(data, "ALLDECAYS", 36, width_A0(param).real());
#endif
#ifdef WIDTH_Hp
    param.Gamma_Hp = get_value(data, "ALLDECAYS", 37, width_Hp(param).real());
#endif

    param.Gamma_c1 = 0.;
    param.Gamma_c2 = 0.;
    param.Gamma_o1 = 0.;
    param.Gamma_o2 = 0.;
    param.Gamma_o3 = 0.;
    param.Gamma_o4 = 0.;

    param.Gamma_t1 = 0.;
    param.Gamma_t2 = 0.;
    param.Gamma_b1 = 0.;
    param.Gamma_b2 = 0.;
    param.Gamma_tau1 = 0.;
    param.Gamma_tau2 = 0.;
    param.Gamma_ul = 0.;
    param.Gamma_ur = 0.;
    param.Gamma_dl = 0.;
    param.Gamma_dr = 0.;
    param.Gamma_el = 0.;
    param.Gamma_er = 0.;
    param.Gamma_cl = 0.;
    param.Gamma_cr = 0.;
    param.Gamma_sl = 0.;
    param.Gamma_sr = 0.;
    param.Gamma_ml = 0.;
    param.Gamma_mr = 0.;
    param.Gamma_nuel = 0.;
    param.Gamma_numl = 0.;
    param.Gamma_nutaul = 0.;

    // ----  Assigning widths ----
    for (auto part = 1; part <= corr::TOTAL_PARTICLES; part++)
    {
      param.widths_vector[part] = 0.;
    }

#ifdef WIDTH_Z
    param.widths_vector[corr::Z] = param.Gamma_Z;
#endif
#ifdef WIDTH_W
    param.widths_vector[corr::W] = param.Gamma_W;
#endif
#ifdef WIDTH_A0
    param.widths_vector[corr::A0] = param.Gamma_A0;
#endif
#ifdef WIDTH_h
    param.widths_vector[corr::h] = param.Gamma_h;
#endif
#ifdef WIDTH_H0
    param.widths_vector[corr::H0] = param.Gamma_H0;
#endif
#ifdef WIDTH_TOP
    param.widths_vector[corr::t] = param.Gamma_top;
#endif
#ifdef WIDTH_Hp
    param.widths_vector[corr::Hp] = param.Gamma_Hp;
#endif

    param.arewidthsassigned = true;

#ifdef DEBUG
    param.Print();
//     scanf("%*c");
#endif
  }

  Param_t ReadLHA(const std::string& name)
  {
    const std::string check_existance = "[ -f " + name + " ]";
    if (std::system(check_existance.c_str()) != 0)
    {
      std::cerr << name << " does not exist\n";
      exit(1);
    }
    Param_t param;
    param.pi = M_PI;
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
    std::cout << "Reading widths\n";
#endif

    ReadWidths(param, name);

#ifdef DEBUG
    param.Print();
    std::cout << "Fixing trilinears\n";
#endif
    input_manip::FixTrilinears(param);

#ifdef DEBUG
    param.Print();
    std::cout << "Exiting read function\n";
#endif
    return param;
  }

  void createwidthsfile(const Param_t& param, const std::string& namefile)
  {
    FILE* fptr = fopen(namefile.c_str(), "w");
    if (fptr == nullptr)
      return;

#ifdef DEBUG
    std::cout << "Printing widths in both formats\n";
    fprintf(fptr, "\nDECAY    %d    %.6E    # %s decay\n", 23, param.widths_vector[corr::Z].get(), "Z");
    fprintf(fptr, "DECAY    %d    %.6E    # %s decay\n", 24, param.widths_vector[corr::W].get(), "W");
    fprintf(fptr, "DECAY    %d    %.6E    # %s decay\n", 37, param.widths_vector[corr::Hp].get(), "H+");
    fprintf(fptr, "DECAY    %d    %.6E    # %s decay\n", 25, param.widths_vector[corr::h].get(), "h");
    fprintf(fptr, "DECAY    %d    %.6E    # %s decay\n", 35, param.widths_vector[corr::H0].get(), "H0");
    fprintf(fptr, "DECAY    %d    %.6E    # %s decay\n", 36, param.widths_vector[corr::A0].get(), "A0");
    fprintf(fptr, "DECAY    %d    %.6E    # %s decay\n", 6, param.widths_vector[corr::t].get(), "top");

    fprintf(fptr, "\nBLOCK ALLDECAYS\n");
    fprintf(fptr, "   %d    %.6E    # %s decay\n", 23, param.widths_vector[corr::Z].get(), "Z");
    fprintf(fptr, "   %d    %.6E    # %s decay\n", 24, param.widths_vector[corr::W].get(), "W");
    fprintf(fptr, "   %d    %.6E    # %s decay\n", 37, param.widths_vector[corr::Hp].get(), "H+");
    fprintf(fptr, "   %d    %.6E    # %s decay\n", 25, param.widths_vector[corr::h].get(), "h");
    fprintf(fptr, "   %d    %.6E    # %s decay\n", 35, param.widths_vector[corr::H0].get(), "H0");
    fprintf(fptr, "   %d    %.6E    # %s decay\n", 36, param.widths_vector[corr::A0].get(), "A0");
    fprintf(fptr, "   %d    %.6E    # %s decay\n", 6, param.widths_vector[corr::t].get(), "top");

    fclose(fptr);
#endif
  }
} // namespace mssm2to2::readmodule
