#include "leshouchesfromsuperiso.hpp"
// #define DEBUG

inline double sq(double& x) { return x * x; }

namespace mssm2to2::superisosupport
{

  int InitInterfaceStruct(const struct parameters* const param, struct param_t& input)
  {
#ifdef DEBUG
    std::cout << "Called InitInterfaceStruct\n  1. Copying values from the SuperIso param structure to the param_t "
                 "structure\n";
#endif
    input.e_em = std::sqrt(4 * M_PI / 137.035999710);
    input.g_s = std::sqrt(4 * M_PI * param->alphas_MZ);

    input.m_e = param->mass_e;
    input.m_mu = param->mass_mu;

    input.m_u = param->mass_u;
    input.m_d = param->mass_d;
    input.m_c = param->mass_c;

    input.m_b = param->mass_b; // mb(mb)
    input.m_s = param->mass_s;
    input.m_t = param->mass_top_pole; // mass_top_pole is going to be fixed with the 1st running
    input.theta_W = std::acos(param->mass_W / param->mass_Z);
    input.delta_wolf = std::atan2(param->CKM_etabar, param->CKM_rhobar); // Angle between eta_bar and rho_bar
    input.V_ud = (static_cast<complex_t>(param->Vud)).real();
    input.V_us = (static_cast<complex_t>(param->Vus)).real();
    input.V_ub_mod = std::abs(static_cast<complex_t>(param->Vub));
    input.V_cd = static_cast<complex_t>(param->Vcd);
    input.V_cs = static_cast<complex_t>(param->Vcs);
    input.V_cb = (static_cast<complex_t>(param->Vcb)).real();
    input.V_td = static_cast<complex_t>(param->Vtd);
    input.V_ts = static_cast<complex_t>(param->Vts);
    input.V_tb = static_cast<complex_t>(param->Vtb);
    input.m_tau = param->mass_tau;

    input.m_sG = param->mass_gluino;
    input.mu_h = param->mu_Q; // it's the mu in the Higgs potential
    input.alpha = param->alpha;
    input.beta = std::atan(param->tan_beta);

    input.m_W = param->mass_W;
    input.m_Z = param->mass_Z;

    // Stau masses at Q scale
    input.M_1 = static_cast<complex_t>(param->M1_Q);
    input.M_2 = static_cast<complex_t>(param->M2_Q);

    input.M_eL = param->MeL_Q;
    input.M_eR = param->MeR_Q;

    //       input.M_eL = param->mass_el;
    //       input.M_eR = param->mass_er;

    input.M_tauL = param->MtauL_Q;
    input.M_tauR = static_cast<complex_t>(param->MtauR_Q).real();

    input.M_q1L = param->MqL1_Q;
    input.M_q3L = static_cast<complex_t>(param->MqL3_Q).real();

    input.M_quR = param->MuR_Q;
    input.M_qdR = param->MdR_Q;

    //       input.M_quR = param->mass_upr;
    //       input.M_qdR = param->mass_dnr;

    input.M_qbR = static_cast<complex_t>(param->MbR_Q).real();
    input.M_qtR = static_cast<complex_t>(param->MtR_Q).real();

    // Trilinear couplings
    input.A_b = param->A_b;
    input.A_t = param->A_t;
    input.A_tau = param->A_tau;
    // Numbered suffixes are for mass eigenstates
    input.m_sb_1 = param->mass_b1;
    input.m_sb_2 = param->mass_b2;
    input.m_st_1 = param->mass_t1;
    input.m_st_2 = param->mass_t2;
    input.m_stau_1 = param->mass_tau1;
    input.m_stau_2 = param->mass_tau2;

    // Mixing matrices
    input.N_B1 = static_cast<complex_t>(param->neut_mix[1][1]);
    input.N_B2 = static_cast<complex_t>(param->neut_mix[2][1]);
    input.N_B3 = static_cast<complex_t>(param->neut_mix[3][1]);
    input.N_B4 = static_cast<complex_t>(param->neut_mix[4][1]);
    input.N_W1 = static_cast<complex_t>(param->neut_mix[1][2]);
    input.N_W2 = static_cast<complex_t>(param->neut_mix[2][2]);
    input.N_W3 = static_cast<complex_t>(param->neut_mix[3][2]);
    input.N_W4 = static_cast<complex_t>(param->neut_mix[4][2]);
    input.N_d1 = static_cast<complex_t>(param->neut_mix[1][3]);
    input.N_d2 = static_cast<complex_t>(param->neut_mix[2][3]);
    input.N_d3 = static_cast<complex_t>(param->neut_mix[3][3]);
    input.N_d4 = static_cast<complex_t>(param->neut_mix[4][3]);
    input.N_u1 = static_cast<complex_t>(param->neut_mix[1][4]);
    input.N_u2 = static_cast<complex_t>(param->neut_mix[2][4]);
    input.N_u3 = static_cast<complex_t>(param->neut_mix[3][4]);
    input.N_u4 = static_cast<complex_t>(param->neut_mix[4][4]);

    input.U_Wm1 = static_cast<complex_t>(param->charg_Umix[1][1]);
    input.U_d1 = static_cast<complex_t>(param->charg_Umix[1][2]);
    input.U_Wm2 = static_cast<complex_t>(param->charg_Umix[2][1]);
    input.U_d2 = static_cast<complex_t>(param->charg_Umix[2][2]);

    input.V_Wp1 = static_cast<complex_t>(param->charg_Vmix[1][1]);
    input.V_u1 = static_cast<complex_t>(param->charg_Vmix[1][2]);
    input.V_Wp2 = static_cast<complex_t>(param->charg_Vmix[2][1]);
    input.V_u2 = static_cast<complex_t>(param->charg_Vmix[2][2]);


    input.U_sb_00 = static_cast<complex_t>(param->sbot_mix[1][1]);
    input.U_sb_01 = static_cast<complex_t>(param->sbot_mix[2][1]);
    input.U_sb_10 = static_cast<complex_t>(param->sbot_mix[1][2]);
    input.U_sb_11 = static_cast<complex_t>(param->sbot_mix[2][2]);

    input.U_st_00 = static_cast<complex_t>(param->stop_mix[1][1]);
    input.U_st_01 = static_cast<complex_t>(param->stop_mix[2][1]);
    input.U_st_10 = static_cast<complex_t>(param->stop_mix[1][2]);
    input.U_st_11 = static_cast<complex_t>(param->stop_mix[2][2]);

    input.U_stau_00 = static_cast<complex_t>(param->stau_mix[1][1]);
    input.U_stau_01 = static_cast<complex_t>(param->stau_mix[2][1]);
    input.U_stau_10 = static_cast<complex_t>(param->stau_mix[1][2]);
    input.U_stau_11 = static_cast<complex_t>(param->stau_mix[2][2]);

    // Charginos and Neutralinos's masses
    input.m_C_1 = param->mass_cha1;
    input.m_C_2 = param->mass_cha2;

    input.m_N_1 = param->mass_neut[1];
    input.m_N_2 = param->mass_neut[2];
    input.m_N_3 = param->mass_neut[3];
    input.m_N_4 = param->mass_neut[4];

    input.reg_prop = 0;

    /*
    input.Gamma_Z = param->width_Z;
    input.Gamma_W = param->width_W;

    input.Gamma_Hp = param->width_H;
    input.Gamma_h = param->width_h0;
    input.Gamma_H0 = param->width_H0;
    input.Gamma_A0 = param->width_A0;
    input.Gamma_top = param->width_top;
    */

    input.Gamma_c1 = param->width_c1 * 0.;
    input.Gamma_c2 = param->width_c2 * 0.;
    input.Gamma_o1 = param->width_o1 * 0.;
    input.Gamma_o2 = param->width_o2 * 0.;
    input.Gamma_o3 = param->width_o3 * 0.;
    input.Gamma_o4 = param->width_o4 * 0.;

    input.Gamma_t1 = param->width_t1 * 0.;
    input.Gamma_t2 = param->width_t2 * 0.;
    input.Gamma_b1 = param->width_b1 * 0.;
    input.Gamma_b2 = param->width_b2 * 0.;
    input.Gamma_tau1 = param->width_tau1 * 0.;
    input.Gamma_tau2 = param->width_tau2 * 0.;
    input.Gamma_ul = param->width_ul * 0.;
    input.Gamma_ur = param->width_ur * 0.;
    input.Gamma_dl = param->width_dl * 0.;
    input.Gamma_dr = param->width_dr * 0.;
    input.Gamma_el = param->width_el * 0.;
    input.Gamma_er = param->width_er * 0.;
    input.Gamma_cl = param->width_cl * 0.;
    input.Gamma_cr = param->width_cr * 0.;
    input.Gamma_sl = param->width_sl * 0.;
    input.Gamma_sr = param->width_sr * 0.;
    input.Gamma_ml = param->width_ml * 0.;
    input.Gamma_mr = param->width_mr * 0.;
    input.Gamma_nuel = param->width_nuel * 0.;
    input.Gamma_numl = param->width_numl * 0.;
    input.Gamma_nutaul = param->width_nutaul * 0.;

    input.m_A0 = param->mass_A0;
    input.m_h = param->mass_h0;
    input.m_H0 = param->mass_H0;
    input.m_Hp = param->mass_H;

    input.m_se_L = param->mass_el;
    input.m_smu_L = param->mass_mul;

    input.m_se_R = param->mass_er;
    input.m_smu_R = param->mass_mur;

    input.m_su_L = param->mass_upl;
    input.m_sc_L = param->mass_chl;

    input.m_su_R = param->mass_upr;
    input.m_sc_R = param->mass_chr;

    input.m_sd_L = param->mass_dnl;
    input.m_ss_L = param->mass_stl;

    input.m_sd_R = param->mass_dnr;
    input.m_ss_R = param->mass_str;

    input.m_snu_e = param->mass_nuel;

    input.m_snu_mu = param->mass_numl;

    input.m_snu_tau = param->mass_nutl;
#ifdef DEBUG
    std::cout << "  2. Fixing MSSM trilinear couplings to make it work in MARTY's convention\n";
#endif
    input_manip::FixTrilinears(input);

    return 0;
  }

  int InitInterfaceStruct_Full(const struct parameters* const param, struct Param_t& input)
  {
    InitInterfaceStruct(param, input);
    // Fixing common parameters
    // They are not manipulated, so they follow the LHA convention
    input.Running_scale = -2.;
    input.alpha_str_Mz = param->alphas_MZ;
    input.m_Z_pole = param->mass_Z;
    input.mbmb = param->mass_b;
    input.mtop_pole = param->mass_top_pole;
    input.m_tau_pole = param->mass_tau;
    input.m_el_pole = param->mass_e;
    input.m_mu_pole = param->mass_mu;
    input.m_d2GeV = param->mass_d;
    input.m_u2GeV = param->mass_u;
    input.m_s2GeV = param->mass_s;
    input.m_c_m_c = param->mass_c;
    input.ckm_lambda = param->CKM_lambda;
    input.ckm_A = param->CKM_A;
    input.ckm_rhobar = param->CKM_rhobar;
    input.ckm_etabar = param->CKM_etabar;

    input.AssignMassesVector();
    input.setLightestBSMparticle();

    // ----  Assigning widths ----
    for (auto part = 1; part <= corr::TOTAL_PARTICLES; part++)
    {
      input.widths_vector[part] = 0.;
    }
    /*
    input.widths_vector[corr::Z] = input.Gamma_Z;
    input.widths_vector[corr::W] = input.Gamma_W;
    input.widths_vector[corr::Hp] = input.Gamma_Hp;
    input.widths_vector[corr::h] = input.Gamma_h;
    input.widths_vector[corr::H0] = input.Gamma_H0;
    input.widths_vector[corr::A0] = input.Gamma_A0;
    input.widths_vector[corr::t] = input.Gamma_top;
    */

    input.arewidthsassigned = true;

    return 0;
  }

  struct parameters ReadLHA(struct Param_t& input, const char* name, int* err)
  {
    // Initialising inputs using SuperIso's algorithm
    struct parameters param; // We need a struct parameters to store data in the SupeIso format
    char* filename = (char*)(malloc((1 + strlen(name)) * sizeof(char)));
    if (filename == nullptr)
    {
      *err = 3;
      return param;
    }
    strcpy(filename, name);
    test_slha(filename); // Testing if the .lha file in input is valid

    Init_param(&param); // Doing a parameter initialisation of the SuperIso structure

    if (!Les_Houches_Reader(filename, &param))
    {
      std::cout << "Invalid LHA !\n";
      free(filename);
      *err = 1;
      return param;
    }
    free(filename);

    if (param.THDM_model > 0)
    {
      std::cout << "THDM problem ?\n";
      *err = 2;
      return param;
    }
    if (param.sign_mu == 0.0)
      param.sign_mu = 1;
#ifdef DEBUG
    std::cout << "Setting values for the masses\n";
    std::cout << "Initialising InterfaceStruct\n";
#endif
    input.reset();
#ifdef pi
  #undef pi
#endif
    input.pi = M_PI;
    InitInterfaceStruct_Full(&param, input);
    input.alpha_em_MZ_inv = pdgValue::alpha_em_MZ_inv;
    input.Gfermi = pdgValue::Gfermi;
    *err = 0;
    return param;
  }
} // namespace mssm2to2::superisosupport