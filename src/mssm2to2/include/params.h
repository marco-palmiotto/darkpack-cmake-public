#ifndef CSL_LIB_PARAM_H_INCLUDED
#define CSL_LIB_PARAM_H_INCLUDED

#include <map>
#include <array>
#include "common.h"
#include "libcomplexop.h"
#include "csl/initSanitizer.h"

namespace mssm2to2 {

struct param_t {

    ///////////////////////////////////////
    // Elementary parameters to be defined 
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> pi { "pi" };
    csl::InitSanitizer<real_t> A_b { "A_b" };
    csl::InitSanitizer<real_t> A_t { "A_t" };
    csl::InitSanitizer<real_t> g_s { "g_s" };
    csl::InitSanitizer<real_t> M_eL { "M_eL" };
    csl::InitSanitizer<real_t> M_eR { "M_eR" };
    csl::InitSanitizer<real_t> V_cb { "V_cb" };
    csl::InitSanitizer<real_t> V_ud { "V_ud" };
    csl::InitSanitizer<real_t> V_us { "V_us" };
    csl::InitSanitizer<real_t> beta { "beta" };
    csl::InitSanitizer<real_t> e_em { "e_em" };
    csl::InitSanitizer<real_t> mu_h { "mu_h" };
    csl::InitSanitizer<real_t> s_12 { "s_12" };
    csl::InitSanitizer<real_t> s_13 { "s_13" };
    csl::InitSanitizer<real_t> s_14 { "s_14" };
    csl::InitSanitizer<real_t> s_23 { "s_23" };
    csl::InitSanitizer<real_t> s_24 { "s_24" };
    csl::InitSanitizer<real_t> s_34 { "s_34" };
    csl::InitSanitizer<real_t> A_tau { "A_tau" };
    csl::InitSanitizer<real_t> M_q1L { "M_q1L" };
    csl::InitSanitizer<real_t> M_q3L { "M_q3L" };
    csl::InitSanitizer<real_t> M_qbR { "M_qbR" };
    csl::InitSanitizer<real_t> M_qdR { "M_qdR" };
    csl::InitSanitizer<real_t> M_quR { "M_quR" };
    csl::InitSanitizer<real_t> alpha { "alpha" };
    csl::InitSanitizer<real_t> M_tauL { "M_tauL" };
    csl::InitSanitizer<real_t> Gamma_W { "Gamma_W" };
    csl::InitSanitizer<real_t> Gamma_Z { "Gamma_Z" };
    csl::InitSanitizer<real_t> Gamma_h { "Gamma_h" };
    csl::InitSanitizer<real_t> theta_W { "theta_W" };
    csl::InitSanitizer<real_t> Gamma_A0 { "Gamma_A0" };
    csl::InitSanitizer<real_t> Gamma_H0 { "Gamma_H0" };
    csl::InitSanitizer<real_t> Gamma_Hp { "Gamma_Hp" };
    csl::InitSanitizer<real_t> Gamma_b1 { "Gamma_b1" };
    csl::InitSanitizer<real_t> Gamma_b2 { "Gamma_b2" };
    csl::InitSanitizer<real_t> Gamma_c1 { "Gamma_c1" };
    csl::InitSanitizer<real_t> Gamma_c2 { "Gamma_c2" };
    csl::InitSanitizer<real_t> Gamma_cl { "Gamma_cl" };
    csl::InitSanitizer<real_t> Gamma_cr { "Gamma_cr" };
    csl::InitSanitizer<real_t> Gamma_dl { "Gamma_dl" };
    csl::InitSanitizer<real_t> Gamma_dr { "Gamma_dr" };
    csl::InitSanitizer<real_t> Gamma_el { "Gamma_el" };
    csl::InitSanitizer<real_t> Gamma_er { "Gamma_er" };
    csl::InitSanitizer<real_t> Gamma_ml { "Gamma_ml" };
    csl::InitSanitizer<real_t> Gamma_mr { "Gamma_mr" };
    csl::InitSanitizer<real_t> Gamma_o1 { "Gamma_o1" };
    csl::InitSanitizer<real_t> Gamma_o2 { "Gamma_o2" };
    csl::InitSanitizer<real_t> Gamma_o3 { "Gamma_o3" };
    csl::InitSanitizer<real_t> Gamma_o4 { "Gamma_o4" };
    csl::InitSanitizer<real_t> Gamma_sl { "Gamma_sl" };
    csl::InitSanitizer<real_t> Gamma_sr { "Gamma_sr" };
    csl::InitSanitizer<real_t> Gamma_t1 { "Gamma_t1" };
    csl::InitSanitizer<real_t> Gamma_t2 { "Gamma_t2" };
    csl::InitSanitizer<real_t> Gamma_ul { "Gamma_ul" };
    csl::InitSanitizer<real_t> Gamma_ur { "Gamma_ur" };
    csl::InitSanitizer<real_t> V_ub_mod { "V_ub_mod" };
    csl::InitSanitizer<real_t> reg_prop { "reg_prop" };
    csl::InitSanitizer<real_t> Gamma_top { "Gamma_top" };
    csl::InitSanitizer<real_t> Gamma_nuel { "Gamma_nuel" };
    csl::InitSanitizer<real_t> Gamma_numl { "Gamma_numl" };
    csl::InitSanitizer<real_t> Gamma_tau1 { "Gamma_tau1" };
    csl::InitSanitizer<real_t> Gamma_tau2 { "Gamma_tau2" };
    csl::InitSanitizer<real_t> delta_wolf { "delta_wolf" };
    csl::InitSanitizer<real_t> Gamma_nutaul { "Gamma_nutaul" };
    csl::InitSanitizer<complex_t> M_1 { "M_1" };
    csl::InitSanitizer<complex_t> M_2 { "M_2" };
    csl::InitSanitizer<complex_t> V_cd { "V_cd" };
    csl::InitSanitizer<complex_t> V_cs { "V_cs" };
    csl::InitSanitizer<complex_t> V_tb { "V_tb" };
    csl::InitSanitizer<complex_t> V_td { "V_td" };
    csl::InitSanitizer<complex_t> V_ts { "V_ts" };
    csl::InitSanitizer<complex_t> M_qtR { "M_qtR" };
    csl::InitSanitizer<complex_t> M_tauR { "M_tauR" };


    ///////////////////////////////////////
    // Parameters functions of others  
    // through diagonalization or mass 
    // expressions, see updateSpectrum()  
    // in global.h or set them by hand  
    // 
    // And other default parameters  
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> m_W { "m_W" };
    csl::InitSanitizer<real_t> m_Z { "m_Z" };
    csl::InitSanitizer<real_t> m_b { "m_b" };
    csl::InitSanitizer<real_t> m_c { "m_c" };
    csl::InitSanitizer<real_t> m_d { "m_d" };
    csl::InitSanitizer<real_t> m_e { "m_e" };
    csl::InitSanitizer<real_t> m_h { "m_h" };
    csl::InitSanitizer<real_t> m_s { "m_s" };
    csl::InitSanitizer<real_t> m_t { "m_t" };
    csl::InitSanitizer<real_t> m_u { "m_u" };
    csl::InitSanitizer<real_t> m_A0 { "m_A0" };
    csl::InitSanitizer<real_t> m_H0 { "m_H0" };
    csl::InitSanitizer<real_t> m_Hp { "m_Hp" };
    csl::InitSanitizer<real_t> m_mu { "m_mu" };
    csl::InitSanitizer<real_t> m_sG { "m_sG" };
    csl::InitSanitizer<real_t> m_C_1 { "m_C_1" };
    csl::InitSanitizer<real_t> m_C_2 { "m_C_2" };
    csl::InitSanitizer<real_t> m_N_1 { "m_N_1" };
    csl::InitSanitizer<real_t> m_N_2 { "m_N_2" };
    csl::InitSanitizer<real_t> m_N_3 { "m_N_3" };
    csl::InitSanitizer<real_t> m_N_4 { "m_N_4" };
    csl::InitSanitizer<real_t> m_tau { "m_tau" };
    csl::InitSanitizer<real_t> m_sb_1 { "m_sb_1" };
    csl::InitSanitizer<real_t> m_sb_2 { "m_sb_2" };
    csl::InitSanitizer<real_t> m_sc_L { "m_sc_L" };
    csl::InitSanitizer<real_t> m_sc_R { "m_sc_R" };
    csl::InitSanitizer<real_t> m_sd_L { "m_sd_L" };
    csl::InitSanitizer<real_t> m_sd_R { "m_sd_R" };
    csl::InitSanitizer<real_t> m_se_L { "m_se_L" };
    csl::InitSanitizer<real_t> m_se_R { "m_se_R" };
    csl::InitSanitizer<real_t> m_ss_L { "m_ss_L" };
    csl::InitSanitizer<real_t> m_ss_R { "m_ss_R" };
    csl::InitSanitizer<real_t> m_st_1 { "m_st_1" };
    csl::InitSanitizer<real_t> m_st_2 { "m_st_2" };
    csl::InitSanitizer<real_t> m_su_L { "m_su_L" };
    csl::InitSanitizer<real_t> m_su_R { "m_su_R" };
    csl::InitSanitizer<real_t> m_smu_L { "m_smu_L" };
    csl::InitSanitizer<real_t> m_smu_R { "m_smu_R" };
    csl::InitSanitizer<real_t> m_snu_e { "m_snu_e" };
    csl::InitSanitizer<real_t> m_snu_mu { "m_snu_mu" };
    csl::InitSanitizer<real_t> m_stau_1 { "m_stau_1" };
    csl::InitSanitizer<real_t> m_stau_2 { "m_stau_2" };
    csl::InitSanitizer<real_t> m_snu_tau { "m_snu_tau" };
    csl::InitSanitizer<complex_t> N_B1 { "N_B1" };
    csl::InitSanitizer<complex_t> N_B2 { "N_B2" };
    csl::InitSanitizer<complex_t> N_B3 { "N_B3" };
    csl::InitSanitizer<complex_t> N_B4 { "N_B4" };
    csl::InitSanitizer<complex_t> N_W1 { "N_W1" };
    csl::InitSanitizer<complex_t> N_W2 { "N_W2" };
    csl::InitSanitizer<complex_t> N_W3 { "N_W3" };
    csl::InitSanitizer<complex_t> N_W4 { "N_W4" };
    csl::InitSanitizer<complex_t> N_d1 { "N_d1" };
    csl::InitSanitizer<complex_t> N_d2 { "N_d2" };
    csl::InitSanitizer<complex_t> N_d3 { "N_d3" };
    csl::InitSanitizer<complex_t> N_d4 { "N_d4" };
    csl::InitSanitizer<complex_t> N_u1 { "N_u1" };
    csl::InitSanitizer<complex_t> N_u2 { "N_u2" };
    csl::InitSanitizer<complex_t> N_u3 { "N_u3" };
    csl::InitSanitizer<complex_t> N_u4 { "N_u4" };
    csl::InitSanitizer<complex_t> U_d1 { "U_d1" };
    csl::InitSanitizer<complex_t> U_d2 { "U_d2" };
    csl::InitSanitizer<complex_t> V_u1 { "V_u1" };
    csl::InitSanitizer<complex_t> V_u2 { "V_u2" };
    csl::InitSanitizer<complex_t> U_Wm1 { "U_Wm1" };
    csl::InitSanitizer<complex_t> U_Wm2 { "U_Wm2" };
    csl::InitSanitizer<complex_t> V_Wp1 { "V_Wp1" };
    csl::InitSanitizer<complex_t> V_Wp2 { "V_Wp2" };
    csl::InitSanitizer<complex_t> U_sb_00 { "U_sb_00" };
    csl::InitSanitizer<complex_t> U_sb_01 { "U_sb_01" };
    csl::InitSanitizer<complex_t> U_sb_10 { "U_sb_10" };
    csl::InitSanitizer<complex_t> U_sb_11 { "U_sb_11" };
    csl::InitSanitizer<complex_t> U_st_00 { "U_st_00" };
    csl::InitSanitizer<complex_t> U_st_01 { "U_st_01" };
    csl::InitSanitizer<complex_t> U_st_10 { "U_st_10" };
    csl::InitSanitizer<complex_t> U_st_11 { "U_st_11" };
    csl::InitSanitizer<complex_t> U_stau_00 { "U_stau_00" };
    csl::InitSanitizer<complex_t> U_stau_01 { "U_stau_01" };
    csl::InitSanitizer<complex_t> U_stau_10 { "U_stau_10" };
    csl::InitSanitizer<complex_t> U_stau_11 { "U_stau_11" };

    void reset()
    {
        using real_params = std::array<csl::InitSanitizer<real_t>*, 107>;
        using complex_params = std::array<csl::InitSanitizer<complex_t>*, 45>;

        for (auto &par : real_params{
                &pi, &A_b, &A_t, &g_s, &M_eL, 
                &M_eR, &V_cb, &V_ud, &V_us, &beta, &e_em, 
                &mu_h, &s_12, &s_13, &s_14, &s_23, &s_24, 
                &s_34, &A_tau, &M_q1L, &M_q3L, &M_qbR, &M_qdR, 
                &M_quR, &alpha, &M_tauL, &Gamma_W, &Gamma_Z, &Gamma_h, 
                &theta_W, &Gamma_A0, &Gamma_H0, &Gamma_Hp, &Gamma_b1, &Gamma_b2, 
                &Gamma_c1, &Gamma_c2, &Gamma_cl, &Gamma_cr, &Gamma_dl, &Gamma_dr, 
                &Gamma_el, &Gamma_er, &Gamma_ml, &Gamma_mr, &Gamma_o1, &Gamma_o2, 
                &Gamma_o3, &Gamma_o4, &Gamma_sl, &Gamma_sr, &Gamma_t1, &Gamma_t2, 
                &Gamma_ul, &Gamma_ur, &V_ub_mod, &reg_prop, &Gamma_top, &Gamma_nuel, 
                &Gamma_numl, &Gamma_tau1, &Gamma_tau2, &delta_wolf, &Gamma_nutaul, &m_W, 
                &m_Z, &m_b, &m_c, &m_d, &m_e, &m_h, 
                &m_s, &m_t, &m_u, &m_A0, &m_H0, &m_Hp, 
                &m_mu, &m_sG, &m_C_1, &m_C_2, &m_N_1, &m_N_2, 
                &m_N_3, &m_N_4, &m_tau, &m_sb_1, &m_sb_2, &m_sc_L, 
                &m_sc_R, &m_sd_L, &m_sd_R, &m_se_L, &m_se_R, &m_ss_L, 
                &m_ss_R, &m_st_1, &m_st_2, &m_su_L, &m_su_R, &m_smu_L, 
                &m_smu_R, &m_snu_e, &m_snu_mu, &m_stau_1, &m_stau_2, &m_snu_tau, })
        {
            par->reset();
        }

        for (auto &par : complex_params{
                &M_1, &M_2, &V_cd, &V_cs, &V_tb, 
                &V_td, &V_ts, &M_qtR, &M_tauR, &N_B1, &N_B2, 
                &N_B3, &N_B4, &N_W1, &N_W2, &N_W3, &N_W4, 
                &N_d1, &N_d2, &N_d3, &N_d4, &N_u1, &N_u2, 
                &N_u3, &N_u4, &U_d1, &U_d2, &V_u1, &V_u2, 
                &U_Wm1, &U_Wm2, &V_Wp1, &V_Wp2, &U_sb_00, &U_sb_01, 
                &U_sb_10, &U_sb_11, &U_st_00, &U_st_01, &U_st_10, &U_st_11, 
                &U_stau_00, &U_stau_01, &U_stau_10, &U_stau_11, })
        {
            par->reset();
        }
    }

    void print(std::ostream &out = std::cout) const
    {
        using real_params = std::array<csl::InitSanitizer<real_t> const*, 107>;
        using complex_params = std::array<csl::InitSanitizer<complex_t> const*, 45>;

        out << "param_t struct:\n";
        out << "Real parameters\n";
        for (auto const &par : real_params{
                &pi, &A_b, &A_t, &g_s, &M_eL, 
                &M_eR, &V_cb, &V_ud, &V_us, &beta, &e_em, 
                &mu_h, &s_12, &s_13, &s_14, &s_23, &s_24, 
                &s_34, &A_tau, &M_q1L, &M_q3L, &M_qbR, &M_qdR, 
                &M_quR, &alpha, &M_tauL, &Gamma_W, &Gamma_Z, &Gamma_h, 
                &theta_W, &Gamma_A0, &Gamma_H0, &Gamma_Hp, &Gamma_b1, &Gamma_b2, 
                &Gamma_c1, &Gamma_c2, &Gamma_cl, &Gamma_cr, &Gamma_dl, &Gamma_dr, 
                &Gamma_el, &Gamma_er, &Gamma_ml, &Gamma_mr, &Gamma_o1, &Gamma_o2, 
                &Gamma_o3, &Gamma_o4, &Gamma_sl, &Gamma_sr, &Gamma_t1, &Gamma_t2, 
                &Gamma_ul, &Gamma_ur, &V_ub_mod, &reg_prop, &Gamma_top, &Gamma_nuel, 
                &Gamma_numl, &Gamma_tau1, &Gamma_tau2, &delta_wolf, &Gamma_nutaul, &m_W, 
                &m_Z, &m_b, &m_c, &m_d, &m_e, &m_h, 
                &m_s, &m_t, &m_u, &m_A0, &m_H0, &m_Hp, 
                &m_mu, &m_sG, &m_C_1, &m_C_2, &m_N_1, &m_N_2, 
                &m_N_3, &m_N_4, &m_tau, &m_sb_1, &m_sb_2, &m_sc_L, 
                &m_sc_R, &m_sd_L, &m_sd_R, &m_se_L, &m_se_R, &m_ss_L, 
                &m_ss_R, &m_st_1, &m_st_2, &m_su_L, &m_su_R, &m_smu_L, 
                &m_smu_R, &m_snu_e, &m_snu_mu, &m_stau_1, &m_stau_2, &m_snu_tau, })
        {
            out << "  -> " << par->getName() << ": ";
            if (par->hasValue()) {
                out << static_cast<double>(par->get()) << '\n';
            }
            else {
                out << "uninitialized" << '\n';
            }
        }

        out << "Complex parameters\n";
        for (auto const &par : complex_params{
                &M_1, &M_2, &V_cd, &V_cs, &V_tb, 
                &V_td, &V_ts, &M_qtR, &M_tauR, &N_B1, &N_B2, 
                &N_B3, &N_B4, &N_W1, &N_W2, &N_W3, &N_W4, 
                &N_d1, &N_d2, &N_d3, &N_d4, &N_u1, &N_u2, 
                &N_u3, &N_u4, &U_d1, &U_d2, &V_u1, &V_u2, 
                &U_Wm1, &U_Wm2, &V_Wp1, &V_Wp2, &U_sb_00, &U_sb_01, 
                &U_sb_10, &U_sb_11, &U_st_00, &U_st_01, &U_st_10, &U_st_11, 
                &U_stau_00, &U_stau_01, &U_stau_10, &U_stau_11, })
        {
            out << "  -> " << par->getName() << ": ";
            if (par->hasValue()) {
                out << static_cast<double>(MTY_REAL(par->get()));
                out << " + i*" << static_cast<double>(MTY_IMAG(par->get())) << '\n';
            }
            else {
                out << "uninitialized" << '\n';
            }
        }
        out << "\n";
    }

    std::map<std::string, csl::InitSanitizer<real_t>*> realParams {
        {"pi", &pi},
        {"A_b", &A_b},
        {"A_t", &A_t},
        {"g_s", &g_s},
        {"M_eL", &M_eL},
        {"M_eR", &M_eR},
        {"V_cb", &V_cb},
        {"V_ud", &V_ud},
        {"V_us", &V_us},
        {"beta", &beta},
        {"e_em", &e_em},
        {"mu_h", &mu_h},
        {"s_12", &s_12},
        {"s_13", &s_13},
        {"s_14", &s_14},
        {"s_23", &s_23},
        {"s_24", &s_24},
        {"s_34", &s_34},
        {"A_tau", &A_tau},
        {"M_q1L", &M_q1L},
        {"M_q3L", &M_q3L},
        {"M_qbR", &M_qbR},
        {"M_qdR", &M_qdR},
        {"M_quR", &M_quR},
        {"alpha", &alpha},
        {"M_tauL", &M_tauL},
        {"Gamma_W", &Gamma_W},
        {"Gamma_Z", &Gamma_Z},
        {"Gamma_h", &Gamma_h},
        {"theta_W", &theta_W},
        {"Gamma_A0", &Gamma_A0},
        {"Gamma_H0", &Gamma_H0},
        {"Gamma_Hp", &Gamma_Hp},
        {"Gamma_b1", &Gamma_b1},
        {"Gamma_b2", &Gamma_b2},
        {"Gamma_c1", &Gamma_c1},
        {"Gamma_c2", &Gamma_c2},
        {"Gamma_cl", &Gamma_cl},
        {"Gamma_cr", &Gamma_cr},
        {"Gamma_dl", &Gamma_dl},
        {"Gamma_dr", &Gamma_dr},
        {"Gamma_el", &Gamma_el},
        {"Gamma_er", &Gamma_er},
        {"Gamma_ml", &Gamma_ml},
        {"Gamma_mr", &Gamma_mr},
        {"Gamma_o1", &Gamma_o1},
        {"Gamma_o2", &Gamma_o2},
        {"Gamma_o3", &Gamma_o3},
        {"Gamma_o4", &Gamma_o4},
        {"Gamma_sl", &Gamma_sl},
        {"Gamma_sr", &Gamma_sr},
        {"Gamma_t1", &Gamma_t1},
        {"Gamma_t2", &Gamma_t2},
        {"Gamma_ul", &Gamma_ul},
        {"Gamma_ur", &Gamma_ur},
        {"V_ub_mod", &V_ub_mod},
        {"reg_prop", &reg_prop},
        {"Gamma_top", &Gamma_top},
        {"Gamma_nuel", &Gamma_nuel},
        {"Gamma_numl", &Gamma_numl},
        {"Gamma_tau1", &Gamma_tau1},
        {"Gamma_tau2", &Gamma_tau2},
        {"delta_wolf", &delta_wolf},
        {"Gamma_nutaul", &Gamma_nutaul},
        {"m_W", &m_W},
        {"m_Z", &m_Z},
        {"m_b", &m_b},
        {"m_c", &m_c},
        {"m_d", &m_d},
        {"m_e", &m_e},
        {"m_h", &m_h},
        {"m_s", &m_s},
        {"m_t", &m_t},
        {"m_u", &m_u},
        {"m_A0", &m_A0},
        {"m_H0", &m_H0},
        {"m_Hp", &m_Hp},
        {"m_mu", &m_mu},
        {"m_sG", &m_sG},
        {"m_C_1", &m_C_1},
        {"m_C_2", &m_C_2},
        {"m_N_1", &m_N_1},
        {"m_N_2", &m_N_2},
        {"m_N_3", &m_N_3},
        {"m_N_4", &m_N_4},
        {"m_tau", &m_tau},
        {"m_sb_1", &m_sb_1},
        {"m_sb_2", &m_sb_2},
        {"m_sc_L", &m_sc_L},
        {"m_sc_R", &m_sc_R},
        {"m_sd_L", &m_sd_L},
        {"m_sd_R", &m_sd_R},
        {"m_se_L", &m_se_L},
        {"m_se_R", &m_se_R},
        {"m_ss_L", &m_ss_L},
        {"m_ss_R", &m_ss_R},
        {"m_st_1", &m_st_1},
        {"m_st_2", &m_st_2},
        {"m_su_L", &m_su_L},
        {"m_su_R", &m_su_R},
        {"m_smu_L", &m_smu_L},
        {"m_smu_R", &m_smu_R},
        {"m_snu_e", &m_snu_e},
        {"m_snu_mu", &m_snu_mu},
        {"m_stau_1", &m_stau_1},
        {"m_stau_2", &m_stau_2},
        {"m_snu_tau", &m_snu_tau},
    };

    std::map<std::string, csl::InitSanitizer<complex_t>*> complexParams {
        {"M_1", &M_1},
        {"M_2", &M_2},
        {"V_cd", &V_cd},
        {"V_cs", &V_cs},
        {"V_tb", &V_tb},
        {"V_td", &V_td},
        {"V_ts", &V_ts},
        {"M_qtR", &M_qtR},
        {"M_tauR", &M_tauR},
        {"N_B1", &N_B1},
        {"N_B2", &N_B2},
        {"N_B3", &N_B3},
        {"N_B4", &N_B4},
        {"N_W1", &N_W1},
        {"N_W2", &N_W2},
        {"N_W3", &N_W3},
        {"N_W4", &N_W4},
        {"N_d1", &N_d1},
        {"N_d2", &N_d2},
        {"N_d3", &N_d3},
        {"N_d4", &N_d4},
        {"N_u1", &N_u1},
        {"N_u2", &N_u2},
        {"N_u3", &N_u3},
        {"N_u4", &N_u4},
        {"U_d1", &U_d1},
        {"U_d2", &U_d2},
        {"V_u1", &V_u1},
        {"V_u2", &V_u2},
        {"U_Wm1", &U_Wm1},
        {"U_Wm2", &U_Wm2},
        {"V_Wp1", &V_Wp1},
        {"V_Wp2", &V_Wp2},
        {"U_sb_00", &U_sb_00},
        {"U_sb_01", &U_sb_01},
        {"U_sb_10", &U_sb_10},
        {"U_sb_11", &U_sb_11},
        {"U_st_00", &U_st_00},
        {"U_st_01", &U_st_01},
        {"U_st_10", &U_st_10},
        {"U_st_11", &U_st_11},
        {"U_stau_00", &U_stau_00},
        {"U_stau_01", &U_stau_01},
        {"U_stau_10", &U_stau_10},
        {"U_stau_11", &U_stau_11},
    };

};


}

#endif
