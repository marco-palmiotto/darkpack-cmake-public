#ifndef CSL_LIB_PARAM_H_INCLUDED
#define CSL_LIB_PARAM_H_INCLUDED

#include <map>
#include <array>
#include "common.h"
#include "libcomplexop.h"
#include "csl/initSanitizer.h"

namespace scalar2to2 {

struct param_t {

    ///////////////////////////////////////
    // Elementary parameters to be defined 
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> pi { "pi" };
    csl::InitSanitizer<real_t> g_d { "g_d" };
    csl::InitSanitizer<real_t> g_l { "g_l" };
    csl::InitSanitizer<real_t> g_u { "g_u" };
    csl::InitSanitizer<real_t> V_cb { "V_cb" };
    csl::InitSanitizer<real_t> V_tb { "V_tb" };
    csl::InitSanitizer<real_t> V_ud { "V_ud" };
    csl::InitSanitizer<real_t> V_us { "V_us" };
    csl::InitSanitizer<real_t> e_em { "e_em" };
    csl::InitSanitizer<real_t> s_12 { "s_12" };
    csl::InitSanitizer<real_t> s_23 { "s_23" };
    csl::InitSanitizer<real_t> s_34 { "s_34" };
    csl::InitSanitizer<real_t> g_chi { "g_chi" };
    csl::InitSanitizer<real_t> theta_W { "theta_W" };
    csl::InitSanitizer<real_t> reg_prop { "reg_prop" };
    csl::InitSanitizer<real_t> Gamma_phi { "Gamma_phi" };
    csl::InitSanitizer<complex_t> V_cd { "V_cd" };
    csl::InitSanitizer<complex_t> V_cs { "V_cs" };
    csl::InitSanitizer<complex_t> V_td { "V_td" };
    csl::InitSanitizer<complex_t> V_ts { "V_ts" };
    csl::InitSanitizer<complex_t> V_ub { "V_ub" };


    ///////////////////////////////////////
    // Parameters functions of others  
    // through diagonalization or mass 
    // expressions, see updateSpectrum()  
    // in global.h or set them by hand  
    // 
    // And other default parameters  
    ///////////////////////////////////////

    csl::InitSanitizer<real_t> g_s { "g_s" };
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
    csl::InitSanitizer<real_t> m_mu { "m_mu" };
    csl::InitSanitizer<real_t> m_chi { "m_chi" };
    csl::InitSanitizer<real_t> m_phi { "m_phi" };
    csl::InitSanitizer<real_t> m_tau { "m_tau" };

    void reset()
    {
        using real_params = std::array<csl::InitSanitizer<real_t>*, 30>;
        using complex_params = std::array<csl::InitSanitizer<complex_t>*, 5>;

        for (auto &par : real_params{
                &pi, &g_d, &g_l, &g_u, &V_cb, 
                &V_tb, &V_ud, &V_us, &e_em, &s_12, &s_34, 
                &g_chi, &theta_W, &reg_prop, &Gamma_phi, &g_s, &m_W, 
                &m_Z, &m_b, &m_c, &m_d, &m_e, &m_h, 
                &m_s, &m_t, &m_u, &m_mu, &m_chi, &m_phi, 
                &m_tau, })
        {
            par->reset();
        }

        for (auto &par : complex_params{&V_cd, &V_cs, &V_td, &V_ts, &V_ub, })
        {
            par->reset();
        }
    }

    void print(std::ostream &out = std::cout) const
    {
        using real_params = std::array<csl::InitSanitizer<real_t> const*, 30>;
        using complex_params = std::array<csl::InitSanitizer<complex_t> const*, 5>;

        out << "param_t struct:\n";
        out << "Real parameters\n";
        for (auto const &par : real_params{
                &pi, &g_d, &g_l, &g_u, &V_cb, 
                &V_tb, &V_ud, &V_us, &e_em, &s_12, &s_34, 
                &g_chi, &theta_W, &reg_prop, &Gamma_phi, &g_s, &m_W, 
                &m_Z, &m_b, &m_c, &m_d, &m_e, &m_h, 
                &m_s, &m_t, &m_u, &m_mu, &m_chi, &m_phi, 
                &m_tau, })
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
        for (auto const &par : complex_params{&V_cd, &V_cs, &V_td, &V_ts, &V_ub, })
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
        {"g_d", &g_d},
        {"g_l", &g_l},
        {"g_u", &g_u},
        {"V_cb", &V_cb},
        {"V_tb", &V_tb},
        {"V_ud", &V_ud},
        {"V_us", &V_us},
        {"e_em", &e_em},
        {"s_12", &s_12},
        {"s_34", &s_34},
        {"g_chi", &g_chi},
        {"theta_W", &theta_W},
        {"reg_prop", &reg_prop},
        {"Gamma_phi", &Gamma_phi},
        {"g_s", &g_s},
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
        {"m_mu", &m_mu},
        {"m_chi", &m_chi},
        {"m_phi", &m_phi},
        {"m_tau", &m_tau},
    };

    std::map<std::string, csl::InitSanitizer<complex_t>*> complexParams {
        {"V_cd", &V_cd},
        {"V_cs", &V_cs},
        {"V_td", &V_td},
        {"V_ts", &V_ts},
        {"V_ub", &V_ub},
    };

};


}

#endif
