#ifndef CSL_LIB_CPARAM_H_INCLUDED
#define CSL_LIB_CPARAM_H_INCLUDED

#include "ccommon.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cparam_s {

    creal_t pi;
    creal_t g_d;
    creal_t g_l;
    creal_t g_u;
    creal_t V_cb;
    creal_t V_tb;
    creal_t V_ud;
    creal_t V_us;
    creal_t e_em;
    creal_t s_12;
    creal_t s_23;
    creal_t s_34;
    creal_t g_chi;
    creal_t theta_W;
    creal_t reg_prop;
    creal_t Gamma_phi;
    ccomplex_t V_cd;
    ccomplex_t V_cs;
    ccomplex_t V_td;
    ccomplex_t V_ts;
    ccomplex_t V_ub;

    creal_t g_s;
    creal_t m_W;
    creal_t m_Z;
    creal_t m_b;
    creal_t m_c;
    creal_t m_d;
    creal_t m_e;
    creal_t m_h;
    creal_t m_s;
    creal_t m_t;
    creal_t m_u;
    creal_t m_mu;
    creal_t m_chi;
    creal_t m_phi;
    creal_t m_tau;
} cparam_t;
#ifdef __cplusplus
} // extern "C" block
#endif


#endif
