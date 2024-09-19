#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_R_to_h_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_R_to_h_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = pow(m_d, 2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = pow(m_c, 2);
    const ccomplex_t IT_0006 = cpow(IT_0004, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*conj(V_cd)
      *e_em*(IT_0000*IT_0001*IT_0003*IT_0004 + (-2)*IT_0002*(m_W*IT_0004 + (-0.5
      )*IT_0000*IT_0005*IT_0006))*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = pow(m_h, 2);
    const ccomplex_t IT_0016 = pow(m_ss_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_s, 2);
    const ccomplex_t IT_0021 = (0 + _Complex_I*(-2.82842712474619))*conj(V_cs)
      *e_em*IT_0008*(IT_0002*(m_W*IT_0004 + (-0.5)*IT_0000*IT_0005*IT_0006) + 
      -1./2*IT_0000*IT_0003*IT_0004*IT_0020);
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = cos(alpha);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0011 + m_s*e_em*mu_h*IT_0000*IT_0003*IT_0008*IT_0023);
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0008*(m_W*IT_0002*IT_0004*conj(U_sb_00) + (-0.5)*IT_0000*(IT_0002
      *IT_0005*IT_0006*conj(U_sb_00) + m_b*(m_b*IT_0003*IT_0004*conj(U_sb_00) +
       mu_h*conj(U_sb_10))));
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_00;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0000*IT_0008;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = cos(theta_W);
    const ccomplex_t IT_0041 = cpow(IT_0040, -2);
    const ccomplex_t IT_0042 = (0 + _Complex_I*-3)*e_em*(IT_0000*IT_0003
      *IT_0008*IT_0011*IT_0020 + -1./3*m_W*IT_0007*(IT_0002*IT_0011 + IT_0004
      *IT_0023)*IT_0041);
    const ccomplex_t IT_0043 = (-0.333333333333333)*IT_0042;
    const ccomplex_t IT_0044 = IT_0039*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0015 + (0 + _Complex_I*1)
      *m_ss_R*Gamma_sr + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0011*U_sb_01;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0008*(m_W*IT_0002*IT_0004*conj(U_sb_01) + (-0.5)*IT_0000*(IT_0002
      *IT_0005*IT_0006*conj(U_sb_01) + m_b*(m_b*IT_0003*IT_0004*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = IT_0048*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = cpow(IT_0004, 3);
    const ccomplex_t IT_0057 = cpow(IT_0002, 2);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_W*e_em*(IT_0008*IT_0011
      *IT_0056 + (-2)*IT_0002*(IT_0004*(IT_0002*IT_0011 + IT_0004*IT_0023)*(0.5
      *IT_0008 + IT_0007*IT_0041) + -1./2*IT_0008*IT_0023*IT_0057));
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0000*IT_0003*IT_0004*IT_0008;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = pow(m_W, 2);
    const ccomplex_t IT_0064 = pow(m_sc_L, 2);
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + -IT_0016 + IT_0063 + -IT_0064 
      + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0062*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0007*IT_0041;
    const ccomplex_t IT_0069 = -IT_0008 + IT_0068;
    const ccomplex_t IT_0070 = 3*IT_0008 + IT_0068;
    const ccomplex_t IT_0071 = cpow(IT_0004, 2);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_W*e_em*(IT_0023*IT_0056
      *IT_0069 + -IT_0002*(IT_0002*(IT_0002*IT_0011*IT_0069 + IT_0004*IT_0023
      *IT_0070) + -IT_0011*IT_0070*IT_0071));
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0039*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_12 + -IT_0016 + -IT_0064 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0000
      *IT_0006*IT_0008*IT_0011;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*conj
      (V_cs)*e_em*IT_0000*IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_13 + IT_0015 + IT_0064 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0005
      *IT_0006*IT_0008*IT_0023 + (-0.5)*m_W*(IT_0002*IT_0011 + IT_0004*IT_0023)*
      (IT_0008 + -1./3*IT_0007*IT_0041));
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0039*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_13 + IT_0015 + (0 + _Complex_I*1)
      *m_sc_L*Gamma_cl + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0019 + IT_0029 + IT_0037 + IT_0047 + IT_0055
       + -IT_0067 + -IT_0077 + IT_0085 + IT_0091;
    return create_ccomplex_return(3*IT_0092*conj(IT_0092));
}

