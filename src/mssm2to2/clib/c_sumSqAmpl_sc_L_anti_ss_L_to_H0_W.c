#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_H0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_H0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0001;
    const ccomplex_t IT_0003 = (-0.5)*IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*m_W*e_em*IT_0001*(IT_0004
      *IT_0005 + IT_0006*IT_0007);
    const ccomplex_t IT_0009 = pow(m_sc_L, 2);
    const ccomplex_t IT_0010 = pow(m_ss_L, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0009 + -IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -2);
    const ccomplex_t IT_0017 = pow(s_14, 2);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0005
      *IT_0006 + -IT_0004*IT_0007);
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = pow(m_s, 2);
    const ccomplex_t IT_0022 = cpow(IT_0004, -1);
    const ccomplex_t IT_0023 = pow(m_c, 2);
    const ccomplex_t IT_0024 = cpow(IT_0006, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0001*(IT_0006*IT_0020*IT_0021*IT_0022 + (-2)*IT_0004*(m_W*IT_0006
       + (-0.5)*IT_0020*IT_0023*IT_0024));
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0019*IT_0029;
    const ccomplex_t IT_0031 = -IT_0002;
    const ccomplex_t IT_0032 = cos(theta_W);
    const ccomplex_t IT_0033 = cpow(IT_0032, -2);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0005
      *IT_0020*IT_0021*IT_0022 + (-0.5)*m_W*(IT_0004*IT_0005 + -IT_0006*IT_0007)
      *(IT_0001 + 1./3*IT_0000*IT_0033));
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = pow(m_H0, 2);
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0036 + (0 + _Complex_I*1)
      *m_ss_L*Gamma_sl + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0031*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0001*conj(U_sb_01);
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0005*U_sb_11;
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0010 + IT_0036 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0042*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0004
      *IT_0005 + IT_0006*IT_0007);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = cpow(IT_0004, 2);
    const ccomplex_t IT_0052 = cpow(IT_0006, 2);
    const ccomplex_t IT_0053 = (0 + _Complex_I*(-1.4142135623731))*conj(V_cs)
      *e_em*IT_0001*(IT_0020*(IT_0021 + -IT_0023) + -m_W*(IT_0051 + -IT_0052));
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = pow(m_W, 2);
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + IT_0055 
      + -reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0001*conj(U_sb_00);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0005*U_sb_10;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0010 + IT_0036 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = IT_0061*IT_0066;
    const ccomplex_t IT_0068 = IT_0030 + -IT_0040 + -IT_0048 + IT_0059 + 
      -IT_0067;
    const ccomplex_t IT_0069 = s_14*s_34*IT_0016;
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = s_13 + IT_0070;
    const ccomplex_t IT_0072 = -IT_0031;
    const ccomplex_t IT_0073 = IT_0003 + IT_0072;
    const ccomplex_t IT_0074 = IT_0013*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0007
      *IT_0020*IT_0023*IT_0024 + 0.5*m_W*(IT_0004*IT_0005 + -IT_0006*IT_0007)*
      (IT_0001 + -1./3*IT_0000*IT_0033));
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_13 + IT_0036 + (0 + _Complex_I*1)
      *m_sc_L*Gamma_cl + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0031*IT_0079;
    const ccomplex_t IT_0081 = IT_0040 + IT_0048 + IT_0067 + -IT_0074 + IT_0080;
    const ccomplex_t IT_0082 = s_14*s_24*IT_0016;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = s_12 + IT_0083;
    const ccomplex_t IT_0085 = (-3)*conj(IT_0015);
    const ccomplex_t IT_0086 = pow(s_34, 2);
    const ccomplex_t IT_0087 = IT_0016*IT_0086;
    const ccomplex_t IT_0088 = s_24*s_34*IT_0016;
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = s_23 + IT_0089;
    const ccomplex_t IT_0091 = pow(s_24, 2);
    const ccomplex_t IT_0092 = IT_0016*IT_0091;
    const ccomplex_t IT_0093 = (-3)*IT_0081;
    const ccomplex_t IT_0094 = (-3)*IT_0015*(conj(IT_0015)*(IT_0009 + -IT_0018
      ) + conj(IT_0068)*IT_0071 + conj(IT_0081)*IT_0084) + IT_0081*IT_0084
      *IT_0085 + (-3)*IT_0068*((-0.333333333333333)*IT_0071*IT_0085 + conj
      (IT_0068)*(IT_0036 + -IT_0087) + conj(IT_0081)*IT_0090) + (conj(IT_0068)
      *IT_0090 + conj(IT_0081)*(IT_0010 + -IT_0092))*IT_0093;
    return create_ccomplex_return(IT_0094);
}

