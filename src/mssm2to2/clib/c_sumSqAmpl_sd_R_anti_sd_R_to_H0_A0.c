#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_anti_sd_R_to_H0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_anti_sd_R_to_H0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + -IT_0004
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0006;
    const ccomplex_t IT_0010 = (-0.666666666666667)*IT_0009;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = s_23*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0008;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = pow(m_sd_R, 2);
    const ccomplex_t IT_0016 = s_12 + IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (-0.333333333333333)*IT_0009;
    const ccomplex_t IT_0019 = IT_0008*IT_0018;
    const ccomplex_t IT_0020 = s_13 + s_23;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0018;
    const ccomplex_t IT_0023 = IT_0016*IT_0022;
    const ccomplex_t IT_0024 = -IT_0017 + -IT_0021 + 2*IT_0023;
    const ccomplex_t IT_0025 = IT_0012 + IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = m_d*e_em*mu_h*IT_0029;
    const ccomplex_t IT_0031 = IT_0007*IT_0030;
    const ccomplex_t IT_0032 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0031 + 0.5*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0004, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0002 + -m_d*e_em*mu_h*IT_0005*IT_0007*IT_0029*IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0035*IT_0038;
    const ccomplex_t IT_0040 = pow(m_H0, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0015 + IT_0040 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*A_b*conj(V_td)*V_ts*IT_0003;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0002;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0015 + IT_0040 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0031 + (-0.5)*IT_0034;
    const ccomplex_t IT_0053 = IT_0038*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_23 + IT_0015 + IT_0040 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*A_b*V_td*conj(V_ts)*IT_0003;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0002;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0058*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0015 + IT_0040 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_23 + IT_0015 + IT_0040 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0066 = 1.4142135623731*A_b*V_tb*V_td*IT_0003*conj
      (U_sb_00);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0002*U_sb_00;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = IT_0065*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0002*conj(U_sb_00);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0003
      *U_sb_00;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_13 + IT_0015 + IT_0040 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_13 + IT_0015 + IT_0040 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0002*conj(U_sb_01);
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0003
      *U_sb_01;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = IT_0081*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0015 + IT_0040 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0002*U_sb_01;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*A_b*V_tb*V_td*IT_0003*conj
      (U_sb_01);
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = IT_0091*IT_0093;
    const ccomplex_t IT_0095 = IT_0089*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = -IT_0028 + IT_0043 + IT_0051 + IT_0056 +
       IT_0064 + IT_0072 + IT_0080 + IT_0088 + IT_0096;
    return create_ccomplex_return(3*IT_0097*conj(IT_0097));
}

