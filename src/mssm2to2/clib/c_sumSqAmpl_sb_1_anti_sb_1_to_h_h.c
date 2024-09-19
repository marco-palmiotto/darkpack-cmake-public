#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_h_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
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
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_h, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0002 = cpow(V_tb, 2);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_b*mu_h*IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0003*IT_0006;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = IT_0008*IT_0012;
    const ccomplex_t IT_0014 = IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_b, 2);
    const ccomplex_t IT_0018 = 2*U_sb_00*(((-0.5)*m_W*IT_0014*(IT_0010 + 1./3
      *IT_0009*IT_0016) + IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*conj(U_sb_00)
       + 1./2*IT_0005*IT_0007*IT_0008*IT_0010*conj(U_sb_10)) + U_sb_10*(IT_0005
      *IT_0007*IT_0008*IT_0010*conj(U_sb_00) + (-0.666666666666667)*(m_W*IT_0009
      *IT_0014*IT_0016 + (-3)*IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*conj
      (U_sb_10));
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(e_em*IT_0018 +
       1.4142135623731*A_b*IT_0002*IT_0003*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = cpow(IT_0020, 2);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0003*U_sb_10;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = pow(m_sb_1, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0000 + IT_0031 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0000 + IT_0031 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0003*U_sb_00;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0000 + IT_0031 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_12 + (-2)*IT_0031 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0009*IT_0016;
    const ccomplex_t IT_0048 = IT_0010 + IT_0047;
    const ccomplex_t IT_0049 = cpow(IT_0003, 3);
    const ccomplex_t IT_0050 = cpow(IT_0008, 2);
    const ccomplex_t IT_0051 = cpow(IT_0003, 2);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_W*e_em*IT_0048*(IT_0012
      *IT_0049 + (-5)*IT_0008*(IT_0003*IT_0008*IT_0012 + (-0.2)*IT_0006*(IT_0050
       + (-5)*IT_0051)));
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = -IT_0003*IT_0012;
    const ccomplex_t IT_0055 = IT_0006*IT_0008;
    const ccomplex_t IT_0056 = IT_0054 + IT_0055;
    const ccomplex_t IT_0057 = (-2)*U_sb_00*((IT_0004*IT_0007*IT_0008*IT_0010
      *IT_0017 + (-0.5)*m_W*(IT_0010 + 1./3*IT_0009*IT_0016)*IT_0056)*conj
      (U_sb_00) + -1./2*IT_0003*IT_0005*IT_0007*IT_0010*conj(U_sb_10)) + U_sb_10
      *(IT_0003*IT_0005*IT_0007*IT_0010*conj(U_sb_00) + (-2)*(IT_0004*IT_0007
      *IT_0008*IT_0010*IT_0017 + -1./3*m_W*IT_0009*IT_0016*IT_0056)*conj(U_sb_10
      ));
    const ccomplex_t IT_0058 = (0 + _Complex_I*-1)*(e_em*IT_0057 + (
      -1.4142135623731)*A_b*IT_0002*IT_0008*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = IT_0053*IT_0059;
    const ccomplex_t IT_0061 = IT_0046*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0031 + -reg_prop, -1);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*m_W*e_em*IT_0048*(IT_0006
      *IT_0049 + -IT_0008*(IT_0008*(IT_0003*IT_0006 + IT_0008*IT_0012) + 
      -IT_0012*IT_0051));
    const ccomplex_t IT_0065 = 1.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0020*IT_0065;
    const ccomplex_t IT_0067 = IT_0063*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_13 + IT_0000 + IT_0031 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0003*U_sb_10;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0071*IT_0073;
    const ccomplex_t IT_0075 = IT_0069*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_23 + IT_0000 + IT_0031 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0074*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0003*U_sb_00;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_13 + IT_0000 + IT_0031 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_23 + IT_0000 + IT_0031 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0084*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_23 + IT_0000 + IT_0031 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0042*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_13 + IT_0000 + IT_0031 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0095 = 2*U_sb_00*(((-0.5)*m_W*IT_0014*(IT_0010 + 1./3
      *IT_0009*IT_0016) + IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*conj(U_sb_01)
       + 1./2*IT_0005*IT_0007*IT_0008*IT_0010*conj(U_sb_11)) + U_sb_10*(IT_0005
      *IT_0007*IT_0008*IT_0010*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0009
      *IT_0014*IT_0016 + (-3)*IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*conj
      (U_sb_11));
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(e_em*IT_0095 +
       1.4142135623731*A_b*IT_0002*IT_0003*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = 2*conj(U_sb_00)*(((-0.5)*m_W*IT_0014*(IT_0010 +
       1./3*IT_0009*IT_0016) + IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*U_sb_01 
      + 1./2*IT_0005*IT_0007*IT_0008*IT_0010*U_sb_11) + conj(U_sb_10)*(IT_0005
      *IT_0007*IT_0008*IT_0010*U_sb_01 + (-0.666666666666667)*(m_W*IT_0009
      *IT_0014*IT_0016 + (-3)*IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*U_sb_11);
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(e_em*IT_0098 +
       1.4142135623731*A_b*IT_0002*IT_0003*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0097*IT_0100;
    const ccomplex_t IT_0102 = IT_0094*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = cpow((-2)*s_23 + IT_0000 + IT_0031 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0101*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = pow(e_em, 2);
    const ccomplex_t IT_0108 = -IT_0051;
    const ccomplex_t IT_0109 = pow(m_W, -2);
    const ccomplex_t IT_0110 = cpow(IT_0006, -2);
    const ccomplex_t IT_0111 = cpow(IT_0009, -2);
    const ccomplex_t IT_0112 = 0.5*(IT_0050*(IT_0016 + 3*IT_0111) + IT_0108*
      (IT_0016 + 3*(1 + (-2)*IT_0017*IT_0109*IT_0110)*IT_0111))*U_sb_00*conj
      (U_sb_00) + (IT_0016*(IT_0050 + IT_0108) + 3*IT_0017*IT_0051*IT_0109
      *IT_0110*IT_0111)*U_sb_10*conj(U_sb_10);
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0107*IT_0112;
    const ccomplex_t IT_0114 = (-0.166666666666667)*IT_0113;
    const ccomplex_t IT_0115 = IT_0023 + IT_0025 + IT_0034 + IT_0037 + IT_0045
       + -IT_0062 + -IT_0068 + IT_0076 + IT_0079 + IT_0087 + IT_0090 + IT_0093 +
       IT_0103 + IT_0106 + IT_0114;
    return create_ccomplex_return(3*IT_0115*conj(IT_0115));
}

