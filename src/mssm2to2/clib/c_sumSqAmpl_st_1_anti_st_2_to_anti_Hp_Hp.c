#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
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
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = pow(m_W, -2);
    const ccomplex_t IT_0010 = pow(m_b, 2);
    const ccomplex_t IT_0011 = cpow(V_tb, 2);
    const ccomplex_t IT_0012 = pow(m_d, 2);
    const ccomplex_t IT_0013 = pow(m_s, 2);
    const ccomplex_t IT_0014 = cpow(IT_0001, -2);
    const ccomplex_t IT_0015 = IT_0002*(IT_0004 + 3*IT_0006) + -IT_0008*
      (IT_0004 + 3*IT_0006*(1 + (-2)*IT_0009*(IT_0010*IT_0011 + V_td*conj(V_td)
      *IT_0012 + V_ts*conj(V_ts)*IT_0013)*IT_0014));
    const ccomplex_t IT_0016 = pow(m_t, 2);
    const ccomplex_t IT_0017 = cpow(IT_0007, -2);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0000*(IT_0015*U_st_00
      *conj(U_st_01) + 4*(IT_0004*IT_0008 + -IT_0002*(IT_0004 + (-1.5)*IT_0006
      *IT_0009*IT_0016*IT_0017))*U_st_10*conj(U_st_11));
    const ccomplex_t IT_0019 = (-0.0833333333333333)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = cpow(IT_0007, -1);
    const ccomplex_t IT_0022 = cpow(IT_0001, -1);
    const ccomplex_t IT_0023 = m_t*conj(U_sb_01)*(m_t*IT_0001*IT_0021*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0007*IT_0022*conj(U_sb_01)*U_st_00 + conj
      (U_sb_11)*(mu_h*U_st_00 + m_t*IT_0021*IT_0022*U_st_10));
    const ccomplex_t IT_0024 = cpow(IT_0005, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*conj
      (U_sb_11)*U_st_00 + (-1.4142135623731)*e_em*IT_0024*((-0.5)*IT_0020
      *IT_0023 + m_W*IT_0001*IT_0007*conj(U_sb_01)*U_st_00) + A_t*IT_0001*conj
      (U_sb_01)*U_st_10);
    const ccomplex_t IT_0026 = m_t*U_sb_01*(m_t*IT_0001*IT_0021*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0007*IT_0022*U_sb_01*conj(U_st_01) +
       U_sb_11*(mu_h*conj(U_st_01) + m_t*IT_0021*IT_0022*conj(U_st_11)));
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*U_sb_11
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0024*((-0.5)*IT_0020*IT_0026 +
       m_W*IT_0001*IT_0007*U_sb_01*conj(U_st_01)) + A_t*IT_0001*U_sb_01*conj
      (U_st_11));
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = pow(m_Hp, 2);
    const ccomplex_t IT_0030 = pow(m_st_2, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0004*IT_0005;
    const ccomplex_t IT_0035 = -IT_0024 + IT_0034;
    const ccomplex_t IT_0036 = cos(alpha);
    const ccomplex_t IT_0037 = cpow(IT_0007, 3);
    const ccomplex_t IT_0038 = 3*IT_0024 + IT_0034;
    const ccomplex_t IT_0039 = sin(alpha);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_W*e_em*(IT_0035*IT_0036
      *IT_0037 + IT_0001*(IT_0008*IT_0038*IT_0039 + -IT_0001*(IT_0007*IT_0036
      *IT_0038 + IT_0001*IT_0035*IT_0039)));
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = m_t*mu_h*IT_0020;
    const ccomplex_t IT_0043 = IT_0001*IT_0039;
    const ccomplex_t IT_0044 = IT_0007*IT_0036;
    const ccomplex_t IT_0045 = IT_0043 + IT_0044;
    const ccomplex_t IT_0046 = 2*U_st_00*((IT_0016*IT_0020*IT_0021*IT_0024
      *IT_0036 + 0.166666666666667*m_W*(IT_0004*IT_0005 + (-3)*IT_0024)*IT_0045)
      *conj(U_st_01) + 1./2*IT_0021*IT_0024*IT_0039*IT_0042*conj(U_st_11)) +
       U_st_10*(IT_0021*IT_0024*IT_0039*IT_0042*conj(U_st_01) + (
      -1.33333333333333)*((-1.5)*IT_0016*IT_0020*IT_0021*IT_0024*IT_0036 + m_W
      *IT_0004*IT_0005*IT_0045)*conj(U_st_11));
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(e_em*IT_0046 +
       1.4142135623731*A_t*IT_0036*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = IT_0041*IT_0048;
    const ccomplex_t IT_0050 = pow(m_st_1, 2);
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0030 + -IT_0050 + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = -IT_0024;
    const ccomplex_t IT_0055 = IT_0034 + IT_0054;
    const ccomplex_t IT_0056 = -IT_0001;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_W*e_em*(IT_0037*IT_0039
      *IT_0055 + IT_0056*(IT_0008*IT_0036*IT_0038 + -(IT_0007*IT_0038*IT_0039 + 
      -IT_0001*IT_0036*IT_0055)*IT_0056));
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0001*IT_0036;
    const ccomplex_t IT_0060 = -IT_0007*IT_0039 + IT_0059;
    const ccomplex_t IT_0061 = (-2)*U_st_00*((IT_0016*IT_0020*IT_0021*IT_0024
      *IT_0039 + (-0.166666666666667)*m_W*(IT_0004*IT_0005 + (-3)*IT_0024)
      *IT_0060)*conj(U_st_01) + -1./2*IT_0021*IT_0024*IT_0036*IT_0042*conj
      (U_st_11)) + U_st_10*(IT_0021*IT_0024*IT_0036*IT_0042*conj(U_st_01) + (
      -1.33333333333333)*(1.5*IT_0016*IT_0020*IT_0021*IT_0024*IT_0039 + m_W
      *IT_0004*IT_0005*IT_0060)*conj(U_st_11));
    const ccomplex_t IT_0062 = (0 + _Complex_I*-1)*(e_em*IT_0061 + (
      -1.4142135623731)*A_t*IT_0039*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0058*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + -IT_0030 + -IT_0050 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0030 + -IT_0050 + -reg_prop, -1);
    const ccomplex_t IT_0069 = cpow(IT_0003, -1);
    const ccomplex_t IT_0070 = (0 + _Complex_I*-1)*e_em*(IT_0003*IT_0024 + 
      -IT_0005*IT_0069);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0003
      *IT_0024 + -1./3*IT_0005*IT_0069)*U_st_00*conj(U_st_01) + IT_0005*IT_0069
      *U_st_10*conj(U_st_11));
    const ccomplex_t IT_0072 = 0.666666666666667*IT_0071;
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = (-0.5)*s_13 + (-0.5)*s_23 + 0.5*IT_0030 + 0.5
      *IT_0050;
    const ccomplex_t IT_0075 = s_12 + IT_0074;
    const ccomplex_t IT_0076 = IT_0073*IT_0075;
    const ccomplex_t IT_0077 = s_12 + IT_0030;
    const ccomplex_t IT_0078 = 0.5*IT_0070;
    const ccomplex_t IT_0079 = 1.33333333333333*IT_0071;
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = IT_0077*IT_0080;
    const ccomplex_t IT_0082 = IT_0072*IT_0078;
    const ccomplex_t IT_0083 = 0.5*IT_0030 + 0.5*IT_0050;
    const ccomplex_t IT_0084 = s_12 + IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = IT_0070*IT_0079;
    const ccomplex_t IT_0087 = -s_23;
    const ccomplex_t IT_0088 = IT_0077 + IT_0087;
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = 0.5*IT_0081 + -IT_0085 + (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0076 + IT_0090;
    const ccomplex_t IT_0092 = IT_0068*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*conj(V_ts)*(A_b*IT_0007
      *U_st_00 + 0.707106781186548*m_s*e_em*IT_0020*IT_0024*(mu_h*U_st_00 + m_t
      *IT_0021*IT_0022*U_st_10));
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*V_ts*(A_b*IT_0007*conj
      (U_st_01) + 0.707106781186548*m_s*e_em*IT_0020*IT_0024*(mu_h*conj(U_st_01)
       + m_t*IT_0021*IT_0022*conj(U_st_11)));
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*conj(V_td)*(A_b*IT_0007
      *U_st_00 + 0.707106781186548*m_d*e_em*IT_0020*IT_0024*(mu_h*U_st_00 + m_t
      *IT_0021*IT_0022*U_st_10));
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*V_td*(A_b*IT_0007*conj
      (U_st_01) + 0.707106781186548*m_d*e_em*IT_0020*IT_0024*(mu_h*conj(U_st_01)
       + m_t*IT_0021*IT_0022*conj(U_st_11)));
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = A_t*IT_0001*U_st_10 + (-1.4142135623731)*e_em
      *IT_0024*(m_W*IT_0001*IT_0007*U_st_00 + (-0.5)*IT_0020*(IT_0007*IT_0013
      *IT_0022*U_st_00 + m_t*(m_t*IT_0001*IT_0021*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*conj(V_ts)*IT_0106;
    const ccomplex_t IT_0108 = A_t*IT_0001*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0024*(m_W*IT_0001*IT_0007*conj(U_st_01) + (-0.5)*IT_0020*(IT_0007
      *IT_0013*IT_0022*conj(U_st_01) + m_t*(m_t*IT_0001*IT_0021*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*V_ts*IT_0108;
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = A_t*IT_0001*U_st_10 + (-1.4142135623731)*e_em
      *IT_0024*(m_W*IT_0001*IT_0007*U_st_00 + (-0.5)*IT_0020*(IT_0007*IT_0012
      *IT_0022*U_st_00 + m_t*(m_t*IT_0001*IT_0021*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*conj(V_td)*IT_0114;
    const ccomplex_t IT_0116 = A_t*IT_0001*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0024*(m_W*IT_0001*IT_0007*conj(U_st_01) + (-0.5)*IT_0020*(IT_0007
      *IT_0012*IT_0022*conj(U_st_01) + m_t*(m_t*IT_0001*IT_0021*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*V_td*IT_0116;
    const ccomplex_t IT_0118 = IT_0115*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = m_t*conj(U_sb_00)*(m_t*IT_0001*IT_0021*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0007*IT_0022*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0021*IT_0022*U_st_10));
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0024*((-0.5)*IT_0020
      *IT_0122 + m_W*IT_0001*IT_0007*conj(U_sb_00)*U_st_00) + A_t*IT_0001*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0124 = m_t*U_sb_00*(m_t*IT_0001*IT_0021*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0007*IT_0022*U_sb_00*conj(U_st_01) +
       U_sb_10*(mu_h*conj(U_st_01) + m_t*IT_0021*IT_0022*conj(U_st_11)));
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*U_sb_10
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0024*((-0.5)*IT_0020*IT_0124 +
       m_W*IT_0001*IT_0007*U_sb_00*conj(U_st_01)) + A_t*IT_0001*U_sb_00*conj
      (U_st_11));
    const ccomplex_t IT_0126 = IT_0123*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_23 + IT_0029 + IT_0030 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0126*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0019 + IT_0033 + -IT_0053 + -IT_0067 + 2
      *IT_0093 + IT_0099 + IT_0105 + IT_0113 + IT_0121 + IT_0129;
    return create_ccomplex_return(3*IT_0130*conj(IT_0130));
}

