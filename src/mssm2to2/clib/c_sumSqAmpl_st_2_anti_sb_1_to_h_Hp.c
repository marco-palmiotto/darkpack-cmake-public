#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_sb_1_to_h_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_sb_1_to_h_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
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
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
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
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = pow(m_W, -2);
    const ccomplex_t IT_0004 = pow(m_t, 2);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = pow(m_b, 2);
    const ccomplex_t IT_0008 = cpow(IT_0001, -2);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = cpow(IT_0001, -1);
    const ccomplex_t IT_0011 = cpow(IT_0005, -1);
    const ccomplex_t IT_0012 = (IT_0001*(IT_0002 + -IT_0002*IT_0003*IT_0004
      *IT_0006) + -IT_0005*(1 + -IT_0003*IT_0007*IT_0008)*IT_0009)*conj(U_sb_00)
      *U_st_01 + -m_b*m_t*IT_0003*(IT_0002*IT_0010 + -IT_0009*IT_0011)*conj
      (U_sb_10)*U_st_11;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -2);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*V_tb*IT_0000
      *IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (-0.25)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_h, 2);
    const ccomplex_t IT_0018 = pow(m_st_2, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = m_t*conj(U_sb_00)*(m_t*IT_0001*IT_0011*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0005*IT_0010*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0010*IT_0011*U_st_10));
    const ccomplex_t IT_0022 = cpow(IT_0013, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0022*((-0.5)*IT_0020
      *IT_0021 + m_W*IT_0001*IT_0005*conj(U_sb_00)*U_st_00) + A_t*IT_0001*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0024 = m_t*mu_h*IT_0020;
    const ccomplex_t IT_0025 = IT_0001*IT_0009;
    const ccomplex_t IT_0026 = IT_0002*IT_0005;
    const ccomplex_t IT_0027 = IT_0025 + IT_0026;
    const ccomplex_t IT_0028 = cos(theta_W);
    const ccomplex_t IT_0029 = cpow(IT_0028, -2);
    const ccomplex_t IT_0030 = 2*conj(U_st_00)*((IT_0002*IT_0004*IT_0011
      *IT_0020*IT_0022 + (-0.5)*m_W*IT_0027*(IT_0022 + -1./3*IT_0013*IT_0029))
      *U_st_01 + 1./2*IT_0009*IT_0011*IT_0022*IT_0024*U_st_11) + conj(U_st_10)*
      (IT_0009*IT_0011*IT_0022*IT_0024*U_st_01 + (-1.33333333333333)*((-1.5)
      *IT_0002*IT_0004*IT_0011*IT_0020*IT_0022 + m_W*IT_0013*IT_0027*IT_0029)
      *U_st_11);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(e_em*IT_0030 +
       1.4142135623731*A_t*IT_0002*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0023*IT_0032;
    const ccomplex_t IT_0034 = IT_0019*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*conj(V_td)*(A_b*IT_0005
      *U_st_01 + 0.707106781186548*m_d*e_em*IT_0020*IT_0022*(mu_h*U_st_01 + m_t
      *IT_0010*IT_0011*U_st_11));
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0009*conj(U_sb_00);
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = pow(m_sb_1, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0017 + IT_0040 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*conj(V_ts)*(A_b*IT_0005
      *U_st_01 + 0.707106781186548*m_s*e_em*IT_0020*IT_0022*(mu_h*U_st_01 + m_t
      *IT_0010*IT_0011*U_st_11));
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0009*conj(U_sb_00);
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0017 + IT_0040 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0017 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0052 = m_t*conj(U_sb_00)*(m_t*IT_0001*IT_0011*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0005*IT_0010*conj(U_sb_00)*U_st_01 + conj
      (U_sb_10)*(mu_h*U_st_01 + m_t*IT_0010*IT_0011*U_st_11));
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_10)*U_st_01 + (-1.4142135623731)*e_em*IT_0022*((-0.5)*IT_0020
      *IT_0052 + m_W*IT_0001*IT_0005*conj(U_sb_00)*U_st_01) + A_t*IT_0001*conj
      (U_sb_00)*U_st_11);
    const ccomplex_t IT_0054 = cpow(V_tb, 2);
    const ccomplex_t IT_0055 = m_b*mu_h*IT_0020;
    const ccomplex_t IT_0056 = U_sb_10*(IT_0002*IT_0010*IT_0022*IT_0055*conj
      (U_sb_00) + 2*(IT_0007*IT_0009*IT_0010*IT_0020*IT_0022 + -1./3*m_W*IT_0013
      *IT_0027*IT_0029)*conj(U_sb_10)) + 2*U_sb_00*((IT_0007*IT_0009*IT_0010
      *IT_0020*IT_0022 + (-0.5)*m_W*IT_0027*(IT_0022 + 1./3*IT_0013*IT_0029))
      *conj(U_sb_00) + 1./2*IT_0002*IT_0010*IT_0022*IT_0055*conj(U_sb_10));
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(e_em*IT_0056 +
       1.4142135623731*A_b*IT_0009*IT_0054*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0053*IT_0058;
    const ccomplex_t IT_0060 = IT_0051*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_12 + -IT_0018 + -IT_0040 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0013*IT_0029;
    const ccomplex_t IT_0064 = -IT_0022 + IT_0063;
    const ccomplex_t IT_0065 = cpow(IT_0005, 3);
    const ccomplex_t IT_0066 = 3*IT_0022 + IT_0063;
    const ccomplex_t IT_0067 = cpow(IT_0005, 2);
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*m_W*e_em*(IT_0002*IT_0064
      *IT_0065 + -IT_0001*(IT_0001*(IT_0001*IT_0009*IT_0064 + IT_0002*IT_0005
      *IT_0066) + -IT_0009*IT_0066*IT_0067));
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = IT_0053*IT_0069;
    const ccomplex_t IT_0071 = IT_0062*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = pow(m_s, 2);
    const ccomplex_t IT_0074 = A_t*IT_0001*U_st_11 + (-1.4142135623731)*e_em
      *IT_0022*(m_W*IT_0001*IT_0005*U_st_01 + (-0.5)*IT_0020*(IT_0005*IT_0010
      *IT_0073*U_st_01 + m_t*(m_t*IT_0001*IT_0011*U_st_01 + mu_h*U_st_11)));
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*conj(V_ts)*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0009*conj(U_sb_10);
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0075*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_23 + IT_0017 + IT_0040 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0018 + -IT_0040 + -reg_prop, -1);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0022*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 +
       IT_0005*IT_0009)*IT_0022;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*s_13 + (-0.5)*s_23 + 0.5*IT_0018 + 0.5
      *IT_0040;
    const ccomplex_t IT_0089 = s_12 + IT_0088;
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = s_12 + IT_0040;
    const ccomplex_t IT_0092 = -IT_0083;
    const ccomplex_t IT_0093 = (-0.5)*IT_0085;
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = IT_0091*IT_0094;
    const ccomplex_t IT_0096 = IT_0086*IT_0092;
    const ccomplex_t IT_0097 = -s_23;
    const ccomplex_t IT_0098 = IT_0091 + IT_0097;
    const ccomplex_t IT_0099 = IT_0096*IT_0098;
    const ccomplex_t IT_0100 = -IT_0090 + (-0.5)*IT_0095 + 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0084*IT_0093;
    const ccomplex_t IT_0102 = 0.5*IT_0018 + 0.5*IT_0040;
    const ccomplex_t IT_0103 = s_12 + IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = IT_0100 + IT_0104;
    const ccomplex_t IT_0106 = IT_0082*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_13 + IT_0017 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0109 = 2*U_st_01*((IT_0002*IT_0004*IT_0011*IT_0020
      *IT_0022 + (-0.5)*m_W*IT_0027*(IT_0022 + -1./3*IT_0013*IT_0029))*conj
      (U_st_01) + 1./2*IT_0009*IT_0011*IT_0022*IT_0024*conj(U_st_11)) + U_st_11*
      (IT_0009*IT_0011*IT_0022*IT_0024*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0002*IT_0004*IT_0011*IT_0020*IT_0022 + m_W*IT_0013*IT_0027
      *IT_0029)*conj(U_st_11));
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(e_em*IT_0109 +
       1.4142135623731*A_t*IT_0002*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = IT_0053*IT_0111;
    const ccomplex_t IT_0113 = IT_0108*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = pow(m_d, 2);
    const ccomplex_t IT_0116 = A_t*IT_0001*U_st_11 + (-1.4142135623731)*e_em
      *IT_0022*(m_W*IT_0001*IT_0005*U_st_01 + (-0.5)*IT_0020*(IT_0005*IT_0010
      *IT_0115*U_st_01 + m_t*(m_t*IT_0001*IT_0011*U_st_01 + mu_h*U_st_11)));
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*conj(V_td)*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0009*conj(U_sb_10);
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0117*IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_23 + IT_0017 + IT_0040 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_23 + IT_0017 + IT_0040 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0125 = m_t*conj(U_sb_01)*(m_t*IT_0001*IT_0011*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0005*IT_0010*conj(U_sb_01)*U_st_01 + conj
      (U_sb_11)*(mu_h*U_st_01 + m_t*IT_0010*IT_0011*U_st_11));
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_11)*U_st_01 + (-1.4142135623731)*e_em*IT_0022*((-0.5)*IT_0020
      *IT_0125 + m_W*IT_0001*IT_0005*conj(U_sb_01)*U_st_01) + A_t*IT_0001*conj
      (U_sb_01)*U_st_11);
    const ccomplex_t IT_0127 = conj(U_sb_10)*(IT_0002*IT_0010*IT_0022*IT_0055
      *U_sb_01 + 2*(IT_0007*IT_0009*IT_0010*IT_0020*IT_0022 + -1./3*m_W*IT_0013
      *IT_0027*IT_0029)*U_sb_11) + 2*conj(U_sb_00)*((IT_0007*IT_0009*IT_0010
      *IT_0020*IT_0022 + (-0.5)*m_W*IT_0027*(IT_0022 + 1./3*IT_0013*IT_0029))
      *U_sb_01 + 1./2*IT_0002*IT_0010*IT_0022*IT_0055*U_sb_11);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(e_em*IT_0127 +
       1.4142135623731*A_b*IT_0009*IT_0054*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = IT_0126*IT_0129;
    const ccomplex_t IT_0131 = IT_0124*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = cpow(IT_0001, 2);
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*m_W*e_em*(IT_0009*IT_0022
      *IT_0065 + (-2)*IT_0001*(IT_0005*(IT_0002*IT_0005 + IT_0001*IT_0009)*(0.5
      *IT_0022 + IT_0013*IT_0029) + -1./2*IT_0002*IT_0022*IT_0133));
    const ccomplex_t IT_0135 = 0.5*IT_0134;
    const ccomplex_t IT_0136 = m_b*mu_h*IT_0005*IT_0010*IT_0020*conj(U_sb_10)
      *U_st_01 + conj(U_sb_00)*(((IT_0004 + -IT_0007)*IT_0020 + -m_W*(IT_0067 + 
      -IT_0133))*U_st_01 + -m_t*mu_h*IT_0001*IT_0011*IT_0020*U_st_11);
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*V_tb*((-0.707106781186548)
      *e_em*IT_0022*IT_0136 + A_b*IT_0001*conj(U_sb_10)*U_st_01 + -A_t*IT_0005
      *conj(U_sb_00)*U_st_11);
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0135*IT_0138;
    const ccomplex_t IT_0140 = pow(m_W, 2);
    const ccomplex_t IT_0141 = cpow((-2)*s_12 + -IT_0018 + -IT_0040 + IT_0140 
      + -reg_prop, -1);
    const ccomplex_t IT_0142 = IT_0139*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0016 + IT_0035 + IT_0043 + IT_0050 + IT_0061
       + -IT_0072 + IT_0081 + (-2)*IT_0107 + IT_0114 + IT_0123 + IT_0132 + 
      -IT_0143;
    return create_ccomplex_return(3*IT_0144*conj(IT_0144));
}

