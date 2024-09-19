#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_sb_1_to_A0_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_sb_1_to_A0_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
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
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = pow(m_W, -2);
    const ccomplex_t IT_0003 = pow(m_b, 2);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -2);
    const ccomplex_t IT_0006 = pow(m_t, 2);
    const ccomplex_t IT_0007 = cpow(IT_0001, -2);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = 1.4142135623731*V_tb*IT_0000*(IT_0001*(IT_0001 
      + -IT_0001*IT_0002*IT_0003*IT_0005) + -IT_0004*(IT_0004 + -IT_0002*IT_0004
      *IT_0006*IT_0007))*IT_0009*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0011 = (-0.25)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = cpow(IT_0001, -1);
    const ccomplex_t IT_0014 = cpow(IT_0004, -1);
    const ccomplex_t IT_0015 = m_t*conj(U_sb_00)*(m_t*IT_0004*IT_0013*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0001*IT_0014*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0013*IT_0014*U_st_10));
    const ccomplex_t IT_0016 = cpow(IT_0008, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0001*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0016*((-0.5)*IT_0012
      *IT_0015 + m_W*IT_0001*IT_0004*conj(U_sb_00)*U_st_00) + A_t*IT_0004*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0018 = m_b*e_em*mu_h*IT_0012*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = m_b*e_em*mu_h*IT_0012*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = cpow(V_tb, 2);
    const ccomplex_t IT_0023 = A_b*IT_0022*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = A_b*IT_0022*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (-0.5)*IT_0019 + 0.5*IT_0021 + (-0.5)*IT_0025 +
       0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0017*IT_0029;
    const ccomplex_t IT_0031 = pow(m_A0, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0031 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = m_t*conj(U_sb_01)*(m_t*IT_0004*IT_0013*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0001*IT_0014*conj(U_sb_01)*U_st_00 + conj
      (U_sb_11)*(mu_h*U_st_00 + m_t*IT_0013*IT_0014*U_st_10));
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0001*conj
      (U_sb_11)*U_st_00 + (-1.4142135623731)*e_em*IT_0016*((-0.5)*IT_0012
      *IT_0035 + m_W*IT_0001*IT_0004*conj(U_sb_01)*U_st_00) + A_t*IT_0004*conj
      (U_sb_01)*U_st_10);
    const ccomplex_t IT_0037 = m_b*e_em*mu_h*IT_0012*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0038 = IT_0016*IT_0037;
    const ccomplex_t IT_0039 = m_b*e_em*mu_h*IT_0012*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0040 = IT_0016*IT_0039;
    const ccomplex_t IT_0041 = A_b*IT_0022*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0042 = IT_0001*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = A_b*IT_0022*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0045 = IT_0001*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0038 + (-0.5)*IT_0040 + 0.5*IT_0043 + (
      -0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0036*IT_0047;
    const ccomplex_t IT_0049 = pow(m_sb_1, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0031 + IT_0049 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*conj(V_ts)*(A_b*IT_0001
      *U_st_00 + 0.707106781186548*m_s*e_em*IT_0012*IT_0016*(mu_h*U_st_00 + m_t
      *IT_0013*IT_0014*U_st_10));
    const ccomplex_t IT_0054 = 1.4142135623731*A_b*V_tb*V_ts*IT_0001*conj
      (U_sb_00);
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0031 + IT_0049 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = pow(m_st_1, 2);
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0049 + -IT_0060 + -reg_prop, -1);
    const ccomplex_t IT_0062 = e_em*IT_0016;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0016*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0049 + 0.5*IT_0060;
    const ccomplex_t IT_0068 = s_12 + IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = s_12 + IT_0049;
    const ccomplex_t IT_0071 = -IT_0064;
    const ccomplex_t IT_0072 = IT_0063*IT_0071;
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = IT_0062*IT_0071;
    const ccomplex_t IT_0075 = -s_23;
    const ccomplex_t IT_0076 = IT_0070 + IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = IT_0062*IT_0065;
    const ccomplex_t IT_0079 = (-0.5)*s_13 + (-0.5)*s_23 + 0.5*IT_0049 + 0.5
      *IT_0060;
    const ccomplex_t IT_0080 = s_12 + IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = (-0.5)*IT_0073 + 0.5*IT_0077 + -IT_0081;
    const ccomplex_t IT_0083 = IT_0069 + IT_0082;
    const ccomplex_t IT_0084 = IT_0061*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = m_t*e_em*mu_h*IT_0012*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0087 = IT_0016*IT_0086;
    const ccomplex_t IT_0088 = m_t*e_em*mu_h*IT_0012*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0089 = IT_0016*IT_0088;
    const ccomplex_t IT_0090 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0091 = IT_0004*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0094 = IT_0004*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (-0.5)*IT_0087 + 0.5*IT_0089 + (-0.5)*IT_0092 +
       0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0017*IT_0096;
    const ccomplex_t IT_0098 = cpow((-2)*s_13 + IT_0031 + (0 + _Complex_I*1)
      *m_st_1*Gamma_t1 + reg_prop, -1);
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = m_t*e_em*mu_h*IT_0012*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0102 = IT_0016*IT_0101;
    const ccomplex_t IT_0103 = m_t*e_em*mu_h*IT_0012*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0104 = IT_0016*IT_0103;
    const ccomplex_t IT_0105 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0106 = IT_0004*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0109 = IT_0004*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (-0.5)*IT_0102 + 0.5*IT_0104 + (-0.5)*IT_0107 +
       0.5*IT_0110;
    const ccomplex_t IT_0112 = m_t*conj(U_sb_00)*(m_t*IT_0004*IT_0013*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0001*IT_0014*conj(U_sb_00)*U_st_01 + conj
      (U_sb_10)*(mu_h*U_st_01 + m_t*IT_0013*IT_0014*U_st_11));
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0001*conj
      (U_sb_10)*U_st_01 + (-1.4142135623731)*e_em*IT_0016*((-0.5)*IT_0012
      *IT_0112 + m_W*IT_0001*IT_0004*conj(U_sb_00)*U_st_01) + A_t*IT_0004*conj
      (U_sb_00)*U_st_11);
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_13 + IT_0031 + IT_0060 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = pow(m_d, 2);
    const ccomplex_t IT_0119 = A_t*IT_0004*U_st_10 + (-1.4142135623731)*e_em
      *IT_0016*(m_W*IT_0001*IT_0004*U_st_00 + (-0.5)*IT_0012*(IT_0001*IT_0014
      *IT_0118*U_st_00 + m_t*(m_t*IT_0004*IT_0013*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*conj(V_td)*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*A_b*V_tb*V_td*IT_0001*conj
      (U_sb_10);
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_23 + IT_0031 + IT_0049 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = pow(m_s, 2);
    const ccomplex_t IT_0128 = A_t*IT_0004*U_st_10 + (-1.4142135623731)*e_em
      *IT_0016*(m_W*IT_0001*IT_0004*U_st_00 + (-0.5)*IT_0012*(IT_0001*IT_0014
      *IT_0127*U_st_00 + m_t*(m_t*IT_0004*IT_0013*U_st_00 + mu_h*U_st_10)));
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*conj(V_ts)*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*A_b*V_tb*V_ts*IT_0001*conj
      (U_sb_10);
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = cpow((-2)*s_23 + IT_0031 + IT_0049 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0134 = IT_0132*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*conj(V_td)*(A_b*IT_0001
      *U_st_00 + 0.707106781186548*m_d*e_em*IT_0012*IT_0016*(mu_h*U_st_00 + m_t
      *IT_0013*IT_0014*U_st_10));
    const ccomplex_t IT_0137 = 1.4142135623731*A_b*V_tb*V_td*IT_0001*conj
      (U_sb_00);
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = IT_0136*IT_0138;
    const ccomplex_t IT_0140 = cpow((-2)*s_23 + IT_0031 + IT_0049 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = pow(m_W, 2);
    const ccomplex_t IT_0144 = cpow((-2)*s_12 + -IT_0049 + -IT_0060 + IT_0143 
      + -reg_prop, -1);
    const ccomplex_t IT_0145 = cpow(IT_0004, 2);
    const ccomplex_t IT_0146 = cpow(IT_0001, 2);
    const ccomplex_t IT_0147 = m_b*mu_h*IT_0001*IT_0012*IT_0014*conj(U_sb_10)
      *U_st_00 + -conj(U_sb_00)*(((IT_0003 + -IT_0006)*IT_0012 + -m_W*(IT_0145 +
       -IT_0146))*U_st_00 + m_t*mu_h*IT_0004*IT_0012*IT_0013*U_st_10);
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*V_tb*((-0.707106781186548)
      *e_em*IT_0016*IT_0147 + A_b*IT_0004*conj(U_sb_10)*U_st_00 + -A_t*IT_0001
      *conj(U_sb_00)*U_st_10);
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = m_W*e_em*IT_0016;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = IT_0149*IT_0151;
    const ccomplex_t IT_0153 = IT_0144*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0011 + IT_0034 + IT_0052 + IT_0059 + (-2)
      *IT_0085 + IT_0100 + IT_0117 + IT_0126 + IT_0135 + IT_0142 + -IT_0154;
    return create_ccomplex_return(3*IT_0155*conj(IT_0155));
}

