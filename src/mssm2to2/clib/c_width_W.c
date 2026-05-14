#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_cb = param->V_cb;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = -IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = pow(m_W, -4);
    const ccomplex_t IT_0009 = pow(m_Z, 2);
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, 2);
    const ccomplex_t IT_0012 = IT_0010 + IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0012, 2);
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = pow(m_Z, -2);
    const ccomplex_t IT_0016 = cpow((-0.5)*IT_0009 + IT_0011, 2);
    const ccomplex_t IT_0017 = -IT_0011 + 0.5*IT_0001*IT_0009*((-0.5)*IT_0009 
      + IT_0011) + IT_0001*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0017, 2);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (1 + 0.5*IT_0014)*(IT_0011 + -IT_0019);
    const ccomplex_t IT_0021 = 2*IT_0005;
    const ccomplex_t IT_0022 = IT_0001*IT_0013;
    const ccomplex_t IT_0023 = IT_0001*IT_0015*IT_0018;
    const ccomplex_t IT_0024 = (IT_0011 + -IT_0022)*(1 + 0.5*IT_0023);
    const ccomplex_t IT_0025 = IT_0001*IT_0015;
    const ccomplex_t IT_0026 = IT_0017*IT_0025;
    const ccomplex_t IT_0027 = IT_0012*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0009;
    const ccomplex_t IT_0029 = IT_0001*IT_0017;
    const ccomplex_t IT_0030 = IT_0012*IT_0029;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = IT_0028 + IT_0031;
    const ccomplex_t IT_0033 = IT_0027*IT_0032;
    const ccomplex_t IT_0034 = IT_0005 + -IT_0006;
    const ccomplex_t IT_0035 = cpow(IT_0028, 2);
    const ccomplex_t IT_0036 = IT_0001*IT_0015*IT_0035;
    const ccomplex_t IT_0037 = (IT_0011 + -IT_0022)*(1 + 0.5*IT_0036);
    const ccomplex_t IT_0038 = (-2)*IT_0005;
    const ccomplex_t IT_0039 = -IT_0006 + IT_0038;
    const ccomplex_t IT_0040 = IT_0025*IT_0028;
    const ccomplex_t IT_0041 = IT_0012*IT_0040;
    const ccomplex_t IT_0042 = IT_0001*IT_0028;
    const ccomplex_t IT_0043 = IT_0012*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0017 + IT_0044;
    const ccomplex_t IT_0046 = IT_0041*IT_0045;
    const ccomplex_t IT_0047 = -IT_0014 + -IT_0023;
    const ccomplex_t IT_0048 = IT_0008*IT_0015;
    const ccomplex_t IT_0049 = IT_0017*IT_0048;
    const ccomplex_t IT_0050 = IT_0028*IT_0049;
    const ccomplex_t IT_0051 = IT_0012*IT_0050;
    const ccomplex_t IT_0052 = IT_0047 + IT_0051;
    const ccomplex_t IT_0053 = 1 + IT_0052;
    const ccomplex_t IT_0054 = IT_0012*IT_0053;
    const ccomplex_t IT_0055 = cpow(IT_0012, 3);
    const ccomplex_t IT_0056 = IT_0008*IT_0055;
    const ccomplex_t IT_0057 = IT_0012*IT_0036;
    const ccomplex_t IT_0058 = IT_0012*IT_0023;
    const ccomplex_t IT_0059 = IT_0015*IT_0017;
    const ccomplex_t IT_0060 = IT_0028*IT_0059;
    const ccomplex_t IT_0061 = IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0062 = IT_0017*IT_0061;
    const ccomplex_t IT_0063 = IT_0028*IT_0062;
    const ccomplex_t IT_0064 = (-0.666666666666667)*IT_0021;
    const ccomplex_t IT_0065 = IT_0015*IT_0035;
    const ccomplex_t IT_0066 = (1 + 0.5*IT_0014)*(IT_0011 + -IT_0065);
    const ccomplex_t IT_0067 = (1 + 0.5*IT_0014)*(IT_0012 + -IT_0060);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = -IT_0014 + -IT_0036;
    const ccomplex_t IT_0071 = IT_0051 + IT_0070;
    const ccomplex_t IT_0072 = 1 + IT_0071;
    const ccomplex_t IT_0073 = IT_0012*IT_0072;
    const ccomplex_t IT_0074 = 0.666666666666667*IT_0007*IT_0020 +
       0.666666666666667*IT_0021*(IT_0021*IT_0024 + IT_0006*IT_0033) +
       0.666666666666667*IT_0034*(IT_0034*IT_0037 + IT_0039*IT_0046 + -IT_0006
      *IT_0054 + 1.5*(IT_0012 + -IT_0056 + -IT_0057 + -IT_0058 + IT_0060 +
       IT_0063)*IT_0064) + 0.666666666666667*IT_0039*(IT_0039*IT_0066 + 2
      *IT_0006*IT_0069 + 1.5*IT_0064*IT_0073);
    const ccomplex_t IT_0075 = cabs(-IT_0011 + IT_0001*IT_0016);
    const ccomplex_t IT_0076 = cpow(IT_0075, 0.5);
    const ccomplex_t IT_0077 = ( ((m_W
) >= (m_W + m_Z
)) ? (0.125*IT_0000
      *IT_0001*IT_0074*IT_0076
) : (0
));
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0004;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = cpow(IT_0079, 2);
    const ccomplex_t IT_0081 = pow(m_c, 2);
    const ccomplex_t IT_0082 = pow(m_b, 2);
    const ccomplex_t IT_0083 = cpow(IT_0011 + IT_0081 + -IT_0082, 2);
    const ccomplex_t IT_0084 = -IT_0081 + 0.25*IT_0001*(IT_0011 + -IT_0081 +
       IT_0082)*(IT_0011 + IT_0081 + -IT_0082) + 0.25*IT_0001*IT_0083;
    const ccomplex_t IT_0085 = 0.5*IT_0011 + 0.5*IT_0081 + (-0.5)*IT_0082;
    const ccomplex_t IT_0086 = 0.5*IT_0011 + (-0.5)*IT_0081 + 0.5*IT_0082;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = IT_0011*IT_0084;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0087 + IT_0089;
    const ccomplex_t IT_0091 = IT_0001*IT_0090;
    const ccomplex_t IT_0092 = cabs(-IT_0081 + 0.25*IT_0001*IT_0083);
    const ccomplex_t IT_0093 = cpow(IT_0092, 0.5);
    const ccomplex_t IT_0094 = ( ((m_W
) >= (m_b + m_c
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0080*(IT_0084 + IT_0091)*IT_0093
) : (0
));
    const ccomplex_t IT_0095 = (0 + _Complex_I*1.4142135623731)*V_cd*e_em
      *IT_0004;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = pow(m_d, 2);
    const ccomplex_t IT_0098 = cpow(IT_0011 + -IT_0081 + IT_0097, 2);
    const ccomplex_t IT_0099 = -IT_0097 + 0.25*IT_0001*(IT_0011 + -IT_0081 +
       IT_0097)*(IT_0011 + IT_0081 + -IT_0097) + 0.25*IT_0001*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0011 + (-0.5)*IT_0081 + 0.5*IT_0097;
    const ccomplex_t IT_0101 = 0.5*IT_0011 + 0.5*IT_0081 + (-0.5)*IT_0097;
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = IT_0011*IT_0099;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0102 + IT_0104;
    const ccomplex_t IT_0106 = IT_0001*IT_0105;
    const ccomplex_t IT_0107 = cabs(-IT_0097 + 0.25*IT_0001*IT_0098);
    const ccomplex_t IT_0108 = cpow(IT_0107, 0.5);
    const ccomplex_t IT_0109 = ( ((m_W
) >= (m_c + m_d
)) ? (0.5*IT_0000
      *IT_0001*IT_0096*conj(IT_0096)*(IT_0099 + IT_0106)*IT_0108
) : (0
));
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004;
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = cpow(IT_0111, 2);
    const ccomplex_t IT_0113 = pow(m_e, 2);
    const ccomplex_t IT_0114 = cpow(IT_0011 + -IT_0113, 2);
    const ccomplex_t IT_0115 = 0.25*IT_0001*(IT_0011 + IT_0113)*(IT_0011 + 
      -IT_0113) + 0.25*IT_0001*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0011 + (-0.5)*IT_0113;
    const ccomplex_t IT_0117 = 0.5*IT_0011 + 0.5*IT_0113;
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = IT_0011*IT_0115;
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0118 + IT_0120;
    const ccomplex_t IT_0122 = IT_0001*IT_0121;
    const ccomplex_t IT_0123 = cabs(0.25*IT_0001*IT_0114);
    const ccomplex_t IT_0124 = cpow(IT_0123, 0.5);
    const ccomplex_t IT_0125 = ( ((m_W
) >= (m_e
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0112*(IT_0115 + IT_0122)*IT_0124
) : (0
));
    const ccomplex_t IT_0126 = cos(alpha);
    const ccomplex_t IT_0127 = sin(beta);
    const ccomplex_t IT_0128 = cos(beta);
    const ccomplex_t IT_0129 = sin(alpha);
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*m_W*e_em*IT_0004*(IT_0126
      *IT_0127 + -IT_0128*IT_0129);
    const ccomplex_t IT_0131 = cpow(IT_0130, 2);
    const ccomplex_t IT_0132 = pow(m_h, 2);
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = IT_0011 + IT_0133;
    const ccomplex_t IT_0135 = cpow(IT_0134, 2);
    const ccomplex_t IT_0136 = IT_0008*IT_0135;
    const ccomplex_t IT_0137 = cpow(IT_0011 + (-0.5)*IT_0132, 2);
    const ccomplex_t IT_0138 = cabs(-IT_0011 + IT_0001*IT_0137);
    const ccomplex_t IT_0139 = cpow(IT_0138, 0.5);
    const ccomplex_t IT_0140 = ( ((m_W
) >= (m_W + m_h
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0131*(1 + 0.5*IT_0136)*IT_0139
) :
       (0
));
    const ccomplex_t IT_0141 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0004;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = pow(m_s, 2);
    const ccomplex_t IT_0144 = cpow(IT_0011 + -IT_0081 + IT_0143, 2);
    const ccomplex_t IT_0145 = -IT_0143 + 0.25*IT_0001*(IT_0011 + -IT_0081 +
       IT_0143)*(IT_0011 + IT_0081 + -IT_0143) + 0.25*IT_0001*IT_0144;
    const ccomplex_t IT_0146 = 0.5*IT_0011 + (-0.5)*IT_0081 + 0.5*IT_0143;
    const ccomplex_t IT_0147 = 0.5*IT_0011 + 0.5*IT_0081 + (-0.5)*IT_0143;
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = IT_0011*IT_0145;
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = IT_0148 + IT_0150;
    const ccomplex_t IT_0152 = IT_0001*IT_0151;
    const ccomplex_t IT_0153 = cabs(-IT_0143 + 0.25*IT_0001*IT_0144);
    const ccomplex_t IT_0154 = cpow(IT_0153, 0.5);
    const ccomplex_t IT_0155 = ( ((m_W
) >= (m_c + m_s
)) ? (0.5*IT_0000
      *IT_0001*IT_0142*conj(IT_0142)*(IT_0145 + IT_0152)*IT_0154
) : (0
));
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004;
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = pow(m_t, 2);
    const ccomplex_t IT_0159 = cpow(IT_0011 + -IT_0082 + IT_0158, 2);
    const ccomplex_t IT_0160 = -IT_0158 + 0.25*IT_0001*(IT_0011 + -IT_0082 +
       IT_0158)*(IT_0011 + IT_0082 + -IT_0158) + 0.25*IT_0001*IT_0159;
    const ccomplex_t IT_0161 = 0.5*IT_0011 + (-0.5)*IT_0082 + 0.5*IT_0158;
    const ccomplex_t IT_0162 = 0.5*IT_0011 + 0.5*IT_0082 + (-0.5)*IT_0158;
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = IT_0011*IT_0160;
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = IT_0163 + IT_0165;
    const ccomplex_t IT_0167 = IT_0001*IT_0166;
    const ccomplex_t IT_0168 = cabs(-IT_0158 + 0.25*IT_0001*IT_0159);
    const ccomplex_t IT_0169 = cpow(IT_0168, 0.5);
    const ccomplex_t IT_0170 = ( ((m_W
) >= (m_b + m_t
)) ? (0.5*IT_0000
      *IT_0001*IT_0157*conj(IT_0157)*(IT_0160 + IT_0167)*IT_0169
) : (0
));
    const ccomplex_t IT_0171 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0004;
    const ccomplex_t IT_0172 = 0.5*IT_0171;
    const ccomplex_t IT_0173 = cpow(IT_0011 + IT_0097 + -IT_0158, 2);
    const ccomplex_t IT_0174 = -IT_0097 + 0.25*IT_0001*(IT_0011 + -IT_0097 +
       IT_0158)*(IT_0011 + IT_0097 + -IT_0158) + 0.25*IT_0001*IT_0173;
    const ccomplex_t IT_0175 = 0.5*IT_0011 + 0.5*IT_0097 + (-0.5)*IT_0158;
    const ccomplex_t IT_0176 = 0.5*IT_0011 + (-0.5)*IT_0097 + 0.5*IT_0158;
    const ccomplex_t IT_0177 = IT_0175*IT_0176;
    const ccomplex_t IT_0178 = IT_0011*IT_0174;
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = IT_0177 + IT_0179;
    const ccomplex_t IT_0181 = IT_0001*IT_0180;
    const ccomplex_t IT_0182 = cabs(-IT_0097 + 0.25*IT_0001*IT_0173);
    const ccomplex_t IT_0183 = cpow(IT_0182, 0.5);
    const ccomplex_t IT_0184 = ( ((m_W
) >= (m_d + m_t
)) ? (0.5*IT_0000
      *IT_0001*IT_0172*conj(IT_0172)*(IT_0174 + IT_0181)*IT_0183
) : (0
));
    const ccomplex_t IT_0185 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0004;
    const ccomplex_t IT_0186 = 0.5*IT_0185;
    const ccomplex_t IT_0187 = cpow(IT_0011 + IT_0143 + -IT_0158, 2);
    const ccomplex_t IT_0188 = -IT_0143 + 0.25*IT_0001*(IT_0011 + -IT_0143 +
       IT_0158)*(IT_0011 + IT_0143 + -IT_0158) + 0.25*IT_0001*IT_0187;
    const ccomplex_t IT_0189 = 0.5*IT_0011 + 0.5*IT_0143 + (-0.5)*IT_0158;
    const ccomplex_t IT_0190 = 0.5*IT_0011 + (-0.5)*IT_0143 + 0.5*IT_0158;
    const ccomplex_t IT_0191 = IT_0189*IT_0190;
    const ccomplex_t IT_0192 = IT_0011*IT_0188;
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0001*IT_0194;
    const ccomplex_t IT_0196 = cabs(-IT_0143 + 0.25*IT_0001*IT_0187);
    const ccomplex_t IT_0197 = cpow(IT_0196, 0.5);
    const ccomplex_t IT_0198 = ( ((m_W
) >= (m_s + m_t
)) ? (0.5*IT_0000
      *IT_0001*IT_0186*conj(IT_0186)*(IT_0188 + IT_0195)*IT_0197
) : (0
));
    const ccomplex_t IT_0199 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0200 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *IT_0199*V_ub_mod;
    const ccomplex_t IT_0201 = 0.5*IT_0200;
    const ccomplex_t IT_0202 = cpow(IT_0201, 2);
    const ccomplex_t IT_0203 = pow(m_u, 2);
    const ccomplex_t IT_0204 = cpow(IT_0011 + -IT_0082 + IT_0203, 2);
    const ccomplex_t IT_0205 = -IT_0203 + 0.25*IT_0001*(IT_0011 + -IT_0082 +
       IT_0203)*(IT_0011 + IT_0082 + -IT_0203) + 0.25*IT_0001*IT_0204;
    const ccomplex_t IT_0206 = 0.5*IT_0011 + (-0.5)*IT_0082 + 0.5*IT_0203;
    const ccomplex_t IT_0207 = 0.5*IT_0011 + 0.5*IT_0082 + (-0.5)*IT_0203;
    const ccomplex_t IT_0208 = IT_0206*IT_0207;
    const ccomplex_t IT_0209 = IT_0011*IT_0205;
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0208 + IT_0210;
    const ccomplex_t IT_0212 = IT_0001*IT_0211;
    const ccomplex_t IT_0213 = cabs(-IT_0203 + 0.25*IT_0001*IT_0204);
    const ccomplex_t IT_0214 = cpow(IT_0213, 0.5);
    const ccomplex_t IT_0215 = ( ((m_W
) >= (m_b + m_u
)) ? (0.5*IT_0000
      *IT_0001*IT_0202*(IT_0205 + IT_0212)*IT_0214
) : (0
));
    const ccomplex_t IT_0216 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0004;
    const ccomplex_t IT_0217 = 0.5*IT_0216;
    const ccomplex_t IT_0218 = cpow(IT_0217, 2);
    const ccomplex_t IT_0219 = cpow(IT_0011 + -IT_0097 + IT_0203, 2);
    const ccomplex_t IT_0220 = -IT_0203 + 0.25*IT_0001*(IT_0011 + -IT_0097 +
       IT_0203)*(IT_0011 + IT_0097 + -IT_0203) + 0.25*IT_0001*IT_0219;
    const ccomplex_t IT_0221 = 0.5*IT_0011 + (-0.5)*IT_0097 + 0.5*IT_0203;
    const ccomplex_t IT_0222 = 0.5*IT_0011 + 0.5*IT_0097 + (-0.5)*IT_0203;
    const ccomplex_t IT_0223 = IT_0221*IT_0222;
    const ccomplex_t IT_0224 = IT_0011*IT_0220;
    const ccomplex_t IT_0225 = (-0.5)*IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = IT_0001*IT_0226;
    const ccomplex_t IT_0228 = cabs(-IT_0203 + 0.25*IT_0001*IT_0219);
    const ccomplex_t IT_0229 = cpow(IT_0228, 0.5);
    const ccomplex_t IT_0230 = ( ((m_W
) >= (m_d + m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0218*(IT_0220 + IT_0227)*IT_0229
) : (0
));
    const ccomplex_t IT_0231 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0004;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = cpow(IT_0232, 2);
    const ccomplex_t IT_0234 = cpow(IT_0011 + -IT_0143 + IT_0203, 2);
    const ccomplex_t IT_0235 = -IT_0203 + 0.25*IT_0001*(IT_0011 + -IT_0143 +
       IT_0203)*(IT_0011 + IT_0143 + -IT_0203) + 0.25*IT_0001*IT_0234;
    const ccomplex_t IT_0236 = 0.5*IT_0011 + (-0.5)*IT_0143 + 0.5*IT_0203;
    const ccomplex_t IT_0237 = 0.5*IT_0011 + 0.5*IT_0143 + (-0.5)*IT_0203;
    const ccomplex_t IT_0238 = IT_0236*IT_0237;
    const ccomplex_t IT_0239 = IT_0011*IT_0235;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = IT_0238 + IT_0240;
    const ccomplex_t IT_0242 = IT_0001*IT_0241;
    const ccomplex_t IT_0243 = cabs(-IT_0203 + 0.25*IT_0001*IT_0234);
    const ccomplex_t IT_0244 = cpow(IT_0243, 0.5);
    const ccomplex_t IT_0245 = ( ((m_W
) >= (m_s + m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0233*(IT_0235 + IT_0242)*IT_0244
) : (0
));
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*m_W*e_em*IT_0004*(IT_0126
      *IT_0128 + IT_0127*IT_0129);
    const ccomplex_t IT_0247 = cpow(IT_0246, 2);
    const ccomplex_t IT_0248 = pow(m_H0, 2);
    const ccomplex_t IT_0249 = (-0.5)*IT_0248;
    const ccomplex_t IT_0250 = IT_0011 + IT_0249;
    const ccomplex_t IT_0251 = cpow(IT_0250, 2);
    const ccomplex_t IT_0252 = IT_0008*IT_0251;
    const ccomplex_t IT_0253 = cpow(IT_0011 + (-0.5)*IT_0248, 2);
    const ccomplex_t IT_0254 = cabs(-IT_0011 + IT_0001*IT_0253);
    const ccomplex_t IT_0255 = cpow(IT_0254, 0.5);
    const ccomplex_t IT_0256 = ( ((m_W
) >= (m_W + m_H0
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0247*(1 + 0.5*IT_0252)*IT_0255
) :
       (0
));
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*e_em*IT_0004*(IT_0126
      *IT_0128 + IT_0127*IT_0129);
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = -IT_0258;
    const ccomplex_t IT_0260 = cpow(IT_0259, 2);
    const ccomplex_t IT_0261 = pow(m_Hp, 2);
    const ccomplex_t IT_0262 = 0.5*IT_0011 + (-0.5)*IT_0132 + 0.5*IT_0261;
    const ccomplex_t IT_0263 = cpow(IT_0262, 2);
    const ccomplex_t IT_0264 = IT_0001*IT_0263;
    const ccomplex_t IT_0265 = cpow(IT_0011 + -IT_0132 + IT_0261, 2);
    const ccomplex_t IT_0266 = cabs(-IT_0261 + 0.25*IT_0001*IT_0265);
    const ccomplex_t IT_0267 = cpow(IT_0266, 0.5);
    const ccomplex_t IT_0268 = ( ((m_W
) >= (m_h + m_Hp
)) ? 
      (0.0416666666666667*IT_0000*IT_0001*IT_0260*(IT_0261 + -IT_0264)*IT_0267
)
       : (0
));
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*e_em*IT_0004*(IT_0126
      *IT_0127 + -IT_0128*IT_0129);
    const ccomplex_t IT_0270 = -IT_0269;
    const ccomplex_t IT_0271 = cpow(IT_0270, 2);
    const ccomplex_t IT_0272 = 0.5*IT_0011 + (-0.5)*IT_0248 + 0.5*IT_0261;
    const ccomplex_t IT_0273 = cpow(IT_0272, 2);
    const ccomplex_t IT_0274 = IT_0001*IT_0273;
    const ccomplex_t IT_0275 = cpow(IT_0011 + -IT_0248 + IT_0261, 2);
    const ccomplex_t IT_0276 = cabs(-IT_0261 + 0.25*IT_0001*IT_0275);
    const ccomplex_t IT_0277 = cpow(IT_0276, 0.5);
    const ccomplex_t IT_0278 = ( ((m_W
) >= (m_H0 + m_Hp
)) ? 
      (0.0416666666666667*IT_0000*IT_0001*IT_0271*(IT_0261 + -IT_0274)*IT_0277
)
       : (0
));
    const ccomplex_t IT_0279 = e_em*IT_0004;
    const ccomplex_t IT_0280 = cpow(IT_0279, 2);
    const ccomplex_t IT_0281 = pow(m_A0, 2);
    const ccomplex_t IT_0282 = 0.5*IT_0011 + (-0.5)*IT_0261 + 0.5*IT_0281;
    const ccomplex_t IT_0283 = cpow(IT_0282, 2);
    const ccomplex_t IT_0284 = IT_0001*IT_0283;
    const ccomplex_t IT_0285 = cpow(IT_0011 + -IT_0261 + IT_0281, 2);
    const ccomplex_t IT_0286 = cabs(-IT_0281 + 0.25*IT_0001*IT_0285);
    const ccomplex_t IT_0287 = cpow(IT_0286, 0.5);
    const ccomplex_t IT_0288 = ( ((m_W
) >= (m_A0 + m_Hp
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0280*(IT_0281 + -IT_0284)*IT_0287

      ) : (0
));
    const ccomplex_t IT_0289 = pow(m_mu, 2);
    const ccomplex_t IT_0290 = cpow(IT_0011 + -IT_0289, 2);
    const ccomplex_t IT_0291 = 0.25*IT_0001*(IT_0011 + IT_0289)*(IT_0011 + 
      -IT_0289) + 0.25*IT_0001*IT_0290;
    const ccomplex_t IT_0292 = 0.5*IT_0011 + (-0.5)*IT_0289;
    const ccomplex_t IT_0293 = 0.5*IT_0011 + 0.5*IT_0289;
    const ccomplex_t IT_0294 = IT_0292*IT_0293;
    const ccomplex_t IT_0295 = IT_0011*IT_0291;
    const ccomplex_t IT_0296 = (-0.5)*IT_0295;
    const ccomplex_t IT_0297 = IT_0294 + IT_0296;
    const ccomplex_t IT_0298 = IT_0001*IT_0297;
    const ccomplex_t IT_0299 = cabs(0.25*IT_0001*IT_0290);
    const ccomplex_t IT_0300 = cpow(IT_0299, 0.5);
    const ccomplex_t IT_0301 = ( ((m_W
) >= (m_mu
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0112*(IT_0291 + IT_0298)*IT_0300
) : (0
));
    const ccomplex_t IT_0302 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0303 = IT_0004*IT_0302;
    const ccomplex_t IT_0304 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0305 = IT_0004*IT_0304;
    const ccomplex_t IT_0306 = 1.4142135623731*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*(IT_0303 + 0.5*IT_0306);
    const ccomplex_t IT_0308 = -IT_0307;
    const ccomplex_t IT_0309 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0310 = IT_0004*IT_0309;
    const ccomplex_t IT_0311 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0312 = IT_0004*IT_0311;
    const ccomplex_t IT_0313 = 1.4142135623731*IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*(IT_0310 + (-0.5)*IT_0313);
    const ccomplex_t IT_0315 = -IT_0314;
    const ccomplex_t IT_0316 = m_C_1*m_N_1;
    const ccomplex_t IT_0317 = pow(m_N_1, 2);
    const ccomplex_t IT_0318 = pow(m_C_1, 2);
    const ccomplex_t IT_0319 = cpow(IT_0011 + IT_0317 + -IT_0318, 2);
    const ccomplex_t IT_0320 = -IT_0317 + 0.25*IT_0001*(IT_0011 + -IT_0317 +
       IT_0318)*(IT_0011 + IT_0317 + -IT_0318) + 0.25*IT_0001*IT_0319;
    const ccomplex_t IT_0321 = 0.5*IT_0011 + 0.5*IT_0317 + (-0.5)*IT_0318;
    const ccomplex_t IT_0322 = 0.5*IT_0011 + (-0.5)*IT_0317 + 0.5*IT_0318;
    const ccomplex_t IT_0323 = IT_0321*IT_0322;
    const ccomplex_t IT_0324 = IT_0011*IT_0320;
    const ccomplex_t IT_0325 = (-0.5)*IT_0324;
    const ccomplex_t IT_0326 = IT_0323 + IT_0325;
    const ccomplex_t IT_0327 = IT_0001*IT_0326;
    const ccomplex_t IT_0328 = IT_0320 + IT_0327;
    const ccomplex_t IT_0329 = cabs(-IT_0317 + 0.25*IT_0001*IT_0319);
    const ccomplex_t IT_0330 = cpow(IT_0329, 0.5);
    const ccomplex_t IT_0331 = ( ((m_W
) >= (m_C_1 + m_N_1
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0308)*(IT_0315*IT_0316 + 0.666666666666667*IT_0308
      *IT_0328) + conj(IT_0315)*(IT_0308*IT_0316 + 0.666666666666667*IT_0315
      *IT_0328))*IT_0330
) : (0
));
    const ccomplex_t IT_0332 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0333 = IT_0004*IT_0332;
    const ccomplex_t IT_0334 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0335 = IT_0004*IT_0334;
    const ccomplex_t IT_0336 = 1.4142135623731*IT_0335;
    const ccomplex_t IT_0337 = (0 + _Complex_I*1)*(IT_0333 + 0.5*IT_0336);
    const ccomplex_t IT_0338 = -IT_0337;
    const ccomplex_t IT_0339 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0340 = IT_0004*IT_0339;
    const ccomplex_t IT_0341 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0342 = IT_0004*IT_0341;
    const ccomplex_t IT_0343 = 1.4142135623731*IT_0342;
    const ccomplex_t IT_0344 = (0 + _Complex_I*1)*(IT_0340 + (-0.5)*IT_0343);
    const ccomplex_t IT_0345 = -IT_0344;
    const ccomplex_t IT_0346 = m_C_2*m_N_1;
    const ccomplex_t IT_0347 = pow(m_C_2, 2);
    const ccomplex_t IT_0348 = cpow(IT_0011 + IT_0317 + -IT_0347, 2);
    const ccomplex_t IT_0349 = -IT_0317 + 0.25*IT_0001*(IT_0011 + -IT_0317 +
       IT_0347)*(IT_0011 + IT_0317 + -IT_0347) + 0.25*IT_0001*IT_0348;
    const ccomplex_t IT_0350 = 0.5*IT_0011 + 0.5*IT_0317 + (-0.5)*IT_0347;
    const ccomplex_t IT_0351 = 0.5*IT_0011 + (-0.5)*IT_0317 + 0.5*IT_0347;
    const ccomplex_t IT_0352 = IT_0350*IT_0351;
    const ccomplex_t IT_0353 = IT_0011*IT_0349;
    const ccomplex_t IT_0354 = (-0.5)*IT_0353;
    const ccomplex_t IT_0355 = IT_0352 + IT_0354;
    const ccomplex_t IT_0356 = IT_0001*IT_0355;
    const ccomplex_t IT_0357 = IT_0349 + IT_0356;
    const ccomplex_t IT_0358 = cabs(-IT_0317 + 0.25*IT_0001*IT_0348);
    const ccomplex_t IT_0359 = cpow(IT_0358, 0.5);
    const ccomplex_t IT_0360 = ( ((m_W
) >= (m_C_2 + m_N_1
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0338)*(IT_0345*IT_0346 + 0.666666666666667*IT_0338
      *IT_0357) + conj(IT_0345)*(IT_0338*IT_0346 + 0.666666666666667*IT_0345
      *IT_0357))*IT_0359
) : (0
));
    const ccomplex_t IT_0361 = conj(N_W2)*e_em*U_Wm1;
    const ccomplex_t IT_0362 = IT_0004*IT_0361;
    const ccomplex_t IT_0363 = conj(N_d2)*U_d1*e_em;
    const ccomplex_t IT_0364 = IT_0004*IT_0363;
    const ccomplex_t IT_0365 = 1.4142135623731*IT_0364;
    const ccomplex_t IT_0366 = (0 + _Complex_I*1)*(IT_0362 + 0.5*IT_0365);
    const ccomplex_t IT_0367 = -IT_0366;
    const ccomplex_t IT_0368 = N_W2*e_em*conj(V_Wp1);
    const ccomplex_t IT_0369 = IT_0004*IT_0368;
    const ccomplex_t IT_0370 = N_u2*conj(V_u1)*e_em;
    const ccomplex_t IT_0371 = IT_0004*IT_0370;
    const ccomplex_t IT_0372 = 1.4142135623731*IT_0371;
    const ccomplex_t IT_0373 = (0 + _Complex_I*1)*(IT_0369 + (-0.5)*IT_0372);
    const ccomplex_t IT_0374 = -IT_0373;
    const ccomplex_t IT_0375 = m_C_1*m_N_2;
    const ccomplex_t IT_0376 = pow(m_N_2, 2);
    const ccomplex_t IT_0377 = cpow(IT_0011 + -IT_0318 + IT_0376, 2);
    const ccomplex_t IT_0378 = -IT_0376 + 0.25*IT_0001*(IT_0011 + -IT_0318 +
       IT_0376)*(IT_0011 + IT_0318 + -IT_0376) + 0.25*IT_0001*IT_0377;
    const ccomplex_t IT_0379 = 0.5*IT_0011 + (-0.5)*IT_0318 + 0.5*IT_0376;
    const ccomplex_t IT_0380 = 0.5*IT_0011 + 0.5*IT_0318 + (-0.5)*IT_0376;
    const ccomplex_t IT_0381 = IT_0379*IT_0380;
    const ccomplex_t IT_0382 = IT_0011*IT_0378;
    const ccomplex_t IT_0383 = (-0.5)*IT_0382;
    const ccomplex_t IT_0384 = IT_0381 + IT_0383;
    const ccomplex_t IT_0385 = IT_0001*IT_0384;
    const ccomplex_t IT_0386 = IT_0378 + IT_0385;
    const ccomplex_t IT_0387 = cabs(-IT_0376 + 0.25*IT_0001*IT_0377);
    const ccomplex_t IT_0388 = cpow(IT_0387, 0.5);
    const ccomplex_t IT_0389 = ( ((m_W
) >= (m_C_1 + m_N_2
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0367)*(IT_0374*IT_0375 + 0.666666666666667*IT_0367
      *IT_0386) + conj(IT_0374)*(IT_0367*IT_0375 + 0.666666666666667*IT_0374
      *IT_0386))*IT_0388
) : (0
));
    const ccomplex_t IT_0390 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0391 = IT_0004*IT_0390;
    const ccomplex_t IT_0392 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0393 = IT_0004*IT_0392;
    const ccomplex_t IT_0394 = 1.4142135623731*IT_0393;
    const ccomplex_t IT_0395 = (0 + _Complex_I*1)*(IT_0391 + 0.5*IT_0394);
    const ccomplex_t IT_0396 = -IT_0395;
    const ccomplex_t IT_0397 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0398 = IT_0004*IT_0397;
    const ccomplex_t IT_0399 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0400 = IT_0004*IT_0399;
    const ccomplex_t IT_0401 = 1.4142135623731*IT_0400;
    const ccomplex_t IT_0402 = (0 + _Complex_I*1)*(IT_0398 + (-0.5)*IT_0401);
    const ccomplex_t IT_0403 = -IT_0402;
    const ccomplex_t IT_0404 = m_C_2*m_N_2;
    const ccomplex_t IT_0405 = cpow(IT_0011 + -IT_0347 + IT_0376, 2);
    const ccomplex_t IT_0406 = -IT_0376 + 0.25*IT_0001*(IT_0011 + -IT_0347 +
       IT_0376)*(IT_0011 + IT_0347 + -IT_0376) + 0.25*IT_0001*IT_0405;
    const ccomplex_t IT_0407 = 0.5*IT_0011 + (-0.5)*IT_0347 + 0.5*IT_0376;
    const ccomplex_t IT_0408 = 0.5*IT_0011 + 0.5*IT_0347 + (-0.5)*IT_0376;
    const ccomplex_t IT_0409 = IT_0407*IT_0408;
    const ccomplex_t IT_0410 = IT_0011*IT_0406;
    const ccomplex_t IT_0411 = (-0.5)*IT_0410;
    const ccomplex_t IT_0412 = IT_0409 + IT_0411;
    const ccomplex_t IT_0413 = IT_0001*IT_0412;
    const ccomplex_t IT_0414 = IT_0406 + IT_0413;
    const ccomplex_t IT_0415 = cabs(-IT_0376 + 0.25*IT_0001*IT_0405);
    const ccomplex_t IT_0416 = cpow(IT_0415, 0.5);
    const ccomplex_t IT_0417 = ( ((m_W
) >= (m_C_2 + m_N_2
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0396)*(IT_0403*IT_0404 + 0.666666666666667*IT_0396
      *IT_0414) + conj(IT_0403)*(IT_0396*IT_0404 + 0.666666666666667*IT_0403
      *IT_0414))*IT_0416
) : (0
));
    const ccomplex_t IT_0418 = conj(N_W3)*e_em*U_Wm1;
    const ccomplex_t IT_0419 = IT_0004*IT_0418;
    const ccomplex_t IT_0420 = conj(N_d3)*U_d1*e_em;
    const ccomplex_t IT_0421 = IT_0004*IT_0420;
    const ccomplex_t IT_0422 = 1.4142135623731*IT_0421;
    const ccomplex_t IT_0423 = (0 + _Complex_I*1)*(IT_0419 + 0.5*IT_0422);
    const ccomplex_t IT_0424 = -IT_0423;
    const ccomplex_t IT_0425 = N_W3*e_em*conj(V_Wp1);
    const ccomplex_t IT_0426 = IT_0004*IT_0425;
    const ccomplex_t IT_0427 = N_u3*conj(V_u1)*e_em;
    const ccomplex_t IT_0428 = IT_0004*IT_0427;
    const ccomplex_t IT_0429 = 1.4142135623731*IT_0428;
    const ccomplex_t IT_0430 = (0 + _Complex_I*1)*(IT_0426 + (-0.5)*IT_0429);
    const ccomplex_t IT_0431 = -IT_0430;
    const ccomplex_t IT_0432 = m_C_1*m_N_3;
    const ccomplex_t IT_0433 = pow(m_N_3, 2);
    const ccomplex_t IT_0434 = cpow(IT_0011 + -IT_0318 + IT_0433, 2);
    const ccomplex_t IT_0435 = -IT_0433 + 0.25*IT_0001*(IT_0011 + -IT_0318 +
       IT_0433)*(IT_0011 + IT_0318 + -IT_0433) + 0.25*IT_0001*IT_0434;
    const ccomplex_t IT_0436 = 0.5*IT_0011 + (-0.5)*IT_0318 + 0.5*IT_0433;
    const ccomplex_t IT_0437 = 0.5*IT_0011 + 0.5*IT_0318 + (-0.5)*IT_0433;
    const ccomplex_t IT_0438 = IT_0436*IT_0437;
    const ccomplex_t IT_0439 = IT_0011*IT_0435;
    const ccomplex_t IT_0440 = (-0.5)*IT_0439;
    const ccomplex_t IT_0441 = IT_0438 + IT_0440;
    const ccomplex_t IT_0442 = IT_0001*IT_0441;
    const ccomplex_t IT_0443 = IT_0435 + IT_0442;
    const ccomplex_t IT_0444 = cabs(-IT_0433 + 0.25*IT_0001*IT_0434);
    const ccomplex_t IT_0445 = cpow(IT_0444, 0.5);
    const ccomplex_t IT_0446 = ( ((m_W
) >= (m_C_1 + m_N_3
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0424)*(IT_0431*IT_0432 + 0.666666666666667*IT_0424
      *IT_0443) + conj(IT_0431)*(IT_0424*IT_0432 + 0.666666666666667*IT_0431
      *IT_0443))*IT_0445
) : (0
));
    const ccomplex_t IT_0447 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0448 = IT_0004*IT_0447;
    const ccomplex_t IT_0449 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0450 = IT_0004*IT_0449;
    const ccomplex_t IT_0451 = 1.4142135623731*IT_0450;
    const ccomplex_t IT_0452 = (0 + _Complex_I*1)*(IT_0448 + 0.5*IT_0451);
    const ccomplex_t IT_0453 = -IT_0452;
    const ccomplex_t IT_0454 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0455 = IT_0004*IT_0454;
    const ccomplex_t IT_0456 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0457 = IT_0004*IT_0456;
    const ccomplex_t IT_0458 = 1.4142135623731*IT_0457;
    const ccomplex_t IT_0459 = (0 + _Complex_I*1)*(IT_0455 + (-0.5)*IT_0458);
    const ccomplex_t IT_0460 = -IT_0459;
    const ccomplex_t IT_0461 = m_C_2*m_N_3;
    const ccomplex_t IT_0462 = cpow(IT_0011 + -IT_0347 + IT_0433, 2);
    const ccomplex_t IT_0463 = -IT_0433 + 0.25*IT_0001*(IT_0011 + -IT_0347 +
       IT_0433)*(IT_0011 + IT_0347 + -IT_0433) + 0.25*IT_0001*IT_0462;
    const ccomplex_t IT_0464 = 0.5*IT_0011 + (-0.5)*IT_0347 + 0.5*IT_0433;
    const ccomplex_t IT_0465 = 0.5*IT_0011 + 0.5*IT_0347 + (-0.5)*IT_0433;
    const ccomplex_t IT_0466 = IT_0464*IT_0465;
    const ccomplex_t IT_0467 = IT_0011*IT_0463;
    const ccomplex_t IT_0468 = (-0.5)*IT_0467;
    const ccomplex_t IT_0469 = IT_0466 + IT_0468;
    const ccomplex_t IT_0470 = IT_0001*IT_0469;
    const ccomplex_t IT_0471 = IT_0463 + IT_0470;
    const ccomplex_t IT_0472 = cabs(-IT_0433 + 0.25*IT_0001*IT_0462);
    const ccomplex_t IT_0473 = cpow(IT_0472, 0.5);
    const ccomplex_t IT_0474 = ( ((m_W
) >= (m_C_2 + m_N_3
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0453)*(IT_0460*IT_0461 + 0.666666666666667*IT_0453
      *IT_0471) + conj(IT_0460)*(IT_0453*IT_0461 + 0.666666666666667*IT_0460
      *IT_0471))*IT_0473
) : (0
));
    const ccomplex_t IT_0475 = conj(N_W4)*e_em*U_Wm1;
    const ccomplex_t IT_0476 = IT_0004*IT_0475;
    const ccomplex_t IT_0477 = conj(N_d4)*U_d1*e_em;
    const ccomplex_t IT_0478 = IT_0004*IT_0477;
    const ccomplex_t IT_0479 = 1.4142135623731*IT_0478;
    const ccomplex_t IT_0480 = (0 + _Complex_I*1)*(IT_0476 + 0.5*IT_0479);
    const ccomplex_t IT_0481 = -IT_0480;
    const ccomplex_t IT_0482 = N_W4*e_em*conj(V_Wp1);
    const ccomplex_t IT_0483 = IT_0004*IT_0482;
    const ccomplex_t IT_0484 = N_u4*conj(V_u1)*e_em;
    const ccomplex_t IT_0485 = IT_0004*IT_0484;
    const ccomplex_t IT_0486 = 1.4142135623731*IT_0485;
    const ccomplex_t IT_0487 = (0 + _Complex_I*1)*(IT_0483 + (-0.5)*IT_0486);
    const ccomplex_t IT_0488 = -IT_0487;
    const ccomplex_t IT_0489 = m_C_1*m_N_4;
    const ccomplex_t IT_0490 = pow(m_N_4, 2);
    const ccomplex_t IT_0491 = cpow(IT_0011 + -IT_0318 + IT_0490, 2);
    const ccomplex_t IT_0492 = -IT_0490 + 0.25*IT_0001*(IT_0011 + -IT_0318 +
       IT_0490)*(IT_0011 + IT_0318 + -IT_0490) + 0.25*IT_0001*IT_0491;
    const ccomplex_t IT_0493 = 0.5*IT_0011 + (-0.5)*IT_0318 + 0.5*IT_0490;
    const ccomplex_t IT_0494 = 0.5*IT_0011 + 0.5*IT_0318 + (-0.5)*IT_0490;
    const ccomplex_t IT_0495 = IT_0493*IT_0494;
    const ccomplex_t IT_0496 = IT_0011*IT_0492;
    const ccomplex_t IT_0497 = (-0.5)*IT_0496;
    const ccomplex_t IT_0498 = IT_0495 + IT_0497;
    const ccomplex_t IT_0499 = IT_0001*IT_0498;
    const ccomplex_t IT_0500 = IT_0492 + IT_0499;
    const ccomplex_t IT_0501 = cabs(-IT_0490 + 0.25*IT_0001*IT_0491);
    const ccomplex_t IT_0502 = cpow(IT_0501, 0.5);
    const ccomplex_t IT_0503 = ( ((m_W
) >= (m_C_1 + m_N_4
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0481)*(IT_0488*IT_0489 + 0.666666666666667*IT_0481
      *IT_0500) + conj(IT_0488)*(IT_0481*IT_0489 + 0.666666666666667*IT_0488
      *IT_0500))*IT_0502
) : (0
));
    const ccomplex_t IT_0504 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0505 = IT_0004*IT_0504;
    const ccomplex_t IT_0506 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0507 = IT_0004*IT_0506;
    const ccomplex_t IT_0508 = 1.4142135623731*IT_0507;
    const ccomplex_t IT_0509 = (0 + _Complex_I*1)*(IT_0505 + 0.5*IT_0508);
    const ccomplex_t IT_0510 = -IT_0509;
    const ccomplex_t IT_0511 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0512 = IT_0004*IT_0511;
    const ccomplex_t IT_0513 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0514 = IT_0004*IT_0513;
    const ccomplex_t IT_0515 = 1.4142135623731*IT_0514;
    const ccomplex_t IT_0516 = (0 + _Complex_I*1)*(IT_0512 + (-0.5)*IT_0515);
    const ccomplex_t IT_0517 = -IT_0516;
    const ccomplex_t IT_0518 = m_C_2*m_N_4;
    const ccomplex_t IT_0519 = cpow(IT_0011 + -IT_0347 + IT_0490, 2);
    const ccomplex_t IT_0520 = -IT_0490 + 0.25*IT_0001*(IT_0011 + -IT_0347 +
       IT_0490)*(IT_0011 + IT_0347 + -IT_0490) + 0.25*IT_0001*IT_0519;
    const ccomplex_t IT_0521 = 0.5*IT_0011 + (-0.5)*IT_0347 + 0.5*IT_0490;
    const ccomplex_t IT_0522 = 0.5*IT_0011 + 0.5*IT_0347 + (-0.5)*IT_0490;
    const ccomplex_t IT_0523 = IT_0521*IT_0522;
    const ccomplex_t IT_0524 = IT_0011*IT_0520;
    const ccomplex_t IT_0525 = (-0.5)*IT_0524;
    const ccomplex_t IT_0526 = IT_0523 + IT_0525;
    const ccomplex_t IT_0527 = IT_0001*IT_0526;
    const ccomplex_t IT_0528 = IT_0520 + IT_0527;
    const ccomplex_t IT_0529 = cabs(-IT_0490 + 0.25*IT_0001*IT_0519);
    const ccomplex_t IT_0530 = cpow(IT_0529, 0.5);
    const ccomplex_t IT_0531 = ( ((m_W
) >= (m_C_2 + m_N_4
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0510)*(IT_0517*IT_0518 + 0.666666666666667*IT_0510
      *IT_0528) + conj(IT_0517)*(IT_0510*IT_0518 + 0.666666666666667*IT_0517
      *IT_0528))*IT_0530
) : (0
));
    const ccomplex_t IT_0532 = pow(m_tau, 2);
    const ccomplex_t IT_0533 = cpow(IT_0011 + -IT_0532, 2);
    const ccomplex_t IT_0534 = 0.25*IT_0001*(IT_0011 + IT_0532)*(IT_0011 + 
      -IT_0532) + 0.25*IT_0001*IT_0533;
    const ccomplex_t IT_0535 = 0.5*IT_0011 + (-0.5)*IT_0532;
    const ccomplex_t IT_0536 = 0.5*IT_0011 + 0.5*IT_0532;
    const ccomplex_t IT_0537 = IT_0535*IT_0536;
    const ccomplex_t IT_0538 = IT_0011*IT_0534;
    const ccomplex_t IT_0539 = (-0.5)*IT_0538;
    const ccomplex_t IT_0540 = IT_0537 + IT_0539;
    const ccomplex_t IT_0541 = IT_0001*IT_0540;
    const ccomplex_t IT_0542 = cabs(0.25*IT_0001*IT_0533);
    const ccomplex_t IT_0543 = cpow(IT_0542, 0.5);
    const ccomplex_t IT_0544 = ( ((m_W
) >= (m_tau
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0112*(IT_0534 + IT_0541)*IT_0543
) : (0
));
    const ccomplex_t IT_0545 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0004*U_sb_00;
    const ccomplex_t IT_0546 = -IT_0545;
    const ccomplex_t IT_0547 = -IT_0546;
    const ccomplex_t IT_0548 = pow(m_sc_L, 2);
    const ccomplex_t IT_0549 = pow(m_sb_1, 2);
    const ccomplex_t IT_0550 = 0.5*IT_0011 + 0.5*IT_0548 + (-0.5)*IT_0549;
    const ccomplex_t IT_0551 = cpow(IT_0550, 2);
    const ccomplex_t IT_0552 = IT_0001*IT_0551;
    const ccomplex_t IT_0553 = cpow(IT_0011 + IT_0548 + -IT_0549, 2);
    const ccomplex_t IT_0554 = cabs(-IT_0548 + 0.25*IT_0001*IT_0553);
    const ccomplex_t IT_0555 = cpow(IT_0554, 0.5);
    const ccomplex_t IT_0556 = ( ((m_W
) >= (m_sb_1 + m_sc_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0547*conj(IT_0547)*(IT_0548 + -IT_0552)*IT_0555
) : 
      (0
));
    const ccomplex_t IT_0557 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0004*U_sb_01;
    const ccomplex_t IT_0558 = -IT_0557;
    const ccomplex_t IT_0559 = -IT_0558;
    const ccomplex_t IT_0560 = pow(m_sb_2, 2);
    const ccomplex_t IT_0561 = 0.5*IT_0011 + 0.5*IT_0548 + (-0.5)*IT_0560;
    const ccomplex_t IT_0562 = cpow(IT_0561, 2);
    const ccomplex_t IT_0563 = IT_0001*IT_0562;
    const ccomplex_t IT_0564 = cpow(IT_0011 + IT_0548 + -IT_0560, 2);
    const ccomplex_t IT_0565 = cabs(-IT_0548 + 0.25*IT_0001*IT_0564);
    const ccomplex_t IT_0566 = cpow(IT_0565, 0.5);
    const ccomplex_t IT_0567 = ( ((m_W
) >= (m_sb_2 + m_sc_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0559*conj(IT_0559)*(IT_0548 + -IT_0563)*IT_0566
) : 
      (0
));
    const ccomplex_t IT_0568 = -IT_0095;
    const ccomplex_t IT_0569 = -IT_0568;
    const ccomplex_t IT_0570 = pow(m_sd_L, 2);
    const ccomplex_t IT_0571 = 0.5*IT_0011 + 0.5*IT_0548 + (-0.5)*IT_0570;
    const ccomplex_t IT_0572 = cpow(IT_0571, 2);
    const ccomplex_t IT_0573 = IT_0001*IT_0572;
    const ccomplex_t IT_0574 = cpow(IT_0011 + IT_0548 + -IT_0570, 2);
    const ccomplex_t IT_0575 = cabs(-IT_0548 + 0.25*IT_0001*IT_0574);
    const ccomplex_t IT_0576 = cpow(IT_0575, 0.5);
    const ccomplex_t IT_0577 = ( ((m_W
) >= (m_sc_L + m_sd_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0569*conj(IT_0569)*(IT_0548 + -IT_0573)*IT_0576
) : 
      (0
));
    const ccomplex_t IT_0578 = -IT_0141;
    const ccomplex_t IT_0579 = -IT_0578;
    const ccomplex_t IT_0580 = pow(m_ss_L, 2);
    const ccomplex_t IT_0581 = 0.5*IT_0011 + 0.5*IT_0548 + (-0.5)*IT_0580;
    const ccomplex_t IT_0582 = cpow(IT_0581, 2);
    const ccomplex_t IT_0583 = IT_0001*IT_0582;
    const ccomplex_t IT_0584 = cpow(IT_0011 + IT_0548 + -IT_0580, 2);
    const ccomplex_t IT_0585 = cabs(-IT_0548 + 0.25*IT_0001*IT_0584);
    const ccomplex_t IT_0586 = cpow(IT_0585, 0.5);
    const ccomplex_t IT_0587 = ( ((m_W
) >= (m_sc_L + m_ss_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0579*conj(IT_0579)*(IT_0548 + -IT_0583)*IT_0586
) : 
      (0
));
    const ccomplex_t IT_0588 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_00*conj(U_st_00);
    const ccomplex_t IT_0589 = -IT_0588;
    const ccomplex_t IT_0590 = -IT_0589;
    const ccomplex_t IT_0591 = pow(m_st_1, 2);
    const ccomplex_t IT_0592 = 0.5*IT_0011 + (-0.5)*IT_0549 + 0.5*IT_0591;
    const ccomplex_t IT_0593 = cpow(IT_0592, 2);
    const ccomplex_t IT_0594 = IT_0001*IT_0593;
    const ccomplex_t IT_0595 = cpow(IT_0011 + -IT_0549 + IT_0591, 2);
    const ccomplex_t IT_0596 = cabs(-IT_0591 + 0.25*IT_0001*IT_0595);
    const ccomplex_t IT_0597 = cpow(IT_0596, 0.5);
    const ccomplex_t IT_0598 = ( ((m_W
) >= (m_sb_1 + m_st_1
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0590*conj(IT_0590)*(IT_0591 + -IT_0594)*IT_0597
) : 
      (0
));
    const ccomplex_t IT_0599 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_01*conj(U_st_00);
    const ccomplex_t IT_0600 = -IT_0599;
    const ccomplex_t IT_0601 = -IT_0600;
    const ccomplex_t IT_0602 = 0.5*IT_0011 + (-0.5)*IT_0560 + 0.5*IT_0591;
    const ccomplex_t IT_0603 = cpow(IT_0602, 2);
    const ccomplex_t IT_0604 = IT_0001*IT_0603;
    const ccomplex_t IT_0605 = cpow(IT_0011 + -IT_0560 + IT_0591, 2);
    const ccomplex_t IT_0606 = cabs(-IT_0591 + 0.25*IT_0001*IT_0605);
    const ccomplex_t IT_0607 = cpow(IT_0606, 0.5);
    const ccomplex_t IT_0608 = ( ((m_W
) >= (m_sb_2 + m_st_1
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0601*conj(IT_0601)*(IT_0591 + -IT_0604)*IT_0607
) : 
      (0
));
    const ccomplex_t IT_0609 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0004*conj(U_st_00);
    const ccomplex_t IT_0610 = -IT_0609;
    const ccomplex_t IT_0611 = 0.5*IT_0011 + 0.5*IT_0570 + (-0.5)*IT_0591;
    const ccomplex_t IT_0612 = cpow(IT_0611, 2);
    const ccomplex_t IT_0613 = IT_0001*IT_0612;
    const ccomplex_t IT_0614 = cpow(IT_0011 + IT_0570 + -IT_0591, 2);
    const ccomplex_t IT_0615 = cabs(-IT_0570 + 0.25*IT_0001*IT_0614);
    const ccomplex_t IT_0616 = cpow(IT_0615, 0.5);
    const ccomplex_t IT_0617 = ( ((m_W
) >= (m_sd_L + m_st_1
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0610*conj(IT_0610)*(IT_0570 + -IT_0613)*IT_0616
) : 
      (0
));
    const ccomplex_t IT_0618 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0004*conj(U_st_00);
    const ccomplex_t IT_0619 = -IT_0618;
    const ccomplex_t IT_0620 = 0.5*IT_0011 + 0.5*IT_0580 + (-0.5)*IT_0591;
    const ccomplex_t IT_0621 = cpow(IT_0620, 2);
    const ccomplex_t IT_0622 = IT_0001*IT_0621;
    const ccomplex_t IT_0623 = cpow(IT_0011 + IT_0580 + -IT_0591, 2);
    const ccomplex_t IT_0624 = cabs(-IT_0580 + 0.25*IT_0001*IT_0623);
    const ccomplex_t IT_0625 = cpow(IT_0624, 0.5);
    const ccomplex_t IT_0626 = ( ((m_W
) >= (m_ss_L + m_st_1
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0619*conj(IT_0619)*(IT_0580 + -IT_0622)*IT_0625
) : 
      (0
));
    const ccomplex_t IT_0627 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_00*conj(U_st_01);
    const ccomplex_t IT_0628 = -IT_0627;
    const ccomplex_t IT_0629 = -IT_0628;
    const ccomplex_t IT_0630 = pow(m_st_2, 2);
    const ccomplex_t IT_0631 = 0.5*IT_0011 + (-0.5)*IT_0549 + 0.5*IT_0630;
    const ccomplex_t IT_0632 = cpow(IT_0631, 2);
    const ccomplex_t IT_0633 = IT_0001*IT_0632;
    const ccomplex_t IT_0634 = cpow(IT_0011 + -IT_0549 + IT_0630, 2);
    const ccomplex_t IT_0635 = cabs(-IT_0630 + 0.25*IT_0001*IT_0634);
    const ccomplex_t IT_0636 = cpow(IT_0635, 0.5);
    const ccomplex_t IT_0637 = ( ((m_W
) >= (m_sb_1 + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0629*conj(IT_0629)*(IT_0630 + -IT_0633)*IT_0636
) : 
      (0
));
    const ccomplex_t IT_0638 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_01*conj(U_st_01);
    const ccomplex_t IT_0639 = -IT_0638;
    const ccomplex_t IT_0640 = -IT_0639;
    const ccomplex_t IT_0641 = 0.5*IT_0011 + (-0.5)*IT_0560 + 0.5*IT_0630;
    const ccomplex_t IT_0642 = cpow(IT_0641, 2);
    const ccomplex_t IT_0643 = IT_0001*IT_0642;
    const ccomplex_t IT_0644 = cpow(IT_0011 + -IT_0560 + IT_0630, 2);
    const ccomplex_t IT_0645 = cabs(-IT_0630 + 0.25*IT_0001*IT_0644);
    const ccomplex_t IT_0646 = cpow(IT_0645, 0.5);
    const ccomplex_t IT_0647 = ( ((m_W
) >= (m_sb_2 + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0640*conj(IT_0640)*(IT_0630 + -IT_0643)*IT_0646
) : 
      (0
));
    const ccomplex_t IT_0648 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0004*conj(U_st_01);
    const ccomplex_t IT_0649 = -IT_0648;
    const ccomplex_t IT_0650 = 0.5*IT_0011 + 0.5*IT_0570 + (-0.5)*IT_0630;
    const ccomplex_t IT_0651 = cpow(IT_0650, 2);
    const ccomplex_t IT_0652 = IT_0001*IT_0651;
    const ccomplex_t IT_0653 = cpow(IT_0011 + IT_0570 + -IT_0630, 2);
    const ccomplex_t IT_0654 = cabs(-IT_0570 + 0.25*IT_0001*IT_0653);
    const ccomplex_t IT_0655 = cpow(IT_0654, 0.5);
    const ccomplex_t IT_0656 = ( ((m_W
) >= (m_sd_L + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0649*conj(IT_0649)*(IT_0570 + -IT_0652)*IT_0655
) : 
      (0
));
    const ccomplex_t IT_0657 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0004*conj(U_st_01);
    const ccomplex_t IT_0658 = -IT_0657;
    const ccomplex_t IT_0659 = 0.5*IT_0011 + 0.5*IT_0580 + (-0.5)*IT_0630;
    const ccomplex_t IT_0660 = cpow(IT_0659, 2);
    const ccomplex_t IT_0661 = IT_0001*IT_0660;
    const ccomplex_t IT_0662 = cpow(IT_0011 + IT_0580 + -IT_0630, 2);
    const ccomplex_t IT_0663 = cabs(-IT_0580 + 0.25*IT_0001*IT_0662);
    const ccomplex_t IT_0664 = cpow(IT_0663, 0.5);
    const ccomplex_t IT_0665 = ( ((m_W
) >= (m_ss_L + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0658*conj(IT_0658)*(IT_0580 + -IT_0661)*IT_0664
) : 
      (0
));
    const ccomplex_t IT_0666 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *IT_0199*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0667 = -IT_0666;
    const ccomplex_t IT_0668 = -IT_0667;
    const ccomplex_t IT_0669 = pow(m_su_L, 2);
    const ccomplex_t IT_0670 = 0.5*IT_0011 + (-0.5)*IT_0549 + 0.5*IT_0669;
    const ccomplex_t IT_0671 = cpow(IT_0670, 2);
    const ccomplex_t IT_0672 = IT_0001*IT_0671;
    const ccomplex_t IT_0673 = cpow(IT_0011 + -IT_0549 + IT_0669, 2);
    const ccomplex_t IT_0674 = cabs(-IT_0669 + 0.25*IT_0001*IT_0673);
    const ccomplex_t IT_0675 = cpow(IT_0674, 0.5);
    const ccomplex_t IT_0676 = ( ((m_W
) >= (m_sb_1 + m_su_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0668*conj(IT_0668)*(IT_0669 + -IT_0672)*IT_0675
) : 
      (0
));
    const ccomplex_t IT_0677 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *IT_0199*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0678 = -IT_0677;
    const ccomplex_t IT_0679 = -IT_0678;
    const ccomplex_t IT_0680 = 0.5*IT_0011 + (-0.5)*IT_0560 + 0.5*IT_0669;
    const ccomplex_t IT_0681 = cpow(IT_0680, 2);
    const ccomplex_t IT_0682 = IT_0001*IT_0681;
    const ccomplex_t IT_0683 = cpow(IT_0011 + -IT_0560 + IT_0669, 2);
    const ccomplex_t IT_0684 = cabs(-IT_0669 + 0.25*IT_0001*IT_0683);
    const ccomplex_t IT_0685 = cpow(IT_0684, 0.5);
    const ccomplex_t IT_0686 = ( ((m_W
) >= (m_sb_2 + m_su_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0679*conj(IT_0679)*(IT_0669 + -IT_0682)*IT_0685
) : 
      (0
));
    const ccomplex_t IT_0687 = -IT_0216;
    const ccomplex_t IT_0688 = -IT_0687;
    const ccomplex_t IT_0689 = cpow(IT_0688, 2);
    const ccomplex_t IT_0690 = 0.5*IT_0011 + (-0.5)*IT_0570 + 0.5*IT_0669;
    const ccomplex_t IT_0691 = cpow(IT_0690, 2);
    const ccomplex_t IT_0692 = IT_0001*IT_0691;
    const ccomplex_t IT_0693 = cpow(IT_0011 + -IT_0570 + IT_0669, 2);
    const ccomplex_t IT_0694 = cabs(-IT_0669 + 0.25*IT_0001*IT_0693);
    const ccomplex_t IT_0695 = cpow(IT_0694, 0.5);
    const ccomplex_t IT_0696 = ( ((m_W
) >= (m_sd_L + m_su_L
)) ? (0.125
      *IT_0000*IT_0001*IT_0689*(IT_0669 + -IT_0692)*IT_0695
) : (0
));
    const ccomplex_t IT_0697 = -IT_0231;
    const ccomplex_t IT_0698 = -IT_0697;
    const ccomplex_t IT_0699 = cpow(IT_0698, 2);
    const ccomplex_t IT_0700 = 0.5*IT_0011 + (-0.5)*IT_0580 + 0.5*IT_0669;
    const ccomplex_t IT_0701 = cpow(IT_0700, 2);
    const ccomplex_t IT_0702 = IT_0001*IT_0701;
    const ccomplex_t IT_0703 = cpow(IT_0011 + -IT_0580 + IT_0669, 2);
    const ccomplex_t IT_0704 = cabs(-IT_0669 + 0.25*IT_0001*IT_0703);
    const ccomplex_t IT_0705 = cpow(IT_0704, 0.5);
    const ccomplex_t IT_0706 = ( ((m_W
) >= (m_ss_L + m_su_L
)) ? (0.125
      *IT_0000*IT_0001*IT_0699*(IT_0669 + -IT_0702)*IT_0705
) : (0
));
    const ccomplex_t IT_0707 = -IT_0110;
    const ccomplex_t IT_0708 = cpow(IT_0707, 2);
    const ccomplex_t IT_0709 = pow(m_se_L, 2);
    const ccomplex_t IT_0710 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0711 = 0.5*IT_0011 + 0.5*IT_0709 + (-0.5)*IT_0710;
    const ccomplex_t IT_0712 = cpow(IT_0711, 2);
    const ccomplex_t IT_0713 = IT_0001*IT_0712;
    const ccomplex_t IT_0714 = cpow(IT_0011 + IT_0709 + -IT_0710, 2);
    const ccomplex_t IT_0715 = cabs(-IT_0709 + 0.25*IT_0001*IT_0714);
    const ccomplex_t IT_0716 = cpow(IT_0715, 0.5);
    const ccomplex_t IT_0717 = ( ((m_W
) >= (m_se_L + m_snu_mu
)) ? 
      (0.0416666666666667*IT_0000*IT_0001*IT_0708*(IT_0709 + -IT_0713)*IT_0716
)
       : (0
));
    const ccomplex_t IT_0718 = pow(m_smu_L, 2);
    const ccomplex_t IT_0719 = 0.5*IT_0011 + (-0.5)*IT_0710 + 0.5*IT_0718;
    const ccomplex_t IT_0720 = cpow(IT_0719, 2);
    const ccomplex_t IT_0721 = IT_0001*IT_0720;
    const ccomplex_t IT_0722 = cpow(IT_0011 + -IT_0710 + IT_0718, 2);
    const ccomplex_t IT_0723 = cabs(-IT_0718 + 0.25*IT_0001*IT_0722);
    const ccomplex_t IT_0724 = cpow(IT_0723, 0.5);
    const ccomplex_t IT_0725 = ( ((m_W
) >= (m_smu_L + m_snu_mu
)) ? 
      (0.0416666666666667*IT_0000*IT_0001*IT_0708*(IT_0718 + -IT_0721)*IT_0724
)
       : (0
));
    const ccomplex_t IT_0726 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *U_stau_00;
    const ccomplex_t IT_0727 = -IT_0726;
    const ccomplex_t IT_0728 = -IT_0727;
    const ccomplex_t IT_0729 = pow(m_stau_1, 2);
    const ccomplex_t IT_0730 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0731 = 0.5*IT_0011 + (-0.5)*IT_0729 + 0.5*IT_0730;
    const ccomplex_t IT_0732 = cpow(IT_0731, 2);
    const ccomplex_t IT_0733 = IT_0001*IT_0732;
    const ccomplex_t IT_0734 = cpow(IT_0011 + -IT_0729 + IT_0730, 2);
    const ccomplex_t IT_0735 = cabs(-IT_0730 + 0.25*IT_0001*IT_0734);
    const ccomplex_t IT_0736 = cpow(IT_0735, 0.5);
    const ccomplex_t IT_0737 = ( ((m_W
) >= (m_stau_1 + m_snu_tau
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0728*conj(IT_0728)*(IT_0730 + 
      -IT_0733)*IT_0736
) : (0
));
    const ccomplex_t IT_0738 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *U_stau_01;
    const ccomplex_t IT_0739 = -IT_0738;
    const ccomplex_t IT_0740 = -IT_0739;
    const ccomplex_t IT_0741 = pow(m_stau_2, 2);
    const ccomplex_t IT_0742 = 0.5*IT_0011 + 0.5*IT_0730 + (-0.5)*IT_0741;
    const ccomplex_t IT_0743 = cpow(IT_0742, 2);
    const ccomplex_t IT_0744 = IT_0001*IT_0743;
    const ccomplex_t IT_0745 = cpow(IT_0011 + IT_0730 + -IT_0741, 2);
    const ccomplex_t IT_0746 = cabs(-IT_0730 + 0.25*IT_0001*IT_0745);
    const ccomplex_t IT_0747 = cpow(IT_0746, 0.5);
    const ccomplex_t IT_0748 = ( ((m_W
) >= (m_stau_2 + m_snu_tau
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0740*conj(IT_0740)*(IT_0730 + 
      -IT_0744)*IT_0747
) : (0
));
    const ccomplex_t IT_0749 = IT_0077 + IT_0094 + IT_0109 + IT_0125 + IT_0140
       + IT_0155 + IT_0170 + IT_0184 + IT_0198 + IT_0215 + IT_0230 + IT_0245 +
       IT_0256 + IT_0268 + IT_0278 + IT_0288 + IT_0301 + IT_0331 + IT_0360 +
       IT_0389 + IT_0417 + IT_0446 + IT_0474 + IT_0503 + IT_0531 + IT_0544 +
       IT_0556 + IT_0567 + IT_0577 + IT_0587 + IT_0598 + IT_0608 + IT_0617 +
       IT_0626 + IT_0637 + IT_0647 + IT_0656 + IT_0665 + IT_0676 + IT_0686 +
       IT_0696 + IT_0706 + IT_0717 + IT_0725 + IT_0737 + IT_0748;
    return create_ccomplex_return(IT_0749);
}

