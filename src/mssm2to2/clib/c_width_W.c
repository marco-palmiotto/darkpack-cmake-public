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
    const ccomplex_t IT_0006 = (-2)*IT_0005;
    const ccomplex_t IT_0007 = -IT_0005;
    const ccomplex_t IT_0008 = IT_0006 + -IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = pow(m_W, -4);
    const ccomplex_t IT_0011 = pow(m_Z, 2);
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, 2);
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    const ccomplex_t IT_0016 = IT_0010*IT_0015;
    const ccomplex_t IT_0017 = pow(m_Z, -2);
    const ccomplex_t IT_0018 = 0.5*IT_0011;
    const ccomplex_t IT_0019 = cpow(IT_0018, 2);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (1 + 0.5*IT_0016)*(IT_0013 + -IT_0020);
    const ccomplex_t IT_0022 = IT_0005 + -IT_0007;
    const ccomplex_t IT_0023 = IT_0001*IT_0015;
    const ccomplex_t IT_0024 = IT_0001*IT_0017*IT_0019;
    const ccomplex_t IT_0025 = (IT_0013 + -IT_0023)*(1 + 0.5*IT_0024);
    const ccomplex_t IT_0026 = 2*IT_0005;
    const ccomplex_t IT_0027 = cpow(IT_0014, 3);
    const ccomplex_t IT_0028 = IT_0010*IT_0027;
    const ccomplex_t IT_0029 = IT_0014*IT_0024;
    const ccomplex_t IT_0030 = cpow((-0.5)*IT_0011 + IT_0013, 2);
    const ccomplex_t IT_0031 = -IT_0013 + 0.5*IT_0001*IT_0011*((-0.5)*IT_0011 
      + IT_0013) + IT_0001*IT_0030;
    const ccomplex_t IT_0032 = cpow(IT_0031, 2);
    const ccomplex_t IT_0033 = IT_0001*IT_0017*IT_0032;
    const ccomplex_t IT_0034 = IT_0014*IT_0033;
    const ccomplex_t IT_0035 = IT_0017*IT_0031;
    const ccomplex_t IT_0036 = IT_0018*IT_0035;
    const ccomplex_t IT_0037 = IT_0010*IT_0015*IT_0017;
    const ccomplex_t IT_0038 = IT_0031*IT_0037;
    const ccomplex_t IT_0039 = IT_0018*IT_0038;
    const ccomplex_t IT_0040 = IT_0001*IT_0017;
    const ccomplex_t IT_0041 = IT_0018*IT_0040;
    const ccomplex_t IT_0042 = IT_0014*IT_0041;
    const ccomplex_t IT_0043 = IT_0001*IT_0018;
    const ccomplex_t IT_0044 = IT_0014*IT_0043;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = IT_0031 + IT_0045;
    const ccomplex_t IT_0047 = IT_0042*IT_0046;
    const ccomplex_t IT_0048 = -IT_0016 + -IT_0033;
    const ccomplex_t IT_0049 = IT_0010*IT_0017;
    const ccomplex_t IT_0050 = IT_0031*IT_0049;
    const ccomplex_t IT_0051 = IT_0018*IT_0050;
    const ccomplex_t IT_0052 = IT_0014*IT_0051;
    const ccomplex_t IT_0053 = IT_0048 + IT_0052;
    const ccomplex_t IT_0054 = 1 + IT_0053;
    const ccomplex_t IT_0055 = IT_0014*IT_0054;
    const ccomplex_t IT_0056 = (IT_0013 + -IT_0023)*(1 + 0.5*IT_0033);
    const ccomplex_t IT_0057 = -IT_0016 + -IT_0024;
    const ccomplex_t IT_0058 = IT_0052 + IT_0057;
    const ccomplex_t IT_0059 = 1 + IT_0058;
    const ccomplex_t IT_0060 = IT_0014*IT_0059;
    const ccomplex_t IT_0061 = IT_0031*IT_0040;
    const ccomplex_t IT_0062 = IT_0014*IT_0061;
    const ccomplex_t IT_0063 = IT_0001*IT_0031;
    const ccomplex_t IT_0064 = IT_0014*IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0018 + IT_0065;
    const ccomplex_t IT_0067 = IT_0062*IT_0066;
    const ccomplex_t IT_0068 = IT_0017*IT_0032;
    const ccomplex_t IT_0069 = (1 + 0.5*IT_0016)*(IT_0013 + -IT_0068);
    const ccomplex_t IT_0070 = (1 + 0.5*IT_0016)*(IT_0014 + -IT_0036);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = 0.666666666666667*IT_0009*IT_0021 +
       0.666666666666667*IT_0022*(IT_0022*IT_0025 + -IT_0026*(IT_0014 + -IT_0028
       + -IT_0029 + -IT_0034 + IT_0036 + IT_0039) + IT_0008*IT_0047 + -IT_0007
      *IT_0055) + 0.666666666666667*IT_0026*(IT_0026*IT_0056 + -IT_0008*IT_0060 
      + IT_0007*IT_0067) + 0.666666666666667*IT_0007*(IT_0007*IT_0069 + 2
      *IT_0008*IT_0072);
    const ccomplex_t IT_0074 = cabs(-IT_0013 + IT_0001*IT_0030);
    const ccomplex_t IT_0075 = cpow(IT_0074, 0.5);
    const ccomplex_t IT_0076 = ( ((m_W
) >= (m_W + m_Z
)) ? (0.125*IT_0000
      *IT_0001*IT_0073*IT_0075
) : (0
));
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0004;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = cpow(IT_0078, 2);
    const ccomplex_t IT_0080 = pow(m_c, 2);
    const ccomplex_t IT_0081 = pow(m_b, 2);
    const ccomplex_t IT_0082 = cpow(IT_0013 + IT_0080 + -IT_0081, 2);
    const ccomplex_t IT_0083 = -IT_0080 + 0.25*IT_0001*(IT_0013 + -IT_0080 +
       IT_0081)*(IT_0013 + IT_0080 + -IT_0081) + 0.25*IT_0001*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0013 + 0.5*IT_0080 + (-0.5)*IT_0081;
    const ccomplex_t IT_0085 = 0.5*IT_0013 + (-0.5)*IT_0080 + 0.5*IT_0081;
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = IT_0013*IT_0083;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = IT_0086 + IT_0088;
    const ccomplex_t IT_0090 = IT_0001*IT_0089;
    const ccomplex_t IT_0091 = cabs(-IT_0080 + 0.25*IT_0001*IT_0082);
    const ccomplex_t IT_0092 = cpow(IT_0091, 0.5);
    const ccomplex_t IT_0093 = ( ((m_W
) >= (m_b + m_c
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0079*(IT_0083 + IT_0090)*IT_0092
) : (0
));
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*V_cd*e_em
      *IT_0004;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = pow(m_d, 2);
    const ccomplex_t IT_0097 = cpow(IT_0013 + -IT_0080 + IT_0096, 2);
    const ccomplex_t IT_0098 = -IT_0096 + 0.25*IT_0001*(IT_0013 + -IT_0080 +
       IT_0096)*(IT_0013 + IT_0080 + -IT_0096) + 0.25*IT_0001*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0013 + (-0.5)*IT_0080 + 0.5*IT_0096;
    const ccomplex_t IT_0100 = 0.5*IT_0013 + 0.5*IT_0080 + (-0.5)*IT_0096;
    const ccomplex_t IT_0101 = IT_0099*IT_0100;
    const ccomplex_t IT_0102 = IT_0013*IT_0098;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0101 + IT_0103;
    const ccomplex_t IT_0105 = IT_0001*IT_0104;
    const ccomplex_t IT_0106 = cabs(-IT_0096 + 0.25*IT_0001*IT_0097);
    const ccomplex_t IT_0107 = cpow(IT_0106, 0.5);
    const ccomplex_t IT_0108 = ( ((m_W
) >= (m_c + m_d
)) ? (0.5*IT_0000
      *IT_0001*IT_0095*conj(IT_0095)*(IT_0098 + IT_0105)*IT_0107
) : (0
));
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = cpow(IT_0110, 2);
    const ccomplex_t IT_0112 = pow(m_e, 2);
    const ccomplex_t IT_0113 = cpow(IT_0013 + -IT_0112, 2);
    const ccomplex_t IT_0114 = 0.25*IT_0001*(IT_0013 + IT_0112)*(IT_0013 + 
      -IT_0112) + 0.25*IT_0001*IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0013 + (-0.5)*IT_0112;
    const ccomplex_t IT_0116 = 0.5*IT_0013 + 0.5*IT_0112;
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = IT_0013*IT_0114;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = IT_0117 + IT_0119;
    const ccomplex_t IT_0121 = IT_0001*IT_0120;
    const ccomplex_t IT_0122 = cabs(0.25*IT_0001*IT_0113);
    const ccomplex_t IT_0123 = cpow(IT_0122, 0.5);
    const ccomplex_t IT_0124 = ( ((m_W
) >= (m_e
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0111*(IT_0114 + IT_0121)*IT_0123
) : (0
));
    const ccomplex_t IT_0125 = cos(alpha);
    const ccomplex_t IT_0126 = sin(beta);
    const ccomplex_t IT_0127 = cos(beta);
    const ccomplex_t IT_0128 = sin(alpha);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_W*e_em*IT_0004*(IT_0125
      *IT_0126 + -IT_0127*IT_0128);
    const ccomplex_t IT_0130 = cpow(IT_0129, 2);
    const ccomplex_t IT_0131 = pow(m_h, 2);
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = IT_0013 + IT_0132;
    const ccomplex_t IT_0134 = cpow(IT_0133, 2);
    const ccomplex_t IT_0135 = IT_0010*IT_0134;
    const ccomplex_t IT_0136 = cpow(IT_0013 + (-0.5)*IT_0131, 2);
    const ccomplex_t IT_0137 = cabs(-IT_0013 + IT_0001*IT_0136);
    const ccomplex_t IT_0138 = cpow(IT_0137, 0.5);
    const ccomplex_t IT_0139 = ( ((m_W
) >= (m_W + m_h
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0130*(1 + 0.5*IT_0135)*IT_0138
) :
       (0
));
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0004;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = pow(m_s, 2);
    const ccomplex_t IT_0143 = cpow(IT_0013 + -IT_0080 + IT_0142, 2);
    const ccomplex_t IT_0144 = -IT_0142 + 0.25*IT_0001*(IT_0013 + -IT_0080 +
       IT_0142)*(IT_0013 + IT_0080 + -IT_0142) + 0.25*IT_0001*IT_0143;
    const ccomplex_t IT_0145 = 0.5*IT_0013 + (-0.5)*IT_0080 + 0.5*IT_0142;
    const ccomplex_t IT_0146 = 0.5*IT_0013 + 0.5*IT_0080 + (-0.5)*IT_0142;
    const ccomplex_t IT_0147 = IT_0145*IT_0146;
    const ccomplex_t IT_0148 = IT_0013*IT_0144;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = IT_0147 + IT_0149;
    const ccomplex_t IT_0151 = IT_0001*IT_0150;
    const ccomplex_t IT_0152 = cabs(-IT_0142 + 0.25*IT_0001*IT_0143);
    const ccomplex_t IT_0153 = cpow(IT_0152, 0.5);
    const ccomplex_t IT_0154 = ( ((m_W
) >= (m_c + m_s
)) ? (0.5*IT_0000
      *IT_0001*IT_0141*conj(IT_0141)*(IT_0144 + IT_0151)*IT_0153
) : (0
));
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004;
    const ccomplex_t IT_0156 = 0.5*IT_0155;
    const ccomplex_t IT_0157 = pow(m_t, 2);
    const ccomplex_t IT_0158 = cpow(IT_0013 + -IT_0081 + IT_0157, 2);
    const ccomplex_t IT_0159 = -IT_0157 + 0.25*IT_0001*(IT_0013 + -IT_0081 +
       IT_0157)*(IT_0013 + IT_0081 + -IT_0157) + 0.25*IT_0001*IT_0158;
    const ccomplex_t IT_0160 = 0.5*IT_0013 + (-0.5)*IT_0081 + 0.5*IT_0157;
    const ccomplex_t IT_0161 = 0.5*IT_0013 + 0.5*IT_0081 + (-0.5)*IT_0157;
    const ccomplex_t IT_0162 = IT_0160*IT_0161;
    const ccomplex_t IT_0163 = IT_0013*IT_0159;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = IT_0162 + IT_0164;
    const ccomplex_t IT_0166 = IT_0001*IT_0165;
    const ccomplex_t IT_0167 = cabs(-IT_0157 + 0.25*IT_0001*IT_0158);
    const ccomplex_t IT_0168 = cpow(IT_0167, 0.5);
    const ccomplex_t IT_0169 = ( ((m_W
) >= (m_b + m_t
)) ? (0.5*IT_0000
      *IT_0001*IT_0156*conj(IT_0156)*(IT_0159 + IT_0166)*IT_0168
) : (0
));
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0004;
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = cpow(IT_0013 + IT_0096 + -IT_0157, 2);
    const ccomplex_t IT_0173 = -IT_0096 + 0.25*IT_0001*(IT_0013 + -IT_0096 +
       IT_0157)*(IT_0013 + IT_0096 + -IT_0157) + 0.25*IT_0001*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0013 + 0.5*IT_0096 + (-0.5)*IT_0157;
    const ccomplex_t IT_0175 = 0.5*IT_0013 + (-0.5)*IT_0096 + 0.5*IT_0157;
    const ccomplex_t IT_0176 = IT_0174*IT_0175;
    const ccomplex_t IT_0177 = IT_0013*IT_0173;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = IT_0001*IT_0179;
    const ccomplex_t IT_0181 = cabs(-IT_0096 + 0.25*IT_0001*IT_0172);
    const ccomplex_t IT_0182 = cpow(IT_0181, 0.5);
    const ccomplex_t IT_0183 = ( ((m_W
) >= (m_d + m_t
)) ? (0.5*IT_0000
      *IT_0001*IT_0171*conj(IT_0171)*(IT_0173 + IT_0180)*IT_0182
) : (0
));
    const ccomplex_t IT_0184 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0004;
    const ccomplex_t IT_0185 = 0.5*IT_0184;
    const ccomplex_t IT_0186 = cpow(IT_0013 + IT_0142 + -IT_0157, 2);
    const ccomplex_t IT_0187 = -IT_0142 + 0.25*IT_0001*(IT_0013 + -IT_0142 +
       IT_0157)*(IT_0013 + IT_0142 + -IT_0157) + 0.25*IT_0001*IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0013 + 0.5*IT_0142 + (-0.5)*IT_0157;
    const ccomplex_t IT_0189 = 0.5*IT_0013 + (-0.5)*IT_0142 + 0.5*IT_0157;
    const ccomplex_t IT_0190 = IT_0188*IT_0189;
    const ccomplex_t IT_0191 = IT_0013*IT_0187;
    const ccomplex_t IT_0192 = (-0.5)*IT_0191;
    const ccomplex_t IT_0193 = IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = IT_0001*IT_0193;
    const ccomplex_t IT_0195 = cabs(-IT_0142 + 0.25*IT_0001*IT_0186);
    const ccomplex_t IT_0196 = cpow(IT_0195, 0.5);
    const ccomplex_t IT_0197 = ( ((m_W
) >= (m_s + m_t
)) ? (0.5*IT_0000
      *IT_0001*IT_0185*conj(IT_0185)*(IT_0187 + IT_0194)*IT_0196
) : (0
));
    const ccomplex_t IT_0198 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0199 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *IT_0198*V_ub_mod;
    const ccomplex_t IT_0200 = 0.5*IT_0199;
    const ccomplex_t IT_0201 = cpow(IT_0200, 2);
    const ccomplex_t IT_0202 = pow(m_u, 2);
    const ccomplex_t IT_0203 = cpow(IT_0013 + -IT_0081 + IT_0202, 2);
    const ccomplex_t IT_0204 = -IT_0202 + 0.25*IT_0001*(IT_0013 + -IT_0081 +
       IT_0202)*(IT_0013 + IT_0081 + -IT_0202) + 0.25*IT_0001*IT_0203;
    const ccomplex_t IT_0205 = 0.5*IT_0013 + (-0.5)*IT_0081 + 0.5*IT_0202;
    const ccomplex_t IT_0206 = 0.5*IT_0013 + 0.5*IT_0081 + (-0.5)*IT_0202;
    const ccomplex_t IT_0207 = IT_0205*IT_0206;
    const ccomplex_t IT_0208 = IT_0013*IT_0204;
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = IT_0207 + IT_0209;
    const ccomplex_t IT_0211 = IT_0001*IT_0210;
    const ccomplex_t IT_0212 = cabs(-IT_0202 + 0.25*IT_0001*IT_0203);
    const ccomplex_t IT_0213 = cpow(IT_0212, 0.5);
    const ccomplex_t IT_0214 = ( ((m_W
) >= (m_b + m_u
)) ? (0.5*IT_0000
      *IT_0001*IT_0201*(IT_0204 + IT_0211)*IT_0213
) : (0
));
    const ccomplex_t IT_0215 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0004;
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = cpow(IT_0216, 2);
    const ccomplex_t IT_0218 = cpow(IT_0013 + -IT_0096 + IT_0202, 2);
    const ccomplex_t IT_0219 = -IT_0202 + 0.25*IT_0001*(IT_0013 + -IT_0096 +
       IT_0202)*(IT_0013 + IT_0096 + -IT_0202) + 0.25*IT_0001*IT_0218;
    const ccomplex_t IT_0220 = 0.5*IT_0013 + (-0.5)*IT_0096 + 0.5*IT_0202;
    const ccomplex_t IT_0221 = 0.5*IT_0013 + 0.5*IT_0096 + (-0.5)*IT_0202;
    const ccomplex_t IT_0222 = IT_0220*IT_0221;
    const ccomplex_t IT_0223 = IT_0013*IT_0219;
    const ccomplex_t IT_0224 = (-0.5)*IT_0223;
    const ccomplex_t IT_0225 = IT_0222 + IT_0224;
    const ccomplex_t IT_0226 = IT_0001*IT_0225;
    const ccomplex_t IT_0227 = cabs(-IT_0202 + 0.25*IT_0001*IT_0218);
    const ccomplex_t IT_0228 = cpow(IT_0227, 0.5);
    const ccomplex_t IT_0229 = ( ((m_W
) >= (m_d + m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0217*(IT_0219 + IT_0226)*IT_0228
) : (0
));
    const ccomplex_t IT_0230 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0004;
    const ccomplex_t IT_0231 = 0.5*IT_0230;
    const ccomplex_t IT_0232 = cpow(IT_0231, 2);
    const ccomplex_t IT_0233 = cpow(IT_0013 + -IT_0142 + IT_0202, 2);
    const ccomplex_t IT_0234 = -IT_0202 + 0.25*IT_0001*(IT_0013 + -IT_0142 +
       IT_0202)*(IT_0013 + IT_0142 + -IT_0202) + 0.25*IT_0001*IT_0233;
    const ccomplex_t IT_0235 = 0.5*IT_0013 + (-0.5)*IT_0142 + 0.5*IT_0202;
    const ccomplex_t IT_0236 = 0.5*IT_0013 + 0.5*IT_0142 + (-0.5)*IT_0202;
    const ccomplex_t IT_0237 = IT_0235*IT_0236;
    const ccomplex_t IT_0238 = IT_0013*IT_0234;
    const ccomplex_t IT_0239 = (-0.5)*IT_0238;
    const ccomplex_t IT_0240 = IT_0237 + IT_0239;
    const ccomplex_t IT_0241 = IT_0001*IT_0240;
    const ccomplex_t IT_0242 = cabs(-IT_0202 + 0.25*IT_0001*IT_0233);
    const ccomplex_t IT_0243 = cpow(IT_0242, 0.5);
    const ccomplex_t IT_0244 = ( ((m_W
) >= (m_s + m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0232*(IT_0234 + IT_0241)*IT_0243
) : (0
));
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*m_W*e_em*IT_0004*(IT_0125
      *IT_0127 + IT_0126*IT_0128);
    const ccomplex_t IT_0246 = cpow(IT_0245, 2);
    const ccomplex_t IT_0247 = pow(m_H0, 2);
    const ccomplex_t IT_0248 = (-0.5)*IT_0247;
    const ccomplex_t IT_0249 = IT_0013 + IT_0248;
    const ccomplex_t IT_0250 = cpow(IT_0249, 2);
    const ccomplex_t IT_0251 = IT_0010*IT_0250;
    const ccomplex_t IT_0252 = cpow(IT_0013 + (-0.5)*IT_0247, 2);
    const ccomplex_t IT_0253 = cabs(-IT_0013 + IT_0001*IT_0252);
    const ccomplex_t IT_0254 = cpow(IT_0253, 0.5);
    const ccomplex_t IT_0255 = ( ((m_W
) >= (m_W + m_H0
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0246*(1 + 0.5*IT_0251)*IT_0254
) :
       (0
));
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*e_em*IT_0004*(IT_0125
      *IT_0127 + IT_0126*IT_0128);
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = cpow(IT_0258, 2);
    const ccomplex_t IT_0260 = pow(m_Hp, 2);
    const ccomplex_t IT_0261 = 0.5*IT_0013 + (-0.5)*IT_0131 + 0.5*IT_0260;
    const ccomplex_t IT_0262 = cpow(IT_0261, 2);
    const ccomplex_t IT_0263 = IT_0001*IT_0262;
    const ccomplex_t IT_0264 = cpow(IT_0013 + -IT_0131 + IT_0260, 2);
    const ccomplex_t IT_0265 = cabs(-IT_0260 + 0.25*IT_0001*IT_0264);
    const ccomplex_t IT_0266 = cpow(IT_0265, 0.5);
    const ccomplex_t IT_0267 = ( ((m_W
) >= (m_h + m_Hp
)) ? 
      (0.0416666666666667*IT_0000*IT_0001*IT_0259*(IT_0260 + -IT_0263)*IT_0266
)
       : (0
));
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*e_em*IT_0004*(IT_0125
      *IT_0126 + -IT_0127*IT_0128);
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = cpow(IT_0269, 2);
    const ccomplex_t IT_0271 = 0.5*IT_0013 + (-0.5)*IT_0247 + 0.5*IT_0260;
    const ccomplex_t IT_0272 = cpow(IT_0271, 2);
    const ccomplex_t IT_0273 = IT_0001*IT_0272;
    const ccomplex_t IT_0274 = cpow(IT_0013 + -IT_0247 + IT_0260, 2);
    const ccomplex_t IT_0275 = cabs(-IT_0260 + 0.25*IT_0001*IT_0274);
    const ccomplex_t IT_0276 = cpow(IT_0275, 0.5);
    const ccomplex_t IT_0277 = ( ((m_W
) >= (m_H0 + m_Hp
)) ? 
      (0.0416666666666667*IT_0000*IT_0001*IT_0270*(IT_0260 + -IT_0273)*IT_0276
)
       : (0
));
    const ccomplex_t IT_0278 = e_em*IT_0004;
    const ccomplex_t IT_0279 = cpow(IT_0278, 2);
    const ccomplex_t IT_0280 = pow(m_A0, 2);
    const ccomplex_t IT_0281 = 0.5*IT_0013 + (-0.5)*IT_0260 + 0.5*IT_0280;
    const ccomplex_t IT_0282 = cpow(IT_0281, 2);
    const ccomplex_t IT_0283 = IT_0001*IT_0282;
    const ccomplex_t IT_0284 = cpow(IT_0013 + -IT_0260 + IT_0280, 2);
    const ccomplex_t IT_0285 = cabs(-IT_0280 + 0.25*IT_0001*IT_0284);
    const ccomplex_t IT_0286 = cpow(IT_0285, 0.5);
    const ccomplex_t IT_0287 = ( ((m_W
) >= (m_A0 + m_Hp
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0279*(IT_0280 + -IT_0283)*IT_0286

      ) : (0
));
    const ccomplex_t IT_0288 = pow(m_mu, 2);
    const ccomplex_t IT_0289 = cpow(IT_0013 + -IT_0288, 2);
    const ccomplex_t IT_0290 = 0.25*IT_0001*(IT_0013 + IT_0288)*(IT_0013 + 
      -IT_0288) + 0.25*IT_0001*IT_0289;
    const ccomplex_t IT_0291 = 0.5*IT_0013 + (-0.5)*IT_0288;
    const ccomplex_t IT_0292 = 0.5*IT_0013 + 0.5*IT_0288;
    const ccomplex_t IT_0293 = IT_0291*IT_0292;
    const ccomplex_t IT_0294 = IT_0013*IT_0290;
    const ccomplex_t IT_0295 = (-0.5)*IT_0294;
    const ccomplex_t IT_0296 = IT_0293 + IT_0295;
    const ccomplex_t IT_0297 = IT_0001*IT_0296;
    const ccomplex_t IT_0298 = cabs(0.25*IT_0001*IT_0289);
    const ccomplex_t IT_0299 = cpow(IT_0298, 0.5);
    const ccomplex_t IT_0300 = ( ((m_W
) >= (m_mu
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0111*(IT_0290 + IT_0297)*IT_0299
) : (0
));
    const ccomplex_t IT_0301 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0302 = IT_0004*IT_0301;
    const ccomplex_t IT_0303 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0304 = IT_0004*IT_0303;
    const ccomplex_t IT_0305 = 1.4142135623731*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*(IT_0302 + 0.5*IT_0305);
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0309 = IT_0004*IT_0308;
    const ccomplex_t IT_0310 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0311 = IT_0004*IT_0310;
    const ccomplex_t IT_0312 = 1.4142135623731*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*(IT_0309 + (-0.5)*IT_0312);
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = m_C_1*m_N_1;
    const ccomplex_t IT_0316 = pow(m_N_1, 2);
    const ccomplex_t IT_0317 = pow(m_C_1, 2);
    const ccomplex_t IT_0318 = cpow(IT_0013 + IT_0316 + -IT_0317, 2);
    const ccomplex_t IT_0319 = -IT_0316 + 0.25*IT_0001*(IT_0013 + -IT_0316 +
       IT_0317)*(IT_0013 + IT_0316 + -IT_0317) + 0.25*IT_0001*IT_0318;
    const ccomplex_t IT_0320 = 0.5*IT_0013 + 0.5*IT_0316 + (-0.5)*IT_0317;
    const ccomplex_t IT_0321 = 0.5*IT_0013 + (-0.5)*IT_0316 + 0.5*IT_0317;
    const ccomplex_t IT_0322 = IT_0320*IT_0321;
    const ccomplex_t IT_0323 = IT_0013*IT_0319;
    const ccomplex_t IT_0324 = (-0.5)*IT_0323;
    const ccomplex_t IT_0325 = IT_0322 + IT_0324;
    const ccomplex_t IT_0326 = IT_0001*IT_0325;
    const ccomplex_t IT_0327 = IT_0319 + IT_0326;
    const ccomplex_t IT_0328 = cabs(-IT_0316 + 0.25*IT_0001*IT_0318);
    const ccomplex_t IT_0329 = cpow(IT_0328, 0.5);
    const ccomplex_t IT_0330 = ( ((m_W
) >= (m_C_1 + m_N_1
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0307)*(IT_0314*IT_0315 + 0.666666666666667*IT_0307
      *IT_0327) + conj(IT_0314)*(IT_0307*IT_0315 + 0.666666666666667*IT_0314
      *IT_0327))*IT_0329
) : (0
));
    const ccomplex_t IT_0331 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0332 = IT_0004*IT_0331;
    const ccomplex_t IT_0333 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0334 = IT_0004*IT_0333;
    const ccomplex_t IT_0335 = 1.4142135623731*IT_0334;
    const ccomplex_t IT_0336 = (0 + _Complex_I*1)*(IT_0332 + 0.5*IT_0335);
    const ccomplex_t IT_0337 = -IT_0336;
    const ccomplex_t IT_0338 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0339 = IT_0004*IT_0338;
    const ccomplex_t IT_0340 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0341 = IT_0004*IT_0340;
    const ccomplex_t IT_0342 = 1.4142135623731*IT_0341;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*(IT_0339 + (-0.5)*IT_0342);
    const ccomplex_t IT_0344 = -IT_0343;
    const ccomplex_t IT_0345 = m_C_2*m_N_1;
    const ccomplex_t IT_0346 = pow(m_C_2, 2);
    const ccomplex_t IT_0347 = cpow(IT_0013 + IT_0316 + -IT_0346, 2);
    const ccomplex_t IT_0348 = -IT_0316 + 0.25*IT_0001*(IT_0013 + -IT_0316 +
       IT_0346)*(IT_0013 + IT_0316 + -IT_0346) + 0.25*IT_0001*IT_0347;
    const ccomplex_t IT_0349 = 0.5*IT_0013 + 0.5*IT_0316 + (-0.5)*IT_0346;
    const ccomplex_t IT_0350 = 0.5*IT_0013 + (-0.5)*IT_0316 + 0.5*IT_0346;
    const ccomplex_t IT_0351 = IT_0349*IT_0350;
    const ccomplex_t IT_0352 = IT_0013*IT_0348;
    const ccomplex_t IT_0353 = (-0.5)*IT_0352;
    const ccomplex_t IT_0354 = IT_0351 + IT_0353;
    const ccomplex_t IT_0355 = IT_0001*IT_0354;
    const ccomplex_t IT_0356 = IT_0348 + IT_0355;
    const ccomplex_t IT_0357 = cabs(-IT_0316 + 0.25*IT_0001*IT_0347);
    const ccomplex_t IT_0358 = cpow(IT_0357, 0.5);
    const ccomplex_t IT_0359 = ( ((m_W
) >= (m_C_2 + m_N_1
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0337)*(IT_0344*IT_0345 + 0.666666666666667*IT_0337
      *IT_0356) + conj(IT_0344)*(IT_0337*IT_0345 + 0.666666666666667*IT_0344
      *IT_0356))*IT_0358
) : (0
));
    const ccomplex_t IT_0360 = conj(N_W2)*e_em*U_Wm1;
    const ccomplex_t IT_0361 = IT_0004*IT_0360;
    const ccomplex_t IT_0362 = conj(N_d2)*U_d1*e_em;
    const ccomplex_t IT_0363 = IT_0004*IT_0362;
    const ccomplex_t IT_0364 = 1.4142135623731*IT_0363;
    const ccomplex_t IT_0365 = (0 + _Complex_I*1)*(IT_0361 + 0.5*IT_0364);
    const ccomplex_t IT_0366 = -IT_0365;
    const ccomplex_t IT_0367 = N_W2*e_em*conj(V_Wp1);
    const ccomplex_t IT_0368 = IT_0004*IT_0367;
    const ccomplex_t IT_0369 = N_u2*conj(V_u1)*e_em;
    const ccomplex_t IT_0370 = IT_0004*IT_0369;
    const ccomplex_t IT_0371 = 1.4142135623731*IT_0370;
    const ccomplex_t IT_0372 = (0 + _Complex_I*1)*(IT_0368 + (-0.5)*IT_0371);
    const ccomplex_t IT_0373 = -IT_0372;
    const ccomplex_t IT_0374 = m_C_1*m_N_2;
    const ccomplex_t IT_0375 = pow(m_N_2, 2);
    const ccomplex_t IT_0376 = cpow(IT_0013 + -IT_0317 + IT_0375, 2);
    const ccomplex_t IT_0377 = -IT_0375 + 0.25*IT_0001*(IT_0013 + -IT_0317 +
       IT_0375)*(IT_0013 + IT_0317 + -IT_0375) + 0.25*IT_0001*IT_0376;
    const ccomplex_t IT_0378 = 0.5*IT_0013 + (-0.5)*IT_0317 + 0.5*IT_0375;
    const ccomplex_t IT_0379 = 0.5*IT_0013 + 0.5*IT_0317 + (-0.5)*IT_0375;
    const ccomplex_t IT_0380 = IT_0378*IT_0379;
    const ccomplex_t IT_0381 = IT_0013*IT_0377;
    const ccomplex_t IT_0382 = (-0.5)*IT_0381;
    const ccomplex_t IT_0383 = IT_0380 + IT_0382;
    const ccomplex_t IT_0384 = IT_0001*IT_0383;
    const ccomplex_t IT_0385 = IT_0377 + IT_0384;
    const ccomplex_t IT_0386 = cabs(-IT_0375 + 0.25*IT_0001*IT_0376);
    const ccomplex_t IT_0387 = cpow(IT_0386, 0.5);
    const ccomplex_t IT_0388 = ( ((m_W
) >= (m_C_1 + m_N_2
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0366)*(IT_0373*IT_0374 + 0.666666666666667*IT_0366
      *IT_0385) + conj(IT_0373)*(IT_0366*IT_0374 + 0.666666666666667*IT_0373
      *IT_0385))*IT_0387
) : (0
));
    const ccomplex_t IT_0389 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0390 = IT_0004*IT_0389;
    const ccomplex_t IT_0391 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0392 = IT_0004*IT_0391;
    const ccomplex_t IT_0393 = 1.4142135623731*IT_0392;
    const ccomplex_t IT_0394 = (0 + _Complex_I*1)*(IT_0390 + 0.5*IT_0393);
    const ccomplex_t IT_0395 = -IT_0394;
    const ccomplex_t IT_0396 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0397 = IT_0004*IT_0396;
    const ccomplex_t IT_0398 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0399 = IT_0004*IT_0398;
    const ccomplex_t IT_0400 = 1.4142135623731*IT_0399;
    const ccomplex_t IT_0401 = (0 + _Complex_I*1)*(IT_0397 + (-0.5)*IT_0400);
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = m_C_2*m_N_2;
    const ccomplex_t IT_0404 = cpow(IT_0013 + -IT_0346 + IT_0375, 2);
    const ccomplex_t IT_0405 = -IT_0375 + 0.25*IT_0001*(IT_0013 + -IT_0346 +
       IT_0375)*(IT_0013 + IT_0346 + -IT_0375) + 0.25*IT_0001*IT_0404;
    const ccomplex_t IT_0406 = 0.5*IT_0013 + (-0.5)*IT_0346 + 0.5*IT_0375;
    const ccomplex_t IT_0407 = 0.5*IT_0013 + 0.5*IT_0346 + (-0.5)*IT_0375;
    const ccomplex_t IT_0408 = IT_0406*IT_0407;
    const ccomplex_t IT_0409 = IT_0013*IT_0405;
    const ccomplex_t IT_0410 = (-0.5)*IT_0409;
    const ccomplex_t IT_0411 = IT_0408 + IT_0410;
    const ccomplex_t IT_0412 = IT_0001*IT_0411;
    const ccomplex_t IT_0413 = IT_0405 + IT_0412;
    const ccomplex_t IT_0414 = cabs(-IT_0375 + 0.25*IT_0001*IT_0404);
    const ccomplex_t IT_0415 = cpow(IT_0414, 0.5);
    const ccomplex_t IT_0416 = ( ((m_W
) >= (m_C_2 + m_N_2
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0395)*(IT_0402*IT_0403 + 0.666666666666667*IT_0395
      *IT_0413) + conj(IT_0402)*(IT_0395*IT_0403 + 0.666666666666667*IT_0402
      *IT_0413))*IT_0415
) : (0
));
    const ccomplex_t IT_0417 = conj(N_W3)*e_em*U_Wm1;
    const ccomplex_t IT_0418 = IT_0004*IT_0417;
    const ccomplex_t IT_0419 = conj(N_d3)*U_d1*e_em;
    const ccomplex_t IT_0420 = IT_0004*IT_0419;
    const ccomplex_t IT_0421 = 1.4142135623731*IT_0420;
    const ccomplex_t IT_0422 = (0 + _Complex_I*1)*(IT_0418 + 0.5*IT_0421);
    const ccomplex_t IT_0423 = -IT_0422;
    const ccomplex_t IT_0424 = N_W3*e_em*conj(V_Wp1);
    const ccomplex_t IT_0425 = IT_0004*IT_0424;
    const ccomplex_t IT_0426 = N_u3*conj(V_u1)*e_em;
    const ccomplex_t IT_0427 = IT_0004*IT_0426;
    const ccomplex_t IT_0428 = 1.4142135623731*IT_0427;
    const ccomplex_t IT_0429 = (0 + _Complex_I*1)*(IT_0425 + (-0.5)*IT_0428);
    const ccomplex_t IT_0430 = -IT_0429;
    const ccomplex_t IT_0431 = m_C_1*m_N_3;
    const ccomplex_t IT_0432 = pow(m_N_3, 2);
    const ccomplex_t IT_0433 = cpow(IT_0013 + -IT_0317 + IT_0432, 2);
    const ccomplex_t IT_0434 = -IT_0432 + 0.25*IT_0001*(IT_0013 + -IT_0317 +
       IT_0432)*(IT_0013 + IT_0317 + -IT_0432) + 0.25*IT_0001*IT_0433;
    const ccomplex_t IT_0435 = 0.5*IT_0013 + (-0.5)*IT_0317 + 0.5*IT_0432;
    const ccomplex_t IT_0436 = 0.5*IT_0013 + 0.5*IT_0317 + (-0.5)*IT_0432;
    const ccomplex_t IT_0437 = IT_0435*IT_0436;
    const ccomplex_t IT_0438 = IT_0013*IT_0434;
    const ccomplex_t IT_0439 = (-0.5)*IT_0438;
    const ccomplex_t IT_0440 = IT_0437 + IT_0439;
    const ccomplex_t IT_0441 = IT_0001*IT_0440;
    const ccomplex_t IT_0442 = IT_0434 + IT_0441;
    const ccomplex_t IT_0443 = cabs(-IT_0432 + 0.25*IT_0001*IT_0433);
    const ccomplex_t IT_0444 = cpow(IT_0443, 0.5);
    const ccomplex_t IT_0445 = ( ((m_W
) >= (m_C_1 + m_N_3
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0423)*(IT_0430*IT_0431 + 0.666666666666667*IT_0423
      *IT_0442) + conj(IT_0430)*(IT_0423*IT_0431 + 0.666666666666667*IT_0430
      *IT_0442))*IT_0444
) : (0
));
    const ccomplex_t IT_0446 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0447 = IT_0004*IT_0446;
    const ccomplex_t IT_0448 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0449 = IT_0004*IT_0448;
    const ccomplex_t IT_0450 = 1.4142135623731*IT_0449;
    const ccomplex_t IT_0451 = (0 + _Complex_I*1)*(IT_0447 + 0.5*IT_0450);
    const ccomplex_t IT_0452 = -IT_0451;
    const ccomplex_t IT_0453 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0454 = IT_0004*IT_0453;
    const ccomplex_t IT_0455 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0456 = IT_0004*IT_0455;
    const ccomplex_t IT_0457 = 1.4142135623731*IT_0456;
    const ccomplex_t IT_0458 = (0 + _Complex_I*1)*(IT_0454 + (-0.5)*IT_0457);
    const ccomplex_t IT_0459 = -IT_0458;
    const ccomplex_t IT_0460 = m_C_2*m_N_3;
    const ccomplex_t IT_0461 = cpow(IT_0013 + -IT_0346 + IT_0432, 2);
    const ccomplex_t IT_0462 = -IT_0432 + 0.25*IT_0001*(IT_0013 + -IT_0346 +
       IT_0432)*(IT_0013 + IT_0346 + -IT_0432) + 0.25*IT_0001*IT_0461;
    const ccomplex_t IT_0463 = 0.5*IT_0013 + (-0.5)*IT_0346 + 0.5*IT_0432;
    const ccomplex_t IT_0464 = 0.5*IT_0013 + 0.5*IT_0346 + (-0.5)*IT_0432;
    const ccomplex_t IT_0465 = IT_0463*IT_0464;
    const ccomplex_t IT_0466 = IT_0013*IT_0462;
    const ccomplex_t IT_0467 = (-0.5)*IT_0466;
    const ccomplex_t IT_0468 = IT_0465 + IT_0467;
    const ccomplex_t IT_0469 = IT_0001*IT_0468;
    const ccomplex_t IT_0470 = IT_0462 + IT_0469;
    const ccomplex_t IT_0471 = cabs(-IT_0432 + 0.25*IT_0001*IT_0461);
    const ccomplex_t IT_0472 = cpow(IT_0471, 0.5);
    const ccomplex_t IT_0473 = ( ((m_W
) >= (m_C_2 + m_N_3
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0452)*(IT_0459*IT_0460 + 0.666666666666667*IT_0452
      *IT_0470) + conj(IT_0459)*(IT_0452*IT_0460 + 0.666666666666667*IT_0459
      *IT_0470))*IT_0472
) : (0
));
    const ccomplex_t IT_0474 = conj(N_W4)*e_em*U_Wm1;
    const ccomplex_t IT_0475 = IT_0004*IT_0474;
    const ccomplex_t IT_0476 = conj(N_d4)*U_d1*e_em;
    const ccomplex_t IT_0477 = IT_0004*IT_0476;
    const ccomplex_t IT_0478 = 1.4142135623731*IT_0477;
    const ccomplex_t IT_0479 = (0 + _Complex_I*1)*(IT_0475 + 0.5*IT_0478);
    const ccomplex_t IT_0480 = -IT_0479;
    const ccomplex_t IT_0481 = N_W4*e_em*conj(V_Wp1);
    const ccomplex_t IT_0482 = IT_0004*IT_0481;
    const ccomplex_t IT_0483 = N_u4*conj(V_u1)*e_em;
    const ccomplex_t IT_0484 = IT_0004*IT_0483;
    const ccomplex_t IT_0485 = 1.4142135623731*IT_0484;
    const ccomplex_t IT_0486 = (0 + _Complex_I*1)*(IT_0482 + (-0.5)*IT_0485);
    const ccomplex_t IT_0487 = -IT_0486;
    const ccomplex_t IT_0488 = m_C_1*m_N_4;
    const ccomplex_t IT_0489 = pow(m_N_4, 2);
    const ccomplex_t IT_0490 = cpow(IT_0013 + -IT_0317 + IT_0489, 2);
    const ccomplex_t IT_0491 = -IT_0489 + 0.25*IT_0001*(IT_0013 + -IT_0317 +
       IT_0489)*(IT_0013 + IT_0317 + -IT_0489) + 0.25*IT_0001*IT_0490;
    const ccomplex_t IT_0492 = 0.5*IT_0013 + (-0.5)*IT_0317 + 0.5*IT_0489;
    const ccomplex_t IT_0493 = 0.5*IT_0013 + 0.5*IT_0317 + (-0.5)*IT_0489;
    const ccomplex_t IT_0494 = IT_0492*IT_0493;
    const ccomplex_t IT_0495 = IT_0013*IT_0491;
    const ccomplex_t IT_0496 = (-0.5)*IT_0495;
    const ccomplex_t IT_0497 = IT_0494 + IT_0496;
    const ccomplex_t IT_0498 = IT_0001*IT_0497;
    const ccomplex_t IT_0499 = IT_0491 + IT_0498;
    const ccomplex_t IT_0500 = cabs(-IT_0489 + 0.25*IT_0001*IT_0490);
    const ccomplex_t IT_0501 = cpow(IT_0500, 0.5);
    const ccomplex_t IT_0502 = ( ((m_W
) >= (m_C_1 + m_N_4
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0480)*(IT_0487*IT_0488 + 0.666666666666667*IT_0480
      *IT_0499) + conj(IT_0487)*(IT_0480*IT_0488 + 0.666666666666667*IT_0487
      *IT_0499))*IT_0501
) : (0
));
    const ccomplex_t IT_0503 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0504 = IT_0004*IT_0503;
    const ccomplex_t IT_0505 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0506 = IT_0004*IT_0505;
    const ccomplex_t IT_0507 = 1.4142135623731*IT_0506;
    const ccomplex_t IT_0508 = (0 + _Complex_I*1)*(IT_0504 + 0.5*IT_0507);
    const ccomplex_t IT_0509 = -IT_0508;
    const ccomplex_t IT_0510 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0511 = IT_0004*IT_0510;
    const ccomplex_t IT_0512 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0513 = IT_0004*IT_0512;
    const ccomplex_t IT_0514 = 1.4142135623731*IT_0513;
    const ccomplex_t IT_0515 = (0 + _Complex_I*1)*(IT_0511 + (-0.5)*IT_0514);
    const ccomplex_t IT_0516 = -IT_0515;
    const ccomplex_t IT_0517 = m_C_2*m_N_4;
    const ccomplex_t IT_0518 = cpow(IT_0013 + -IT_0346 + IT_0489, 2);
    const ccomplex_t IT_0519 = -IT_0489 + 0.25*IT_0001*(IT_0013 + -IT_0346 +
       IT_0489)*(IT_0013 + IT_0346 + -IT_0489) + 0.25*IT_0001*IT_0518;
    const ccomplex_t IT_0520 = 0.5*IT_0013 + (-0.5)*IT_0346 + 0.5*IT_0489;
    const ccomplex_t IT_0521 = 0.5*IT_0013 + 0.5*IT_0346 + (-0.5)*IT_0489;
    const ccomplex_t IT_0522 = IT_0520*IT_0521;
    const ccomplex_t IT_0523 = IT_0013*IT_0519;
    const ccomplex_t IT_0524 = (-0.5)*IT_0523;
    const ccomplex_t IT_0525 = IT_0522 + IT_0524;
    const ccomplex_t IT_0526 = IT_0001*IT_0525;
    const ccomplex_t IT_0527 = IT_0519 + IT_0526;
    const ccomplex_t IT_0528 = cabs(-IT_0489 + 0.25*IT_0001*IT_0518);
    const ccomplex_t IT_0529 = cpow(IT_0528, 0.5);
    const ccomplex_t IT_0530 = ( ((m_W
) >= (m_C_2 + m_N_4
)) ? (0.25*IT_0000
      *IT_0001*(conj(IT_0509)*(IT_0516*IT_0517 + 0.666666666666667*IT_0509
      *IT_0527) + conj(IT_0516)*(IT_0509*IT_0517 + 0.666666666666667*IT_0516
      *IT_0527))*IT_0529
) : (0
));
    const ccomplex_t IT_0531 = pow(m_tau, 2);
    const ccomplex_t IT_0532 = cpow(IT_0013 + -IT_0531, 2);
    const ccomplex_t IT_0533 = 0.25*IT_0001*(IT_0013 + IT_0531)*(IT_0013 + 
      -IT_0531) + 0.25*IT_0001*IT_0532;
    const ccomplex_t IT_0534 = 0.5*IT_0013 + (-0.5)*IT_0531;
    const ccomplex_t IT_0535 = 0.5*IT_0013 + 0.5*IT_0531;
    const ccomplex_t IT_0536 = IT_0534*IT_0535;
    const ccomplex_t IT_0537 = IT_0013*IT_0533;
    const ccomplex_t IT_0538 = (-0.5)*IT_0537;
    const ccomplex_t IT_0539 = IT_0536 + IT_0538;
    const ccomplex_t IT_0540 = IT_0001*IT_0539;
    const ccomplex_t IT_0541 = cabs(0.25*IT_0001*IT_0532);
    const ccomplex_t IT_0542 = cpow(IT_0541, 0.5);
    const ccomplex_t IT_0543 = ( ((m_W
) >= (m_tau
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0111*(IT_0533 + IT_0540)*IT_0542
) : (0
));
    const ccomplex_t IT_0544 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0004*U_sb_00;
    const ccomplex_t IT_0545 = -IT_0544;
    const ccomplex_t IT_0546 = -IT_0545;
    const ccomplex_t IT_0547 = pow(m_sc_L, 2);
    const ccomplex_t IT_0548 = pow(m_sb_1, 2);
    const ccomplex_t IT_0549 = 0.5*IT_0013 + 0.5*IT_0547 + (-0.5)*IT_0548;
    const ccomplex_t IT_0550 = cpow(IT_0549, 2);
    const ccomplex_t IT_0551 = IT_0001*IT_0550;
    const ccomplex_t IT_0552 = cpow(IT_0013 + IT_0547 + -IT_0548, 2);
    const ccomplex_t IT_0553 = cabs(-IT_0547 + 0.25*IT_0001*IT_0552);
    const ccomplex_t IT_0554 = cpow(IT_0553, 0.5);
    const ccomplex_t IT_0555 = ( ((m_W
) >= (m_sb_1 + m_sc_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0546*conj(IT_0546)*(IT_0547 + -IT_0551)*IT_0554
) : 
      (0
));
    const ccomplex_t IT_0556 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0004*U_sb_01;
    const ccomplex_t IT_0557 = -IT_0556;
    const ccomplex_t IT_0558 = -IT_0557;
    const ccomplex_t IT_0559 = pow(m_sb_2, 2);
    const ccomplex_t IT_0560 = 0.5*IT_0013 + 0.5*IT_0547 + (-0.5)*IT_0559;
    const ccomplex_t IT_0561 = cpow(IT_0560, 2);
    const ccomplex_t IT_0562 = IT_0001*IT_0561;
    const ccomplex_t IT_0563 = cpow(IT_0013 + IT_0547 + -IT_0559, 2);
    const ccomplex_t IT_0564 = cabs(-IT_0547 + 0.25*IT_0001*IT_0563);
    const ccomplex_t IT_0565 = cpow(IT_0564, 0.5);
    const ccomplex_t IT_0566 = ( ((m_W
) >= (m_sb_2 + m_sc_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0558*conj(IT_0558)*(IT_0547 + -IT_0562)*IT_0565
) : 
      (0
));
    const ccomplex_t IT_0567 = -IT_0094;
    const ccomplex_t IT_0568 = -IT_0567;
    const ccomplex_t IT_0569 = pow(m_sd_L, 2);
    const ccomplex_t IT_0570 = 0.5*IT_0013 + 0.5*IT_0547 + (-0.5)*IT_0569;
    const ccomplex_t IT_0571 = cpow(IT_0570, 2);
    const ccomplex_t IT_0572 = IT_0001*IT_0571;
    const ccomplex_t IT_0573 = cpow(IT_0013 + IT_0547 + -IT_0569, 2);
    const ccomplex_t IT_0574 = cabs(-IT_0547 + 0.25*IT_0001*IT_0573);
    const ccomplex_t IT_0575 = cpow(IT_0574, 0.5);
    const ccomplex_t IT_0576 = ( ((m_W
) >= (m_sc_L + m_sd_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0568*conj(IT_0568)*(IT_0547 + -IT_0572)*IT_0575
) : 
      (0
));
    const ccomplex_t IT_0577 = -IT_0140;
    const ccomplex_t IT_0578 = -IT_0577;
    const ccomplex_t IT_0579 = pow(m_ss_L, 2);
    const ccomplex_t IT_0580 = 0.5*IT_0013 + 0.5*IT_0547 + (-0.5)*IT_0579;
    const ccomplex_t IT_0581 = cpow(IT_0580, 2);
    const ccomplex_t IT_0582 = IT_0001*IT_0581;
    const ccomplex_t IT_0583 = cpow(IT_0013 + IT_0547 + -IT_0579, 2);
    const ccomplex_t IT_0584 = cabs(-IT_0547 + 0.25*IT_0001*IT_0583);
    const ccomplex_t IT_0585 = cpow(IT_0584, 0.5);
    const ccomplex_t IT_0586 = ( ((m_W
) >= (m_sc_L + m_ss_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0578*conj(IT_0578)*(IT_0547 + -IT_0582)*IT_0585
) : 
      (0
));
    const ccomplex_t IT_0587 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_00*conj(U_st_00);
    const ccomplex_t IT_0588 = -IT_0587;
    const ccomplex_t IT_0589 = -IT_0588;
    const ccomplex_t IT_0590 = pow(m_st_1, 2);
    const ccomplex_t IT_0591 = 0.5*IT_0013 + (-0.5)*IT_0548 + 0.5*IT_0590;
    const ccomplex_t IT_0592 = cpow(IT_0591, 2);
    const ccomplex_t IT_0593 = IT_0001*IT_0592;
    const ccomplex_t IT_0594 = cpow(IT_0013 + -IT_0548 + IT_0590, 2);
    const ccomplex_t IT_0595 = cabs(-IT_0590 + 0.25*IT_0001*IT_0594);
    const ccomplex_t IT_0596 = cpow(IT_0595, 0.5);
    const ccomplex_t IT_0597 = ( ((m_W
) >= (m_sb_1 + m_st_1
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0589*conj(IT_0589)*(IT_0590 + -IT_0593)*IT_0596
) : 
      (0
));
    const ccomplex_t IT_0598 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_01*conj(U_st_00);
    const ccomplex_t IT_0599 = -IT_0598;
    const ccomplex_t IT_0600 = -IT_0599;
    const ccomplex_t IT_0601 = 0.5*IT_0013 + (-0.5)*IT_0559 + 0.5*IT_0590;
    const ccomplex_t IT_0602 = cpow(IT_0601, 2);
    const ccomplex_t IT_0603 = IT_0001*IT_0602;
    const ccomplex_t IT_0604 = cpow(IT_0013 + -IT_0559 + IT_0590, 2);
    const ccomplex_t IT_0605 = cabs(-IT_0590 + 0.25*IT_0001*IT_0604);
    const ccomplex_t IT_0606 = cpow(IT_0605, 0.5);
    const ccomplex_t IT_0607 = ( ((m_W
) >= (m_sb_2 + m_st_1
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0600*conj(IT_0600)*(IT_0590 + -IT_0603)*IT_0606
) : 
      (0
));
    const ccomplex_t IT_0608 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0004*conj(U_st_00);
    const ccomplex_t IT_0609 = -IT_0608;
    const ccomplex_t IT_0610 = 0.5*IT_0013 + 0.5*IT_0569 + (-0.5)*IT_0590;
    const ccomplex_t IT_0611 = cpow(IT_0610, 2);
    const ccomplex_t IT_0612 = IT_0001*IT_0611;
    const ccomplex_t IT_0613 = cpow(IT_0013 + IT_0569 + -IT_0590, 2);
    const ccomplex_t IT_0614 = cabs(-IT_0569 + 0.25*IT_0001*IT_0613);
    const ccomplex_t IT_0615 = cpow(IT_0614, 0.5);
    const ccomplex_t IT_0616 = ( ((m_W
) >= (m_sd_L + m_st_1
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0609*conj(IT_0609)*(IT_0569 + -IT_0612)*IT_0615
) : 
      (0
));
    const ccomplex_t IT_0617 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0004*conj(U_st_00);
    const ccomplex_t IT_0618 = -IT_0617;
    const ccomplex_t IT_0619 = 0.5*IT_0013 + 0.5*IT_0579 + (-0.5)*IT_0590;
    const ccomplex_t IT_0620 = cpow(IT_0619, 2);
    const ccomplex_t IT_0621 = IT_0001*IT_0620;
    const ccomplex_t IT_0622 = cpow(IT_0013 + IT_0579 + -IT_0590, 2);
    const ccomplex_t IT_0623 = cabs(-IT_0579 + 0.25*IT_0001*IT_0622);
    const ccomplex_t IT_0624 = cpow(IT_0623, 0.5);
    const ccomplex_t IT_0625 = ( ((m_W
) >= (m_ss_L + m_st_1
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0618*conj(IT_0618)*(IT_0579 + -IT_0621)*IT_0624
) : 
      (0
));
    const ccomplex_t IT_0626 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_00*conj(U_st_01);
    const ccomplex_t IT_0627 = -IT_0626;
    const ccomplex_t IT_0628 = -IT_0627;
    const ccomplex_t IT_0629 = pow(m_st_2, 2);
    const ccomplex_t IT_0630 = 0.5*IT_0013 + (-0.5)*IT_0548 + 0.5*IT_0629;
    const ccomplex_t IT_0631 = cpow(IT_0630, 2);
    const ccomplex_t IT_0632 = IT_0001*IT_0631;
    const ccomplex_t IT_0633 = cpow(IT_0013 + -IT_0548 + IT_0629, 2);
    const ccomplex_t IT_0634 = cabs(-IT_0629 + 0.25*IT_0001*IT_0633);
    const ccomplex_t IT_0635 = cpow(IT_0634, 0.5);
    const ccomplex_t IT_0636 = ( ((m_W
) >= (m_sb_1 + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0628*conj(IT_0628)*(IT_0629 + -IT_0632)*IT_0635
) : 
      (0
));
    const ccomplex_t IT_0637 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*U_sb_01*conj(U_st_01);
    const ccomplex_t IT_0638 = -IT_0637;
    const ccomplex_t IT_0639 = -IT_0638;
    const ccomplex_t IT_0640 = 0.5*IT_0013 + (-0.5)*IT_0559 + 0.5*IT_0629;
    const ccomplex_t IT_0641 = cpow(IT_0640, 2);
    const ccomplex_t IT_0642 = IT_0001*IT_0641;
    const ccomplex_t IT_0643 = cpow(IT_0013 + -IT_0559 + IT_0629, 2);
    const ccomplex_t IT_0644 = cabs(-IT_0629 + 0.25*IT_0001*IT_0643);
    const ccomplex_t IT_0645 = cpow(IT_0644, 0.5);
    const ccomplex_t IT_0646 = ( ((m_W
) >= (m_sb_2 + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0639*conj(IT_0639)*(IT_0629 + -IT_0642)*IT_0645
) : 
      (0
));
    const ccomplex_t IT_0647 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0004*conj(U_st_01);
    const ccomplex_t IT_0648 = -IT_0647;
    const ccomplex_t IT_0649 = 0.5*IT_0013 + 0.5*IT_0569 + (-0.5)*IT_0629;
    const ccomplex_t IT_0650 = cpow(IT_0649, 2);
    const ccomplex_t IT_0651 = IT_0001*IT_0650;
    const ccomplex_t IT_0652 = cpow(IT_0013 + IT_0569 + -IT_0629, 2);
    const ccomplex_t IT_0653 = cabs(-IT_0569 + 0.25*IT_0001*IT_0652);
    const ccomplex_t IT_0654 = cpow(IT_0653, 0.5);
    const ccomplex_t IT_0655 = ( ((m_W
) >= (m_sd_L + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0648*conj(IT_0648)*(IT_0569 + -IT_0651)*IT_0654
) : 
      (0
));
    const ccomplex_t IT_0656 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0004*conj(U_st_01);
    const ccomplex_t IT_0657 = -IT_0656;
    const ccomplex_t IT_0658 = 0.5*IT_0013 + 0.5*IT_0579 + (-0.5)*IT_0629;
    const ccomplex_t IT_0659 = cpow(IT_0658, 2);
    const ccomplex_t IT_0660 = IT_0001*IT_0659;
    const ccomplex_t IT_0661 = cpow(IT_0013 + IT_0579 + -IT_0629, 2);
    const ccomplex_t IT_0662 = cabs(-IT_0579 + 0.25*IT_0001*IT_0661);
    const ccomplex_t IT_0663 = cpow(IT_0662, 0.5);
    const ccomplex_t IT_0664 = ( ((m_W
) >= (m_ss_L + m_st_2
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0657*conj(IT_0657)*(IT_0579 + -IT_0660)*IT_0663
) : 
      (0
));
    const ccomplex_t IT_0665 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *IT_0198*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0666 = -IT_0665;
    const ccomplex_t IT_0667 = -IT_0666;
    const ccomplex_t IT_0668 = pow(m_su_L, 2);
    const ccomplex_t IT_0669 = 0.5*IT_0013 + (-0.5)*IT_0548 + 0.5*IT_0668;
    const ccomplex_t IT_0670 = cpow(IT_0669, 2);
    const ccomplex_t IT_0671 = IT_0001*IT_0670;
    const ccomplex_t IT_0672 = cpow(IT_0013 + -IT_0548 + IT_0668, 2);
    const ccomplex_t IT_0673 = cabs(-IT_0668 + 0.25*IT_0001*IT_0672);
    const ccomplex_t IT_0674 = cpow(IT_0673, 0.5);
    const ccomplex_t IT_0675 = ( ((m_W
) >= (m_sb_1 + m_su_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0667*conj(IT_0667)*(IT_0668 + -IT_0671)*IT_0674
) : 
      (0
));
    const ccomplex_t IT_0676 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *IT_0198*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0677 = -IT_0676;
    const ccomplex_t IT_0678 = -IT_0677;
    const ccomplex_t IT_0679 = 0.5*IT_0013 + (-0.5)*IT_0559 + 0.5*IT_0668;
    const ccomplex_t IT_0680 = cpow(IT_0679, 2);
    const ccomplex_t IT_0681 = IT_0001*IT_0680;
    const ccomplex_t IT_0682 = cpow(IT_0013 + -IT_0559 + IT_0668, 2);
    const ccomplex_t IT_0683 = cabs(-IT_0668 + 0.25*IT_0001*IT_0682);
    const ccomplex_t IT_0684 = cpow(IT_0683, 0.5);
    const ccomplex_t IT_0685 = ( ((m_W
) >= (m_sb_2 + m_su_L
)) ? ((-0.125)
      *IT_0000*IT_0001*IT_0678*conj(IT_0678)*(IT_0668 + -IT_0681)*IT_0684
) : 
      (0
));
    const ccomplex_t IT_0686 = -IT_0215;
    const ccomplex_t IT_0687 = -IT_0686;
    const ccomplex_t IT_0688 = cpow(IT_0687, 2);
    const ccomplex_t IT_0689 = 0.5*IT_0013 + (-0.5)*IT_0569 + 0.5*IT_0668;
    const ccomplex_t IT_0690 = cpow(IT_0689, 2);
    const ccomplex_t IT_0691 = IT_0001*IT_0690;
    const ccomplex_t IT_0692 = cpow(IT_0013 + -IT_0569 + IT_0668, 2);
    const ccomplex_t IT_0693 = cabs(-IT_0668 + 0.25*IT_0001*IT_0692);
    const ccomplex_t IT_0694 = cpow(IT_0693, 0.5);
    const ccomplex_t IT_0695 = ( ((m_W
) >= (m_sd_L + m_su_L
)) ? (0.125
      *IT_0000*IT_0001*IT_0688*(IT_0668 + -IT_0691)*IT_0694
) : (0
));
    const ccomplex_t IT_0696 = -IT_0230;
    const ccomplex_t IT_0697 = -IT_0696;
    const ccomplex_t IT_0698 = cpow(IT_0697, 2);
    const ccomplex_t IT_0699 = 0.5*IT_0013 + (-0.5)*IT_0579 + 0.5*IT_0668;
    const ccomplex_t IT_0700 = cpow(IT_0699, 2);
    const ccomplex_t IT_0701 = IT_0001*IT_0700;
    const ccomplex_t IT_0702 = cpow(IT_0013 + -IT_0579 + IT_0668, 2);
    const ccomplex_t IT_0703 = cabs(-IT_0668 + 0.25*IT_0001*IT_0702);
    const ccomplex_t IT_0704 = cpow(IT_0703, 0.5);
    const ccomplex_t IT_0705 = ( ((m_W
) >= (m_ss_L + m_su_L
)) ? (0.125
      *IT_0000*IT_0001*IT_0698*(IT_0668 + -IT_0701)*IT_0704
) : (0
));
    const ccomplex_t IT_0706 = -IT_0109;
    const ccomplex_t IT_0707 = cpow(IT_0706, 2);
    const ccomplex_t IT_0708 = pow(m_se_L, 2);
    const ccomplex_t IT_0709 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0710 = 0.5*IT_0013 + 0.5*IT_0708 + (-0.5)*IT_0709;
    const ccomplex_t IT_0711 = cpow(IT_0710, 2);
    const ccomplex_t IT_0712 = IT_0001*IT_0711;
    const ccomplex_t IT_0713 = cpow(IT_0013 + IT_0708 + -IT_0709, 2);
    const ccomplex_t IT_0714 = cabs(-IT_0708 + 0.25*IT_0001*IT_0713);
    const ccomplex_t IT_0715 = cpow(IT_0714, 0.5);
    const ccomplex_t IT_0716 = ( ((m_W
) >= (m_se_L + m_snu_mu
)) ? 
      (0.0416666666666667*IT_0000*IT_0001*IT_0707*(IT_0708 + -IT_0712)*IT_0715
)
       : (0
));
    const ccomplex_t IT_0717 = pow(m_smu_L, 2);
    const ccomplex_t IT_0718 = 0.5*IT_0013 + (-0.5)*IT_0709 + 0.5*IT_0717;
    const ccomplex_t IT_0719 = cpow(IT_0718, 2);
    const ccomplex_t IT_0720 = IT_0001*IT_0719;
    const ccomplex_t IT_0721 = cpow(IT_0013 + -IT_0709 + IT_0717, 2);
    const ccomplex_t IT_0722 = cabs(-IT_0717 + 0.25*IT_0001*IT_0721);
    const ccomplex_t IT_0723 = cpow(IT_0722, 0.5);
    const ccomplex_t IT_0724 = ( ((m_W
) >= (m_smu_L + m_snu_mu
)) ? 
      (0.0416666666666667*IT_0000*IT_0001*IT_0707*(IT_0717 + -IT_0720)*IT_0723
)
       : (0
));
    const ccomplex_t IT_0725 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *U_stau_00;
    const ccomplex_t IT_0726 = -IT_0725;
    const ccomplex_t IT_0727 = -IT_0726;
    const ccomplex_t IT_0728 = pow(m_stau_1, 2);
    const ccomplex_t IT_0729 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0730 = 0.5*IT_0013 + (-0.5)*IT_0728 + 0.5*IT_0729;
    const ccomplex_t IT_0731 = cpow(IT_0730, 2);
    const ccomplex_t IT_0732 = IT_0001*IT_0731;
    const ccomplex_t IT_0733 = cpow(IT_0013 + -IT_0728 + IT_0729, 2);
    const ccomplex_t IT_0734 = cabs(-IT_0729 + 0.25*IT_0001*IT_0733);
    const ccomplex_t IT_0735 = cpow(IT_0734, 0.5);
    const ccomplex_t IT_0736 = ( ((m_W
) >= (m_stau_1 + m_snu_tau
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0727*conj(IT_0727)*(IT_0729 + 
      -IT_0732)*IT_0735
) : (0
));
    const ccomplex_t IT_0737 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *U_stau_01;
    const ccomplex_t IT_0738 = -IT_0737;
    const ccomplex_t IT_0739 = -IT_0738;
    const ccomplex_t IT_0740 = pow(m_stau_2, 2);
    const ccomplex_t IT_0741 = 0.5*IT_0013 + 0.5*IT_0729 + (-0.5)*IT_0740;
    const ccomplex_t IT_0742 = cpow(IT_0741, 2);
    const ccomplex_t IT_0743 = IT_0001*IT_0742;
    const ccomplex_t IT_0744 = cpow(IT_0013 + IT_0729 + -IT_0740, 2);
    const ccomplex_t IT_0745 = cabs(-IT_0729 + 0.25*IT_0001*IT_0744);
    const ccomplex_t IT_0746 = cpow(IT_0745, 0.5);
    const ccomplex_t IT_0747 = ( ((m_W
) >= (m_stau_2 + m_snu_tau
)) ? ((
      -0.0416666666666667)*IT_0000*IT_0001*IT_0739*conj(IT_0739)*(IT_0729 + 
      -IT_0743)*IT_0746
) : (0
));
    const ccomplex_t IT_0748 = IT_0076 + IT_0093 + IT_0108 + IT_0124 + IT_0139
       + IT_0154 + IT_0169 + IT_0183 + IT_0197 + IT_0214 + IT_0229 + IT_0244 +
       IT_0255 + IT_0267 + IT_0277 + IT_0287 + IT_0300 + IT_0330 + IT_0359 +
       IT_0388 + IT_0416 + IT_0445 + IT_0473 + IT_0502 + IT_0530 + IT_0543 +
       IT_0555 + IT_0566 + IT_0576 + IT_0586 + IT_0597 + IT_0607 + IT_0616 +
       IT_0625 + IT_0636 + IT_0646 + IT_0655 + IT_0664 + IT_0675 + IT_0685 +
       IT_0695 + IT_0705 + IT_0716 + IT_0724 + IT_0736 + IT_0747;
    return create_ccomplex_return(IT_0748);
}

