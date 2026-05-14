#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_t, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0006 = (-0.666666666666667)*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = pow(m_Z, 2);
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = pow(m_t, 2);
    const ccomplex_t IT_0011 = IT_0009 + IT_0010;
    const ccomplex_t IT_0012 = pow(m_Z, -2);
    const ccomplex_t IT_0013 = IT_0008*IT_0011;
    const ccomplex_t IT_0014 = 0.5*IT_0008;
    const ccomplex_t IT_0015 = pow(m_Z, 4);
    const ccomplex_t IT_0016 = -IT_0008 + (-0.25)*IT_0001*IT_0008*(IT_0008 + (
      -2)*IT_0010) + 0.25*IT_0001*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (-2)*IT_0017;
    const ccomplex_t IT_0019 = IT_0013 + IT_0018;
    const ccomplex_t IT_0020 = IT_0012*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0011 + IT_0021;
    const ccomplex_t IT_0023 = IT_0003*IT_0004;
    const ccomplex_t IT_0024 = e_em*IT_0023;
    const ccomplex_t IT_0025 = cpow(IT_0004, -1);
    const ccomplex_t IT_0026 = IT_0002*IT_0025;
    const ccomplex_t IT_0027 = e_em*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + (-3)*IT_0027);
    const ccomplex_t IT_0029 = (-0.166666666666667)*IT_0028;
    const ccomplex_t IT_0030 = cabs(-IT_0008 + 0.25*IT_0001*IT_0015);
    const ccomplex_t IT_0031 = cpow(IT_0030, 0.5);
    const ccomplex_t IT_0032 = ( ((m_t
) >= (m_Z + m_t
)) ? ((-0.25)*IT_0000
      *IT_0001*(IT_0007*IT_0022 + (-3)*IT_0029*(IT_0006*IT_0010 + (
      -0.333333333333333)*IT_0022*IT_0029))*IT_0031
) : (0
));
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = cos(alpha);
    const ccomplex_t IT_0035 = sin(beta);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_t*e_em*IT_0025*IT_0033
      *IT_0034*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = cpow(IT_0038, 2);
    const ccomplex_t IT_0040 = pow(m_h, 2);
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = IT_0010 + IT_0041;
    const ccomplex_t IT_0043 = pow(m_h, 4);
    const ccomplex_t IT_0044 = cabs(-IT_0040 + 0.25*IT_0001*IT_0043);
    const ccomplex_t IT_0045 = cpow(IT_0044, 0.5);
    const ccomplex_t IT_0046 = ( ((m_t
) >= (m_h + m_t
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0039*(IT_0010 + IT_0042)*IT_0045
) : (0
));
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0025;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = pow(m_W, 2);
    const ccomplex_t IT_0050 = pow(m_b, 2);
    const ccomplex_t IT_0051 = 0.5*IT_0010 + (-0.5)*IT_0049 + 0.5*IT_0050;
    const ccomplex_t IT_0052 = pow(m_W, -2);
    const ccomplex_t IT_0053 = IT_0049*IT_0051;
    const ccomplex_t IT_0054 = 0.5*IT_0010 + 0.5*IT_0049 + (-0.5)*IT_0050;
    const ccomplex_t IT_0055 = cpow(IT_0010 + IT_0049 + -IT_0050, 2);
    const ccomplex_t IT_0056 = -IT_0049 + 0.25*IT_0001*(IT_0010 + -IT_0049 +
       IT_0050)*(IT_0010 + IT_0049 + -IT_0050) + 0.25*IT_0001*IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = (-2)*IT_0057;
    const ccomplex_t IT_0059 = IT_0053 + IT_0058;
    const ccomplex_t IT_0060 = IT_0052*IT_0059;
    const ccomplex_t IT_0061 = cabs(-IT_0049 + 0.25*IT_0001*IT_0055);
    const ccomplex_t IT_0062 = cpow(IT_0061, 0.5);
    const ccomplex_t IT_0063 = ( ((m_t
) >= (m_W + m_b
)) ? (0.25*IT_0000
      *IT_0001*IT_0048*conj(IT_0048)*(IT_0051 + (-0.5)*IT_0060)*IT_0062
) : (0
));
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0025;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = pow(m_d, 2);
    const ccomplex_t IT_0067 = 0.5*IT_0010 + (-0.5)*IT_0049 + 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0049*IT_0067;
    const ccomplex_t IT_0069 = 0.5*IT_0010 + 0.5*IT_0049 + (-0.5)*IT_0066;
    const ccomplex_t IT_0070 = cpow(IT_0010 + IT_0049 + -IT_0066, 2);
    const ccomplex_t IT_0071 = -IT_0049 + 0.25*IT_0001*(IT_0010 + -IT_0049 +
       IT_0066)*(IT_0010 + IT_0049 + -IT_0066) + 0.25*IT_0001*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = (-2)*IT_0072;
    const ccomplex_t IT_0074 = IT_0068 + IT_0073;
    const ccomplex_t IT_0075 = IT_0052*IT_0074;
    const ccomplex_t IT_0076 = cabs(-IT_0049 + 0.25*IT_0001*IT_0070);
    const ccomplex_t IT_0077 = cpow(IT_0076, 0.5);
    const ccomplex_t IT_0078 = ( ((m_t
) >= (m_W + m_d
)) ? (0.25*IT_0000
      *IT_0001*IT_0065*conj(IT_0065)*(IT_0067 + (-0.5)*IT_0075)*IT_0077
) : (0
));
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0025;
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = pow(m_s, 2);
    const ccomplex_t IT_0082 = 0.5*IT_0010 + (-0.5)*IT_0049 + 0.5*IT_0081;
    const ccomplex_t IT_0083 = IT_0049*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0010 + 0.5*IT_0049 + (-0.5)*IT_0081;
    const ccomplex_t IT_0085 = cpow(IT_0010 + IT_0049 + -IT_0081, 2);
    const ccomplex_t IT_0086 = -IT_0049 + 0.25*IT_0001*(IT_0010 + -IT_0049 +
       IT_0081)*(IT_0010 + IT_0049 + -IT_0081) + 0.25*IT_0001*IT_0085;
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = (-2)*IT_0087;
    const ccomplex_t IT_0089 = IT_0083 + IT_0088;
    const ccomplex_t IT_0090 = IT_0052*IT_0089;
    const ccomplex_t IT_0091 = cabs(-IT_0049 + 0.25*IT_0001*IT_0085);
    const ccomplex_t IT_0092 = cpow(IT_0091, 0.5);
    const ccomplex_t IT_0093 = ( ((m_t
) >= (m_W + m_s
)) ? (0.25*IT_0000
      *IT_0001*IT_0080*conj(IT_0080)*(IT_0082 + (-0.5)*IT_0090)*IT_0092
) : (0
));
    const ccomplex_t IT_0094 = cos(beta);
    const ccomplex_t IT_0095 = m_t*e_em*IT_0025*IT_0033*IT_0036*IT_0094;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = cpow(IT_0096, 2);
    const ccomplex_t IT_0098 = pow(m_A0, 2);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0010 + IT_0099;
    const ccomplex_t IT_0101 = pow(m_A0, 4);
    const ccomplex_t IT_0102 = cabs(-IT_0098 + 0.25*IT_0001*IT_0101);
    const ccomplex_t IT_0103 = cpow(IT_0102, 0.5);
    const ccomplex_t IT_0104 = ( ((m_t
) >= (m_t + m_A0
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0097*(IT_0010 + -IT_0100)*IT_0103
) : (0
));
    const ccomplex_t IT_0105 = sin(alpha);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_t*e_em*IT_0025*IT_0033
      *IT_0036*IT_0105;
    const ccomplex_t IT_0107 = (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = cpow(IT_0107, 2);
    const ccomplex_t IT_0109 = pow(m_H0, 2);
    const ccomplex_t IT_0110 = (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = IT_0010 + IT_0110;
    const ccomplex_t IT_0112 = pow(m_H0, 4);
    const ccomplex_t IT_0113 = cabs(-IT_0109 + 0.25*IT_0001*IT_0112);
    const ccomplex_t IT_0114 = cpow(IT_0113, 0.5);
    const ccomplex_t IT_0115 = ( ((m_t
) >= (m_t + m_H0
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0108*(IT_0010 + IT_0111)*IT_0114
) : (0
));
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0025*IT_0033*IT_0036*IT_0094;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = cpow(IT_0094, -1);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0025*IT_0033*IT_0035*IT_0118;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = m_b*m_t;
    const ccomplex_t IT_0122 = pow(m_Hp, 2);
    const ccomplex_t IT_0123 = 0.5*IT_0010 + 0.5*IT_0050 + (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = cpow(-IT_0010 + IT_0050 + -IT_0122, 2);
    const ccomplex_t IT_0125 = cabs(-IT_0122 + 0.25*IT_0001*IT_0124);
    const ccomplex_t IT_0126 = cpow(IT_0125, 0.5);
    const ccomplex_t IT_0127 = ( ((m_t
) >= (m_b + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*((conj(IT_0117)*IT_0120 + IT_0117*conj(IT_0120))*IT_0121 + 
      (IT_0117*conj(IT_0117) + IT_0120*conj(IT_0120))*IT_0123)*IT_0126
) : (0
));
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *e_em*IT_0025*IT_0033*IT_0036*IT_0094;
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_td)
      *e_em*IT_0025*IT_0033*IT_0035*IT_0118;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = m_d*m_t;
    const ccomplex_t IT_0133 = 0.5*IT_0010 + 0.5*IT_0066 + (-0.5)*IT_0122;
    const ccomplex_t IT_0134 = cpow(-IT_0010 + IT_0066 + -IT_0122, 2);
    const ccomplex_t IT_0135 = cabs(-IT_0122 + 0.25*IT_0001*IT_0134);
    const ccomplex_t IT_0136 = cpow(IT_0135, 0.5);
    const ccomplex_t IT_0137 = ( ((m_t
) >= (m_d + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*((conj(IT_0129)*IT_0131 + IT_0129*conj(IT_0131))*IT_0132 + 
      (IT_0129*conj(IT_0129) + IT_0131*conj(IT_0131))*IT_0133)*IT_0136
) : (0
));
    const ccomplex_t IT_0138 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0025*IT_0033*IT_0036*IT_0094;
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0025*IT_0033*IT_0035*IT_0118;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = m_s*m_t;
    const ccomplex_t IT_0143 = 0.5*IT_0010 + 0.5*IT_0081 + (-0.5)*IT_0122;
    const ccomplex_t IT_0144 = cpow(-IT_0010 + IT_0081 + -IT_0122, 2);
    const ccomplex_t IT_0145 = cabs(-IT_0122 + 0.25*IT_0001*IT_0144);
    const ccomplex_t IT_0146 = cpow(IT_0145, 0.5);
    const ccomplex_t IT_0147 = ( ((m_t
) >= (m_s + m_Hp
)) ? (0.125*IT_0000
      *IT_0001*((conj(IT_0139)*IT_0141 + IT_0139*conj(IT_0141))*IT_0142 + 
      (IT_0139*conj(IT_0139) + IT_0141*conj(IT_0141))*IT_0143)*IT_0146
) : (0
));
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u1)*e_em*IT_0025*IT_0033*IT_0036*conj(U_sb_00);
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = V_tb*e_em*U_Wm1*conj(U_sb_00);
    const ccomplex_t IT_0151 = IT_0025*IT_0150;
    const ccomplex_t IT_0152 = IT_0025*IT_0118;
    const ccomplex_t IT_0153 = m_b*U_d1*V_tb*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0154 = IT_0152*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0151 + (-0.5)*IT_0155);
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = m_t*m_C_1;
    const ccomplex_t IT_0159 = pow(m_C_1, 2);
    const ccomplex_t IT_0160 = pow(m_sb_1, 2);
    const ccomplex_t IT_0161 = 0.5*IT_0010 + 0.5*IT_0159 + (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = cpow(IT_0010 + -IT_0159 + IT_0160, 2);
    const ccomplex_t IT_0163 = cabs(-IT_0160 + 0.25*IT_0001*IT_0162);
    const ccomplex_t IT_0164 = cpow(IT_0163, 0.5);
    const ccomplex_t IT_0165 = ( ((m_t
) >= (m_C_1 + m_sb_1
)) ? (0.125
      *IT_0000*IT_0001*((conj(IT_0149)*IT_0157 + IT_0149*conj(IT_0157))*IT_0158 
      + (IT_0149*conj(IT_0149) + IT_0157*conj(IT_0157))*IT_0161)*IT_0164
) : (0

      ));
    const ccomplex_t IT_0166 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0025*IT_0033*IT_0036*conj(U_sb_00);
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = V_tb*e_em*U_Wm2*conj(U_sb_00);
    const ccomplex_t IT_0169 = IT_0025*IT_0168;
    const ccomplex_t IT_0170 = m_b*U_d2*V_tb*e_em*IT_0033*conj(U_sb_10);
    const ccomplex_t IT_0171 = IT_0152*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*(IT_0169 + (-0.5)*IT_0172);
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = m_t*m_C_2;
    const ccomplex_t IT_0176 = pow(m_C_2, 2);
    const ccomplex_t IT_0177 = 0.5*IT_0010 + (-0.5)*IT_0160 + 0.5*IT_0176;
    const ccomplex_t IT_0178 = cpow(IT_0010 + IT_0160 + -IT_0176, 2);
    const ccomplex_t IT_0179 = cabs(-IT_0160 + 0.25*IT_0001*IT_0178);
    const ccomplex_t IT_0180 = cpow(IT_0179, 0.5);
    const ccomplex_t IT_0181 = ( ((m_t
) >= (m_C_2 + m_sb_1
)) ? (0.125
      *IT_0000*IT_0001*((conj(IT_0167)*IT_0174 + IT_0167*conj(IT_0174))*IT_0175 
      + (IT_0167*conj(IT_0167) + IT_0174*conj(IT_0174))*IT_0177)*IT_0180
) : (0

      ));
    const ccomplex_t IT_0182 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u1)*e_em*IT_0025*IT_0033*IT_0036*conj(U_sb_01);
    const ccomplex_t IT_0183 = 0.5*IT_0182;
    const ccomplex_t IT_0184 = V_tb*e_em*U_Wm1*conj(U_sb_01);
    const ccomplex_t IT_0185 = IT_0025*IT_0184;
    const ccomplex_t IT_0186 = m_b*U_d1*V_tb*e_em*IT_0033*conj(U_sb_11);
    const ccomplex_t IT_0187 = IT_0152*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0185 + (-0.5)*IT_0188);
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = pow(m_sb_2, 2);
    const ccomplex_t IT_0192 = 0.5*IT_0010 + 0.5*IT_0159 + (-0.5)*IT_0191;
    const ccomplex_t IT_0193 = cpow(IT_0010 + -IT_0159 + IT_0191, 2);
    const ccomplex_t IT_0194 = cabs(-IT_0191 + 0.25*IT_0001*IT_0193);
    const ccomplex_t IT_0195 = cpow(IT_0194, 0.5);
    const ccomplex_t IT_0196 = ( ((m_t
) >= (m_C_1 + m_sb_2
)) ? (0.125
      *IT_0000*IT_0001*(IT_0158*(conj(IT_0183)*IT_0190 + IT_0183*conj(IT_0190)) 
      + (IT_0183*conj(IT_0183) + IT_0190*conj(IT_0190))*IT_0192)*IT_0195
) : (0

      ));
    const ccomplex_t IT_0197 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*conj
      (V_u2)*e_em*IT_0025*IT_0033*IT_0036*conj(U_sb_01);
    const ccomplex_t IT_0198 = 0.5*IT_0197;
    const ccomplex_t IT_0199 = V_tb*e_em*U_Wm2*conj(U_sb_01);
    const ccomplex_t IT_0200 = IT_0025*IT_0199;
    const ccomplex_t IT_0201 = m_b*U_d2*V_tb*e_em*IT_0033*conj(U_sb_11);
    const ccomplex_t IT_0202 = IT_0152*IT_0201;
    const ccomplex_t IT_0203 = 1.4142135623731*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*(IT_0200 + (-0.5)*IT_0203);
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = 0.5*IT_0010 + 0.5*IT_0176 + (-0.5)*IT_0191;
    const ccomplex_t IT_0207 = cpow(IT_0010 + -IT_0176 + IT_0191, 2);
    const ccomplex_t IT_0208 = cabs(-IT_0191 + 0.25*IT_0001*IT_0207);
    const ccomplex_t IT_0209 = cpow(IT_0208, 0.5);
    const ccomplex_t IT_0210 = ( ((m_t
) >= (m_C_2 + m_sb_2
)) ? (0.125
      *IT_0000*IT_0001*(IT_0175*(conj(IT_0198)*IT_0205 + IT_0198*conj(IT_0205)) 
      + (IT_0198*conj(IT_0198) + IT_0205*conj(IT_0205))*IT_0206)*IT_0209
) : (0

      ));
    const ccomplex_t IT_0211 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *conj(V_u1)*e_em*IT_0025*IT_0033*IT_0036;
    const ccomplex_t IT_0212 = 0.5*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*conj(V_td)*e_em*U_Wm1*IT_0025;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = pow(m_sd_L, 2);
    const ccomplex_t IT_0216 = 0.5*IT_0010 + 0.5*IT_0159 + (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = cpow(IT_0010 + -IT_0159 + IT_0215, 2);
    const ccomplex_t IT_0218 = cabs(-IT_0215 + 0.25*IT_0001*IT_0217);
    const ccomplex_t IT_0219 = cpow(IT_0218, 0.5);
    const ccomplex_t IT_0220 = ( ((m_t
) >= (m_C_1 + m_sd_L
)) ? (0.125
      *IT_0000*IT_0001*(IT_0158*(conj(IT_0212)*IT_0214 + IT_0212*conj(IT_0214)) 
      + (IT_0212*conj(IT_0212) + IT_0214*conj(IT_0214))*IT_0216)*IT_0219
) : (0

      ));
    const ccomplex_t IT_0221 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_td)
      *conj(V_u2)*e_em*IT_0025*IT_0033*IT_0036;
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*conj(V_td)*e_em*U_Wm2*IT_0025;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = 0.5*IT_0010 + 0.5*IT_0176 + (-0.5)*IT_0215;
    const ccomplex_t IT_0226 = cpow(IT_0010 + -IT_0176 + IT_0215, 2);
    const ccomplex_t IT_0227 = cabs(-IT_0215 + 0.25*IT_0001*IT_0226);
    const ccomplex_t IT_0228 = cpow(IT_0227, 0.5);
    const ccomplex_t IT_0229 = ( ((m_t
) >= (m_C_2 + m_sd_L
)) ? (0.125
      *IT_0000*IT_0001*(IT_0175*(conj(IT_0222)*IT_0224 + IT_0222*conj(IT_0224)) 
      + (IT_0222*conj(IT_0222) + IT_0224*conj(IT_0224))*IT_0225)*IT_0228
) : (0

      ));
    const ccomplex_t IT_0230 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*conj
      (V_td)*e_em*IT_0025*IT_0033*IT_0118;
    const ccomplex_t IT_0231 = (-0.5)*IT_0230;
    const ccomplex_t IT_0232 = -IT_0231;
    const ccomplex_t IT_0233 = pow(m_sd_R, 2);
    const ccomplex_t IT_0234 = 0.5*IT_0010 + 0.5*IT_0159 + (-0.5)*IT_0233;
    const ccomplex_t IT_0235 = cpow(IT_0010 + -IT_0159 + IT_0233, 2);
    const ccomplex_t IT_0236 = cabs(-IT_0233 + 0.25*IT_0001*IT_0235);
    const ccomplex_t IT_0237 = cpow(IT_0236, 0.5);
    const ccomplex_t IT_0238 = ( ((m_t
) >= (m_C_1 + m_sd_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0232*conj(IT_0232)*IT_0234*IT_0237
) : (0
));
    const ccomplex_t IT_0239 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0025*IT_0033*IT_0118;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = -IT_0240;
    const ccomplex_t IT_0242 = 0.5*IT_0010 + 0.5*IT_0176 + (-0.5)*IT_0233;
    const ccomplex_t IT_0243 = cpow(IT_0010 + -IT_0176 + IT_0233, 2);
    const ccomplex_t IT_0244 = cabs(-IT_0233 + 0.25*IT_0001*IT_0243);
    const ccomplex_t IT_0245 = cpow(IT_0244, 0.5);
    const ccomplex_t IT_0246 = ( ((m_t
) >= (m_C_2 + m_sd_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0241*conj(IT_0241)*IT_0242*IT_0245
) : (0
));
    const ccomplex_t IT_0247 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *conj(V_u1)*e_em*IT_0025*IT_0033*IT_0036;
    const ccomplex_t IT_0248 = 0.5*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*conj(V_ts)*e_em*U_Wm1*IT_0025;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = pow(m_ss_L, 2);
    const ccomplex_t IT_0252 = 0.5*IT_0010 + 0.5*IT_0159 + (-0.5)*IT_0251;
    const ccomplex_t IT_0253 = cpow(IT_0010 + -IT_0159 + IT_0251, 2);
    const ccomplex_t IT_0254 = cabs(-IT_0251 + 0.25*IT_0001*IT_0253);
    const ccomplex_t IT_0255 = cpow(IT_0254, 0.5);
    const ccomplex_t IT_0256 = ( ((m_t
) >= (m_C_1 + m_ss_L
)) ? (0.125
      *IT_0000*IT_0001*(IT_0158*(conj(IT_0248)*IT_0250 + IT_0248*conj(IT_0250)) 
      + (IT_0248*conj(IT_0248) + IT_0250*conj(IT_0250))*IT_0252)*IT_0255
) : (0

      ));
    const ccomplex_t IT_0257 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *conj(V_u2)*e_em*IT_0025*IT_0033*IT_0036;
    const ccomplex_t IT_0258 = 0.5*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*conj(V_ts)*e_em*U_Wm2*IT_0025;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = 0.5*IT_0010 + 0.5*IT_0176 + (-0.5)*IT_0251;
    const ccomplex_t IT_0262 = cpow(IT_0010 + -IT_0176 + IT_0251, 2);
    const ccomplex_t IT_0263 = cabs(-IT_0251 + 0.25*IT_0001*IT_0262);
    const ccomplex_t IT_0264 = cpow(IT_0263, 0.5);
    const ccomplex_t IT_0265 = ( ((m_t
) >= (m_C_2 + m_ss_L
)) ? (0.125
      *IT_0000*IT_0001*(IT_0175*(conj(IT_0258)*IT_0260 + IT_0258*conj(IT_0260)) 
      + (IT_0258*conj(IT_0258) + IT_0260*conj(IT_0260))*IT_0261)*IT_0264
) : (0

      ));
    const ccomplex_t IT_0266 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0025*IT_0033*IT_0118;
    const ccomplex_t IT_0267 = (-0.5)*IT_0266;
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = pow(m_ss_R, 2);
    const ccomplex_t IT_0270 = 0.5*IT_0010 + 0.5*IT_0159 + (-0.5)*IT_0269;
    const ccomplex_t IT_0271 = cpow(IT_0010 + -IT_0159 + IT_0269, 2);
    const ccomplex_t IT_0272 = cabs(-IT_0269 + 0.25*IT_0001*IT_0271);
    const ccomplex_t IT_0273 = cpow(IT_0272, 0.5);
    const ccomplex_t IT_0274 = ( ((m_t
) >= (m_C_1 + m_ss_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0268*conj(IT_0268)*IT_0270*IT_0273
) : (0
));
    const ccomplex_t IT_0275 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0025*IT_0033*IT_0118;
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = -IT_0276;
    const ccomplex_t IT_0278 = 0.5*IT_0010 + 0.5*IT_0176 + (-0.5)*IT_0269;
    const ccomplex_t IT_0279 = cpow(IT_0010 + -IT_0176 + IT_0269, 2);
    const ccomplex_t IT_0280 = cabs(-IT_0269 + 0.25*IT_0001*IT_0279);
    const ccomplex_t IT_0281 = cpow(IT_0280, 0.5);
    const ccomplex_t IT_0282 = ( ((m_t
) >= (m_C_2 + m_ss_R
)) ? (0.125
      *IT_0000*IT_0001*IT_0277*conj(IT_0277)*IT_0278*IT_0281
) : (0
));
    const ccomplex_t IT_0283 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0284 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0285 = -IT_0284;
    const ccomplex_t IT_0286 = m_t*m_sG;
    const ccomplex_t IT_0287 = pow(m_sG, 2);
    const ccomplex_t IT_0288 = pow(m_st_1, 2);
    const ccomplex_t IT_0289 = 0.5*IT_0010 + 0.5*IT_0287 + (-0.5)*IT_0288;
    const ccomplex_t IT_0290 = 1.33333333333333*IT_0289;
    const ccomplex_t IT_0291 = cpow(IT_0010 + -IT_0287 + IT_0288, 2);
    const ccomplex_t IT_0292 = cabs(-IT_0288 + 0.25*IT_0001*IT_0291);
    const ccomplex_t IT_0293 = cpow(IT_0292, 0.5);
    const ccomplex_t IT_0294 = ( ((m_t
) >= (m_sG + m_st_1
)) ? (0.125*IT_0000
      *IT_0001*(1.33333333333333*conj(IT_0283)*(IT_0285*IT_0286 + 0.75*IT_0283
      *IT_0290) + 1.33333333333333*conj(IT_0285)*(IT_0283*IT_0286 + 0.75*IT_0285
      *IT_0290))*IT_0293
) : (0
));
    const ccomplex_t IT_0295 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0296 = IT_0003*IT_0295;
    const ccomplex_t IT_0297 = 1.4142135623731*IT_0296;
    const ccomplex_t IT_0298 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0299 = IT_0025*IT_0298;
    const ccomplex_t IT_0300 = 1.4142135623731*IT_0299;
    const ccomplex_t IT_0301 = IT_0025*IT_0036;
    const ccomplex_t IT_0302 = m_t*N_u1*e_em*IT_0033*conj(U_st_10);
    const ccomplex_t IT_0303 = IT_0301*IT_0302;
    const ccomplex_t IT_0304 = 1.4142135623731*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*(IT_0297 + 3*IT_0300 + 3
      *IT_0304);
    const ccomplex_t IT_0306 = 0.166666666666667*IT_0305;
    const ccomplex_t IT_0307 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0308 = IT_0003*IT_0307;
    const ccomplex_t IT_0309 = 1.4142135623731*IT_0308;
    const ccomplex_t IT_0310 = m_t*conj(N_u1)*e_em*IT_0033*conj(U_st_00);
    const ccomplex_t IT_0311 = IT_0301*IT_0310;
    const ccomplex_t IT_0312 = 1.4142135623731*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*(IT_0309 + (-0.75)*IT_0312);
    const ccomplex_t IT_0314 = 0.666666666666667*IT_0313;
    const ccomplex_t IT_0315 = -IT_0314;
    const ccomplex_t IT_0316 = m_t*m_N_1;
    const ccomplex_t IT_0317 = pow(m_N_1, 2);
    const ccomplex_t IT_0318 = 0.5*IT_0010 + (-0.5)*IT_0288 + 0.5*IT_0317;
    const ccomplex_t IT_0319 = cpow(IT_0010 + IT_0288 + -IT_0317, 2);
    const ccomplex_t IT_0320 = cabs(-IT_0288 + 0.25*IT_0001*IT_0319);
    const ccomplex_t IT_0321 = cpow(IT_0320, 0.5);
    const ccomplex_t IT_0322 = ( ((m_t
) >= (m_N_1 + m_st_1
)) ? (0.125
      *IT_0000*IT_0001*((conj(IT_0306)*IT_0315 + IT_0306*conj(IT_0315))*IT_0316 
      + (IT_0306*conj(IT_0306) + IT_0315*conj(IT_0315))*IT_0318)*IT_0321
) : (0

      ));
    const ccomplex_t IT_0323 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0324 = IT_0003*IT_0323;
    const ccomplex_t IT_0325 = 1.4142135623731*IT_0324;
    const ccomplex_t IT_0326 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0327 = IT_0025*IT_0326;
    const ccomplex_t IT_0328 = 1.4142135623731*IT_0327;
    const ccomplex_t IT_0329 = m_t*N_u2*e_em*IT_0033*conj(U_st_10);
    const ccomplex_t IT_0330 = IT_0301*IT_0329;
    const ccomplex_t IT_0331 = 1.4142135623731*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*(IT_0325 + 3*IT_0328 + 3
      *IT_0331);
    const ccomplex_t IT_0333 = 0.166666666666667*IT_0332;
    const ccomplex_t IT_0334 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0335 = IT_0003*IT_0334;
    const ccomplex_t IT_0336 = 1.4142135623731*IT_0335;
    const ccomplex_t IT_0337 = m_t*conj(N_u2)*e_em*IT_0033*conj(U_st_00);
    const ccomplex_t IT_0338 = IT_0301*IT_0337;
    const ccomplex_t IT_0339 = 1.4142135623731*IT_0338;
    const ccomplex_t IT_0340 = (0 + _Complex_I*1)*(IT_0336 + (-0.75)*IT_0339);
    const ccomplex_t IT_0341 = 0.666666666666667*IT_0340;
    const ccomplex_t IT_0342 = -IT_0341;
    const ccomplex_t IT_0343 = m_t*m_N_2;
    const ccomplex_t IT_0344 = pow(m_N_2, 2);
    const ccomplex_t IT_0345 = 0.5*IT_0010 + (-0.5)*IT_0288 + 0.5*IT_0344;
    const ccomplex_t IT_0346 = cpow(IT_0010 + IT_0288 + -IT_0344, 2);
    const ccomplex_t IT_0347 = cabs(-IT_0288 + 0.25*IT_0001*IT_0346);
    const ccomplex_t IT_0348 = cpow(IT_0347, 0.5);
    const ccomplex_t IT_0349 = ( ((m_t
) >= (m_N_2 + m_st_1
)) ? (0.125
      *IT_0000*IT_0001*((conj(IT_0333)*IT_0342 + IT_0333*conj(IT_0342))*IT_0343 
      + (IT_0333*conj(IT_0333) + IT_0342*conj(IT_0342))*IT_0345)*IT_0348
) : (0

      ));
    const ccomplex_t IT_0350 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0351 = IT_0003*IT_0350;
    const ccomplex_t IT_0352 = 1.4142135623731*IT_0351;
    const ccomplex_t IT_0353 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0354 = IT_0025*IT_0353;
    const ccomplex_t IT_0355 = 1.4142135623731*IT_0354;
    const ccomplex_t IT_0356 = m_t*N_u3*e_em*IT_0033*conj(U_st_10);
    const ccomplex_t IT_0357 = IT_0301*IT_0356;
    const ccomplex_t IT_0358 = 1.4142135623731*IT_0357;
    const ccomplex_t IT_0359 = (0 + _Complex_I*1)*(IT_0352 + 3*IT_0355 + 3
      *IT_0358);
    const ccomplex_t IT_0360 = 0.166666666666667*IT_0359;
    const ccomplex_t IT_0361 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0362 = IT_0003*IT_0361;
    const ccomplex_t IT_0363 = 1.4142135623731*IT_0362;
    const ccomplex_t IT_0364 = m_t*conj(N_u3)*e_em*IT_0033*conj(U_st_00);
    const ccomplex_t IT_0365 = IT_0301*IT_0364;
    const ccomplex_t IT_0366 = 1.4142135623731*IT_0365;
    const ccomplex_t IT_0367 = (0 + _Complex_I*1)*(IT_0363 + (-0.75)*IT_0366);
    const ccomplex_t IT_0368 = 0.666666666666667*IT_0367;
    const ccomplex_t IT_0369 = -IT_0368;
    const ccomplex_t IT_0370 = m_t*m_N_3;
    const ccomplex_t IT_0371 = pow(m_N_3, 2);
    const ccomplex_t IT_0372 = 0.5*IT_0010 + (-0.5)*IT_0288 + 0.5*IT_0371;
    const ccomplex_t IT_0373 = cpow(IT_0010 + IT_0288 + -IT_0371, 2);
    const ccomplex_t IT_0374 = cabs(-IT_0288 + 0.25*IT_0001*IT_0373);
    const ccomplex_t IT_0375 = cpow(IT_0374, 0.5);
    const ccomplex_t IT_0376 = ( ((m_t
) >= (m_N_3 + m_st_1
)) ? (0.125
      *IT_0000*IT_0001*((conj(IT_0360)*IT_0369 + IT_0360*conj(IT_0369))*IT_0370 
      + (IT_0360*conj(IT_0360) + IT_0369*conj(IT_0369))*IT_0372)*IT_0375
) : (0

      ));
    const ccomplex_t IT_0377 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0378 = IT_0003*IT_0377;
    const ccomplex_t IT_0379 = 1.4142135623731*IT_0378;
    const ccomplex_t IT_0380 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0381 = IT_0025*IT_0380;
    const ccomplex_t IT_0382 = 1.4142135623731*IT_0381;
    const ccomplex_t IT_0383 = m_t*N_u4*e_em*IT_0033*conj(U_st_10);
    const ccomplex_t IT_0384 = IT_0301*IT_0383;
    const ccomplex_t IT_0385 = 1.4142135623731*IT_0384;
    const ccomplex_t IT_0386 = (0 + _Complex_I*1)*(IT_0379 + 3*IT_0382 + 3
      *IT_0385);
    const ccomplex_t IT_0387 = 0.166666666666667*IT_0386;
    const ccomplex_t IT_0388 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0389 = IT_0003*IT_0388;
    const ccomplex_t IT_0390 = 1.4142135623731*IT_0389;
    const ccomplex_t IT_0391 = m_t*conj(N_u4)*e_em*IT_0033*conj(U_st_00);
    const ccomplex_t IT_0392 = IT_0301*IT_0391;
    const ccomplex_t IT_0393 = 1.4142135623731*IT_0392;
    const ccomplex_t IT_0394 = (0 + _Complex_I*1)*(IT_0390 + (-0.75)*IT_0393);
    const ccomplex_t IT_0395 = 0.666666666666667*IT_0394;
    const ccomplex_t IT_0396 = -IT_0395;
    const ccomplex_t IT_0397 = m_t*m_N_4;
    const ccomplex_t IT_0398 = pow(m_N_4, 2);
    const ccomplex_t IT_0399 = 0.5*IT_0010 + (-0.5)*IT_0288 + 0.5*IT_0398;
    const ccomplex_t IT_0400 = cpow(IT_0010 + IT_0288 + -IT_0398, 2);
    const ccomplex_t IT_0401 = cabs(-IT_0288 + 0.25*IT_0001*IT_0400);
    const ccomplex_t IT_0402 = cpow(IT_0401, 0.5);
    const ccomplex_t IT_0403 = ( ((m_t
) >= (m_N_4 + m_st_1
)) ? (0.125
      *IT_0000*IT_0001*((conj(IT_0387)*IT_0396 + IT_0387*conj(IT_0396))*IT_0397 
      + (IT_0387*conj(IT_0387) + IT_0396*conj(IT_0396))*IT_0399)*IT_0402
) : (0

      ));
    const ccomplex_t IT_0404 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_01);
    const ccomplex_t IT_0405 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_11);
    const ccomplex_t IT_0406 = -IT_0405;
    const ccomplex_t IT_0407 = pow(m_st_2, 2);
    const ccomplex_t IT_0408 = 0.5*IT_0010 + 0.5*IT_0287 + (-0.5)*IT_0407;
    const ccomplex_t IT_0409 = 1.33333333333333*IT_0408;
    const ccomplex_t IT_0410 = cpow(IT_0010 + -IT_0287 + IT_0407, 2);
    const ccomplex_t IT_0411 = cabs(-IT_0407 + 0.25*IT_0001*IT_0410);
    const ccomplex_t IT_0412 = cpow(IT_0411, 0.5);
    const ccomplex_t IT_0413 = ( ((m_t
) >= (m_sG + m_st_2
)) ? (0.125*IT_0000
      *IT_0001*(1.33333333333333*conj(IT_0404)*(IT_0286*IT_0406 + 0.75*IT_0404
      *IT_0409) + 1.33333333333333*conj(IT_0406)*(IT_0286*IT_0404 + 0.75*IT_0406
      *IT_0409))*IT_0412
) : (0
));
    const ccomplex_t IT_0414 = N_B1*e_em*conj(U_st_01);
    const ccomplex_t IT_0415 = IT_0003*IT_0414;
    const ccomplex_t IT_0416 = 1.4142135623731*IT_0415;
    const ccomplex_t IT_0417 = N_W1*e_em*conj(U_st_01);
    const ccomplex_t IT_0418 = IT_0025*IT_0417;
    const ccomplex_t IT_0419 = 1.4142135623731*IT_0418;
    const ccomplex_t IT_0420 = m_t*N_u1*e_em*IT_0033*conj(U_st_11);
    const ccomplex_t IT_0421 = IT_0301*IT_0420;
    const ccomplex_t IT_0422 = 1.4142135623731*IT_0421;
    const ccomplex_t IT_0423 = (0 + _Complex_I*1)*(IT_0416 + 3*IT_0419 + 3
      *IT_0422);
    const ccomplex_t IT_0424 = 0.166666666666667*IT_0423;
    const ccomplex_t IT_0425 = conj(N_B1)*e_em*conj(U_st_11);
    const ccomplex_t IT_0426 = IT_0003*IT_0425;
    const ccomplex_t IT_0427 = 1.4142135623731*IT_0426;
    const ccomplex_t IT_0428 = m_t*conj(N_u1)*e_em*IT_0033*conj(U_st_01);
    const ccomplex_t IT_0429 = IT_0301*IT_0428;
    const ccomplex_t IT_0430 = 1.4142135623731*IT_0429;
    const ccomplex_t IT_0431 = (0 + _Complex_I*1)*(IT_0427 + (-0.75)*IT_0430);
    const ccomplex_t IT_0432 = 0.666666666666667*IT_0431;
    const ccomplex_t IT_0433 = -IT_0432;
    const ccomplex_t IT_0434 = 0.5*IT_0010 + 0.5*IT_0317 + (-0.5)*IT_0407;
    const ccomplex_t IT_0435 = cpow(IT_0010 + -IT_0317 + IT_0407, 2);
    const ccomplex_t IT_0436 = cabs(-IT_0407 + 0.25*IT_0001*IT_0435);
    const ccomplex_t IT_0437 = cpow(IT_0436, 0.5);
    const ccomplex_t IT_0438 = ( ((m_t
) >= (m_N_1 + m_st_2
)) ? (0.125
      *IT_0000*IT_0001*(IT_0316*(conj(IT_0424)*IT_0433 + IT_0424*conj(IT_0433)) 
      + (IT_0424*conj(IT_0424) + IT_0433*conj(IT_0433))*IT_0434)*IT_0437
) : (0

      ));
    const ccomplex_t IT_0439 = N_B2*e_em*conj(U_st_01);
    const ccomplex_t IT_0440 = IT_0003*IT_0439;
    const ccomplex_t IT_0441 = 1.4142135623731*IT_0440;
    const ccomplex_t IT_0442 = N_W2*e_em*conj(U_st_01);
    const ccomplex_t IT_0443 = IT_0025*IT_0442;
    const ccomplex_t IT_0444 = 1.4142135623731*IT_0443;
    const ccomplex_t IT_0445 = m_t*N_u2*e_em*IT_0033*conj(U_st_11);
    const ccomplex_t IT_0446 = IT_0301*IT_0445;
    const ccomplex_t IT_0447 = 1.4142135623731*IT_0446;
    const ccomplex_t IT_0448 = (0 + _Complex_I*1)*(IT_0441 + 3*IT_0444 + 3
      *IT_0447);
    const ccomplex_t IT_0449 = 0.166666666666667*IT_0448;
    const ccomplex_t IT_0450 = conj(N_B2)*e_em*conj(U_st_11);
    const ccomplex_t IT_0451 = IT_0003*IT_0450;
    const ccomplex_t IT_0452 = 1.4142135623731*IT_0451;
    const ccomplex_t IT_0453 = m_t*conj(N_u2)*e_em*IT_0033*conj(U_st_01);
    const ccomplex_t IT_0454 = IT_0301*IT_0453;
    const ccomplex_t IT_0455 = 1.4142135623731*IT_0454;
    const ccomplex_t IT_0456 = (0 + _Complex_I*1)*(IT_0452 + (-0.75)*IT_0455);
    const ccomplex_t IT_0457 = 0.666666666666667*IT_0456;
    const ccomplex_t IT_0458 = -IT_0457;
    const ccomplex_t IT_0459 = 0.5*IT_0010 + 0.5*IT_0344 + (-0.5)*IT_0407;
    const ccomplex_t IT_0460 = cpow(IT_0010 + -IT_0344 + IT_0407, 2);
    const ccomplex_t IT_0461 = cabs(-IT_0407 + 0.25*IT_0001*IT_0460);
    const ccomplex_t IT_0462 = cpow(IT_0461, 0.5);
    const ccomplex_t IT_0463 = ( ((m_t
) >= (m_N_2 + m_st_2
)) ? (0.125
      *IT_0000*IT_0001*(IT_0343*(conj(IT_0449)*IT_0458 + IT_0449*conj(IT_0458)) 
      + (IT_0449*conj(IT_0449) + IT_0458*conj(IT_0458))*IT_0459)*IT_0462
) : (0

      ));
    const ccomplex_t IT_0464 = N_B3*e_em*conj(U_st_01);
    const ccomplex_t IT_0465 = IT_0003*IT_0464;
    const ccomplex_t IT_0466 = 1.4142135623731*IT_0465;
    const ccomplex_t IT_0467 = N_W3*e_em*conj(U_st_01);
    const ccomplex_t IT_0468 = IT_0025*IT_0467;
    const ccomplex_t IT_0469 = 1.4142135623731*IT_0468;
    const ccomplex_t IT_0470 = m_t*N_u3*e_em*IT_0033*conj(U_st_11);
    const ccomplex_t IT_0471 = IT_0301*IT_0470;
    const ccomplex_t IT_0472 = 1.4142135623731*IT_0471;
    const ccomplex_t IT_0473 = (0 + _Complex_I*1)*(IT_0466 + 3*IT_0469 + 3
      *IT_0472);
    const ccomplex_t IT_0474 = 0.166666666666667*IT_0473;
    const ccomplex_t IT_0475 = conj(N_B3)*e_em*conj(U_st_11);
    const ccomplex_t IT_0476 = IT_0003*IT_0475;
    const ccomplex_t IT_0477 = 1.4142135623731*IT_0476;
    const ccomplex_t IT_0478 = m_t*conj(N_u3)*e_em*IT_0033*conj(U_st_01);
    const ccomplex_t IT_0479 = IT_0301*IT_0478;
    const ccomplex_t IT_0480 = 1.4142135623731*IT_0479;
    const ccomplex_t IT_0481 = (0 + _Complex_I*1)*(IT_0477 + (-0.75)*IT_0480);
    const ccomplex_t IT_0482 = 0.666666666666667*IT_0481;
    const ccomplex_t IT_0483 = -IT_0482;
    const ccomplex_t IT_0484 = 0.5*IT_0010 + 0.5*IT_0371 + (-0.5)*IT_0407;
    const ccomplex_t IT_0485 = cpow(IT_0010 + -IT_0371 + IT_0407, 2);
    const ccomplex_t IT_0486 = cabs(-IT_0407 + 0.25*IT_0001*IT_0485);
    const ccomplex_t IT_0487 = cpow(IT_0486, 0.5);
    const ccomplex_t IT_0488 = ( ((m_t
) >= (m_N_3 + m_st_2
)) ? (0.125
      *IT_0000*IT_0001*(IT_0370*(conj(IT_0474)*IT_0483 + IT_0474*conj(IT_0483)) 
      + (IT_0474*conj(IT_0474) + IT_0483*conj(IT_0483))*IT_0484)*IT_0487
) : (0

      ));
    const ccomplex_t IT_0489 = N_B4*e_em*conj(U_st_01);
    const ccomplex_t IT_0490 = IT_0003*IT_0489;
    const ccomplex_t IT_0491 = 1.4142135623731*IT_0490;
    const ccomplex_t IT_0492 = N_W4*e_em*conj(U_st_01);
    const ccomplex_t IT_0493 = IT_0025*IT_0492;
    const ccomplex_t IT_0494 = 1.4142135623731*IT_0493;
    const ccomplex_t IT_0495 = m_t*N_u4*e_em*IT_0033*conj(U_st_11);
    const ccomplex_t IT_0496 = IT_0301*IT_0495;
    const ccomplex_t IT_0497 = 1.4142135623731*IT_0496;
    const ccomplex_t IT_0498 = (0 + _Complex_I*1)*(IT_0491 + 3*IT_0494 + 3
      *IT_0497);
    const ccomplex_t IT_0499 = 0.166666666666667*IT_0498;
    const ccomplex_t IT_0500 = conj(N_B4)*e_em*conj(U_st_11);
    const ccomplex_t IT_0501 = IT_0003*IT_0500;
    const ccomplex_t IT_0502 = 1.4142135623731*IT_0501;
    const ccomplex_t IT_0503 = m_t*conj(N_u4)*e_em*IT_0033*conj(U_st_01);
    const ccomplex_t IT_0504 = IT_0301*IT_0503;
    const ccomplex_t IT_0505 = 1.4142135623731*IT_0504;
    const ccomplex_t IT_0506 = (0 + _Complex_I*1)*(IT_0502 + (-0.75)*IT_0505);
    const ccomplex_t IT_0507 = 0.666666666666667*IT_0506;
    const ccomplex_t IT_0508 = -IT_0507;
    const ccomplex_t IT_0509 = 0.5*IT_0010 + 0.5*IT_0398 + (-0.5)*IT_0407;
    const ccomplex_t IT_0510 = cpow(IT_0010 + -IT_0398 + IT_0407, 2);
    const ccomplex_t IT_0511 = cabs(-IT_0407 + 0.25*IT_0001*IT_0510);
    const ccomplex_t IT_0512 = cpow(IT_0511, 0.5);
    const ccomplex_t IT_0513 = ( ((m_t
) >= (m_N_4 + m_st_2
)) ? (0.125
      *IT_0000*IT_0001*(IT_0397*(conj(IT_0499)*IT_0508 + IT_0499*conj(IT_0508)) 
      + (IT_0499*conj(IT_0499) + IT_0508*conj(IT_0508))*IT_0509)*IT_0512
) : (0

      ));
    const ccomplex_t IT_0514 = IT_0032 + IT_0046 + IT_0063 + IT_0078 + IT_0093
       + IT_0104 + IT_0115 + IT_0127 + IT_0137 + IT_0147 + IT_0165 + IT_0181 +
       IT_0196 + IT_0210 + IT_0220 + IT_0229 + IT_0238 + IT_0246 + IT_0256 +
       IT_0265 + IT_0274 + IT_0282 + IT_0294 + IT_0322 + IT_0349 + IT_0376 +
       IT_0403 + IT_0413 + IT_0438 + IT_0463 + IT_0488 + IT_0513;
    return create_ccomplex_return(IT_0514);
}

