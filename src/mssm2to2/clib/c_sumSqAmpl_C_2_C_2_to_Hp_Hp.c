#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_C_2_to_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_C_2_to_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(beta);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = (-2)*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0000*IT_0008;
    const ccomplex_t IT_0010 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = IT_0012 + IT_0015;
    const ccomplex_t IT_0017 = IT_0006 + IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = IT_0002*IT_0020;
    const ccomplex_t IT_0022 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 2*IT_0023;
    const ccomplex_t IT_0025 = IT_0008*IT_0020;
    const ccomplex_t IT_0026 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = IT_0028 + IT_0031;
    const ccomplex_t IT_0033 = IT_0024 + IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0019*IT_0035;
    const ccomplex_t IT_0037 = pow(m_Hp, 2);
    const ccomplex_t IT_0038 = pow(m_C_2, 2);
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0037 + IT_0038 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0043 = IT_0003*IT_0042;
    const ccomplex_t IT_0044 = (-2)*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0046 = IT_0009*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0049 = IT_0003*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = IT_0047 + IT_0050;
    const ccomplex_t IT_0052 = IT_0044 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0056 = IT_0021*IT_0055;
    const ccomplex_t IT_0057 = 2*IT_0056;
    const ccomplex_t IT_0058 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0059 = IT_0025*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0062 = IT_0021*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = IT_0060 + IT_0063;
    const ccomplex_t IT_0065 = IT_0057 + IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0054*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0037 + IT_0038 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0037 + IT_0038 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0073 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0074 = IT_0021*IT_0073;
    const ccomplex_t IT_0075 = 2*IT_0074;
    const ccomplex_t IT_0076 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0077 = IT_0025*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0080 = IT_0021*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = IT_0078 + IT_0081;
    const ccomplex_t IT_0083 = IT_0075 + IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0090 = IT_0009*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = IT_0088 + IT_0091;
    const ccomplex_t IT_0093 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0094 = IT_0003*IT_0093;
    const ccomplex_t IT_0095 = (-2)*IT_0094;
    const ccomplex_t IT_0096 = IT_0092 + IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0085*IT_0098;
    const ccomplex_t IT_0100 = IT_0072*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_23 + IT_0037 + IT_0038 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0099*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_23 + IT_0037 + IT_0038 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0036*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_13 + IT_0037 + IT_0038 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0068*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = cpow((-2)*s_13 + IT_0037 + IT_0038 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0112 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0113 = IT_0003*IT_0112;
    const ccomplex_t IT_0114 = (-2)*IT_0113;
    const ccomplex_t IT_0115 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0116 = IT_0009*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0119 = IT_0003*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = IT_0117 + IT_0120;
    const ccomplex_t IT_0122 = IT_0114 + IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0126 = IT_0021*IT_0125;
    const ccomplex_t IT_0127 = 2*IT_0126;
    const ccomplex_t IT_0128 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0129 = IT_0025*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0132 = IT_0021*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = IT_0130 + IT_0133;
    const ccomplex_t IT_0135 = IT_0127 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0124*IT_0137;
    const ccomplex_t IT_0139 = IT_0111*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = cpow((-2)*s_23 + IT_0037 + IT_0038 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0142 = IT_0138*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0041 + -IT_0071 + IT_0101 + -IT_0104 + 
      -IT_0107 + IT_0110 + IT_0140 + -IT_0143;
    const ccomplex_t IT_0145 = cpow(IT_0054, 2);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_N_3*IT_0145;
    const ccomplex_t IT_0147 = IT_0108*IT_0146;
    const ccomplex_t IT_0148 = IT_0069*IT_0146;
    const ccomplex_t IT_0149 = cpow(IT_0019, 2);
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*m_N_2*IT_0149;
    const ccomplex_t IT_0151 = IT_0105*IT_0150;
    const ccomplex_t IT_0152 = cpow(IT_0124, 2);
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*m_N_4*IT_0152;
    const ccomplex_t IT_0154 = IT_0111*IT_0153;
    const ccomplex_t IT_0155 = IT_0141*IT_0153;
    const ccomplex_t IT_0156 = m_C_2*IT_0110;
    const ccomplex_t IT_0157 = m_C_2*IT_0071;
    const ccomplex_t IT_0158 = m_C_2*IT_0101;
    const ccomplex_t IT_0159 = m_C_2*IT_0104;
    const ccomplex_t IT_0160 = m_C_2*IT_0041;
    const ccomplex_t IT_0161 = m_C_2*IT_0107;
    const ccomplex_t IT_0162 = cpow(IT_0098, 2);
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_N_1*IT_0162;
    const ccomplex_t IT_0164 = IT_0102*IT_0163;
    const ccomplex_t IT_0165 = IT_0039*IT_0150;
    const ccomplex_t IT_0166 = IT_0072*IT_0163;
    const ccomplex_t IT_0167 = m_C_2*IT_0140;
    const ccomplex_t IT_0168 = m_C_2*IT_0143;
    const ccomplex_t IT_0169 = IT_0147 + IT_0148 + IT_0151 + IT_0154 + IT_0155
       + -IT_0156 + -IT_0157 + -IT_0158 + -IT_0159 + -IT_0160 + -IT_0161 +
       IT_0164 + IT_0165 + IT_0166 + -IT_0167 + -IT_0168;
    const ccomplex_t IT_0170 = s_24*m_C_2;
    const ccomplex_t IT_0171 = cpow(IT_0035, 2);
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*m_N_2*IT_0171;
    const ccomplex_t IT_0173 = IT_0039*IT_0172;
    const ccomplex_t IT_0174 = IT_0105*IT_0172;
    const ccomplex_t IT_0175 = cpow(IT_0067, 2);
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*m_N_3*IT_0175;
    const ccomplex_t IT_0177 = IT_0108*IT_0176;
    const ccomplex_t IT_0178 = IT_0069*IT_0176;
    const ccomplex_t IT_0179 = cpow(IT_0137, 2);
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_N_4*IT_0179;
    const ccomplex_t IT_0181 = IT_0111*IT_0180;
    const ccomplex_t IT_0182 = IT_0141*IT_0180;
    const ccomplex_t IT_0183 = cpow(IT_0085, 2);
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*m_N_1*IT_0183;
    const ccomplex_t IT_0185 = IT_0072*IT_0184;
    const ccomplex_t IT_0186 = IT_0102*IT_0184;
    const ccomplex_t IT_0187 = -IT_0156 + -IT_0157 + -IT_0158 + -IT_0159 + 
      -IT_0160 + -IT_0161 + -IT_0167 + -IT_0168 + IT_0173 + IT_0174 + IT_0177 +
       IT_0178 + IT_0181 + IT_0182 + IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = s_14*m_C_2;
    const ccomplex_t IT_0189 = IT_0037*IT_0038;
    const ccomplex_t IT_0190 = s_12*IT_0037;
    const ccomplex_t IT_0191 = s_14*s_24;
    const ccomplex_t IT_0192 = (-2)*IT_0191;
    const ccomplex_t IT_0193 = IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = 2*s_12;
    const ccomplex_t IT_0195 = (-2)*conj(IT_0187);
    const ccomplex_t IT_0196 = (-2)*IT_0144;
    const ccomplex_t IT_0197 = (-2)*IT_0144*(conj(IT_0169)*IT_0170 + -conj
      (IT_0187)*IT_0188 + conj(IT_0144)*IT_0189 + conj(IT_0144)*IT_0193) + (-2)
      *IT_0187*(IT_0038*conj(IT_0169) + conj(IT_0144)*IT_0170 + -conj(IT_0144)
      *IT_0188 + (-0.5)*conj(IT_0187)*IT_0194) + (-2)*IT_0169*(conj(IT_0144)
      *IT_0170 + -conj(IT_0144)*IT_0188 + (-0.5)*conj(IT_0169)*IT_0194 + (-0.5)
      *IT_0038*IT_0195) + 2*IT_0144*(conj(IT_0169)*IT_0188 + 0.5*IT_0170*IT_0195
      ) + (conj(IT_0144)*IT_0189 + conj(IT_0144)*IT_0193)*IT_0196;
    return create_ccomplex_return(IT_0197);
}

