#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_sb_1_to_c_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_sb_1_to_c_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sc_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = m_b*N_d3*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0014 + 1.5*IT_0020);
    const ccomplex_t IT_0022 = (-0.333333333333333)*IT_0021;
    const ccomplex_t IT_0023 = IT_0009*IT_0022;
    const ccomplex_t IT_0024 = IT_0002*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0030 = IT_0011*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = m_b*N_d4*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + 1.5*IT_0034);
    const ccomplex_t IT_0036 = (-0.333333333333333)*IT_0035;
    const ccomplex_t IT_0037 = IT_0028*IT_0036;
    const ccomplex_t IT_0038 = IT_0026*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0044 = IT_0011*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = m_b*N_d1*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0047 = IT_0017*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + 1.5*IT_0048);
    const ccomplex_t IT_0050 = (-0.333333333333333)*IT_0049;
    const ccomplex_t IT_0051 = IT_0042*IT_0050;
    const ccomplex_t IT_0052 = IT_0040*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0058 = IT_0011*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = m_b*N_d2*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0061 = IT_0017*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + 1.5*IT_0062);
    const ccomplex_t IT_0064 = (-0.333333333333333)*IT_0063;
    const ccomplex_t IT_0065 = IT_0056*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0025 + -IT_0039 + -IT_0053 + -IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0011;
    const ccomplex_t IT_0070 = 0.666666666666667*IT_0069;
    const ccomplex_t IT_0071 = m_b*conj(N_d2)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0072 = IT_0017*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0075 = IT_0011*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0078 = IT_0007*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0073 + 0.333333333333333
      *IT_0076 + -IT_0079);
    const ccomplex_t IT_0081 = 3*IT_0080;
    const ccomplex_t IT_0082 = 0.166666666666667*IT_0081;
    const ccomplex_t IT_0083 = IT_0070*IT_0082;
    const ccomplex_t IT_0084 = IT_0054*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0087 = 0.666666666666667*IT_0086;
    const ccomplex_t IT_0088 = m_b*conj(N_d4)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0089 = IT_0017*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0092 = IT_0007*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0095 = IT_0011*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0090 + -IT_0093 +
       0.333333333333333*IT_0096);
    const ccomplex_t IT_0098 = 3*IT_0097;
    const ccomplex_t IT_0099 = 0.166666666666667*IT_0098;
    const ccomplex_t IT_0100 = IT_0087*IT_0099;
    const ccomplex_t IT_0101 = IT_0026*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = pow(m_sG, 2);
    const ccomplex_t IT_0104 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0103 +
       reg_prop, -1);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = IT_0104*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0011;
    const ccomplex_t IT_0111 = 0.666666666666667*IT_0110;
    const ccomplex_t IT_0112 = m_b*conj(N_d1)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0113 = IT_0017*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0116 = IT_0007*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0119 = IT_0011*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0114 + -IT_0117 +
       0.333333333333333*IT_0120);
    const ccomplex_t IT_0122 = 3*IT_0121;
    const ccomplex_t IT_0123 = 0.166666666666667*IT_0122;
    const ccomplex_t IT_0124 = IT_0111*IT_0123;
    const ccomplex_t IT_0125 = IT_0040*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0011;
    const ccomplex_t IT_0128 = 0.666666666666667*IT_0127;
    const ccomplex_t IT_0129 = m_b*conj(N_d3)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0130 = IT_0017*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0133 = IT_0007*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0136 = IT_0011*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0131 + -IT_0134 +
       0.333333333333333*IT_0137);
    const ccomplex_t IT_0139 = 3*IT_0138;
    const ccomplex_t IT_0140 = 0.166666666666667*IT_0139;
    const ccomplex_t IT_0141 = IT_0128*IT_0140;
    const ccomplex_t IT_0142 = IT_0002*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = -IT_0085 + -IT_0102 + 0.166666666666667*IT_0109
       + -IT_0126 + -IT_0143;
    const ccomplex_t IT_0145 = (-0.5)*IT_0109;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*conj
      (V_u2)*e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0147 = (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0003*IT_0005*IT_0007*U_sb_00;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = IT_0147*IT_0149;
    const ccomplex_t IT_0151 = pow(m_sb_1, 2);
    const ccomplex_t IT_0152 = cpow((-2)*s_23 + IT_0000 + IT_0151 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0153 = IT_0150*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = cpow((-2)*s_23 + IT_0000 + IT_0151 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*conj
      (V_u1)*e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0003*IT_0005*IT_0007*U_sb_00;
    const ccomplex_t IT_0159 = 0.5*IT_0158;
    const ccomplex_t IT_0160 = IT_0157*IT_0159;
    const ccomplex_t IT_0161 = IT_0155*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = IT_0154 + IT_0162;
    const ccomplex_t IT_0164 = m_b*m_c*IT_0001;
    const ccomplex_t IT_0165 = s_34*IT_0001;
    const ccomplex_t IT_0166 = (-6)*IT_0165;
    const ccomplex_t IT_0167 = s_13*s_14;
    const ccomplex_t IT_0168 = 12*IT_0167;
    const ccomplex_t IT_0169 = IT_0166 + IT_0168;
    const ccomplex_t IT_0170 = IT_0145 + IT_0163;
    const ccomplex_t IT_0171 = conj(IT_0145) + conj(IT_0163);
    const ccomplex_t IT_0172 = (-18)*IT_0165;
    const ccomplex_t IT_0173 = 36*IT_0167;
    const ccomplex_t IT_0174 = IT_0172 + IT_0173;
    const ccomplex_t IT_0175 = m_b*m_c;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0177 = IT_0064*IT_0070;
    const ccomplex_t IT_0178 = IT_0054*IT_0176*IT_0177;
    const ccomplex_t IT_0179 = m_c*IT_0025;
    const ccomplex_t IT_0180 = m_c*IT_0039;
    const ccomplex_t IT_0181 = m_c*IT_0053;
    const ccomplex_t IT_0182 = m_c*IT_0067;
    const ccomplex_t IT_0183 = IT_0178 + -IT_0179 + -IT_0180 + -IT_0181 + 
      -IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0185 = IT_0036*IT_0087;
    const ccomplex_t IT_0186 = IT_0026*IT_0184*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0188 = IT_0050*IT_0111;
    const ccomplex_t IT_0189 = IT_0040*IT_0187*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0191 = IT_0022*IT_0128;
    const ccomplex_t IT_0192 = IT_0002*IT_0190*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0194 = IT_0105*IT_0193;
    const ccomplex_t IT_0195 = IT_0104*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = m_sG*IT_0196;
    const ccomplex_t IT_0198 = -IT_0186 + -IT_0189 + -IT_0192 +
       0.166666666666667*IT_0197;
    const ccomplex_t IT_0199 = m_b*IT_0162;
    const ccomplex_t IT_0200 = m_b*IT_0154;
    const ccomplex_t IT_0201 = (-0.5)*IT_0197 + IT_0199 + IT_0200;
    const ccomplex_t IT_0202 = m_c*IT_0109;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0204 = m_b*conj(U_d1)*V_cb*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0205 = IT_0017*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = V_cb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0208 = IT_0007*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*(IT_0206 + (-2)*IT_0208);
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0157*IT_0210;
    const ccomplex_t IT_0212 = IT_0155*IT_0203*IT_0211;
    const ccomplex_t IT_0213 = (-0.5)*IT_0202 + IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0215 = m_b*conj(U_d2)*V_cb*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0216 = IT_0017*IT_0215;
    const ccomplex_t IT_0217 = 1.4142135623731*IT_0216;
    const ccomplex_t IT_0218 = V_cb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0219 = IT_0007*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*(IT_0217 + (-2)*IT_0219);
    const ccomplex_t IT_0221 = (-0.5)*IT_0220;
    const ccomplex_t IT_0222 = IT_0147*IT_0221;
    const ccomplex_t IT_0223 = IT_0152*IT_0214*IT_0222;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = 6*IT_0183;
    const ccomplex_t IT_0226 = (-18)*IT_0201;
    const ccomplex_t IT_0227 = (-6)*IT_0198;
    const ccomplex_t IT_0228 = IT_0225 + IT_0226 + IT_0227;
    const ccomplex_t IT_0229 = 6*conj(IT_0183);
    const ccomplex_t IT_0230 = (-18)*conj(IT_0201);
    const ccomplex_t IT_0231 = (-6)*conj(IT_0198);
    const ccomplex_t IT_0232 = IT_0229 + IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = IT_0056*IT_0082;
    const ccomplex_t IT_0234 = IT_0054*IT_0176*IT_0233;
    const ccomplex_t IT_0235 = m_c*IT_0085;
    const ccomplex_t IT_0236 = m_c*IT_0126;
    const ccomplex_t IT_0237 = m_c*IT_0143;
    const ccomplex_t IT_0238 = m_c*IT_0102;
    const ccomplex_t IT_0239 = 0.166666666666667*IT_0202 + IT_0234 + -IT_0235 
      + -IT_0236 + -IT_0237 + -IT_0238;
    const ccomplex_t IT_0240 = 6*IT_0201;
    const ccomplex_t IT_0241 = 18*IT_0198;
    const ccomplex_t IT_0242 = (-18)*IT_0183;
    const ccomplex_t IT_0243 = IT_0240 + IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = 6*conj(IT_0201);
    const ccomplex_t IT_0245 = 18*conj(IT_0198);
    const ccomplex_t IT_0246 = (-18)*conj(IT_0183);
    const ccomplex_t IT_0247 = IT_0244 + IT_0245 + IT_0246;
    const ccomplex_t IT_0248 = IT_0042*IT_0123;
    const ccomplex_t IT_0249 = IT_0040*IT_0187*IT_0248;
    const ccomplex_t IT_0250 = IT_0009*IT_0140;
    const ccomplex_t IT_0251 = IT_0002*IT_0190*IT_0250;
    const ccomplex_t IT_0252 = IT_0028*IT_0099;
    const ccomplex_t IT_0253 = IT_0026*IT_0184*IT_0252;
    const ccomplex_t IT_0254 = -IT_0249 + -IT_0251 + -IT_0253;
    const ccomplex_t IT_0255 = 18*IT_0183;
    const ccomplex_t IT_0256 = (-18)*IT_0198;
    const ccomplex_t IT_0257 = (-6)*IT_0201;
    const ccomplex_t IT_0258 = IT_0255 + IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = 18*conj(IT_0183);
    const ccomplex_t IT_0260 = (-18)*conj(IT_0198);
    const ccomplex_t IT_0261 = (-6)*conj(IT_0201);
    const ccomplex_t IT_0262 = IT_0259 + IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = (conj(IT_0183) + -conj(IT_0198) + (-3)*conj
      (IT_0201))*IT_0213 + (IT_0183 + -IT_0198 + (-3)*IT_0201)*conj(IT_0213) + (
      -0.166666666666667)*conj(IT_0224)*IT_0228 + (-0.166666666666667)*IT_0224
      *IT_0232 + (-0.166666666666667)*conj(IT_0239)*IT_0243 + (
      -0.166666666666667)*IT_0239*IT_0247 + (-0.166666666666667)*conj(IT_0254)
      *IT_0258 + (-0.166666666666667)*IT_0254*IT_0262;
    const ccomplex_t IT_0264 = 6*s_34;
    const ccomplex_t IT_0265 = 18*s_34;
    const ccomplex_t IT_0266 = (-18)*s_34;
    const ccomplex_t IT_0267 = (-6)*s_34;
    const ccomplex_t IT_0268 = m_c*s_14;
    const ccomplex_t IT_0269 = -conj(IT_0213);
    const ccomplex_t IT_0270 = (-3)*conj(IT_0239);
    const ccomplex_t IT_0271 = conj(IT_0224) + 3*conj(IT_0254) + IT_0269 +
       IT_0270;
    const ccomplex_t IT_0272 = 6*IT_0254;
    const ccomplex_t IT_0273 = 18*IT_0224;
    const ccomplex_t IT_0274 = (-18)*IT_0213;
    const ccomplex_t IT_0275 = (-6)*IT_0239;
    const ccomplex_t IT_0276 = 6*conj(IT_0254);
    const ccomplex_t IT_0277 = 18*conj(IT_0224);
    const ccomplex_t IT_0278 = (-18)*conj(IT_0213);
    const ccomplex_t IT_0279 = (-6)*conj(IT_0239);
    const ccomplex_t IT_0280 = (-3)*IT_0239;
    const ccomplex_t IT_0281 = -IT_0213;
    const ccomplex_t IT_0282 = 0.166666666666667*conj(IT_0068)*IT_0243 +
       0.166666666666667*IT_0068*IT_0247 + IT_0144*IT_0271 + 0.166666666666667
      *IT_0171*(IT_0272 + IT_0273 + IT_0274 + IT_0275) + 0.166666666666667
      *IT_0170*(IT_0276 + IT_0277 + IT_0278 + IT_0279) + conj(IT_0144)*(IT_0224 
      + 3*IT_0254 + IT_0280 + IT_0281);
    const ccomplex_t IT_0283 = m_b*s_13;
    const ccomplex_t IT_0284 = (-6)*((-0.166666666666667)*IT_0171*IT_0228 + (
      -0.166666666666667)*IT_0170*IT_0232 + (-0.166666666666667)*conj(IT_0144)
      *IT_0258 + (-0.166666666666667)*IT_0144*IT_0262 + IT_0068*IT_0271 + conj
      (IT_0068)*(IT_0224 + 3*IT_0254 + IT_0280 + IT_0281))*IT_0283;
    const ccomplex_t IT_0285 = (conj(IT_0068)*((-18)*IT_0144 + (-6)*IT_0145 + 
      (-6)*IT_0163) + IT_0068*((-18)*conj(IT_0144) + (-6)*conj(IT_0145) + (-6)
      *conj(IT_0163)))*IT_0164 + IT_0169*(conj(IT_0144)*IT_0170 + IT_0144
      *IT_0171) + (IT_0068*conj(IT_0068) + IT_0144*conj(IT_0144) + IT_0170
      *IT_0171)*IT_0174 + (-6)*IT_0175*IT_0263 + (conj(IT_0198)*IT_0201 +
       IT_0198*conj(IT_0201) + conj(IT_0213)*IT_0239 + IT_0213*conj(IT_0239) +
       conj(IT_0224)*IT_0254 + IT_0224*conj(IT_0254))*IT_0264 + (IT_0183*conj
      (IT_0183) + IT_0198*conj(IT_0198) + IT_0201*conj(IT_0201) + IT_0213*conj
      (IT_0213) + IT_0224*conj(IT_0224) + IT_0239*conj(IT_0239) + IT_0254*conj
      (IT_0254))*IT_0265 + (conj(IT_0183)*IT_0198 + IT_0183*conj(IT_0198) + conj
      (IT_0213)*IT_0224 + IT_0213*conj(IT_0224) + conj(IT_0239)*IT_0254 +
       IT_0239*conj(IT_0254))*IT_0266 + (conj(IT_0183)*IT_0201 + IT_0183*conj
      (IT_0201) + conj(IT_0224)*IT_0239 + IT_0224*conj(IT_0239) + conj(IT_0213)
      *IT_0254 + IT_0213*conj(IT_0254))*IT_0267 + 6*IT_0268*IT_0282 + IT_0284;
    return create_ccomplex_return(IT_0285);
}

