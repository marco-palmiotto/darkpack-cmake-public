#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_sb_2_to_u_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_sb_2_to_u_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_u = param->m_u;
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
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_su_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = m_b*N_d3*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0014 + 1.5*IT_0020);
    const ccomplex_t IT_0022 = 0.666666666666667*IT_0021;
    const ccomplex_t IT_0023 = 1.5*IT_0022;
    const ccomplex_t IT_0024 = 0.666666666666667*IT_0023;
    const ccomplex_t IT_0025 = 1.5*IT_0024;
    const ccomplex_t IT_0026 = 0.666666666666667*IT_0025;
    const ccomplex_t IT_0027 = 1.5*IT_0026;
    const ccomplex_t IT_0028 = (-0.333333333333333)*IT_0027;
    const ccomplex_t IT_0029 = IT_0009*IT_0028;
    const ccomplex_t IT_0030 = IT_0002*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0036 = IT_0011*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = m_b*N_d2*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0039 = IT_0017*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + 1.5*IT_0040);
    const ccomplex_t IT_0042 = 0.666666666666667*IT_0041;
    const ccomplex_t IT_0043 = 1.5*IT_0042;
    const ccomplex_t IT_0044 = (-0.333333333333333)*IT_0043;
    const ccomplex_t IT_0045 = IT_0034*IT_0044;
    const ccomplex_t IT_0046 = IT_0032*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0052 = IT_0011*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = m_b*N_d4*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0055 = IT_0017*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 1.5*IT_0056);
    const ccomplex_t IT_0058 = 0.666666666666667*IT_0057;
    const ccomplex_t IT_0059 = 1.5*IT_0058;
    const ccomplex_t IT_0060 = (-0.333333333333333)*IT_0059;
    const ccomplex_t IT_0061 = IT_0050*IT_0060;
    const ccomplex_t IT_0062 = IT_0048*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0068 = IT_0011*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_b*N_d1*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0071 = IT_0017*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 1.5*IT_0072);
    const ccomplex_t IT_0074 = (-0.333333333333333)*IT_0073;
    const ccomplex_t IT_0075 = IT_0066*IT_0074;
    const ccomplex_t IT_0076 = IT_0064*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0031 + -IT_0047 + -IT_0063 + -IT_0077;
    const ccomplex_t IT_0079 = pow(m_sG, 2);
    const ccomplex_t IT_0080 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0079 +
       reg_prop, -1);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = IT_0080*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0011;
    const ccomplex_t IT_0088 = 0.666666666666667*IT_0087;
    const ccomplex_t IT_0089 = m_b*conj(N_d1)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0090 = IT_0017*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0093 = IT_0007*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0096 = IT_0011*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0091 + -IT_0094 +
       0.333333333333333*IT_0097);
    const ccomplex_t IT_0099 = 3*IT_0098;
    const ccomplex_t IT_0100 = 0.166666666666667*IT_0099;
    const ccomplex_t IT_0101 = IT_0088*IT_0100;
    const ccomplex_t IT_0102 = IT_0064*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0011;
    const ccomplex_t IT_0105 = 0.666666666666667*IT_0104;
    const ccomplex_t IT_0106 = m_b*conj(N_d2)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0107 = IT_0017*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0110 = IT_0007*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0113 = IT_0011*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0108 + -IT_0111 +
       0.333333333333333*IT_0114);
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = (-3)*IT_0116;
    const ccomplex_t IT_0118 = 0.166666666666667*IT_0117;
    const ccomplex_t IT_0119 = IT_0105*IT_0118;
    const ccomplex_t IT_0120 = IT_0032*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0011;
    const ccomplex_t IT_0123 = 0.666666666666667*IT_0122;
    const ccomplex_t IT_0124 = m_b*conj(N_d3)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0125 = IT_0017*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0128 = IT_0007*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0131 = IT_0011*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0126 + -IT_0129 +
       0.333333333333333*IT_0132);
    const ccomplex_t IT_0134 = 3*IT_0133;
    const ccomplex_t IT_0135 = 0.166666666666667*IT_0134;
    const ccomplex_t IT_0136 = IT_0123*IT_0135;
    const ccomplex_t IT_0137 = IT_0002*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0140 = 0.666666666666667*IT_0139;
    const ccomplex_t IT_0141 = m_b*conj(N_d4)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0142 = IT_0017*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0145 = IT_0007*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0148 = IT_0011*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0143 + -IT_0146 +
       0.333333333333333*IT_0149);
    const ccomplex_t IT_0151 = 3*IT_0150;
    const ccomplex_t IT_0152 = 0.166666666666667*IT_0151;
    const ccomplex_t IT_0153 = IT_0140*IT_0152;
    const ccomplex_t IT_0154 = IT_0048*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = 0.166666666666667*IT_0085 + -IT_0103 + -IT_0121
       + -IT_0138 + -IT_0155;
    const ccomplex_t IT_0157 = pow(m_sb_2, 2);
    const ccomplex_t IT_0158 = cpow((-2)*s_23 + IT_0000 + IT_0157 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0159 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0003*IT_0005*IT_0007*IT_0159*V_ub_mod;
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0003*IT_0005*IT_0007*IT_0162*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = IT_0161*IT_0164;
    const ccomplex_t IT_0166 = IT_0158*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = cpow((-2)*s_23 + IT_0000 + IT_0157 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0169 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *e_em*IT_0003*IT_0005*IT_0007*IT_0159*V_ub_mod;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0003*IT_0005*IT_0007*IT_0162*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0172 = 0.5*IT_0171;
    const ccomplex_t IT_0173 = IT_0170*IT_0172;
    const ccomplex_t IT_0174 = IT_0168*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = IT_0167 + IT_0175;
    const ccomplex_t IT_0177 = m_b*m_u*IT_0001;
    const ccomplex_t IT_0178 = s_34*IT_0001;
    const ccomplex_t IT_0179 = (-6)*IT_0178;
    const ccomplex_t IT_0180 = s_13*s_14;
    const ccomplex_t IT_0181 = 12*IT_0180;
    const ccomplex_t IT_0182 = IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = IT_0086 + IT_0176;
    const ccomplex_t IT_0184 = conj(IT_0086) + conj(IT_0176);
    const ccomplex_t IT_0185 = (-18)*IT_0178;
    const ccomplex_t IT_0186 = 36*IT_0180;
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = m_u*IT_0047;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0190 = IT_0074*IT_0088;
    const ccomplex_t IT_0191 = IT_0064*IT_0189*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0193 = IT_0028*IT_0123;
    const ccomplex_t IT_0194 = IT_0002*IT_0192*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0196 = IT_0060*IT_0140;
    const ccomplex_t IT_0197 = IT_0048*IT_0195*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0199 = IT_0081*IT_0198;
    const ccomplex_t IT_0200 = IT_0080*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = m_sG*IT_0201;
    const ccomplex_t IT_0203 = m_u*IT_0077;
    const ccomplex_t IT_0204 = m_u*IT_0031;
    const ccomplex_t IT_0205 = m_u*IT_0063;
    const ccomplex_t IT_0206 = -IT_0188 + IT_0191 + IT_0194 + IT_0197 + (
      -0.166666666666667)*IT_0202 + -IT_0203 + -IT_0204 + -IT_0205;
    const ccomplex_t IT_0207 = 0.5*IT_0202;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0209 = IT_0044*IT_0105;
    const ccomplex_t IT_0210 = IT_0032*IT_0208*IT_0209;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = m_b*IT_0176;
    const ccomplex_t IT_0213 = 6*s_34;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0215 = IT_0007*IT_0016*IT_0162;
    const ccomplex_t IT_0216 = m_b*conj(U_d2)*e_em*IT_0003*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0217 = IT_0215*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = IT_0007*IT_0162;
    const ccomplex_t IT_0220 = e_em*conj(U_Wm2)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0221 = IT_0219*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*(IT_0218 + (-2)*IT_0221);
    const ccomplex_t IT_0223 = (-0.5)*IT_0222;
    const ccomplex_t IT_0224 = IT_0161*IT_0223;
    const ccomplex_t IT_0225 = IT_0158*IT_0214*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0227 = e_em*conj(U_Wm1)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0228 = IT_0219*IT_0227;
    const ccomplex_t IT_0229 = m_b*conj(U_d1)*e_em*IT_0003*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0230 = IT_0215*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0228 + (-0.5)*IT_0231);
    const ccomplex_t IT_0233 = IT_0170*IT_0232;
    const ccomplex_t IT_0234 = IT_0168*IT_0226*IT_0233;
    const ccomplex_t IT_0235 = IT_0225 + IT_0234;
    const ccomplex_t IT_0236 = m_u*IT_0085;
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = IT_0066*IT_0100;
    const ccomplex_t IT_0240 = IT_0064*IT_0189*IT_0239;
    const ccomplex_t IT_0241 = IT_0009*IT_0135;
    const ccomplex_t IT_0242 = IT_0002*IT_0192*IT_0241;
    const ccomplex_t IT_0243 = IT_0050*IT_0152;
    const ccomplex_t IT_0244 = IT_0048*IT_0195*IT_0243;
    const ccomplex_t IT_0245 = m_u*IT_0103;
    const ccomplex_t IT_0246 = m_u*IT_0121;
    const ccomplex_t IT_0247 = m_u*IT_0138;
    const ccomplex_t IT_0248 = m_u*IT_0155;
    const ccomplex_t IT_0249 = 0.166666666666667*IT_0236 + IT_0240 + IT_0242 +
       IT_0244 + -IT_0245 + -IT_0246 + -IT_0247 + -IT_0248;
    const ccomplex_t IT_0250 = IT_0213*conj(IT_0249);
    const ccomplex_t IT_0251 = IT_0034*IT_0118;
    const ccomplex_t IT_0252 = IT_0032*IT_0208*IT_0251;
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = conj(IT_0235) + conj(IT_0237);
    const ccomplex_t IT_0255 = 18*s_34;
    const ccomplex_t IT_0256 = IT_0213*IT_0249;
    const ccomplex_t IT_0257 = IT_0235*IT_0255;
    const ccomplex_t IT_0258 = (-18)*s_34;
    const ccomplex_t IT_0259 = (-6)*s_34;
    const ccomplex_t IT_0260 = m_b*m_u;
    const ccomplex_t IT_0261 = 6*IT_0207;
    const ccomplex_t IT_0262 = (-18)*IT_0211;
    const ccomplex_t IT_0263 = (-0.166666666666667)*IT_0262;
    const ccomplex_t IT_0264 = 18*IT_0206;
    const ccomplex_t IT_0265 = (-0.166666666666667)*IT_0264;
    const ccomplex_t IT_0266 = IT_0212 + (-0.166666666666667)*IT_0261 +
       IT_0263 + IT_0265;
    const ccomplex_t IT_0267 = 18*IT_0211;
    const ccomplex_t IT_0268 = (-6)*IT_0207;
    const ccomplex_t IT_0269 = 0.166666666666667*IT_0268;
    const ccomplex_t IT_0270 = (-18)*IT_0206;
    const ccomplex_t IT_0271 = 0.166666666666667*IT_0270;
    const ccomplex_t IT_0272 = IT_0212 + 0.166666666666667*IT_0267 + IT_0269 +
       IT_0271;
    const ccomplex_t IT_0273 = (-18)*conj(IT_0211);
    const ccomplex_t IT_0274 = (-0.166666666666667)*IT_0273;
    const ccomplex_t IT_0275 = 18*conj(IT_0206);
    const ccomplex_t IT_0276 = (-0.166666666666667)*IT_0275;
    const ccomplex_t IT_0277 = 6*conj(IT_0207);
    const ccomplex_t IT_0278 = (-0.166666666666667)*IT_0277;
    const ccomplex_t IT_0279 = conj(IT_0212) + IT_0274 + IT_0276 + IT_0278;
    const ccomplex_t IT_0280 = (-6)*conj(IT_0207);
    const ccomplex_t IT_0281 = 0.166666666666667*IT_0280;
    const ccomplex_t IT_0282 = (-18)*conj(IT_0206);
    const ccomplex_t IT_0283 = 0.166666666666667*IT_0282;
    const ccomplex_t IT_0284 = 18*conj(IT_0211);
    const ccomplex_t IT_0285 = 0.166666666666667*IT_0284;
    const ccomplex_t IT_0286 = conj(IT_0212) + IT_0281 + IT_0283 + IT_0285;
    const ccomplex_t IT_0287 = 6*IT_0211;
    const ccomplex_t IT_0288 = (-18)*IT_0207;
    const ccomplex_t IT_0289 = (-6)*IT_0206;
    const ccomplex_t IT_0290 = 6*conj(IT_0211);
    const ccomplex_t IT_0291 = (-18)*conj(IT_0207);
    const ccomplex_t IT_0292 = (-6)*conj(IT_0206);
    const ccomplex_t IT_0293 = conj(IT_0253)*IT_0266 + -conj(IT_0249)*IT_0272 
      + IT_0253*IT_0279 + -IT_0249*IT_0286 + (-0.166666666666667)*IT_0254*(18
      *IT_0212 + IT_0287 + IT_0288 + IT_0289) + (-0.166666666666667)*IT_0238*(18
      *conj(IT_0212) + IT_0290 + IT_0291 + IT_0292);
    const ccomplex_t IT_0294 = m_b*s_13;
    const ccomplex_t IT_0295 = 6*IT_0206;
    const ccomplex_t IT_0296 = 18*IT_0207;
    const ccomplex_t IT_0297 = (-6)*IT_0211;
    const ccomplex_t IT_0298 = 6*conj(IT_0206);
    const ccomplex_t IT_0299 = 18*conj(IT_0207);
    const ccomplex_t IT_0300 = (-6)*conj(IT_0211);
    const ccomplex_t IT_0301 = 6*conj(IT_0078)*(IT_0235 + IT_0237 + 3*IT_0249 
      + (-3)*IT_0253) + 6*IT_0078*(conj(IT_0235) + conj(IT_0237) + 3*conj
      (IT_0249) + (-3)*conj(IT_0253)) + (-6)*conj(IT_0156)*IT_0266 + (-6)
      *IT_0156*IT_0279 + IT_0184*((-18)*IT_0212 + IT_0295 + IT_0296 + IT_0297) +
       IT_0183*((-18)*conj(IT_0212) + IT_0298 + IT_0299 + IT_0300);
    const ccomplex_t IT_0302 = m_u*s_14;
    const ccomplex_t IT_0303 = 6*IT_0253;
    const ccomplex_t IT_0304 = (-18)*IT_0235;
    const ccomplex_t IT_0305 = (-18)*IT_0237;
    const ccomplex_t IT_0306 = (-6)*IT_0249;
    const ccomplex_t IT_0307 = 6*conj(IT_0253);
    const ccomplex_t IT_0308 = (-18)*conj(IT_0235);
    const ccomplex_t IT_0309 = (-18)*conj(IT_0237);
    const ccomplex_t IT_0310 = (-6)*conj(IT_0249);
    const ccomplex_t IT_0311 = (-6)*conj(IT_0156)*(IT_0235 + IT_0237 + 3
      *IT_0249 + (-3)*IT_0253) + (-6)*IT_0156*(conj(IT_0235) + conj(IT_0237) + 3
      *conj(IT_0249) + (-3)*conj(IT_0253)) + 6*conj(IT_0078)*IT_0272 + 6*IT_0078
      *IT_0286 + IT_0184*(IT_0303 + IT_0304 + IT_0305 + IT_0306) + IT_0183*
      (IT_0307 + IT_0308 + IT_0309 + IT_0310);
    const ccomplex_t IT_0312 = (-6)*(conj(IT_0078)*(IT_0086 + 3*IT_0156 +
       IT_0176) + IT_0078*(conj(IT_0086) + 3*conj(IT_0156) + conj(IT_0176)))
      *IT_0177 + IT_0182*(conj(IT_0156)*IT_0183 + IT_0156*IT_0184) + (IT_0078
      *conj(IT_0078) + IT_0156*conj(IT_0156) + IT_0183*IT_0184)*IT_0187 + (conj
      (IT_0206)*IT_0207 + IT_0206*conj(IT_0207) + conj(IT_0211)*IT_0212 +
       IT_0211*conj(IT_0212))*IT_0213 + IT_0238*IT_0250 + (IT_0206*conj(IT_0206)
       + IT_0207*conj(IT_0207) + IT_0211*conj(IT_0211) + IT_0212*conj(IT_0212) +
       IT_0249*conj(IT_0249) + IT_0253*conj(IT_0253) + IT_0237*IT_0254)*IT_0255 
      + IT_0254*(IT_0256 + IT_0257) + (conj(IT_0206)*IT_0211 + IT_0206*conj
      (IT_0211) + conj(IT_0207)*IT_0212 + IT_0207*conj(IT_0212) + conj(IT_0249)
      *IT_0253 + IT_0249*conj(IT_0253))*IT_0258 + (conj(IT_0207)*IT_0211 +
       IT_0207*conj(IT_0211) + conj(IT_0206)*IT_0212 + IT_0206*conj(IT_0212) +
       IT_0238*conj(IT_0253) + IT_0253*IT_0254)*IT_0259 + (-6)*IT_0260*IT_0293 +
       IT_0294*IT_0301 + IT_0302*IT_0311;
    return create_ccomplex_return(IT_0312);
}

