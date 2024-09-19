#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_sb_1_to_u_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_sb_1_to_u_b(
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
    const creal_t m_sb_1 = param->m_sb_1;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_su_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = m_b*N_d1*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0014 + 1.5*IT_0020);
    const ccomplex_t IT_0022 = (-0.333333333333333)*IT_0021;
    const ccomplex_t IT_0023 = IT_0009*IT_0022;
    const ccomplex_t IT_0024 = IT_0002*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0030 = IT_0011*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = m_b*N_d2*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + 1.5*IT_0034);
    const ccomplex_t IT_0036 = (-0.333333333333333)*IT_0035;
    const ccomplex_t IT_0037 = IT_0028*IT_0036;
    const ccomplex_t IT_0038 = IT_0026*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0044 = IT_0011*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = m_b*N_d3*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0047 = IT_0017*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + 1.5*IT_0048);
    const ccomplex_t IT_0050 = (-0.333333333333333)*IT_0049;
    const ccomplex_t IT_0051 = IT_0042*IT_0050;
    const ccomplex_t IT_0052 = IT_0040*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0058 = IT_0011*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = m_b*N_d4*e_em*IT_0003*U_sb_00;
    const ccomplex_t IT_0061 = IT_0017*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + 1.5*IT_0062);
    const ccomplex_t IT_0064 = 0.666666666666667*IT_0063;
    const ccomplex_t IT_0065 = 1.5*IT_0064;
    const ccomplex_t IT_0066 = (-0.333333333333333)*IT_0065;
    const ccomplex_t IT_0067 = IT_0056*IT_0066;
    const ccomplex_t IT_0068 = IT_0054*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = -IT_0025 + -IT_0039 + -IT_0053 + -IT_0069;
    const ccomplex_t IT_0071 = pow(m_sG, 2);
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0071 +
       reg_prop, -1);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = IT_0072*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0011;
    const ccomplex_t IT_0079 = 0.666666666666667*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0081 = IT_0011*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = m_b*conj(N_d1)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0084 = IT_0017*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0087 = IT_0007*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0082 + 3*IT_0085 + (-3)
      *IT_0088);
    const ccomplex_t IT_0090 = (-0.333333333333333)*IT_0089;
    const ccomplex_t IT_0091 = (-3)*IT_0090;
    const ccomplex_t IT_0092 = (-0.333333333333333)*IT_0091;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = 3*IT_0093;
    const ccomplex_t IT_0095 = 0.166666666666667*IT_0094;
    const ccomplex_t IT_0096 = IT_0079*IT_0095;
    const ccomplex_t IT_0097 = IT_0002*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0011;
    const ccomplex_t IT_0100 = 0.666666666666667*IT_0099;
    const ccomplex_t IT_0101 = m_b*conj(N_d2)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0102 = IT_0017*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0105 = IT_0007*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0108 = IT_0011*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + -IT_0106 +
       0.333333333333333*IT_0109);
    const ccomplex_t IT_0111 = 3*IT_0110;
    const ccomplex_t IT_0112 = 0.166666666666667*IT_0111;
    const ccomplex_t IT_0113 = IT_0100*IT_0112;
    const ccomplex_t IT_0114 = IT_0026*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0011;
    const ccomplex_t IT_0117 = 0.666666666666667*IT_0116;
    const ccomplex_t IT_0118 = m_b*conj(N_d3)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0119 = IT_0017*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0122 = IT_0011*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0125 = IT_0007*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0120 + 0.333333333333333
      *IT_0123 + -IT_0126);
    const ccomplex_t IT_0128 = 3*IT_0127;
    const ccomplex_t IT_0129 = 0.166666666666667*IT_0128;
    const ccomplex_t IT_0130 = IT_0117*IT_0129;
    const ccomplex_t IT_0131 = IT_0040*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0134 = 0.666666666666667*IT_0133;
    const ccomplex_t IT_0135 = m_b*conj(N_d4)*e_em*IT_0003*U_sb_10;
    const ccomplex_t IT_0136 = IT_0017*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0139 = IT_0011*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0142 = IT_0007*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0137 + 0.333333333333333
      *IT_0140 + -IT_0143);
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = 3*IT_0146;
    const ccomplex_t IT_0148 = 0.166666666666667*IT_0147;
    const ccomplex_t IT_0149 = IT_0134*IT_0148;
    const ccomplex_t IT_0150 = IT_0054*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = 0.166666666666667*IT_0077 + -IT_0098 + -IT_0115
       + -IT_0132 + -IT_0151;
    const ccomplex_t IT_0153 = (-0.5)*IT_0077;
    const ccomplex_t IT_0154 = pow(m_sb_1, 2);
    const ccomplex_t IT_0155 = cpow((-2)*s_23 + IT_0000 + IT_0154 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0156 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0157 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *e_em*IT_0003*IT_0005*IT_0007*IT_0156*V_ub_mod;
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0003*IT_0005*IT_0007*IT_0159*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = IT_0158*IT_0161;
    const ccomplex_t IT_0163 = IT_0155*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = cpow((-2)*s_23 + IT_0000 + IT_0154 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0166 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0003*IT_0005*IT_0007*IT_0156*V_ub_mod;
    const ccomplex_t IT_0167 = (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0003*IT_0005*IT_0007*IT_0159*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = IT_0167*IT_0169;
    const ccomplex_t IT_0171 = IT_0165*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = IT_0164 + IT_0172;
    const ccomplex_t IT_0174 = m_b*m_u*IT_0001;
    const ccomplex_t IT_0175 = s_34*IT_0001;
    const ccomplex_t IT_0176 = (-6)*IT_0175;
    const ccomplex_t IT_0177 = s_13*s_14;
    const ccomplex_t IT_0178 = 12*IT_0177;
    const ccomplex_t IT_0179 = IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = IT_0153 + IT_0173;
    const ccomplex_t IT_0181 = conj(IT_0153) + conj(IT_0173);
    const ccomplex_t IT_0182 = (-18)*IT_0175;
    const ccomplex_t IT_0183 = 36*IT_0177;
    const ccomplex_t IT_0184 = IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0186 = IT_0022*IT_0079;
    const ccomplex_t IT_0187 = IT_0002*IT_0185*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0189 = IT_0073*IT_0188;
    const ccomplex_t IT_0190 = IT_0072*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = m_sG*IT_0191;
    const ccomplex_t IT_0193 = -IT_0187 + 0.166666666666667*IT_0192;
    const ccomplex_t IT_0194 = (-0.5)*IT_0192;
    const ccomplex_t IT_0195 = m_b*IT_0173;
    const ccomplex_t IT_0196 = IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = conj(IT_0194) + conj(IT_0195);
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0199 = IT_0007*IT_0159;
    const ccomplex_t IT_0200 = e_em*conj(U_Wm2)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0201 = IT_0199*IT_0200;
    const ccomplex_t IT_0202 = IT_0007*IT_0016*IT_0159;
    const ccomplex_t IT_0203 = m_b*conj(U_d2)*e_em*IT_0003*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0204 = IT_0202*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*(IT_0201 + (-0.5)*IT_0205);
    const ccomplex_t IT_0207 = IT_0167*IT_0206;
    const ccomplex_t IT_0208 = IT_0165*IT_0198*IT_0207;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*(0.333333333333333*IT_0082 +
       IT_0085 + -IT_0088);
    const ccomplex_t IT_0211 = 3*IT_0210;
    const ccomplex_t IT_0212 = (-0.333333333333333)*IT_0211;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = 3*IT_0213;
    const ccomplex_t IT_0215 = 0.166666666666667*IT_0214;
    const ccomplex_t IT_0216 = IT_0009*IT_0215;
    const ccomplex_t IT_0217 = IT_0002*IT_0185*IT_0216;
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0220 = IT_0042*IT_0129;
    const ccomplex_t IT_0221 = IT_0040*IT_0219*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0223 = -IT_0146;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = 3*IT_0224;
    const ccomplex_t IT_0226 = 0.166666666666667*IT_0225;
    const ccomplex_t IT_0227 = IT_0056*IT_0226;
    const ccomplex_t IT_0228 = IT_0054*IT_0222*IT_0227;
    const ccomplex_t IT_0229 = -IT_0221 + -IT_0228;
    const ccomplex_t IT_0230 = IT_0218 + IT_0229;
    const ccomplex_t IT_0231 = conj(IT_0218) + conj(IT_0229);
    const ccomplex_t IT_0232 = 6*s_34;
    const ccomplex_t IT_0233 = m_b*s_13;
    const ccomplex_t IT_0234 = IT_0098 + IT_0132 + IT_0151;
    const ccomplex_t IT_0235 = m_u*IT_0234;
    const ccomplex_t IT_0236 = -IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0238 = IT_0028*IT_0112;
    const ccomplex_t IT_0239 = IT_0026*IT_0237*IT_0238;
    const ccomplex_t IT_0240 = m_u*IT_0077;
    const ccomplex_t IT_0241 = m_u*IT_0115;
    const ccomplex_t IT_0242 = IT_0239 + 0.166666666666667*IT_0240 + -IT_0241;
    const ccomplex_t IT_0243 = (-0.5)*IT_0240;
    const ccomplex_t IT_0244 = m_u*IT_0069;
    const ccomplex_t IT_0245 = IT_0036*IT_0100;
    const ccomplex_t IT_0246 = IT_0026*IT_0237*IT_0245;
    const ccomplex_t IT_0247 = m_u*IT_0039;
    const ccomplex_t IT_0248 = m_u*IT_0053;
    const ccomplex_t IT_0249 = -IT_0244 + IT_0246 + -IT_0247 + -IT_0248;
    const ccomplex_t IT_0250 = 6*IT_0249;
    const ccomplex_t IT_0251 = m_u*IT_0025;
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = 6*IT_0252;
    const ccomplex_t IT_0254 = (-18)*IT_0194;
    const ccomplex_t IT_0255 = (-6)*IT_0193;
    const ccomplex_t IT_0256 = IT_0050*IT_0117;
    const ccomplex_t IT_0257 = IT_0040*IT_0219*IT_0256;
    const ccomplex_t IT_0258 = IT_0066*IT_0134;
    const ccomplex_t IT_0259 = IT_0054*IT_0222*IT_0258;
    const ccomplex_t IT_0260 = -IT_0257 + -IT_0259;
    const ccomplex_t IT_0261 = (-6)*IT_0260;
    const ccomplex_t IT_0262 = IT_0250 + IT_0253 + IT_0254 + IT_0255 + IT_0261;
    const ccomplex_t IT_0263 = 6*conj(IT_0249);
    const ccomplex_t IT_0264 = 6*conj(IT_0252);
    const ccomplex_t IT_0265 = (-18)*conj(IT_0194);
    const ccomplex_t IT_0266 = (-6)*conj(IT_0193);
    const ccomplex_t IT_0267 = (-6)*conj(IT_0260);
    const ccomplex_t IT_0268 = IT_0263 + IT_0264 + IT_0265 + IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = 18*IT_0249;
    const ccomplex_t IT_0270 = 18*IT_0252;
    const ccomplex_t IT_0271 = (-18)*IT_0193;
    const ccomplex_t IT_0272 = (-18)*IT_0260;
    const ccomplex_t IT_0273 = (-6)*IT_0194;
    const ccomplex_t IT_0274 = IT_0269 + IT_0270 + IT_0271 + IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = 18*conj(IT_0249);
    const ccomplex_t IT_0276 = 18*conj(IT_0252);
    const ccomplex_t IT_0277 = (-18)*conj(IT_0193);
    const ccomplex_t IT_0278 = (-18)*conj(IT_0260);
    const ccomplex_t IT_0279 = (-6)*conj(IT_0194);
    const ccomplex_t IT_0280 = IT_0275 + IT_0276 + IT_0277 + IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = (-18)*IT_0153;
    const ccomplex_t IT_0282 = (-18)*IT_0173;
    const ccomplex_t IT_0283 = (-6)*IT_0152;
    const ccomplex_t IT_0284 = IT_0281 + IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = (-18)*conj(IT_0153);
    const ccomplex_t IT_0286 = (-18)*conj(IT_0173);
    const ccomplex_t IT_0287 = (-6)*conj(IT_0152);
    const ccomplex_t IT_0288 = IT_0285 + IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = e_em*conj(U_Wm1)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0290 = IT_0199*IT_0289;
    const ccomplex_t IT_0291 = m_b*conj(U_d1)*e_em*IT_0003*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0292 = IT_0202*IT_0291;
    const ccomplex_t IT_0293 = 1.4142135623731*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*(IT_0290 + (-0.5)*IT_0293);
    const ccomplex_t IT_0295 = IT_0158*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0297 = IT_0155*IT_0295*IT_0296;
    const ccomplex_t IT_0298 = 6*IT_0070;
    const ccomplex_t IT_0299 = 6*conj(IT_0070);
    const ccomplex_t IT_0300 = (-6)*conj(IT_0070)*(IT_0209 + 3*IT_0218 + 3
      *IT_0229 + (-3)*IT_0236 + (-3)*IT_0242 + -IT_0243) + (-6)*IT_0070*(conj
      (IT_0209) + 3*conj(IT_0218) + 3*conj(IT_0229) + (-3)*conj(IT_0236) + (-3)
      *conj(IT_0242) + -conj(IT_0243)) + IT_0181*IT_0262 + IT_0180*IT_0268 +
       conj(IT_0152)*IT_0274 + IT_0152*IT_0280 + conj(IT_0195)*IT_0284 + IT_0195
      *IT_0288 + conj(IT_0297)*IT_0298 + IT_0297*IT_0299;
    const ccomplex_t IT_0301 = IT_0249 + IT_0252;
    const ccomplex_t IT_0302 = conj(IT_0249) + conj(IT_0252);
    const ccomplex_t IT_0303 = conj(IT_0236) + conj(IT_0242);
    const ccomplex_t IT_0304 = conj(IT_0243) + conj(IT_0297);
    const ccomplex_t IT_0305 = 18*s_34;
    const ccomplex_t IT_0306 = IT_0236*IT_0305;
    const ccomplex_t IT_0307 = IT_0243 + IT_0297;
    const ccomplex_t IT_0308 = (-18)*s_34;
    const ccomplex_t IT_0309 = conj(IT_0193)*IT_0308;
    const ccomplex_t IT_0310 = IT_0218*IT_0305;
    const ccomplex_t IT_0311 = IT_0236*IT_0308;
    const ccomplex_t IT_0312 = IT_0242*IT_0308;
    const ccomplex_t IT_0313 = conj(IT_0236)*IT_0308;
    const ccomplex_t IT_0314 = conj(IT_0242)*IT_0308;
    const ccomplex_t IT_0315 = IT_0236 + IT_0242;
    const ccomplex_t IT_0316 = (-6)*s_34;
    const ccomplex_t IT_0317 = IT_0232*IT_0260;
    const ccomplex_t IT_0318 = IT_0194*IT_0305;
    const ccomplex_t IT_0319 = IT_0252*IT_0316;
    const ccomplex_t IT_0320 = IT_0232*conj(IT_0260);
    const ccomplex_t IT_0321 = conj(IT_0249)*IT_0316;
    const ccomplex_t IT_0322 = conj(IT_0252)*IT_0316;
    const ccomplex_t IT_0323 = IT_0232*IT_0236;
    const ccomplex_t IT_0324 = IT_0232*IT_0242;
    const ccomplex_t IT_0325 = IT_0243*IT_0305;
    const ccomplex_t IT_0326 = IT_0218*IT_0316;
    const ccomplex_t IT_0327 = IT_0232*conj(IT_0236);
    const ccomplex_t IT_0328 = IT_0232*conj(IT_0242);
    const ccomplex_t IT_0329 = conj(IT_0218)*IT_0316;
    const ccomplex_t IT_0330 = m_b*m_u;
    const ccomplex_t IT_0331 = 6*IT_0194;
    const ccomplex_t IT_0332 = 18*IT_0193;
    const ccomplex_t IT_0333 = 18*IT_0260;
    const ccomplex_t IT_0334 = (-18)*IT_0249;
    const ccomplex_t IT_0335 = (-18)*IT_0252;
    const ccomplex_t IT_0336 = IT_0331 + IT_0332 + IT_0333 + IT_0334 + IT_0335;
    const ccomplex_t IT_0337 = 6*conj(IT_0194);
    const ccomplex_t IT_0338 = 18*conj(IT_0193);
    const ccomplex_t IT_0339 = 18*conj(IT_0260);
    const ccomplex_t IT_0340 = (-18)*conj(IT_0249);
    const ccomplex_t IT_0341 = (-18)*conj(IT_0252);
    const ccomplex_t IT_0342 = IT_0337 + IT_0338 + IT_0339 + IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = 6*IT_0193;
    const ccomplex_t IT_0344 = 6*IT_0260;
    const ccomplex_t IT_0345 = 18*IT_0194;
    const ccomplex_t IT_0346 = (-6)*IT_0249;
    const ccomplex_t IT_0347 = (-6)*IT_0252;
    const ccomplex_t IT_0348 = 6*conj(IT_0193);
    const ccomplex_t IT_0349 = 6*conj(IT_0260);
    const ccomplex_t IT_0350 = 18*conj(IT_0194);
    const ccomplex_t IT_0351 = (-6)*conj(IT_0249);
    const ccomplex_t IT_0352 = (-6)*conj(IT_0252);
    const ccomplex_t IT_0353 = (-0.166666666666667)*conj(IT_0209)*IT_0262 + (
      -0.166666666666667)*IT_0209*IT_0268 + (-0.166666666666667)*IT_0231*IT_0274
       + (-0.166666666666667)*IT_0230*IT_0280 + 3*IT_0195*(conj(IT_0209) + 1./3
      *conj(IT_0218) + 1./3*conj(IT_0229) + -conj(IT_0243) + -conj(IT_0297) + 
      -1./3*IT_0303) + 3*conj(IT_0195)*(IT_0209 + 1./3*IT_0218 + 1./3*IT_0229 + 
      -IT_0243 + -IT_0297 + -1./3*IT_0315) + (-0.166666666666667)*IT_0303
      *IT_0336 + (-0.166666666666667)*IT_0315*IT_0342 + (-0.166666666666667)
      *IT_0304*(IT_0343 + IT_0344 + IT_0345 + IT_0346 + IT_0347) + (
      -0.166666666666667)*IT_0307*(IT_0348 + IT_0349 + IT_0350 + IT_0351 +
       IT_0352);
    const ccomplex_t IT_0354 = m_u*s_14;
    const ccomplex_t IT_0355 = 6*IT_0218;
    const ccomplex_t IT_0356 = 6*IT_0229;
    const ccomplex_t IT_0357 = 18*IT_0209;
    const ccomplex_t IT_0358 = (-18)*IT_0243;
    const ccomplex_t IT_0359 = (-6)*IT_0236;
    const ccomplex_t IT_0360 = (-6)*IT_0242;
    const ccomplex_t IT_0361 = 6*conj(IT_0218);
    const ccomplex_t IT_0362 = 6*conj(IT_0229);
    const ccomplex_t IT_0363 = 18*conj(IT_0209);
    const ccomplex_t IT_0364 = (-18)*conj(IT_0243);
    const ccomplex_t IT_0365 = (-6)*conj(IT_0236);
    const ccomplex_t IT_0366 = (-6)*conj(IT_0242);
    const ccomplex_t IT_0367 = 6*conj(IT_0152)*(IT_0209 + 3*IT_0218 + 3
      *IT_0229 + (-3)*IT_0236 + (-3)*IT_0242 + -IT_0243) + 6*IT_0152*(conj
      (IT_0209) + 3*conj(IT_0218) + 3*conj(IT_0229) + (-3)*conj(IT_0236) + (-3)
      *conj(IT_0242) + -conj(IT_0243)) + IT_0288*IT_0297 + IT_0284*conj(IT_0297)
       + conj(IT_0195)*IT_0298 + IT_0195*IT_0299 + conj(IT_0070)*IT_0336 +
       IT_0070*IT_0342 + IT_0181*(IT_0355 + IT_0356 + IT_0357 + IT_0358 +
       IT_0359 + IT_0360) + IT_0180*(IT_0361 + IT_0362 + IT_0363 + IT_0364 +
       IT_0365 + IT_0366);
    const ccomplex_t IT_0368 = (conj(IT_0070)*((-18)*IT_0152 + (-6)*IT_0153 + 
      (-6)*IT_0173) + IT_0070*((-18)*conj(IT_0152) + (-6)*conj(IT_0153) + (-6)
      *conj(IT_0173)))*IT_0174 + IT_0179*(conj(IT_0152)*IT_0180 + IT_0152
      *IT_0181) + (IT_0070*conj(IT_0070) + IT_0152*conj(IT_0152) + IT_0180
      *IT_0181)*IT_0184 + (conj(IT_0193)*IT_0196 + IT_0193*IT_0197 + conj
      (IT_0209)*IT_0230 + IT_0209*IT_0231)*IT_0232 + IT_0233*IT_0300 + (IT_0195
      *IT_0197 + IT_0209*conj(IT_0209) + IT_0229*IT_0231 + (IT_0193 + IT_0260)*
      (conj(IT_0193) + conj(IT_0260)) + IT_0301*IT_0302 + IT_0242*IT_0303 +
       IT_0297*IT_0304)*IT_0305 + IT_0303*IT_0306 + (conj(IT_0260)*IT_0301 + 
      (IT_0193 + IT_0260)*IT_0302 + IT_0209*IT_0304 + conj(IT_0209)*IT_0307)
      *IT_0308 + IT_0301*IT_0309 + IT_0231*(IT_0310 + IT_0311 + IT_0312) +
       IT_0230*(IT_0313 + IT_0314) + (IT_0197*IT_0249 + IT_0209*IT_0303 +
       IT_0229*IT_0304 + conj(IT_0229)*IT_0307 + conj(IT_0209)*IT_0315)*IT_0316 
      + IT_0197*(IT_0317 + IT_0318 + IT_0319) + IT_0196*(IT_0320 + IT_0321 +
       IT_0322) + IT_0304*(IT_0323 + IT_0324 + IT_0325 + IT_0326) + IT_0307*
      (IT_0327 + IT_0328 + IT_0329) + (-6)*IT_0330*IT_0353 + IT_0354*IT_0367;
    return create_ccomplex_return(IT_0368);
}

