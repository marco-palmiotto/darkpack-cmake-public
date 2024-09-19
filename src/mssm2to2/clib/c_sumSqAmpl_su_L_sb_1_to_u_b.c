#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_sb_1_to_u_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_sb_1_to_u_b(
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
    const creal_t m_su_L = param->m_su_L;
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
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_su_L, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_13*s_14;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_u, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = m_b*conj(N_d4)*e_em*IT_0008*U_sb_10;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0022 = IT_0012*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = cos(theta_W);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0020 + -IT_0023 +
       0.333333333333333*IT_0028);
    const ccomplex_t IT_0030 = 3*IT_0029;
    const ccomplex_t IT_0031 = 0.333333333333333*IT_0030;
    const ccomplex_t IT_0032 = 3*IT_0031;
    const ccomplex_t IT_0033 = 0.166666666666667*IT_0032;
    const ccomplex_t IT_0034 = IT_0014*IT_0033;
    const ccomplex_t IT_0035 = IT_0007*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = m_b*conj(N_d1)*e_em*IT_0008*U_sb_10;
    const ccomplex_t IT_0041 = IT_0017*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0044 = IT_0025*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0047 = IT_0012*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0042 + 0.333333333333333
      *IT_0045 + -IT_0048);
    const ccomplex_t IT_0050 = 3*IT_0049;
    const ccomplex_t IT_0051 = (-0.333333333333333)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = 3*IT_0052;
    const ccomplex_t IT_0054 = (-0.333333333333333)*IT_0053;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = (-3)*IT_0056;
    const ccomplex_t IT_0058 = (-0.333333333333333)*IT_0057;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = 3*IT_0059;
    const ccomplex_t IT_0061 = 0.166666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0039*IT_0061;
    const ccomplex_t IT_0063 = IT_0037*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = m_b*conj(N_d2)*e_em*IT_0008*U_sb_10;
    const ccomplex_t IT_0069 = IT_0017*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0072 = IT_0025*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0075 = IT_0012*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0070 + 0.333333333333333
      *IT_0073 + -IT_0076);
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = 3*IT_0079;
    const ccomplex_t IT_0081 = 0.333333333333333*IT_0080;
    const ccomplex_t IT_0082 = 3*IT_0081;
    const ccomplex_t IT_0083 = 0.333333333333333*IT_0082;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = 3*IT_0087;
    const ccomplex_t IT_0089 = 0.333333333333333*IT_0088;
    const ccomplex_t IT_0090 = 3*IT_0089;
    const ccomplex_t IT_0091 = (-0.333333333333333)*IT_0090;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = 3*IT_0092;
    const ccomplex_t IT_0094 = 0.166666666666667*IT_0093;
    const ccomplex_t IT_0095 = IT_0067*IT_0094;
    const ccomplex_t IT_0096 = IT_0065*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0008*IT_0010*IT_0012;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = m_b*conj(N_d3)*e_em*IT_0008*U_sb_10;
    const ccomplex_t IT_0102 = IT_0017*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0105 = IT_0025*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0108 = IT_0012*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + 0.333333333333333
      *IT_0106 + -IT_0109);
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = 3*IT_0114;
    const ccomplex_t IT_0116 = (-0.333333333333333)*IT_0115;
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = 3*IT_0117;
    const ccomplex_t IT_0119 = (-0.333333333333333)*IT_0118;
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = -IT_0121;
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = 3*IT_0124;
    const ccomplex_t IT_0126 = (-0.333333333333333)*IT_0125;
    const ccomplex_t IT_0127 = (-3)*IT_0126;
    const ccomplex_t IT_0128 = 0.333333333333333*IT_0127;
    const ccomplex_t IT_0129 = 3*IT_0128;
    const ccomplex_t IT_0130 = 0.166666666666667*IT_0129;
    const ccomplex_t IT_0131 = IT_0100*IT_0130;
    const ccomplex_t IT_0132 = IT_0098*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = -IT_0036 + -IT_0064 + -IT_0097 + -IT_0133;
    const ccomplex_t IT_0135 = pow(m_sb_1, 2);
    const ccomplex_t IT_0136 = cpow((-2)*s_23 + IT_0006 + IT_0135 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0137 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0012
      *IT_0137*V_ub_mod;
    const ccomplex_t IT_0139 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0008*IT_0010*IT_0012*IT_0139*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0138*IT_0141;
    const ccomplex_t IT_0143 = IT_0136*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_23 + IT_0006 + IT_0135 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0012
      *IT_0137*V_ub_mod;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0008*IT_0010*IT_0012*IT_0139*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = IT_0146*IT_0148;
    const ccomplex_t IT_0150 = IT_0145*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0144 + IT_0151;
    const ccomplex_t IT_0153 = conj(N_B1)*e_em;
    const ccomplex_t IT_0154 = IT_0025*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W1)*e_em;
    const ccomplex_t IT_0157 = IT_0012*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0155 + 3*IT_0158);
    const ccomplex_t IT_0160 = 0.166666666666667*IT_0159;
    const ccomplex_t IT_0161 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0162 = IT_0025*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = m_b*N_d1*e_em*IT_0008*U_sb_00;
    const ccomplex_t IT_0165 = IT_0017*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0163 + 1.5*IT_0166);
    const ccomplex_t IT_0168 = (-0.333333333333333)*IT_0167;
    const ccomplex_t IT_0169 = IT_0160*IT_0168;
    const ccomplex_t IT_0170 = IT_0037*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = conj(N_B2)*e_em;
    const ccomplex_t IT_0173 = IT_0025*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = conj(N_W2)*e_em;
    const ccomplex_t IT_0176 = IT_0012*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + 3*IT_0177);
    const ccomplex_t IT_0179 = 0.166666666666667*IT_0178;
    const ccomplex_t IT_0180 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0181 = IT_0025*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = m_b*N_d2*e_em*IT_0008*U_sb_00;
    const ccomplex_t IT_0184 = IT_0017*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0182 + 1.5*IT_0185);
    const ccomplex_t IT_0187 = 0.666666666666667*IT_0186;
    const ccomplex_t IT_0188 = 1.5*IT_0187;
    const ccomplex_t IT_0189 = 0.666666666666667*IT_0188;
    const ccomplex_t IT_0190 = 1.5*IT_0189;
    const ccomplex_t IT_0191 = 0.666666666666667*IT_0190;
    const ccomplex_t IT_0192 = 1.5*IT_0191;
    const ccomplex_t IT_0193 = 0.666666666666667*IT_0192;
    const ccomplex_t IT_0194 = 1.5*IT_0193;
    const ccomplex_t IT_0195 = 0.666666666666667*IT_0194;
    const ccomplex_t IT_0196 = 1.5*IT_0195;
    const ccomplex_t IT_0197 = 0.666666666666667*IT_0196;
    const ccomplex_t IT_0198 = 1.5*IT_0197;
    const ccomplex_t IT_0199 = 0.666666666666667*IT_0198;
    const ccomplex_t IT_0200 = 1.5*IT_0199;
    const ccomplex_t IT_0201 = (-0.333333333333333)*IT_0200;
    const ccomplex_t IT_0202 = IT_0179*IT_0201;
    const ccomplex_t IT_0203 = IT_0065*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = conj(N_B3)*e_em;
    const ccomplex_t IT_0206 = IT_0025*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = conj(N_W3)*e_em;
    const ccomplex_t IT_0209 = IT_0012*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0207 + 3*IT_0210);
    const ccomplex_t IT_0212 = 0.166666666666667*IT_0211;
    const ccomplex_t IT_0213 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0214 = IT_0025*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = m_b*N_d3*e_em*IT_0008*U_sb_00;
    const ccomplex_t IT_0217 = IT_0017*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*(IT_0215 + 1.5*IT_0218);
    const ccomplex_t IT_0220 = 0.666666666666667*IT_0219;
    const ccomplex_t IT_0221 = 1.5*IT_0220;
    const ccomplex_t IT_0222 = (-0.333333333333333)*IT_0221;
    const ccomplex_t IT_0223 = IT_0212*IT_0222;
    const ccomplex_t IT_0224 = IT_0098*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = conj(N_B4)*e_em;
    const ccomplex_t IT_0227 = IT_0025*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = conj(N_W4)*e_em;
    const ccomplex_t IT_0230 = IT_0012*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0228 + 3*IT_0231);
    const ccomplex_t IT_0233 = 0.166666666666667*IT_0232;
    const ccomplex_t IT_0234 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0235 = IT_0025*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = m_b*N_d4*e_em*IT_0008*U_sb_00;
    const ccomplex_t IT_0238 = IT_0017*IT_0237;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*(IT_0236 + 1.5*IT_0239);
    const ccomplex_t IT_0241 = 0.666666666666667*IT_0240;
    const ccomplex_t IT_0242 = 1.5*IT_0241;
    const ccomplex_t IT_0243 = (-0.333333333333333)*IT_0242;
    const ccomplex_t IT_0244 = IT_0233*IT_0243;
    const ccomplex_t IT_0245 = IT_0007*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*IT_0245;
    const ccomplex_t IT_0247 = -IT_0171 + -IT_0204 + -IT_0225 + -IT_0246;
    const ccomplex_t IT_0248 = pow(m_sG, 2);
    const ccomplex_t IT_0249 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -IT_0248 +
       reg_prop, -1);
    const ccomplex_t IT_0250 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0252 = IT_0250*IT_0251;
    const ccomplex_t IT_0253 = IT_0249*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = 0.166666666666667*IT_0254;
    const ccomplex_t IT_0256 = conj(IT_0247) + conj(IT_0255);
    const ccomplex_t IT_0257 = (-0.5)*IT_0254;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*e_em
      *IT_0008*IT_0012*IT_0016*IT_0137*V_ub_mod;
    const ccomplex_t IT_0259 = 0.5*IT_0258;
    const ccomplex_t IT_0260 = IT_0012*IT_0139;
    const ccomplex_t IT_0261 = e_em*conj(U_Wm2)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0262 = IT_0260*IT_0261;
    const ccomplex_t IT_0263 = IT_0012*IT_0016*IT_0139;
    const ccomplex_t IT_0264 = m_b*conj(U_d2)*e_em*IT_0008*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0265 = IT_0263*IT_0264;
    const ccomplex_t IT_0266 = 1.4142135623731*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*(IT_0262 + (-0.5)*IT_0266);
    const ccomplex_t IT_0268 = IT_0259*IT_0267;
    const ccomplex_t IT_0269 = IT_0145*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = e_em*conj(U_Wm1)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0272 = IT_0260*IT_0271;
    const ccomplex_t IT_0273 = m_b*conj(U_d1)*e_em*IT_0008*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0274 = IT_0263*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*(IT_0272 + (-0.5)*IT_0275);
    const ccomplex_t IT_0277 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*e_em
      *IT_0008*IT_0012*IT_0016*IT_0137*V_ub_mod;
    const ccomplex_t IT_0278 = 0.5*IT_0277;
    const ccomplex_t IT_0279 = IT_0276*IT_0278;
    const ccomplex_t IT_0280 = IT_0136*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = IT_0270 + IT_0281;
    const ccomplex_t IT_0283 = IT_0257 + IT_0282;
    const ccomplex_t IT_0284 = IT_0247 + IT_0255;
    const ccomplex_t IT_0285 = conj(IT_0257) + conj(IT_0282);
    const ccomplex_t IT_0286 = (-18)*IT_0001;
    const ccomplex_t IT_0287 = 36*IT_0003;
    const ccomplex_t IT_0288 = IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = m_u*IT_0254;
    const ccomplex_t IT_0290 = (-0.5)*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0292 = IT_0141*IT_0278;
    const ccomplex_t IT_0293 = IT_0136*IT_0291*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0295 = IT_0039*IT_0168;
    const ccomplex_t IT_0296 = IT_0037*IT_0294*IT_0295;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0298 = IT_0014*IT_0243;
    const ccomplex_t IT_0299 = IT_0007*IT_0297*IT_0298;
    const ccomplex_t IT_0300 = m_u*IT_0171;
    const ccomplex_t IT_0301 = m_u*IT_0204;
    const ccomplex_t IT_0302 = m_u*IT_0246;
    const ccomplex_t IT_0303 = IT_0296 + IT_0299 + -IT_0300 + -IT_0301 + 
      -IT_0302;
    const ccomplex_t IT_0304 = (-6)*IT_0225;
    const ccomplex_t IT_0305 = IT_0254 + IT_0304;
    const ccomplex_t IT_0306 = m_u*IT_0305;
    const ccomplex_t IT_0307 = 0.166666666666667*IT_0306;
    const ccomplex_t IT_0308 = conj(IT_0303) + conj(IT_0307);
    const ccomplex_t IT_0309 = conj(IT_0290) + conj(IT_0293);
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0311 = IT_0067*IT_0201;
    const ccomplex_t IT_0312 = IT_0065*IT_0310*IT_0311;
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0315 = IT_0148*IT_0259;
    const ccomplex_t IT_0316 = IT_0145*IT_0314*IT_0315;
    const ccomplex_t IT_0317 = -IT_0316;
    const ccomplex_t IT_0318 = m_b*IT_0152;
    const ccomplex_t IT_0319 = IT_0317 + IT_0318;
    const ccomplex_t IT_0320 = conj(IT_0317) + conj(IT_0318);
    const ccomplex_t IT_0321 = IT_0146*IT_0267;
    const ccomplex_t IT_0322 = IT_0145*IT_0314*IT_0321;
    const ccomplex_t IT_0323 = IT_0061*IT_0160;
    const ccomplex_t IT_0324 = IT_0037*IT_0294*IT_0323;
    const ccomplex_t IT_0325 = m_u*IT_0133;
    const ccomplex_t IT_0326 = m_u*IT_0036;
    const ccomplex_t IT_0327 = IT_0324 + -IT_0325 + -IT_0326;
    const ccomplex_t IT_0328 = IT_0064 + IT_0097;
    const ccomplex_t IT_0329 = m_u*IT_0328;
    const ccomplex_t IT_0330 = -IT_0329;
    const ccomplex_t IT_0331 = IT_0033*IT_0233;
    const ccomplex_t IT_0332 = IT_0007*IT_0297*IT_0331;
    const ccomplex_t IT_0333 = IT_0327 + IT_0330 + IT_0332;
    const ccomplex_t IT_0334 = conj(IT_0327) + conj(IT_0330) + conj(IT_0332);
    const ccomplex_t IT_0335 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0336 = IT_0250*IT_0335;
    const ccomplex_t IT_0337 = IT_0249*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*IT_0337;
    const ccomplex_t IT_0339 = m_sG*IT_0338;
    const ccomplex_t IT_0340 = 0.166666666666667*IT_0339;
    const ccomplex_t IT_0341 = IT_0138*IT_0276;
    const ccomplex_t IT_0342 = IT_0136*IT_0291*IT_0341;
    const ccomplex_t IT_0343 = (-0.5)*IT_0339 + -IT_0342;
    const ccomplex_t IT_0344 = m_b*IT_0282;
    const ccomplex_t IT_0345 = IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = conj(IT_0343) + conj(IT_0344);
    const ccomplex_t IT_0347 = 6*s_34;
    const ccomplex_t IT_0348 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0349 = IT_0100*IT_0222;
    const ccomplex_t IT_0350 = IT_0098*IT_0348*IT_0349;
    const ccomplex_t IT_0351 = -IT_0350;
    const ccomplex_t IT_0352 = conj(IT_0313) + conj(IT_0351);
    const ccomplex_t IT_0353 = IT_0130*IT_0212;
    const ccomplex_t IT_0354 = IT_0098*IT_0348*IT_0353;
    const ccomplex_t IT_0355 = -IT_0354;
    const ccomplex_t IT_0356 = conj(IT_0340) + conj(IT_0355);
    const ccomplex_t IT_0357 = conj(IT_0330) + conj(IT_0332);
    const ccomplex_t IT_0358 = 18*s_34;
    const ccomplex_t IT_0359 = IT_0290*IT_0358;
    const ccomplex_t IT_0360 = IT_0327*IT_0358;
    const ccomplex_t IT_0361 = IT_0094*IT_0179;
    const ccomplex_t IT_0362 = IT_0065*IT_0310*IT_0361;
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = IT_0358*conj(IT_0363);
    const ccomplex_t IT_0365 = IT_0355*IT_0358;
    const ccomplex_t IT_0366 = IT_0330 + IT_0332;
    const ccomplex_t IT_0367 = (-18)*s_34;
    const ccomplex_t IT_0368 = conj(IT_0340) + conj(IT_0355) + conj(IT_0363);
    const ccomplex_t IT_0369 = IT_0358*IT_0363;
    const ccomplex_t IT_0370 = IT_0327*IT_0367;
    const ccomplex_t IT_0371 = IT_0340 + IT_0355 + IT_0363;
    const ccomplex_t IT_0372 = conj(IT_0327)*IT_0367;
    const ccomplex_t IT_0373 = IT_0313 + IT_0351;
    const ccomplex_t IT_0374 = (-6)*s_34;
    const ccomplex_t IT_0375 = IT_0351*IT_0358;
    const ccomplex_t IT_0376 = IT_0307*IT_0367;
    const ccomplex_t IT_0377 = IT_0290*IT_0374;
    const ccomplex_t IT_0378 = conj(IT_0303)*IT_0367;
    const ccomplex_t IT_0379 = conj(IT_0307)*IT_0367;
    const ccomplex_t IT_0380 = conj(IT_0290)*IT_0374;
    const ccomplex_t IT_0381 = IT_0347*IT_0351;
    const ccomplex_t IT_0382 = IT_0317*IT_0358;
    const ccomplex_t IT_0383 = IT_0290*IT_0367;
    const ccomplex_t IT_0384 = IT_0307*IT_0374;
    const ccomplex_t IT_0385 = IT_0347*conj(IT_0351);
    const ccomplex_t IT_0386 = conj(IT_0290)*IT_0367;
    const ccomplex_t IT_0387 = conj(IT_0303)*IT_0374;
    const ccomplex_t IT_0388 = conj(IT_0307)*IT_0374;
    const ccomplex_t IT_0389 = IT_0347*IT_0363;
    const ccomplex_t IT_0390 = IT_0347*IT_0355;
    const ccomplex_t IT_0391 = IT_0343*IT_0358;
    const ccomplex_t IT_0392 = IT_0327*IT_0374;
    const ccomplex_t IT_0393 = IT_0347*conj(IT_0363);
    const ccomplex_t IT_0394 = IT_0347*conj(IT_0355);
    const ccomplex_t IT_0395 = conj(IT_0327)*IT_0374;
    const ccomplex_t IT_0396 = IT_0330*IT_0358;
    const ccomplex_t IT_0397 = IT_0363*IT_0367;
    const ccomplex_t IT_0398 = IT_0355*IT_0367;
    const ccomplex_t IT_0399 = IT_0343*IT_0374;
    const ccomplex_t IT_0400 = conj(IT_0327)*IT_0358;
    const ccomplex_t IT_0401 = conj(IT_0363)*IT_0367;
    const ccomplex_t IT_0402 = conj(IT_0355)*IT_0367;
    const ccomplex_t IT_0403 = conj(IT_0343)*IT_0374;
    const ccomplex_t IT_0404 = m_u*s_14;
    const ccomplex_t IT_0405 = 6*IT_0351;
    const ccomplex_t IT_0406 = 6*IT_0313;
    const ccomplex_t IT_0407 = 18*IT_0317;
    const ccomplex_t IT_0408 = (-18)*IT_0290;
    const ccomplex_t IT_0409 = (-6)*IT_0303;
    const ccomplex_t IT_0410 = (-6)*IT_0307;
    const ccomplex_t IT_0411 = IT_0405 + IT_0406 + IT_0407 + IT_0408 + IT_0409
       + IT_0410;
    const ccomplex_t IT_0412 = 6*conj(IT_0351);
    const ccomplex_t IT_0413 = 6*conj(IT_0313);
    const ccomplex_t IT_0414 = 18*conj(IT_0317);
    const ccomplex_t IT_0415 = (-18)*conj(IT_0290);
    const ccomplex_t IT_0416 = (-6)*conj(IT_0303);
    const ccomplex_t IT_0417 = (-6)*conj(IT_0307);
    const ccomplex_t IT_0418 = IT_0412 + IT_0413 + IT_0414 + IT_0415 + IT_0416
       + IT_0417;
    const ccomplex_t IT_0419 = 6*IT_0317;
    const ccomplex_t IT_0420 = 18*IT_0351;
    const ccomplex_t IT_0421 = 18*IT_0313;
    const ccomplex_t IT_0422 = (-18)*IT_0303;
    const ccomplex_t IT_0423 = (-18)*IT_0307;
    const ccomplex_t IT_0424 = (-6)*IT_0290;
    const ccomplex_t IT_0425 = IT_0419 + IT_0420 + IT_0421 + IT_0422 + IT_0423
       + IT_0424;
    const ccomplex_t IT_0426 = 6*conj(IT_0317);
    const ccomplex_t IT_0427 = 18*conj(IT_0351);
    const ccomplex_t IT_0428 = 18*conj(IT_0313);
    const ccomplex_t IT_0429 = (-18)*conj(IT_0303);
    const ccomplex_t IT_0430 = (-18)*conj(IT_0307);
    const ccomplex_t IT_0431 = (-6)*conj(IT_0290);
    const ccomplex_t IT_0432 = IT_0426 + IT_0427 + IT_0428 + IT_0429 + IT_0430
       + IT_0431;
    const ccomplex_t IT_0433 = 6*IT_0134;
    const ccomplex_t IT_0434 = 18*IT_0152;
    const ccomplex_t IT_0435 = IT_0433 + IT_0434;
    const ccomplex_t IT_0436 = 6*conj(IT_0134);
    const ccomplex_t IT_0437 = 18*conj(IT_0152);
    const ccomplex_t IT_0438 = IT_0436 + IT_0437;
    const ccomplex_t IT_0439 = 6*IT_0247;
    const ccomplex_t IT_0440 = 6*IT_0255;
    const ccomplex_t IT_0441 = 18*IT_0257;
    const ccomplex_t IT_0442 = 18*IT_0282;
    const ccomplex_t IT_0443 = IT_0439 + IT_0440 + IT_0441 + IT_0442;
    const ccomplex_t IT_0444 = 6*conj(IT_0247);
    const ccomplex_t IT_0445 = 6*conj(IT_0255);
    const ccomplex_t IT_0446 = 18*conj(IT_0257);
    const ccomplex_t IT_0447 = 18*conj(IT_0282);
    const ccomplex_t IT_0448 = IT_0444 + IT_0445 + IT_0446 + IT_0447;
    const ccomplex_t IT_0449 = (-18)*IT_0257;
    const ccomplex_t IT_0450 = (-18)*IT_0282;
    const ccomplex_t IT_0451 = (-6)*IT_0247;
    const ccomplex_t IT_0452 = (-6)*IT_0255;
    const ccomplex_t IT_0453 = IT_0449 + IT_0450 + IT_0451 + IT_0452;
    const ccomplex_t IT_0454 = (-18)*conj(IT_0257);
    const ccomplex_t IT_0455 = (-18)*conj(IT_0282);
    const ccomplex_t IT_0456 = (-6)*conj(IT_0247);
    const ccomplex_t IT_0457 = (-6)*conj(IT_0255);
    const ccomplex_t IT_0458 = IT_0454 + IT_0455 + IT_0456 + IT_0457;
    const ccomplex_t IT_0459 = (-18)*IT_0152;
    const ccomplex_t IT_0460 = (-6)*IT_0134;
    const ccomplex_t IT_0461 = IT_0459 + IT_0460;
    const ccomplex_t IT_0462 = (-18)*conj(IT_0152);
    const ccomplex_t IT_0463 = (-6)*conj(IT_0134);
    const ccomplex_t IT_0464 = IT_0462 + IT_0463;
    const ccomplex_t IT_0465 = (-18)*IT_0134;
    const ccomplex_t IT_0466 = (-18)*conj(IT_0134);
    const ccomplex_t IT_0467 = 6*IT_0363;
    const ccomplex_t IT_0468 = 6*IT_0355;
    const ccomplex_t IT_0469 = 6*IT_0340;
    const ccomplex_t IT_0470 = 18*IT_0343;
    const ccomplex_t IT_0471 = (-6)*IT_0327;
    const ccomplex_t IT_0472 = (-6)*IT_0330;
    const ccomplex_t IT_0473 = 6*conj(IT_0363);
    const ccomplex_t IT_0474 = 6*conj(IT_0355);
    const ccomplex_t IT_0475 = 6*conj(IT_0340);
    const ccomplex_t IT_0476 = 18*conj(IT_0343);
    const ccomplex_t IT_0477 = (-6)*conj(IT_0327);
    const ccomplex_t IT_0478 = (-6)*conj(IT_0330);
    const ccomplex_t IT_0479 = conj(IT_0134)*((-18)*IT_0327 + (-18)*IT_0330 +
       18*IT_0340 + 6*IT_0343 + 18*IT_0355 + 18*IT_0363) + IT_0134*((-18)*conj
      (IT_0327) + (-18)*conj(IT_0330) + 18*conj(IT_0340) + 6*conj(IT_0343) + 18
      *conj(IT_0355) + 18*conj(IT_0363)) + IT_0285*IT_0411 + IT_0283*IT_0418 +
       IT_0256*IT_0425 + IT_0284*IT_0432 + conj(IT_0344)*IT_0435 + IT_0344
      *IT_0438 + conj(IT_0318)*IT_0443 + IT_0318*IT_0448 + conj(IT_0293)*IT_0453
       + IT_0293*IT_0458 + conj(IT_0322)*IT_0461 + IT_0322*IT_0464 + conj
      (IT_0332)*IT_0465 + IT_0332*IT_0466 + (-6)*conj(IT_0152)*(IT_0332 + (
      -0.166666666666667)*IT_0467 + (-0.166666666666667)*IT_0468 + (
      -0.166666666666667)*IT_0469 + (-0.166666666666667)*IT_0470 + (
      -0.166666666666667)*IT_0471 + (-0.166666666666667)*IT_0472) + (-6)*IT_0152
      *(conj(IT_0332) + (-0.166666666666667)*IT_0473 + (-0.166666666666667)
      *IT_0474 + (-0.166666666666667)*IT_0475 + (-0.166666666666667)*IT_0476 + (
      -0.166666666666667)*IT_0477 + (-0.166666666666667)*IT_0478);
    const ccomplex_t IT_0480 = m_b*m_u;
    const ccomplex_t IT_0481 = 6*IT_0303;
    const ccomplex_t IT_0482 = 6*IT_0307;
    const ccomplex_t IT_0483 = 18*IT_0290;
    const ccomplex_t IT_0484 = (-18)*IT_0317;
    const ccomplex_t IT_0485 = (-6)*IT_0351;
    const ccomplex_t IT_0486 = (-6)*IT_0313;
    const ccomplex_t IT_0487 = IT_0481 + IT_0482 + IT_0483 + IT_0484 + IT_0485
       + IT_0486;
    const ccomplex_t IT_0488 = 6*conj(IT_0303);
    const ccomplex_t IT_0489 = 6*conj(IT_0307);
    const ccomplex_t IT_0490 = 18*conj(IT_0290);
    const ccomplex_t IT_0491 = (-18)*conj(IT_0317);
    const ccomplex_t IT_0492 = (-6)*conj(IT_0351);
    const ccomplex_t IT_0493 = (-6)*conj(IT_0313);
    const ccomplex_t IT_0494 = IT_0488 + IT_0489 + IT_0490 + IT_0491 + IT_0492
       + IT_0493;
    const ccomplex_t IT_0495 = 6*IT_0290;
    const ccomplex_t IT_0496 = 18*IT_0303;
    const ccomplex_t IT_0497 = 18*IT_0307;
    const ccomplex_t IT_0498 = (-18)*IT_0351;
    const ccomplex_t IT_0499 = (-18)*IT_0313;
    const ccomplex_t IT_0500 = (-6)*IT_0317;
    const ccomplex_t IT_0501 = IT_0495 + IT_0496 + IT_0497 + IT_0498 + IT_0499
       + IT_0500;
    const ccomplex_t IT_0502 = 6*conj(IT_0290);
    const ccomplex_t IT_0503 = 18*conj(IT_0303);
    const ccomplex_t IT_0504 = 18*conj(IT_0307);
    const ccomplex_t IT_0505 = (-18)*conj(IT_0351);
    const ccomplex_t IT_0506 = (-18)*conj(IT_0313);
    const ccomplex_t IT_0507 = (-6)*conj(IT_0317);
    const ccomplex_t IT_0508 = IT_0502 + IT_0503 + IT_0504 + IT_0505 + IT_0506
       + IT_0507;
    const ccomplex_t IT_0509 = 6*IT_0327;
    const ccomplex_t IT_0510 = 6*IT_0330;
    const ccomplex_t IT_0511 = (-18)*IT_0343;
    const ccomplex_t IT_0512 = (-6)*IT_0363;
    const ccomplex_t IT_0513 = (-6)*IT_0355;
    const ccomplex_t IT_0514 = (-6)*IT_0340;
    const ccomplex_t IT_0515 = 6*conj(IT_0327);
    const ccomplex_t IT_0516 = 6*conj(IT_0330);
    const ccomplex_t IT_0517 = (-18)*conj(IT_0343);
    const ccomplex_t IT_0518 = (-6)*conj(IT_0363);
    const ccomplex_t IT_0519 = (-6)*conj(IT_0355);
    const ccomplex_t IT_0520 = (-6)*conj(IT_0340);
    const ccomplex_t IT_0521 = conj(IT_0322)*IT_0411 + IT_0322*IT_0418 +
       IT_0334*IT_0425 + IT_0333*IT_0432 + conj(IT_0293)*((-18)*IT_0322 + (-6)
      *IT_0332 + 18*IT_0344 + IT_0467 + IT_0468 + IT_0469 + IT_0470 + IT_0471 +
       IT_0472) + IT_0293*((-18)*conj(IT_0322) + (-6)*conj(IT_0332) + 18*conj
      (IT_0344) + IT_0473 + IT_0474 + IT_0475 + IT_0476 + IT_0477 + IT_0478) +
       IT_0346*IT_0487 + IT_0345*IT_0494 + IT_0368*IT_0501 + IT_0371*IT_0508 +
       conj(IT_0318)*(18*IT_0322 + 6*IT_0332 + (-18)*IT_0344 + IT_0509 + IT_0510
       + IT_0511 + IT_0512 + IT_0513 + IT_0514) + IT_0318*(18*conj(IT_0322) + 6
      *conj(IT_0332) + (-18)*conj(IT_0344) + IT_0515 + IT_0516 + IT_0517 +
       IT_0518 + IT_0519 + IT_0520);
    const ccomplex_t IT_0522 = m_b*s_13;
    const ccomplex_t IT_0523 = 18*IT_0327;
    const ccomplex_t IT_0524 = 18*IT_0330;
    const ccomplex_t IT_0525 = (-18)*IT_0363;
    const ccomplex_t IT_0526 = (-18)*IT_0355;
    const ccomplex_t IT_0527 = (-18)*IT_0340;
    const ccomplex_t IT_0528 = (-6)*IT_0343;
    const ccomplex_t IT_0529 = 18*conj(IT_0327);
    const ccomplex_t IT_0530 = 18*conj(IT_0330);
    const ccomplex_t IT_0531 = (-18)*conj(IT_0363);
    const ccomplex_t IT_0532 = (-18)*conj(IT_0355);
    const ccomplex_t IT_0533 = (-18)*conj(IT_0340);
    const ccomplex_t IT_0534 = (-6)*conj(IT_0343);
    const ccomplex_t IT_0535 = conj(IT_0293)*IT_0435 + IT_0293*IT_0438 + conj
      (IT_0322)*IT_0443 + IT_0322*IT_0448 + conj(IT_0344)*IT_0453 + IT_0344
      *IT_0458 + conj(IT_0318)*IT_0461 + IT_0318*IT_0464 + conj(IT_0152)*IT_0487
       + IT_0152*IT_0494 + conj(IT_0134)*IT_0501 + IT_0134*IT_0508 + 6*IT_0285*
      (IT_0332 + 0.166666666666667*IT_0509 + 0.166666666666667*IT_0510 +
       0.166666666666667*IT_0511 + 0.166666666666667*IT_0512 + 0.166666666666667
      *IT_0513 + 0.166666666666667*IT_0514) + 6*IT_0283*(conj(IT_0332) +
       0.166666666666667*IT_0515 + 0.166666666666667*IT_0516 + 0.166666666666667
      *IT_0517 + 0.166666666666667*IT_0518 + 0.166666666666667*IT_0519 +
       0.166666666666667*IT_0520) + IT_0256*(18*IT_0332 + IT_0523 + IT_0524 +
       IT_0525 + IT_0526 + IT_0527 + IT_0528) + IT_0284*(18*conj(IT_0332) +
       IT_0529 + IT_0530 + IT_0531 + IT_0532 + IT_0533 + IT_0534);
    const ccomplex_t IT_0536 = m_b*m_u*IT_0000;
    const ccomplex_t IT_0537 = IT_0005*(conj(IT_0134)*IT_0152 + IT_0134*conj
      (IT_0152) + IT_0256*IT_0283 + IT_0284*IT_0285) + (IT_0134*conj(IT_0134) +
       IT_0152*conj(IT_0152) + IT_0256*IT_0284 + IT_0283*IT_0285)*IT_0288 + (
      (IT_0290 + IT_0293)*IT_0308 + (IT_0303 + IT_0307)*IT_0309 + conj(IT_0313)
      *IT_0319 + IT_0313*IT_0320 + conj(IT_0322)*IT_0333 + IT_0322*IT_0334 +
       conj(IT_0340)*IT_0345 + IT_0340*IT_0346)*IT_0347 + ((IT_0303 + IT_0307)
      *IT_0308 + IT_0293*IT_0309 + IT_0318*IT_0320 + IT_0322*conj(IT_0322) +
       IT_0344*IT_0346 + IT_0313*IT_0352 + IT_0340*IT_0356 + IT_0332*IT_0357)
      *IT_0358 + IT_0309*IT_0359 + conj(IT_0327)*IT_0360 + (IT_0340 + IT_0355)
      *IT_0364 + IT_0356*IT_0365 + (conj(IT_0293)*IT_0319 + IT_0293*IT_0320 +
       conj(IT_0322)*IT_0345 + IT_0322*IT_0346 + IT_0303*IT_0352 + IT_0340
      *IT_0357 + conj(IT_0340)*IT_0366)*IT_0367 + IT_0368*(IT_0369 + IT_0370) +
       IT_0371*IT_0372 + (IT_0303*IT_0320 + IT_0293*IT_0352 + IT_0344*IT_0357 +
       conj(IT_0344)*IT_0366 + IT_0322*IT_0368 + conj(IT_0322)*IT_0371 + conj
      (IT_0293)*IT_0373)*IT_0374 + IT_0352*(IT_0375 + IT_0376 + IT_0377) +
       IT_0373*(IT_0378 + IT_0379 + IT_0380) + IT_0320*(IT_0381 + IT_0382 +
       IT_0383 + IT_0384) + IT_0319*(IT_0385 + IT_0386 + IT_0387 + IT_0388) +
       IT_0346*(IT_0389 + IT_0390 + IT_0391 + IT_0392) + IT_0345*(IT_0393 +
       IT_0394 + IT_0395) + IT_0357*(IT_0360 + IT_0396 + IT_0397 + IT_0398 +
       IT_0399) + IT_0366*(IT_0400 + IT_0401 + IT_0402 + IT_0403) + IT_0404
      *IT_0479 + IT_0480*IT_0521 + IT_0522*IT_0535 + (conj(IT_0152)*IT_0453 +
       IT_0152*IT_0458 + IT_0285*IT_0460 + IT_0283*IT_0463 + IT_0256*IT_0465 +
       IT_0284*IT_0466)*IT_0536;
    return create_ccomplex_return(IT_0537);
}

