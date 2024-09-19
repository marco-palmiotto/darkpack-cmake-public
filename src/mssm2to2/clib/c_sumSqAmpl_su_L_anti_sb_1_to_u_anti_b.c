#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sb_1_to_u_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sb_1_to_u_anti_b(
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
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_su_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_B1)*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_W1)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0007 + 3*IT_0012);
    const ccomplex_t IT_0014 = 0.166666666666667*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = m_b*N_d1*e_em*IT_0018*conj(U_sb_10);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + 0.333333333333333
      *IT_0024 + -IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = (-3)*IT_0029;
    const ccomplex_t IT_0031 = 0.166666666666667*IT_0030;
    const ccomplex_t IT_0032 = IT_0014*IT_0031;
    const ccomplex_t IT_0033 = IT_0002*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0036 = conj(N_B2)*e_em;
    const ccomplex_t IT_0037 = IT_0004*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W2)*e_em;
    const ccomplex_t IT_0040 = IT_0009*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + 3*IT_0041);
    const ccomplex_t IT_0043 = 0.166666666666667*IT_0042;
    const ccomplex_t IT_0044 = m_b*N_d2*e_em*IT_0018*conj(U_sb_10);
    const ccomplex_t IT_0045 = IT_0017*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0048 = IT_0009*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0046 + -IT_0049 +
       0.333333333333333*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = (-3)*IT_0054;
    const ccomplex_t IT_0056 = 0.166666666666667*IT_0055;
    const ccomplex_t IT_0057 = IT_0043*IT_0056;
    const ccomplex_t IT_0058 = IT_0035*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0061 = conj(N_B3)*e_em;
    const ccomplex_t IT_0062 = IT_0004*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = conj(N_W3)*e_em;
    const ccomplex_t IT_0065 = IT_0009*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + 3*IT_0066);
    const ccomplex_t IT_0068 = 0.166666666666667*IT_0067;
    const ccomplex_t IT_0069 = m_b*N_d3*e_em*IT_0018*conj(U_sb_10);
    const ccomplex_t IT_0070 = IT_0017*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0073 = IT_0009*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0076 = IT_0004*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0071 + -IT_0074 +
       0.333333333333333*IT_0077);
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = (-3)*IT_0079;
    const ccomplex_t IT_0081 = 0.166666666666667*IT_0080;
    const ccomplex_t IT_0082 = IT_0068*IT_0081;
    const ccomplex_t IT_0083 = IT_0060*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0086 = conj(N_B4)*e_em;
    const ccomplex_t IT_0087 = IT_0004*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W4)*e_em;
    const ccomplex_t IT_0090 = IT_0009*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + 3*IT_0091);
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0092;
    const ccomplex_t IT_0094 = m_b*N_d4*e_em*IT_0018*conj(U_sb_10);
    const ccomplex_t IT_0095 = IT_0017*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0098 = IT_0009*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0101 = IT_0004*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0096 + -IT_0099 +
       0.333333333333333*IT_0102);
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = (-3)*IT_0104;
    const ccomplex_t IT_0106 = 0.166666666666667*IT_0105;
    const ccomplex_t IT_0107 = IT_0093*IT_0106;
    const ccomplex_t IT_0108 = IT_0085*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = -IT_0034 + -IT_0059 + -IT_0084 + -IT_0109;
    const ccomplex_t IT_0111 = pow(m_sG, 2);
    const ccomplex_t IT_0112 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0111 +
       reg_prop, -1);
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = IT_0112*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = 0.166666666666667*IT_0117;
    const ccomplex_t IT_0120 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0009
      *IT_0120*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0122 = -IT_0121;
    const ccomplex_t IT_0123 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0009
      *IT_0123*V_ub_mod;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = pow(m_sb_1, 2);
    const ccomplex_t IT_0127 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0001 + -IT_0126 + -reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0125*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0122*IT_0129;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = sin(beta);
    const ccomplex_t IT_0133 = cpow(IT_0132, -1);
    const ccomplex_t IT_0134 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0009*IT_0018*IT_0133;
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0137 = IT_0004*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = m_b*conj(N_d1)*e_em*IT_0018*conj(U_sb_00);
    const ccomplex_t IT_0140 = IT_0017*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + 1.5*IT_0141);
    const ccomplex_t IT_0143 = (-0.333333333333333)*IT_0142;
    const ccomplex_t IT_0144 = IT_0135*IT_0143;
    const ccomplex_t IT_0145 = IT_0002*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0009*IT_0018*IT_0133;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0150 = IT_0004*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = m_b*conj(N_d3)*e_em*IT_0018*conj(U_sb_00);
    const ccomplex_t IT_0153 = IT_0017*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0151 + 1.5*IT_0154);
    const ccomplex_t IT_0156 = 0.666666666666667*IT_0155;
    const ccomplex_t IT_0157 = 1.5*IT_0156;
    const ccomplex_t IT_0158 = (-0.333333333333333)*IT_0157;
    const ccomplex_t IT_0159 = IT_0148*IT_0158;
    const ccomplex_t IT_0160 = IT_0060*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0009*IT_0018*IT_0133;
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0165 = IT_0004*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = m_b*conj(N_d2)*e_em*IT_0018*conj(U_sb_00);
    const ccomplex_t IT_0168 = IT_0017*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0166 + 1.5*IT_0169);
    const ccomplex_t IT_0171 = 0.666666666666667*IT_0170;
    const ccomplex_t IT_0172 = 1.5*IT_0171;
    const ccomplex_t IT_0173 = 0.666666666666667*IT_0172;
    const ccomplex_t IT_0174 = 1.5*IT_0173;
    const ccomplex_t IT_0175 = 0.666666666666667*IT_0174;
    const ccomplex_t IT_0176 = 1.5*IT_0175;
    const ccomplex_t IT_0177 = 0.666666666666667*IT_0176;
    const ccomplex_t IT_0178 = 1.5*IT_0177;
    const ccomplex_t IT_0179 = (-0.333333333333333)*IT_0178;
    const ccomplex_t IT_0180 = IT_0163*IT_0179;
    const ccomplex_t IT_0181 = IT_0035*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0009*IT_0018*IT_0133;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0186 = IT_0004*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = m_b*conj(N_d4)*e_em*IT_0018*conj(U_sb_00);
    const ccomplex_t IT_0189 = IT_0017*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0187 + 1.5*IT_0190);
    const ccomplex_t IT_0192 = 0.666666666666667*IT_0191;
    const ccomplex_t IT_0193 = 1.5*IT_0192;
    const ccomplex_t IT_0194 = (-0.333333333333333)*IT_0193;
    const ccomplex_t IT_0195 = IT_0184*IT_0194;
    const ccomplex_t IT_0196 = IT_0085*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = -IT_0146 + -IT_0161 + -IT_0182 + -IT_0197;
    const ccomplex_t IT_0199 = m_b*m_u*IT_0001;
    const ccomplex_t IT_0200 = s_34*IT_0001;
    const ccomplex_t IT_0201 = (-18)*IT_0200;
    const ccomplex_t IT_0202 = s_13*s_14;
    const ccomplex_t IT_0203 = 36*IT_0202;
    const ccomplex_t IT_0204 = IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = IT_0110 + IT_0119;
    const ccomplex_t IT_0206 = conj(IT_0110) + conj(IT_0119);
    const ccomplex_t IT_0207 = conj(IT_0118) + conj(IT_0131);
    const ccomplex_t IT_0208 = (-6)*IT_0200;
    const ccomplex_t IT_0209 = 12*IT_0202;
    const ccomplex_t IT_0210 = IT_0208 + IT_0209;
    const ccomplex_t IT_0211 = IT_0118 + IT_0131;
    const ccomplex_t IT_0212 = m_u*IT_0117;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0214 = IT_0031*IT_0135;
    const ccomplex_t IT_0215 = IT_0002*IT_0213*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0217 = IT_0056*IT_0163;
    const ccomplex_t IT_0218 = IT_0035*IT_0216*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0220 = IT_0081*IT_0148;
    const ccomplex_t IT_0221 = IT_0060*IT_0219*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0223 = IT_0106*IT_0184;
    const ccomplex_t IT_0224 = IT_0085*IT_0222*IT_0223;
    const ccomplex_t IT_0225 = m_u*IT_0034;
    const ccomplex_t IT_0226 = m_u*IT_0059;
    const ccomplex_t IT_0227 = m_u*IT_0084;
    const ccomplex_t IT_0228 = m_u*IT_0109;
    const ccomplex_t IT_0229 = (-0.166666666666667)*IT_0212 + -IT_0215 + 
      -IT_0218 + -IT_0221 + -IT_0224 + IT_0225 + IT_0226 + IT_0227 + IT_0228;
    const ccomplex_t IT_0230 = pow(m_W, 2);
    const ccomplex_t IT_0231 = cpow((-2)*s_12 + -IT_0001 + -IT_0126 + IT_0230 
      + -reg_prop, -1);
    const ccomplex_t IT_0232 = pow(m_b, 2);
    const ccomplex_t IT_0233 = cpow(IT_0015, 2);
    const ccomplex_t IT_0234 = cpow(IT_0132, 2);
    const ccomplex_t IT_0235 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0009
      *IT_0120*((IT_0018*(IT_0000 + -IT_0232) + m_W*(IT_0233 + -IT_0234))*conj
      (U_sb_00) + m_b*mu_h*IT_0016*IT_0018*IT_0132*conj(U_sb_10))*V_ub_mod;
    const ccomplex_t IT_0236 = 0.5*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0009*IT_0018*IT_0123*V_ub_mod;
    const ccomplex_t IT_0238 = 0.5*IT_0237;
    const ccomplex_t IT_0239 = IT_0236*IT_0238;
    const ccomplex_t IT_0240 = IT_0231*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = cpow((-2)*s_12 + -IT_0001 + -IT_0126 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0243 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0009
      *IT_0120*(m_W*IT_0015*IT_0132*conj(U_sb_00) + (-0.5)*IT_0018*(IT_0000
      *IT_0015*IT_0133*conj(U_sb_00) + m_b*(m_b*IT_0016*IT_0132*conj(U_sb_00) +
       mu_h*conj(U_sb_10))))*V_ub_mod;
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0009*IT_0015*IT_0018*IT_0123*IT_0133*V_ub_mod;
    const ccomplex_t IT_0246 = 0.5*IT_0245;
    const ccomplex_t IT_0247 = IT_0244*IT_0246;
    const ccomplex_t IT_0248 = IT_0242*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = (-0.5)*IT_0121;
    const ccomplex_t IT_0251 = m_u*IT_0250;
    const ccomplex_t IT_0252 = m_u*IT_0122;
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = IT_0251 + IT_0253;
    const ccomplex_t IT_0255 = IT_0129*IT_0254;
    const ccomplex_t IT_0256 = 0.5*IT_0212 + -IT_0241 + -IT_0249 + -IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0259 = IT_0113*IT_0258;
    const ccomplex_t IT_0260 = IT_0112*IT_0257*IT_0259;
    const ccomplex_t IT_0261 = IT_0068*IT_0158;
    const ccomplex_t IT_0262 = IT_0060*IT_0219*IT_0261;
    const ccomplex_t IT_0263 = IT_0014*IT_0143;
    const ccomplex_t IT_0264 = IT_0002*IT_0213*IT_0263;
    const ccomplex_t IT_0265 = IT_0043*IT_0179;
    const ccomplex_t IT_0266 = IT_0035*IT_0216*IT_0265;
    const ccomplex_t IT_0267 = IT_0093*IT_0194;
    const ccomplex_t IT_0268 = IT_0085*IT_0222*IT_0267;
    const ccomplex_t IT_0269 = m_u*IT_0146;
    const ccomplex_t IT_0270 = m_u*IT_0182;
    const ccomplex_t IT_0271 = m_u*IT_0161;
    const ccomplex_t IT_0272 = m_u*IT_0197;
    const ccomplex_t IT_0273 = 0.166666666666667*IT_0260 + -IT_0262 + -IT_0264
       + -IT_0266 + -IT_0268 + IT_0269 + IT_0270 + IT_0271 + IT_0272;
    const ccomplex_t IT_0274 = (-0.5)*IT_0260;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0009*IT_0018*IT_0123*V_ub_mod;
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = IT_0236*IT_0276;
    const ccomplex_t IT_0278 = IT_0231*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0009*IT_0016*IT_0018*IT_0123*IT_0132*V_ub_mod;
    const ccomplex_t IT_0281 = 0.5*IT_0280;
    const ccomplex_t IT_0282 = IT_0244*IT_0281;
    const ccomplex_t IT_0283 = IT_0242*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = m_b*IT_0250;
    const ccomplex_t IT_0286 = m_b*IT_0122;
    const ccomplex_t IT_0287 = -IT_0286;
    const ccomplex_t IT_0288 = IT_0285 + IT_0287;
    const ccomplex_t IT_0289 = IT_0129*IT_0288;
    const ccomplex_t IT_0290 = -IT_0279 + -IT_0284 + IT_0289;
    const ccomplex_t IT_0291 = IT_0274 + IT_0290;
    const ccomplex_t IT_0292 = conj(IT_0274) + conj(IT_0290);
    const ccomplex_t IT_0293 = 6*s_34;
    const ccomplex_t IT_0294 = 18*s_34;
    const ccomplex_t IT_0295 = IT_0274*IT_0294;
    const ccomplex_t IT_0296 = m_u*s_14;
    const ccomplex_t IT_0297 = 6*IT_0256;
    const ccomplex_t IT_0298 = 18*IT_0229;
    const ccomplex_t IT_0299 = 6*conj(IT_0256);
    const ccomplex_t IT_0300 = 18*conj(IT_0229);
    const ccomplex_t IT_0301 = 6*IT_0229;
    const ccomplex_t IT_0302 = 18*IT_0256;
    const ccomplex_t IT_0303 = 6*conj(IT_0229);
    const ccomplex_t IT_0304 = 18*conj(IT_0256);
    const ccomplex_t IT_0305 = 3*conj(IT_0198)*(IT_0273 + 1./3*IT_0274 + 1./3
      *IT_0290) + 3*IT_0198*(conj(IT_0273) + 1./3*conj(IT_0274) + 1./3*conj
      (IT_0290)) + 0.166666666666667*IT_0206*(IT_0297 + IT_0298) +
       0.166666666666667*IT_0205*(IT_0299 + IT_0300) + 0.166666666666667*IT_0207
      *(IT_0301 + IT_0302) + 0.166666666666667*IT_0211*(IT_0303 + IT_0304);
    const ccomplex_t IT_0306 = m_b*m_u;
    const ccomplex_t IT_0307 = (-18)*conj(IT_0229);
    const ccomplex_t IT_0308 = (-6)*conj(IT_0256);
    const ccomplex_t IT_0309 = IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = (-18)*IT_0256;
    const ccomplex_t IT_0311 = (-18)*conj(IT_0256);
    const ccomplex_t IT_0312 = (-6)*conj(IT_0229);
    const ccomplex_t IT_0313 = (-6)*IT_0256;
    const ccomplex_t IT_0314 = m_b*s_13;
    const ccomplex_t IT_0315 = (-18)*IT_0274;
    const ccomplex_t IT_0316 = (-18)*IT_0290;
    const ccomplex_t IT_0317 = (-6)*IT_0273;
    const ccomplex_t IT_0318 = (-18)*conj(IT_0274);
    const ccomplex_t IT_0319 = (-18)*conj(IT_0290);
    const ccomplex_t IT_0320 = (-6)*conj(IT_0273);
    const ccomplex_t IT_0321 = (-18)*IT_0273;
    const ccomplex_t IT_0322 = (-6)*IT_0274;
    const ccomplex_t IT_0323 = (-6)*IT_0290;
    const ccomplex_t IT_0324 = (-18)*conj(IT_0273);
    const ccomplex_t IT_0325 = (-6)*conj(IT_0274);
    const ccomplex_t IT_0326 = (-6)*conj(IT_0290);
    const ccomplex_t IT_0327 = IT_0198*IT_0309 + conj(IT_0198)*((-18)*IT_0229 
      + IT_0313) + IT_0207*(IT_0315 + IT_0316 + IT_0317) + IT_0211*(IT_0318 +
       IT_0319 + IT_0320) + IT_0206*(IT_0321 + IT_0322 + IT_0323) + IT_0205*
      (IT_0324 + IT_0325 + IT_0326);
    const ccomplex_t IT_0328 = (((-18)*conj(IT_0110) + (-6)*conj(IT_0118) + (
      -18)*conj(IT_0119) + (-6)*conj(IT_0131))*IT_0198 + ((-18)*IT_0110 + (-6)
      *IT_0118 + (-18)*IT_0119 + (-6)*IT_0131)*conj(IT_0198))*IT_0199 + IT_0204*
      (IT_0198*conj(IT_0198) + IT_0205*IT_0206 + (IT_0118 + IT_0131)*IT_0207) +
       IT_0210*(IT_0205*IT_0207 + IT_0206*IT_0211) + (conj(IT_0229)*IT_0256 +
       IT_0229*conj(IT_0256) + conj(IT_0273)*IT_0291 + IT_0273*IT_0292)*IT_0293 
      + (IT_0229*conj(IT_0229) + IT_0256*conj(IT_0256) + IT_0273*conj(IT_0273) +
       IT_0290*IT_0292)*IT_0294 + IT_0292*IT_0295 + 6*IT_0296*IT_0305 + (-18)
      *IT_0306*(IT_0229*(conj(IT_0273) + 1./3*conj(IT_0274) + 1./3*conj(IT_0290)
      ) + (-0.0555555555555555)*IT_0273*IT_0309 + (-0.0555555555555555)*IT_0292
      *IT_0310 + (-0.0555555555555555)*IT_0291*(IT_0311 + IT_0312) + (
      -0.0555555555555555)*conj(IT_0273)*IT_0313) + IT_0314*IT_0327;
    return create_ccomplex_return(IT_0328);
}

