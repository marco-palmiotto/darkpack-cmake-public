#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sb_2_to_u_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sb_2_to_u_anti_b(
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
    const creal_t m_sb_2 = param->m_sb_2;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_su_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_B1)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = conj(N_W1)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0014 + 3*IT_0019);
    const ccomplex_t IT_0021 = 0.166666666666667*IT_0020;
    const ccomplex_t IT_0022 = cos(beta);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = IT_0016*IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = m_b*N_d1*e_em*IT_0025*conj(U_sb_11);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0030 = IT_0011*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0033 = IT_0016*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0028 + 0.333333333333333
      *IT_0031 + -IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = (-3)*IT_0036;
    const ccomplex_t IT_0038 = 0.166666666666667*IT_0037;
    const ccomplex_t IT_0039 = IT_0021*IT_0038;
    const ccomplex_t IT_0040 = IT_0009*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0043 = conj(N_B2)*e_em;
    const ccomplex_t IT_0044 = IT_0011*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = conj(N_W2)*e_em;
    const ccomplex_t IT_0047 = IT_0016*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + 3*IT_0048);
    const ccomplex_t IT_0050 = 0.166666666666667*IT_0049;
    const ccomplex_t IT_0051 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0052 = IT_0016*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0055 = IT_0011*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = m_b*N_d2*e_em*IT_0025*conj(U_sb_11);
    const ccomplex_t IT_0058 = IT_0024*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0053 + (
      -0.333333333333333)*IT_0056 + -IT_0059);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = 3*IT_0061;
    const ccomplex_t IT_0063 = 0.166666666666667*IT_0062;
    const ccomplex_t IT_0064 = IT_0050*IT_0063;
    const ccomplex_t IT_0065 = IT_0042*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0068 = conj(N_B3)*e_em;
    const ccomplex_t IT_0069 = IT_0011*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = conj(N_W3)*e_em;
    const ccomplex_t IT_0072 = IT_0016*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + 3*IT_0073);
    const ccomplex_t IT_0075 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0076 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0077 = IT_0016*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = m_b*N_d3*e_em*IT_0025*conj(U_sb_11);
    const ccomplex_t IT_0080 = IT_0024*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0083 = IT_0011*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0078 + -IT_0081 + (
      -0.333333333333333)*IT_0084);
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = 3*IT_0086;
    const ccomplex_t IT_0088 = 0.166666666666667*IT_0087;
    const ccomplex_t IT_0089 = IT_0075*IT_0088;
    const ccomplex_t IT_0090 = IT_0067*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0093 = conj(N_B4)*e_em;
    const ccomplex_t IT_0094 = IT_0011*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = conj(N_W4)*e_em;
    const ccomplex_t IT_0097 = IT_0016*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + 3*IT_0098);
    const ccomplex_t IT_0100 = 0.166666666666667*IT_0099;
    const ccomplex_t IT_0101 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0102 = IT_0011*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0105 = IT_0016*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = m_b*N_d4*e_em*IT_0025*conj(U_sb_11);
    const ccomplex_t IT_0108 = IT_0024*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + (-3)*IT_0106 + 3
      *IT_0109);
    const ccomplex_t IT_0111 = (-0.333333333333333)*IT_0110;
    const ccomplex_t IT_0112 = (-3)*IT_0111;
    const ccomplex_t IT_0113 = 0.166666666666667*IT_0112;
    const ccomplex_t IT_0114 = IT_0100*IT_0113;
    const ccomplex_t IT_0115 = IT_0092*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 0.166666666666667*IT_0008 + -IT_0041 + -IT_0066
       + -IT_0091 + -IT_0116;
    const ccomplex_t IT_0118 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0016
      *IT_0118*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0016
      *IT_0121*V_ub_mod;
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = pow(m_sb_2, 2);
    const ccomplex_t IT_0125 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0002 + -IT_0124 + -reg_prop, -1);
    const ccomplex_t IT_0126 = IT_0123*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = IT_0120*IT_0127;
    const ccomplex_t IT_0129 = (-0.5)*IT_0008 + -IT_0128;
    const ccomplex_t IT_0130 = sin(beta);
    const ccomplex_t IT_0131 = cpow(IT_0130, -1);
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0016*IT_0025*IT_0131;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = m_b*conj(N_d4)*e_em*IT_0025*conj(U_sb_01);
    const ccomplex_t IT_0135 = IT_0024*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0138 = IT_0011*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0136 + 0.666666666666667
      *IT_0139);
    const ccomplex_t IT_0141 = 1.5*IT_0140;
    const ccomplex_t IT_0142 = (-0.333333333333333)*IT_0141;
    const ccomplex_t IT_0143 = IT_0133*IT_0142;
    const ccomplex_t IT_0144 = IT_0092*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0016*IT_0025*IT_0131;
    const ccomplex_t IT_0147 = (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0149 = IT_0011*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = m_b*conj(N_d1)*e_em*IT_0025*conj(U_sb_01);
    const ccomplex_t IT_0152 = IT_0024*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0150 + 1.5*IT_0153);
    const ccomplex_t IT_0155 = (-0.333333333333333)*IT_0154;
    const ccomplex_t IT_0156 = IT_0147*IT_0155;
    const ccomplex_t IT_0157 = IT_0009*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0016*IT_0025*IT_0131;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0162 = IT_0011*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = m_b*conj(N_d2)*e_em*IT_0025*conj(U_sb_01);
    const ccomplex_t IT_0165 = IT_0024*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0163 + 1.5*IT_0166);
    const ccomplex_t IT_0168 = (-0.333333333333333)*IT_0167;
    const ccomplex_t IT_0169 = IT_0160*IT_0168;
    const ccomplex_t IT_0170 = IT_0042*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0016*IT_0025*IT_0131;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0175 = IT_0011*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = m_b*conj(N_d3)*e_em*IT_0025*conj(U_sb_01);
    const ccomplex_t IT_0178 = IT_0024*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(IT_0176 + 1.5*IT_0179);
    const ccomplex_t IT_0181 = (-0.333333333333333)*IT_0180;
    const ccomplex_t IT_0182 = IT_0173*IT_0181;
    const ccomplex_t IT_0183 = IT_0067*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = -IT_0145 + -IT_0158 + -IT_0171 + -IT_0184;
    const ccomplex_t IT_0186 = m_b*m_u*IT_0002;
    const ccomplex_t IT_0187 = s_34*IT_0002;
    const ccomplex_t IT_0188 = (-18)*IT_0187;
    const ccomplex_t IT_0189 = s_13*s_14;
    const ccomplex_t IT_0190 = 36*IT_0189;
    const ccomplex_t IT_0191 = IT_0188 + IT_0190;
    const ccomplex_t IT_0192 = m_u*s_14;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0016
      *IT_0118*(m_W*IT_0022*IT_0130*conj(U_sb_01) + (-0.5)*IT_0025*(IT_0000
      *IT_0022*IT_0131*conj(U_sb_01) + m_b*(m_b*IT_0023*IT_0130*conj(U_sb_01) +
       mu_h*conj(U_sb_11))))*V_ub_mod;
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0016*IT_0022*IT_0025*IT_0121*IT_0131*V_ub_mod;
    const ccomplex_t IT_0196 = 0.5*IT_0195;
    const ccomplex_t IT_0197 = IT_0194*IT_0196;
    const ccomplex_t IT_0198 = cpow((-2)*s_12 + -IT_0002 + -IT_0124 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0199 = IT_0197*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = (-0.5)*IT_0119;
    const ccomplex_t IT_0202 = m_u*IT_0201;
    const ccomplex_t IT_0203 = m_u*IT_0120;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = IT_0127*IT_0205;
    const ccomplex_t IT_0207 = pow(m_W, 2);
    const ccomplex_t IT_0208 = cpow((-2)*s_12 + -IT_0002 + -IT_0124 + IT_0207 
      + -reg_prop, -1);
    const ccomplex_t IT_0209 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0016*IT_0025*IT_0121*V_ub_mod;
    const ccomplex_t IT_0210 = 0.5*IT_0209;
    const ccomplex_t IT_0211 = pow(m_b, 2);
    const ccomplex_t IT_0212 = cpow(IT_0022, 2);
    const ccomplex_t IT_0213 = cpow(IT_0130, 2);
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0016
      *IT_0118*((IT_0025*(IT_0000 + -IT_0211) + m_W*(IT_0212 + -IT_0213))*conj
      (U_sb_01) + m_b*mu_h*IT_0023*IT_0025*IT_0130*conj(U_sb_11))*V_ub_mod;
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = IT_0210*IT_0215;
    const ccomplex_t IT_0217 = IT_0208*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = m_u*IT_0008;
    const ccomplex_t IT_0220 = -IT_0200 + -IT_0206 + -IT_0218 + 0.5*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0222 = IT_0088*IT_0173;
    const ccomplex_t IT_0223 = IT_0067*IT_0221*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0225 = IT_0038*IT_0147;
    const ccomplex_t IT_0226 = IT_0009*IT_0224*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0228 = IT_0063*IT_0160;
    const ccomplex_t IT_0229 = IT_0042*IT_0227*IT_0228;
    const ccomplex_t IT_0230 = m_u*IT_0066;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0232 = IT_0113*IT_0133;
    const ccomplex_t IT_0233 = IT_0092*IT_0231*IT_0232;
    const ccomplex_t IT_0234 = m_u*IT_0041;
    const ccomplex_t IT_0235 = m_u*IT_0091;
    const ccomplex_t IT_0236 = m_u*IT_0116;
    const ccomplex_t IT_0237 = (-0.166666666666667)*IT_0219 + -IT_0223 + 
      -IT_0226 + -IT_0229 + IT_0230 + -IT_0233 + IT_0234 + IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0239 = IT_0004*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0241 = IT_0003*IT_0239*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0016*IT_0025*IT_0121*V_ub_mod;
    const ccomplex_t IT_0243 = (-0.5)*IT_0242;
    const ccomplex_t IT_0244 = IT_0215*IT_0243;
    const ccomplex_t IT_0245 = IT_0208*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0016*IT_0023*IT_0025*IT_0121*IT_0130*V_ub_mod;
    const ccomplex_t IT_0248 = 0.5*IT_0247;
    const ccomplex_t IT_0249 = IT_0194*IT_0248;
    const ccomplex_t IT_0250 = IT_0198*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = m_b*IT_0201;
    const ccomplex_t IT_0253 = m_b*IT_0120;
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = IT_0252 + IT_0254;
    const ccomplex_t IT_0256 = IT_0127*IT_0255;
    const ccomplex_t IT_0257 = (-0.5)*IT_0241 + -IT_0246 + -IT_0251 + IT_0256;
    const ccomplex_t IT_0258 = m_u*IT_0158;
    const ccomplex_t IT_0259 = IT_0075*IT_0181;
    const ccomplex_t IT_0260 = IT_0067*IT_0221*IT_0259;
    const ccomplex_t IT_0261 = IT_0021*IT_0155;
    const ccomplex_t IT_0262 = IT_0009*IT_0224*IT_0261;
    const ccomplex_t IT_0263 = IT_0050*IT_0168;
    const ccomplex_t IT_0264 = IT_0042*IT_0227*IT_0263;
    const ccomplex_t IT_0265 = IT_0100*IT_0142;
    const ccomplex_t IT_0266 = IT_0092*IT_0231*IT_0265;
    const ccomplex_t IT_0267 = m_u*IT_0171;
    const ccomplex_t IT_0268 = m_u*IT_0184;
    const ccomplex_t IT_0269 = m_u*IT_0145;
    const ccomplex_t IT_0270 = 0.166666666666667*IT_0241 + IT_0258 + -IT_0260 
      + -IT_0262 + -IT_0264 + -IT_0266 + IT_0267 + IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = 3*IT_0270;
    const ccomplex_t IT_0272 = IT_0257 + IT_0271;
    const ccomplex_t IT_0273 = 3*conj(IT_0270);
    const ccomplex_t IT_0274 = conj(IT_0257) + IT_0273;
    const ccomplex_t IT_0275 = conj(IT_0117)*(IT_0220 + 3*IT_0237) + 3*conj
      (IT_0129)*(IT_0220 + 0.333333333333333*IT_0237) + IT_0117*(conj(IT_0220) +
       3*conj(IT_0237)) + 3*IT_0129*(conj(IT_0220) + 0.333333333333333*conj
      (IT_0237)) + conj(IT_0185)*IT_0272 + IT_0185*IT_0274;
    const ccomplex_t IT_0276 = 6*s_34;
    const ccomplex_t IT_0277 = 18*s_34;
    const ccomplex_t IT_0278 = m_b*m_u;
    const ccomplex_t IT_0279 = (-18)*conj(IT_0237);
    const ccomplex_t IT_0280 = (-6)*conj(IT_0220);
    const ccomplex_t IT_0281 = IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = (-18)*IT_0237;
    const ccomplex_t IT_0283 = m_b*s_13;
    const ccomplex_t IT_0284 = conj(IT_0129)*(IT_0257 + 0.333333333333333
      *IT_0270) + IT_0129*(conj(IT_0257) + 0.333333333333333*conj(IT_0270)) +
       0.333333333333333*conj(IT_0117)*IT_0272 + 0.333333333333333*IT_0117
      *IT_0274 + (-0.0555555555555556)*IT_0185*IT_0281 + 0.333333333333333*conj
      (IT_0185)*(IT_0220 + (-0.166666666666667)*IT_0282);
    const ccomplex_t IT_0285 = 12*IT_0189;
    const ccomplex_t IT_0286 = (-6)*IT_0187;
    const ccomplex_t IT_0287 = (((-18)*conj(IT_0117) + (-6)*conj(IT_0129))
      *IT_0185 + ((-18)*IT_0117 + (-6)*IT_0129)*conj(IT_0185))*IT_0186 + 
      (IT_0117*conj(IT_0117) + IT_0129*conj(IT_0129) + IT_0185*conj(IT_0185))
      *IT_0191 + 6*IT_0192*IT_0275 + (conj(IT_0220)*IT_0237 + IT_0220*conj
      (IT_0237) + conj(IT_0257)*IT_0270 + IT_0257*conj(IT_0270))*IT_0276 + 
      (IT_0220*conj(IT_0220) + IT_0237*conj(IT_0237) + IT_0257*conj(IT_0257) +
       IT_0270*conj(IT_0270))*IT_0277 + IT_0278*(((-18)*conj(IT_0220) + (-6)
      *conj(IT_0237))*IT_0257 + (-6)*IT_0237*conj(IT_0257) + IT_0220*((-18)*conj
      (IT_0257) + (-6)*conj(IT_0270)) + IT_0270*IT_0281 + conj(IT_0270)*IT_0282)
       + (-18)*IT_0283*IT_0284 + (conj(IT_0117)*IT_0129 + IT_0117*conj(IT_0129))
      *(IT_0285 + IT_0286);
    return create_ccomplex_return(IT_0287);
}

