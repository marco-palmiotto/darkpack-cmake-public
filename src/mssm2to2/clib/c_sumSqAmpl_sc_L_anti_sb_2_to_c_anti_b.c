#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_sb_2_to_c_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_sb_2_to_c_anti_b(
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sc_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0010*conj(U_sb_01);
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0010;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = pow(m_sb_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0002 + -IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0012*IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0008 + -IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0022 = cos(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = conj(N_B2)*e_em;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W2)*e_em;
    const ccomplex_t IT_0028 = IT_0010*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0026 + 3*IT_0029);
    const ccomplex_t IT_0031 = 0.166666666666667*IT_0030;
    const ccomplex_t IT_0032 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0033 = IT_0023*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0036 = IT_0010*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = cos(beta);
    const ccomplex_t IT_0039 = cpow(IT_0038, -1);
    const ccomplex_t IT_0040 = IT_0010*IT_0039;
    const ccomplex_t IT_0041 = pow(m_W, -1);
    const ccomplex_t IT_0042 = m_b*N_d2*e_em*IT_0041*conj(U_sb_11);
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0034 + (-3)*IT_0037 + 3
      *IT_0044);
    const ccomplex_t IT_0046 = 0.166666666666667*IT_0045;
    const ccomplex_t IT_0047 = IT_0031*IT_0046;
    const ccomplex_t IT_0048 = IT_0021*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0051 = conj(N_B1)*e_em;
    const ccomplex_t IT_0052 = IT_0023*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W1)*e_em;
    const ccomplex_t IT_0055 = IT_0010*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 3*IT_0056);
    const ccomplex_t IT_0058 = 0.166666666666667*IT_0057;
    const ccomplex_t IT_0059 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0060 = IT_0023*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0063 = IT_0010*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = m_b*N_d1*e_em*IT_0041*conj(U_sb_11);
    const ccomplex_t IT_0066 = IT_0040*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0061 + (-3)*IT_0064 + 3
      *IT_0067);
    const ccomplex_t IT_0069 = 0.166666666666667*IT_0068;
    const ccomplex_t IT_0070 = IT_0058*IT_0069;
    const ccomplex_t IT_0071 = IT_0050*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0074 = conj(N_B3)*e_em;
    const ccomplex_t IT_0075 = IT_0023*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W3)*e_em;
    const ccomplex_t IT_0078 = IT_0010*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0076 + 3*IT_0079);
    const ccomplex_t IT_0081 = 0.166666666666667*IT_0080;
    const ccomplex_t IT_0082 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0083 = IT_0023*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0086 = IT_0010*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = m_b*N_d3*e_em*IT_0041*conj(U_sb_11);
    const ccomplex_t IT_0089 = IT_0040*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0084 + (-3)*IT_0087 + 3
      *IT_0090);
    const ccomplex_t IT_0092 = 0.166666666666667*IT_0091;
    const ccomplex_t IT_0093 = IT_0081*IT_0092;
    const ccomplex_t IT_0094 = IT_0073*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0097 = conj(N_B4)*e_em;
    const ccomplex_t IT_0098 = IT_0023*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = conj(N_W4)*e_em;
    const ccomplex_t IT_0101 = IT_0010*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0099 + 3*IT_0102);
    const ccomplex_t IT_0104 = 0.166666666666667*IT_0103;
    const ccomplex_t IT_0105 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0106 = IT_0023*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0109 = IT_0010*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = m_b*N_d4*e_em*IT_0041*conj(U_sb_11);
    const ccomplex_t IT_0112 = IT_0040*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0107 + (-3)*IT_0110 + 3
      *IT_0113);
    const ccomplex_t IT_0115 = 0.166666666666667*IT_0114;
    const ccomplex_t IT_0116 = IT_0104*IT_0115;
    const ccomplex_t IT_0117 = IT_0096*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = 0.166666666666667*IT_0008 + -IT_0049 + -IT_0072
       + -IT_0095 + -IT_0118;
    const ccomplex_t IT_0120 = sin(beta);
    const ccomplex_t IT_0121 = cpow(IT_0120, -1);
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0010*IT_0041*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0125 = IT_0023*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = m_b*conj(N_d1)*e_em*IT_0041*conj(U_sb_01);
    const ccomplex_t IT_0128 = IT_0040*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + 1.5*IT_0129);
    const ccomplex_t IT_0131 = (-0.333333333333333)*IT_0130;
    const ccomplex_t IT_0132 = IT_0123*IT_0131;
    const ccomplex_t IT_0133 = IT_0050*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0010*IT_0041*IT_0121;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0138 = IT_0023*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = m_b*conj(N_d2)*e_em*IT_0041*conj(U_sb_01);
    const ccomplex_t IT_0141 = IT_0040*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + 1.5*IT_0142);
    const ccomplex_t IT_0144 = (-0.333333333333333)*IT_0143;
    const ccomplex_t IT_0145 = IT_0136*IT_0144;
    const ccomplex_t IT_0146 = IT_0021*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0010*IT_0041*IT_0121;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0151 = IT_0023*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = m_b*conj(N_d3)*e_em*IT_0041*conj(U_sb_01);
    const ccomplex_t IT_0154 = IT_0040*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0152 + 1.5*IT_0155);
    const ccomplex_t IT_0157 = (-0.333333333333333)*IT_0156;
    const ccomplex_t IT_0158 = IT_0149*IT_0157;
    const ccomplex_t IT_0159 = IT_0073*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0010*IT_0041*IT_0121;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0164 = IT_0023*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = m_b*conj(N_d4)*e_em*IT_0041*conj(U_sb_01);
    const ccomplex_t IT_0167 = IT_0040*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0165 + 1.5*IT_0168);
    const ccomplex_t IT_0170 = (-0.333333333333333)*IT_0169;
    const ccomplex_t IT_0171 = IT_0162*IT_0170;
    const ccomplex_t IT_0172 = IT_0096*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = -IT_0134 + -IT_0147 + -IT_0160 + -IT_0173;
    const ccomplex_t IT_0175 = s_34*IT_0002;
    const ccomplex_t IT_0176 = (-18)*IT_0175;
    const ccomplex_t IT_0177 = s_13*s_14;
    const ccomplex_t IT_0178 = 36*IT_0177;
    const ccomplex_t IT_0179 = IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = m_c*s_14;
    const ccomplex_t IT_0181 = cpow((-2)*s_12 + -IT_0002 + -IT_0015 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0182 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0010*IT_0038*IT_0041*IT_0121;
    const ccomplex_t IT_0183 = 0.5*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0010*(m_W*IT_0038*IT_0120*conj(U_sb_01) + (-0.5)*IT_0041*(IT_0000
      *IT_0038*IT_0121*conj(U_sb_01) + m_b*(m_b*IT_0039*IT_0120*conj(U_sb_01) +
       mu_h*conj(U_sb_11))));
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = IT_0183*IT_0185;
    const ccomplex_t IT_0187 = IT_0181*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0010*IT_0041;
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = pow(m_b, 2);
    const ccomplex_t IT_0192 = cpow(IT_0038, 2);
    const ccomplex_t IT_0193 = cpow(IT_0120, 2);
    const ccomplex_t IT_0194 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0010*((IT_0041*(IT_0000 + -IT_0191) + m_W*(IT_0192 + -IT_0193))*conj
      (U_sb_01) + m_b*mu_h*IT_0039*IT_0041*IT_0120*conj(U_sb_11));
    const ccomplex_t IT_0195 = 0.5*IT_0194;
    const ccomplex_t IT_0196 = IT_0190*IT_0195;
    const ccomplex_t IT_0197 = pow(m_W, 2);
    const ccomplex_t IT_0198 = cpow((-2)*s_12 + -IT_0002 + -IT_0015 + IT_0197 
      + -reg_prop, -1);
    const ccomplex_t IT_0199 = IT_0196*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = m_c*IT_0008;
    const ccomplex_t IT_0202 = (-0.5)*IT_0011;
    const ccomplex_t IT_0203 = m_c*IT_0202;
    const ccomplex_t IT_0204 = m_c*IT_0012;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = IT_0018*IT_0206;
    const ccomplex_t IT_0208 = -IT_0188 + -IT_0200 + 0.5*IT_0201 + -IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0010*IT_0041;
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0195*IT_0210;
    const ccomplex_t IT_0212 = IT_0198*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0010*IT_0039*IT_0041*IT_0120;
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = IT_0185*IT_0215;
    const ccomplex_t IT_0217 = IT_0181*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0221 = IT_0004*IT_0220;
    const ccomplex_t IT_0222 = IT_0003*IT_0219*IT_0221;
    const ccomplex_t IT_0223 = m_b*IT_0202;
    const ccomplex_t IT_0224 = m_b*IT_0012;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = IT_0018*IT_0226;
    const ccomplex_t IT_0228 = -IT_0213 + -IT_0218 + (-0.5)*IT_0222 + IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0230 = IT_0031*IT_0144;
    const ccomplex_t IT_0231 = IT_0021*IT_0229*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0233 = IT_0104*IT_0170;
    const ccomplex_t IT_0234 = IT_0096*IT_0232*IT_0233;
    const ccomplex_t IT_0235 = m_c*IT_0134;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0237 = IT_0081*IT_0157;
    const ccomplex_t IT_0238 = IT_0073*IT_0236*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0240 = IT_0058*IT_0131;
    const ccomplex_t IT_0241 = IT_0050*IT_0239*IT_0240;
    const ccomplex_t IT_0242 = m_c*IT_0147;
    const ccomplex_t IT_0243 = m_c*IT_0160;
    const ccomplex_t IT_0244 = m_c*IT_0173;
    const ccomplex_t IT_0245 = 0.166666666666667*IT_0222 + -IT_0231 + -IT_0234
       + IT_0235 + -IT_0238 + -IT_0241 + IT_0242 + IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = IT_0069*IT_0123;
    const ccomplex_t IT_0247 = IT_0050*IT_0239*IT_0246;
    const ccomplex_t IT_0248 = IT_0046*IT_0136;
    const ccomplex_t IT_0249 = IT_0021*IT_0229*IT_0248;
    const ccomplex_t IT_0250 = IT_0092*IT_0149;
    const ccomplex_t IT_0251 = IT_0073*IT_0236*IT_0250;
    const ccomplex_t IT_0252 = IT_0115*IT_0162;
    const ccomplex_t IT_0253 = IT_0096*IT_0232*IT_0252;
    const ccomplex_t IT_0254 = m_c*IT_0072;
    const ccomplex_t IT_0255 = m_c*IT_0049;
    const ccomplex_t IT_0256 = m_c*IT_0095;
    const ccomplex_t IT_0257 = m_c*IT_0118;
    const ccomplex_t IT_0258 = (-0.166666666666667)*IT_0201 + -IT_0247 + 
      -IT_0249 + -IT_0251 + -IT_0253 + IT_0254 + IT_0255 + IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = 3*IT_0119;
    const ccomplex_t IT_0260 = IT_0020 + IT_0259;
    const ccomplex_t IT_0261 = 3*conj(IT_0119);
    const ccomplex_t IT_0262 = conj(IT_0020) + IT_0261;
    const ccomplex_t IT_0263 = (18*conj(IT_0020) + 6*conj(IT_0119))*IT_0208 + 
      (18*IT_0020 + 6*IT_0119)*conj(IT_0208) + 6*conj(IT_0174)*(IT_0228 + 3
      *IT_0245) + 6*IT_0174*(conj(IT_0228) + 3*conj(IT_0245)) + 6*conj(IT_0258)
      *IT_0260 + 6*IT_0258*IT_0262;
    const ccomplex_t IT_0264 = m_b*m_c*IT_0002;
    const ccomplex_t IT_0265 = 6*s_34;
    const ccomplex_t IT_0266 = 18*s_34;
    const ccomplex_t IT_0267 = 12*IT_0177;
    const ccomplex_t IT_0268 = (-6)*IT_0175;
    const ccomplex_t IT_0269 = m_b*s_13;
    const ccomplex_t IT_0270 = (-18)*conj(IT_0258);
    const ccomplex_t IT_0271 = (-6)*conj(IT_0208);
    const ccomplex_t IT_0272 = IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = (-18)*IT_0258;
    const ccomplex_t IT_0274 = (-6)*conj(IT_0245);
    const ccomplex_t IT_0275 = ((-18)*conj(IT_0020) + (-6)*conj(IT_0119))
      *IT_0228 + ((-18)*IT_0020 + (-6)*IT_0119)*conj(IT_0228) + (-6)*IT_0245
      *IT_0262 + IT_0174*IT_0272 + (-6)*conj(IT_0174)*(IT_0208 + (
      -0.166666666666667)*IT_0273) + IT_0260*IT_0274;
    const ccomplex_t IT_0276 = m_b*m_c;
    const ccomplex_t IT_0277 = (IT_0020*conj(IT_0020) + IT_0119*conj(IT_0119) 
      + IT_0174*conj(IT_0174))*IT_0179 + IT_0180*IT_0263 + (-6)*(conj(IT_0174)
      *IT_0260 + IT_0174*IT_0262)*IT_0264 + (conj(IT_0228)*IT_0245 + IT_0228
      *conj(IT_0245) + conj(IT_0208)*IT_0258 + IT_0208*conj(IT_0258))*IT_0265 + 
      (IT_0208*conj(IT_0208) + IT_0228*conj(IT_0228) + IT_0245*conj(IT_0245) +
       IT_0258*conj(IT_0258))*IT_0266 + (conj(IT_0020)*IT_0119 + IT_0020*conj
      (IT_0119))*(IT_0267 + IT_0268) + IT_0269*IT_0275 + (-6)*(conj(IT_0228)
      *IT_0258 + -1./6*IT_0228*((-18)*conj(IT_0208) + (-6)*conj(IT_0258)) + -1.
      /6*IT_0245*IT_0272 + -1./6*conj(IT_0245)*IT_0273 + -1./6*IT_0208*((-18)
      *conj(IT_0228) + IT_0274))*IT_0276;
    return create_ccomplex_return(IT_0277);
}

