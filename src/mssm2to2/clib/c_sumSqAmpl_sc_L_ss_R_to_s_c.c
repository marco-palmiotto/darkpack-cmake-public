#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_ss_R_to_s_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_ss_R_to_s_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_R = param->m_ss_R;
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
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_ss_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*IT_0005;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = pow(m_sc_L, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0000 + IT_0008 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_cs)*e_em*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = IT_0009*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0000 + IT_0008 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_cs*e_em*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_cs)*e_em*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = IT_0022*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0007 + IT_0021 + IT_0029;
    const ccomplex_t IT_0031 = cos(theta_W);
    const ccomplex_t IT_0032 = cpow(IT_0031, -1);
    const ccomplex_t IT_0033 = conj(N_B2)*e_em;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W2)*e_em;
    const ccomplex_t IT_0037 = IT_0014*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + 3*IT_0038);
    const ccomplex_t IT_0040 = 0.166666666666667*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0032;
    const ccomplex_t IT_0042 = (-0.333333333333333)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0048 = conj(N_B3)*e_em;
    const ccomplex_t IT_0049 = IT_0032*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W3)*e_em;
    const ccomplex_t IT_0052 = IT_0014*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + 3*IT_0053);
    const ccomplex_t IT_0055 = 0.166666666666667*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0032;
    const ccomplex_t IT_0057 = (-0.333333333333333)*IT_0056;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = conj(N_B1)*e_em;
    const ccomplex_t IT_0062 = IT_0032*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = conj(N_W1)*e_em;
    const ccomplex_t IT_0065 = IT_0014*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + 3*IT_0066);
    const ccomplex_t IT_0068 = 0.166666666666667*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0032;
    const ccomplex_t IT_0070 = (-0.333333333333333)*IT_0069;
    const ccomplex_t IT_0071 = IT_0068*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B4)*e_em;
    const ccomplex_t IT_0076 = IT_0032*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W4)*e_em;
    const ccomplex_t IT_0079 = IT_0014*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 3*IT_0080);
    const ccomplex_t IT_0082 = 0.166666666666667*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0032;
    const ccomplex_t IT_0084 = (-0.333333333333333)*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = 0.166666666666667*IT_0007 + -IT_0046 + -IT_0060
       + -IT_0074 + -IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = sin(beta);
    const ccomplex_t IT_0093 = cpow(IT_0092, -1);
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0010*IT_0014*IT_0093;
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0091*IT_0095;
    const ccomplex_t IT_0097 = IT_0072*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0010*IT_0014*IT_0093;
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = IT_0100*IT_0102;
    const ccomplex_t IT_0104 = IT_0086*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0010*IT_0014*IT_0093;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = IT_0044*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0010*IT_0012*IT_0014;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0010*IT_0014*IT_0093;
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = IT_0114*IT_0116;
    const ccomplex_t IT_0118 = IT_0047*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = -IT_0098 + -IT_0105 + -IT_0112 + -IT_0119;
    const ccomplex_t IT_0121 = m_c*m_s*IT_0002;
    const ccomplex_t IT_0122 = s_34*IT_0002;
    const ccomplex_t IT_0123 = (-18)*IT_0122;
    const ccomplex_t IT_0124 = s_23*s_24;
    const ccomplex_t IT_0125 = 36*IT_0124;
    const ccomplex_t IT_0126 = IT_0123 + IT_0125;
    const ccomplex_t IT_0127 = IT_0021 + IT_0029;
    const ccomplex_t IT_0128 = m_c*IT_0127;
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = IT_0070*IT_0095;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0132 = IT_0072*IT_0130*IT_0131;
    const ccomplex_t IT_0133 = IT_0084*IT_0102;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0135 = IT_0086*IT_0133*IT_0134;
    const ccomplex_t IT_0136 = m_s*IT_0098;
    const ccomplex_t IT_0137 = m_s*IT_0112;
    const ccomplex_t IT_0138 = m_s*IT_0119;
    const ccomplex_t IT_0139 = m_s*IT_0105;
    const ccomplex_t IT_0140 = IT_0132 + IT_0135 + -IT_0136 + -IT_0137 + 
      -IT_0138 + -IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0142 = IT_0040*IT_0107;
    const ccomplex_t IT_0143 = IT_0044*IT_0141*IT_0142;
    const ccomplex_t IT_0144 = IT_0082*IT_0100;
    const ccomplex_t IT_0145 = IT_0086*IT_0134*IT_0144;
    const ccomplex_t IT_0146 = IT_0068*IT_0091;
    const ccomplex_t IT_0147 = IT_0072*IT_0131*IT_0146;
    const ccomplex_t IT_0148 = -IT_0143 + -IT_0145 + -IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp2)
      *IT_0014;
    const ccomplex_t IT_0151 = IT_0024*IT_0150;
    const ccomplex_t IT_0152 = IT_0022*IT_0149*IT_0151;
    const ccomplex_t IT_0153 = -IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp1)
      *IT_0014;
    const ccomplex_t IT_0156 = IT_0016*IT_0155;
    const ccomplex_t IT_0157 = IT_0009*IT_0154*IT_0156;
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = IT_0153 + IT_0158;
    const ccomplex_t IT_0160 = conj(IT_0153) + conj(IT_0158);
    const ccomplex_t IT_0161 = 6*s_34;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0163 = IT_0055*IT_0114;
    const ccomplex_t IT_0164 = IT_0047*IT_0162*IT_0163;
    const ccomplex_t IT_0165 = m_s*IT_0088;
    const ccomplex_t IT_0166 = m_s*IT_0074;
    const ccomplex_t IT_0167 = m_s*IT_0060;
    const ccomplex_t IT_0168 = IT_0164 + -IT_0165 + -IT_0166 + -IT_0167;
    const ccomplex_t IT_0169 = (-6)*IT_0046;
    const ccomplex_t IT_0170 = IT_0007 + IT_0169;
    const ccomplex_t IT_0171 = m_s*IT_0170;
    const ccomplex_t IT_0172 = 0.166666666666667*IT_0171;
    const ccomplex_t IT_0173 = m_s*IT_0007;
    const ccomplex_t IT_0174 = (-0.5)*IT_0173;
    const ccomplex_t IT_0175 = IT_0161*IT_0174;
    const ccomplex_t IT_0176 = IT_0042*IT_0109;
    const ccomplex_t IT_0177 = IT_0044*IT_0141*IT_0176;
    const ccomplex_t IT_0178 = IT_0057*IT_0116;
    const ccomplex_t IT_0179 = IT_0047*IT_0162*IT_0178;
    const ccomplex_t IT_0180 = -IT_0177 + -IT_0179;
    const ccomplex_t IT_0181 = conj(IT_0168) + conj(IT_0172);
    const ccomplex_t IT_0182 = 18*s_34;
    const ccomplex_t IT_0183 = IT_0168*IT_0182;
    const ccomplex_t IT_0184 = IT_0168 + IT_0172;
    const ccomplex_t IT_0185 = (-18)*s_34;
    const ccomplex_t IT_0186 = (-6)*s_34;
    const ccomplex_t IT_0187 = IT_0174*IT_0186;
    const ccomplex_t IT_0188 = IT_0153*IT_0182;
    const ccomplex_t IT_0189 = IT_0168*IT_0186;
    const ccomplex_t IT_0190 = IT_0174*IT_0182;
    const ccomplex_t IT_0191 = conj(IT_0168)*IT_0186;
    const ccomplex_t IT_0192 = m_c*m_s;
    const ccomplex_t IT_0193 = 6*IT_0129;
    const ccomplex_t IT_0194 = 18*IT_0140;
    const ccomplex_t IT_0195 = (-18)*IT_0180;
    const ccomplex_t IT_0196 = IT_0193 + IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = 6*conj(IT_0129);
    const ccomplex_t IT_0198 = 18*conj(IT_0140);
    const ccomplex_t IT_0199 = (-18)*conj(IT_0180);
    const ccomplex_t IT_0200 = IT_0197 + IT_0198 + IT_0199;
    const ccomplex_t IT_0201 = 6*IT_0140;
    const ccomplex_t IT_0202 = 18*IT_0129;
    const ccomplex_t IT_0203 = (-6)*IT_0180;
    const ccomplex_t IT_0204 = IT_0201 + IT_0202 + IT_0203;
    const ccomplex_t IT_0205 = 6*conj(IT_0140);
    const ccomplex_t IT_0206 = 18*conj(IT_0129);
    const ccomplex_t IT_0207 = (-6)*conj(IT_0180);
    const ccomplex_t IT_0208 = IT_0205 + IT_0206 + IT_0207;
    const ccomplex_t IT_0209 = 18*IT_0180;
    const ccomplex_t IT_0210 = (-18)*IT_0140;
    const ccomplex_t IT_0211 = (-6)*IT_0129;
    const ccomplex_t IT_0212 = IT_0209 + IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = 18*conj(IT_0180);
    const ccomplex_t IT_0214 = (-18)*conj(IT_0140);
    const ccomplex_t IT_0215 = (-6)*conj(IT_0129);
    const ccomplex_t IT_0216 = IT_0213 + IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = conj(IT_0148)*IT_0196 + IT_0148*IT_0200 +
       IT_0174*((-18)*conj(IT_0129) + (-6)*conj(IT_0140) + 6*conj(IT_0180) +
       IT_0201 + IT_0202 + IT_0203) + IT_0160*IT_0204 + IT_0159*IT_0208 +
       IT_0181*IT_0212 + IT_0184*IT_0216;
    const ccomplex_t IT_0218 = m_s*s_24;
    const ccomplex_t IT_0219 = 6*IT_0174;
    const ccomplex_t IT_0220 = (-6)*IT_0174;
    const ccomplex_t IT_0221 = 6*IT_0030*(conj(IT_0148) + 3*conj(IT_0153) + 3
      *conj(IT_0158) + -conj(IT_0168) + -conj(IT_0172) + 3*IT_0174) + 6*conj
      (IT_0030)*(IT_0148 + 3*IT_0153 + 3*IT_0158 + -IT_0168 + -IT_0172 + (-3)
      *IT_0174) + conj(IT_0120)*IT_0212 + IT_0120*IT_0216 + IT_0089*(18*conj
      (IT_0148) + 6*conj(IT_0153) + 6*conj(IT_0158) + (-18)*conj(IT_0168) + (-18
      )*conj(IT_0172) + IT_0219) + conj(IT_0089)*(18*IT_0148 + 6*IT_0153 + 6
      *IT_0158 + (-18)*IT_0168 + (-18)*IT_0172 + IT_0220);
    const ccomplex_t IT_0222 = m_c*s_23;
    const ccomplex_t IT_0223 = conj(IT_0089)*IT_0196 + IT_0089*IT_0200 + conj
      (IT_0030)*IT_0204 + IT_0030*IT_0208 + conj(IT_0120)*((-18)*IT_0148 + (-6)
      *IT_0153 + (-6)*IT_0158 + 18*IT_0168 + 18*IT_0172 + IT_0219) + IT_0120*((
      -18)*conj(IT_0148) + (-6)*conj(IT_0153) + (-6)*conj(IT_0158) + 18*conj
      (IT_0168) + 18*conj(IT_0172) + IT_0220);
    const ccomplex_t IT_0224 = 12*IT_0124;
    const ccomplex_t IT_0225 = (-6)*IT_0122;
    const ccomplex_t IT_0226 = (-6)*((conj(IT_0030) + 3*conj(IT_0089))*IT_0120
       + (IT_0030 + 3*IT_0089)*conj(IT_0120))*IT_0121 + (IT_0030*conj(IT_0030) +
       IT_0089*conj(IT_0089) + IT_0120*conj(IT_0120))*IT_0126 + (conj(IT_0129)
      *IT_0140 + IT_0129*conj(IT_0140) + conj(IT_0148)*IT_0159 + IT_0148*IT_0160
      )*IT_0161 + (IT_0148 + conj(IT_0168) + conj(IT_0172))*IT_0175 + (IT_0129
      *conj(IT_0129) + IT_0140*conj(IT_0140) + IT_0148*conj(IT_0148) + IT_0158
      *IT_0160 + IT_0180*conj(IT_0180) + IT_0172*IT_0181)*IT_0182 + IT_0181
      *IT_0183 + ((conj(IT_0153) + conj(IT_0158) + IT_0174)*IT_0174 + conj
      (IT_0140)*IT_0180 + IT_0140*conj(IT_0180) + IT_0148*IT_0181 + conj(IT_0148
      )*IT_0184)*IT_0185 + (IT_0160*IT_0172 + IT_0159*conj(IT_0172) + conj
      (IT_0129)*IT_0180 + IT_0129*conj(IT_0180))*IT_0186 + (conj(IT_0148) +
       IT_0168 + IT_0172)*IT_0187 + IT_0160*(IT_0188 + IT_0189) + IT_0159*
      (IT_0190 + IT_0191) + IT_0192*IT_0217 + IT_0218*IT_0221 + IT_0222*IT_0223 
      + (conj(IT_0030)*IT_0089 + IT_0030*conj(IT_0089))*(IT_0224 + IT_0225);
    return create_ccomplex_return(IT_0226);
}

