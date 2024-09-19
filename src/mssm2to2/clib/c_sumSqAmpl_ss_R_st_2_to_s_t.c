#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_st_2_to_s_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_st_2_to_s_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
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
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_st_2 = param->m_st_2;
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
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_st_2, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_23*s_24;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_s, 2);
    const ccomplex_t IT_0007 = pow(m_sG, 2);
    const ccomplex_t IT_0008 = pow(m_ss_R, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0006 + -IT_0007 + IT_0008 +
       reg_prop, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0006 + IT_0008 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = sin(theta_W);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = sin(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = IT_0022*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = m_t*conj(N_u1)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0020 + 3*IT_0025 + 3
      *IT_0032);
    const ccomplex_t IT_0034 = 0.166666666666667*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0017;
    const ccomplex_t IT_0036 = (-0.333333333333333)*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = IT_0015*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0006 + IT_0008 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0041 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0042 = IT_0017*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0045 = IT_0022*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = m_t*conj(N_u2)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0048 = IT_0028*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0043 + 3*IT_0046 + 3
      *IT_0049);
    const ccomplex_t IT_0051 = 0.166666666666667*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0017;
    const ccomplex_t IT_0053 = (-0.333333333333333)*IT_0052;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = IT_0040*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0017;
    const ccomplex_t IT_0058 = (-0.333333333333333)*IT_0057;
    const ccomplex_t IT_0059 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0060 = IT_0017*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0063 = IT_0022*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = m_t*conj(N_u3)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0066 = IT_0028*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0061 + 3*IT_0064 + 3
      *IT_0067);
    const ccomplex_t IT_0069 = 0.166666666666667*IT_0068;
    const ccomplex_t IT_0070 = IT_0058*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0006 + IT_0008 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0017;
    const ccomplex_t IT_0075 = (-0.333333333333333)*IT_0074;
    const ccomplex_t IT_0076 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0077 = IT_0017*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0080 = IT_0022*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = m_t*conj(N_u4)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0083 = IT_0028*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0078 + 3*IT_0081 + 3
      *IT_0084);
    const ccomplex_t IT_0086 = 0.166666666666667*IT_0085;
    const ccomplex_t IT_0087 = IT_0075*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0006 + IT_0008 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = 0.166666666666667*IT_0014 + -IT_0039 + -IT_0056
       + -IT_0073 + -IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0014;
    const ccomplex_t IT_0093 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0094 = cos(beta);
    const ccomplex_t IT_0095 = cpow(IT_0094, -1);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0022*IT_0029*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0022*IT_0029*IT_0095*U_st_01;
    const ccomplex_t IT_0099 = 0.5*IT_0098;
    const ccomplex_t IT_0100 = IT_0097*IT_0099;
    const ccomplex_t IT_0101 = IT_0093*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0022*IT_0029*IT_0095;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0022*IT_0029*IT_0095*U_st_01;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0102 + IT_0110;
    const ccomplex_t IT_0112 = IT_0092 + IT_0111;
    const ccomplex_t IT_0113 = conj(IT_0092) + conj(IT_0111);
    const ccomplex_t IT_0114 = (-18)*IT_0001;
    const ccomplex_t IT_0115 = 36*IT_0003;
    const ccomplex_t IT_0116 = IT_0114 + IT_0115;
    const ccomplex_t IT_0117 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0118 = IT_0017*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = m_t*N_u1*e_em*IT_0029*U_st_01;
    const ccomplex_t IT_0121 = IT_0028*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0119 + (-0.75)*IT_0122);
    const ccomplex_t IT_0124 = 0.666666666666667*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0022*IT_0029*IT_0095;
    const ccomplex_t IT_0126 = 0.5*IT_0125;
    const ccomplex_t IT_0127 = IT_0124*IT_0126;
    const ccomplex_t IT_0128 = IT_0015*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0131 = IT_0017*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = m_t*N_u2*e_em*IT_0029*U_st_01;
    const ccomplex_t IT_0134 = IT_0028*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + (-0.75)*IT_0135);
    const ccomplex_t IT_0137 = 0.666666666666667*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0022*IT_0029*IT_0095;
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = IT_0137*IT_0139;
    const ccomplex_t IT_0141 = IT_0040*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0022*IT_0029*IT_0095;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0146 = IT_0017*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = m_t*N_u3*e_em*IT_0029*U_st_01;
    const ccomplex_t IT_0149 = IT_0028*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + (-0.75)*IT_0150);
    const ccomplex_t IT_0152 = 0.666666666666667*IT_0151;
    const ccomplex_t IT_0153 = IT_0144*IT_0152;
    const ccomplex_t IT_0154 = IT_0071*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0022*IT_0029*IT_0095;
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = m_t*N_u4*e_em*IT_0029*U_st_01;
    const ccomplex_t IT_0159 = IT_0028*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0162 = IT_0017*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + (
      -1.33333333333333)*IT_0163);
    const ccomplex_t IT_0165 = (-0.75)*IT_0164;
    const ccomplex_t IT_0166 = 0.666666666666667*IT_0165;
    const ccomplex_t IT_0167 = IT_0157*IT_0166;
    const ccomplex_t IT_0168 = IT_0088*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = -IT_0129 + -IT_0142 + -IT_0155 + -IT_0169;
    const ccomplex_t IT_0171 = m_s*m_t;
    const ccomplex_t IT_0172 = m_t*IT_0014;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0175 = IT_0010*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0177 = IT_0009*IT_0175*IT_0176;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = IT_0058*IT_0152;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0181 = IT_0071*IT_0179*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0183 = IT_0053*IT_0137;
    const ccomplex_t IT_0184 = IT_0040*IT_0182*IT_0183;
    const ccomplex_t IT_0185 = m_t*IT_0039;
    const ccomplex_t IT_0186 = m_t*IT_0056;
    const ccomplex_t IT_0187 = m_t*IT_0073;
    const ccomplex_t IT_0188 = m_t*IT_0090;
    const ccomplex_t IT_0189 = 0.166666666666667*IT_0172 + IT_0181 + IT_0184 +
       -IT_0185 + -IT_0186 + -IT_0187 + -IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0191 = IT_0036*IT_0124;
    const ccomplex_t IT_0192 = IT_0015*IT_0190*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0194 = IT_0075*IT_0166;
    const ccomplex_t IT_0195 = IT_0088*IT_0193*IT_0194;
    const ccomplex_t IT_0196 = 0.166666666666667*IT_0177 + -IT_0192 + -IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0198 = conj(V_ts)*e_em*conj(V_Wp1)*U_st_01;
    const ccomplex_t IT_0199 = IT_0022*IT_0198;
    const ccomplex_t IT_0200 = m_t*conj(V_ts)*conj(V_u1)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0201 = IT_0028*IT_0200;
    const ccomplex_t IT_0202 = 1.4142135623731*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*(IT_0199 + (-0.5)*IT_0202);
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = IT_0097*IT_0204;
    const ccomplex_t IT_0206 = IT_0093*IT_0197*IT_0205;
    const ccomplex_t IT_0207 = m_s*IT_0102;
    const ccomplex_t IT_0208 = m_s*IT_0110;
    const ccomplex_t IT_0209 = IT_0206 + IT_0207 + IT_0208;
    const ccomplex_t IT_0210 = IT_0069*IT_0144;
    const ccomplex_t IT_0211 = IT_0071*IT_0180*IT_0210;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = 6*IT_0173;
    const ccomplex_t IT_0214 = 18*IT_0189;
    const ccomplex_t IT_0215 = (-18)*IT_0196;
    const ccomplex_t IT_0216 = (-6)*IT_0178;
    const ccomplex_t IT_0217 = IT_0213 + IT_0214 + IT_0215 + IT_0216;
    const ccomplex_t IT_0218 = 6*conj(IT_0173);
    const ccomplex_t IT_0219 = 18*conj(IT_0189);
    const ccomplex_t IT_0220 = (-18)*conj(IT_0196);
    const ccomplex_t IT_0221 = (-6)*conj(IT_0178);
    const ccomplex_t IT_0222 = IT_0218 + IT_0219 + IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = IT_0051*IT_0139;
    const ccomplex_t IT_0224 = IT_0040*IT_0182*IT_0223;
    const ccomplex_t IT_0225 = IT_0034*IT_0126;
    const ccomplex_t IT_0226 = IT_0015*IT_0190*IT_0225;
    const ccomplex_t IT_0227 = IT_0086*IT_0157;
    const ccomplex_t IT_0228 = IT_0088*IT_0193*IT_0227;
    const ccomplex_t IT_0229 = m_t*IT_0129;
    const ccomplex_t IT_0230 = m_t*IT_0142;
    const ccomplex_t IT_0231 = m_t*IT_0155;
    const ccomplex_t IT_0232 = m_t*IT_0169;
    const ccomplex_t IT_0233 = IT_0224 + IT_0226 + IT_0228 + -IT_0229 + 
      -IT_0230 + -IT_0231 + -IT_0232;
    const ccomplex_t IT_0234 = 6*IT_0178;
    const ccomplex_t IT_0235 = 18*IT_0196;
    const ccomplex_t IT_0236 = (-18)*IT_0189;
    const ccomplex_t IT_0237 = (-6)*IT_0173;
    const ccomplex_t IT_0238 = IT_0234 + IT_0235 + IT_0236 + IT_0237;
    const ccomplex_t IT_0239 = 6*conj(IT_0178);
    const ccomplex_t IT_0240 = 18*conj(IT_0196);
    const ccomplex_t IT_0241 = (-18)*conj(IT_0189);
    const ccomplex_t IT_0242 = (-6)*conj(IT_0173);
    const ccomplex_t IT_0243 = IT_0239 + IT_0240 + IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = conj(V_ts)*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0245 = IT_0022*IT_0244;
    const ccomplex_t IT_0246 = m_t*conj(V_ts)*conj(V_u2)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0247 = IT_0028*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*(IT_0245 + (-0.5)*IT_0248);
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = IT_0104*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0253 = IT_0108*IT_0251*IT_0252;
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = 6*IT_0196;
    const ccomplex_t IT_0256 = 18*IT_0178;
    const ccomplex_t IT_0257 = (-18)*IT_0173;
    const ccomplex_t IT_0258 = (-6)*IT_0189;
    const ccomplex_t IT_0259 = IT_0255 + IT_0256 + IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = 6*conj(IT_0196);
    const ccomplex_t IT_0261 = 18*conj(IT_0178);
    const ccomplex_t IT_0262 = (-18)*conj(IT_0173);
    const ccomplex_t IT_0263 = (-6)*conj(IT_0189);
    const ccomplex_t IT_0264 = IT_0260 + IT_0261 + IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = (18*conj(IT_0173) + (-18)*conj(IT_0178) + 6
      *conj(IT_0189) + (-6)*conj(IT_0196))*IT_0209 + (18*IT_0173 + (-18)*IT_0178
       + 6*IT_0189 + (-6)*IT_0196)*conj(IT_0209) + conj(IT_0212)*IT_0217 +
       IT_0212*IT_0222 + conj(IT_0233)*IT_0238 + IT_0233*IT_0243 + conj(IT_0254)
      *IT_0259 + IT_0254*IT_0264;
    const ccomplex_t IT_0266 = m_t*s_23;
    const ccomplex_t IT_0267 = -IT_0209 + (-3)*IT_0212 + 3*IT_0233 + IT_0254;
    const ccomplex_t IT_0268 = (-3)*conj(IT_0212);
    const ccomplex_t IT_0269 = -conj(IT_0209);
    const ccomplex_t IT_0270 = 3*conj(IT_0233) + conj(IT_0254) + IT_0268 +
       IT_0269;
    const ccomplex_t IT_0271 = 6*s_34;
    const ccomplex_t IT_0272 = 18*s_34;
    const ccomplex_t IT_0273 = (-18)*s_34;
    const ccomplex_t IT_0274 = (-6)*s_34;
    const ccomplex_t IT_0275 = m_s*s_24;
    const ccomplex_t IT_0276 = 6*IT_0233;
    const ccomplex_t IT_0277 = 18*IT_0254;
    const ccomplex_t IT_0278 = (-18)*IT_0209;
    const ccomplex_t IT_0279 = (-6)*IT_0212;
    const ccomplex_t IT_0280 = 6*conj(IT_0233);
    const ccomplex_t IT_0281 = 18*conj(IT_0254);
    const ccomplex_t IT_0282 = (-18)*conj(IT_0209);
    const ccomplex_t IT_0283 = (-6)*conj(IT_0212);
    const ccomplex_t IT_0284 = 0.166666666666667*conj(IT_0170)*IT_0217 +
       0.166666666666667*IT_0170*IT_0222 + conj(IT_0091)*IT_0267 + IT_0091
      *IT_0270 + 0.166666666666667*IT_0113*(IT_0276 + IT_0277 + IT_0278 +
       IT_0279) + 0.166666666666667*IT_0112*(IT_0280 + IT_0281 + IT_0282 +
       IT_0283);
    const ccomplex_t IT_0285 = m_s*m_t*IT_0000;
    const ccomplex_t IT_0286 = (-6)*IT_0170;
    const ccomplex_t IT_0287 = (-6)*conj(IT_0170);
    const ccomplex_t IT_0288 = IT_0005*(conj(IT_0091)*IT_0112 + IT_0091
      *IT_0113) + IT_0116*(IT_0091*conj(IT_0091) + IT_0112*IT_0113 + IT_0170
      *conj(IT_0170)) + IT_0171*IT_0265 + (-6)*IT_0266*((-0.166666666666667)
      *conj(IT_0091)*IT_0238 + (-0.166666666666667)*IT_0091*IT_0243 + (
      -0.166666666666667)*IT_0113*IT_0259 + (-0.166666666666667)*IT_0112*IT_0264
       + conj(IT_0170)*IT_0267 + IT_0170*IT_0270) + (conj(IT_0173)*IT_0189 +
       IT_0173*conj(IT_0189) + conj(IT_0178)*IT_0196 + IT_0178*conj(IT_0196) +
       conj(IT_0209)*IT_0212 + IT_0209*conj(IT_0212) + conj(IT_0233)*IT_0254 +
       IT_0233*conj(IT_0254))*IT_0271 + (IT_0173*conj(IT_0173) + IT_0178*conj
      (IT_0178) + IT_0189*conj(IT_0189) + IT_0196*conj(IT_0196) + IT_0209*conj
      (IT_0209) + IT_0212*conj(IT_0212) + IT_0233*conj(IT_0233) + IT_0254*conj
      (IT_0254))*IT_0272 + (conj(IT_0173)*IT_0178 + IT_0173*conj(IT_0178) + conj
      (IT_0189)*IT_0196 + IT_0189*conj(IT_0196) + conj(IT_0212)*IT_0233 +
       IT_0212*conj(IT_0233) + conj(IT_0209)*IT_0254 + IT_0209*conj(IT_0254))
      *IT_0273 + (conj(IT_0178)*IT_0189 + IT_0178*conj(IT_0189) + conj(IT_0173)
      *IT_0196 + IT_0173*conj(IT_0196) + conj(IT_0209)*IT_0233 + IT_0209*conj
      (IT_0233) + conj(IT_0212)*IT_0254 + IT_0212*conj(IT_0254))*IT_0274 + 6
      *IT_0275*IT_0284 + IT_0285*((-18)*conj(IT_0091)*IT_0170 + (-18)*IT_0091
      *conj(IT_0170) + IT_0113*IT_0286 + IT_0112*IT_0287);
    return create_ccomplex_return(IT_0288);
}

