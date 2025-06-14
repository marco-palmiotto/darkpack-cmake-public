#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_st_2_to_c_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_st_2_to_c_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_st_2 = param->m_st_2;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_c, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sc_L, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_13 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = conj(N_B1)*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = sin(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = conj(N_W1)*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0015 + 3*IT_0020);
    const ccomplex_t IT_0022 = 0.166666666666667*IT_0021;
    const ccomplex_t IT_0023 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0024 = IT_0012*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = sin(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = IT_0017*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -1);
    const ccomplex_t IT_0030 = m_t*N_u1*e_em*IT_0029*U_st_01;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0025 + (-0.75)*IT_0032);
    const ccomplex_t IT_0034 = 0.666666666666667*IT_0033;
    const ccomplex_t IT_0035 = IT_0022*IT_0034;
    const ccomplex_t IT_0036 = IT_0010*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0039 = conj(N_B2)*e_em;
    const ccomplex_t IT_0040 = IT_0012*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = conj(N_W2)*e_em;
    const ccomplex_t IT_0043 = IT_0017*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0041 + 3*IT_0044);
    const ccomplex_t IT_0046 = 0.166666666666667*IT_0045;
    const ccomplex_t IT_0047 = N_B2*e_em*U_st_11;
    const ccomplex_t IT_0048 = IT_0012*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = m_t*N_u2*e_em*IT_0029*U_st_01;
    const ccomplex_t IT_0051 = IT_0028*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-0.75)*IT_0052);
    const ccomplex_t IT_0054 = 0.666666666666667*IT_0053;
    const ccomplex_t IT_0055 = IT_0046*IT_0054;
    const ccomplex_t IT_0056 = IT_0038*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0059 = conj(N_B3)*e_em;
    const ccomplex_t IT_0060 = IT_0012*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W3)*e_em;
    const ccomplex_t IT_0063 = IT_0017*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + 3*IT_0064);
    const ccomplex_t IT_0066 = 0.166666666666667*IT_0065;
    const ccomplex_t IT_0067 = N_B3*e_em*U_st_11;
    const ccomplex_t IT_0068 = IT_0012*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_t*N_u3*e_em*IT_0029*U_st_01;
    const ccomplex_t IT_0071 = IT_0028*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + (-0.75)*IT_0072);
    const ccomplex_t IT_0074 = 0.666666666666667*IT_0073;
    const ccomplex_t IT_0075 = IT_0066*IT_0074;
    const ccomplex_t IT_0076 = IT_0058*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_13 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0079 = conj(N_B4)*e_em;
    const ccomplex_t IT_0080 = IT_0012*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W4)*e_em;
    const ccomplex_t IT_0083 = IT_0017*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + 3*IT_0084);
    const ccomplex_t IT_0086 = 0.166666666666667*IT_0085;
    const ccomplex_t IT_0087 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0088 = IT_0012*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = m_t*N_u4*e_em*IT_0029*U_st_01;
    const ccomplex_t IT_0091 = IT_0028*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0089 + (-0.75)*IT_0092);
    const ccomplex_t IT_0094 = 0.666666666666667*IT_0093;
    const ccomplex_t IT_0095 = IT_0086*IT_0094;
    const ccomplex_t IT_0096 = IT_0078*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0008 + -IT_0037 + -IT_0057
       + -IT_0077 + -IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0017*IT_0027*IT_0029;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = conj(N_B3)*e_em*U_st_01;
    const ccomplex_t IT_0102 = IT_0012*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = conj(N_W3)*e_em*U_st_01;
    const ccomplex_t IT_0105 = IT_0017*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = m_t*conj(N_u3)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0108 = IT_0028*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + 3*IT_0106 + 3
      *IT_0109);
    const ccomplex_t IT_0111 = 0.166666666666667*IT_0110;
    const ccomplex_t IT_0112 = IT_0100*IT_0111;
    const ccomplex_t IT_0113 = IT_0058*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0017*IT_0027*IT_0029;
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0118 = IT_0012*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = m_t*conj(N_u4)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0124 = IT_0028*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0119 + 3*IT_0122 + 3
      *IT_0125);
    const ccomplex_t IT_0127 = 0.166666666666667*IT_0126;
    const ccomplex_t IT_0128 = IT_0116*IT_0127;
    const ccomplex_t IT_0129 = IT_0078*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0017*IT_0027*IT_0029;
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0134 = IT_0012*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0137 = IT_0017*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = m_t*conj(N_u1)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0140 = IT_0028*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0135 + 3*IT_0138 + 3
      *IT_0141);
    const ccomplex_t IT_0143 = 0.166666666666667*IT_0142;
    const ccomplex_t IT_0144 = IT_0132*IT_0143;
    const ccomplex_t IT_0145 = IT_0010*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0017*IT_0027*IT_0029;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = conj(N_B2)*e_em*U_st_01;
    const ccomplex_t IT_0150 = IT_0012*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W2)*e_em*U_st_01;
    const ccomplex_t IT_0153 = IT_0017*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = m_t*conj(N_u2)*e_em*IT_0029*U_st_11;
    const ccomplex_t IT_0156 = IT_0028*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0151 + 3*IT_0154 + 3
      *IT_0157);
    const ccomplex_t IT_0159 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0160 = IT_0148*IT_0159;
    const ccomplex_t IT_0161 = IT_0038*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = -IT_0114 + -IT_0130 + -IT_0146 + -IT_0162;
    const ccomplex_t IT_0164 = pow(m_st_2, 2);
    const ccomplex_t IT_0165 = m_c*m_t*IT_0164;
    const ccomplex_t IT_0166 = s_34*IT_0164;
    const ccomplex_t IT_0167 = (-18)*IT_0166;
    const ccomplex_t IT_0168 = s_23*s_24;
    const ccomplex_t IT_0169 = 36*IT_0168;
    const ccomplex_t IT_0170 = IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = m_c*s_24;
    const ccomplex_t IT_0172 = m_t*IT_0008;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0175 = IT_0000*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0177 = IT_0006*IT_0175*IT_0176;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = m_t*IT_0037;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0181 = IT_0066*IT_0111;
    const ccomplex_t IT_0182 = IT_0058*IT_0180*IT_0181;
    const ccomplex_t IT_0183 = m_t*IT_0057;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0185 = IT_0022*IT_0143;
    const ccomplex_t IT_0186 = IT_0010*IT_0184*IT_0185;
    const ccomplex_t IT_0187 = m_t*IT_0097;
    const ccomplex_t IT_0188 = m_t*IT_0077;
    const ccomplex_t IT_0189 = 0.166666666666667*IT_0172 + -IT_0179 + IT_0182 
      + -IT_0183 + IT_0186 + -IT_0187 + -IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0191 = IT_0046*IT_0159;
    const ccomplex_t IT_0192 = IT_0038*IT_0190*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0194 = IT_0086*IT_0127;
    const ccomplex_t IT_0195 = IT_0078*IT_0193*IT_0194;
    const ccomplex_t IT_0196 = 0.166666666666667*IT_0177 + -IT_0192 + -IT_0195;
    const ccomplex_t IT_0197 = IT_0173 + -IT_0178 + 3*IT_0189 + (-3)*IT_0196;
    const ccomplex_t IT_0198 = (-3)*conj(IT_0196);
    const ccomplex_t IT_0199 = 3*conj(IT_0189);
    const ccomplex_t IT_0200 = -conj(IT_0178);
    const ccomplex_t IT_0201 = conj(IT_0173) + IT_0198 + IT_0199 + IT_0200;
    const ccomplex_t IT_0202 = m_t*IT_0130;
    const ccomplex_t IT_0203 = IT_0054*IT_0148;
    const ccomplex_t IT_0204 = IT_0038*IT_0190*IT_0203;
    const ccomplex_t IT_0205 = m_t*IT_0146;
    const ccomplex_t IT_0206 = m_t*IT_0162;
    const ccomplex_t IT_0207 = m_t*IT_0114;
    const ccomplex_t IT_0208 = -IT_0202 + IT_0204 + -IT_0205 + -IT_0206 + 
      -IT_0207;
    const ccomplex_t IT_0209 = 6*IT_0208;
    const ccomplex_t IT_0210 = IT_0034*IT_0132;
    const ccomplex_t IT_0211 = IT_0010*IT_0184*IT_0210;
    const ccomplex_t IT_0212 = IT_0094*IT_0116;
    const ccomplex_t IT_0213 = IT_0078*IT_0193*IT_0212;
    const ccomplex_t IT_0214 = IT_0074*IT_0100;
    const ccomplex_t IT_0215 = IT_0058*IT_0180*IT_0214;
    const ccomplex_t IT_0216 = -IT_0211 + -IT_0213 + -IT_0215;
    const ccomplex_t IT_0217 = (-6)*IT_0216;
    const ccomplex_t IT_0218 = IT_0209 + IT_0217;
    const ccomplex_t IT_0219 = 6*conj(IT_0208);
    const ccomplex_t IT_0220 = (-6)*conj(IT_0216);
    const ccomplex_t IT_0221 = IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = 18*IT_0208;
    const ccomplex_t IT_0223 = (-18)*IT_0216;
    const ccomplex_t IT_0224 = IT_0222 + IT_0223;
    const ccomplex_t IT_0225 = 18*conj(IT_0208);
    const ccomplex_t IT_0226 = (-18)*conj(IT_0216);
    const ccomplex_t IT_0227 = IT_0225 + IT_0226;
    const ccomplex_t IT_0228 = m_c*m_t;
    const ccomplex_t IT_0229 = 18*IT_0216;
    const ccomplex_t IT_0230 = (-18)*IT_0208;
    const ccomplex_t IT_0231 = IT_0229 + IT_0230;
    const ccomplex_t IT_0232 = 18*conj(IT_0216);
    const ccomplex_t IT_0233 = (-18)*conj(IT_0208);
    const ccomplex_t IT_0234 = IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = conj(IT_0173)*(IT_0208 + -IT_0216) + IT_0173*
      (conj(IT_0208) + -conj(IT_0216)) + (-0.166666666666667)*conj(IT_0178)
      *IT_0218 + (-0.166666666666667)*IT_0178*IT_0221 + (-0.166666666666667)
      *conj(IT_0196)*IT_0224 + (-0.166666666666667)*IT_0196*IT_0227 + (
      -0.166666666666667)*conj(IT_0189)*IT_0231 + (-0.166666666666667)*IT_0189
      *IT_0234;
    const ccomplex_t IT_0236 = m_t*s_23;
    const ccomplex_t IT_0237 = conj(IT_0009)*((-18)*IT_0173 + 18*IT_0178 + (-6
      )*IT_0189 + 6*IT_0196) + IT_0009*((-18)*conj(IT_0173) + 18*conj(IT_0178) +
       (-6)*conj(IT_0189) + 6*conj(IT_0196)) + (-6)*conj(IT_0098)*IT_0197 + (-6)
      *IT_0098*IT_0201 + conj(IT_0163)*IT_0231 + IT_0163*IT_0234;
    const ccomplex_t IT_0238 = 6*s_34;
    const ccomplex_t IT_0239 = 18*s_34;
    const ccomplex_t IT_0240 = (-18)*s_34;
    const ccomplex_t IT_0241 = (-6)*s_34;
    const ccomplex_t IT_0242 = 12*IT_0168;
    const ccomplex_t IT_0243 = (-6)*IT_0166;
    const ccomplex_t IT_0244 = (-6)*((conj(IT_0009) + 3*conj(IT_0098))*IT_0163
       + (IT_0009 + 3*IT_0098)*conj(IT_0163))*IT_0165 + (IT_0009*conj(IT_0009) +
       IT_0098*conj(IT_0098) + IT_0163*conj(IT_0163))*IT_0170 + 6*IT_0171*(conj
      (IT_0163)*IT_0197 + IT_0163*IT_0201 + 0.166666666666667*conj(IT_0009)
      *IT_0218 + 0.166666666666667*IT_0009*IT_0221 + 0.166666666666667*conj
      (IT_0098)*IT_0224 + 0.166666666666667*IT_0098*IT_0227) + (-6)*IT_0228
      *IT_0235 + IT_0236*IT_0237 + (conj(IT_0173)*IT_0189 + IT_0173*conj(IT_0189
      ) + conj(IT_0178)*IT_0196 + IT_0178*conj(IT_0196))*IT_0238 + (IT_0173*conj
      (IT_0173) + IT_0178*conj(IT_0178) + IT_0189*conj(IT_0189) + IT_0196*conj
      (IT_0196) + IT_0208*conj(IT_0208) + IT_0216*conj(IT_0216))*IT_0239 + (conj
      (IT_0173)*IT_0178 + IT_0173*conj(IT_0178) + conj(IT_0189)*IT_0196 +
       IT_0189*conj(IT_0196) + conj(IT_0208)*IT_0216 + IT_0208*conj(IT_0216))
      *IT_0240 + (conj(IT_0178)*IT_0189 + IT_0178*conj(IT_0189) + conj(IT_0173)
      *IT_0196 + IT_0173*conj(IT_0196))*IT_0241 + (conj(IT_0009)*IT_0098 +
       IT_0009*conj(IT_0098))*(IT_0242 + IT_0243);
    return create_ccomplex_return(IT_0244);
}

