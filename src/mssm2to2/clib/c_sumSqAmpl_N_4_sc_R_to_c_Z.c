#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_sc_R_to_c_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_sc_R_to_c_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0015 = 0.666666666666667*IT_0014;
    const ccomplex_t IT_0016 = pow(m_c, 2);
    const ccomplex_t IT_0017 = pow(m_sc_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0012;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0026 = IT_0003*IT_0025;
    const ccomplex_t IT_0027 = IT_0007*IT_0025;
    const ccomplex_t IT_0028 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0029 = IT_0007*IT_0028;
    const ccomplex_t IT_0030 = IT_0003*IT_0028;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0026 + IT_0027 + 
      -IT_0029 + -IT_0030);
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0037 = 0.666666666666667*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0035*IT_0040;
    const ccomplex_t IT_0042 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0043 = IT_0003*IT_0042;
    const ccomplex_t IT_0044 = IT_0007*IT_0042;
    const ccomplex_t IT_0045 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0046 = IT_0003*IT_0045;
    const ccomplex_t IT_0047 = IT_0007*IT_0045;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0043 + IT_0044 + 
      -IT_0046 + -IT_0047);
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0051 = 0.666666666666667*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0049*IT_0054;
    const ccomplex_t IT_0056 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0057 = IT_0003*IT_0056;
    const ccomplex_t IT_0058 = IT_0007*IT_0056;
    const ccomplex_t IT_0059 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0060 = IT_0003*IT_0059;
    const ccomplex_t IT_0061 = IT_0007*IT_0059;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0057 + IT_0058 + 
      -IT_0060 + -IT_0061);
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0065 = 0.666666666666667*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0063*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0024 + -IT_0041 + -IT_0055 + -IT_0069;
    const ccomplex_t IT_0071 = pow(m_W, -1);
    const ccomplex_t IT_0072 = sin(beta);
    const ccomplex_t IT_0073 = cpow(IT_0072, -1);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0006*IT_0071*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0018*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = IT_0023*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = IT_0013*IT_0077;
    const ccomplex_t IT_0081 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0082 = IT_0003*IT_0081;
    const ccomplex_t IT_0083 = IT_0007*IT_0081;
    const ccomplex_t IT_0084 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0085 = IT_0003*IT_0084;
    const ccomplex_t IT_0086 = IT_0007*IT_0084;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0082 + IT_0083 + 
      -IT_0085 + -IT_0086);
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0006*IT_0071*IT_0073;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = IT_0066*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0088*IT_0092;
    const ccomplex_t IT_0094 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = IT_0007*IT_0094;
    const ccomplex_t IT_0097 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0098 = IT_0003*IT_0097;
    const ccomplex_t IT_0099 = IT_0007*IT_0097;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0095 + IT_0096 + 
      -IT_0098 + -IT_0099);
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0006*IT_0071*IT_0073;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0052*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0101*IT_0105;
    const ccomplex_t IT_0107 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0108 = IT_0003*IT_0107;
    const ccomplex_t IT_0109 = IT_0007*IT_0107;
    const ccomplex_t IT_0110 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0111 = IT_0007*IT_0110;
    const ccomplex_t IT_0112 = IT_0003*IT_0110;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0108 + IT_0109 + 
      -IT_0111 + -IT_0112);
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0006*IT_0071*IT_0073;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = IT_0038*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = IT_0115*IT_0119;
    const ccomplex_t IT_0121 = 0.5*IT_0080 + IT_0093 + IT_0106 + IT_0120;
    const ccomplex_t IT_0122 = pow(m_Z, -2);
    const ccomplex_t IT_0123 = pow(m_Z, 4);
    const ccomplex_t IT_0124 = s_13*IT_0123;
    const ccomplex_t IT_0125 = IT_0122*IT_0124;
    const ccomplex_t IT_0126 = (-6)*IT_0125;
    const ccomplex_t IT_0127 = s_14*s_34;
    const ccomplex_t IT_0128 = (-24)*IT_0127;
    const ccomplex_t IT_0129 = pow(m_Z, 2);
    const ccomplex_t IT_0130 = s_13*IT_0129;
    const ccomplex_t IT_0131 = 12*IT_0130;
    const ccomplex_t IT_0132 = IT_0126 + IT_0128 + IT_0131;
    const ccomplex_t IT_0133 = 6*IT_0125;
    const ccomplex_t IT_0134 = 24*IT_0127;
    const ccomplex_t IT_0135 = (-12)*IT_0130;
    const ccomplex_t IT_0136 = IT_0133 + IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = (-2)*IT_0127;
    const ccomplex_t IT_0138 = IT_0130 + IT_0137;
    const ccomplex_t IT_0139 = IT_0122*IT_0138;
    const ccomplex_t IT_0140 = (-6)*IT_0139;
    const ccomplex_t IT_0141 = 12*s_13;
    const ccomplex_t IT_0142 = IT_0140 + IT_0141;
    const ccomplex_t IT_0143 = m_N_4*IT_0013;
    const ccomplex_t IT_0144 = m_N_4*IT_0023;
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = IT_0143 + IT_0145;
    const ccomplex_t IT_0147 = IT_0077*IT_0146;
    const ccomplex_t IT_0148 = m_N_4*IT_0101;
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = m_N_1*IT_0049;
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = IT_0105*IT_0151;
    const ccomplex_t IT_0153 = m_N_2*IT_0035;
    const ccomplex_t IT_0154 = m_N_4*IT_0115;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = IT_0119*IT_0156;
    const ccomplex_t IT_0158 = m_N_3*IT_0063;
    const ccomplex_t IT_0159 = m_N_4*IT_0088;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0158 + IT_0160;
    const ccomplex_t IT_0162 = IT_0092*IT_0161;
    const ccomplex_t IT_0163 = (-0.5)*IT_0147 + IT_0152 + IT_0157 + IT_0162;
    const ccomplex_t IT_0164 = pow(m_N_4, 2);
    const ccomplex_t IT_0165 = cpow((-2)*s_12 + IT_0016 + -IT_0017 + -IT_0164 
      + -reg_prop, -1);
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0167 = e_em*IT_0003;
    const ccomplex_t IT_0168 = e_em*IT_0007;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0167 + (-3)*IT_0168);
    const ccomplex_t IT_0170 = (-0.166666666666667)*IT_0169;
    const ccomplex_t IT_0171 = IT_0015*IT_0170;
    const ccomplex_t IT_0172 = IT_0165*IT_0166*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0174 = (-0.666666666666667)*IT_0173;
    const ccomplex_t IT_0175 = IT_0015*IT_0174;
    const ccomplex_t IT_0176 = IT_0165*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = m_c*IT_0177;
    const ccomplex_t IT_0179 = (-0.5)*IT_0147 + -IT_0172 + IT_0178;
    const ccomplex_t IT_0180 = IT_0163 + IT_0179;
    const ccomplex_t IT_0181 = conj(IT_0163) + conj(IT_0179);
    const ccomplex_t IT_0182 = IT_0075*IT_0174;
    const ccomplex_t IT_0183 = IT_0165*IT_0166*IT_0182;
    const ccomplex_t IT_0184 = IT_0075*IT_0170;
    const ccomplex_t IT_0185 = IT_0165*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = m_c*IT_0186;
    const ccomplex_t IT_0188 = IT_0020*IT_0146;
    const ccomplex_t IT_0189 = IT_0183 + -IT_0187 + (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = m_N_4*IT_0049;
    const ccomplex_t IT_0191 = -IT_0190;
    const ccomplex_t IT_0192 = m_N_1*IT_0101;
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = IT_0054*IT_0193;
    const ccomplex_t IT_0195 = m_N_2*IT_0115;
    const ccomplex_t IT_0196 = m_N_4*IT_0035;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = IT_0195 + IT_0197;
    const ccomplex_t IT_0199 = IT_0040*IT_0198;
    const ccomplex_t IT_0200 = m_N_3*IT_0088;
    const ccomplex_t IT_0201 = m_N_4*IT_0063;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = IT_0068*IT_0203;
    const ccomplex_t IT_0205 = (-0.5)*IT_0188 + -IT_0194 + -IT_0199 + -IT_0204;
    const ccomplex_t IT_0206 = IT_0189 + IT_0205;
    const ccomplex_t IT_0207 = conj(IT_0189) + conj(IT_0205);
    const ccomplex_t IT_0208 = m_c*m_N_4;
    const ccomplex_t IT_0209 = pow(s_14, 2);
    const ccomplex_t IT_0210 = IT_0122*IT_0209;
    const ccomplex_t IT_0211 = -IT_0164;
    const ccomplex_t IT_0212 = IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = IT_0208*IT_0212;
    const ccomplex_t IT_0214 = 1.33333333333333*IT_0173;
    const ccomplex_t IT_0215 = cpow((-2)*s_13 + IT_0016 + IT_0164 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0216 = IT_0075*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = IT_0214*IT_0217;
    const ccomplex_t IT_0219 = -IT_0078 + -IT_0218;
    const ccomplex_t IT_0220 = IT_0015*IT_0215;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = IT_0214*IT_0221;
    const ccomplex_t IT_0223 = IT_0021 + IT_0222;
    const ccomplex_t IT_0224 = 6*IT_0223;
    const ccomplex_t IT_0225 = IT_0024 + 2*IT_0041 + 2*IT_0055 + 2*IT_0069;
    const ccomplex_t IT_0226 = (-6)*IT_0225;
    const ccomplex_t IT_0227 = IT_0224 + IT_0226;
    const ccomplex_t IT_0228 = 6*conj(IT_0223);
    const ccomplex_t IT_0229 = (-6)*conj(IT_0225);
    const ccomplex_t IT_0230 = IT_0228 + IT_0229;
    const ccomplex_t IT_0231 = -IT_0080 + (-2)*IT_0093 + (-2)*IT_0106 + (-2)
      *IT_0120;
    const ccomplex_t IT_0232 = 6*IT_0225;
    const ccomplex_t IT_0233 = (-6)*IT_0223;
    const ccomplex_t IT_0234 = IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = 6*conj(IT_0225);
    const ccomplex_t IT_0236 = (-6)*conj(IT_0223);
    const ccomplex_t IT_0237 = IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = m_c*IT_0212;
    const ccomplex_t IT_0239 = 6*IT_0219;
    const ccomplex_t IT_0240 = (-6)*IT_0231;
    const ccomplex_t IT_0241 = IT_0239 + IT_0240;
    const ccomplex_t IT_0242 = 6*conj(IT_0219);
    const ccomplex_t IT_0243 = (-6)*conj(IT_0231);
    const ccomplex_t IT_0244 = IT_0242 + IT_0243;
    const ccomplex_t IT_0245 = s_14*s_34*IT_0122;
    const ccomplex_t IT_0246 = -IT_0245;
    const ccomplex_t IT_0247 = s_13 + IT_0246;
    const ccomplex_t IT_0248 = s_13*IT_0247;
    const ccomplex_t IT_0249 = (-2)*IT_0177 + -IT_0222;
    const ccomplex_t IT_0250 = 2*IT_0186 + IT_0218;
    const ccomplex_t IT_0251 = 6*IT_0231;
    const ccomplex_t IT_0252 = (-6)*IT_0219;
    const ccomplex_t IT_0253 = IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 6*conj(IT_0231);
    const ccomplex_t IT_0255 = (-6)*conj(IT_0219);
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = s_13*s_14;
    const ccomplex_t IT_0258 = s_34*IT_0164;
    const ccomplex_t IT_0259 = -IT_0258;
    const ccomplex_t IT_0260 = IT_0257 + IT_0259;
    const ccomplex_t IT_0261 = -IT_0177;
    const ccomplex_t IT_0262 = conj(IT_0070) + conj(IT_0261);
    const ccomplex_t IT_0263 = IT_0070 + IT_0261;
    const ccomplex_t IT_0264 = conj(IT_0121) + conj(IT_0186);
    const ccomplex_t IT_0265 = IT_0121 + IT_0186;
    const ccomplex_t IT_0266 = IT_0129*IT_0208;
    const ccomplex_t IT_0267 = 3*IT_0266;
    const ccomplex_t IT_0268 = 6*IT_0261;
    const ccomplex_t IT_0269 = 6*conj(IT_0261);
    const ccomplex_t IT_0270 = (-18)*conj(IT_0179);
    const ccomplex_t IT_0271 = s_13*IT_0212;
    const ccomplex_t IT_0272 = IT_0208*IT_0247;
    const ccomplex_t IT_0273 = 6*IT_0249;
    const ccomplex_t IT_0274 = 6*conj(IT_0249);
    const ccomplex_t IT_0275 = (-6)*IT_0249;
    const ccomplex_t IT_0276 = (-6)*conj(IT_0249);
    const ccomplex_t IT_0277 = pow(s_34, 2);
    const ccomplex_t IT_0278 = IT_0122*IT_0277;
    const ccomplex_t IT_0279 = -IT_0278;
    const ccomplex_t IT_0280 = IT_0016 + IT_0279;
    const ccomplex_t IT_0281 = IT_0208*IT_0280;
    const ccomplex_t IT_0282 = s_13*IT_0280;
    const ccomplex_t IT_0283 = (-6)*IT_0250;
    const ccomplex_t IT_0284 = m_c*IT_0247;
    const ccomplex_t IT_0285 = (-6)*conj(IT_0250);
    const ccomplex_t IT_0286 = s_14*IT_0016;
    const ccomplex_t IT_0287 = s_13*s_34;
    const ccomplex_t IT_0288 = -IT_0287;
    const ccomplex_t IT_0289 = IT_0286 + IT_0288;
    const ccomplex_t IT_0290 = conj(IT_0079)*IT_0250 + IT_0079*conj(IT_0250) +
       0.166666666666667*conj(IT_0022)*IT_0273 + 0.166666666666667*IT_0022
      *IT_0274 + 0.166666666666667*IT_0262*IT_0275 + 0.166666666666667*IT_0263
      *IT_0276 + 0.166666666666667*IT_0264*IT_0283 + 0.166666666666667*IT_0265
      *IT_0285;
    const ccomplex_t IT_0291 = m_N_4*IT_0247;
    const ccomplex_t IT_0292 = 6*IT_0179;
    const ccomplex_t IT_0293 = 6*conj(IT_0179);
    const ccomplex_t IT_0294 = (-6)*IT_0179;
    const ccomplex_t IT_0295 = (-6)*conj(IT_0179);
    const ccomplex_t IT_0296 = m_N_4*IT_0280;
    const ccomplex_t IT_0297 = s_34*IT_0129;
    const ccomplex_t IT_0298 = IT_0122*IT_0297;
    const ccomplex_t IT_0299 = (-0.25)*IT_0298;
    const ccomplex_t IT_0300 = s_34 + IT_0299;
    const ccomplex_t IT_0301 = m_N_4*IT_0300;
    const ccomplex_t IT_0302 = (-4)*IT_0301;
    const ccomplex_t IT_0303 = (-0.166666666666667)*IT_0294;
    const ccomplex_t IT_0304 = IT_0163 + IT_0303;
    const ccomplex_t IT_0305 = 6*conj(IT_0189);
    const ccomplex_t IT_0306 = 6*conj(IT_0205);
    const ccomplex_t IT_0307 = IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = 6*IT_0189;
    const ccomplex_t IT_0309 = 6*IT_0205;
    const ccomplex_t IT_0310 = IT_0308 + IT_0309;
    const ccomplex_t IT_0311 = (-0.166666666666667)*conj(IT_0022);
    const ccomplex_t IT_0312 = 6*conj(IT_0163);
    const ccomplex_t IT_0313 = IT_0293 + IT_0312;
    const ccomplex_t IT_0314 = (-0.166666666666667)*IT_0079;
    const ccomplex_t IT_0315 = 6*IT_0163;
    const ccomplex_t IT_0316 = IT_0292 + IT_0315;
    const ccomplex_t IT_0317 = (-0.166666666666667)*conj(IT_0079);
    const ccomplex_t IT_0318 = (-6)*conj(IT_0189);
    const ccomplex_t IT_0319 = (-6)*conj(IT_0205);
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = (-0.166666666666667)*IT_0070;
    const ccomplex_t IT_0322 = (-6)*IT_0189;
    const ccomplex_t IT_0323 = (-6)*IT_0205;
    const ccomplex_t IT_0324 = IT_0322 + IT_0323;
    const ccomplex_t IT_0325 = (-0.166666666666667)*conj(IT_0070);
    const ccomplex_t IT_0326 = (-6)*conj(IT_0163);
    const ccomplex_t IT_0327 = IT_0295 + IT_0326;
    const ccomplex_t IT_0328 = (-0.166666666666667)*IT_0121;
    const ccomplex_t IT_0329 = s_14*IT_0129;
    const ccomplex_t IT_0330 = IT_0122*IT_0329;
    const ccomplex_t IT_0331 = 0.5*IT_0330;
    const ccomplex_t IT_0332 = s_14 + IT_0331;
    const ccomplex_t IT_0333 = m_c*IT_0332;
    const ccomplex_t IT_0334 = 2*IT_0333;
    const ccomplex_t IT_0335 = (-0.166666666666667)*conj(IT_0121);
    const ccomplex_t IT_0336 = (-3)*IT_0266;
    const ccomplex_t IT_0337 = 6*IT_0022;
    const ccomplex_t IT_0338 = IT_0070 + (-0.166666666666667)*IT_0337;
    const ccomplex_t IT_0339 = 6*conj(IT_0022);
    const ccomplex_t IT_0340 = (-0.166666666666667)*IT_0339;
    const ccomplex_t IT_0341 = conj(IT_0070) + IT_0340;
    const ccomplex_t IT_0342 = (-0.166666666666667)*conj(IT_0186);
    const ccomplex_t IT_0343 = IT_0128 + IT_0133;
    const ccomplex_t IT_0344 = IT_0126 + IT_0134;
    const ccomplex_t IT_0345 = (-0.25)*IT_0330;
    const ccomplex_t IT_0346 = s_14 + IT_0345;
    const ccomplex_t IT_0347 = m_c*IT_0346;
    const ccomplex_t IT_0348 = (-4)*IT_0347;
    const ccomplex_t IT_0349 = (-0.166666666666667)*IT_0348;
    const ccomplex_t IT_0350 = 0.5*IT_0298;
    const ccomplex_t IT_0351 = s_34 + IT_0350;
    const ccomplex_t IT_0352 = m_N_4*IT_0351;
    const ccomplex_t IT_0353 = 2*IT_0352;
    const ccomplex_t IT_0354 = (-0.166666666666667)*IT_0353;
    const ccomplex_t IT_0355 = (conj(IT_0022)*IT_0070 + IT_0022*conj(IT_0070) 
      + conj(IT_0079)*IT_0121 + IT_0079*conj(IT_0121))*IT_0132 + (IT_0022*conj
      (IT_0022) + IT_0070*conj(IT_0070) + IT_0079*conj(IT_0079) + IT_0121*conj
      (IT_0121))*IT_0136 + IT_0142*(IT_0180*IT_0181 + IT_0206*IT_0207) + IT_0213
      *(conj(IT_0219)*IT_0227 + IT_0219*IT_0230 + conj(IT_0231)*IT_0234 +
       IT_0231*IT_0237) + IT_0238*(IT_0181*IT_0227 + IT_0180*IT_0230 + IT_0207
      *IT_0241 + IT_0206*IT_0244) + IT_0248*(IT_0237*IT_0249 + IT_0234*conj
      (IT_0249) + conj(IT_0250)*IT_0253 + IT_0250*IT_0256) + IT_0260*(conj
      (IT_0022)*IT_0227 + IT_0022*IT_0230 + conj(IT_0079)*IT_0241 + IT_0079
      *IT_0244 + IT_0234*IT_0262 + IT_0237*IT_0263 + IT_0253*IT_0264 + IT_0256
      *IT_0265) + IT_0267*(conj(IT_0079)*IT_0268 + IT_0079*IT_0269) + IT_0208*
      (conj(IT_0163)*((-18)*IT_0189 + (-18)*IT_0205) + IT_0180*((-18)*conj
      (IT_0189) + (-18)*conj(IT_0205)) + IT_0206*IT_0270) + (conj(IT_0223)
      *IT_0224 + conj(IT_0225)*IT_0234 + IT_0225*IT_0236 + conj(IT_0219)*IT_0239
       + conj(IT_0231)*IT_0253 + IT_0231*IT_0255)*IT_0271 + IT_0272*(IT_0237
      *IT_0250 + IT_0234*conj(IT_0250) + conj(IT_0231)*IT_0273 + IT_0231*IT_0274
       + conj(IT_0219)*IT_0275 + IT_0219*IT_0276) + (conj(IT_0250)*IT_0275 +
       IT_0250*IT_0276)*IT_0281 + IT_0282*(conj(IT_0249)*IT_0275 + conj(IT_0250)
      *IT_0283) + IT_0284*(IT_0181*IT_0275 + IT_0180*IT_0276 + IT_0207*IT_0283 +
       IT_0206*IT_0285) + 6*IT_0289*IT_0290 + IT_0291*(IT_0207*IT_0234 + IT_0206
      *IT_0237 + conj(IT_0163)*IT_0253 + IT_0163*IT_0256 + conj(IT_0231)*IT_0292
       + IT_0231*IT_0293 + conj(IT_0219)*IT_0294 + IT_0219*IT_0295) + (IT_0207
      *IT_0275 + IT_0206*IT_0276 + conj(IT_0163)*IT_0283 + IT_0163*IT_0285 +
       conj(IT_0250)*IT_0294 + IT_0250*IT_0295)*IT_0296 + (-6)*IT_0302*(conj
      (IT_0121)*IT_0304 + (-0.166666666666667)*IT_0022*IT_0307 + IT_0310*IT_0311
       + IT_0313*IT_0314 + IT_0316*IT_0317 + IT_0320*IT_0321 + IT_0324*IT_0325 +
       IT_0327*IT_0328) + (-6)*IT_0334*(conj(IT_0070)*IT_0304 + (
      -0.166666666666667)*IT_0022*IT_0313 + IT_0307*IT_0314 + IT_0311*IT_0316 +
       IT_0310*IT_0317 + IT_0321*IT_0327 + IT_0320*IT_0328 + IT_0324*IT_0335) + 
      (-6)*IT_0336*((conj(IT_0022) + -conj(IT_0070))*IT_0121 + (IT_0022 + 
      -IT_0070)*conj(IT_0121) + (-0.166666666666667)*IT_0186*IT_0269 + conj
      (IT_0079)*IT_0338 + IT_0079*IT_0341 + IT_0268*IT_0342) + (-6)*IT_0261*
      (conj(IT_0121)*IT_0267 + IT_0311*IT_0343 + IT_0325*IT_0344 + IT_0313
      *IT_0349 + IT_0307*IT_0354) + (-6)*conj(IT_0261)*((-0.166666666666667)
      *IT_0136*IT_0261 + IT_0121*IT_0267 + (-0.166666666666667)*IT_0022*IT_0343 
      + IT_0321*IT_0344 + IT_0316*IT_0349 + IT_0310*IT_0354) + (-6)*IT_0186*
      (IT_0267*IT_0341 + IT_0136*IT_0342 + IT_0317*IT_0343 + IT_0335*IT_0344 +
       IT_0307*IT_0349 + IT_0313*IT_0354) + (-6)*conj(IT_0186)*(IT_0267*IT_0338 
      + IT_0314*IT_0343 + IT_0328*IT_0344 + IT_0310*IT_0349 + IT_0316*IT_0354);
    return create_ccomplex_return(IT_0355);
}

