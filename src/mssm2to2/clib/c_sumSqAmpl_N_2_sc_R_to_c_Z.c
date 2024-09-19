#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sc_R_to_c_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sc_R_to_c_Z(
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
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.33333333333333*IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0006 = 0.666666666666667*IT_0005;
    const ccomplex_t IT_0007 = pow(m_c, 2);
    const ccomplex_t IT_0008 = pow(m_N_2, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0004*IT_0011;
    const ccomplex_t IT_0013 = IT_0001*IT_0002;
    const ccomplex_t IT_0014 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = cpow(IT_0002, -1);
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0020 = IT_0013*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0019;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0015 + IT_0018 + 
      -IT_0020 + -IT_0021);
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = pow(m_sc_R, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0007 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0006*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0023*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0022;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0032 = IT_0013*IT_0031;
    const ccomplex_t IT_0033 = IT_0017*IT_0031;
    const ccomplex_t IT_0034 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0035 = IT_0013*IT_0034;
    const ccomplex_t IT_0036 = IT_0017*IT_0034;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0032 + IT_0033 + 
      -IT_0035 + -IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0040 = 0.666666666666667*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0007 + IT_0024 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0038*IT_0043;
    const ccomplex_t IT_0045 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0046 = IT_0013*IT_0045;
    const ccomplex_t IT_0047 = IT_0017*IT_0045;
    const ccomplex_t IT_0048 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0049 = IT_0013*IT_0048;
    const ccomplex_t IT_0050 = IT_0017*IT_0048;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0046 + IT_0047 + 
      -IT_0049 + -IT_0050);
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0054 = 0.666666666666667*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0007 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0052*IT_0057;
    const ccomplex_t IT_0059 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0060 = IT_0013*IT_0059;
    const ccomplex_t IT_0061 = IT_0017*IT_0059;
    const ccomplex_t IT_0062 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0063 = IT_0013*IT_0062;
    const ccomplex_t IT_0064 = IT_0017*IT_0062;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0060 + IT_0061 + 
      -IT_0063 + -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0068 = 0.666666666666667*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0007 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0066*IT_0071;
    const ccomplex_t IT_0073 = IT_0012 + IT_0028 + -IT_0030 + (-2)*IT_0044 + 2
      *IT_0058 + 2*IT_0072;
    const ccomplex_t IT_0074 = pow(m_Z, -2);
    const ccomplex_t IT_0075 = pow(s_14, 2);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = -IT_0008;
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = s_13*IT_0078;
    const ccomplex_t IT_0080 = pow(m_W, -1);
    const ccomplex_t IT_0081 = sin(beta);
    const ccomplex_t IT_0082 = cpow(IT_0081, -1);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0016*IT_0080*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0025*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0029*IT_0086;
    const ccomplex_t IT_0088 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0089 = IT_0017*IT_0088;
    const ccomplex_t IT_0090 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0091 = IT_0013*IT_0090;
    const ccomplex_t IT_0092 = IT_0017*IT_0090;
    const ccomplex_t IT_0093 = IT_0013*IT_0088;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0089 + -IT_0091 + 
      -IT_0092 + IT_0093);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0016*IT_0080*IT_0082;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0069*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0096*IT_0100;
    const ccomplex_t IT_0102 = IT_0009*IT_0084;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0004*IT_0103;
    const ccomplex_t IT_0105 = IT_0023*IT_0086;
    const ccomplex_t IT_0106 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0107 = IT_0013*IT_0106;
    const ccomplex_t IT_0108 = IT_0017*IT_0106;
    const ccomplex_t IT_0109 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0110 = IT_0013*IT_0109;
    const ccomplex_t IT_0111 = IT_0017*IT_0109;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0107 + IT_0108 + 
      -IT_0110 + -IT_0111);
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0016*IT_0080*IT_0082;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = IT_0041*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = IT_0113*IT_0117;
    const ccomplex_t IT_0119 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0120 = IT_0013*IT_0119;
    const ccomplex_t IT_0121 = IT_0017*IT_0119;
    const ccomplex_t IT_0122 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0123 = IT_0013*IT_0122;
    const ccomplex_t IT_0124 = IT_0017*IT_0122;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0120 + IT_0121 + 
      -IT_0123 + -IT_0124);
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0016*IT_0080*IT_0082;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0055*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = IT_0126*IT_0130;
    const ccomplex_t IT_0132 = -IT_0087 + (-2)*IT_0101 + -IT_0104 + IT_0105 +
       2*IT_0118 + (-2)*IT_0131;
    const ccomplex_t IT_0133 = m_c*m_N_2;
    const ccomplex_t IT_0134 = IT_0078*IT_0133;
    const ccomplex_t IT_0135 = (-0.5)*IT_0028 + 0.5*IT_0030 + IT_0044 + 
      -IT_0058 + -IT_0072;
    const ccomplex_t IT_0136 = s_13*s_14;
    const ccomplex_t IT_0137 = s_34*IT_0008;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = m_c*IT_0078;
    const ccomplex_t IT_0141 = e_em*IT_0013;
    const ccomplex_t IT_0142 = e_em*IT_0017;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0141 + (-3)*IT_0142);
    const ccomplex_t IT_0144 = (-0.166666666666667)*IT_0143;
    const ccomplex_t IT_0145 = IT_0006*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_12 + IT_0007 + -IT_0008 + -IT_0024 
      + -reg_prop, -1);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0148 = IT_0145*IT_0146*IT_0147;
    const ccomplex_t IT_0149 = m_N_2*IT_0023;
    const ccomplex_t IT_0150 = IT_0086*IT_0149;
    const ccomplex_t IT_0151 = m_N_2*IT_0029;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = IT_0149 + IT_0152;
    const ccomplex_t IT_0154 = IT_0086*IT_0153;
    const ccomplex_t IT_0155 = m_N_2*IT_0113;
    const ccomplex_t IT_0156 = IT_0117*IT_0155;
    const ccomplex_t IT_0157 = m_N_3*IT_0052;
    const ccomplex_t IT_0158 = m_N_2*IT_0126;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0157 + IT_0159;
    const ccomplex_t IT_0161 = IT_0130*IT_0160;
    const ccomplex_t IT_0162 = m_N_4*IT_0066;
    const ccomplex_t IT_0163 = m_N_2*IT_0096;
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = IT_0162 + IT_0164;
    const ccomplex_t IT_0166 = IT_0100*IT_0165;
    const ccomplex_t IT_0167 = (-0.666666666666667)*IT_0003;
    const ccomplex_t IT_0168 = IT_0006*IT_0167;
    const ccomplex_t IT_0169 = IT_0146*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = m_c*IT_0170;
    const ccomplex_t IT_0172 = -IT_0148 + (-0.5)*IT_0150 + (-0.5)*IT_0154 + 
      -IT_0156 + -IT_0161 + -IT_0166 + IT_0171;
    const ccomplex_t IT_0173 = IT_0086*IT_0151;
    const ccomplex_t IT_0174 = m_N_1*IT_0038;
    const ccomplex_t IT_0175 = IT_0117*IT_0174;
    const ccomplex_t IT_0176 = 0.5*IT_0173 + IT_0175;
    const ccomplex_t IT_0177 = IT_0172 + IT_0176;
    const ccomplex_t IT_0178 = IT_0084*IT_0167;
    const ccomplex_t IT_0179 = IT_0146*IT_0147*IT_0178;
    const ccomplex_t IT_0180 = IT_0084*IT_0144;
    const ccomplex_t IT_0181 = IT_0146*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = m_c*IT_0182;
    const ccomplex_t IT_0184 = IT_0027*IT_0151;
    const ccomplex_t IT_0185 = IT_0027*IT_0153;
    const ccomplex_t IT_0186 = m_N_2*IT_0038;
    const ccomplex_t IT_0187 = IT_0043*IT_0186;
    const ccomplex_t IT_0188 = m_N_3*IT_0126;
    const ccomplex_t IT_0189 = m_N_2*IT_0052;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = IT_0188 + IT_0190;
    const ccomplex_t IT_0192 = IT_0057*IT_0191;
    const ccomplex_t IT_0193 = m_N_4*IT_0096;
    const ccomplex_t IT_0194 = m_N_2*IT_0066;
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = IT_0193 + IT_0195;
    const ccomplex_t IT_0197 = IT_0071*IT_0196;
    const ccomplex_t IT_0198 = IT_0179 + -IT_0183 + 0.5*IT_0184 + (-0.5)
      *IT_0185 + IT_0187 + IT_0192 + IT_0197;
    const ccomplex_t IT_0199 = IT_0027*IT_0149;
    const ccomplex_t IT_0200 = m_N_1*IT_0113;
    const ccomplex_t IT_0201 = IT_0043*IT_0200;
    const ccomplex_t IT_0202 = (-0.5)*IT_0199 + -IT_0201;
    const ccomplex_t IT_0203 = IT_0198 + IT_0202;
    const ccomplex_t IT_0204 = (-18)*conj(IT_0172);
    const ccomplex_t IT_0205 = -IT_0012 + (-2)*IT_0170;
    const ccomplex_t IT_0206 = s_14*s_34*IT_0074;
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = s_13 + IT_0207;
    const ccomplex_t IT_0209 = s_13*IT_0208;
    const ccomplex_t IT_0210 = IT_0104 + 2*IT_0182;
    const ccomplex_t IT_0211 = IT_0133*IT_0208;
    const ccomplex_t IT_0212 = -IT_0170;
    const ccomplex_t IT_0213 = m_N_2*IT_0208;
    const ccomplex_t IT_0214 = (-6)*conj(IT_0073);
    const ccomplex_t IT_0215 = (-6)*conj(IT_0132);
    const ccomplex_t IT_0216 = pow(m_Z, 2);
    const ccomplex_t IT_0217 = s_14*IT_0216;
    const ccomplex_t IT_0218 = IT_0074*IT_0217;
    const ccomplex_t IT_0219 = (-0.25)*IT_0218;
    const ccomplex_t IT_0220 = s_14 + IT_0219;
    const ccomplex_t IT_0221 = m_c*IT_0220;
    const ccomplex_t IT_0222 = (-4)*IT_0221;
    const ccomplex_t IT_0223 = 0.5*IT_0218;
    const ccomplex_t IT_0224 = s_14 + IT_0223;
    const ccomplex_t IT_0225 = m_c*IT_0224;
    const ccomplex_t IT_0226 = 2*IT_0225;
    const ccomplex_t IT_0227 = 0.5*IT_0087 + IT_0101 + (-0.5)*IT_0105 + 
      -IT_0118 + IT_0131;
    const ccomplex_t IT_0228 = s_34*IT_0216;
    const ccomplex_t IT_0229 = IT_0074*IT_0228;
    const ccomplex_t IT_0230 = (-0.25)*IT_0229;
    const ccomplex_t IT_0231 = s_34 + IT_0230;
    const ccomplex_t IT_0232 = m_N_2*IT_0231;
    const ccomplex_t IT_0233 = (-4)*IT_0232;
    const ccomplex_t IT_0234 = 0.5*IT_0229;
    const ccomplex_t IT_0235 = s_34 + IT_0234;
    const ccomplex_t IT_0236 = m_N_2*IT_0235;
    const ccomplex_t IT_0237 = 2*IT_0236;
    const ccomplex_t IT_0238 = s_13*IT_0216;
    const ccomplex_t IT_0239 = s_14*s_34;
    const ccomplex_t IT_0240 = (-2)*IT_0239;
    const ccomplex_t IT_0241 = IT_0238 + IT_0240;
    const ccomplex_t IT_0242 = IT_0074*IT_0241;
    const ccomplex_t IT_0243 = (-6)*IT_0242;
    const ccomplex_t IT_0244 = 12*s_13;
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = m_c*IT_0208;
    const ccomplex_t IT_0247 = (-6)*IT_0205;
    const ccomplex_t IT_0248 = pow(s_34, 2);
    const ccomplex_t IT_0249 = IT_0074*IT_0248;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = IT_0007 + IT_0250;
    const ccomplex_t IT_0252 = IT_0133*IT_0251;
    const ccomplex_t IT_0253 = s_14*IT_0007;
    const ccomplex_t IT_0254 = s_13*s_34;
    const ccomplex_t IT_0255 = -IT_0254;
    const ccomplex_t IT_0256 = IT_0253 + IT_0255;
    const ccomplex_t IT_0257 = m_N_2*IT_0251;
    const ccomplex_t IT_0258 = (-6)*conj(IT_0205);
    const ccomplex_t IT_0259 = (-6)*conj(IT_0210);
    const ccomplex_t IT_0260 = (-6)*conj(IT_0172);
    const ccomplex_t IT_0261 = (-6)*IT_0073;
    const ccomplex_t IT_0262 = (-6)*IT_0210;
    const ccomplex_t IT_0263 = 6*IT_0212;
    const ccomplex_t IT_0264 = 6*IT_0135;
    const ccomplex_t IT_0265 = 6*IT_0227;
    const ccomplex_t IT_0266 = 6*IT_0182;
    const ccomplex_t IT_0267 = 6*IT_0132*IT_0140 + IT_0133*((-18)*IT_0172 + (
      -18)*IT_0176) + IT_0203*IT_0245 + IT_0247*IT_0257 + IT_0213*IT_0261 +
       IT_0246*IT_0262 + IT_0237*IT_0263 + IT_0233*IT_0264 + IT_0226*IT_0265 +
       IT_0222*IT_0266;
    const ccomplex_t IT_0268 = (-6)*IT_0132;
    const ccomplex_t IT_0269 = 6*IT_0073*IT_0140 + IT_0133*((-18)*IT_0198 + (
      -18)*IT_0202) + IT_0177*IT_0245 + IT_0246*IT_0247 + IT_0257*IT_0262 +
       IT_0222*IT_0263 + IT_0226*IT_0264 + IT_0233*IT_0265 + IT_0237*IT_0266 +
       IT_0213*IT_0268;
    const ccomplex_t IT_0270 = IT_0133*IT_0216;
    const ccomplex_t IT_0271 = (-3)*IT_0270;
    const ccomplex_t IT_0272 = 3*IT_0270;
    const ccomplex_t IT_0273 = pow(m_Z, 4);
    const ccomplex_t IT_0274 = s_13*IT_0273;
    const ccomplex_t IT_0275 = IT_0074*IT_0274;
    const ccomplex_t IT_0276 = (-24)*IT_0239;
    const ccomplex_t IT_0277 = 6*IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = 0.166666666666667*IT_0212;
    const ccomplex_t IT_0279 = 24*IT_0239;
    const ccomplex_t IT_0280 = (-12)*IT_0238;
    const ccomplex_t IT_0281 = 6*IT_0275 + IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = 0.166666666666667*IT_0135;
    const ccomplex_t IT_0283 = 0.166666666666667*IT_0227;
    const ccomplex_t IT_0284 = 0.166666666666667*IT_0182;
    const ccomplex_t IT_0285 = (-6)*IT_0172;
    const ccomplex_t IT_0286 = 0.166666666666667*IT_0211;
    const ccomplex_t IT_0287 = 0.166666666666667*IT_0209;
    const ccomplex_t IT_0288 = s_13*IT_0251;
    const ccomplex_t IT_0289 = 0.166666666666667*IT_0288;
    const ccomplex_t IT_0290 = 0.166666666666667*IT_0139;
    const ccomplex_t IT_0291 = 0.166666666666667*IT_0256;
    const ccomplex_t IT_0292 = 6*conj(IT_0073)*(IT_0073*IT_0079 + IT_0132
      *IT_0134 + IT_0135*IT_0139 + IT_0140*IT_0177) + IT_0133*IT_0203*IT_0204 + 
      (IT_0205*IT_0209 + IT_0210*IT_0211 + IT_0139*IT_0212 + IT_0203*IT_0213)
      *IT_0214 + (IT_0139*IT_0182 + IT_0209*IT_0210 + IT_0176*IT_0213)*IT_0215 +
       6*conj(IT_0172)*(IT_0073*IT_0140 + IT_0212*IT_0222 + IT_0135*IT_0226 +
       IT_0227*IT_0233 + IT_0182*IT_0237 + 0.166666666666667*IT_0177*IT_0245 +
       0.166666666666667*IT_0246*IT_0247) + (IT_0132*IT_0211 + IT_0177*IT_0246 +
       IT_0210*IT_0252 + IT_0212*IT_0256 + IT_0203*IT_0257)*IT_0258 + (IT_0203
      *IT_0246 + IT_0182*IT_0256 + IT_0176*IT_0257)*IT_0259 + (IT_0132*IT_0213 +
       IT_0210*IT_0257)*IT_0260 + (conj(IT_0198) + conj(IT_0202))*IT_0267 + conj
      (IT_0176)*IT_0269 + 6*conj(IT_0135)*(IT_0073*IT_0139 + IT_0177*IT_0226 +
       IT_0203*IT_0233 + IT_0205*IT_0256 + IT_0227*IT_0271 + IT_0182*IT_0272 +
       IT_0277*IT_0278 + IT_0281*IT_0282) + 6*conj(IT_0227)*(IT_0132*IT_0139 +
       IT_0203*IT_0226 + IT_0177*IT_0233 + IT_0210*IT_0256 + IT_0135*IT_0271 +
       IT_0212*IT_0272 + IT_0281*IT_0283 + IT_0277*IT_0284) + 6*conj(IT_0132)*
      (IT_0079*IT_0132 + IT_0073*IT_0134 + IT_0140*IT_0203 + IT_0139*IT_0227 +
       0.166666666666667*IT_0213*IT_0285 + IT_0247*IT_0286) + 6*conj(IT_0205)*
      (IT_0135*IT_0256 + IT_0261*IT_0287 + IT_0247*IT_0289) + 6*conj(IT_0210)*
      (0.166666666666667*IT_0247*IT_0252 + IT_0227*IT_0256 + 0.166666666666667
      *IT_0257*IT_0285 + IT_0261*IT_0286 + IT_0268*IT_0287 + IT_0262*IT_0289) +
       6*conj(IT_0212)*(IT_0177*IT_0222 + IT_0203*IT_0237 + IT_0182*IT_0271 +
       IT_0227*IT_0272 + IT_0278*IT_0281 + IT_0277*IT_0282 + IT_0261*IT_0290 +
       IT_0247*IT_0291) + 6*conj(IT_0182)*(IT_0203*IT_0222 + IT_0177*IT_0237 +
       IT_0212*IT_0271 + IT_0135*IT_0272 + IT_0277*IT_0283 + IT_0281*IT_0284 +
       IT_0268*IT_0290 + IT_0262*IT_0291);
    return create_ccomplex_return(IT_0292);
}

