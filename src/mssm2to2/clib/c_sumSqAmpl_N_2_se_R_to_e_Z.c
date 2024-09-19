#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_se_R_to_e_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_se_R_to_e_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
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
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_er = param->Gamma_er;
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
    const ccomplex_t IT_0000 = s_14*s_34;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(m_Z, 4);
    const ccomplex_t IT_0003 = s_13*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = (-0.25)*IT_0004;
    const ccomplex_t IT_0006 = IT_0000 + IT_0005;
    const ccomplex_t IT_0007 = pow(m_e, 2);
    const ccomplex_t IT_0008 = pow(m_N_2, 2);
    const ccomplex_t IT_0009 = pow(m_se_R, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + IT_0007 + -IT_0008 + -IT_0009 
      + -reg_prop, -1);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*e_em*IT_0012*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = IT_0012*IT_0015;
    const ccomplex_t IT_0022 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0015, -1);
    const ccomplex_t IT_0025 = IT_0011*IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = IT_0025*IT_0027;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0023 + IT_0026 + 
      -IT_0028 + -IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0012;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_23 + IT_0007 + IT_0009 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0031*IT_0036;
    const ccomplex_t IT_0038 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0039 = IT_0021*IT_0038;
    const ccomplex_t IT_0040 = IT_0025*IT_0038;
    const ccomplex_t IT_0041 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0042 = IT_0021*IT_0041;
    const ccomplex_t IT_0043 = IT_0025*IT_0041;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0039 + IT_0040 + 
      -IT_0042 + -IT_0043);
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_23 + IT_0007 + IT_0009 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0045*IT_0048;
    const ccomplex_t IT_0050 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0051 = IT_0021*IT_0050;
    const ccomplex_t IT_0052 = IT_0025*IT_0050;
    const ccomplex_t IT_0053 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0054 = IT_0021*IT_0053;
    const ccomplex_t IT_0055 = IT_0025*IT_0053;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0051 + IT_0052 + 
      -IT_0054 + -IT_0055);
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0012;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_23 + IT_0007 + IT_0009 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0057*IT_0062;
    const ccomplex_t IT_0064 = -IT_0037 + (-0.5)*IT_0049 + -IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0044;
    const ccomplex_t IT_0066 = IT_0048*IT_0065;
    const ccomplex_t IT_0067 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0068 = IT_0021*IT_0067;
    const ccomplex_t IT_0069 = IT_0025*IT_0067;
    const ccomplex_t IT_0070 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0071 = IT_0021*IT_0070;
    const ccomplex_t IT_0072 = IT_0025*IT_0070;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + IT_0069 + 
      -IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0012;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_23 + IT_0007 + IT_0009 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0074*IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0066 + -IT_0080;
    const ccomplex_t IT_0082 = pow(m_W, -1);
    const ccomplex_t IT_0083 = cos(beta);
    const ccomplex_t IT_0084 = cpow(IT_0083, -1);
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d2)
      *e_em*IT_0024*IT_0082*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = e_em*IT_0021;
    const ccomplex_t IT_0088 = e_em*IT_0025;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0087 + -IT_0088);
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = IT_0086*IT_0090;
    const ccomplex_t IT_0092 = IT_0010*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = IT_0046*IT_0086;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0045*IT_0095;
    const ccomplex_t IT_0097 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0098 = IT_0021*IT_0097;
    const ccomplex_t IT_0099 = IT_0025*IT_0097;
    const ccomplex_t IT_0100 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0101 = IT_0021*IT_0100;
    const ccomplex_t IT_0102 = IT_0025*IT_0100;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0098 + IT_0099 + 
      -IT_0101 + -IT_0102);
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d1)
      *e_em*IT_0024*IT_0082*IT_0084;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0077*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0104*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0096 + IT_0109;
    const ccomplex_t IT_0111 = IT_0065*IT_0095;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0114 = IT_0021*IT_0113;
    const ccomplex_t IT_0115 = IT_0025*IT_0113;
    const ccomplex_t IT_0116 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0117 = IT_0021*IT_0116;
    const ccomplex_t IT_0118 = IT_0025*IT_0116;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0114 + IT_0115 + 
      -IT_0117 + -IT_0118);
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0024*IT_0082*IT_0084;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0060*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = IT_0120*IT_0124;
    const ccomplex_t IT_0126 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0127 = IT_0025*IT_0126;
    const ccomplex_t IT_0128 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0129 = IT_0021*IT_0128;
    const ccomplex_t IT_0130 = IT_0025*IT_0128;
    const ccomplex_t IT_0131 = IT_0021*IT_0126;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0127 + -IT_0129 + 
      -IT_0130 + IT_0131);
    const ccomplex_t IT_0133 = -IT_0132;
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0024*IT_0082*IT_0084;
    const ccomplex_t IT_0136 = 0.5*IT_0135;
    const ccomplex_t IT_0137 = IT_0034*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = IT_0134*IT_0138;
    const ccomplex_t IT_0140 = IT_0125 + IT_0139;
    const ccomplex_t IT_0141 = IT_0112 + IT_0140;
    const ccomplex_t IT_0142 = -IT_0110 + IT_0141;
    const ccomplex_t IT_0143 = conj(IT_0112) + conj(IT_0140);
    const ccomplex_t IT_0144 = -conj(IT_0110) + IT_0143;
    const ccomplex_t IT_0145 = conj(IT_0064) + -conj(IT_0081);
    const ccomplex_t IT_0146 = pow(m_Z, 2);
    const ccomplex_t IT_0147 = s_13*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0004;
    const ccomplex_t IT_0149 = (-2)*IT_0000;
    const ccomplex_t IT_0150 = IT_0147 + IT_0148 + IT_0149;
    const ccomplex_t IT_0151 = m_e*m_N_2;
    const ccomplex_t IT_0152 = pow(s_14, 2);
    const ccomplex_t IT_0153 = IT_0001*IT_0152;
    const ccomplex_t IT_0154 = -IT_0008;
    const ccomplex_t IT_0155 = IT_0153 + IT_0154;
    const ccomplex_t IT_0156 = IT_0151*IT_0155;
    const ccomplex_t IT_0157 = (-2)*IT_0016;
    const ccomplex_t IT_0158 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0159 = IT_0086*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0157*IT_0160;
    const ccomplex_t IT_0162 = -IT_0111 + (-2)*IT_0125 + (-2)*IT_0139 + 
      -IT_0161;
    const ccomplex_t IT_0163 = IT_0014*IT_0158;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0157*IT_0164;
    const ccomplex_t IT_0166 = 2*IT_0037 + IT_0049 + 2*IT_0063 + IT_0165;
    const ccomplex_t IT_0167 = 2*IT_0166;
    const ccomplex_t IT_0168 = IT_0066 + 2*IT_0080;
    const ccomplex_t IT_0169 = (-2)*IT_0168;
    const ccomplex_t IT_0170 = IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = 2*conj(IT_0166);
    const ccomplex_t IT_0172 = (-2)*conj(IT_0168);
    const ccomplex_t IT_0173 = IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = -IT_0096 + (-2)*IT_0109;
    const ccomplex_t IT_0175 = 2*IT_0168;
    const ccomplex_t IT_0176 = (-2)*IT_0166;
    const ccomplex_t IT_0177 = IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = 2*conj(IT_0168);
    const ccomplex_t IT_0179 = (-2)*conj(IT_0166);
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = m_e*IT_0155;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_e;
    const ccomplex_t IT_0183 = IT_0014*IT_0090;
    const ccomplex_t IT_0184 = IT_0010*IT_0182*IT_0183;
    const ccomplex_t IT_0185 = m_N_2*IT_0045;
    const ccomplex_t IT_0186 = m_N_2*IT_0065;
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = IT_0185 + IT_0187;
    const ccomplex_t IT_0189 = IT_0095*IT_0188;
    const ccomplex_t IT_0190 = m_e*IT_0019;
    const ccomplex_t IT_0191 = m_N_2*IT_0120;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = m_N_3*IT_0057;
    const ccomplex_t IT_0194 = IT_0192 + IT_0193;
    const ccomplex_t IT_0195 = IT_0124*IT_0194;
    const ccomplex_t IT_0196 = m_N_2*IT_0134;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = m_N_4*IT_0031;
    const ccomplex_t IT_0199 = IT_0197 + IT_0198;
    const ccomplex_t IT_0200 = IT_0138*IT_0199;
    const ccomplex_t IT_0201 = -IT_0184 + (-0.5)*IT_0189 + IT_0190 + -IT_0195 
      + -IT_0200;
    const ccomplex_t IT_0202 = m_N_1*IT_0074;
    const ccomplex_t IT_0203 = m_N_2*IT_0104;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = IT_0108*IT_0205;
    const ccomplex_t IT_0207 = (-0.5)*IT_0189 + IT_0206;
    const ccomplex_t IT_0208 = conj(IT_0201) + conj(IT_0207);
    const ccomplex_t IT_0209 = IT_0201 + IT_0207;
    const ccomplex_t IT_0210 = IT_0016*IT_0086;
    const ccomplex_t IT_0211 = IT_0010*IT_0182*IT_0210;
    const ccomplex_t IT_0212 = IT_0048*IT_0188;
    const ccomplex_t IT_0213 = m_e*IT_0093;
    const ccomplex_t IT_0214 = m_N_3*IT_0120;
    const ccomplex_t IT_0215 = m_N_2*IT_0057;
    const ccomplex_t IT_0216 = -IT_0215;
    const ccomplex_t IT_0217 = IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = IT_0062*IT_0217;
    const ccomplex_t IT_0219 = m_N_4*IT_0134;
    const ccomplex_t IT_0220 = m_N_2*IT_0031;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = IT_0219 + IT_0221;
    const ccomplex_t IT_0223 = IT_0036*IT_0222;
    const ccomplex_t IT_0224 = IT_0211 + (-0.5)*IT_0212 + -IT_0213 + IT_0218 +
       IT_0223;
    const ccomplex_t IT_0225 = m_N_1*IT_0104;
    const ccomplex_t IT_0226 = m_N_2*IT_0074;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = IT_0225 + IT_0227;
    const ccomplex_t IT_0229 = IT_0079*IT_0228;
    const ccomplex_t IT_0230 = (-0.5)*IT_0212 + -IT_0229;
    const ccomplex_t IT_0231 = conj(IT_0224) + conj(IT_0230);
    const ccomplex_t IT_0232 = 2*IT_0162;
    const ccomplex_t IT_0233 = (-2)*IT_0174;
    const ccomplex_t IT_0234 = IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = IT_0224 + IT_0230;
    const ccomplex_t IT_0236 = 2*conj(IT_0162);
    const ccomplex_t IT_0237 = (-2)*conj(IT_0174);
    const ccomplex_t IT_0238 = IT_0236 + IT_0237;
    const ccomplex_t IT_0239 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0240 = -IT_0239;
    const ccomplex_t IT_0241 = s_13 + IT_0240;
    const ccomplex_t IT_0242 = s_13*IT_0241;
    const ccomplex_t IT_0243 = (-2)*IT_0019 + -IT_0165;
    const ccomplex_t IT_0244 = 2*IT_0093 + IT_0161;
    const ccomplex_t IT_0245 = 2*IT_0174;
    const ccomplex_t IT_0246 = (-2)*IT_0162;
    const ccomplex_t IT_0247 = IT_0245 + IT_0246;
    const ccomplex_t IT_0248 = 2*conj(IT_0174);
    const ccomplex_t IT_0249 = (-2)*conj(IT_0162);
    const ccomplex_t IT_0250 = IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = s_13*s_14;
    const ccomplex_t IT_0252 = s_34*IT_0008;
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = IT_0251 + IT_0253;
    const ccomplex_t IT_0255 = conj(IT_0020) + conj(IT_0081);
    const ccomplex_t IT_0256 = IT_0020 + IT_0081;
    const ccomplex_t IT_0257 = conj(IT_0093) + conj(IT_0110);
    const ccomplex_t IT_0258 = IT_0093 + IT_0110;
    const ccomplex_t IT_0259 = s_34*IT_0146;
    const ccomplex_t IT_0260 = IT_0001*IT_0259;
    const ccomplex_t IT_0261 = 0.5*IT_0260;
    const ccomplex_t IT_0262 = s_34 + IT_0261;
    const ccomplex_t IT_0263 = m_N_2*IT_0262;
    const ccomplex_t IT_0264 = 2*IT_0263;
    const ccomplex_t IT_0265 = 2*IT_0201;
    const ccomplex_t IT_0266 = 2*IT_0207;
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = 2*conj(IT_0201);
    const ccomplex_t IT_0269 = 2*conj(IT_0207);
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = 2*IT_0224;
    const ccomplex_t IT_0272 = 2*IT_0230;
    const ccomplex_t IT_0273 = IT_0271 + IT_0272;
    const ccomplex_t IT_0274 = 2*conj(IT_0224);
    const ccomplex_t IT_0275 = 2*conj(IT_0230);
    const ccomplex_t IT_0276 = IT_0274 + IT_0275;
    const ccomplex_t IT_0277 = s_14*IT_0146;
    const ccomplex_t IT_0278 = IT_0001*IT_0277;
    const ccomplex_t IT_0279 = (-0.25)*IT_0278;
    const ccomplex_t IT_0280 = s_14 + IT_0279;
    const ccomplex_t IT_0281 = m_e*IT_0280;
    const ccomplex_t IT_0282 = (-4)*IT_0281;
    const ccomplex_t IT_0283 = (-0.25)*IT_0260;
    const ccomplex_t IT_0284 = s_34 + IT_0283;
    const ccomplex_t IT_0285 = m_N_2*IT_0284;
    const ccomplex_t IT_0286 = (-4)*IT_0285;
    const ccomplex_t IT_0287 = (-2)*IT_0201;
    const ccomplex_t IT_0288 = (-2)*IT_0207;
    const ccomplex_t IT_0289 = IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = (-2)*conj(IT_0201);
    const ccomplex_t IT_0291 = (-2)*conj(IT_0207);
    const ccomplex_t IT_0292 = IT_0290 + IT_0291;
    const ccomplex_t IT_0293 = (-2)*IT_0224;
    const ccomplex_t IT_0294 = (-2)*IT_0230;
    const ccomplex_t IT_0295 = IT_0293 + IT_0294;
    const ccomplex_t IT_0296 = (-2)*conj(IT_0224);
    const ccomplex_t IT_0297 = (-2)*conj(IT_0230);
    const ccomplex_t IT_0298 = IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = 0.5*IT_0278;
    const ccomplex_t IT_0300 = s_14 + IT_0299;
    const ccomplex_t IT_0301 = m_e*IT_0300;
    const ccomplex_t IT_0302 = 2*IT_0301;
    const ccomplex_t IT_0303 = IT_0147 + IT_0149;
    const ccomplex_t IT_0304 = IT_0001*IT_0303;
    const ccomplex_t IT_0305 = (-0.5)*IT_0304;
    const ccomplex_t IT_0306 = s_13 + IT_0305;
    const ccomplex_t IT_0307 = 4*IT_0224;
    const ccomplex_t IT_0308 = (-6)*conj(IT_0201);
    const ccomplex_t IT_0309 = s_13*IT_0155;
    const ccomplex_t IT_0310 = IT_0151*IT_0241;
    const ccomplex_t IT_0311 = 2*IT_0243;
    const ccomplex_t IT_0312 = 2*conj(IT_0243);
    const ccomplex_t IT_0313 = (-2)*IT_0243;
    const ccomplex_t IT_0314 = (-2)*conj(IT_0243);
    const ccomplex_t IT_0315 = pow(s_34, 2);
    const ccomplex_t IT_0316 = IT_0001*IT_0315;
    const ccomplex_t IT_0317 = -IT_0316;
    const ccomplex_t IT_0318 = IT_0007 + IT_0317;
    const ccomplex_t IT_0319 = IT_0151*IT_0318;
    const ccomplex_t IT_0320 = s_13*IT_0318;
    const ccomplex_t IT_0321 = (-2)*IT_0244;
    const ccomplex_t IT_0322 = m_e*IT_0241;
    const ccomplex_t IT_0323 = (-2)*conj(IT_0244);
    const ccomplex_t IT_0324 = s_14*IT_0007;
    const ccomplex_t IT_0325 = s_13*s_34;
    const ccomplex_t IT_0326 = -IT_0325;
    const ccomplex_t IT_0327 = IT_0324 + IT_0326;
    const ccomplex_t IT_0328 = 2*IT_0244;
    const ccomplex_t IT_0329 = 2*conj(IT_0244);
    const ccomplex_t IT_0330 = m_N_2*IT_0241;
    const ccomplex_t IT_0331 = m_N_2*IT_0318;
    const ccomplex_t IT_0332 = IT_0146*IT_0151;
    const ccomplex_t IT_0333 = 3*IT_0332;
    const ccomplex_t IT_0334 = 2*conj(IT_0020);
    const ccomplex_t IT_0335 = 2*IT_0064;
    const ccomplex_t IT_0336 = 2*conj(IT_0064);
    const ccomplex_t IT_0337 = 2*IT_0020;
    const ccomplex_t IT_0338 = (-0.5)*IT_0143;
    const ccomplex_t IT_0339 = (-3)*IT_0332;
    const ccomplex_t IT_0340 = (-2)*(conj(IT_0064)*IT_0110 + IT_0064*conj
      (IT_0110) + conj(IT_0081)*IT_0142 + IT_0081*IT_0144 + (-0.5)*IT_0093
      *IT_0334 + (-0.5)*IT_0141*IT_0336 + (-0.5)*conj(IT_0093)*IT_0337 + IT_0335
      *IT_0338)*IT_0339;
    const ccomplex_t IT_0341 = (-8)*IT_0006*(conj(IT_0020)*(IT_0064 + -IT_0081
      ) + conj(IT_0093)*IT_0142 + IT_0093*IT_0144 + IT_0020*IT_0145) + (-4)*
      (IT_0020*conj(IT_0020) + IT_0093*conj(IT_0093) + IT_0142*IT_0144 + 
      (IT_0064 + -IT_0081)*IT_0145)*IT_0150 + IT_0156*(conj(IT_0162)*IT_0170 +
       IT_0162*IT_0173 + conj(IT_0174)*IT_0177 + IT_0174*IT_0180) + IT_0181*
      (IT_0170*IT_0208 + IT_0173*IT_0209 + IT_0231*IT_0234 + IT_0235*IT_0238) +
       IT_0242*(IT_0180*IT_0243 + IT_0177*conj(IT_0243) + conj(IT_0244)*IT_0247 
      + IT_0244*IT_0250) + IT_0254*(conj(IT_0064)*IT_0170 + IT_0064*IT_0173 +
       IT_0143*IT_0234 + IT_0141*IT_0238 + IT_0177*IT_0255 + IT_0180*IT_0256 +
       IT_0247*IT_0257 + IT_0250*IT_0258) + IT_0264*(conj(IT_0093)*IT_0267 +
       IT_0093*IT_0270 + conj(IT_0020)*IT_0273 + IT_0020*IT_0276) + (conj
      (IT_0020)*IT_0267 + IT_0020*IT_0270 + conj(IT_0093)*IT_0273 + IT_0093
      *IT_0276)*IT_0282 + IT_0286*(IT_0143*IT_0267 + IT_0141*IT_0270 + conj
      (IT_0064)*IT_0273 + IT_0064*IT_0276 + conj(IT_0110)*IT_0289 + IT_0110
      *IT_0292 + conj(IT_0081)*IT_0295 + IT_0081*IT_0298) + (conj(IT_0064)
      *IT_0267 + IT_0064*IT_0270 + IT_0143*IT_0273 + IT_0141*IT_0276 + conj
      (IT_0081)*IT_0289 + IT_0081*IT_0292 + conj(IT_0110)*IT_0295 + IT_0110
      *IT_0298)*IT_0302 + 4*IT_0306*(IT_0208*IT_0209 + IT_0231*(IT_0230 + 0.25
      *IT_0307)) + (-6)*IT_0151*(IT_0209*IT_0231 + IT_0235*(conj(IT_0207) + (
      -0.166666666666667)*IT_0308)) + (conj(IT_0166)*IT_0167 + conj(IT_0168)
      *IT_0177 + IT_0168*IT_0179 + conj(IT_0162)*IT_0232 + conj(IT_0174)*IT_0247
       + IT_0174*IT_0249)*IT_0309 + IT_0310*(IT_0180*IT_0244 + IT_0177*conj
      (IT_0244) + conj(IT_0174)*IT_0311 + IT_0174*IT_0312 + conj(IT_0162)
      *IT_0313 + IT_0162*IT_0314) + (conj(IT_0244)*IT_0313 + IT_0244*IT_0314)
      *IT_0319 + IT_0320*(conj(IT_0243)*IT_0313 + conj(IT_0244)*IT_0321) +
       IT_0322*(IT_0208*IT_0313 + IT_0209*IT_0314 + IT_0231*IT_0321 + IT_0235
      *IT_0323) + IT_0327*(conj(IT_0064)*IT_0311 + IT_0064*IT_0312 + IT_0255
      *IT_0313 + IT_0256*IT_0314 + IT_0257*IT_0321 + IT_0258*IT_0323 + IT_0143
      *IT_0328 + IT_0141*IT_0329) + (IT_0177*IT_0231 + IT_0180*IT_0235 + conj
      (IT_0207)*IT_0247 + IT_0207*IT_0250 + conj(IT_0174)*IT_0265 + IT_0174
      *IT_0268 + conj(IT_0162)*IT_0287 + IT_0162*IT_0290)*IT_0330 + (conj
      (IT_0244)*IT_0287 + IT_0244*IT_0290 + IT_0231*IT_0313 + IT_0235*IT_0314 +
       conj(IT_0207)*IT_0321 + IT_0207*IT_0323)*IT_0331 + (-2)*IT_0333*(conj
      (IT_0020)*IT_0110 + IT_0020*conj(IT_0110) + (-0.5)*IT_0141*IT_0334 + conj
      (IT_0093)*(IT_0081 + (-0.5)*IT_0335) + IT_0093*(conj(IT_0081) + (-0.5)
      *IT_0336) + IT_0337*IT_0338) + IT_0340;
    return create_ccomplex_return(IT_0341);
}

