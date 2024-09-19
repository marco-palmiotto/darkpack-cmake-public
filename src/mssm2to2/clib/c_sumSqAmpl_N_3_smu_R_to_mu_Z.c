#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_smu_R_to_mu_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_smu_R_to_mu_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_mr = param->Gamma_mr;
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
    const ccomplex_t IT_0000 = pow(m_Z, 2);
    const ccomplex_t IT_0001 = s_13*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = pow(m_Z, 4);
    const ccomplex_t IT_0004 = s_13*IT_0003;
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = s_14*s_34;
    const ccomplex_t IT_0008 = (-2)*IT_0007;
    const ccomplex_t IT_0009 = IT_0001 + IT_0006 + IT_0008;
    const ccomplex_t IT_0010 = pow(m_mu, 2);
    const ccomplex_t IT_0011 = pow(m_N_3, 2);
    const ccomplex_t IT_0012 = pow(m_smu_R, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*IT_0015*IT_0018;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0015*IT_0018;
    const ccomplex_t IT_0025 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0018, -1);
    const ccomplex_t IT_0028 = IT_0014*IT_0027;
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0031 = IT_0024*IT_0030;
    const ccomplex_t IT_0032 = IT_0028*IT_0030;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0026 + IT_0029 + 
      -IT_0031 + -IT_0032);
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0010 + IT_0012 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0017*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0034*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0033;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0042 = IT_0024*IT_0041;
    const ccomplex_t IT_0043 = IT_0028*IT_0041;
    const ccomplex_t IT_0044 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0045 = IT_0024*IT_0044;
    const ccomplex_t IT_0046 = IT_0028*IT_0044;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0042 + IT_0043 + 
      -IT_0045 + -IT_0046);
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0015;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0010 + IT_0012 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0048*IT_0053;
    const ccomplex_t IT_0055 = (-0.5)*IT_0038 + 0.5*IT_0040 + -IT_0054;
    const ccomplex_t IT_0056 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0057 = IT_0024*IT_0056;
    const ccomplex_t IT_0058 = IT_0028*IT_0056;
    const ccomplex_t IT_0059 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0060 = IT_0024*IT_0059;
    const ccomplex_t IT_0061 = IT_0028*IT_0059;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0057 + IT_0058 + 
      -IT_0060 + -IT_0061);
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0015;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_23 + IT_0010 + IT_0012 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0063*IT_0068;
    const ccomplex_t IT_0070 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0071 = IT_0024*IT_0070;
    const ccomplex_t IT_0072 = IT_0028*IT_0070;
    const ccomplex_t IT_0073 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0074 = IT_0024*IT_0073;
    const ccomplex_t IT_0075 = IT_0028*IT_0073;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0071 + IT_0072 + 
      -IT_0074 + -IT_0075);
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0015;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_23 + IT_0010 + IT_0012 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = IT_0077*IT_0082;
    const ccomplex_t IT_0084 = -IT_0069 + -IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0055 + IT_0085;
    const ccomplex_t IT_0087 = -conj(IT_0084);
    const ccomplex_t IT_0088 = conj(IT_0055) + IT_0087;
    const ccomplex_t IT_0089 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0090 = IT_0024*IT_0089;
    const ccomplex_t IT_0091 = IT_0028*IT_0089;
    const ccomplex_t IT_0092 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0093 = IT_0024*IT_0092;
    const ccomplex_t IT_0094 = IT_0028*IT_0092;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0090 + IT_0091 + 
      -IT_0093 + -IT_0094);
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = pow(m_W, -1);
    const ccomplex_t IT_0098 = cos(beta);
    const ccomplex_t IT_0099 = cpow(IT_0098, -1);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0027*IT_0097*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0051*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0096*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0027*IT_0097*IT_0099;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0035*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0034*IT_0108;
    const ccomplex_t IT_0110 = IT_0039*IT_0108;
    const ccomplex_t IT_0111 = IT_0104 + (-0.5)*IT_0109 + 0.5*IT_0110;
    const ccomplex_t IT_0112 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0113 = IT_0024*IT_0112;
    const ccomplex_t IT_0114 = IT_0028*IT_0112;
    const ccomplex_t IT_0115 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0116 = IT_0024*IT_0115;
    const ccomplex_t IT_0117 = IT_0028*IT_0115;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0113 + IT_0114 + 
      -IT_0116 + -IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0027*IT_0097*IT_0099;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = IT_0066*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = IT_0119*IT_0123;
    const ccomplex_t IT_0125 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0126 = IT_0024*IT_0125;
    const ccomplex_t IT_0127 = IT_0028*IT_0125;
    const ccomplex_t IT_0128 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0129 = IT_0024*IT_0128;
    const ccomplex_t IT_0130 = IT_0028*IT_0128;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0126 + IT_0127 + 
      -IT_0129 + -IT_0130);
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0027*IT_0097*IT_0099;
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = IT_0080*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = IT_0132*IT_0136;
    const ccomplex_t IT_0138 = IT_0124 + IT_0137;
    const ccomplex_t IT_0139 = conj(IT_0111) + -conj(IT_0138);
    const ccomplex_t IT_0140 = m_mu*m_N_3;
    const ccomplex_t IT_0141 = IT_0000*IT_0140;
    const ccomplex_t IT_0142 = (-3)*IT_0141;
    const ccomplex_t IT_0143 = 2*IT_0055;
    const ccomplex_t IT_0144 = (-0.5)*IT_0143;
    const ccomplex_t IT_0145 = IT_0084 + IT_0144;
    const ccomplex_t IT_0146 = 2*conj(IT_0055);
    const ccomplex_t IT_0147 = (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = conj(IT_0084) + IT_0147;
    const ccomplex_t IT_0149 = pow(s_14, 2);
    const ccomplex_t IT_0150 = IT_0002*IT_0149;
    const ccomplex_t IT_0151 = -IT_0011;
    const ccomplex_t IT_0152 = IT_0150 + IT_0151;
    const ccomplex_t IT_0153 = IT_0140*IT_0152;
    const ccomplex_t IT_0154 = (-2)*IT_0019;
    const ccomplex_t IT_0155 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0156 = IT_0106*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = IT_0154*IT_0157;
    const ccomplex_t IT_0159 = (-2)*IT_0104 + IT_0109 + -IT_0110 + -IT_0158;
    const ccomplex_t IT_0160 = IT_0017*IT_0155;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = IT_0154*IT_0161;
    const ccomplex_t IT_0163 = IT_0038 + -IT_0040 + 2*IT_0054 + IT_0162;
    const ccomplex_t IT_0164 = 2*IT_0163;
    const ccomplex_t IT_0165 = 2*IT_0069 + 2*IT_0083;
    const ccomplex_t IT_0166 = (-2)*IT_0165;
    const ccomplex_t IT_0167 = IT_0164 + IT_0166;
    const ccomplex_t IT_0168 = 2*conj(IT_0163);
    const ccomplex_t IT_0169 = (-2)*conj(IT_0165);
    const ccomplex_t IT_0170 = IT_0168 + IT_0169;
    const ccomplex_t IT_0171 = (-2)*IT_0124 + (-2)*IT_0137;
    const ccomplex_t IT_0172 = 2*IT_0165;
    const ccomplex_t IT_0173 = (-2)*IT_0163;
    const ccomplex_t IT_0174 = IT_0172 + IT_0173;
    const ccomplex_t IT_0175 = 2*conj(IT_0165);
    const ccomplex_t IT_0176 = (-2)*conj(IT_0163);
    const ccomplex_t IT_0177 = IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = m_mu*IT_0152;
    const ccomplex_t IT_0179 = e_em*IT_0024;
    const ccomplex_t IT_0180 = e_em*IT_0028;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*(IT_0179 + -IT_0180);
    const ccomplex_t IT_0182 = 0.5*IT_0181;
    const ccomplex_t IT_0183 = IT_0017*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*m_mu;
    const ccomplex_t IT_0185 = IT_0013*IT_0183*IT_0184;
    const ccomplex_t IT_0186 = m_N_3*IT_0034;
    const ccomplex_t IT_0187 = IT_0108*IT_0186;
    const ccomplex_t IT_0188 = m_N_3*IT_0039;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0186 + IT_0189;
    const ccomplex_t IT_0191 = IT_0108*IT_0190;
    const ccomplex_t IT_0192 = m_N_4*IT_0048;
    const ccomplex_t IT_0193 = m_N_3*IT_0096;
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = IT_0192 + IT_0194;
    const ccomplex_t IT_0196 = IT_0103*IT_0195;
    const ccomplex_t IT_0197 = m_mu*IT_0022;
    const ccomplex_t IT_0198 = -IT_0185 + (-0.5)*IT_0187 + (-0.5)*IT_0191 + 
      -IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = IT_0108*IT_0188;
    const ccomplex_t IT_0200 = m_N_1*IT_0063;
    const ccomplex_t IT_0201 = m_N_3*IT_0119;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = IT_0123*IT_0203;
    const ccomplex_t IT_0205 = m_N_2*IT_0077;
    const ccomplex_t IT_0206 = m_N_3*IT_0132;
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = IT_0136*IT_0208;
    const ccomplex_t IT_0210 = 0.5*IT_0199 + IT_0204 + IT_0209;
    const ccomplex_t IT_0211 = conj(IT_0198) + conj(IT_0210);
    const ccomplex_t IT_0212 = IT_0198 + IT_0210;
    const ccomplex_t IT_0213 = m_N_4*IT_0096;
    const ccomplex_t IT_0214 = m_N_3*IT_0048;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = IT_0213 + IT_0215;
    const ccomplex_t IT_0217 = IT_0053*IT_0216;
    const ccomplex_t IT_0218 = IT_0019*IT_0106;
    const ccomplex_t IT_0219 = IT_0013*IT_0184*IT_0218;
    const ccomplex_t IT_0220 = IT_0106*IT_0182;
    const ccomplex_t IT_0221 = IT_0013*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = m_mu*IT_0222;
    const ccomplex_t IT_0224 = IT_0037*IT_0188;
    const ccomplex_t IT_0225 = IT_0037*IT_0190;
    const ccomplex_t IT_0226 = IT_0217 + IT_0219 + -IT_0223 + 0.5*IT_0224 + (
      -0.5)*IT_0225;
    const ccomplex_t IT_0227 = IT_0037*IT_0186;
    const ccomplex_t IT_0228 = m_N_1*IT_0119;
    const ccomplex_t IT_0229 = m_N_3*IT_0063;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = IT_0068*IT_0231;
    const ccomplex_t IT_0233 = m_N_2*IT_0132;
    const ccomplex_t IT_0234 = m_N_3*IT_0077;
    const ccomplex_t IT_0235 = -IT_0234;
    const ccomplex_t IT_0236 = IT_0233 + IT_0235;
    const ccomplex_t IT_0237 = IT_0082*IT_0236;
    const ccomplex_t IT_0238 = (-0.5)*IT_0227 + -IT_0232 + -IT_0237;
    const ccomplex_t IT_0239 = conj(IT_0226) + conj(IT_0238);
    const ccomplex_t IT_0240 = 2*IT_0159;
    const ccomplex_t IT_0241 = (-2)*IT_0171;
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = IT_0226 + IT_0238;
    const ccomplex_t IT_0244 = 2*conj(IT_0159);
    const ccomplex_t IT_0245 = (-2)*conj(IT_0171);
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = s_14*s_34*IT_0002;
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = s_13 + IT_0248;
    const ccomplex_t IT_0250 = s_13*IT_0249;
    const ccomplex_t IT_0251 = (-2)*IT_0022 + -IT_0162;
    const ccomplex_t IT_0252 = IT_0158 + 2*IT_0222;
    const ccomplex_t IT_0253 = 2*IT_0171;
    const ccomplex_t IT_0254 = (-2)*IT_0159;
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = 2*conj(IT_0171);
    const ccomplex_t IT_0257 = (-2)*conj(IT_0159);
    const ccomplex_t IT_0258 = IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = s_13*s_14;
    const ccomplex_t IT_0260 = s_34*IT_0011;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0259 + IT_0261;
    const ccomplex_t IT_0263 = conj(IT_0023) + conj(IT_0084);
    const ccomplex_t IT_0264 = IT_0023 + IT_0084;
    const ccomplex_t IT_0265 = conj(IT_0138) + conj(IT_0222);
    const ccomplex_t IT_0266 = IT_0138 + IT_0222;
    const ccomplex_t IT_0267 = s_34*IT_0000;
    const ccomplex_t IT_0268 = IT_0002*IT_0267;
    const ccomplex_t IT_0269 = (-0.25)*IT_0268;
    const ccomplex_t IT_0270 = s_34 + IT_0269;
    const ccomplex_t IT_0271 = m_N_3*IT_0270;
    const ccomplex_t IT_0272 = (-4)*IT_0271;
    const ccomplex_t IT_0273 = 2*IT_0198;
    const ccomplex_t IT_0274 = 2*IT_0210;
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = 2*conj(IT_0198);
    const ccomplex_t IT_0277 = 2*conj(IT_0210);
    const ccomplex_t IT_0278 = IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = 2*IT_0226;
    const ccomplex_t IT_0280 = 2*IT_0238;
    const ccomplex_t IT_0281 = IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = 2*conj(IT_0226);
    const ccomplex_t IT_0283 = 2*conj(IT_0238);
    const ccomplex_t IT_0284 = IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = (-2)*IT_0198;
    const ccomplex_t IT_0286 = (-2)*IT_0210;
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = (-2)*conj(IT_0198);
    const ccomplex_t IT_0289 = (-2)*conj(IT_0210);
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = (-2)*IT_0226;
    const ccomplex_t IT_0292 = (-2)*IT_0238;
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = (-2)*conj(IT_0226);
    const ccomplex_t IT_0295 = (-2)*conj(IT_0238);
    const ccomplex_t IT_0296 = IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = s_14*IT_0000;
    const ccomplex_t IT_0298 = IT_0002*IT_0297;
    const ccomplex_t IT_0299 = 0.5*IT_0298;
    const ccomplex_t IT_0300 = s_14 + IT_0299;
    const ccomplex_t IT_0301 = m_mu*IT_0300;
    const ccomplex_t IT_0302 = 2*IT_0301;
    const ccomplex_t IT_0303 = 3*IT_0141;
    const ccomplex_t IT_0304 = 2*IT_0023;
    const ccomplex_t IT_0305 = 2*conj(IT_0023);
    const ccomplex_t IT_0306 = IT_0001 + IT_0008;
    const ccomplex_t IT_0307 = IT_0002*IT_0306;
    const ccomplex_t IT_0308 = (-0.5)*IT_0307;
    const ccomplex_t IT_0309 = s_13 + IT_0308;
    const ccomplex_t IT_0310 = 4*IT_0226;
    const ccomplex_t IT_0311 = (-6)*conj(IT_0198);
    const ccomplex_t IT_0312 = s_13*IT_0152;
    const ccomplex_t IT_0313 = IT_0140*IT_0249;
    const ccomplex_t IT_0314 = 2*IT_0251;
    const ccomplex_t IT_0315 = 2*conj(IT_0251);
    const ccomplex_t IT_0316 = (-2)*IT_0251;
    const ccomplex_t IT_0317 = (-2)*conj(IT_0251);
    const ccomplex_t IT_0318 = pow(s_34, 2);
    const ccomplex_t IT_0319 = IT_0002*IT_0318;
    const ccomplex_t IT_0320 = -IT_0319;
    const ccomplex_t IT_0321 = IT_0010 + IT_0320;
    const ccomplex_t IT_0322 = IT_0140*IT_0321;
    const ccomplex_t IT_0323 = s_13*IT_0321;
    const ccomplex_t IT_0324 = (-2)*IT_0252;
    const ccomplex_t IT_0325 = m_mu*IT_0249;
    const ccomplex_t IT_0326 = (-2)*conj(IT_0252);
    const ccomplex_t IT_0327 = s_14*IT_0010;
    const ccomplex_t IT_0328 = s_13*s_34;
    const ccomplex_t IT_0329 = -IT_0328;
    const ccomplex_t IT_0330 = IT_0327 + IT_0329;
    const ccomplex_t IT_0331 = conj(IT_0111)*IT_0252 + IT_0111*conj(IT_0252) +
       0.5*conj(IT_0055)*IT_0314 + 0.5*IT_0055*IT_0315 + 0.5*IT_0263*IT_0316 +
       0.5*IT_0264*IT_0317 + 0.5*IT_0265*IT_0324 + 0.5*IT_0266*IT_0326;
    const ccomplex_t IT_0332 = m_N_3*IT_0249;
    const ccomplex_t IT_0333 = m_N_3*IT_0321;
    const ccomplex_t IT_0334 = (-0.25)*IT_0005;
    const ccomplex_t IT_0335 = IT_0007 + IT_0334;
    const ccomplex_t IT_0336 = (-0.25)*IT_0298;
    const ccomplex_t IT_0337 = s_14 + IT_0336;
    const ccomplex_t IT_0338 = m_mu*IT_0337;
    const ccomplex_t IT_0339 = (-4)*IT_0338;
    const ccomplex_t IT_0340 = (-0.5)*IT_0339;
    const ccomplex_t IT_0341 = 0.5*IT_0268;
    const ccomplex_t IT_0342 = s_34 + IT_0341;
    const ccomplex_t IT_0343 = m_N_3*IT_0342;
    const ccomplex_t IT_0344 = 2*IT_0343;
    const ccomplex_t IT_0345 = (-0.5)*IT_0344;
    const ccomplex_t IT_0346 = (-0.5)*IT_0142;
    const ccomplex_t IT_0347 = (-4)*IT_0009*(IT_0023*conj(IT_0023) + IT_0086
      *IT_0088 + (IT_0111 + -IT_0138)*IT_0139) + (-2)*IT_0142*(IT_0088*IT_0138 +
       IT_0086*conj(IT_0138) + conj(IT_0111)*IT_0145 + IT_0111*IT_0148) +
       IT_0153*(conj(IT_0159)*IT_0167 + IT_0159*IT_0170 + conj(IT_0171)*IT_0174 
      + IT_0171*IT_0177) + IT_0178*(IT_0167*IT_0211 + IT_0170*IT_0212 + IT_0239
      *IT_0242 + IT_0243*IT_0246) + IT_0250*(IT_0177*IT_0251 + IT_0174*conj
      (IT_0251) + conj(IT_0252)*IT_0255 + IT_0252*IT_0258) + IT_0262*(conj
      (IT_0055)*IT_0167 + IT_0055*IT_0170 + conj(IT_0111)*IT_0242 + IT_0111
      *IT_0246 + IT_0174*IT_0263 + IT_0177*IT_0264 + IT_0255*IT_0265 + IT_0258
      *IT_0266) + IT_0272*(conj(IT_0111)*IT_0275 + IT_0111*IT_0278 + conj
      (IT_0055)*IT_0281 + IT_0055*IT_0284 + conj(IT_0138)*IT_0287 + IT_0138
      *IT_0290 + conj(IT_0084)*IT_0293 + IT_0084*IT_0296) + (conj(IT_0055)
      *IT_0275 + IT_0055*IT_0278 + conj(IT_0111)*IT_0281 + IT_0111*IT_0284 +
       conj(IT_0084)*IT_0287 + IT_0084*IT_0290 + conj(IT_0138)*IT_0293 + IT_0138
      *IT_0296)*IT_0302 + IT_0303*(conj(IT_0111)*IT_0304 + IT_0111*IT_0305) + 4
      *IT_0309*(IT_0211*IT_0212 + IT_0239*(IT_0238 + 0.25*IT_0310)) + (-6)
      *IT_0140*(IT_0212*IT_0239 + IT_0243*(conj(IT_0210) + (-0.166666666666667)
      *IT_0311)) + (conj(IT_0163)*IT_0164 + conj(IT_0165)*IT_0174 + IT_0165
      *IT_0176 + conj(IT_0159)*IT_0240 + conj(IT_0171)*IT_0255 + IT_0171*IT_0257
      )*IT_0312 + IT_0313*(IT_0177*IT_0252 + IT_0174*conj(IT_0252) + conj
      (IT_0171)*IT_0314 + IT_0171*IT_0315 + conj(IT_0159)*IT_0316 + IT_0159
      *IT_0317) + (conj(IT_0252)*IT_0316 + IT_0252*IT_0317)*IT_0322 + IT_0323*
      (conj(IT_0251)*IT_0316 + conj(IT_0252)*IT_0324) + IT_0325*(IT_0211*IT_0316
       + IT_0212*IT_0317 + IT_0239*IT_0324 + IT_0243*IT_0326) + 2*IT_0330
      *IT_0331 + (IT_0174*IT_0239 + IT_0177*IT_0243 + conj(IT_0210)*IT_0255 +
       IT_0210*IT_0258 + conj(IT_0171)*IT_0273 + IT_0171*IT_0276 + conj(IT_0159)
      *IT_0285 + IT_0159*IT_0288)*IT_0332 + (conj(IT_0252)*IT_0285 + IT_0252
      *IT_0288 + IT_0239*IT_0316 + IT_0243*IT_0317 + conj(IT_0210)*IT_0324 +
       IT_0210*IT_0326)*IT_0333 + (-2)*conj(IT_0023)*(IT_0138*IT_0303 + 4
      *IT_0086*IT_0335 + IT_0275*IT_0340 + IT_0281*IT_0345) + (-2)*IT_0023*(conj
      (IT_0138)*IT_0303 + 4*IT_0088*IT_0335 + IT_0278*IT_0340 + IT_0284*IT_0345)
       + (-8)*conj(IT_0222)*(0.25*IT_0145*IT_0303 + (IT_0111 + -IT_0138)*IT_0335
       + 0.25*IT_0281*IT_0340 + 0.25*IT_0275*IT_0345 + 0.25*IT_0304*IT_0346) + (
      -4)*IT_0222*(IT_0009*conj(IT_0222) + 0.5*IT_0148*IT_0303 + 2*IT_0139
      *IT_0335 + 0.5*IT_0284*IT_0340 + 0.5*IT_0278*IT_0345 + 0.5*IT_0305*IT_0346);
    return create_ccomplex_return(IT_0347);
}

