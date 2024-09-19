#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_ss_L_to_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_ss_L_to_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
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
    const creal_t m_ss_L = param->m_ss_L;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_ss_L, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-18)*IT_0001;
    const ccomplex_t IT_0003 = s_23*s_24;
    const ccomplex_t IT_0004 = 36*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_s, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_B2)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = conj(N_W2)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0012 + (-3)*IT_0017);
    const ccomplex_t IT_0019 = 0.166666666666667*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = cos(beta);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0014*IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = IT_0007*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0030 = conj(N_B4)*e_em;
    const ccomplex_t IT_0031 = IT_0009*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = conj(N_W4)*e_em;
    const ccomplex_t IT_0034 = IT_0014*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0032 + (-3)*IT_0035);
    const ccomplex_t IT_0037 = 0.166666666666667*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0014*IT_0020*IT_0022;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = IT_0029*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0014*IT_0020*IT_0022;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = conj(N_B3)*e_em;
    const ccomplex_t IT_0048 = IT_0009*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = conj(N_W3)*e_em;
    const ccomplex_t IT_0051 = IT_0014*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-3)*IT_0052);
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = IT_0046*IT_0054;
    const ccomplex_t IT_0056 = IT_0044*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = conj(N_B1)*e_em;
    const ccomplex_t IT_0059 = IT_0009*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = conj(N_W1)*e_em;
    const ccomplex_t IT_0062 = IT_0014*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0060 + (-3)*IT_0063);
    const ccomplex_t IT_0065 = 0.166666666666667*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0014*IT_0020*IT_0022;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0057 + -IT_0071;
    const ccomplex_t IT_0073 = IT_0028 + IT_0043 + IT_0072;
    const ccomplex_t IT_0074 = conj(IT_0028) + conj(IT_0043) + conj(IT_0072);
    const ccomplex_t IT_0075 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0040*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0055*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0025*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0080 + IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0068*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = IT_0077 + IT_0088;
    const ccomplex_t IT_0090 = conj(IT_0084) + conj(IT_0087);
    const ccomplex_t IT_0091 = conj(IT_0077) + IT_0090;
    const ccomplex_t IT_0092 = -IT_0027 + -IT_0042;
    const ccomplex_t IT_0093 = -IT_0057;
    const ccomplex_t IT_0094 = -IT_0071;
    const ccomplex_t IT_0095 = IT_0092 + IT_0093 + IT_0094;
    const ccomplex_t IT_0096 = conj(IT_0092) + conj(IT_0093) + conj(IT_0094);
    const ccomplex_t IT_0097 = IT_0077 + IT_0080;
    const ccomplex_t IT_0098 = IT_0083 + IT_0087;
    const ccomplex_t IT_0099 = conj(IT_0097) + conj(IT_0098);
    const ccomplex_t IT_0100 = (-6)*IT_0001;
    const ccomplex_t IT_0101 = 12*IT_0003;
    const ccomplex_t IT_0102 = IT_0100 + IT_0101;
    const ccomplex_t IT_0103 = IT_0097 + IT_0098;
    const ccomplex_t IT_0104 = m_s*IT_0077;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = m_s*IT_0071;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = 6*s_34;
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = m_s*IT_0027;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = IT_0108*IT_0111;
    const ccomplex_t IT_0113 = IT_0109 + IT_0112;
    const ccomplex_t IT_0114 = m_s*IT_0057;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = IT_0108*IT_0115;
    const ccomplex_t IT_0117 = m_s*IT_0042;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = IT_0108*IT_0118;
    const ccomplex_t IT_0120 = IT_0113 + IT_0116 + IT_0119;
    const ccomplex_t IT_0121 = m_s*IT_0083;
    const ccomplex_t IT_0122 = -IT_0121;
    const ccomplex_t IT_0123 = 18*s_34;
    const ccomplex_t IT_0124 = IT_0122*IT_0123;
    const ccomplex_t IT_0125 = cpow(IT_0046, 2);
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*m_N_3*IT_0125;
    const ccomplex_t IT_0127 = IT_0078*IT_0126;
    const ccomplex_t IT_0128 = cpow(IT_0039, 2);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*m_N_4*IT_0128;
    const ccomplex_t IT_0130 = IT_0075*IT_0129;
    const ccomplex_t IT_0131 = IT_0127 + IT_0130;
    const ccomplex_t IT_0132 = IT_0123*IT_0131;
    const ccomplex_t IT_0133 = IT_0105*IT_0123;
    const ccomplex_t IT_0134 = IT_0120 + IT_0124 + IT_0132 + IT_0133;
    const ccomplex_t IT_0135 = conj(IT_0107)*IT_0108;
    const ccomplex_t IT_0136 = IT_0108*conj(IT_0111);
    const ccomplex_t IT_0137 = IT_0135 + IT_0136;
    const ccomplex_t IT_0138 = IT_0077 + IT_0080 + IT_0083;
    const ccomplex_t IT_0139 = m_s*IT_0138;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0123*IT_0140;
    const ccomplex_t IT_0142 = IT_0119 + IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0144 = cpow(IT_0143, 2);
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*m_sG*IT_0144;
    const ccomplex_t IT_0146 = pow(m_sG, 2);
    const ccomplex_t IT_0147 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -IT_0146 +
       reg_prop, -1);
    const ccomplex_t IT_0148 = IT_0145*IT_0147;
    const ccomplex_t IT_0149 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -IT_0146 +
       reg_prop, -1);
    const ccomplex_t IT_0150 = IT_0145*IT_0149;
    const ccomplex_t IT_0151 = cpow(IT_0065, 2);
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*m_N_1*IT_0151;
    const ccomplex_t IT_0153 = IT_0085*IT_0152;
    const ccomplex_t IT_0154 = 0.5*IT_0148 + (-0.166666666666667)*IT_0150 +
       IT_0153;
    const ccomplex_t IT_0155 = IT_0123*IT_0154;
    const ccomplex_t IT_0156 = IT_0142 + IT_0155;
    const ccomplex_t IT_0157 = cpow(IT_0037, 2);
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*m_N_4*IT_0157;
    const ccomplex_t IT_0159 = IT_0029*IT_0158;
    const ccomplex_t IT_0160 = cpow(IT_0054, 2);
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*m_N_3*IT_0160;
    const ccomplex_t IT_0162 = IT_0044*IT_0161;
    const ccomplex_t IT_0163 = IT_0159 + IT_0162;
    const ccomplex_t IT_0164 = IT_0108*IT_0140;
    const ccomplex_t IT_0165 = IT_0108*IT_0154;
    const ccomplex_t IT_0166 = IT_0075*IT_0158;
    const ccomplex_t IT_0167 = IT_0078*IT_0161;
    const ccomplex_t IT_0168 = IT_0166 + IT_0167;
    const ccomplex_t IT_0169 = IT_0108*IT_0168;
    const ccomplex_t IT_0170 = IT_0118*IT_0123;
    const ccomplex_t IT_0171 = IT_0123*IT_0163;
    const ccomplex_t IT_0172 = IT_0164 + IT_0165 + IT_0169 + IT_0170 + IT_0171;
    const ccomplex_t IT_0173 = IT_0108*conj(IT_0140);
    const ccomplex_t IT_0174 = IT_0108*conj(IT_0154);
    const ccomplex_t IT_0175 = IT_0108*conj(IT_0168);
    const ccomplex_t IT_0176 = conj(IT_0118)*IT_0123;
    const ccomplex_t IT_0177 = IT_0173 + IT_0174 + IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = IT_0000*IT_0006;
    const ccomplex_t IT_0179 = (-18)*IT_0028;
    const ccomplex_t IT_0180 = (-18)*IT_0043;
    const ccomplex_t IT_0181 = (-18)*IT_0072;
    const ccomplex_t IT_0182 = (-6)*IT_0084;
    const ccomplex_t IT_0183 = IT_0179 + IT_0180 + IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = (-18)*conj(IT_0028);
    const ccomplex_t IT_0185 = (-18)*conj(IT_0043);
    const ccomplex_t IT_0186 = (-18)*conj(IT_0072);
    const ccomplex_t IT_0187 = (-6)*conj(IT_0084);
    const ccomplex_t IT_0188 = IT_0184 + IT_0185 + IT_0186 + IT_0187;
    const ccomplex_t IT_0189 = (-18)*IT_0084;
    const ccomplex_t IT_0190 = (-6)*IT_0028;
    const ccomplex_t IT_0191 = (-6)*IT_0043;
    const ccomplex_t IT_0192 = (-6)*IT_0072;
    const ccomplex_t IT_0193 = IT_0189 + IT_0190 + IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = (-18)*conj(IT_0084);
    const ccomplex_t IT_0195 = (-6)*conj(IT_0028);
    const ccomplex_t IT_0196 = (-6)*conj(IT_0043);
    const ccomplex_t IT_0197 = (-6)*conj(IT_0072);
    const ccomplex_t IT_0198 = IT_0194 + IT_0195 + IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = conj(IT_0077) + conj(IT_0087);
    const ccomplex_t IT_0200 = (-18)*IT_0097;
    const ccomplex_t IT_0201 = (-18)*IT_0098;
    const ccomplex_t IT_0202 = (-6)*IT_0092;
    const ccomplex_t IT_0203 = (-6)*IT_0093;
    const ccomplex_t IT_0204 = (-6)*IT_0094;
    const ccomplex_t IT_0205 = IT_0200 + IT_0201 + IT_0202 + IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = IT_0077 + IT_0087;
    const ccomplex_t IT_0207 = (-18)*conj(IT_0097);
    const ccomplex_t IT_0208 = (-18)*conj(IT_0098);
    const ccomplex_t IT_0209 = (-6)*conj(IT_0092);
    const ccomplex_t IT_0210 = (-6)*conj(IT_0093);
    const ccomplex_t IT_0211 = (-6)*conj(IT_0094);
    const ccomplex_t IT_0212 = IT_0207 + IT_0208 + IT_0209 + IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = IT_0140 + IT_0154 + IT_0168;
    const ccomplex_t IT_0214 = IT_0108*conj(IT_0118);
    const ccomplex_t IT_0215 = IT_0107*IT_0123;
    const ccomplex_t IT_0216 = conj(IT_0107)*IT_0123;
    const ccomplex_t IT_0217 = IT_0111*IT_0123;
    const ccomplex_t IT_0218 = IT_0108*conj(IT_0122);
    const ccomplex_t IT_0219 = conj(IT_0111)*IT_0123;
    const ccomplex_t IT_0220 = IT_0108*conj(IT_0115);
    const ccomplex_t IT_0221 = IT_0137 + IT_0220;
    const ccomplex_t IT_0222 = conj(IT_0122)*IT_0123;
    const ccomplex_t IT_0223 = IT_0108*IT_0122;
    const ccomplex_t IT_0224 = IT_0115*IT_0123;
    const ccomplex_t IT_0225 = IT_0108*IT_0131;
    const ccomplex_t IT_0226 = IT_0108*conj(IT_0131);
    const ccomplex_t IT_0227 = conj(IT_0115)*IT_0123;
    const ccomplex_t IT_0228 = IT_0214 + IT_0221;
    const ccomplex_t IT_0229 = IT_0123*conj(IT_0131);
    const ccomplex_t IT_0230 = cpow(IT_0024, 2);
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*m_N_2*IT_0230;
    const ccomplex_t IT_0232 = IT_0007*IT_0231;
    const ccomplex_t IT_0233 = IT_0044*IT_0126;
    const ccomplex_t IT_0234 = IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = conj(IT_0105)*IT_0108;
    const ccomplex_t IT_0236 = m_s*IT_0080;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = IT_0108*conj(IT_0237);
    const ccomplex_t IT_0239 = conj(IT_0105)*IT_0123;
    const ccomplex_t IT_0240 = IT_0123*IT_0237;
    const ccomplex_t IT_0241 = m_s*IT_0087;
    const ccomplex_t IT_0242 = -IT_0241;
    const ccomplex_t IT_0243 = IT_0108*conj(IT_0234);
    const ccomplex_t IT_0244 = IT_0123*conj(IT_0237);
    const ccomplex_t IT_0245 = IT_0105*IT_0108;
    const ccomplex_t IT_0246 = IT_0108*IT_0237;
    const ccomplex_t IT_0247 = IT_0123*IT_0234;
    const ccomplex_t IT_0248 = IT_0029*IT_0129;
    const ccomplex_t IT_0249 = IT_0108*conj(IT_0242);
    const ccomplex_t IT_0250 = cpow(IT_0067, 2);
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*m_N_1*IT_0250;
    const ccomplex_t IT_0252 = IT_0085*IT_0251;
    const ccomplex_t IT_0253 = IT_0108*conj(IT_0252);
    const ccomplex_t IT_0254 = IT_0081*IT_0231;
    const ccomplex_t IT_0255 = IT_0108*conj(IT_0254);
    const ccomplex_t IT_0256 = IT_0123*conj(IT_0234);
    const ccomplex_t IT_0257 = IT_0069*IT_0251;
    const ccomplex_t IT_0258 = conj(IT_0248) + conj(IT_0257);
    const ccomplex_t IT_0259 = IT_0108*IT_0234;
    const ccomplex_t IT_0260 = IT_0123*IT_0242;
    const ccomplex_t IT_0261 = IT_0123*conj(IT_0242);
    const ccomplex_t IT_0262 = conj(IT_0252) + conj(IT_0254);
    const ccomplex_t IT_0263 = IT_0123*conj(IT_0140);
    const ccomplex_t IT_0264 = IT_0123*conj(IT_0154);
    const ccomplex_t IT_0265 = IT_0123*IT_0168;
    const ccomplex_t IT_0266 = IT_0108*conj(IT_0163);
    const ccomplex_t IT_0267 = IT_0069*IT_0152;
    const ccomplex_t IT_0268 = -IT_0106 + -IT_0110 + -IT_0114 + (
      -0.166666666666667)*IT_0148 + 0.5*IT_0150 + IT_0267;
    const ccomplex_t IT_0269 = IT_0108*conj(IT_0268);
    const ccomplex_t IT_0270 = cpow(IT_0019, 2);
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*m_N_2*IT_0270;
    const ccomplex_t IT_0272 = IT_0007*IT_0271;
    const ccomplex_t IT_0273 = IT_0108*conj(IT_0272);
    const ccomplex_t IT_0274 = IT_0214 + IT_0266 + IT_0269 + IT_0273;
    const ccomplex_t IT_0275 = IT_0123*conj(IT_0168);
    const ccomplex_t IT_0276 = IT_0081*IT_0271;
    const ccomplex_t IT_0277 = conj(IT_0242) + conj(IT_0276);
    const ccomplex_t IT_0278 = IT_0123*conj(IT_0163);
    const ccomplex_t IT_0279 = conj(IT_0268) + conj(IT_0272);
    const ccomplex_t IT_0280 = IT_0123*IT_0268;
    const ccomplex_t IT_0281 = IT_0108*IT_0163;
    const ccomplex_t IT_0282 = IT_0108*IT_0268;
    const ccomplex_t IT_0283 = IT_0123*IT_0252;
    const ccomplex_t IT_0284 = IT_0108*IT_0242;
    const ccomplex_t IT_0285 = IT_0108*IT_0252;
    const ccomplex_t IT_0286 = IT_0108*IT_0254;
    const ccomplex_t IT_0287 = IT_0123*IT_0248;
    const ccomplex_t IT_0288 = m_s*s_24;
    const ccomplex_t IT_0289 = 6*IT_0122;
    const ccomplex_t IT_0290 = 6*IT_0131;
    const ccomplex_t IT_0291 = 6*IT_0105;
    const ccomplex_t IT_0292 = 6*IT_0237;
    const ccomplex_t IT_0293 = 6*IT_0242;
    const ccomplex_t IT_0294 = 18*IT_0107;
    const ccomplex_t IT_0295 = 18*IT_0111;
    const ccomplex_t IT_0296 = 18*IT_0115;
    const ccomplex_t IT_0297 = 18*IT_0118;
    const ccomplex_t IT_0298 = 18*IT_0234;
    const ccomplex_t IT_0299 = 6*conj(IT_0122);
    const ccomplex_t IT_0300 = 6*conj(IT_0131);
    const ccomplex_t IT_0301 = 6*conj(IT_0105);
    const ccomplex_t IT_0302 = 6*conj(IT_0237);
    const ccomplex_t IT_0303 = 6*conj(IT_0242);
    const ccomplex_t IT_0304 = 18*conj(IT_0107);
    const ccomplex_t IT_0305 = 18*conj(IT_0111);
    const ccomplex_t IT_0306 = 18*conj(IT_0115);
    const ccomplex_t IT_0307 = 18*conj(IT_0118);
    const ccomplex_t IT_0308 = 18*conj(IT_0234);
    const ccomplex_t IT_0309 = 6*IT_0107;
    const ccomplex_t IT_0310 = 6*IT_0111;
    const ccomplex_t IT_0311 = 6*IT_0115;
    const ccomplex_t IT_0312 = 6*IT_0118;
    const ccomplex_t IT_0313 = 6*IT_0234;
    const ccomplex_t IT_0314 = 18*IT_0122;
    const ccomplex_t IT_0315 = 18*IT_0131;
    const ccomplex_t IT_0316 = 18*IT_0105;
    const ccomplex_t IT_0317 = 18*IT_0237;
    const ccomplex_t IT_0318 = 18*IT_0242;
    const ccomplex_t IT_0319 = 6*conj(IT_0107);
    const ccomplex_t IT_0320 = 6*conj(IT_0111);
    const ccomplex_t IT_0321 = 6*conj(IT_0115);
    const ccomplex_t IT_0322 = 6*conj(IT_0118);
    const ccomplex_t IT_0323 = 6*conj(IT_0234);
    const ccomplex_t IT_0324 = 18*conj(IT_0122);
    const ccomplex_t IT_0325 = 18*conj(IT_0131);
    const ccomplex_t IT_0326 = 18*conj(IT_0105);
    const ccomplex_t IT_0327 = 18*conj(IT_0237);
    const ccomplex_t IT_0328 = 18*conj(IT_0242);
    const ccomplex_t IT_0329 = 6*IT_0140;
    const ccomplex_t IT_0330 = 6*IT_0154;
    const ccomplex_t IT_0331 = 6*IT_0168;
    const ccomplex_t IT_0332 = 18*IT_0163;
    const ccomplex_t IT_0333 = 18*IT_0268;
    const ccomplex_t IT_0334 = 6*conj(IT_0140);
    const ccomplex_t IT_0335 = 6*conj(IT_0154);
    const ccomplex_t IT_0336 = 6*conj(IT_0168);
    const ccomplex_t IT_0337 = 18*conj(IT_0163);
    const ccomplex_t IT_0338 = 18*conj(IT_0268);
    const ccomplex_t IT_0339 = IT_0272 + 3*IT_0276;
    const ccomplex_t IT_0340 = 6*IT_0163;
    const ccomplex_t IT_0341 = 6*IT_0268;
    const ccomplex_t IT_0342 = 18*IT_0140;
    const ccomplex_t IT_0343 = 18*IT_0154;
    const ccomplex_t IT_0344 = 18*IT_0168;
    const ccomplex_t IT_0345 = 3*conj(IT_0276);
    const ccomplex_t IT_0346 = conj(IT_0272) + IT_0345;
    const ccomplex_t IT_0347 = 6*conj(IT_0163);
    const ccomplex_t IT_0348 = 6*conj(IT_0268);
    const ccomplex_t IT_0349 = 18*conj(IT_0140);
    const ccomplex_t IT_0350 = 18*conj(IT_0154);
    const ccomplex_t IT_0351 = 18*conj(IT_0168);
    const ccomplex_t IT_0352 = 6*IT_0084;
    const ccomplex_t IT_0353 = 6*IT_0087;
    const ccomplex_t IT_0354 = 18*IT_0028;
    const ccomplex_t IT_0355 = 18*IT_0043;
    const ccomplex_t IT_0356 = 18*IT_0072;
    const ccomplex_t IT_0357 = IT_0248 + IT_0257;
    const ccomplex_t IT_0358 = 6*conj(IT_0084);
    const ccomplex_t IT_0359 = 6*conj(IT_0087);
    const ccomplex_t IT_0360 = 18*conj(IT_0028);
    const ccomplex_t IT_0361 = 18*conj(IT_0043);
    const ccomplex_t IT_0362 = 18*conj(IT_0072);
    const ccomplex_t IT_0363 = 6*IT_0028;
    const ccomplex_t IT_0364 = 6*IT_0043;
    const ccomplex_t IT_0365 = 6*IT_0072;
    const ccomplex_t IT_0366 = 18*IT_0084;
    const ccomplex_t IT_0367 = 18*IT_0087;
    const ccomplex_t IT_0368 = IT_0252 + IT_0254;
    const ccomplex_t IT_0369 = 6*conj(IT_0028);
    const ccomplex_t IT_0370 = 6*conj(IT_0043);
    const ccomplex_t IT_0371 = 6*conj(IT_0072);
    const ccomplex_t IT_0372 = 18*conj(IT_0084);
    const ccomplex_t IT_0373 = 18*conj(IT_0087);
    const ccomplex_t IT_0374 = IT_0074*(IT_0289 + IT_0290 + IT_0291 + IT_0292 
      + IT_0293 + IT_0294 + IT_0295 + IT_0296 + IT_0297 + IT_0298) + IT_0073*
      (IT_0299 + IT_0300 + IT_0301 + IT_0302 + IT_0303 + IT_0304 + IT_0305 +
       IT_0306 + IT_0307 + IT_0308) + IT_0090*(IT_0309 + IT_0310 + IT_0311 +
       IT_0312 + IT_0313 + IT_0314 + IT_0315 + IT_0316 + IT_0317 + IT_0318) + 6
      *conj(IT_0077)*(IT_0248 + 3*IT_0252 + 3*IT_0254 + IT_0257 + 1./6*IT_0309 +
       1./6*IT_0310 + 1./6*IT_0311 + 1./6*IT_0312 + 1./6*IT_0313 + 1./6*IT_0314 
      + 1./6*IT_0315 + 1./6*IT_0316 + 1./6*IT_0317 + 1./6*IT_0318) + IT_0088*
      (IT_0319 + IT_0320 + IT_0321 + IT_0322 + IT_0323 + IT_0324 + IT_0325 +
       IT_0326 + IT_0327 + IT_0328) + 6*IT_0077*(conj(IT_0248) + 3*conj(IT_0252)
       + 3*conj(IT_0254) + conj(IT_0257) + 1./6*IT_0319 + 1./6*IT_0320 + 1./6
      *IT_0321 + 1./6*IT_0322 + 1./6*IT_0323 + 1./6*IT_0324 + 1./6*IT_0325 + 1.
      /6*IT_0326 + 1./6*IT_0327 + 1./6*IT_0328) + IT_0096*(18*IT_0272 + 6
      *IT_0276 + IT_0293 + IT_0297 + IT_0329 + IT_0330 + IT_0331 + IT_0332 +
       IT_0333) + IT_0095*(18*conj(IT_0272) + 6*conj(IT_0276) + IT_0303 +
       IT_0307 + IT_0334 + IT_0335 + IT_0336 + IT_0337 + IT_0338) + 6*IT_0099*
      (0.166666666666667*IT_0312 + 0.166666666666667*IT_0318 + IT_0339 +
       0.166666666666667*IT_0340 + 0.166666666666667*IT_0341 + 0.166666666666667
      *IT_0342 + 0.166666666666667*IT_0343 + 0.166666666666667*IT_0344) + 6
      *IT_0103*(0.166666666666667*IT_0322 + 0.166666666666667*IT_0328 + IT_0346 
      + 0.166666666666667*IT_0347 + 0.166666666666667*IT_0348 +
       0.166666666666667*IT_0349 + 0.166666666666667*IT_0350 + 0.166666666666667
      *IT_0351) + IT_0258*(IT_0352 + IT_0353 + IT_0354 + IT_0355 + IT_0356) +
       IT_0357*(IT_0358 + IT_0359 + IT_0360 + IT_0361 + IT_0362) + IT_0262*
      (IT_0363 + IT_0364 + IT_0365 + IT_0366 + IT_0367) + IT_0368*(IT_0369 +
       IT_0370 + IT_0371 + IT_0372 + IT_0373);
    const ccomplex_t IT_0375 = (-18)*conj(IT_0107);
    const ccomplex_t IT_0376 = (-18)*conj(IT_0111);
    const ccomplex_t IT_0377 = (-18)*conj(IT_0115);
    const ccomplex_t IT_0378 = (-18)*conj(IT_0118);
    const ccomplex_t IT_0379 = (-18)*conj(IT_0234);
    const ccomplex_t IT_0380 = (-6)*conj(IT_0122);
    const ccomplex_t IT_0381 = (-6)*conj(IT_0131);
    const ccomplex_t IT_0382 = (-6)*conj(IT_0105);
    const ccomplex_t IT_0383 = (-6)*conj(IT_0237);
    const ccomplex_t IT_0384 = (-6)*conj(IT_0242);
    const ccomplex_t IT_0385 = IT_0375 + IT_0376 + IT_0377 + IT_0378 + IT_0379
       + IT_0380 + IT_0381 + IT_0382 + IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = (-18)*IT_0122;
    const ccomplex_t IT_0387 = (-18)*IT_0131;
    const ccomplex_t IT_0388 = (-18)*IT_0105;
    const ccomplex_t IT_0389 = (-18)*IT_0237;
    const ccomplex_t IT_0390 = (-18)*IT_0242;
    const ccomplex_t IT_0391 = (-6)*IT_0107;
    const ccomplex_t IT_0392 = (-6)*IT_0111;
    const ccomplex_t IT_0393 = (-6)*IT_0115;
    const ccomplex_t IT_0394 = (-6)*IT_0118;
    const ccomplex_t IT_0395 = (-6)*IT_0234;
    const ccomplex_t IT_0396 = IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390
       + IT_0391 + IT_0392 + IT_0393 + IT_0394 + IT_0395;
    const ccomplex_t IT_0397 = (-18)*conj(IT_0122);
    const ccomplex_t IT_0398 = (-18)*conj(IT_0131);
    const ccomplex_t IT_0399 = (-18)*conj(IT_0105);
    const ccomplex_t IT_0400 = (-18)*conj(IT_0237);
    const ccomplex_t IT_0401 = (-18)*conj(IT_0242);
    const ccomplex_t IT_0402 = (-6)*conj(IT_0107);
    const ccomplex_t IT_0403 = (-6)*conj(IT_0111);
    const ccomplex_t IT_0404 = (-6)*conj(IT_0115);
    const ccomplex_t IT_0405 = (-6)*conj(IT_0118);
    const ccomplex_t IT_0406 = (-6)*conj(IT_0234);
    const ccomplex_t IT_0407 = IT_0397 + IT_0398 + IT_0399 + IT_0400 + IT_0401
       + IT_0402 + IT_0403 + IT_0404 + IT_0405 + IT_0406;
    const ccomplex_t IT_0408 = (-18)*IT_0118;
    const ccomplex_t IT_0409 = (-18)*IT_0163;
    const ccomplex_t IT_0410 = (-18)*IT_0268;
    const ccomplex_t IT_0411 = (-6)*IT_0140;
    const ccomplex_t IT_0412 = (-6)*IT_0154;
    const ccomplex_t IT_0413 = (-6)*IT_0168;
    const ccomplex_t IT_0414 = (-6)*IT_0242;
    const ccomplex_t IT_0415 = IT_0408 + IT_0409 + IT_0410 + IT_0411 + IT_0412
       + IT_0413 + IT_0414;
    const ccomplex_t IT_0416 = (-18)*conj(IT_0163);
    const ccomplex_t IT_0417 = (-18)*conj(IT_0268);
    const ccomplex_t IT_0418 = (-6)*conj(IT_0140);
    const ccomplex_t IT_0419 = (-6)*conj(IT_0154);
    const ccomplex_t IT_0420 = (-6)*conj(IT_0168);
    const ccomplex_t IT_0421 = IT_0378 + IT_0384 + IT_0416 + IT_0417 + IT_0418
       + IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = (-18)*IT_0140;
    const ccomplex_t IT_0423 = (-18)*IT_0154;
    const ccomplex_t IT_0424 = (-18)*IT_0168;
    const ccomplex_t IT_0425 = (-6)*IT_0163;
    const ccomplex_t IT_0426 = (-6)*IT_0268;
    const ccomplex_t IT_0427 = IT_0390 + IT_0394 + IT_0422 + IT_0423 + IT_0424
       + IT_0425 + IT_0426;
    const ccomplex_t IT_0428 = (-18)*conj(IT_0140);
    const ccomplex_t IT_0429 = (-18)*conj(IT_0154);
    const ccomplex_t IT_0430 = (-18)*conj(IT_0168);
    const ccomplex_t IT_0431 = (-6)*conj(IT_0163);
    const ccomplex_t IT_0432 = (-6)*conj(IT_0268);
    const ccomplex_t IT_0433 = IT_0401 + IT_0405 + IT_0428 + IT_0429 + IT_0430
       + IT_0431 + IT_0432;
    const ccomplex_t IT_0434 = (-18)*IT_0107;
    const ccomplex_t IT_0435 = (-18)*IT_0111;
    const ccomplex_t IT_0436 = (-18)*IT_0115;
    const ccomplex_t IT_0437 = IT_0434 + IT_0435 + IT_0436;
    const ccomplex_t IT_0438 = (-6)*IT_0122;
    const ccomplex_t IT_0439 = (-6)*IT_0131;
    const ccomplex_t IT_0440 = (-18)*IT_0234;
    const ccomplex_t IT_0441 = (-6)*IT_0105;
    const ccomplex_t IT_0442 = (-6)*IT_0237;
    const ccomplex_t IT_0443 = IT_0248 + IT_0257 + 0.333333333333333*IT_0368 +
       (-0.0555555555555556)*IT_0408 + (-0.0555555555555556)*IT_0414 + (
      -0.0555555555555556)*IT_0434 + (-0.0555555555555556)*IT_0435 + (
      -0.0555555555555556)*IT_0436 + (-0.0555555555555556)*IT_0438 + (
      -0.0555555555555556)*IT_0439 + (-0.0555555555555556)*IT_0440 + (
      -0.0555555555555556)*IT_0441 + (-0.0555555555555556)*IT_0442;
    const ccomplex_t IT_0444 = conj(IT_0248) + conj(IT_0257) +
       0.333333333333333*IT_0262 + (-0.0555555555555556)*IT_0375 + (
      -0.0555555555555556)*IT_0376 + (-0.0555555555555556)*IT_0377 + (
      -0.0555555555555556)*IT_0378 + (-0.0555555555555556)*IT_0379 + (
      -0.0555555555555556)*IT_0380 + (-0.0555555555555556)*IT_0381 + (
      -0.0555555555555556)*IT_0382 + (-0.0555555555555556)*IT_0383 + (
      -0.0555555555555556)*IT_0384;
    const ccomplex_t IT_0445 = IT_0234*IT_0378 + IT_0118*((-36)*conj(IT_0118) 
      + IT_0375 + IT_0376 + IT_0377 + IT_0380 + IT_0381) + (IT_0163 + IT_0268)
      *IT_0385 + (conj(IT_0105) + conj(IT_0237) + conj(IT_0242))*IT_0394 + (-6)
      *conj(IT_0276)*(IT_0248 + 3*IT_0252 + 3*IT_0254 + IT_0257 + -1./6*IT_0386 
      + -1./6*IT_0387 + -1./6*IT_0388 + -1./6*IT_0389 + -1./6*IT_0390 + -1./6
      *IT_0391 + -1./6*IT_0392 + -1./6*IT_0393 + -1./6*IT_0394 + -1./6*IT_0395) 
      + (conj(IT_0140) + conj(IT_0154) + conj(IT_0168) + conj(IT_0242))*IT_0396 
      + (IT_0105 + IT_0237 + IT_0242)*IT_0405 + (-6)*IT_0276*(conj(IT_0248) + 3
      *conj(IT_0252) + 3*conj(IT_0254) + conj(IT_0257) + -1./6*IT_0397 + -1./6
      *IT_0398 + -1./6*IT_0399 + -1./6*IT_0400 + -1./6*IT_0401 + -1./6*IT_0402 +
       -1./6*IT_0403 + -1./6*IT_0404 + -1./6*IT_0405 + -1./6*IT_0406) + (IT_0213
       + IT_0242)*IT_0407 + conj(IT_0234)*IT_0408 + IT_0258*IT_0415 + IT_0357
      *IT_0421 + IT_0262*IT_0427 + IT_0368*IT_0433 + conj(IT_0118)*(IT_0437 +
       IT_0438 + IT_0439) + (conj(IT_0163) + conj(IT_0268))*(IT_0408 + IT_0414 +
       IT_0437 + IT_0438 + IT_0439 + IT_0440 + IT_0441 + IT_0442) + (-18)*conj
      (IT_0272)*IT_0443 + (-18)*IT_0272*IT_0444;
    const ccomplex_t IT_0446 = m_s*s_23;
    const ccomplex_t IT_0447 = (-18)*IT_0092;
    const ccomplex_t IT_0448 = (-18)*IT_0093;
    const ccomplex_t IT_0449 = (-18)*IT_0094;
    const ccomplex_t IT_0450 = (-6)*IT_0097;
    const ccomplex_t IT_0451 = (-6)*IT_0098;
    const ccomplex_t IT_0452 = (-18)*conj(IT_0092);
    const ccomplex_t IT_0453 = (-18)*conj(IT_0093);
    const ccomplex_t IT_0454 = (-18)*conj(IT_0094);
    const ccomplex_t IT_0455 = (-6)*conj(IT_0097);
    const ccomplex_t IT_0456 = (-6)*conj(IT_0098);
    const ccomplex_t IT_0457 = IT_0205*IT_0262 + IT_0188*IT_0272 + IT_0183
      *conj(IT_0272) + IT_0198*IT_0276 + IT_0193*conj(IT_0276) + IT_0212*IT_0368
       + IT_0095*IT_0385 + IT_0099*IT_0396 + IT_0103*IT_0407 + IT_0074*IT_0415 +
       IT_0073*IT_0421 + (-6)*IT_0199*(IT_0339 + (-0.166666666666667)*IT_0390 + 
      (-0.166666666666667)*IT_0394 + (-0.166666666666667)*IT_0422 + (
      -0.166666666666667)*IT_0423 + (-0.166666666666667)*IT_0424 + (
      -0.166666666666667)*IT_0425 + (-0.166666666666667)*IT_0426) + conj(IT_0084
      )*IT_0427 + (-6)*IT_0206*(IT_0346 + (-0.166666666666667)*IT_0401 + (
      -0.166666666666667)*IT_0405 + (-0.166666666666667)*IT_0428 + (
      -0.166666666666667)*IT_0429 + (-0.166666666666667)*IT_0430 + (
      -0.166666666666667)*IT_0431 + (-0.166666666666667)*IT_0432) + IT_0084
      *IT_0433 + IT_0096*(IT_0408 + IT_0414 + IT_0437 + IT_0438 + IT_0439 +
       IT_0440 + IT_0441 + IT_0442) + IT_0258*(IT_0447 + IT_0448 + IT_0449 +
       IT_0450 + IT_0451) + IT_0357*(IT_0452 + IT_0453 + IT_0454 + IT_0455 +
       IT_0456);
    const ccomplex_t IT_0458 = IT_0005*(IT_0073*IT_0074 + IT_0089*IT_0091 +
       IT_0095*IT_0096 + (IT_0097 + IT_0098)*IT_0099) + IT_0102*(IT_0074*IT_0089
       + IT_0073*IT_0091 + IT_0095*IT_0099 + IT_0096*IT_0103) + conj(IT_0122)*
      (IT_0113 + IT_0124) + conj(IT_0131)*(IT_0113 + IT_0116 + IT_0124 + IT_0132
      ) + conj(IT_0105)*IT_0134 + IT_0122*IT_0137 + conj(IT_0140)*IT_0142 + conj
      (IT_0154)*IT_0156 + conj(IT_0163)*IT_0172 + IT_0163*IT_0177 + IT_0178*
      (IT_0096*IT_0183 + IT_0095*IT_0188 + IT_0099*IT_0193 + IT_0103*IT_0198 +
       IT_0199*IT_0205 + IT_0206*IT_0212) + IT_0213*IT_0214 + (conj(IT_0107) +
       conj(IT_0111))*IT_0215 + IT_0111*IT_0216 + conj(IT_0111)*IT_0217 +
       IT_0115*(IT_0216 + IT_0218 + IT_0219) + IT_0131*(IT_0221 + IT_0222) +
       conj(IT_0115)*(IT_0215 + IT_0217 + IT_0223 + IT_0224) + conj(IT_0118)*
      (IT_0215 + IT_0217 + IT_0223 + IT_0224 + IT_0225) + IT_0118*(36*s_34*conj
      (IT_0118) + IT_0216 + IT_0218 + IT_0219 + IT_0226 + IT_0227) + IT_0105*
      (IT_0222 + IT_0228 + IT_0229) + IT_0234*(IT_0176 + IT_0216 + IT_0218 +
       IT_0219 + IT_0226 + IT_0227 + IT_0235 + IT_0238) + IT_0237*(IT_0222 +
       IT_0228 + IT_0229 + IT_0239) + conj(IT_0237)*(IT_0134 + IT_0240) +
       IT_0242*(IT_0222 + IT_0228 + IT_0229 + IT_0239 + IT_0243 + IT_0244) +
       conj(IT_0234)*(IT_0170 + IT_0215 + IT_0217 + IT_0223 + IT_0224 + IT_0225 
      + IT_0245 + IT_0246 + IT_0247) + IT_0248*(IT_0176 + IT_0216 + IT_0218 +
       IT_0219 + IT_0226 + IT_0227 + IT_0235 + IT_0238 + IT_0249 + IT_0253 +
       IT_0255 + IT_0256) + IT_0257*(IT_0176 + IT_0216 + IT_0218 + IT_0219 +
       IT_0226 + IT_0227 + IT_0235 + IT_0238 + IT_0249 + IT_0253 + IT_0255 +
       IT_0256 + IT_0123*IT_0258) + conj(IT_0242)*(IT_0120 + IT_0124 + IT_0132 +
       IT_0133 + IT_0240 + IT_0259 + IT_0260) + IT_0252*(IT_0222 + IT_0228 +
       IT_0229 + IT_0239 + IT_0243 + IT_0244 + IT_0261) + IT_0254*(IT_0222 +
       IT_0228 + IT_0229 + IT_0239 + IT_0243 + IT_0244 + IT_0261 + IT_0123
      *IT_0262) + (IT_0154 + IT_0168)*IT_0263 + IT_0168*IT_0264 + conj(IT_0168)*
      (IT_0156 + IT_0265) + IT_0242*(IT_0263 + IT_0264 + IT_0274 + IT_0275) +
       IT_0276*(IT_0263 + IT_0264 + IT_0274 + IT_0275 + IT_0123*IT_0277) +
       IT_0268*(IT_0177 + IT_0278) + IT_0272*(IT_0173 + IT_0174 + IT_0175 +
       IT_0176 + IT_0108*IT_0277 + IT_0278 + IT_0123*IT_0279) + IT_0279*(IT_0172
       + IT_0280) + IT_0277*(IT_0119 + IT_0141 + IT_0155 + IT_0260 + IT_0265 +
       IT_0281 + IT_0282) + IT_0262*(IT_0120 + IT_0124 + IT_0132 + IT_0133 +
       IT_0240 + IT_0259 + IT_0260 + IT_0283) + IT_0258*(IT_0170 + IT_0215 +
       IT_0217 + IT_0223 + IT_0224 + IT_0225 + IT_0245 + IT_0246 + IT_0247 +
       IT_0284 + IT_0285 + IT_0286 + IT_0287) + IT_0288*IT_0374 + IT_0006
      *IT_0445 + IT_0446*IT_0457;
    return create_ccomplex_return(IT_0458);
}

