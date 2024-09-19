#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_sd_R_to_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_sd_R_to_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sd_R = param->m_sd_R;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_sd_R, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_23*s_24;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_d, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0009;
    const ccomplex_t IT_0011 = (-0.333333333333333)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0012*IT_0014*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0012*IT_0014*IT_0016;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0009;
    const ccomplex_t IT_0026 = (-0.333333333333333)*IT_0025;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = IT_0022*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0012*IT_0014*IT_0016;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0009;
    const ccomplex_t IT_0033 = (-0.333333333333333)*IT_0032;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0021 + IT_0029 + IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0019*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0034*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0027*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = -IT_0045 + -IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0012*IT_0014*IT_0016;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0009;
    const ccomplex_t IT_0054 = (-0.333333333333333)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = IT_0050*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = conj(IT_0042) + conj(IT_0049) + conj(IT_0058);
    const ccomplex_t IT_0060 = -IT_0045;
    const ccomplex_t IT_0061 = -IT_0048;
    const ccomplex_t IT_0062 = IT_0021 + IT_0029;
    const ccomplex_t IT_0063 = IT_0037 + IT_0062;
    const ccomplex_t IT_0064 = conj(IT_0037) + conj(IT_0062);
    const ccomplex_t IT_0065 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0055*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0058 + IT_0060 + IT_0061;
    const ccomplex_t IT_0069 = IT_0049 + IT_0058 + IT_0068;
    const ccomplex_t IT_0070 = conj(IT_0058) + conj(IT_0060) + conj(IT_0061);
    const ccomplex_t IT_0071 = conj(IT_0038)*(IT_0042 + IT_0049 + IT_0058) +
       IT_0038*IT_0059 + (conj(IT_0042) + conj(IT_0058) + conj(IT_0060) + conj
      (IT_0061))*IT_0063 + (IT_0042 + IT_0058 + IT_0060 + IT_0061)*IT_0064 +
       conj(IT_0067)*IT_0069 + IT_0067*(conj(IT_0049) + conj(IT_0058) + IT_0070);
    const ccomplex_t IT_0072 = (-18)*IT_0001;
    const ccomplex_t IT_0073 = 36*IT_0003;
    const ccomplex_t IT_0074 = IT_0072 + IT_0073;
    const ccomplex_t IT_0075 = IT_0049 + IT_0058;
    const ccomplex_t IT_0076 = conj(IT_0049) + conj(IT_0058);
    const ccomplex_t IT_0077 = IT_0042*(conj(IT_0042) + 0.5*conj(IT_0049) +
       conj(IT_0058) + 0.5*conj(IT_0060) + 0.5*conj(IT_0061)) + 0.5*(conj
      (IT_0037) + conj(IT_0038) + conj(IT_0062))*IT_0067 + 0.5*IT_0038*(conj
      (IT_0038) + conj(IT_0067)) + 0.5*IT_0063*(conj(IT_0037) + conj(IT_0062) +
       conj(IT_0067)) + 0.5*conj(IT_0042)*IT_0069 + 0.5*IT_0068*IT_0070 + 0.5
      *IT_0075*IT_0076;
    const ccomplex_t IT_0078 = m_d*IT_0021;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = m_d*IT_0045;
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = 12*s_34;
    const ccomplex_t IT_0083 = m_d*IT_0057;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = pow(m_sG, 2);
    const ccomplex_t IT_0087 = cpow((-2)*s_13 + IT_0000 + IT_0006 + -IT_0086 +
       reg_prop, -1);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0089 = cpow(IT_0088, 2);
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_sG*IT_0089;
    const ccomplex_t IT_0091 = IT_0087*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = 18*s_34;
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = 36*s_34;
    const ccomplex_t IT_0096 = IT_0079*IT_0095;
    const ccomplex_t IT_0097 = IT_0081*IT_0082 + IT_0085 + IT_0094 + IT_0096;
    const ccomplex_t IT_0098 = m_d*s_24;
    const ccomplex_t IT_0099 = 6*IT_0084;
    const ccomplex_t IT_0100 = conj(IT_0067)*IT_0098*IT_0099;
    const ccomplex_t IT_0101 = 6*conj(IT_0084);
    const ccomplex_t IT_0102 = IT_0067*IT_0098*IT_0101;
    const ccomplex_t IT_0103 = 6*IT_0079;
    const ccomplex_t IT_0104 = conj(IT_0042)*IT_0098*IT_0103;
    const ccomplex_t IT_0105 = 6*conj(IT_0079);
    const ccomplex_t IT_0106 = IT_0042*IT_0098*IT_0105;
    const ccomplex_t IT_0107 = 6*IT_0081;
    const ccomplex_t IT_0108 = conj(IT_0067)*IT_0098*IT_0107;
    const ccomplex_t IT_0109 = 6*conj(IT_0081);
    const ccomplex_t IT_0110 = IT_0067*IT_0098*IT_0109;
    const ccomplex_t IT_0111 = m_d*IT_0037;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = 6*IT_0112;
    const ccomplex_t IT_0114 = conj(IT_0042)*IT_0098*IT_0113;
    const ccomplex_t IT_0115 = 6*conj(IT_0112);
    const ccomplex_t IT_0116 = IT_0042*IT_0098*IT_0115;
    const ccomplex_t IT_0117 = 18*IT_0084;
    const ccomplex_t IT_0118 = conj(IT_0042)*IT_0098*IT_0117;
    const ccomplex_t IT_0119 = 18*conj(IT_0084);
    const ccomplex_t IT_0120 = IT_0042*IT_0098*IT_0119;
    const ccomplex_t IT_0121 = 18*IT_0079;
    const ccomplex_t IT_0122 = conj(IT_0067)*IT_0098*IT_0121;
    const ccomplex_t IT_0123 = 18*conj(IT_0079);
    const ccomplex_t IT_0124 = IT_0067*IT_0098*IT_0123;
    const ccomplex_t IT_0125 = 18*IT_0081;
    const ccomplex_t IT_0126 = conj(IT_0042)*IT_0098*IT_0125;
    const ccomplex_t IT_0127 = 18*conj(IT_0081);
    const ccomplex_t IT_0128 = IT_0042*IT_0098*IT_0127;
    const ccomplex_t IT_0129 = 18*IT_0112;
    const ccomplex_t IT_0130 = conj(IT_0067)*IT_0098*IT_0129;
    const ccomplex_t IT_0131 = 18*conj(IT_0112);
    const ccomplex_t IT_0132 = IT_0067*IT_0098*IT_0131;
    const ccomplex_t IT_0133 = cpow(IT_0092, 2);
    const ccomplex_t IT_0134 = (-18)*s_34;
    const ccomplex_t IT_0135 = IT_0082*conj(IT_0084);
    const ccomplex_t IT_0136 = IT_0092*IT_0134;
    const ccomplex_t IT_0137 = m_d*s_23;
    const ccomplex_t IT_0138 = (-18)*IT_0084;
    const ccomplex_t IT_0139 = conj(IT_0042)*IT_0137*IT_0138;
    const ccomplex_t IT_0140 = (-18)*conj(IT_0084);
    const ccomplex_t IT_0141 = IT_0042*IT_0137*IT_0140;
    const ccomplex_t IT_0142 = (-18)*IT_0079;
    const ccomplex_t IT_0143 = conj(IT_0067)*IT_0137*IT_0142;
    const ccomplex_t IT_0144 = (-18)*conj(IT_0079);
    const ccomplex_t IT_0145 = IT_0067*IT_0137*IT_0144;
    const ccomplex_t IT_0146 = (-18)*IT_0081;
    const ccomplex_t IT_0147 = conj(IT_0042)*IT_0137*IT_0146;
    const ccomplex_t IT_0148 = (-18)*conj(IT_0081);
    const ccomplex_t IT_0149 = IT_0042*IT_0137*IT_0148;
    const ccomplex_t IT_0150 = (-18)*IT_0112;
    const ccomplex_t IT_0151 = conj(IT_0067)*IT_0137*IT_0150;
    const ccomplex_t IT_0152 = (-18)*conj(IT_0112);
    const ccomplex_t IT_0153 = IT_0067*IT_0137*IT_0152;
    const ccomplex_t IT_0154 = (-6)*s_34;
    const ccomplex_t IT_0155 = IT_0092*IT_0154;
    const ccomplex_t IT_0156 = cpow((-2)*s_23 + IT_0000 + IT_0006 + -IT_0086 +
       reg_prop, -1);
    const ccomplex_t IT_0157 = IT_0090*IT_0156;
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = IT_0093*IT_0158;
    const ccomplex_t IT_0160 = cpow(IT_0054, 2);
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*m_N_1*IT_0160;
    const ccomplex_t IT_0162 = IT_0050*IT_0161;
    const ccomplex_t IT_0163 = cpow(IT_0026, 2);
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*m_N_3*IT_0163;
    const ccomplex_t IT_0165 = IT_0046*IT_0164;
    const ccomplex_t IT_0166 = -IT_0162 + -IT_0165;
    const ccomplex_t IT_0167 = IT_0065*IT_0161;
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = 6*s_34;
    const ccomplex_t IT_0170 = IT_0168*IT_0169;
    const ccomplex_t IT_0171 = IT_0093*IT_0166;
    const ccomplex_t IT_0172 = IT_0084*IT_0134;
    const ccomplex_t IT_0173 = IT_0081*IT_0134;
    const ccomplex_t IT_0174 = IT_0079*IT_0154;
    const ccomplex_t IT_0175 = IT_0079*IT_0169;
    const ccomplex_t IT_0176 = IT_0084*IT_0093;
    const ccomplex_t IT_0177 = conj(IT_0084)*IT_0134;
    const ccomplex_t IT_0178 = conj(IT_0079)*IT_0154;
    const ccomplex_t IT_0179 = IT_0092*IT_0169;
    const ccomplex_t IT_0180 = conj(IT_0168)*IT_0169;
    const ccomplex_t IT_0181 = IT_0134*IT_0158;
    const ccomplex_t IT_0182 = conj(IT_0081)*IT_0134;
    const ccomplex_t IT_0183 = IT_0158*IT_0169;
    const ccomplex_t IT_0184 = IT_0093*IT_0168;
    const ccomplex_t IT_0185 = IT_0079*IT_0134;
    const ccomplex_t IT_0186 = IT_0084*IT_0154;
    const ccomplex_t IT_0187 = IT_0081*IT_0154;
    const ccomplex_t IT_0188 = cpow(IT_0033, 2);
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*m_N_4*IT_0188;
    const ccomplex_t IT_0190 = IT_0035*IT_0189;
    const ccomplex_t IT_0191 = 0.166666666666667*IT_0157 + -IT_0190;
    const ccomplex_t IT_0192 = IT_0166*IT_0169;
    const ccomplex_t IT_0193 = IT_0093*IT_0191;
    const ccomplex_t IT_0194 = cpow(IT_0018, 2);
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*m_N_2*IT_0194;
    const ccomplex_t IT_0196 = IT_0007*IT_0195;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = cpow(IT_0052, 2);
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*m_N_1*IT_0198;
    const ccomplex_t IT_0200 = IT_0050*IT_0199;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0112*IT_0134;
    const ccomplex_t IT_0203 = cpow(IT_0011, 2);
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*m_N_2*IT_0203;
    const ccomplex_t IT_0205 = IT_0007*IT_0204;
    const ccomplex_t IT_0206 = -IT_0205;
    const ccomplex_t IT_0207 = conj(IT_0079)*IT_0134;
    const ccomplex_t IT_0208 = conj(IT_0112)*IT_0134;
    const ccomplex_t IT_0209 = conj(IT_0084)*IT_0154;
    const ccomplex_t IT_0210 = conj(IT_0081)*IT_0154;
    const ccomplex_t IT_0211 = IT_0154*IT_0158;
    const ccomplex_t IT_0212 = conj(IT_0166)*IT_0169;
    const ccomplex_t IT_0213 = IT_0093*conj(IT_0168);
    const ccomplex_t IT_0214 = IT_0043*IT_0189;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = conj(IT_0191) + conj(IT_0206);
    const ccomplex_t IT_0217 = IT_0029 + IT_0067;
    const ccomplex_t IT_0218 = m_d*IT_0217;
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = IT_0022*IT_0164;
    const ccomplex_t IT_0221 = conj(IT_0219) + conj(IT_0220);
    const ccomplex_t IT_0222 = IT_0134*IT_0168;
    const ccomplex_t IT_0223 = IT_0134*IT_0191;
    const ccomplex_t IT_0224 = IT_0154*IT_0166;
    const ccomplex_t IT_0225 = IT_0134*conj(IT_0168);
    const ccomplex_t IT_0226 = IT_0134*conj(IT_0191);
    const ccomplex_t IT_0227 = IT_0154*conj(IT_0166);
    const ccomplex_t IT_0228 = IT_0093*IT_0201;
    const ccomplex_t IT_0229 = IT_0112*IT_0154;
    const ccomplex_t IT_0230 = cpow(IT_0031, 2);
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*m_N_4*IT_0230;
    const ccomplex_t IT_0232 = IT_0043*IT_0231;
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0039*IT_0195;
    const ccomplex_t IT_0235 = cpow(IT_0024, 2);
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*m_N_3*IT_0235;
    const ccomplex_t IT_0237 = IT_0046*IT_0236;
    const ccomplex_t IT_0238 = -IT_0234 + -IT_0237;
    const ccomplex_t IT_0239 = conj(IT_0233) + conj(IT_0238);
    const ccomplex_t IT_0240 = IT_0169*IT_0197;
    const ccomplex_t IT_0241 = IT_0093*IT_0233;
    const ccomplex_t IT_0242 = conj(IT_0112)*IT_0154;
    const ccomplex_t IT_0243 = IT_0169*conj(IT_0197);
    const ccomplex_t IT_0244 = IT_0093*conj(IT_0201);
    const ccomplex_t IT_0245 = IT_0084*IT_0169;
    const ccomplex_t IT_0246 = IT_0081*IT_0169;
    const ccomplex_t IT_0247 = IT_0079*IT_0093;
    const ccomplex_t IT_0248 = IT_0093*IT_0112;
    const ccomplex_t IT_0249 = IT_0093*IT_0219;
    const ccomplex_t IT_0250 = IT_0154*IT_0215;
    const ccomplex_t IT_0251 = conj(IT_0084)*IT_0169;
    const ccomplex_t IT_0252 = conj(IT_0081)*IT_0169;
    const ccomplex_t IT_0253 = conj(IT_0079)*IT_0093;
    const ccomplex_t IT_0254 = IT_0093*conj(IT_0112);
    const ccomplex_t IT_0255 = IT_0154*conj(IT_0215);
    const ccomplex_t IT_0256 = IT_0039*IT_0204;
    const ccomplex_t IT_0257 = m_d*IT_0041;
    const ccomplex_t IT_0258 = IT_0256 + -IT_0257;
    const ccomplex_t IT_0259 = m_d*IT_0048;
    const ccomplex_t IT_0260 = (-0.166666666666667)*IT_0091 + -IT_0259;
    const ccomplex_t IT_0261 = conj(IT_0258) + conj(IT_0260);
    const ccomplex_t IT_0262 = IT_0112*IT_0169;
    const ccomplex_t IT_0263 = IT_0169*IT_0219;
    const ccomplex_t IT_0264 = IT_0081*IT_0093;
    const ccomplex_t IT_0265 = IT_0093*IT_0258;
    const ccomplex_t IT_0266 = IT_0134*IT_0166;
    const ccomplex_t IT_0267 = IT_0134*IT_0215;
    const ccomplex_t IT_0268 = IT_0154*IT_0168;
    const ccomplex_t IT_0269 = IT_0154*IT_0191;
    const ccomplex_t IT_0270 = IT_0154*IT_0206;
    const ccomplex_t IT_0271 = -IT_0259;
    const ccomplex_t IT_0272 = -IT_0257;
    const ccomplex_t IT_0273 = conj(IT_0271) + conj(IT_0272);
    const ccomplex_t IT_0274 = IT_0093*IT_0271;
    const ccomplex_t IT_0275 = IT_0134*IT_0201;
    const ccomplex_t IT_0276 = IT_0134*IT_0233;
    const ccomplex_t IT_0277 = IT_0134*IT_0238;
    const ccomplex_t IT_0278 = IT_0154*IT_0197;
    const ccomplex_t IT_0279 = conj(IT_0079)*IT_0169;
    const ccomplex_t IT_0280 = conj(IT_0112)*IT_0169;
    const ccomplex_t IT_0281 = IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = conj(IT_0084)*IT_0093;
    const ccomplex_t IT_0283 = conj(IT_0081)*IT_0093;
    const ccomplex_t IT_0284 = IT_0134*conj(IT_0201);
    const ccomplex_t IT_0285 = IT_0134*conj(IT_0233);
    const ccomplex_t IT_0286 = IT_0134*conj(IT_0238);
    const ccomplex_t IT_0287 = IT_0154*conj(IT_0197);
    const ccomplex_t IT_0288 = m_d*IT_0029;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = IT_0169*IT_0271;
    const ccomplex_t IT_0291 = IT_0169*IT_0272;
    const ccomplex_t IT_0292 = IT_0093*IT_0289;
    const ccomplex_t IT_0293 = IT_0134*IT_0197;
    const ccomplex_t IT_0294 = IT_0154*IT_0201;
    const ccomplex_t IT_0295 = IT_0154*IT_0233;
    const ccomplex_t IT_0296 = IT_0154*IT_0238;
    const ccomplex_t IT_0297 = m_d*IT_0067;
    const ccomplex_t IT_0298 = -IT_0297;
    const ccomplex_t IT_0299 = IT_0093*IT_0298;
    const ccomplex_t IT_0300 = IT_0022*IT_0236;
    const ccomplex_t IT_0301 = IT_0093*IT_0300;
    const ccomplex_t IT_0302 = IT_0035*IT_0231;
    const ccomplex_t IT_0303 = IT_0065*IT_0199;
    const ccomplex_t IT_0304 = conj(IT_0302) + conj(IT_0303);
    const ccomplex_t IT_0305 = IT_0093*IT_0302;
    const ccomplex_t IT_0306 = IT_0169*conj(IT_0271);
    const ccomplex_t IT_0307 = IT_0169*conj(IT_0272);
    const ccomplex_t IT_0308 = IT_0134*conj(IT_0197);
    const ccomplex_t IT_0309 = IT_0154*conj(IT_0201);
    const ccomplex_t IT_0310 = IT_0154*conj(IT_0233);
    const ccomplex_t IT_0311 = IT_0154*conj(IT_0238);
    const ccomplex_t IT_0312 = IT_0093*conj(IT_0289);
    const ccomplex_t IT_0313 = IT_0093*conj(IT_0298);
    const ccomplex_t IT_0314 = IT_0093*conj(IT_0300);
    const ccomplex_t IT_0315 = (-6)*IT_0084;
    const ccomplex_t IT_0316 = conj(IT_0067)*IT_0137*IT_0315;
    const ccomplex_t IT_0317 = (-6)*conj(IT_0084);
    const ccomplex_t IT_0318 = IT_0067*IT_0137*IT_0317;
    const ccomplex_t IT_0319 = (-6)*IT_0079;
    const ccomplex_t IT_0320 = conj(IT_0042)*IT_0137*IT_0319;
    const ccomplex_t IT_0321 = (-6)*conj(IT_0079);
    const ccomplex_t IT_0322 = IT_0042*IT_0137*IT_0321;
    const ccomplex_t IT_0323 = (-6)*IT_0081;
    const ccomplex_t IT_0324 = conj(IT_0067)*IT_0137*IT_0323;
    const ccomplex_t IT_0325 = (-6)*conj(IT_0081);
    const ccomplex_t IT_0326 = IT_0067*IT_0137*IT_0325;
    const ccomplex_t IT_0327 = (-6)*IT_0112;
    const ccomplex_t IT_0328 = conj(IT_0042)*IT_0137*IT_0327;
    const ccomplex_t IT_0329 = (-6)*conj(IT_0112);
    const ccomplex_t IT_0330 = IT_0042*IT_0137*IT_0329;
    const ccomplex_t IT_0331 = IT_0000*IT_0006;
    const ccomplex_t IT_0332 = (-6)*IT_0042;
    const ccomplex_t IT_0333 = conj(IT_0067)*IT_0331*IT_0332;
    const ccomplex_t IT_0334 = (-6)*conj(IT_0042);
    const ccomplex_t IT_0335 = IT_0067*IT_0331*IT_0334;
    const ccomplex_t IT_0336 = (-18)*IT_0038;
    const ccomplex_t IT_0337 = (-6)*IT_0049;
    const ccomplex_t IT_0338 = (-6)*IT_0058;
    const ccomplex_t IT_0339 = IT_0332 + IT_0336 + IT_0337 + IT_0338;
    const ccomplex_t IT_0340 = (-18)*conj(IT_0038);
    const ccomplex_t IT_0341 = (-6)*conj(IT_0049);
    const ccomplex_t IT_0342 = (-6)*conj(IT_0058);
    const ccomplex_t IT_0343 = IT_0334 + IT_0340 + IT_0341 + IT_0342;
    const ccomplex_t IT_0344 = 6*IT_0168;
    const ccomplex_t IT_0345 = 6*IT_0191;
    const ccomplex_t IT_0346 = 6*IT_0206;
    const ccomplex_t IT_0347 = 18*IT_0158;
    const ccomplex_t IT_0348 = 18*IT_0166;
    const ccomplex_t IT_0349 = 18*IT_0215;
    const ccomplex_t IT_0350 = (-18)*IT_0258;
    const ccomplex_t IT_0351 = (-18)*IT_0260;
    const ccomplex_t IT_0352 = (-6)*IT_0092;
    const ccomplex_t IT_0353 = (-6)*IT_0219;
    const ccomplex_t IT_0354 = 6*IT_0092;
    const ccomplex_t IT_0355 = 6*conj(IT_0168);
    const ccomplex_t IT_0356 = 6*conj(IT_0191);
    const ccomplex_t IT_0357 = 6*conj(IT_0206);
    const ccomplex_t IT_0358 = 18*conj(IT_0166);
    const ccomplex_t IT_0359 = 18*conj(IT_0215);
    const ccomplex_t IT_0360 = (-18)*IT_0158;
    const ccomplex_t IT_0361 = (-18)*conj(IT_0258);
    const ccomplex_t IT_0362 = (-18)*conj(IT_0260);
    const ccomplex_t IT_0363 = (-6)*conj(IT_0219);
    const ccomplex_t IT_0364 = 6*IT_0158;
    const ccomplex_t IT_0365 = 6*IT_0166;
    const ccomplex_t IT_0366 = 6*IT_0215;
    const ccomplex_t IT_0367 = 18*IT_0168;
    const ccomplex_t IT_0368 = 18*IT_0191;
    const ccomplex_t IT_0369 = 18*IT_0206;
    const ccomplex_t IT_0370 = (-18)*IT_0092;
    const ccomplex_t IT_0371 = (-18)*IT_0219;
    const ccomplex_t IT_0372 = (-6)*IT_0258;
    const ccomplex_t IT_0373 = (-6)*IT_0260;
    const ccomplex_t IT_0374 = 6*conj(IT_0166);
    const ccomplex_t IT_0375 = 6*conj(IT_0215);
    const ccomplex_t IT_0376 = 18*IT_0092;
    const ccomplex_t IT_0377 = 18*conj(IT_0168);
    const ccomplex_t IT_0378 = 18*conj(IT_0191);
    const ccomplex_t IT_0379 = 18*conj(IT_0206);
    const ccomplex_t IT_0380 = (-18)*conj(IT_0219);
    const ccomplex_t IT_0381 = (-6)*IT_0158;
    const ccomplex_t IT_0382 = (-6)*conj(IT_0258);
    const ccomplex_t IT_0383 = (-6)*conj(IT_0260);
    const ccomplex_t IT_0384 = 6*IT_0201;
    const ccomplex_t IT_0385 = 6*IT_0233;
    const ccomplex_t IT_0386 = 6*IT_0238;
    const ccomplex_t IT_0387 = 18*IT_0197;
    const ccomplex_t IT_0388 = (-18)*IT_0289;
    const ccomplex_t IT_0389 = (-18)*IT_0298;
    const ccomplex_t IT_0390 = (-6)*IT_0271;
    const ccomplex_t IT_0391 = (-6)*IT_0272;
    const ccomplex_t IT_0392 = (-18)*IT_0220;
    const ccomplex_t IT_0393 = 6*conj(IT_0201);
    const ccomplex_t IT_0394 = 6*conj(IT_0233);
    const ccomplex_t IT_0395 = 6*conj(IT_0238);
    const ccomplex_t IT_0396 = 18*conj(IT_0197);
    const ccomplex_t IT_0397 = (-18)*conj(IT_0289);
    const ccomplex_t IT_0398 = (-18)*conj(IT_0298);
    const ccomplex_t IT_0399 = (-6)*conj(IT_0271);
    const ccomplex_t IT_0400 = (-6)*conj(IT_0272);
    const ccomplex_t IT_0401 = (-18)*conj(IT_0220);
    const ccomplex_t IT_0402 = 6*IT_0197;
    const ccomplex_t IT_0403 = 18*IT_0201;
    const ccomplex_t IT_0404 = 18*IT_0233;
    const ccomplex_t IT_0405 = 18*IT_0238;
    const ccomplex_t IT_0406 = (-18)*IT_0271;
    const ccomplex_t IT_0407 = (-18)*IT_0272;
    const ccomplex_t IT_0408 = (-6)*IT_0289;
    const ccomplex_t IT_0409 = (-6)*IT_0298;
    const ccomplex_t IT_0410 = 6*conj(IT_0197);
    const ccomplex_t IT_0411 = 18*conj(IT_0201);
    const ccomplex_t IT_0412 = 18*conj(IT_0233);
    const ccomplex_t IT_0413 = 18*conj(IT_0238);
    const ccomplex_t IT_0414 = (-18)*conj(IT_0271);
    const ccomplex_t IT_0415 = (-18)*conj(IT_0272);
    const ccomplex_t IT_0416 = (-6)*conj(IT_0289);
    const ccomplex_t IT_0417 = (-6)*conj(IT_0298);
    const ccomplex_t IT_0418 = conj(IT_0300) + conj(IT_0302) + conj(IT_0303);
    const ccomplex_t IT_0419 = (-18)*IT_0062;
    const ccomplex_t IT_0420 = (-18)*IT_0067;
    const ccomplex_t IT_0421 = (-18)*IT_0037;
    const ccomplex_t IT_0422 = (-6)*IT_0060;
    const ccomplex_t IT_0423 = (-6)*IT_0061;
    const ccomplex_t IT_0424 = IT_0300 + IT_0302 + IT_0303;
    const ccomplex_t IT_0425 = (-18)*conj(IT_0062);
    const ccomplex_t IT_0426 = (-18)*conj(IT_0067);
    const ccomplex_t IT_0427 = (-18)*conj(IT_0037);
    const ccomplex_t IT_0428 = (-6)*conj(IT_0060);
    const ccomplex_t IT_0429 = (-6)*conj(IT_0061);
    const ccomplex_t IT_0430 = conj(IT_0220)*IT_0339 + IT_0220*IT_0343 +
       IT_0076*(IT_0138 + IT_0146 + IT_0319 + IT_0327 + IT_0344 + IT_0345 +
       IT_0346 + IT_0347 + IT_0348 + IT_0349 + IT_0350 + IT_0351 + IT_0352 +
       IT_0353) + IT_0075*(IT_0140 + IT_0148 + IT_0321 + IT_0329 + IT_0354 +
       IT_0355 + IT_0356 + IT_0357 + IT_0358 + IT_0359 + IT_0360 + IT_0361 +
       IT_0362 + IT_0363) + conj(IT_0038)*(IT_0142 + IT_0150 + IT_0315 + IT_0323
       + IT_0364 + IT_0365 + IT_0366 + IT_0367 + IT_0368 + IT_0369 + IT_0370 +
       IT_0371 + IT_0372 + IT_0373) + IT_0038*(IT_0144 + IT_0152 + IT_0317 +
       IT_0325 + IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378 + IT_0379 +
       IT_0380 + IT_0381 + IT_0382 + IT_0383) + IT_0064*(IT_0142 + IT_0150 +
       IT_0315 + IT_0323 + IT_0384 + IT_0385 + IT_0386 + IT_0387 + IT_0388 +
       IT_0389 + IT_0390 + IT_0391) + conj(IT_0067)*(IT_0364 + IT_0365 + IT_0366
       + IT_0367 + IT_0368 + IT_0369 + IT_0370 + IT_0371 + IT_0372 + IT_0373 +
       IT_0384 + IT_0385 + IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390 +
       IT_0391 + IT_0392) + IT_0063*(IT_0144 + IT_0152 + IT_0317 + IT_0325 +
       IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397 + IT_0398 + IT_0399 +
       IT_0400) + IT_0067*(IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378 +
       IT_0379 + IT_0380 + IT_0381 + IT_0382 + IT_0383 + IT_0393 + IT_0394 +
       IT_0395 + IT_0396 + IT_0397 + IT_0398 + IT_0399 + IT_0400 + IT_0401) +
       IT_0070*(IT_0138 + IT_0146 + IT_0319 + IT_0327 + IT_0402 + IT_0403 +
       IT_0404 + IT_0405 + IT_0406 + IT_0407 + IT_0408 + IT_0409) + conj(IT_0042
      )*(IT_0344 + IT_0345 + IT_0346 + IT_0347 + IT_0348 + IT_0349 + IT_0350 +
       IT_0351 + IT_0352 + IT_0353 + IT_0402 + IT_0403 + IT_0404 + IT_0405 +
       IT_0406 + IT_0407 + IT_0408 + IT_0409) + IT_0068*(IT_0140 + IT_0148 +
       IT_0321 + IT_0329 + IT_0410 + IT_0411 + IT_0412 + IT_0413 + IT_0414 +
       IT_0415 + IT_0416 + IT_0417) + IT_0042*(IT_0354 + IT_0355 + IT_0356 +
       IT_0357 + IT_0358 + IT_0359 + IT_0360 + IT_0361 + IT_0362 + IT_0363 +
       IT_0410 + IT_0411 + IT_0412 + IT_0413 + IT_0414 + IT_0415 + IT_0416 +
       IT_0417) + IT_0418*(IT_0332 + IT_0338 + IT_0419 + IT_0420 + IT_0421 +
       IT_0422 + IT_0423) + IT_0424*(IT_0334 + IT_0342 + IT_0425 + IT_0426 +
       IT_0427 + IT_0428 + IT_0429);
    const ccomplex_t IT_0431 = (-18)*IT_0042;
    const ccomplex_t IT_0432 = (-18)*IT_0060;
    const ccomplex_t IT_0433 = (-18)*IT_0061;
    const ccomplex_t IT_0434 = (-18)*IT_0058;
    const ccomplex_t IT_0435 = (-6)*IT_0062;
    const ccomplex_t IT_0436 = (-18)*conj(IT_0042);
    const ccomplex_t IT_0437 = (-18)*conj(IT_0060);
    const ccomplex_t IT_0438 = (-18)*conj(IT_0061);
    const ccomplex_t IT_0439 = (-18)*conj(IT_0058);
    const ccomplex_t IT_0440 = (-6)*conj(IT_0062);
    const ccomplex_t IT_0441 = conj(IT_0037)*IT_0339 + IT_0037*IT_0343 + conj
      (IT_0038)*(IT_0332 + IT_0338 + IT_0419 + IT_0422 + IT_0423) + conj(IT_0067
      )*(IT_0336 + IT_0337 + 2*IT_0338 + IT_0419 + IT_0420 + IT_0421 + IT_0422 +
       IT_0423) + IT_0038*(IT_0334 + IT_0342 + IT_0425 + IT_0428 + IT_0429) +
       IT_0067*(IT_0340 + IT_0341 + 2*IT_0342 + IT_0425 + IT_0426 + IT_0427 +
       IT_0428 + IT_0429) + IT_0059*(IT_0431 + IT_0432 + IT_0433 + IT_0434 +
       IT_0435) + (IT_0042 + IT_0049 + IT_0058)*(IT_0436 + IT_0437 + IT_0438 +
       IT_0439 + IT_0440);
    const ccomplex_t IT_0442 = IT_0005*IT_0042*conj(IT_0067);
    const ccomplex_t IT_0443 = IT_0005*conj(IT_0042)*IT_0067;
    const ccomplex_t IT_0444 = IT_0067*conj(IT_0067)*IT_0074;
    const ccomplex_t IT_0445 = 6*IT_0042;
    const ccomplex_t IT_0446 = 6*conj(IT_0042);
    const ccomplex_t IT_0447 = 6*IT_0049;
    const ccomplex_t IT_0448 = 6*IT_0058;
    const ccomplex_t IT_0449 = 18*IT_0038;
    const ccomplex_t IT_0450 = 6*conj(IT_0049);
    const ccomplex_t IT_0451 = 6*conj(IT_0058);
    const ccomplex_t IT_0452 = 18*conj(IT_0038);
    const ccomplex_t IT_0453 = 6*IT_0258;
    const ccomplex_t IT_0454 = 6*IT_0260;
    const ccomplex_t IT_0455 = 18*IT_0219;
    const ccomplex_t IT_0456 = (-18)*IT_0168;
    const ccomplex_t IT_0457 = (-18)*IT_0191;
    const ccomplex_t IT_0458 = (-18)*IT_0206;
    const ccomplex_t IT_0459 = (-6)*IT_0166;
    const ccomplex_t IT_0460 = (-6)*IT_0215;
    const ccomplex_t IT_0461 = 6*IT_0219;
    const ccomplex_t IT_0462 = 18*IT_0258;
    const ccomplex_t IT_0463 = 18*IT_0260;
    const ccomplex_t IT_0464 = (-18)*IT_0166;
    const ccomplex_t IT_0465 = (-18)*IT_0215;
    const ccomplex_t IT_0466 = (-6)*IT_0168;
    const ccomplex_t IT_0467 = (-6)*IT_0191;
    const ccomplex_t IT_0468 = (-6)*IT_0206;
    const ccomplex_t IT_0469 = 0.166666666666667*IT_0103 + 0.166666666666667
      *IT_0113 + 0.166666666666667*IT_0117 + 0.166666666666667*IT_0125 + IT_0220
       + 0.166666666666667*IT_0354 + 0.166666666666667*IT_0360 +
       0.166666666666667*IT_0461 + 0.166666666666667*IT_0462 + 0.166666666666667
      *IT_0463 + 0.166666666666667*IT_0464 + 0.166666666666667*IT_0465 +
       0.166666666666667*IT_0466 + 0.166666666666667*IT_0467 + 0.166666666666667
      *IT_0468;
    const ccomplex_t IT_0470 = 6*conj(IT_0219);
    const ccomplex_t IT_0471 = 18*conj(IT_0258);
    const ccomplex_t IT_0472 = 18*conj(IT_0260);
    const ccomplex_t IT_0473 = (-18)*conj(IT_0166);
    const ccomplex_t IT_0474 = (-18)*conj(IT_0215);
    const ccomplex_t IT_0475 = (-6)*conj(IT_0168);
    const ccomplex_t IT_0476 = (-6)*conj(IT_0191);
    const ccomplex_t IT_0477 = (-6)*conj(IT_0206);
    const ccomplex_t IT_0478 = 0.166666666666667*IT_0105 + 0.166666666666667
      *IT_0115 + 0.166666666666667*IT_0119 + 0.166666666666667*IT_0127 + conj
      (IT_0220) + 0.166666666666667*IT_0347 + 0.166666666666667*IT_0352 +
       0.166666666666667*IT_0470 + 0.166666666666667*IT_0471 + 0.166666666666667
      *IT_0472 + 0.166666666666667*IT_0473 + 0.166666666666667*IT_0474 +
       0.166666666666667*IT_0475 + 0.166666666666667*IT_0476 + 0.166666666666667
      *IT_0477;
    const ccomplex_t IT_0479 = 6*IT_0289;
    const ccomplex_t IT_0480 = 6*IT_0298;
    const ccomplex_t IT_0481 = 18*IT_0271;
    const ccomplex_t IT_0482 = 18*IT_0272;
    const ccomplex_t IT_0483 = (-18)*IT_0201;
    const ccomplex_t IT_0484 = (-18)*IT_0233;
    const ccomplex_t IT_0485 = (-18)*IT_0238;
    const ccomplex_t IT_0486 = (-6)*IT_0197;
    const ccomplex_t IT_0487 = 6*conj(IT_0289);
    const ccomplex_t IT_0488 = 6*conj(IT_0298);
    const ccomplex_t IT_0489 = 18*conj(IT_0271);
    const ccomplex_t IT_0490 = 18*conj(IT_0272);
    const ccomplex_t IT_0491 = (-18)*conj(IT_0201);
    const ccomplex_t IT_0492 = (-18)*conj(IT_0233);
    const ccomplex_t IT_0493 = (-18)*conj(IT_0238);
    const ccomplex_t IT_0494 = (-6)*conj(IT_0197);
    const ccomplex_t IT_0495 = 6*IT_0271;
    const ccomplex_t IT_0496 = 6*IT_0272;
    const ccomplex_t IT_0497 = 18*IT_0289;
    const ccomplex_t IT_0498 = 18*IT_0298;
    const ccomplex_t IT_0499 = (-18)*IT_0197;
    const ccomplex_t IT_0500 = (-6)*IT_0201;
    const ccomplex_t IT_0501 = (-6)*IT_0233;
    const ccomplex_t IT_0502 = (-6)*IT_0238;
    const ccomplex_t IT_0503 = 18*IT_0220 + 18*IT_0300 + 18*IT_0302 + 18
      *IT_0303 + IT_0376 + IT_0381 + IT_0453 + IT_0454 + IT_0455 + IT_0456 +
       IT_0457 + IT_0458 + IT_0459 + IT_0460 + IT_0495 + IT_0496 + IT_0497 +
       IT_0498 + IT_0499 + IT_0500 + IT_0501 + IT_0502;
    const ccomplex_t IT_0504 = 6*conj(IT_0271);
    const ccomplex_t IT_0505 = 6*conj(IT_0272);
    const ccomplex_t IT_0506 = 18*conj(IT_0289);
    const ccomplex_t IT_0507 = 18*conj(IT_0298);
    const ccomplex_t IT_0508 = (-18)*conj(IT_0197);
    const ccomplex_t IT_0509 = (-6)*conj(IT_0201);
    const ccomplex_t IT_0510 = (-6)*conj(IT_0233);
    const ccomplex_t IT_0511 = (-6)*conj(IT_0238);
    const ccomplex_t IT_0512 = 6*conj(IT_0258);
    const ccomplex_t IT_0513 = 6*conj(IT_0260);
    const ccomplex_t IT_0514 = 18*conj(IT_0219);
    const ccomplex_t IT_0515 = (-18)*conj(IT_0168);
    const ccomplex_t IT_0516 = (-18)*conj(IT_0191);
    const ccomplex_t IT_0517 = (-18)*conj(IT_0206);
    const ccomplex_t IT_0518 = (-6)*conj(IT_0166);
    const ccomplex_t IT_0519 = (-6)*conj(IT_0215);
    const ccomplex_t IT_0520 = 18*conj(IT_0220);
    const ccomplex_t IT_0521 = 18*conj(IT_0300) + 18*conj(IT_0302) + 18*conj
      (IT_0303) + IT_0364 + IT_0370 + IT_0504 + IT_0505 + IT_0506 + IT_0507 +
       IT_0508 + IT_0509 + IT_0510 + IT_0511 + IT_0512 + IT_0513 + IT_0514 +
       IT_0515 + IT_0516 + IT_0517 + IT_0518 + IT_0519 + IT_0520;
    const ccomplex_t IT_0522 = (conj(IT_0220) + conj(IT_0300) + conj(IT_0302) 
      + conj(IT_0303))*IT_0445 + (IT_0220 + IT_0300 + IT_0302 + IT_0303)*IT_0446
       + IT_0418*(IT_0447 + IT_0448 + IT_0449) + IT_0424*(IT_0450 + IT_0451 +
       IT_0452) + IT_0064*(IT_0099 + IT_0107 + IT_0121 + IT_0129 + 18*IT_0220 +
       IT_0376 + IT_0381 + IT_0453 + IT_0454 + IT_0455 + IT_0456 + IT_0457 +
       IT_0458 + IT_0459 + IT_0460) + 6*IT_0070*IT_0469 + 6*IT_0068*IT_0478 +
       IT_0076*(IT_0103 + IT_0113 + IT_0117 + IT_0125 + IT_0479 + IT_0480 +
       IT_0481 + IT_0482 + IT_0483 + IT_0484 + IT_0485 + IT_0486) + conj(IT_0042
      )*(IT_0354 + IT_0360 + IT_0461 + IT_0462 + IT_0463 + IT_0464 + IT_0465 +
       IT_0466 + IT_0467 + IT_0468 + IT_0479 + IT_0480 + IT_0481 + IT_0482 +
       IT_0483 + IT_0484 + IT_0485 + IT_0486) + IT_0075*(IT_0105 + IT_0115 +
       IT_0119 + IT_0127 + IT_0487 + IT_0488 + IT_0489 + IT_0490 + IT_0491 +
       IT_0492 + IT_0493 + IT_0494) + IT_0042*(IT_0347 + IT_0352 + IT_0470 +
       IT_0471 + IT_0472 + IT_0473 + IT_0474 + IT_0475 + IT_0476 + IT_0477 +
       IT_0487 + IT_0488 + IT_0489 + IT_0490 + IT_0491 + IT_0492 + IT_0493 +
       IT_0494) + conj(IT_0038)*(IT_0099 + IT_0107 + IT_0121 + IT_0129 + IT_0495
       + IT_0496 + IT_0497 + IT_0498 + IT_0499 + IT_0500 + IT_0501 + IT_0502) +
       conj(IT_0067)*IT_0503 + IT_0038*(IT_0101 + IT_0109 + IT_0123 + IT_0131 +
       IT_0504 + IT_0505 + IT_0506 + IT_0507 + IT_0508 + IT_0509 + IT_0510 +
       IT_0511) + IT_0063*(IT_0101 + IT_0109 + IT_0123 + IT_0131 + IT_0364 +
       IT_0370 + IT_0512 + IT_0513 + IT_0514 + IT_0515 + IT_0516 + IT_0517 +
       IT_0518 + IT_0519 + IT_0520) + IT_0067*IT_0521;
    const ccomplex_t IT_0523 = (-36)*conj(IT_0081);
    const ccomplex_t IT_0524 = (-12)*conj(IT_0079);
    const ccomplex_t IT_0525 = (-36)*conj(IT_0084) + IT_0354 + IT_0523 +
       IT_0524;
    const ccomplex_t IT_0526 = IT_0099 + IT_0107 + IT_0121;
    const ccomplex_t IT_0527 = IT_0101 + IT_0109 + IT_0123;
    const ccomplex_t IT_0528 = IT_0105 + IT_0119;
    const ccomplex_t IT_0529 = (-36)*IT_0079;
    const ccomplex_t IT_0530 = (-12)*conj(IT_0084);
    const ccomplex_t IT_0531 = (-12)*conj(IT_0081);
    const ccomplex_t IT_0532 = (IT_0103 + IT_0117 + IT_0125)*conj(IT_0166) + 
      (IT_0105 + IT_0115 + IT_0119 + IT_0127)*IT_0215 + (IT_0103 + IT_0113 +
       IT_0117 + IT_0125)*conj(IT_0215) + conj(IT_0219)*(IT_0142 + IT_0150 +
       IT_0315 + IT_0323) + IT_0219*(IT_0144 + IT_0152 + IT_0317 + IT_0325) +
       IT_0261*(IT_0138 + IT_0146 + IT_0319 + IT_0327) + (IT_0258 + IT_0260)*
      (IT_0140 + IT_0148 + IT_0321 + IT_0329) + conj(IT_0084)*IT_0352 + conj
      (IT_0081)*(IT_0347 + IT_0352) + IT_0273*(IT_0138 + IT_0146 + IT_0319 +
       IT_0327 + IT_0344 + IT_0345 + IT_0346 + IT_0347 + IT_0348 + IT_0349 +
       IT_0350 + IT_0351 + IT_0352 + IT_0353) + (IT_0271 + IT_0272)*(IT_0140 +
       IT_0148 + IT_0321 + IT_0329 + IT_0354 + IT_0355 + IT_0356 + IT_0357 +
       IT_0358 + IT_0359 + IT_0360 + IT_0361 + IT_0362 + IT_0363) + (conj
      (IT_0289) + conj(IT_0298))*(IT_0142 + IT_0150 + IT_0315 + IT_0323 +
       IT_0364 + IT_0365 + IT_0366 + IT_0367 + IT_0368 + IT_0369 + IT_0370 +
       IT_0371 + IT_0372 + IT_0373) + (IT_0289 + IT_0298)*(IT_0144 + IT_0152 +
       IT_0317 + IT_0325 + IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378 +
       IT_0379 + IT_0380 + IT_0381 + IT_0382 + IT_0383) + conj(IT_0220)*(IT_0142
       + IT_0150 + IT_0315 + IT_0323 + IT_0384 + IT_0385 + IT_0386 + IT_0387 +
       IT_0388 + IT_0389 + IT_0390 + IT_0391) + IT_0220*(IT_0144 + IT_0152 +
       IT_0317 + IT_0325 + IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397 +
       IT_0398 + IT_0399 + IT_0400) + (IT_0142 + IT_0150 + IT_0315 + IT_0323 +
       IT_0364 + IT_0365 + IT_0366 + IT_0367 + IT_0368 + IT_0369 + IT_0370 +
       IT_0371 + IT_0372 + IT_0373 + IT_0392)*IT_0418 + (IT_0144 + IT_0152 +
       IT_0317 + IT_0325 + IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378 +
       IT_0379 + IT_0380 + IT_0381 + IT_0382 + IT_0383 + IT_0401)*IT_0424 + conj
      (IT_0197)*(IT_0099 + IT_0107 + IT_0121 + IT_0129 + IT_0376 + IT_0381 +
       IT_0453 + IT_0454 + IT_0455 + IT_0456 + IT_0457 + IT_0458 + IT_0459 +
       IT_0460) + (conj(IT_0201) + conj(IT_0233) + conj(IT_0238))*(IT_0103 +
       IT_0113 + IT_0117 + IT_0125 + IT_0354 + IT_0360 + IT_0461 + IT_0462 +
       IT_0463 + IT_0464 + IT_0465 + IT_0466 + IT_0467 + IT_0468) + (IT_0201 +
       IT_0233 + IT_0238)*(IT_0105 + IT_0115 + IT_0119 + IT_0127 + IT_0347 +
       IT_0352 + IT_0470 + IT_0471 + IT_0472 + IT_0473 + IT_0474 + IT_0475 +
       IT_0476 + IT_0477) + IT_0197*(IT_0101 + IT_0109 + IT_0123 + IT_0131 +
       IT_0364 + IT_0370 + IT_0512 + IT_0513 + IT_0514 + IT_0515 + IT_0516 +
       IT_0517 + IT_0518 + IT_0519) + IT_0084*IT_0525 + IT_0081*(IT_0360 +
       IT_0525) + (conj(IT_0168) + conj(IT_0191))*IT_0526 + conj(IT_0206)*
      (IT_0129 + IT_0526) + (IT_0168 + IT_0191)*IT_0527 + IT_0206*(IT_0131 +
       IT_0527) + IT_0166*(IT_0127 + IT_0528) + IT_0158*(IT_0138 + IT_0319 +
       IT_0528) + conj(IT_0079)*(IT_0370 + IT_0529) + conj(IT_0112)*((-12)
      *IT_0081 + (-12)*IT_0084 + IT_0364 + IT_0365 + IT_0367 + IT_0368 + IT_0370
       + IT_0529) + IT_0079*(IT_0376 + IT_0530 + IT_0531) + IT_0112*((-36)*conj
      (IT_0079) + (-36)*conj(IT_0112) + IT_0374 + IT_0376 + IT_0377 + IT_0378 +
       IT_0381 + IT_0530 + IT_0531);
    const ccomplex_t IT_0533 = (IT_0081 + IT_0084)*IT_0095;
    const ccomplex_t IT_0534 = IT_0169*conj(IT_0219);
    const ccomplex_t IT_0535 = IT_0134*conj(IT_0166);
    const ccomplex_t IT_0536 = IT_0134*conj(IT_0215);
    const ccomplex_t IT_0537 = IT_0154*conj(IT_0168);
    const ccomplex_t IT_0538 = IT_0154*conj(IT_0191);
    const ccomplex_t IT_0539 = IT_0154*conj(IT_0206);
    const ccomplex_t IT_0540 = IT_0169*conj(IT_0220);
    const ccomplex_t IT_0541 = IT_0005*IT_0071 + 2*IT_0074*IT_0077 + conj
      (IT_0079)*IT_0097 + 2*IT_0100 + 2*IT_0102 + 2*IT_0104 + 2*IT_0106 + 2
      *IT_0108 + 2*IT_0110 + 2*IT_0114 + 2*IT_0116 + 2*IT_0118 + 2*IT_0120 + 2
      *IT_0122 + 2*IT_0124 + 2*IT_0126 + 2*IT_0128 + 2*IT_0130 + 2*IT_0132 +
       IT_0133*IT_0134 + IT_0079*(IT_0135 + IT_0136) + 2*IT_0139 + 2*IT_0141 + 2
      *IT_0143 + 2*IT_0145 + 2*IT_0147 + 2*IT_0149 + 2*IT_0151 + 2*IT_0153 +
       IT_0084*IT_0155 + IT_0081*(IT_0155 + IT_0159) + conj(IT_0166)*(IT_0155 +
       IT_0159 + IT_0170 + IT_0171 + IT_0172 + IT_0173 + IT_0174) + IT_0158*
      (IT_0082*IT_0092 + IT_0134*IT_0158 + IT_0175 + IT_0176 + IT_0177 + IT_0178
      ) + IT_0166*(IT_0177 + IT_0178 + IT_0179 + IT_0180 + IT_0181 + IT_0182) +
       conj(IT_0168)*(IT_0136 + IT_0183 + IT_0184 + IT_0185 + IT_0186 + IT_0187)
       + conj(IT_0191)*(IT_0136 + IT_0183 + IT_0184 + IT_0185 + IT_0186 +
       IT_0187 + IT_0192 + IT_0193) + conj(IT_0197)*(IT_0185 + IT_0186 + IT_0187
       + IT_0093*IT_0197 + IT_0169*IT_0201 + IT_0202) + (IT_0136 + IT_0183 +
       IT_0184 + IT_0185 + IT_0186 + IT_0187 + IT_0192 + IT_0193 + IT_0202)*conj
      (IT_0206) + IT_0197*(IT_0207 + IT_0208 + IT_0209 + IT_0210) + IT_0168*
      (IT_0094 + IT_0207 + IT_0209 + IT_0210 + IT_0211) + IT_0191*(IT_0094 +
       IT_0207 + IT_0209 + IT_0210 + IT_0211 + IT_0212 + IT_0213) + IT_0206*
      (IT_0094 + IT_0207 + IT_0208 + IT_0209 + IT_0210 + IT_0211 + IT_0212 +
       IT_0213 + IT_0169*conj(IT_0215) + IT_0093*IT_0216 + IT_0134*IT_0221) +
       conj(IT_0112)*(IT_0097 + IT_0211 + IT_0222 + IT_0223 + IT_0224) + IT_0112
      *(conj(IT_0081)*IT_0082 + IT_0095*(conj(IT_0079) + conj(IT_0112)) +
       IT_0135 + IT_0136 + IT_0183 + IT_0225 + IT_0226 + IT_0227) + conj(IT_0215
      )*(IT_0155 + IT_0159 + IT_0170 + IT_0171 + IT_0172 + IT_0173 + IT_0174 +
       IT_0169*IT_0191 + IT_0229) + conj(IT_0201)*(IT_0172 + IT_0173 + IT_0174 +
       IT_0169*IT_0197 + IT_0228 + IT_0229) + IT_0239*(IT_0172 + IT_0173 +
       IT_0174 + IT_0228 + IT_0229 + IT_0240 + IT_0241) + IT_0201*(IT_0177 +
       IT_0178 + IT_0182 + IT_0242) + IT_0215*(IT_0177 + IT_0178 + IT_0179 +
       IT_0180 + IT_0181 + IT_0182 + IT_0093*(conj(IT_0166) + conj(IT_0215)) +
       IT_0169*IT_0216 + IT_0242) + IT_0233*(IT_0177 + IT_0178 + IT_0182 +
       IT_0242 + IT_0243 + IT_0244) + IT_0238*(IT_0177 + IT_0178 + IT_0182 +
       IT_0093*IT_0239 + IT_0242 + IT_0243 + IT_0244) + IT_0221*(IT_0094 +
       IT_0211 + IT_0222 + IT_0223 + IT_0224 + IT_0245 + IT_0246 + IT_0247 +
       IT_0248 + IT_0249 + IT_0250) + IT_0219*(IT_0136 + IT_0183 + IT_0134*conj
      (IT_0206) + IT_0225 + IT_0226 + IT_0227 + IT_0251 + IT_0252 + IT_0253 +
       IT_0254 + IT_0255) + IT_0220*(IT_0136 + IT_0183 + IT_0134*conj(IT_0206) +
       IT_0093*IT_0221 + IT_0225 + IT_0226 + IT_0227 + IT_0251 + IT_0252 +
       IT_0253 + IT_0254 + IT_0255 + IT_0169*IT_0261) + IT_0261*(IT_0175 +
       IT_0176 + IT_0179 + IT_0181 + IT_0262 + IT_0263 + IT_0264 + IT_0265 +
       IT_0266 + IT_0267 + IT_0268 + IT_0269 + IT_0270) + IT_0273*(IT_0175 +
       IT_0176 + IT_0262 + IT_0264 + IT_0274 + IT_0275 + IT_0276 + IT_0277 +
       IT_0278) + IT_0271*(IT_0281 + IT_0282 + IT_0283 + IT_0284 + IT_0285 +
       IT_0286 + IT_0287) + IT_0272*(IT_0093*IT_0273 + IT_0281 + IT_0282 +
       IT_0283 + IT_0284 + IT_0285 + IT_0286 + IT_0287) + conj(IT_0289)*(IT_0245
       + IT_0246 + IT_0247 + IT_0248 + IT_0290 + IT_0291 + IT_0292 + IT_0293 +
       IT_0294 + IT_0295 + IT_0296) + conj(IT_0298)*(IT_0245 + IT_0246 + IT_0247
       + IT_0248 + IT_0290 + IT_0291 + IT_0292 + IT_0293 + IT_0294 + IT_0295 +
       IT_0296 + IT_0299) + conj(IT_0300)*(IT_0245 + IT_0246 + IT_0247 + IT_0248
       + IT_0290 + IT_0291 + IT_0292 + IT_0293 + IT_0294 + IT_0295 + IT_0296 +
       IT_0299 + IT_0301) + IT_0304*(IT_0245 + IT_0246 + IT_0247 + IT_0248 +
       IT_0290 + IT_0291 + IT_0292 + IT_0293 + IT_0294 + IT_0295 + IT_0296 +
       IT_0299 + IT_0301 + IT_0305) + IT_0289*(IT_0251 + IT_0252 + IT_0253 +
       IT_0254 + IT_0306 + IT_0307 + IT_0308 + IT_0309 + IT_0310 + IT_0311) +
       IT_0298*(IT_0251 + IT_0252 + IT_0253 + IT_0254 + IT_0306 + IT_0307 +
       IT_0308 + IT_0309 + IT_0310 + IT_0311 + IT_0312) + IT_0300*(IT_0251 +
       IT_0252 + IT_0253 + IT_0254 + IT_0306 + IT_0307 + IT_0308 + IT_0309 +
       IT_0310 + IT_0311 + IT_0312 + IT_0313) + IT_0302*(IT_0251 + IT_0252 +
       IT_0253 + IT_0254 + IT_0306 + IT_0307 + IT_0308 + IT_0309 + IT_0310 +
       IT_0311 + IT_0312 + IT_0313 + IT_0314) + IT_0303*(IT_0251 + IT_0252 +
       IT_0253 + IT_0254 + IT_0093*IT_0304 + IT_0306 + IT_0307 + IT_0308 +
       IT_0309 + IT_0310 + IT_0311 + IT_0312 + IT_0313 + IT_0314) + 2*IT_0316 +
       2*IT_0318 + 2*IT_0320 + 2*IT_0322 + 2*IT_0324 + 2*IT_0326 + 2*IT_0328 + 2
      *IT_0330 + 2*IT_0333 + 2*IT_0335 + IT_0137*IT_0430 + IT_0331*IT_0441 + 2
      *IT_0442 + 2*IT_0443 + 2*IT_0444 + IT_0098*IT_0522 + IT_0006*IT_0532 +
       conj(IT_0084)*(IT_0179 + IT_0533) + conj(IT_0081)*(IT_0079*IT_0082 +
       IT_0179 + IT_0181 + IT_0533) + IT_0258*(IT_0155 + IT_0159 + IT_0279 +
       IT_0280 + IT_0282 + IT_0283 + IT_0534 + IT_0535 + IT_0536 + IT_0537 +
       IT_0538 + IT_0539 + IT_0540) + IT_0260*(IT_0155 + IT_0159 + IT_0093
      *IT_0261 + IT_0279 + IT_0280 + IT_0282 + IT_0283 + IT_0534 + IT_0535 +
       IT_0536 + IT_0537 + IT_0538 + IT_0539 + IT_0540);
    return create_ccomplex_return(IT_0541);
}

