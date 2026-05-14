#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_smu_L_to_anti_mu_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_smu_L_to_anti_mu_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
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
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0008 + 0.5*IT_0011);
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = pow(m_mu, 2);
    const ccomplex_t IT_0016 = pow(m_smu_L, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0025 = IT_0004*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0036 = IT_0004*IT_0035;
    const ccomplex_t IT_0037 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0038 = IT_0004*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 0.5*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = IT_0034*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0049 = IT_0004*IT_0048;
    const ccomplex_t IT_0050 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + 0.5*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0019 + IT_0032 + IT_0045 + IT_0058;
    const ccomplex_t IT_0060 = m_C_1*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0062 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0063 = IT_0004*IT_0062;
    const ccomplex_t IT_0064 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0065 = IT_0004*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + (-0.5)*IT_0066);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0006*IT_0068;
    const ccomplex_t IT_0070 = IT_0017*IT_0061*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0072 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0073 = IT_0004*IT_0072;
    const ccomplex_t IT_0074 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0075 = IT_0004*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + (-0.5)*IT_0076);
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = IT_0021*IT_0078;
    const ccomplex_t IT_0080 = IT_0030*IT_0071*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0082 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0083 = IT_0004*IT_0082;
    const ccomplex_t IT_0084 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0085 = IT_0004*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + (-0.5)*IT_0086);
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0034*IT_0088;
    const ccomplex_t IT_0090 = IT_0043*IT_0081*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0092 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0093 = IT_0004*IT_0092;
    const ccomplex_t IT_0094 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0095 = IT_0004*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + (-0.5)*IT_0096);
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = IT_0047*IT_0098;
    const ccomplex_t IT_0100 = IT_0056*IT_0091*IT_0099;
    const ccomplex_t IT_0101 = -IT_0070 + -IT_0080 + -IT_0090 + -IT_0100;
    const ccomplex_t IT_0102 = pow(m_W, -2);
    const ccomplex_t IT_0103 = pow(m_W, 2);
    const ccomplex_t IT_0104 = s_13*IT_0103;
    const ccomplex_t IT_0105 = s_14*s_34;
    const ccomplex_t IT_0106 = (-2)*IT_0105;
    const ccomplex_t IT_0107 = IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = IT_0102*IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = s_13 + IT_0109;
    const ccomplex_t IT_0111 = cos(theta_W);
    const ccomplex_t IT_0112 = cpow(IT_0111, -1);
    const ccomplex_t IT_0113 = N_B1*e_em;
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = N_W1*e_em;
    const ccomplex_t IT_0117 = IT_0004*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + IT_0118);
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0068*IT_0120;
    const ccomplex_t IT_0122 = IT_0017*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = N_B2*e_em;
    const ccomplex_t IT_0125 = IT_0112*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = N_W2*e_em;
    const ccomplex_t IT_0128 = IT_0004*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + IT_0129);
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = IT_0078*IT_0131;
    const ccomplex_t IT_0133 = IT_0030*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = N_B3*e_em;
    const ccomplex_t IT_0136 = IT_0112*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = N_W3*e_em;
    const ccomplex_t IT_0139 = IT_0004*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0137 + IT_0140);
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = IT_0088*IT_0142;
    const ccomplex_t IT_0144 = IT_0043*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = N_B4*e_em;
    const ccomplex_t IT_0147 = IT_0112*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = N_W4*e_em;
    const ccomplex_t IT_0150 = IT_0004*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0148 + IT_0151);
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = IT_0098*IT_0153;
    const ccomplex_t IT_0155 = IT_0056*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = IT_0123 + IT_0134 + IT_0145 + IT_0156;
    const ccomplex_t IT_0158 = m_C_1*IT_0157;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0013*IT_0120;
    const ccomplex_t IT_0161 = IT_0017*IT_0061*IT_0160;
    const ccomplex_t IT_0162 = IT_0028*IT_0131;
    const ccomplex_t IT_0163 = IT_0030*IT_0071*IT_0162;
    const ccomplex_t IT_0164 = IT_0041*IT_0142;
    const ccomplex_t IT_0165 = IT_0043*IT_0081*IT_0164;
    const ccomplex_t IT_0166 = IT_0054*IT_0153;
    const ccomplex_t IT_0167 = IT_0056*IT_0091*IT_0166;
    const ccomplex_t IT_0168 = IT_0161 + IT_0163 + IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = m_mu*m_C_1;
    const ccomplex_t IT_0170 = 2*IT_0123 + 2*IT_0134 + 2*IT_0145 + 2*IT_0156;
    const ccomplex_t IT_0171 = s_14*s_34*IT_0102;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = s_13 + IT_0172;
    const ccomplex_t IT_0174 = m_C_1*IT_0173;
    const ccomplex_t IT_0175 = (-2)*IT_0174;
    const ccomplex_t IT_0176 = -IT_0019 + -IT_0032 + -IT_0045 + -IT_0058;
    const ccomplex_t IT_0177 = s_14*IT_0103;
    const ccomplex_t IT_0178 = IT_0102*IT_0177;
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = s_14 + IT_0179;
    const ccomplex_t IT_0181 = m_mu*IT_0180;
    const ccomplex_t IT_0182 = 2*IT_0181;
    const ccomplex_t IT_0183 = (-2)*IT_0182;
    const ccomplex_t IT_0184 = (-0.166666666666667)*IT_0176*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004;
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0188 = (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = pow(m_C_1, 2);
    const ccomplex_t IT_0190 = cpow((-2)*s_13 + IT_0015 + IT_0189 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0191 = IT_0188*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = IT_0186*IT_0192;
    const ccomplex_t IT_0194 = 2*IT_0019 + 2*IT_0032 + 2*IT_0045 + 2*IT_0058 +
       IT_0193;
    const ccomplex_t IT_0195 = pow(s_14, 2);
    const ccomplex_t IT_0196 = IT_0102*IT_0195;
    const ccomplex_t IT_0197 = -IT_0189;
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = m_mu*IT_0198;
    const ccomplex_t IT_0200 = (-2)*IT_0199;
    const ccomplex_t IT_0201 = IT_0194*IT_0200;
    const ccomplex_t IT_0202 = (-0.166666666666667)*IT_0201;
    const ccomplex_t IT_0203 = (-6)*IT_0169;
    const ccomplex_t IT_0204 = (-0.166666666666667)*IT_0101*IT_0203;
    const ccomplex_t IT_0205 = IT_0159 + IT_0168;
    const ccomplex_t IT_0206 = 4*IT_0110;
    const ccomplex_t IT_0207 = (-0.166666666666667)*IT_0205*IT_0206;
    const ccomplex_t IT_0208 = -IT_0123 + -IT_0134 + -IT_0145 + -IT_0156;
    const ccomplex_t IT_0209 = s_34*IT_0103;
    const ccomplex_t IT_0210 = IT_0102*IT_0209;
    const ccomplex_t IT_0211 = (-0.25)*IT_0210;
    const ccomplex_t IT_0212 = s_34 + IT_0211;
    const ccomplex_t IT_0213 = m_C_1*IT_0212;
    const ccomplex_t IT_0214 = (-4)*IT_0213;
    const ccomplex_t IT_0215 = 2*IT_0214;
    const ccomplex_t IT_0216 = (-0.166666666666667)*IT_0208*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0004;
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = IT_0190*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = IT_0186*IT_0220;
    const ccomplex_t IT_0222 = pow(s_34, 2);
    const ccomplex_t IT_0223 = IT_0102*IT_0222;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0015 + IT_0224;
    const ccomplex_t IT_0226 = m_C_1*IT_0225;
    const ccomplex_t IT_0227 = 2*IT_0226;
    const ccomplex_t IT_0228 = (-0.166666666666667)*IT_0221*IT_0227;
    const ccomplex_t IT_0229 = -IT_0221;
    const ccomplex_t IT_0230 = 2*IT_0174;
    const ccomplex_t IT_0231 = (-0.166666666666667)*IT_0229*IT_0230;
    const ccomplex_t IT_0232 = -IT_0193;
    const ccomplex_t IT_0233 = m_mu*IT_0173;
    const ccomplex_t IT_0234 = 2*IT_0233;
    const ccomplex_t IT_0235 = (-0.166666666666667)*IT_0232*IT_0234;
    const ccomplex_t IT_0236 = IT_0060*IT_0169 + (-0.166666666666667)*IT_0170
      *IT_0175 + IT_0184 + IT_0202 + IT_0204 + IT_0207 + IT_0216 + IT_0228 +
       IT_0231 + IT_0235;
    const ccomplex_t IT_0237 = 2*conj(IT_0060);
    const ccomplex_t IT_0238 = (-2)*conj(IT_0060);
    const ccomplex_t IT_0239 = conj(IT_0232)*IT_0234;
    const ccomplex_t IT_0240 = conj(IT_0229)*IT_0230;
    const ccomplex_t IT_0241 = conj(IT_0221)*IT_0227;
    const ccomplex_t IT_0242 = (-6)*conj(IT_0060)*IT_0169;
    const ccomplex_t IT_0243 = conj(IT_0101)*IT_0203;
    const ccomplex_t IT_0244 = conj(IT_0194)*IT_0200;
    const ccomplex_t IT_0245 = conj(IT_0170)*IT_0175;
    const ccomplex_t IT_0246 = 2*IT_0060;
    const ccomplex_t IT_0247 = s_13*s_14;
    const ccomplex_t IT_0248 = s_34*IT_0189;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = IT_0247 + IT_0249;
    const ccomplex_t IT_0251 = 2*IT_0250;
    const ccomplex_t IT_0252 = 2*IT_0199;
    const ccomplex_t IT_0253 = s_13*IT_0198;
    const ccomplex_t IT_0254 = 2*IT_0253;
    const ccomplex_t IT_0255 = s_13*IT_0173;
    const ccomplex_t IT_0256 = 2*IT_0255;
    const ccomplex_t IT_0257 = IT_0169*IT_0173;
    const ccomplex_t IT_0258 = 2*IT_0257;
    const ccomplex_t IT_0259 = (-2)*IT_0253;
    const ccomplex_t IT_0260 = IT_0169*IT_0198;
    const ccomplex_t IT_0261 = (-2)*IT_0260;
    const ccomplex_t IT_0262 = 2*IT_0260;
    const ccomplex_t IT_0263 = (-2)*IT_0255;
    const ccomplex_t IT_0264 = (-2)*IT_0257;
    const ccomplex_t IT_0265 = s_14*IT_0015;
    const ccomplex_t IT_0266 = s_13*s_34;
    const ccomplex_t IT_0267 = -IT_0266;
    const ccomplex_t IT_0268 = IT_0265 + IT_0267;
    const ccomplex_t IT_0269 = (-2)*IT_0268;
    const ccomplex_t IT_0270 = s_13*IT_0225;
    const ccomplex_t IT_0271 = (-2)*IT_0270;
    const ccomplex_t IT_0272 = IT_0169*IT_0225;
    const ccomplex_t IT_0273 = (-2)*IT_0272;
    const ccomplex_t IT_0274 = 2*IT_0268;
    const ccomplex_t IT_0275 = 2*IT_0182;
    const ccomplex_t IT_0276 = (-2)*IT_0214;
    const ccomplex_t IT_0277 = IT_0060*IT_0110 + 0.25*IT_0101*IT_0206 + 0.25
      *IT_0200*IT_0229 + 0.25*IT_0194*IT_0230 + 0.25*IT_0227*IT_0232 + 0.25
      *IT_0221*IT_0234 + 0.25*IT_0170*IT_0252 + 0.25*IT_0208*IT_0275 + 0.25
      *IT_0176*IT_0276;
    const ccomplex_t IT_0278 = pow(m_W, 4);
    const ccomplex_t IT_0279 = s_13*IT_0278;
    const ccomplex_t IT_0280 = IT_0102*IT_0279;
    const ccomplex_t IT_0281 = (-0.5)*IT_0280;
    const ccomplex_t IT_0282 = IT_0104 + IT_0106 + IT_0281;
    const ccomplex_t IT_0283 = (-4)*IT_0282;
    const ccomplex_t IT_0284 = (-2)*IT_0250;
    const ccomplex_t IT_0285 = IT_0103*IT_0169;
    const ccomplex_t IT_0286 = (-3)*IT_0285;
    const ccomplex_t IT_0287 = (-2)*IT_0286;
    const ccomplex_t IT_0288 = (-0.5)*IT_0194;
    const ccomplex_t IT_0289 = (-0.5)*IT_0232;
    const ccomplex_t IT_0290 = (-0.5)*IT_0208;
    const ccomplex_t IT_0291 = (-0.5)*IT_0101;
    const ccomplex_t IT_0292 = (-2)*conj(IT_0229)*(IT_0060*IT_0199 + (-0.5)
      *IT_0229*IT_0254 + (-0.5)*IT_0170*IT_0259 + (-0.5)*IT_0221*IT_0263 +
       IT_0262*IT_0288 + IT_0264*IT_0289 + IT_0284*IT_0290 + IT_0200*IT_0291);
    const ccomplex_t IT_0293 = 4*conj(IT_0060)*(IT_0060 + IT_0101)*IT_0110 + (
      -6)*(conj(IT_0159) + conj(IT_0168))*IT_0236 + (IT_0174*IT_0194 + IT_0170
      *IT_0199 + IT_0182*IT_0208 + IT_0226*IT_0232 + IT_0221*IT_0233)*IT_0237 + 
      (IT_0176*IT_0214 + IT_0199*IT_0229)*IT_0238 + IT_0205*(IT_0239 + IT_0240 +
       IT_0241 + IT_0242 + IT_0243 + IT_0244 + IT_0245) + conj(IT_0170)*(IT_0199
      *IT_0246 + IT_0208*IT_0251 + IT_0101*IT_0252 + IT_0170*IT_0254 + IT_0221
      *IT_0256 + IT_0232*IT_0258 + IT_0229*IT_0259 + IT_0194*IT_0261) + conj
      (IT_0194)*(IT_0101*IT_0230 + IT_0174*IT_0246 + IT_0176*IT_0251 + IT_0194
      *IT_0254 + IT_0170*IT_0261 + IT_0229*IT_0262 + IT_0232*IT_0263 + IT_0221
      *IT_0264) + conj(IT_0221)*(IT_0101*IT_0234 + IT_0233*IT_0246 + IT_0170
      *IT_0256 + IT_0229*IT_0263 + IT_0194*IT_0264 + IT_0208*IT_0269 + IT_0221
      *IT_0271 + IT_0232*IT_0273) + conj(IT_0232)*(IT_0101*IT_0227 + IT_0226
      *IT_0246 + IT_0170*IT_0258 + IT_0194*IT_0263 + IT_0229*IT_0264 + IT_0232
      *IT_0271 + IT_0221*IT_0273 + IT_0176*IT_0274) + 4*conj(IT_0101)*IT_0277 +
       conj(IT_0208)*(IT_0205*IT_0215 + IT_0182*IT_0246 + IT_0170*IT_0251 +
       IT_0221*IT_0269 + IT_0101*IT_0275 + IT_0208*IT_0283 + IT_0229*IT_0284 +
       IT_0176*IT_0287) + (-2)*conj(IT_0176)*((-0.5)*IT_0183*IT_0205 + IT_0060
      *IT_0214 + (-0.5)*IT_0176*IT_0283 + IT_0251*IT_0288 + IT_0274*IT_0289 +
       IT_0287*IT_0290 + IT_0276*IT_0291) + IT_0292;
    return create_ccomplex_return(IT_0293);
}

