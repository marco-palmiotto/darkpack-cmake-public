#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_sd_L_to_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_sd_L_to_d_d(
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
    const creal_t m_sd_L = param->m_sd_L;
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
    const ccomplex_t IT_0000 = pow(m_sd_L, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = (-18)*IT_0001;
    const ccomplex_t IT_0003 = s_23*s_24;
    const ccomplex_t IT_0004 = 36*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = cos(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = conj(N_B4)*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = conj(N_W4)*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0010 + (-3)*IT_0015);
    const ccomplex_t IT_0017 = 0.166666666666667*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = cos(beta);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0012*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = pow(m_d, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = conj(N_B2)*e_em;
    const ccomplex_t IT_0030 = IT_0007*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = conj(N_W2)*e_em;
    const ccomplex_t IT_0033 = IT_0012*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + (-3)*IT_0034);
    const ccomplex_t IT_0036 = 0.166666666666667*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0012*IT_0018*IT_0020;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = conj(N_B3)*e_em;
    const ccomplex_t IT_0045 = IT_0007*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W3)*e_em;
    const ccomplex_t IT_0048 = IT_0012*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + (-3)*IT_0049);
    const ccomplex_t IT_0051 = 0.166666666666667*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0012*IT_0018*IT_0020;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0059 = conj(N_B1)*e_em;
    const ccomplex_t IT_0060 = IT_0007*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W1)*e_em;
    const ccomplex_t IT_0063 = IT_0012*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + (-3)*IT_0064);
    const ccomplex_t IT_0066 = 0.166666666666667*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0012*IT_0018*IT_0020;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = IT_0058*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0057 + -IT_0071;
    const ccomplex_t IT_0073 = IT_0028 + IT_0043 + IT_0072;
    const ccomplex_t IT_0074 = conj(IT_0028) + conj(IT_0043) + conj(IT_0072);
    const ccomplex_t IT_0075 = -IT_0071;
    const ccomplex_t IT_0076 = -IT_0027 + -IT_0042;
    const ccomplex_t IT_0077 = -IT_0057;
    const ccomplex_t IT_0078 = IT_0075 + IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = conj(IT_0075) + conj(IT_0076) + conj(IT_0077);
    const ccomplex_t IT_0080 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0069*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0039*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0054*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0023*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0089 + IT_0092;
    const ccomplex_t IT_0094 = IT_0082 + IT_0093;
    const ccomplex_t IT_0095 = conj(IT_0089) + conj(IT_0092);
    const ccomplex_t IT_0096 = conj(IT_0082) + IT_0095;
    const ccomplex_t IT_0097 = IT_0088 + IT_0092;
    const ccomplex_t IT_0098 = IT_0082 + IT_0085;
    const ccomplex_t IT_0099 = conj(IT_0097) + conj(IT_0098);
    const ccomplex_t IT_0100 = (-6)*IT_0001;
    const ccomplex_t IT_0101 = 12*IT_0003;
    const ccomplex_t IT_0102 = IT_0100 + IT_0101;
    const ccomplex_t IT_0103 = IT_0097 + IT_0098;
    const ccomplex_t IT_0104 = m_d*IT_0092;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = m_d*IT_0042;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = 6*s_34;
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = m_d*IT_0071;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = IT_0108*IT_0111;
    const ccomplex_t IT_0113 = IT_0109 + IT_0112;
    const ccomplex_t IT_0114 = m_d*IT_0057;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = IT_0108*IT_0115;
    const ccomplex_t IT_0117 = m_d*IT_0027;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = IT_0108*IT_0118;
    const ccomplex_t IT_0120 = m_d*IT_0085;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = 18*s_34;
    const ccomplex_t IT_0123 = IT_0121*IT_0122;
    const ccomplex_t IT_0124 = cpow(IT_0053, 2);
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*m_N_3*IT_0124;
    const ccomplex_t IT_0126 = IT_0086*IT_0125;
    const ccomplex_t IT_0127 = cpow(IT_0022, 2);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*m_N_4*IT_0127;
    const ccomplex_t IT_0129 = IT_0090*IT_0128;
    const ccomplex_t IT_0130 = IT_0126 + IT_0129;
    const ccomplex_t IT_0131 = IT_0122*IT_0130;
    const ccomplex_t IT_0132 = IT_0105*IT_0122;
    const ccomplex_t IT_0133 = IT_0113 + IT_0116 + IT_0119 + IT_0123 + IT_0131
       + IT_0132;
    const ccomplex_t IT_0134 = conj(IT_0107)*IT_0108;
    const ccomplex_t IT_0135 = IT_0108*conj(IT_0111);
    const ccomplex_t IT_0136 = IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = cpow(IT_0036, 2);
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*m_N_2*IT_0137;
    const ccomplex_t IT_0139 = IT_0040*IT_0138;
    const ccomplex_t IT_0140 = cpow(IT_0017, 2);
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_N_4*IT_0140;
    const ccomplex_t IT_0142 = IT_0025*IT_0141;
    const ccomplex_t IT_0143 = cpow(IT_0051, 2);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*m_N_3*IT_0143;
    const ccomplex_t IT_0145 = IT_0055*IT_0144;
    const ccomplex_t IT_0146 = IT_0139 + IT_0142 + IT_0145;
    const ccomplex_t IT_0147 = pow(m_sG, 2);
    const ccomplex_t IT_0148 = cpow((-2)*s_23 + IT_0000 + IT_0024 + -IT_0147 +
       reg_prop, -1);
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0150 = cpow(IT_0149, 2);
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*m_sG*IT_0150;
    const ccomplex_t IT_0152 = IT_0148*IT_0151;
    const ccomplex_t IT_0153 = cpow((-2)*s_13 + IT_0000 + IT_0024 + -IT_0147 +
       reg_prop, -1);
    const ccomplex_t IT_0154 = IT_0151*IT_0153;
    const ccomplex_t IT_0155 = (-0.166666666666667)*IT_0152 + 0.5*IT_0154;
    const ccomplex_t IT_0156 = IT_0108*IT_0155;
    const ccomplex_t IT_0157 = IT_0085 + IT_0088 + IT_0092;
    const ccomplex_t IT_0158 = m_d*IT_0157;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0108*IT_0159;
    const ccomplex_t IT_0161 = IT_0083*IT_0138;
    const ccomplex_t IT_0162 = IT_0090*IT_0141;
    const ccomplex_t IT_0163 = IT_0086*IT_0144;
    const ccomplex_t IT_0164 = IT_0161 + IT_0162 + IT_0163;
    const ccomplex_t IT_0165 = IT_0108*IT_0164;
    const ccomplex_t IT_0166 = IT_0118*IT_0122;
    const ccomplex_t IT_0167 = IT_0122*IT_0146;
    const ccomplex_t IT_0168 = IT_0156 + IT_0160 + IT_0165 + IT_0166 + IT_0167;
    const ccomplex_t IT_0169 = IT_0000*IT_0024;
    const ccomplex_t IT_0170 = (-18)*IT_0028;
    const ccomplex_t IT_0171 = (-18)*IT_0043;
    const ccomplex_t IT_0172 = (-18)*IT_0072;
    const ccomplex_t IT_0173 = (-6)*IT_0089;
    const ccomplex_t IT_0174 = IT_0170 + IT_0171 + IT_0172 + IT_0173;
    const ccomplex_t IT_0175 = (-18)*conj(IT_0028);
    const ccomplex_t IT_0176 = (-18)*conj(IT_0043);
    const ccomplex_t IT_0177 = (-18)*conj(IT_0072);
    const ccomplex_t IT_0178 = (-6)*conj(IT_0089);
    const ccomplex_t IT_0179 = IT_0175 + IT_0176 + IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = (-18)*IT_0089;
    const ccomplex_t IT_0181 = (-6)*IT_0028;
    const ccomplex_t IT_0182 = (-6)*IT_0043;
    const ccomplex_t IT_0183 = (-6)*IT_0072;
    const ccomplex_t IT_0184 = IT_0180 + IT_0181 + IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = (-18)*conj(IT_0089);
    const ccomplex_t IT_0186 = (-6)*conj(IT_0028);
    const ccomplex_t IT_0187 = (-6)*conj(IT_0043);
    const ccomplex_t IT_0188 = (-6)*conj(IT_0072);
    const ccomplex_t IT_0189 = IT_0185 + IT_0186 + IT_0187 + IT_0188;
    const ccomplex_t IT_0190 = conj(IT_0082) + conj(IT_0092);
    const ccomplex_t IT_0191 = (-18)*IT_0097;
    const ccomplex_t IT_0192 = (-18)*IT_0098;
    const ccomplex_t IT_0193 = (-6)*IT_0075;
    const ccomplex_t IT_0194 = (-6)*IT_0076;
    const ccomplex_t IT_0195 = (-6)*IT_0077;
    const ccomplex_t IT_0196 = IT_0191 + IT_0192 + IT_0193 + IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = IT_0082 + IT_0092;
    const ccomplex_t IT_0198 = (-18)*conj(IT_0097);
    const ccomplex_t IT_0199 = (-18)*conj(IT_0098);
    const ccomplex_t IT_0200 = (-6)*conj(IT_0075);
    const ccomplex_t IT_0201 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0202 = (-6)*conj(IT_0077);
    const ccomplex_t IT_0203 = IT_0198 + IT_0199 + IT_0200 + IT_0201 + IT_0202;
    const ccomplex_t IT_0204 = m_d*IT_0082;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0119*conj(IT_0205);
    const ccomplex_t IT_0207 = IT_0108*conj(IT_0118);
    const ccomplex_t IT_0208 = IT_0205*IT_0207;
    const ccomplex_t IT_0209 = IT_0108*conj(IT_0115);
    const ccomplex_t IT_0210 = IT_0136 + IT_0209;
    const ccomplex_t IT_0211 = conj(IT_0121)*IT_0122;
    const ccomplex_t IT_0212 = IT_0207 + IT_0210;
    const ccomplex_t IT_0213 = IT_0122*conj(IT_0130);
    const ccomplex_t IT_0214 = IT_0111*IT_0122;
    const ccomplex_t IT_0215 = IT_0107*IT_0122;
    const ccomplex_t IT_0216 = IT_0108*conj(IT_0121);
    const ccomplex_t IT_0217 = conj(IT_0107)*IT_0122;
    const ccomplex_t IT_0218 = conj(IT_0111)*IT_0122;
    const ccomplex_t IT_0219 = IT_0108*IT_0121;
    const ccomplex_t IT_0220 = IT_0115*IT_0122;
    const ccomplex_t IT_0221 = IT_0108*IT_0130;
    const ccomplex_t IT_0222 = IT_0055*IT_0125;
    const ccomplex_t IT_0223 = cpow(IT_0038, 2);
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*m_N_2*IT_0223;
    const ccomplex_t IT_0225 = IT_0040*IT_0224;
    const ccomplex_t IT_0226 = IT_0222 + IT_0225;
    const ccomplex_t IT_0227 = IT_0108*conj(IT_0130);
    const ccomplex_t IT_0228 = conj(IT_0105)*IT_0108;
    const ccomplex_t IT_0229 = m_d*IT_0088;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0108*conj(IT_0230);
    const ccomplex_t IT_0232 = conj(IT_0115)*IT_0122;
    const ccomplex_t IT_0233 = conj(IT_0118)*IT_0122;
    const ccomplex_t IT_0234 = conj(IT_0105)*IT_0122;
    const ccomplex_t IT_0235 = IT_0122*IT_0230;
    const ccomplex_t IT_0236 = IT_0105*IT_0108;
    const ccomplex_t IT_0237 = IT_0108*IT_0230;
    const ccomplex_t IT_0238 = IT_0122*IT_0226;
    const ccomplex_t IT_0239 = IT_0025*IT_0128;
    const ccomplex_t IT_0240 = IT_0108*conj(IT_0205);
    const ccomplex_t IT_0241 = IT_0083*IT_0224;
    const ccomplex_t IT_0242 = IT_0108*conj(IT_0241);
    const ccomplex_t IT_0243 = cpow(IT_0068, 2);
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*m_N_1*IT_0243;
    const ccomplex_t IT_0245 = IT_0080*IT_0244;
    const ccomplex_t IT_0246 = IT_0108*conj(IT_0245);
    const ccomplex_t IT_0247 = IT_0122*conj(IT_0226);
    const ccomplex_t IT_0248 = IT_0058*IT_0244;
    const ccomplex_t IT_0249 = conj(IT_0239) + conj(IT_0248);
    const ccomplex_t IT_0250 = IT_0122*IT_0205;
    const ccomplex_t IT_0251 = conj(IT_0205)*IT_0250;
    const ccomplex_t IT_0252 = IT_0108*conj(IT_0226);
    const ccomplex_t IT_0253 = IT_0122*conj(IT_0230);
    const ccomplex_t IT_0254 = IT_0122*conj(IT_0205);
    const ccomplex_t IT_0255 = conj(IT_0241) + conj(IT_0245);
    const ccomplex_t IT_0256 = IT_0122*IT_0159;
    const ccomplex_t IT_0257 = IT_0122*IT_0155;
    const ccomplex_t IT_0258 = IT_0122*IT_0164;
    const ccomplex_t IT_0259 = IT_0108*IT_0226;
    const ccomplex_t IT_0260 = IT_0108*IT_0146;
    const ccomplex_t IT_0261 = -IT_0106 + -IT_0110 + -IT_0114 + 0.5*IT_0152 + 
      (-0.166666666666667)*IT_0154;
    const ccomplex_t IT_0262 = IT_0108*IT_0261;
    const ccomplex_t IT_0263 = cpow(IT_0066, 2);
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*m_N_1*IT_0263;
    const ccomplex_t IT_0265 = IT_0080*IT_0264;
    const ccomplex_t IT_0266 = IT_0122*IT_0261;
    const ccomplex_t IT_0267 = IT_0058*IT_0264;
    const ccomplex_t IT_0268 = IT_0108*IT_0205;
    const ccomplex_t IT_0269 = IT_0108*IT_0241;
    const ccomplex_t IT_0270 = IT_0108*IT_0245;
    const ccomplex_t IT_0271 = IT_0122*IT_0239;
    const ccomplex_t IT_0272 = IT_0122*IT_0245;
    const ccomplex_t IT_0273 = IT_0122*conj(IT_0159);
    const ccomplex_t IT_0274 = IT_0207 + IT_0273;
    const ccomplex_t IT_0275 = (-18)*s_34;
    const ccomplex_t IT_0276 = IT_0155*IT_0275;
    const ccomplex_t IT_0277 = IT_0108*conj(IT_0146);
    const ccomplex_t IT_0278 = IT_0108*conj(IT_0261);
    const ccomplex_t IT_0279 = IT_0122*conj(IT_0164);
    const ccomplex_t IT_0280 = (-6)*s_34;
    const ccomplex_t IT_0281 = IT_0108*conj(IT_0159);
    const ccomplex_t IT_0282 = IT_0108*conj(IT_0164);
    const ccomplex_t IT_0283 = IT_0233 + IT_0281 + IT_0282;
    const ccomplex_t IT_0284 = IT_0155*IT_0280;
    const ccomplex_t IT_0285 = IT_0122*conj(IT_0146);
    const ccomplex_t IT_0286 = m_d*s_24;
    const ccomplex_t IT_0287 = 6*IT_0121;
    const ccomplex_t IT_0288 = 6*IT_0130;
    const ccomplex_t IT_0289 = 6*IT_0105;
    const ccomplex_t IT_0290 = 6*IT_0230;
    const ccomplex_t IT_0291 = 6*IT_0205;
    const ccomplex_t IT_0292 = 18*IT_0107;
    const ccomplex_t IT_0293 = 18*IT_0111;
    const ccomplex_t IT_0294 = 18*IT_0115;
    const ccomplex_t IT_0295 = 18*IT_0118;
    const ccomplex_t IT_0296 = 18*IT_0226;
    const ccomplex_t IT_0297 = 6*conj(IT_0121);
    const ccomplex_t IT_0298 = 6*conj(IT_0130);
    const ccomplex_t IT_0299 = 6*conj(IT_0105);
    const ccomplex_t IT_0300 = 6*conj(IT_0230);
    const ccomplex_t IT_0301 = 6*conj(IT_0205);
    const ccomplex_t IT_0302 = 18*conj(IT_0107);
    const ccomplex_t IT_0303 = 18*conj(IT_0111);
    const ccomplex_t IT_0304 = 18*conj(IT_0115);
    const ccomplex_t IT_0305 = 18*conj(IT_0118);
    const ccomplex_t IT_0306 = 18*conj(IT_0226);
    const ccomplex_t IT_0307 = 6*IT_0107;
    const ccomplex_t IT_0308 = 6*IT_0111;
    const ccomplex_t IT_0309 = 6*IT_0115;
    const ccomplex_t IT_0310 = 6*IT_0118;
    const ccomplex_t IT_0311 = 6*IT_0226;
    const ccomplex_t IT_0312 = 18*IT_0121;
    const ccomplex_t IT_0313 = 18*IT_0130;
    const ccomplex_t IT_0314 = 18*IT_0105;
    const ccomplex_t IT_0315 = 18*IT_0230;
    const ccomplex_t IT_0316 = 18*IT_0205;
    const ccomplex_t IT_0317 = 6*conj(IT_0107);
    const ccomplex_t IT_0318 = 6*conj(IT_0111);
    const ccomplex_t IT_0319 = 6*conj(IT_0115);
    const ccomplex_t IT_0320 = 6*conj(IT_0118);
    const ccomplex_t IT_0321 = 6*conj(IT_0226);
    const ccomplex_t IT_0322 = 18*conj(IT_0121);
    const ccomplex_t IT_0323 = 18*conj(IT_0130);
    const ccomplex_t IT_0324 = 18*conj(IT_0105);
    const ccomplex_t IT_0325 = 18*conj(IT_0230);
    const ccomplex_t IT_0326 = 18*conj(IT_0205);
    const ccomplex_t IT_0327 = 3*IT_0265;
    const ccomplex_t IT_0328 = IT_0267 + IT_0327;
    const ccomplex_t IT_0329 = 6*IT_0146;
    const ccomplex_t IT_0330 = 6*IT_0261;
    const ccomplex_t IT_0331 = 18*IT_0155;
    const ccomplex_t IT_0332 = 18*IT_0159;
    const ccomplex_t IT_0333 = 18*IT_0164;
    const ccomplex_t IT_0334 = 6*IT_0155;
    const ccomplex_t IT_0335 = 6*IT_0159;
    const ccomplex_t IT_0336 = 6*IT_0164;
    const ccomplex_t IT_0337 = 18*IT_0146;
    const ccomplex_t IT_0338 = 18*IT_0261;
    const ccomplex_t IT_0339 = 6*IT_0089;
    const ccomplex_t IT_0340 = 6*IT_0092;
    const ccomplex_t IT_0341 = 18*IT_0028;
    const ccomplex_t IT_0342 = 18*IT_0043;
    const ccomplex_t IT_0343 = 18*IT_0072;
    const ccomplex_t IT_0344 = IT_0239 + IT_0248;
    const ccomplex_t IT_0345 = 6*conj(IT_0089);
    const ccomplex_t IT_0346 = 6*conj(IT_0092);
    const ccomplex_t IT_0347 = 18*conj(IT_0028);
    const ccomplex_t IT_0348 = 18*conj(IT_0043);
    const ccomplex_t IT_0349 = 18*conj(IT_0072);
    const ccomplex_t IT_0350 = 6*IT_0028;
    const ccomplex_t IT_0351 = 6*IT_0043;
    const ccomplex_t IT_0352 = 6*IT_0072;
    const ccomplex_t IT_0353 = 18*IT_0089;
    const ccomplex_t IT_0354 = 18*IT_0092;
    const ccomplex_t IT_0355 = IT_0241 + IT_0245;
    const ccomplex_t IT_0356 = 6*conj(IT_0028);
    const ccomplex_t IT_0357 = 6*conj(IT_0043);
    const ccomplex_t IT_0358 = 6*conj(IT_0072);
    const ccomplex_t IT_0359 = 18*conj(IT_0089);
    const ccomplex_t IT_0360 = 18*conj(IT_0092);
    const ccomplex_t IT_0361 = 3*conj(IT_0265);
    const ccomplex_t IT_0362 = conj(IT_0267) + IT_0361;
    const ccomplex_t IT_0363 = 6*conj(IT_0146);
    const ccomplex_t IT_0364 = 6*conj(IT_0261);
    const ccomplex_t IT_0365 = 18*conj(IT_0159);
    const ccomplex_t IT_0366 = 18*conj(IT_0164);
    const ccomplex_t IT_0367 = (-18)*IT_0155;
    const ccomplex_t IT_0368 = 6*conj(IT_0159);
    const ccomplex_t IT_0369 = 6*conj(IT_0164);
    const ccomplex_t IT_0370 = 18*conj(IT_0146);
    const ccomplex_t IT_0371 = 18*conj(IT_0261);
    const ccomplex_t IT_0372 = (-6)*IT_0155;
    const ccomplex_t IT_0373 = IT_0074*(IT_0287 + IT_0288 + IT_0289 + IT_0290 
      + IT_0291 + IT_0292 + IT_0293 + IT_0294 + IT_0295 + IT_0296) + IT_0073*
      (IT_0297 + IT_0298 + IT_0299 + IT_0300 + IT_0301 + IT_0302 + IT_0303 +
       IT_0304 + IT_0305 + IT_0306) + IT_0095*(IT_0307 + IT_0308 + IT_0309 +
       IT_0310 + IT_0311 + IT_0312 + IT_0313 + IT_0314 + IT_0315 + IT_0316) + 6
      *conj(IT_0082)*(IT_0239 + 3*IT_0241 + 3*IT_0245 + IT_0248 + 1./6*IT_0307 +
       1./6*IT_0308 + 1./6*IT_0309 + 1./6*IT_0310 + 1./6*IT_0311 + 1./6*IT_0312 
      + 1./6*IT_0313 + 1./6*IT_0314 + 1./6*IT_0315 + 1./6*IT_0316) + IT_0093*
      (IT_0317 + IT_0318 + IT_0319 + IT_0320 + IT_0321 + IT_0322 + IT_0323 +
       IT_0324 + IT_0325 + IT_0326) + 6*IT_0082*(conj(IT_0239) + 3*conj(IT_0241)
       + 3*conj(IT_0245) + conj(IT_0248) + 1./6*IT_0317 + 1./6*IT_0318 + 1./6
      *IT_0319 + 1./6*IT_0320 + 1./6*IT_0321 + 1./6*IT_0322 + 1./6*IT_0323 + 1.
      /6*IT_0324 + 1./6*IT_0325 + 1./6*IT_0326) + 6*IT_0099*(0.166666666666667
      *IT_0310 + 0.166666666666667*IT_0316 + IT_0328 + 0.166666666666667*IT_0329
       + 0.166666666666667*IT_0330 + 0.166666666666667*IT_0331 +
       0.166666666666667*IT_0332 + 0.166666666666667*IT_0333) + 6*IT_0079*
      (IT_0265 + 3*IT_0267 + 1./6*IT_0291 + 1./6*IT_0295 + 1./6*IT_0334 + 1./6
      *IT_0335 + 1./6*IT_0336 + 1./6*IT_0337 + 1./6*IT_0338) + IT_0249*(IT_0339 
      + IT_0340 + IT_0341 + IT_0342 + IT_0343) + IT_0344*(IT_0345 + IT_0346 +
       IT_0347 + IT_0348 + IT_0349) + IT_0255*(IT_0350 + IT_0351 + IT_0352 +
       IT_0353 + IT_0354) + IT_0355*(IT_0356 + IT_0357 + IT_0358 + IT_0359 +
       IT_0360) + 6*IT_0103*(0.166666666666667*IT_0320 + 0.166666666666667
      *IT_0326 + IT_0362 + 0.166666666666667*IT_0363 + 0.166666666666667*IT_0364
       + 0.166666666666667*IT_0365 + 0.166666666666667*IT_0366 +
       0.166666666666667*IT_0367) + 6*IT_0078*(conj(IT_0265) + 3*conj(IT_0267) +
       1./6*IT_0301 + 1./6*IT_0305 + 1./6*IT_0368 + 1./6*IT_0369 + 1./6*IT_0370 
      + 1./6*IT_0371 + 1./6*IT_0372);
    const ccomplex_t IT_0374 = (-6)*IT_0118;
    const ccomplex_t IT_0375 = IT_0024*conj(IT_0205)*IT_0374;
    const ccomplex_t IT_0376 = (-6)*conj(IT_0118);
    const ccomplex_t IT_0377 = IT_0024*IT_0205*IT_0376;
    const ccomplex_t IT_0378 = (-18)*conj(IT_0118);
    const ccomplex_t IT_0379 = (-18)*conj(IT_0146);
    const ccomplex_t IT_0380 = (-18)*conj(IT_0261);
    const ccomplex_t IT_0381 = (-6)*conj(IT_0205);
    const ccomplex_t IT_0382 = (-6)*conj(IT_0159);
    const ccomplex_t IT_0383 = (-6)*conj(IT_0164);
    const ccomplex_t IT_0384 = IT_0334 + IT_0378 + IT_0379 + IT_0380 + IT_0381
       + IT_0382 + IT_0383;
    const ccomplex_t IT_0385 = (-18)*conj(IT_0205);
    const ccomplex_t IT_0386 = (-18)*conj(IT_0159);
    const ccomplex_t IT_0387 = (-18)*conj(IT_0164);
    const ccomplex_t IT_0388 = (-6)*conj(IT_0146);
    const ccomplex_t IT_0389 = (-6)*conj(IT_0261);
    const ccomplex_t IT_0390 = IT_0331 + IT_0376 + IT_0385 + IT_0386 + IT_0387
       + IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = (-18)*conj(IT_0107);
    const ccomplex_t IT_0392 = (-18)*conj(IT_0111);
    const ccomplex_t IT_0393 = (-18)*conj(IT_0115);
    const ccomplex_t IT_0394 = (-18)*conj(IT_0226);
    const ccomplex_t IT_0395 = (-6)*conj(IT_0121);
    const ccomplex_t IT_0396 = (-6)*conj(IT_0130);
    const ccomplex_t IT_0397 = (-6)*conj(IT_0105);
    const ccomplex_t IT_0398 = (-6)*conj(IT_0230);
    const ccomplex_t IT_0399 = IT_0378 + IT_0381 + IT_0391 + IT_0392 + IT_0393
       + IT_0394 + IT_0395 + IT_0396 + IT_0397 + IT_0398;
    const ccomplex_t IT_0400 = (-18)*IT_0205;
    const ccomplex_t IT_0401 = (-18)*IT_0159;
    const ccomplex_t IT_0402 = (-18)*IT_0164;
    const ccomplex_t IT_0403 = (-6)*IT_0146;
    const ccomplex_t IT_0404 = (-6)*IT_0261;
    const ccomplex_t IT_0405 = IT_0367 + IT_0374 + IT_0400 + IT_0401 + IT_0402
       + IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = (-18)*IT_0118;
    const ccomplex_t IT_0407 = (-18)*IT_0146;
    const ccomplex_t IT_0408 = (-18)*IT_0261;
    const ccomplex_t IT_0409 = (-6)*IT_0205;
    const ccomplex_t IT_0410 = (-6)*IT_0159;
    const ccomplex_t IT_0411 = (-6)*IT_0164;
    const ccomplex_t IT_0412 = IT_0372 + IT_0406 + IT_0407 + IT_0408 + IT_0409
       + IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = (-18)*IT_0107;
    const ccomplex_t IT_0414 = (-18)*IT_0111;
    const ccomplex_t IT_0415 = (-18)*IT_0115;
    const ccomplex_t IT_0416 = IT_0413 + IT_0414 + IT_0415;
    const ccomplex_t IT_0417 = (-6)*IT_0121;
    const ccomplex_t IT_0418 = (-6)*IT_0130;
    const ccomplex_t IT_0419 = (-18)*IT_0121;
    const ccomplex_t IT_0420 = (-18)*IT_0130;
    const ccomplex_t IT_0421 = (-18)*IT_0105;
    const ccomplex_t IT_0422 = (-18)*IT_0230;
    const ccomplex_t IT_0423 = (-6)*IT_0107;
    const ccomplex_t IT_0424 = (-6)*IT_0111;
    const ccomplex_t IT_0425 = (-6)*IT_0115;
    const ccomplex_t IT_0426 = IT_0400 + IT_0419 + IT_0420 + IT_0421 + IT_0422
       + IT_0423 + IT_0424 + IT_0425;
    const ccomplex_t IT_0427 = (-6)*IT_0226;
    const ccomplex_t IT_0428 = (-18)*conj(IT_0121);
    const ccomplex_t IT_0429 = (-18)*conj(IT_0130);
    const ccomplex_t IT_0430 = (-18)*conj(IT_0105);
    const ccomplex_t IT_0431 = (-18)*conj(IT_0230);
    const ccomplex_t IT_0432 = (-6)*conj(IT_0107);
    const ccomplex_t IT_0433 = (-6)*conj(IT_0111);
    const ccomplex_t IT_0434 = (-6)*conj(IT_0115);
    const ccomplex_t IT_0435 = IT_0385 + IT_0428 + IT_0429 + IT_0430 + IT_0431
       + IT_0432 + IT_0433 + IT_0434;
    const ccomplex_t IT_0436 = (-6)*conj(IT_0226);
    const ccomplex_t IT_0437 = IT_0155*(IT_0307 + IT_0308 + IT_0309 + IT_0310 
      + IT_0311 + IT_0312 + IT_0313 + IT_0314 + IT_0315 + IT_0316 + IT_0376 +
       IT_0385 + IT_0428 + IT_0429 + IT_0430 + IT_0431 + IT_0432 + IT_0433 +
       IT_0434 + IT_0436);
    const ccomplex_t IT_0438 = (-18)*IT_0226;
    const ccomplex_t IT_0439 = (-6)*IT_0105;
    const ccomplex_t IT_0440 = (-6)*IT_0230;
    const ccomplex_t IT_0441 = 3*IT_0239 + IT_0241 + IT_0245 + 3*IT_0248 + (
      -0.166666666666667)*IT_0406 + (-0.166666666666667)*IT_0409 + (
      -0.166666666666667)*IT_0413 + (-0.166666666666667)*IT_0414 + (
      -0.166666666666667)*IT_0415 + (-0.166666666666667)*IT_0417 + (
      -0.166666666666667)*IT_0418 + (-0.166666666666667)*IT_0438 + (
      -0.166666666666667)*IT_0439 + (-0.166666666666667)*IT_0440;
    const ccomplex_t IT_0442 = 3*conj(IT_0239) + conj(IT_0241) + conj(IT_0245)
       + 3*conj(IT_0248) + (-0.166666666666667)*IT_0378 + (-0.166666666666667)
      *IT_0381 + (-0.166666666666667)*IT_0391 + (-0.166666666666667)*IT_0392 + (
      -0.166666666666667)*IT_0393 + (-0.166666666666667)*IT_0394 + (
      -0.166666666666667)*IT_0395 + (-0.166666666666667)*IT_0396 + (
      -0.166666666666667)*IT_0397 + (-0.166666666666667)*IT_0398;
    const ccomplex_t IT_0443 = (conj(IT_0105) + conj(IT_0230))*IT_0374 + 
      (IT_0105 + IT_0230)*IT_0376 + IT_0226*IT_0378 + IT_0344*IT_0384 + IT_0355
      *IT_0390 + IT_0118*((-36)*conj(IT_0118) + IT_0391 + IT_0392 + IT_0393 +
       IT_0395 + IT_0396) + (IT_0146 + IT_0261)*IT_0399 + IT_0255*IT_0405 + conj
      (IT_0226)*IT_0406 + IT_0249*IT_0412 + conj(IT_0118)*(IT_0416 + IT_0417 +
       IT_0418) + conj(IT_0205)*(IT_0426 + IT_0427) + (conj(IT_0159) + conj
      (IT_0164))*(IT_0374 + IT_0426 + IT_0427) + (-6)*conj(IT_0265)*(IT_0239 + 3
      *IT_0241 + 3*IT_0245 + IT_0248 + -1./6*IT_0374 + -1./6*IT_0400 + -1./6
      *IT_0419 + -1./6*IT_0420 + -1./6*IT_0421 + -1./6*IT_0422 + -1./6*IT_0423 +
       -1./6*IT_0424 + -1./6*IT_0425 + -1./6*IT_0427) + IT_0205*(IT_0435 +
       IT_0436) + (IT_0159 + IT_0164)*(IT_0376 + IT_0435 + IT_0436) + (-6)
      *IT_0265*(conj(IT_0239) + 3*conj(IT_0241) + 3*conj(IT_0245) + conj(IT_0248
      ) + -1./6*IT_0376 + -1./6*IT_0385 + -1./6*IT_0428 + -1./6*IT_0429 + -1./6
      *IT_0430 + -1./6*IT_0431 + -1./6*IT_0432 + -1./6*IT_0433 + -1./6*IT_0434 +
       -1./6*IT_0436) + IT_0437 + (conj(IT_0146) + conj(IT_0261))*(IT_0406 +
       IT_0409 + IT_0416 + IT_0417 + IT_0418 + IT_0438 + IT_0439 + IT_0440) + (
      -6)*conj(IT_0267)*IT_0441 + (-6)*IT_0267*IT_0442;
    const ccomplex_t IT_0444 = m_d*s_23;
    const ccomplex_t IT_0445 = (-18)*IT_0075;
    const ccomplex_t IT_0446 = (-18)*IT_0076;
    const ccomplex_t IT_0447 = (-18)*IT_0077;
    const ccomplex_t IT_0448 = (-6)*IT_0097;
    const ccomplex_t IT_0449 = (-6)*IT_0098;
    const ccomplex_t IT_0450 = (-18)*conj(IT_0075);
    const ccomplex_t IT_0451 = (-18)*conj(IT_0076);
    const ccomplex_t IT_0452 = (-18)*conj(IT_0077);
    const ccomplex_t IT_0453 = (-6)*conj(IT_0097);
    const ccomplex_t IT_0454 = (-6)*conj(IT_0098);
    const ccomplex_t IT_0455 = IT_0196*IT_0255 + IT_0189*IT_0265 + IT_0184
      *conj(IT_0265) + IT_0179*IT_0267 + IT_0174*conj(IT_0267) + IT_0203*IT_0355
       + IT_0073*IT_0384 + (-6)*IT_0197*((-0.166666666666667)*IT_0331 + IT_0362 
      + (-0.166666666666667)*IT_0376 + (-0.166666666666667)*IT_0385 + (
      -0.166666666666667)*IT_0386 + (-0.166666666666667)*IT_0387 + (
      -0.166666666666667)*IT_0388 + (-0.166666666666667)*IT_0389) + IT_0089
      *IT_0390 + IT_0078*IT_0399 + (-6)*IT_0190*(IT_0328 + (-0.166666666666667)
      *IT_0367 + (-0.166666666666667)*IT_0374 + (-0.166666666666667)*IT_0400 + (
      -0.166666666666667)*IT_0401 + (-0.166666666666667)*IT_0402 + (
      -0.166666666666667)*IT_0403 + (-0.166666666666667)*IT_0404) + conj(IT_0089
      )*IT_0405 + IT_0074*IT_0412 + IT_0099*(IT_0374 + IT_0426 + IT_0427) +
       IT_0103*(IT_0376 + IT_0435 + IT_0436) + IT_0079*(IT_0406 + IT_0409 +
       IT_0416 + IT_0417 + IT_0418 + IT_0438 + IT_0439 + IT_0440) + IT_0249*
      (IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449) + IT_0344*(IT_0450 +
       IT_0451 + IT_0452 + IT_0453 + IT_0454);
    const ccomplex_t IT_0456 = IT_0005*(IT_0073*IT_0074 + IT_0078*IT_0079 +
       IT_0094*IT_0096 + (IT_0097 + IT_0098)*IT_0099) + IT_0102*(IT_0074*IT_0094
       + IT_0073*IT_0096 + IT_0078*IT_0099 + IT_0079*IT_0103) + conj(IT_0121)*
      (IT_0113 + IT_0123) + conj(IT_0130)*(IT_0113 + IT_0116 + IT_0123 + IT_0131
      ) + conj(IT_0105)*IT_0133 + IT_0121*IT_0136 + conj(IT_0146)*IT_0168 +
       IT_0169*(IT_0079*IT_0174 + IT_0078*IT_0179 + IT_0099*IT_0184 + IT_0103
      *IT_0189 + IT_0190*IT_0196 + IT_0197*IT_0203) + 2*IT_0206 + IT_0159
      *IT_0207 + 2*IT_0208 + IT_0130*(IT_0210 + IT_0211) + IT_0105*(IT_0211 +
       IT_0212 + IT_0213) + conj(IT_0111)*IT_0214 + conj(IT_0107)*(IT_0214 +
       IT_0215) + IT_0115*(IT_0216 + IT_0217 + IT_0218) + IT_0107*IT_0218 + conj
      (IT_0115)*(IT_0214 + IT_0215 + IT_0219 + IT_0220) + conj(IT_0118)*(IT_0214
       + IT_0215 + IT_0219 + IT_0220 + IT_0221) + IT_0226*(IT_0216 + IT_0217 +
       IT_0218 + IT_0227 + IT_0228 + IT_0231 + IT_0232 + IT_0233) + IT_0230*
      (IT_0211 + IT_0212 + IT_0213 + IT_0234) + conj(IT_0230)*(IT_0133 + IT_0235
      ) + conj(IT_0226)*(IT_0166 + IT_0214 + IT_0215 + IT_0219 + IT_0220 +
       IT_0221 + IT_0236 + IT_0237 + IT_0238) + IT_0239*(IT_0216 + IT_0217 +
       IT_0218 + IT_0227 + IT_0228 + IT_0231 + IT_0232 + IT_0233 + IT_0240 +
       IT_0242 + IT_0246 + IT_0247) + IT_0248*(IT_0216 + IT_0217 + IT_0218 +
       IT_0227 + IT_0228 + IT_0231 + IT_0232 + IT_0233 + IT_0240 + IT_0242 +
       IT_0246 + IT_0247 + IT_0122*IT_0249) + 2*IT_0251 + IT_0245*(IT_0211 +
       IT_0212 + IT_0213 + IT_0234 + IT_0252 + IT_0253 + IT_0254) + IT_0241*
      (IT_0211 + IT_0212 + IT_0213 + IT_0234 + IT_0252 + IT_0253 + IT_0254 +
       IT_0122*IT_0255) + conj(IT_0159)*(IT_0119 + IT_0256) + conj(IT_0164)*
      (IT_0119 + IT_0256 + IT_0257 + IT_0258) + conj(IT_0205)*(IT_0113 + IT_0116
       + IT_0123 + IT_0131 + IT_0132 + IT_0235 + IT_0256 + IT_0257 + IT_0258 +
       IT_0259 + IT_0260 + IT_0262) + (IT_0119 + IT_0250 + IT_0256 + IT_0257 +
       IT_0258 + IT_0260 + IT_0262)*conj(IT_0265) + conj(IT_0261)*(IT_0168 +
       IT_0266) + conj(IT_0267)*(IT_0156 + IT_0160 + IT_0165 + IT_0166 + IT_0167
       + IT_0266 + IT_0268) + IT_0249*(IT_0166 + IT_0214 + IT_0215 + IT_0219 +
       IT_0220 + IT_0221 + IT_0236 + IT_0237 + IT_0238 + IT_0268 + IT_0269 +
       IT_0270 + IT_0271) + IT_0255*(IT_0113 + IT_0116 + IT_0119 + IT_0123 +
       IT_0131 + IT_0132 + IT_0235 + IT_0250 + IT_0259 + IT_0272) + IT_0155*
      (IT_0274 + (IT_0155 + IT_0159)*IT_0275) + IT_0164*(IT_0274 + IT_0276) +
       IT_0205*(IT_0210 + IT_0211 + IT_0213 + IT_0234 + IT_0252 + IT_0253 +
       IT_0273 + IT_0276 + IT_0277 + IT_0278 + IT_0279) + IT_0265*(IT_0207 +
       IT_0254 + IT_0122*conj(IT_0265) + IT_0108*conj(IT_0267) + IT_0273 +
       IT_0276 + IT_0277 + IT_0278 + IT_0279) + IT_0118*(36*s_34*conj(IT_0118) +
       IT_0216 + IT_0217 + IT_0218 + IT_0227 + IT_0232 + IT_0155*IT_0280) +
       IT_0146*(IT_0283 + IT_0284) + IT_0267*(IT_0233 + IT_0240 + IT_0108*conj
      (IT_0265) + IT_0122*(conj(IT_0261) + conj(IT_0267)) + IT_0281 + IT_0282 +
       IT_0284 + IT_0285) + IT_0261*(IT_0283 + IT_0284 + IT_0285) + IT_0286
      *IT_0373 + 2*IT_0375 + 2*IT_0377 + IT_0024*IT_0443 + IT_0444*IT_0455;
    return create_ccomplex_return(IT_0456);
}

