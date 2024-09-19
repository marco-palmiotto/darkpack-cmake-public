#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_mu_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t IT_0000 = m_mu*m_N_1;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_B1)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_W1)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0009 + IT_0014);
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = N_B1*e_em;
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = N_W1*e_em;
    const ccomplex_t IT_0021 = IT_0011*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + IT_0022);
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0016*IT_0024;
    const ccomplex_t IT_0026 = pow(m_mu, 2);
    const ccomplex_t IT_0027 = pow(m_N_1, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_13 + IT_0026 + IT_0027 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0011*IT_0031*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0011*IT_0031*IT_0033;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0026 + IT_0027 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0030 + -IT_0041;
    const ccomplex_t IT_0043 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0044 = cos(alpha);
    const ccomplex_t IT_0045 = IT_0006*IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0048 = sin(alpha);
    const ccomplex_t IT_0049 = IT_0006*IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0052 = IT_0011*IT_0044;
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0055 = IT_0011*IT_0048;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0046 + IT_0050 + 
      -IT_0053 + -IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*e_em*m_mu*IT_0011*IT_0031
      *IT_0033*IT_0048;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0058*IT_0063;
    const ccomplex_t IT_0065 = IT_0045*IT_0047;
    const ccomplex_t IT_0066 = IT_0043*IT_0049;
    const ccomplex_t IT_0067 = IT_0052*IT_0054;
    const ccomplex_t IT_0068 = IT_0051*IT_0055;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + -IT_0066 + 
      -IT_0067 + IT_0068);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*e_em*m_mu*IT_0011*IT_0031
      *IT_0033*IT_0044;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_12 + (-2)*IT_0027 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0069*IT_0074;
    const ccomplex_t IT_0076 = (-0.5)*IT_0064 + (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0064 + 0.5*IT_0075;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = IT_0076 + IT_0078;
    const ccomplex_t IT_0080 = -conj(IT_0077);
    const ccomplex_t IT_0081 = conj(IT_0076) + IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0006;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0006;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = IT_0039*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0034;
    const ccomplex_t IT_0090 = (-0.5)*IT_0036;
    const ccomplex_t IT_0091 = IT_0089*IT_0090;
    const ccomplex_t IT_0092 = IT_0028*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = -IT_0088 + -IT_0093;
    const ccomplex_t IT_0095 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0096 = IT_0045*IT_0095;
    const ccomplex_t IT_0097 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0098 = IT_0049*IT_0097;
    const ccomplex_t IT_0099 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0100 = IT_0052*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0102 = IT_0055*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0096 + IT_0098 + 
      -IT_0100 + -IT_0102);
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = IT_0063*IT_0104;
    const ccomplex_t IT_0106 = IT_0045*IT_0097;
    const ccomplex_t IT_0107 = IT_0049*IT_0095;
    const ccomplex_t IT_0108 = IT_0052*IT_0101;
    const ccomplex_t IT_0109 = IT_0055*IT_0099;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + -IT_0107 + 
      -IT_0108 + IT_0109);
    const ccomplex_t IT_0111 = IT_0074*IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0105 + (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = 0.5*IT_0105 + 0.5*IT_0111;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0112 + IT_0114;
    const ccomplex_t IT_0116 = -conj(IT_0113);
    const ccomplex_t IT_0117 = conj(IT_0112) + IT_0116;
    const ccomplex_t IT_0118 = -s_24;
    const ccomplex_t IT_0119 = s_23 + IT_0118;
    const ccomplex_t IT_0120 = IT_0000*IT_0119;
    const ccomplex_t IT_0121 = (-2)*IT_0120;
    const ccomplex_t IT_0122 = 4*IT_0003;
    const ccomplex_t IT_0123 = IT_0006*IT_0010;
    const ccomplex_t IT_0124 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = IT_0005*IT_0011;
    const ccomplex_t IT_0127 = IT_0124*IT_0126;
    const ccomplex_t IT_0128 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0129 = IT_0123*IT_0128;
    const ccomplex_t IT_0130 = IT_0126*IT_0128;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0125 + IT_0127 + 
      -IT_0129 + -IT_0130);
    const ccomplex_t IT_0132 = 0.5*IT_0131;
    const ccomplex_t IT_0133 = e_em*IT_0123;
    const ccomplex_t IT_0134 = e_em*IT_0126;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0133 + -IT_0134);
    const ccomplex_t IT_0136 = 0.5*IT_0135;
    const ccomplex_t IT_0137 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0138 = IT_0136*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = IT_0132*IT_0139;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0010;
    const ccomplex_t IT_0143 = IT_0137*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0132*IT_0144;
    const ccomplex_t IT_0146 = 0.5*IT_0145;
    const ccomplex_t IT_0147 = (-0.5)*IT_0131;
    const ccomplex_t IT_0148 = IT_0139*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = IT_0144*IT_0147;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0149 + IT_0151;
    const ccomplex_t IT_0153 = 2*IT_0120;
    const ccomplex_t IT_0154 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0155 = IT_0086*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = cpow((-2)*s_23 + IT_0026 + IT_0027 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0158 = IT_0091*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = IT_0156 + IT_0159;
    const ccomplex_t IT_0161 = -IT_0112;
    const ccomplex_t IT_0162 = IT_0113 + IT_0161;
    const ccomplex_t IT_0163 = -conj(IT_0112);
    const ccomplex_t IT_0164 = conj(IT_0113) + IT_0163;
    const ccomplex_t IT_0165 = IT_0025*IT_0157;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = IT_0038*IT_0154;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0166 + IT_0168;
    const ccomplex_t IT_0170 = -IT_0076;
    const ccomplex_t IT_0171 = IT_0077 + IT_0170;
    const ccomplex_t IT_0172 = -conj(IT_0076);
    const ccomplex_t IT_0173 = conj(IT_0077) + IT_0172;
    const ccomplex_t IT_0174 = (-2)*IT_0003;
    const ccomplex_t IT_0175 = 4*IT_0120;
    const ccomplex_t IT_0176 = (-0.5)*IT_0148;
    const ccomplex_t IT_0177 = (-0.5)*IT_0150;
    const ccomplex_t IT_0178 = (-0.5)*IT_0140;
    const ccomplex_t IT_0179 = (-0.5)*IT_0145;
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = s_12*IT_0026;
    const ccomplex_t IT_0182 = s_13*s_24;
    const ccomplex_t IT_0183 = s_12*s_34;
    const ccomplex_t IT_0184 = s_14*s_23;
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = IT_0181 + IT_0182 + IT_0183 + IT_0185;
    const ccomplex_t IT_0187 = IT_0024*IT_0090;
    const ccomplex_t IT_0188 = IT_0028*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0037*IT_0085;
    const ccomplex_t IT_0191 = IT_0039*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = IT_0189 + IT_0192;
    const ccomplex_t IT_0194 = IT_0006*IT_0032;
    const ccomplex_t IT_0195 = IT_0047*IT_0194;
    const ccomplex_t IT_0196 = sin(beta);
    const ccomplex_t IT_0197 = IT_0006*IT_0196;
    const ccomplex_t IT_0198 = IT_0043*IT_0197;
    const ccomplex_t IT_0199 = IT_0011*IT_0032;
    const ccomplex_t IT_0200 = IT_0054*IT_0199;
    const ccomplex_t IT_0201 = IT_0011*IT_0196;
    const ccomplex_t IT_0202 = IT_0051*IT_0201;
    const ccomplex_t IT_0203 = -IT_0195 + -IT_0198 + IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = e_em*m_mu*IT_0011*IT_0031;
    const ccomplex_t IT_0205 = 0.5*IT_0204;
    const ccomplex_t IT_0206 = pow(m_Z, 2);
    const ccomplex_t IT_0207 = cpow((-2)*s_12 + (-2)*IT_0027 + IT_0206 + 
      -reg_prop, -1);
    const ccomplex_t IT_0208 = IT_0205*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = IT_0203*IT_0209;
    const ccomplex_t IT_0211 = IT_0043*IT_0194;
    const ccomplex_t IT_0212 = IT_0047*IT_0197;
    const ccomplex_t IT_0213 = IT_0051*IT_0199;
    const ccomplex_t IT_0214 = IT_0054*IT_0201;
    const ccomplex_t IT_0215 = -IT_0211 + IT_0212 + IT_0213 + -IT_0214;
    const ccomplex_t IT_0216 = e_em*m_mu*IT_0011*IT_0031*IT_0033*IT_0196;
    const ccomplex_t IT_0217 = (-0.5)*IT_0216;
    const ccomplex_t IT_0218 = cpow((-2)*s_12 + (-2)*IT_0027 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0219 = IT_0217*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = IT_0215*IT_0220;
    const ccomplex_t IT_0222 = 0.5*IT_0210 + 0.5*IT_0221;
    const ccomplex_t IT_0223 = 2*IT_0222;
    const ccomplex_t IT_0224 = (-0.5)*IT_0210 + (-0.5)*IT_0221;
    const ccomplex_t IT_0225 = (-2)*IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = 2*conj(IT_0222);
    const ccomplex_t IT_0228 = (-2)*conj(IT_0224);
    const ccomplex_t IT_0229 = IT_0227 + IT_0228;
    const ccomplex_t IT_0230 = IT_0016*IT_0089;
    const ccomplex_t IT_0231 = IT_0028*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = IT_0035*IT_0083;
    const ccomplex_t IT_0234 = IT_0039*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*IT_0234;
    const ccomplex_t IT_0236 = IT_0232 + IT_0235;
    const ccomplex_t IT_0237 = IT_0097*IT_0194;
    const ccomplex_t IT_0238 = IT_0095*IT_0197;
    const ccomplex_t IT_0239 = IT_0101*IT_0199;
    const ccomplex_t IT_0240 = IT_0099*IT_0201;
    const ccomplex_t IT_0241 = IT_0237 + IT_0238 + -IT_0239 + -IT_0240;
    const ccomplex_t IT_0242 = IT_0209*IT_0241;
    const ccomplex_t IT_0243 = IT_0095*IT_0194;
    const ccomplex_t IT_0244 = IT_0097*IT_0197;
    const ccomplex_t IT_0245 = IT_0099*IT_0199;
    const ccomplex_t IT_0246 = IT_0101*IT_0201;
    const ccomplex_t IT_0247 = IT_0243 + -IT_0244 + -IT_0245 + IT_0246;
    const ccomplex_t IT_0248 = IT_0220*IT_0247;
    const ccomplex_t IT_0249 = (-0.5)*IT_0242 + (-0.5)*IT_0248;
    const ccomplex_t IT_0250 = 2*IT_0249;
    const ccomplex_t IT_0251 = 0.5*IT_0242 + 0.5*IT_0248;
    const ccomplex_t IT_0252 = (-2)*IT_0251;
    const ccomplex_t IT_0253 = IT_0250 + IT_0252;
    const ccomplex_t IT_0254 = 2*conj(IT_0249);
    const ccomplex_t IT_0255 = (-2)*conj(IT_0251);
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = -IT_0182;
    const ccomplex_t IT_0258 = IT_0181 + IT_0183 + IT_0184 + IT_0257;
    const ccomplex_t IT_0259 = IT_0157*IT_0187;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = IT_0154*IT_0190;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = -IT_0260 + -IT_0262;
    const ccomplex_t IT_0264 = 2*IT_0224;
    const ccomplex_t IT_0265 = (-2)*IT_0222;
    const ccomplex_t IT_0266 = IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = 2*conj(IT_0224);
    const ccomplex_t IT_0268 = (-2)*conj(IT_0222);
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = IT_0157*IT_0230;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = IT_0154*IT_0233;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = -IT_0271 + -IT_0273;
    const ccomplex_t IT_0275 = 2*IT_0251;
    const ccomplex_t IT_0276 = (-2)*IT_0249;
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = 2*conj(IT_0251);
    const ccomplex_t IT_0279 = (-2)*conj(IT_0249);
    const ccomplex_t IT_0280 = IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = s_23 + s_24;
    const ccomplex_t IT_0282 = IT_0000*IT_0281;
    const ccomplex_t IT_0283 = conj(IT_0149) + conj(IT_0178);
    const ccomplex_t IT_0284 = 4*IT_0224;
    const ccomplex_t IT_0285 = (-4)*IT_0222;
    const ccomplex_t IT_0286 = IT_0284 + IT_0285;
    const ccomplex_t IT_0287 = IT_0149 + IT_0178;
    const ccomplex_t IT_0288 = 4*conj(IT_0224);
    const ccomplex_t IT_0289 = (-4)*conj(IT_0222);
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = conj(IT_0151) + conj(IT_0179);
    const ccomplex_t IT_0292 = 4*IT_0222;
    const ccomplex_t IT_0293 = (-4)*IT_0224;
    const ccomplex_t IT_0294 = IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = IT_0151 + IT_0179;
    const ccomplex_t IT_0296 = 4*conj(IT_0222);
    const ccomplex_t IT_0297 = (-4)*conj(IT_0224);
    const ccomplex_t IT_0298 = IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = conj(IT_0141) + conj(IT_0176);
    const ccomplex_t IT_0300 = 4*IT_0249;
    const ccomplex_t IT_0301 = (-4)*IT_0251;
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = IT_0141 + IT_0176;
    const ccomplex_t IT_0304 = 4*conj(IT_0249);
    const ccomplex_t IT_0305 = (-4)*conj(IT_0251);
    const ccomplex_t IT_0306 = IT_0304 + IT_0305;
    const ccomplex_t IT_0307 = conj(IT_0146) + conj(IT_0177);
    const ccomplex_t IT_0308 = 4*IT_0251;
    const ccomplex_t IT_0309 = (-4)*IT_0249;
    const ccomplex_t IT_0310 = IT_0308 + IT_0309;
    const ccomplex_t IT_0311 = IT_0146 + IT_0177;
    const ccomplex_t IT_0312 = 4*conj(IT_0251);
    const ccomplex_t IT_0313 = (-4)*conj(IT_0249);
    const ccomplex_t IT_0314 = IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = conj(IT_0094)*IT_0226 + IT_0094*IT_0229 + conj
      (IT_0042)*IT_0253 + IT_0042*IT_0256 + conj(IT_0169)*IT_0266 + IT_0169
      *IT_0269 + conj(IT_0160)*IT_0277 + IT_0160*IT_0280 + IT_0283*IT_0286 +
       IT_0287*IT_0290 + IT_0291*IT_0294 + IT_0295*IT_0298 + IT_0299*IT_0302 +
       IT_0303*IT_0306 + IT_0307*IT_0310 + IT_0311*IT_0314;
    const ccomplex_t IT_0316 = s_13 + s_14;
    const ccomplex_t IT_0317 = IT_0000*IT_0316;
    const ccomplex_t IT_0318 = 2*IT_0094;
    const ccomplex_t IT_0319 = 2*conj(IT_0094);
    const ccomplex_t IT_0320 = 2*IT_0169;
    const ccomplex_t IT_0321 = 2*conj(IT_0169);
    const ccomplex_t IT_0322 = conj(IT_0042)*IT_0226 + IT_0042*IT_0229 + conj
      (IT_0160)*IT_0266 + IT_0160*IT_0269 + IT_0294*IT_0299 + IT_0291*IT_0302 +
       IT_0298*IT_0303 + IT_0295*IT_0306 + IT_0286*IT_0307 + IT_0283*IT_0310 +
       IT_0290*IT_0311 + IT_0287*IT_0314 + (-2)*conj(IT_0249)*(IT_0169 + (-0.5)
      *IT_0318) + (-2)*IT_0249*(conj(IT_0169) + (-0.5)*IT_0319) + (-2)*conj
      (IT_0251)*(IT_0094 + (-0.5)*IT_0320) + (-2)*IT_0251*(conj(IT_0094) + (-0.5
      )*IT_0321);
    const ccomplex_t IT_0323 = m_mu*s_13*m_N_1;
    const ccomplex_t IT_0324 = 4*IT_0094;
    const ccomplex_t IT_0325 = 0.125*IT_0324;
    const ccomplex_t IT_0326 = IT_0295 + IT_0325;
    const ccomplex_t IT_0327 = 4*conj(IT_0094);
    const ccomplex_t IT_0328 = 0.125*IT_0327;
    const ccomplex_t IT_0329 = IT_0291 + IT_0328;
    const ccomplex_t IT_0330 = 4*IT_0263;
    const ccomplex_t IT_0331 = 0.125*IT_0330;
    const ccomplex_t IT_0332 = IT_0193 + IT_0331;
    const ccomplex_t IT_0333 = 4*conj(IT_0263);
    const ccomplex_t IT_0334 = 0.125*IT_0333;
    const ccomplex_t IT_0335 = conj(IT_0193) + IT_0334;
    const ccomplex_t IT_0336 = 2*IT_0263;
    const ccomplex_t IT_0337 = 4*IT_0193;
    const ccomplex_t IT_0338 = IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = 2*conj(IT_0263);
    const ccomplex_t IT_0340 = 4*conj(IT_0193);
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = 4*IT_0274;
    const ccomplex_t IT_0343 = 4*conj(IT_0274);
    const ccomplex_t IT_0344 = 0.125*conj(IT_0274)*IT_0318 + 0.125*IT_0274
      *IT_0319 + conj(IT_0236)*IT_0326 + IT_0236*IT_0329 + IT_0299*IT_0332 +
       IT_0303*IT_0335 + 0.125*conj(IT_0042)*IT_0338 + 0.125*IT_0042*IT_0341 +
       0.125*IT_0291*IT_0342 + 0.125*IT_0295*IT_0343;
    const ccomplex_t IT_0345 = m_mu*s_23*m_N_1;
    const ccomplex_t IT_0346 = 4*IT_0160;
    const ccomplex_t IT_0347 = 0.125*IT_0346;
    const ccomplex_t IT_0348 = IT_0311 + IT_0347;
    const ccomplex_t IT_0349 = 4*conj(IT_0160);
    const ccomplex_t IT_0350 = 0.125*IT_0349;
    const ccomplex_t IT_0351 = IT_0307 + IT_0350;
    const ccomplex_t IT_0352 = 0.125*IT_0337;
    const ccomplex_t IT_0353 = IT_0263 + IT_0352;
    const ccomplex_t IT_0354 = 0.125*IT_0340;
    const ccomplex_t IT_0355 = conj(IT_0263) + IT_0354;
    const ccomplex_t IT_0356 = 2*IT_0193;
    const ccomplex_t IT_0357 = IT_0330 + IT_0356;
    const ccomplex_t IT_0358 = 2*conj(IT_0193);
    const ccomplex_t IT_0359 = IT_0333 + IT_0358;
    const ccomplex_t IT_0360 = 4*IT_0236;
    const ccomplex_t IT_0361 = 4*conj(IT_0236);
    const ccomplex_t IT_0362 = conj(IT_0160)*IT_0236 + IT_0160*conj(IT_0236) +
       4*conj(IT_0274)*IT_0348 + 4*IT_0274*IT_0351 + 4*IT_0283*IT_0353 + 4
      *IT_0287*IT_0355 + 0.5*conj(IT_0169)*IT_0357 + 0.5*IT_0169*IT_0359 + 0.5
      *IT_0307*IT_0360 + 0.5*IT_0311*IT_0361;
    const ccomplex_t IT_0363 = s_34 + IT_0026;
    const ccomplex_t IT_0364 = s_12*IT_0363;
    const ccomplex_t IT_0365 = (-8)*IT_0224;
    const ccomplex_t IT_0366 = (-8)*conj(IT_0224);
    const ccomplex_t IT_0367 = (-4)*IT_0003;
    const ccomplex_t IT_0368 = -s_34;
    const ccomplex_t IT_0369 = IT_0026 + IT_0368;
    const ccomplex_t IT_0370 = s_12*IT_0369;
    const ccomplex_t IT_0371 = (-8)*IT_0077;
    const ccomplex_t IT_0372 = (-4)*IT_0120;
    const ccomplex_t IT_0373 = m_mu*s_14*m_N_1;
    const ccomplex_t IT_0374 = 4*conj(IT_0169);
    const ccomplex_t IT_0375 = 0.125*IT_0374;
    const ccomplex_t IT_0376 = IT_0283 + IT_0375;
    const ccomplex_t IT_0377 = 4*IT_0169;
    const ccomplex_t IT_0378 = 0.125*IT_0377;
    const ccomplex_t IT_0379 = 8*IT_0311;
    const ccomplex_t IT_0380 = 8*IT_0307;
    const ccomplex_t IT_0381 = 0.125*conj(IT_0236)*IT_0320 + 0.125*IT_0236
      *IT_0321 + 0.125*conj(IT_0160)*IT_0357 + 0.125*IT_0160*IT_0359 + 0.125
      *IT_0283*IT_0360 + 0.125*IT_0287*IT_0361 + IT_0274*IT_0376 + conj(IT_0274)
      *(IT_0287 + IT_0378) + 0.125*IT_0355*IT_0379 + 0.125*IT_0353*IT_0380;
    const ccomplex_t IT_0382 = s_34*IT_0027;
    const ccomplex_t IT_0383 = 8*IT_0149;
    const ccomplex_t IT_0384 = 0.125*IT_0383;
    const ccomplex_t IT_0385 = IT_0178 + IT_0378 + IT_0384;
    const ccomplex_t IT_0386 = 8*conj(IT_0149);
    const ccomplex_t IT_0387 = 0.125*IT_0386;
    const ccomplex_t IT_0388 = conj(IT_0178) + IT_0375 + IT_0387;
    const ccomplex_t IT_0389 = 4*IT_0042;
    const ccomplex_t IT_0390 = 4*conj(IT_0042);
    const ccomplex_t IT_0391 = 8*IT_0141;
    const ccomplex_t IT_0392 = 8*conj(IT_0141);
    const ccomplex_t IT_0393 = 8*IT_0295;
    const ccomplex_t IT_0394 = 8*IT_0291;
    const ccomplex_t IT_0395 = 0.125*conj(IT_0160)*IT_0318 + 0.125*IT_0160
      *IT_0319 + 0.125*conj(IT_0042)*IT_0320 + 0.125*IT_0042*IT_0321 + 0.125
      *IT_0307*IT_0324 + 0.125*IT_0311*IT_0327 + 0.125*conj(IT_0141)*(IT_0377 +
       IT_0383) + conj(IT_0176)*IT_0385 + 0.125*IT_0141*(IT_0374 + IT_0386) +
       IT_0176*IT_0388 + 0.125*IT_0283*IT_0389 + 0.125*IT_0287*IT_0390 + 0.125
      *conj(IT_0178)*IT_0391 + 0.125*IT_0178*IT_0392 + 0.125*IT_0351*IT_0393 +
       0.125*IT_0348*IT_0394;
    const ccomplex_t IT_0396 = 16*IT_0149;
    const ccomplex_t IT_0397 = 0.125*IT_0396;
    const ccomplex_t IT_0398 = IT_0169 + IT_0397;
    const ccomplex_t IT_0399 = 16*IT_0146;
    const ccomplex_t IT_0400 = 0.125*IT_0399;
    const ccomplex_t IT_0401 = IT_0160 + IT_0400;
    const ccomplex_t IT_0402 = 16*conj(IT_0146);
    const ccomplex_t IT_0403 = 16*IT_0177;
    const ccomplex_t IT_0404 = 0.125*conj(IT_0274);
    const ccomplex_t IT_0405 = 16*conj(IT_0149);
    const ccomplex_t IT_0406 = 0.125*IT_0178;
    const ccomplex_t IT_0407 = 0.125*conj(IT_0263);
    const ccomplex_t IT_0408 = conj(IT_0160)*IT_0348 + conj(IT_0169)*(IT_0287 
      + IT_0378) + conj(IT_0178)*(IT_0169 + 2*IT_0178 + IT_0397) + conj(IT_0149)
      *IT_0398 + conj(IT_0146)*IT_0401 + 0.125*IT_0177*IT_0402 + conj(IT_0177)*
      (IT_0401 + 0.125*IT_0403) + IT_0342*IT_0404 + IT_0405*IT_0406 + IT_0330
      *IT_0407;
    const ccomplex_t IT_0409 = IT_0027*IT_0369;
    const ccomplex_t IT_0410 = 2*conj(IT_0113);
    const ccomplex_t IT_0411 = (-2)*conj(IT_0112);
    const ccomplex_t IT_0412 = IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = (-2)*conj(IT_0113);
    const ccomplex_t IT_0414 = 2*conj(IT_0077);
    const ccomplex_t IT_0415 = (-2)*conj(IT_0076);
    const ccomplex_t IT_0416 = IT_0414 + IT_0415;
    const ccomplex_t IT_0417 = 0.125*IT_0274;
    const ccomplex_t IT_0418 = 2*IT_0077;
    const ccomplex_t IT_0419 = 2*IT_0112;
    const ccomplex_t IT_0420 = (-2)*IT_0113;
    const ccomplex_t IT_0421 = IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = 0.125*conj(IT_0193);
    const ccomplex_t IT_0423 = 2*IT_0113;
    const ccomplex_t IT_0424 = (-2)*IT_0112;
    const ccomplex_t IT_0425 = IT_0423 + IT_0424;
    const ccomplex_t IT_0426 = 2*conj(IT_0076);
    const ccomplex_t IT_0427 = (-2)*conj(IT_0077);
    const ccomplex_t IT_0428 = IT_0426 + IT_0427;
    const ccomplex_t IT_0429 = 0.125*IT_0236;
    const ccomplex_t IT_0430 = (-2)*IT_0077;
    const ccomplex_t IT_0431 = 0.125*conj(IT_0236);
    const ccomplex_t IT_0432 = IT_0077*conj(IT_0113) + IT_0081*(IT_0112 +
       IT_0114) + IT_0076*(conj(IT_0112) + IT_0116 + 0.25*conj(IT_0236) + (-0.25
      )*conj(IT_0274)) + 0.25*conj(IT_0112)*(IT_0193 + 0.5*IT_0371) + 0.125
      *IT_0263*IT_0412 + 0.125*IT_0193*IT_0413 + IT_0416*IT_0417 + IT_0404
      *IT_0418 + IT_0421*IT_0422 + IT_0407*IT_0425 + IT_0428*IT_0429 + IT_0430
      *IT_0431;
    const ccomplex_t IT_0433 = m_mu*s_24*m_N_1;
    const ccomplex_t IT_0434 = 2*conj(IT_0042);
    const ccomplex_t IT_0435 = 2*IT_0042;
    const ccomplex_t IT_0436 = 0.125*IT_0389;
    const ccomplex_t IT_0437 = 0.125*IT_0390;
    const ccomplex_t IT_0438 = IT_0291*IT_0332 + IT_0295*IT_0335 + 0.125*conj
      (IT_0094)*IT_0338 + 0.125*IT_0094*IT_0341 + 0.125*IT_0299*IT_0342 + 0.125
      *IT_0303*IT_0343 + IT_0417*IT_0434 + IT_0404*IT_0435 + conj(IT_0236)*
      (IT_0303 + IT_0436) + IT_0236*(IT_0299 + IT_0437);
    const ccomplex_t IT_0439 = conj(IT_0169)*IT_0318 + IT_0169*IT_0319 +
       IT_0283*IT_0324 + IT_0287*IT_0327 + IT_0299*IT_0346 + IT_0303*IT_0349 +
       IT_0388*IT_0393 + IT_0385*IT_0394 + IT_0160*IT_0434 + conj(IT_0160)
      *IT_0435 + IT_0380*(IT_0176 + 0.125*IT_0391 + IT_0436) + IT_0379*(conj
      (IT_0176) + 0.125*IT_0392 + IT_0437);
    const ccomplex_t IT_0440 = IT_0026*IT_0027;
    const ccomplex_t IT_0441 = 16*conj(IT_0151);
    const ccomplex_t IT_0442 = 16*IT_0141;
    const ccomplex_t IT_0443 = 16*IT_0176;
    const ccomplex_t IT_0444 = 16*conj(IT_0141);
    const ccomplex_t IT_0445 = 0.125*IT_0179;
    const ccomplex_t IT_0446 = 16*IT_0151;
    const ccomplex_t IT_0447 = 0.125*IT_0446;
    const ccomplex_t IT_0448 = 2*IT_0179;
    const ccomplex_t IT_0449 = conj(IT_0042)*IT_0326 + IT_0042*IT_0329 +
       IT_0160*IT_0376 + conj(IT_0160)*(IT_0287 + IT_0378) + conj(IT_0177)*
      (IT_0169 + 2*IT_0178 + IT_0397) + conj(IT_0146)*IT_0398 + 0.125*conj
      (IT_0178)*(IT_0399 + IT_0403) + IT_0357*IT_0404 + IT_0311*(conj(IT_0169) +
       0.125*IT_0405) + IT_0402*IT_0406 + IT_0359*IT_0417 + IT_0341*IT_0429 +
       IT_0338*IT_0431 + IT_0303*(conj(IT_0094) + 0.125*IT_0441) + 0.125*conj
      (IT_0179)*(IT_0442 + IT_0443) + IT_0444*IT_0445 + conj(IT_0141)*(IT_0094 +
       IT_0447) + conj(IT_0176)*(IT_0094 + IT_0447 + IT_0448);
    const ccomplex_t IT_0450 = 0.125*IT_0442;
    const ccomplex_t IT_0451 = IT_0042 + IT_0450;
    const ccomplex_t IT_0452 = conj(IT_0094)*IT_0326 + IT_0337*IT_0422 +
       IT_0360*IT_0431 + conj(IT_0042)*(IT_0303 + IT_0436) + 0.125*IT_0176
      *IT_0444 + IT_0441*IT_0445 + conj(IT_0151)*(IT_0094 + IT_0447) + conj
      (IT_0179)*(IT_0094 + IT_0447 + IT_0448) + conj(IT_0141)*IT_0451 + conj
      (IT_0176)*(0.125*IT_0443 + IT_0451);
    const ccomplex_t IT_0453 = -IT_0183;
    const ccomplex_t IT_0454 = 0.5*IT_0263;
    const ccomplex_t IT_0455 = IT_0182 + IT_0184 + IT_0453;
    const ccomplex_t IT_0456 = (-0.5)*IT_0236;
    const ccomplex_t IT_0457 = (-0.5)*conj(IT_0236);
    const ccomplex_t IT_0458 = IT_0027*IT_0363;
    const ccomplex_t IT_0459 = IT_0269*IT_0274 + IT_0266*conj(IT_0274) + (-2)
      *conj(IT_0251)*(IT_0193 + 4*IT_0222 + (-4)*IT_0224 + (-0.5)*IT_0336) + (-2
      )*IT_0251*(conj(IT_0193) + 4*conj(IT_0222) + (-4)*conj(IT_0224) + (-0.5)
      *IT_0339) + 8*conj(IT_0249)*(IT_0222 + -1./4*IT_0263 + 0.125*IT_0356 +
       0.125*IT_0365) + 8*IT_0249*(conj(IT_0222) + -1./4*conj(IT_0263) + 0.125
      *IT_0358 + 0.125*IT_0366) + (-2)*IT_0229*IT_0456 + (-2)*IT_0226*IT_0457;
    const ccomplex_t IT_0460 = IT_0004*(conj(IT_0042)*IT_0079 + IT_0042
      *IT_0081 + conj(IT_0094)*IT_0115 + IT_0094*IT_0117) + (IT_0081*IT_0094 +
       IT_0079*conj(IT_0094) + conj(IT_0042)*IT_0115 + IT_0042*IT_0117)*IT_0121 
      + IT_0122*(IT_0081*(IT_0141 + IT_0146) + IT_0079*(conj(IT_0141) + conj
      (IT_0146)) + IT_0115*(conj(IT_0149) + conj(IT_0151)) + IT_0117*IT_0152) +
       IT_0153*(conj(IT_0160)*IT_0162 + IT_0160*IT_0164 + conj(IT_0169)*IT_0171 
      + IT_0169*IT_0173) + (IT_0164*IT_0169 + IT_0162*conj(IT_0169) + conj
      (IT_0160)*IT_0171 + IT_0160*IT_0173)*IT_0174 + IT_0175*(IT_0164*(IT_0176 +
       IT_0177) + IT_0162*(conj(IT_0176) + conj(IT_0177)) + IT_0171*(conj
      (IT_0178) + conj(IT_0179)) + IT_0173*IT_0180) + IT_0186*(conj(IT_0193)
      *IT_0226 + IT_0193*IT_0229 + conj(IT_0236)*IT_0253 + IT_0236*IT_0256) +
       IT_0258*(conj(IT_0263)*IT_0266 + IT_0263*IT_0269 + conj(IT_0274)*IT_0277 
      + IT_0274*IT_0280) + IT_0282*IT_0315 + IT_0317*IT_0322 + 8*IT_0323*IT_0344
       + 2*IT_0345*IT_0362 + 8*IT_0364*(IT_0224*conj(IT_0224) + (IT_0249 + 
      -IT_0251)*(conj(IT_0249) + -conj(IT_0251)) + 0.125*conj(IT_0222)*IT_0365 +
       IT_0222*(conj(IT_0222) + 0.125*IT_0366)) + (IT_0173*(IT_0176 + IT_0177) +
       IT_0171*(conj(IT_0176) + conj(IT_0177)) + IT_0162*(conj(IT_0178) + conj
      (IT_0179)) + (conj(IT_0113) + IT_0163)*IT_0180)*IT_0367 + 8*IT_0370*(conj
      (IT_0076)*IT_0077 + -IT_0076*IT_0081 + -(IT_0112 + IT_0114)*(conj(IT_0112)
       + IT_0116) + 0.125*conj(IT_0077)*IT_0371) + (IT_0117*(IT_0141 + IT_0146) 
      + IT_0115*(conj(IT_0141) + conj(IT_0146)) + IT_0079*(conj(IT_0149) + conj
      (IT_0151)) + (conj(IT_0076) + IT_0080)*IT_0152)*IT_0372 + 8*IT_0373
      *IT_0381 + 8*IT_0382*IT_0395 + 8*IT_0184*IT_0408 + 8*IT_0409*IT_0432 + 8
      *IT_0433*IT_0438 + IT_0181*IT_0439 + 8*IT_0440*IT_0449 + 8*IT_0182*IT_0452
       + 2*(IT_0181 + IT_0182 + IT_0185 + IT_0453)*(IT_0274*(conj(IT_0112) + 0.5
      *IT_0413) + 1./2*conj(IT_0274)*IT_0421 + conj(IT_0263)*(IT_0076 + 0.5
      *IT_0430) + IT_0428*IT_0454) + 2*(conj(IT_0236)*IT_0274 + IT_0236*conj
      (IT_0274) + 0.5*conj(IT_0263)*IT_0356 + IT_0358*IT_0454)*IT_0455 + 
      (IT_0181 + IT_0184 + IT_0257 + IT_0453)*(IT_0193*IT_0416 + (-2)*conj
      (IT_0193)*(IT_0076 + (-0.5)*IT_0418) + (-2)*IT_0412*IT_0456 + (-2)*IT_0425
      *IT_0457) + IT_0458*IT_0459;
    return create_ccomplex_return(IT_0460);
}

