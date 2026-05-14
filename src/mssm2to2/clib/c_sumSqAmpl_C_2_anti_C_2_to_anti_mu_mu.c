#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_anti_mu_mu(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = e_em*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + -IT_0007);
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0010;
    const ccomplex_t IT_0013 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + -IT_0012 + (-2)
      *IT_0014);
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_2, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = s_14*s_23;
    const ccomplex_t IT_0023 = m_mu*m_C_2;
    const ccomplex_t IT_0024 = s_13 + s_14;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + (-2)*IT_0018 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = cos(beta);
    const ccomplex_t IT_0029 = cpow(IT_0028, -1);
    const ccomplex_t IT_0030 = sin(beta);
    const ccomplex_t IT_0031 = e_em*m_mu*IT_0005*IT_0027*IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = IT_0005*(IT_0034 + IT_0036);
    const ccomplex_t IT_0038 = 1.4142135623731*e_em*IT_0037;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0032*IT_0039;
    const ccomplex_t IT_0041 = IT_0026*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = pow(m_Z, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + (-2)*IT_0018 + IT_0043 + 
      -reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0028*IT_0035;
    const ccomplex_t IT_0046 = IT_0030*IT_0033;
    const ccomplex_t IT_0047 = IT_0005*(IT_0045 + -IT_0046);
    const ccomplex_t IT_0048 = 1.4142135623731*e_em*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = e_em*m_mu*IT_0005*IT_0027;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = IT_0044*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = -IT_0042 + -IT_0054;
    const ccomplex_t IT_0056 = -s_24;
    const ccomplex_t IT_0057 = s_23 + IT_0056;
    const ccomplex_t IT_0058 = IT_0023*IT_0057;
    const ccomplex_t IT_0059 = 4*IT_0058;
    const ccomplex_t IT_0060 = cos(alpha);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*e_em*m_mu*IT_0005*IT_0027
      *IT_0029*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = U_d2*V_Wp2;
    const ccomplex_t IT_0064 = IT_0060*IT_0063;
    const ccomplex_t IT_0065 = sin(alpha);
    const ccomplex_t IT_0066 = V_u2*U_Wm2;
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = IT_0005*(IT_0064 + IT_0067);
    const ccomplex_t IT_0069 = 1.4142135623731*e_em*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0062*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_12 + (-2)*IT_0018 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*e_em*m_mu*IT_0005*IT_0027
      *IT_0029*IT_0065;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0060*IT_0066;
    const ccomplex_t IT_0079 = IT_0063*IT_0065;
    const ccomplex_t IT_0080 = IT_0005*(IT_0078 + -IT_0079);
    const ccomplex_t IT_0081 = 1.4142135623731*e_em*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0077*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = -IT_0075 + -IT_0087;
    const ccomplex_t IT_0089 = -s_14;
    const ccomplex_t IT_0090 = s_13 + IT_0089;
    const ccomplex_t IT_0091 = IT_0023*IT_0090;
    const ccomplex_t IT_0092 = (-4)*IT_0091;
    const ccomplex_t IT_0093 = IT_0033*IT_0060;
    const ccomplex_t IT_0094 = IT_0035*IT_0065;
    const ccomplex_t IT_0095 = IT_0005*(IT_0093 + -IT_0094);
    const ccomplex_t IT_0096 = 1.4142135623731*e_em*IT_0095;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0077*IT_0098;
    const ccomplex_t IT_0100 = IT_0085*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = IT_0035*IT_0060;
    const ccomplex_t IT_0103 = IT_0033*IT_0065;
    const ccomplex_t IT_0104 = IT_0005*(IT_0102 + IT_0103);
    const ccomplex_t IT_0105 = 1.4142135623731*e_em*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0062*IT_0107;
    const ccomplex_t IT_0109 = IT_0073*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = -IT_0101 + -IT_0110;
    const ccomplex_t IT_0112 = m_mu*s_23*m_C_2;
    const ccomplex_t IT_0113 = pow(m_mu, 2);
    const ccomplex_t IT_0114 = cpow((-2)*s_13 + IT_0018 + IT_0113 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em*m_mu
      *IT_0005*IT_0027*IT_0029;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = IT_0114*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = -IT_0121;
    const ccomplex_t IT_0123 = m_mu*s_14*m_C_2;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0005*IT_0027*IT_0029;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = IT_0114*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0018*IT_0113;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0133 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0134 = IT_0003*IT_0133;
    const ccomplex_t IT_0135 = IT_0006*IT_0133;
    const ccomplex_t IT_0136 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0137 = IT_0006*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + -IT_0135 + (-2)
      *IT_0137);
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = IT_0132*IT_0139;
    const ccomplex_t IT_0141 = IT_0019*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = s_12*IT_0113;
    const ccomplex_t IT_0145 = IT_0116*IT_0125;
    const ccomplex_t IT_0146 = IT_0114*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = 4*conj(IT_0147);
    const ccomplex_t IT_0149 = 16*conj(IT_0021)*IT_0022 + 4*IT_0025*conj
      (IT_0055) + IT_0059*conj(IT_0088) + IT_0092*conj(IT_0111) + 4*IT_0112*conj
      (IT_0122) + 4*IT_0123*conj(IT_0130) + 16*IT_0131*conj(IT_0143) + IT_0144
      *IT_0148;
    const ccomplex_t IT_0150 = IT_0016*IT_0132;
    const ccomplex_t IT_0151 = IT_0019*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = s_23 + s_24;
    const ccomplex_t IT_0154 = IT_0023*IT_0153;
    const ccomplex_t IT_0155 = IT_0028*IT_0066;
    const ccomplex_t IT_0156 = IT_0030*IT_0063;
    const ccomplex_t IT_0157 = IT_0005*(IT_0155 + IT_0156);
    const ccomplex_t IT_0158 = 1.4142135623731*e_em*IT_0157;
    const ccomplex_t IT_0159 = 0.5*IT_0158;
    const ccomplex_t IT_0160 = IT_0032*IT_0159;
    const ccomplex_t IT_0161 = IT_0026*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = IT_0028*IT_0063;
    const ccomplex_t IT_0164 = IT_0030*IT_0066;
    const ccomplex_t IT_0165 = IT_0005*(IT_0163 + -IT_0164);
    const ccomplex_t IT_0166 = 1.4142135623731*e_em*IT_0165;
    const ccomplex_t IT_0167 = (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = IT_0051*IT_0167;
    const ccomplex_t IT_0169 = IT_0044*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = -IT_0162 + -IT_0170;
    const ccomplex_t IT_0172 = s_13*s_24;
    const ccomplex_t IT_0173 = 8*conj(IT_0147);
    const ccomplex_t IT_0174 = 8*conj(IT_0021);
    const ccomplex_t IT_0175 = 8*IT_0021;
    const ccomplex_t IT_0176 = 16*IT_0152;
    const ccomplex_t IT_0177 = IT_0009*IT_0139;
    const ccomplex_t IT_0178 = IT_0019*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0118*IT_0126;
    const ccomplex_t IT_0182 = IT_0114*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = (-8)*IT_0152;
    const ccomplex_t IT_0186 = (-4)*IT_0021;
    const ccomplex_t IT_0187 = s_34*IT_0018;
    const ccomplex_t IT_0188 = 8*IT_0152;
    const ccomplex_t IT_0189 = s_12*s_34;
    const ccomplex_t IT_0190 = -IT_0022;
    const ccomplex_t IT_0191 = IT_0144 + IT_0172 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = 2*IT_0191;
    const ccomplex_t IT_0193 = s_34 + IT_0113;
    const ccomplex_t IT_0194 = IT_0018*IT_0193;
    const ccomplex_t IT_0195 = 2*IT_0194;
    const ccomplex_t IT_0196 = 4*IT_0025;
    const ccomplex_t IT_0197 = s_12*IT_0193;
    const ccomplex_t IT_0198 = 8*IT_0197;
    const ccomplex_t IT_0199 = (-8)*IT_0194;
    const ccomplex_t IT_0200 = 2*conj(IT_0152);
    const ccomplex_t IT_0201 = conj(IT_0147) + IT_0200;
    const ccomplex_t IT_0202 = (-4)*conj(IT_0021);
    const ccomplex_t IT_0203 = (-4)*IT_0025;
    const ccomplex_t IT_0204 = (-2)*IT_0025;
    const ccomplex_t IT_0205 = conj(IT_0021) + conj(IT_0152);
    const ccomplex_t IT_0206 = conj(IT_0143) + conj(IT_0180);
    const ccomplex_t IT_0207 = -IT_0189;
    const ccomplex_t IT_0208 = -IT_0172;
    const ccomplex_t IT_0209 = IT_0022 + IT_0144 + IT_0207 + IT_0208;
    const ccomplex_t IT_0210 = 2*IT_0209;
    const ccomplex_t IT_0211 = -s_34;
    const ccomplex_t IT_0212 = IT_0113 + IT_0211;
    const ccomplex_t IT_0213 = IT_0018*IT_0212;
    const ccomplex_t IT_0214 = 8*IT_0213;
    const ccomplex_t IT_0215 = s_12*IT_0212;
    const ccomplex_t IT_0216 = (-8)*IT_0215;
    const ccomplex_t IT_0217 = (-2)*IT_0213;
    const ccomplex_t IT_0218 = (-2)*IT_0058;
    const ccomplex_t IT_0219 = -conj(IT_0147);
    const ccomplex_t IT_0220 = (-2)*IT_0091;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = 2*IT_0154;
    const ccomplex_t IT_0223 = (-8)*IT_0172;
    const ccomplex_t IT_0224 = (-4)*IT_0144;
    const ccomplex_t IT_0225 = m_mu*s_24*m_C_2;
    const ccomplex_t IT_0226 = (-4)*IT_0225;
    const ccomplex_t IT_0227 = m_mu*s_13*m_C_2;
    const ccomplex_t IT_0228 = (-4)*IT_0227;
    const ccomplex_t IT_0229 = 2*IT_0058;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0172*conj(IT_0184) + -IT_0131*IT_0201 + 0.25
      *IT_0187*IT_0202 + 0.25*conj(IT_0171)*IT_0204 + 0.25*conj(IT_0088)*IT_0221
       + 0.25*conj(IT_0055)*IT_0222 + 0.25*conj(IT_0180)*IT_0223 + 0.25*conj
      (IT_0143)*IT_0224 + 0.25*conj(IT_0130)*IT_0226 + 0.25*conj(IT_0122)
      *IT_0228 + 0.25*conj(IT_0111)*IT_0230;
    const ccomplex_t IT_0232 = 2*IT_0091;
    const ccomplex_t IT_0233 = 8*IT_0144;
    const ccomplex_t IT_0234 = 8*IT_0225;
    const ccomplex_t IT_0235 = 8*IT_0227;
    const ccomplex_t IT_0236 = (-4)*IT_0154;
    const ccomplex_t IT_0237 = IT_0059*conj(IT_0111);
    const ccomplex_t IT_0238 = conj(IT_0088)*IT_0092;
    const ccomplex_t IT_0239 = IT_0131*(16*conj(IT_0152) + IT_0173) + 16
      *IT_0172*conj(IT_0180) + IT_0174*IT_0187 + conj(IT_0171)*IT_0196 + conj
      (IT_0184)*IT_0223 + conj(IT_0143)*IT_0233 + conj(IT_0130)*IT_0234 + conj
      (IT_0122)*IT_0235 + conj(IT_0055)*IT_0236 + IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = 8*conj(IT_0152);
    const ccomplex_t IT_0241 = IT_0148 + IT_0240;
    const ccomplex_t IT_0242 = 4*IT_0154;
    const ccomplex_t IT_0243 = 4*IT_0123;
    const ccomplex_t IT_0244 = 4*IT_0112;
    const ccomplex_t IT_0245 = 16*conj(IT_0021)*IT_0131 + 16*IT_0022*conj
      (IT_0143) + conj(IT_0171)*IT_0203 + conj(IT_0184)*IT_0224 + conj(IT_0180)
      *IT_0233 + IT_0237 + IT_0238 + IT_0187*IT_0241 + conj(IT_0055)*IT_0242 +
       conj(IT_0122)*IT_0243 + conj(IT_0130)*IT_0244;
    const ccomplex_t IT_0246 = (-2)*IT_0191;
    const ccomplex_t IT_0247 = (-2)*IT_0194;
    const ccomplex_t IT_0248 = IT_0055*(conj(IT_0055)*IT_0198 + conj(IT_0171)
      *IT_0199 + 4*IT_0025*(conj(IT_0021) + -1./2*conj(IT_0147) + -1./2*IT_0200)
       + conj(IT_0184)*IT_0222 + conj(IT_0180)*IT_0236 + conj(IT_0143)*IT_0242 +
       conj(IT_0130)*IT_0246 + conj(IT_0122)*IT_0247);
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = IT_0249*IT_0250;
    const ccomplex_t IT_0252 = cpow(s_12 + IT_0018 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0253 = IT_0251*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = 0.5*IT_0254;
    const ccomplex_t IT_0256 = -s_14 + -s_23;
    const ccomplex_t IT_0257 = s_13 + s_24;
    const ccomplex_t IT_0258 = IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = IT_0023*IT_0258;
    const ccomplex_t IT_0260 = (-8)*IT_0259;
    const ccomplex_t IT_0261 = 2*IT_0131;
    const ccomplex_t IT_0262 = 2*IT_0022;
    const ccomplex_t IT_0263 = IT_0144 + IT_0187 + IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0265 = IT_0257 + IT_0264;
    const ccomplex_t IT_0266 = IT_0023*IT_0265;
    const ccomplex_t IT_0267 = 8*IT_0266;
    const ccomplex_t IT_0268 = (-8)*IT_0266;
    const ccomplex_t IT_0269 = -IT_0260;
    const ccomplex_t IT_0270 = 2*IT_0172;
    const ccomplex_t IT_0271 = IT_0144 + IT_0187 + IT_0261 + IT_0270;
    const ccomplex_t IT_0272 = (-2)*conj(IT_0180);
    const ccomplex_t IT_0273 = -(32*IT_0022 + 64*IT_0131 + 32*IT_0144 + 32
      *IT_0172 + 32*IT_0187)*IT_0255 + (conj(IT_0088) + conj(IT_0111))*IT_0260 +
       (-8)*(IT_0021 + IT_0143 + -conj(IT_0143) + (-0.125)*IT_0174)*IT_0263 + 
      (conj(IT_0122) + conj(IT_0130))*IT_0267 + (IT_0122 + IT_0130)*IT_0268 + 
      (IT_0088 + IT_0111)*IT_0269 + (-4)*IT_0271*(IT_0147 + (-0.25)*IT_0148 + 2
      *IT_0180 + -IT_0184 + conj(IT_0184) + (-0.25)*IT_0185 + (-0.25)*IT_0240 +
       IT_0272);
    const ccomplex_t IT_0274 = 0.25*conj(IT_0055);
    const ccomplex_t IT_0275 = 0.25*conj(IT_0088);
    const ccomplex_t IT_0276 = 0.25*conj(IT_0111);
    const ccomplex_t IT_0277 = 0.25*conj(IT_0171);
    const ccomplex_t IT_0278 = 0.25*conj(IT_0184);
    const ccomplex_t IT_0279 = 0.25*conj(IT_0180);
    const ccomplex_t IT_0280 = 0.25*conj(IT_0143);
    const ccomplex_t IT_0281 = conj(IT_0021)*IT_0112 + conj(IT_0130)*IT_0131 +
       conj(IT_0122)*IT_0172 + 0.25*IT_0225*IT_0241 + IT_0247*IT_0274 + IT_0210
      *IT_0275 + IT_0217*IT_0276 + IT_0192*IT_0277 + IT_0228*IT_0278 + IT_0235
      *IT_0279 + IT_0243*IT_0280;
    const ccomplex_t IT_0282 = conj(IT_0021)*IT_0123 + conj(IT_0122)*IT_0131 +
       conj(IT_0130)*IT_0172 + 0.25*IT_0227*IT_0241 + IT_0246*IT_0274 + IT_0217
      *IT_0275 + IT_0210*IT_0276 + IT_0195*IT_0277 + IT_0226*IT_0278 + IT_0234
      *IT_0279 + IT_0244*IT_0280;
    const ccomplex_t IT_0283 = conj(IT_0122)*IT_0225;
    const ccomplex_t IT_0284 = conj(IT_0130)*IT_0227;
    const ccomplex_t IT_0285 = (-0.5)*IT_0025*conj(IT_0055) + conj(IT_0021)
      *IT_0144 + 0.5*IT_0154*conj(IT_0171) + conj(IT_0143)*IT_0187 + IT_0172
      *IT_0201 + (-0.25)*conj(IT_0088)*IT_0218 + (-0.25)*conj(IT_0111)*IT_0232 +
       -IT_0131*(conj(IT_0184) + IT_0272) + IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = IT_0021*IT_0149 + (-4)*IT_0152*(IT_0025*conj
      (IT_0055) + (-0.25)*IT_0059*conj(IT_0088) + (-0.25)*IT_0092*conj(IT_0111) 
      + -IT_0154*conj(IT_0171) + (-0.25)*IT_0172*IT_0173 + (-0.25)*IT_0144
      *IT_0174) + conj(IT_0152)*(IT_0144*IT_0175 + IT_0172*IT_0176) + IT_0131*
      (IT_0176*conj(IT_0180) + conj(IT_0184)*IT_0185) + IT_0154*conj(IT_0171)
      *IT_0186 + IT_0187*(IT_0175*conj(IT_0180) + conj(IT_0184)*IT_0186 + conj
      (IT_0143)*IT_0188) + IT_0171*(conj(IT_0122)*IT_0192 + conj(IT_0130)
      *IT_0195 + conj(IT_0180)*IT_0196 + conj(IT_0171)*IT_0198 + conj(IT_0055)
      *IT_0199 + 2*IT_0154*(IT_0201 + 0.5*IT_0202) + conj(IT_0143)*IT_0203 +
       conj(IT_0184)*IT_0204) + IT_0088*(IT_0059*IT_0205 + IT_0092*IT_0206 +
       conj(IT_0122)*IT_0210 + conj(IT_0111)*IT_0214 + conj(IT_0088)*IT_0216 +
       conj(IT_0130)*IT_0217 + IT_0218*IT_0219 + conj(IT_0184)*IT_0221) + 4
      *IT_0184*IT_0231 + IT_0111*(IT_0092*IT_0205 + IT_0059*IT_0206 + conj
      (IT_0130)*IT_0210 + conj(IT_0088)*IT_0214 + conj(IT_0111)*IT_0216 + conj
      (IT_0122)*IT_0217 + conj(IT_0184)*IT_0230 + IT_0219*IT_0232) + IT_0180
      *IT_0239 + IT_0143*IT_0245 + IT_0248 + IT_0255*IT_0273 + 4*IT_0122*IT_0281
       + 4*IT_0130*IT_0282 + IT_0188*(IT_0283 + IT_0284) + 4*IT_0147*IT_0285;
    return create_ccomplex_return(IT_0286);
}

