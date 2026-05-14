#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_tau_to_anti_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_tau_to_anti_tau_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = m_C_2*m_tau;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_C_2, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0011;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0011*IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_tau, 2);
    const ccomplex_t IT_0020 = cpow(s_13 + (-0.5)*IT_0004 + (-0.5)*IT_0019 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0013*IT_0022;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0011;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = IT_0020*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0013*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = 2*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0011;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = IT_0009*IT_0010;
    const ccomplex_t IT_0035 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0008*IT_0011;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0036 + -IT_0038 + (-2)
      *IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0033*IT_0042;
    const ccomplex_t IT_0044 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0045 = cpow(s_23 + (-0.5)*IT_0019 + (-0.5)*IT_0044 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0049 = IT_0034*IT_0048;
    const ccomplex_t IT_0050 = IT_0037*IT_0048;
    const ccomplex_t IT_0051 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0052 = IT_0037*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + -IT_0050 + (-2)
      *IT_0052);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0026*IT_0054;
    const ccomplex_t IT_0056 = cpow(s_23 + (-0.5)*IT_0019 + (-0.5)*IT_0044 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0047 + IT_0058;
    const ccomplex_t IT_0060 = (-2)*IT_0059;
    const ccomplex_t IT_0061 = IT_0031 + IT_0060;
    const ccomplex_t IT_0062 = 2*conj(IT_0030);
    const ccomplex_t IT_0063 = (-2)*conj(IT_0059);
    const ccomplex_t IT_0064 = IT_0062 + IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0011*IT_0014*IT_0016;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0068 = IT_0034*IT_0067;
    const ccomplex_t IT_0069 = IT_0037*IT_0067;
    const ccomplex_t IT_0070 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0071 = IT_0037*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0068 + -IT_0069 + (-2)
      *IT_0071);
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0066*IT_0073;
    const ccomplex_t IT_0075 = IT_0045*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0078 = IT_0034*IT_0077;
    const ccomplex_t IT_0079 = IT_0037*IT_0077;
    const ccomplex_t IT_0080 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0081 = IT_0037*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0078 + -IT_0079 + (-2)
      *IT_0081);
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0018*IT_0083;
    const ccomplex_t IT_0085 = IT_0056*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0076 + IT_0086;
    const ccomplex_t IT_0088 = 2*IT_0059;
    const ccomplex_t IT_0089 = (-2)*IT_0030;
    const ccomplex_t IT_0090 = IT_0088 + IT_0089;
    const ccomplex_t IT_0091 = 2*conj(IT_0059);
    const ccomplex_t IT_0092 = (-2)*conj(IT_0030);
    const ccomplex_t IT_0093 = IT_0091 + IT_0092;
    const ccomplex_t IT_0094 = s_13*s_14;
    const ccomplex_t IT_0095 = s_34*IT_0004;
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = IT_0094 + IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0076 + (-0.5)*IT_0086;
    const ccomplex_t IT_0099 = cpow(s_12 + 0.5*IT_0004 + (-0.5)*IT_0019 + 0.5
      *IT_0044 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0010;
    const ccomplex_t IT_0101 = IT_0018*IT_0100;
    const ccomplex_t IT_0102 = IT_0099*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = 2*IT_0087;
    const ccomplex_t IT_0106 = (-2)*IT_0024;
    const ccomplex_t IT_0107 = IT_0105 + IT_0106;
    const ccomplex_t IT_0108 = 2*conj(IT_0087);
    const ccomplex_t IT_0109 = (-2)*conj(IT_0024);
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = e_em*IT_0034;
    const ccomplex_t IT_0112 = e_em*IT_0037;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0111 + -IT_0112);
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0026*IT_0114;
    const ccomplex_t IT_0116 = IT_0099*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = (-0.5)*IT_0047 + (-0.5)*IT_0058;
    const ccomplex_t IT_0120 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = s_13 + IT_0121;
    const ccomplex_t IT_0123 = s_13*IT_0122;
    const ccomplex_t IT_0124 = 0.5*IT_0023 + IT_0103;
    const ccomplex_t IT_0125 = (-0.5)*IT_0029 + -IT_0117;
    const ccomplex_t IT_0126 = pow(s_34, 2);
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = -IT_0019;
    const ccomplex_t IT_0129 = IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = IT_0000*IT_0129;
    const ccomplex_t IT_0131 = 2*IT_0125;
    const ccomplex_t IT_0132 = 2*conj(IT_0125);
    const ccomplex_t IT_0133 = pow(m_Z, 2);
    const ccomplex_t IT_0134 = s_14*IT_0133;
    const ccomplex_t IT_0135 = IT_0001*IT_0134;
    const ccomplex_t IT_0136 = 0.5*IT_0135;
    const ccomplex_t IT_0137 = s_14 + IT_0136;
    const ccomplex_t IT_0138 = m_tau*IT_0137;
    const ccomplex_t IT_0139 = 2*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0141 = IT_0018*IT_0054;
    const ccomplex_t IT_0142 = IT_0056*IT_0140*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0144 = IT_0042*IT_0066;
    const ccomplex_t IT_0145 = IT_0045*IT_0143*IT_0144;
    const ccomplex_t IT_0146 = m_C_2*IT_0076;
    const ccomplex_t IT_0147 = m_C_2*IT_0086;
    const ccomplex_t IT_0148 = (-0.5)*IT_0142 + (-0.5)*IT_0145 + (-0.5)
      *IT_0146 + (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = 2*conj(IT_0148);
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0151 = IT_0026*IT_0100;
    const ccomplex_t IT_0152 = IT_0099*IT_0150*IT_0151;
    const ccomplex_t IT_0153 = m_tau*IT_0117;
    const ccomplex_t IT_0154 = 0.5*IT_0152 + (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = 2*conj(IT_0154);
    const ccomplex_t IT_0156 = IT_0149 + IT_0155;
    const ccomplex_t IT_0157 = s_34*IT_0133;
    const ccomplex_t IT_0158 = IT_0001*IT_0157;
    const ccomplex_t IT_0159 = (-0.25)*IT_0158;
    const ccomplex_t IT_0160 = s_34 + IT_0159;
    const ccomplex_t IT_0161 = m_C_2*IT_0160;
    const ccomplex_t IT_0162 = (-4)*IT_0161;
    const ccomplex_t IT_0163 = IT_0018*IT_0114;
    const ccomplex_t IT_0164 = IT_0099*IT_0150*IT_0163;
    const ccomplex_t IT_0165 = m_tau*IT_0103;
    const ccomplex_t IT_0166 = (-0.5)*IT_0164 + 0.5*IT_0165;
    const ccomplex_t IT_0167 = 2*conj(IT_0166);
    const ccomplex_t IT_0168 = IT_0033*IT_0073;
    const ccomplex_t IT_0169 = IT_0045*IT_0143*IT_0168;
    const ccomplex_t IT_0170 = m_C_2*IT_0058;
    const ccomplex_t IT_0171 = IT_0026*IT_0083;
    const ccomplex_t IT_0172 = IT_0056*IT_0140*IT_0171;
    const ccomplex_t IT_0173 = m_C_2*IT_0047;
    const ccomplex_t IT_0174 = (-0.5)*IT_0169 + (-0.5)*IT_0170 + (-0.5)
      *IT_0172 + (-0.5)*IT_0173;
    const ccomplex_t IT_0175 = 2*conj(IT_0174);
    const ccomplex_t IT_0176 = IT_0167 + IT_0175;
    const ccomplex_t IT_0177 = s_13*s_34;
    const ccomplex_t IT_0178 = s_14*IT_0019;
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = IT_0177 + IT_0179;
    const ccomplex_t IT_0181 = (-0.25)*IT_0135;
    const ccomplex_t IT_0182 = s_14 + IT_0181;
    const ccomplex_t IT_0183 = m_tau*IT_0182;
    const ccomplex_t IT_0184 = (-4)*IT_0183;
    const ccomplex_t IT_0185 = (-2)*conj(IT_0148);
    const ccomplex_t IT_0186 = (-2)*conj(IT_0154);
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0158;
    const ccomplex_t IT_0189 = s_34 + IT_0188;
    const ccomplex_t IT_0190 = m_C_2*IT_0189;
    const ccomplex_t IT_0191 = 2*IT_0190;
    const ccomplex_t IT_0192 = (-2)*conj(IT_0166);
    const ccomplex_t IT_0193 = (-2)*conj(IT_0174);
    const ccomplex_t IT_0194 = IT_0192 + IT_0193;
    const ccomplex_t IT_0195 = s_13*IT_0129;
    const ccomplex_t IT_0196 = 2*IT_0124;
    const ccomplex_t IT_0197 = m_tau*IT_0122;
    const ccomplex_t IT_0198 = conj(IT_0148) + conj(IT_0154);
    const ccomplex_t IT_0199 = IT_0148 + IT_0154;
    const ccomplex_t IT_0200 = 2*conj(IT_0124);
    const ccomplex_t IT_0201 = 2*IT_0166;
    const ccomplex_t IT_0202 = s_13*IT_0133;
    const ccomplex_t IT_0203 = s_14*s_34;
    const ccomplex_t IT_0204 = (-2)*IT_0203;
    const ccomplex_t IT_0205 = IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = IT_0001*IT_0205;
    const ccomplex_t IT_0207 = (-0.5)*IT_0206;
    const ccomplex_t IT_0208 = s_13 + IT_0207;
    const ccomplex_t IT_0209 = conj(IT_0166) + conj(IT_0174);
    const ccomplex_t IT_0210 = 4*IT_0166;
    const ccomplex_t IT_0211 = 4*IT_0154;
    const ccomplex_t IT_0212 = IT_0166 + IT_0174;
    const ccomplex_t IT_0213 = (-6)*conj(IT_0166);
    const ccomplex_t IT_0214 = m_C_2*IT_0122;
    const ccomplex_t IT_0215 = 2*IT_0024;
    const ccomplex_t IT_0216 = (-2)*IT_0087;
    const ccomplex_t IT_0217 = 2*conj(IT_0024);
    const ccomplex_t IT_0218 = (-2)*conj(IT_0087);
    const ccomplex_t IT_0219 = s_13*IT_0006;
    const ccomplex_t IT_0220 = IT_0000*IT_0122;
    const ccomplex_t IT_0221 = (-2)*IT_0125;
    const ccomplex_t IT_0222 = (-2)*conj(IT_0125);
    const ccomplex_t IT_0223 = m_C_2*IT_0129;
    const ccomplex_t IT_0224 = (-2)*IT_0124;
    const ccomplex_t IT_0225 = (-2)*conj(IT_0124);
    const ccomplex_t IT_0226 = m_tau*IT_0006;
    const ccomplex_t IT_0227 = (-2)*IT_0166;
    const ccomplex_t IT_0228 = IT_0180*IT_0200;
    const ccomplex_t IT_0229 = IT_0000*IT_0133;
    const ccomplex_t IT_0230 = (-3)*IT_0229;
    const ccomplex_t IT_0231 = 3*IT_0229;
    const ccomplex_t IT_0232 = pow(m_Z, 4);
    const ccomplex_t IT_0233 = s_13*IT_0232;
    const ccomplex_t IT_0234 = IT_0001*IT_0233;
    const ccomplex_t IT_0235 = (-0.5)*IT_0234;
    const ccomplex_t IT_0236 = IT_0202 + IT_0204 + IT_0235;
    const ccomplex_t IT_0237 = (-0.25)*IT_0234;
    const ccomplex_t IT_0238 = IT_0203 + IT_0237;
    const ccomplex_t IT_0239 = (-2)*IT_0174;
    const ccomplex_t IT_0240 = IT_0227 + IT_0239;
    const ccomplex_t IT_0241 = (-2)*IT_0148;
    const ccomplex_t IT_0242 = (-2)*IT_0154;
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = 0.5*IT_0184;
    const ccomplex_t IT_0245 = 0.5*IT_0180;
    const ccomplex_t IT_0246 = 2*IT_0174;
    const ccomplex_t IT_0247 = IT_0201 + IT_0246;
    const ccomplex_t IT_0248 = 2*IT_0148;
    const ccomplex_t IT_0249 = 2*IT_0154;
    const ccomplex_t IT_0250 = IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = 0.5*IT_0139;
    const ccomplex_t IT_0252 = IT_0196*IT_0245;
    const ccomplex_t IT_0253 = IT_0007*(conj(IT_0024)*IT_0061 + IT_0024
      *IT_0064 + conj(IT_0087)*IT_0090 + IT_0087*IT_0093) + IT_0097*((conj
      (IT_0098) + conj(IT_0104))*IT_0107 + (IT_0098 + IT_0104)*IT_0110 + IT_0093
      *(IT_0118 + IT_0119) + IT_0090*(conj(IT_0118) + conj(IT_0119))) + IT_0123*
      (IT_0110*IT_0124 + IT_0107*conj(IT_0124) + IT_0093*IT_0125 + IT_0090*conj
      (IT_0125)) + IT_0130*(conj(IT_0124)*IT_0131 + IT_0124*IT_0132) + IT_0119*
      (IT_0139*IT_0156 + IT_0162*IT_0176 + IT_0132*IT_0180) + IT_0118*(IT_0132
      *IT_0180 + IT_0184*IT_0187 + IT_0191*IT_0194) + IT_0195*(conj(IT_0125)
      *IT_0131 + conj(IT_0124)*IT_0196) + IT_0197*(IT_0124*IT_0167 + conj
      (IT_0174)*IT_0196 + IT_0131*IT_0198 + IT_0132*IT_0199 + IT_0174*IT_0200 +
       conj(IT_0124)*IT_0201) + 4*IT_0208*(IT_0209*(IT_0174 + 0.25*IT_0210) +
       IT_0198*(IT_0148 + 0.25*IT_0211)) + (-6)*IT_0000*(IT_0198*IT_0212 +
       IT_0199*(conj(IT_0174) + (-0.166666666666667)*IT_0213)) + IT_0214*
      (IT_0061*IT_0209 + IT_0064*IT_0212 + IT_0198*(IT_0215 + IT_0216) + IT_0199
      *(IT_0217 + IT_0218)) + (conj(IT_0030)*IT_0031 + conj(IT_0059)*IT_0090 +
       IT_0059*IT_0092 + conj(IT_0087)*IT_0107 + IT_0087*IT_0109 + conj(IT_0024)
      *IT_0215)*IT_0219 + IT_0220*(IT_0093*IT_0124 + IT_0090*conj(IT_0124) +
       conj(IT_0087)*IT_0131 + IT_0087*IT_0132 + conj(IT_0024)*IT_0221 + IT_0024
      *IT_0222) + IT_0223*(IT_0209*IT_0221 + IT_0212*IT_0222 + IT_0198*IT_0224 +
       IT_0199*IT_0225) + IT_0226*(IT_0087*IT_0167 + IT_0110*IT_0174 + IT_0107
      *conj(IT_0174) + IT_0024*IT_0192 + IT_0090*IT_0198 + IT_0093*IT_0199 +
       conj(IT_0087)*IT_0201 + conj(IT_0024)*IT_0227) + IT_0098*(IT_0156*IT_0162
       + IT_0139*IT_0176 + IT_0228) + IT_0104*(IT_0187*IT_0191 + IT_0184*IT_0194
       + IT_0228) + 2*conj(IT_0118)*(IT_0104*IT_0230 + -IT_0098*IT_0231 + (-2)
      *IT_0118*IT_0236 + 4*IT_0119*IT_0238 + 0.5*IT_0191*IT_0240 + IT_0243
      *IT_0244 + IT_0131*IT_0245) + 2*conj(IT_0119)*(IT_0098*IT_0230 + -IT_0104
      *IT_0231 + (-2)*IT_0119*IT_0236 + 4*IT_0118*IT_0238 + IT_0131*IT_0245 +
       0.5*IT_0162*IT_0247 + IT_0250*IT_0251) + 2*conj(IT_0104)*(IT_0118*IT_0230
       + -IT_0119*IT_0231 + (-2)*IT_0104*IT_0236 + 4*IT_0098*IT_0238 + 0.5
      *IT_0191*IT_0243 + IT_0240*IT_0244 + IT_0252) + 2*conj(IT_0098)*(IT_0119
      *IT_0230 + -IT_0118*IT_0231 + (-2)*IT_0098*IT_0236 + 4*IT_0104*IT_0238 +
       0.5*IT_0162*IT_0250 + IT_0247*IT_0251 + IT_0252);
    return create_ccomplex_return(IT_0253);
}

