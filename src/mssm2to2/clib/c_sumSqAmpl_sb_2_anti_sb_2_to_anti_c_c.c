#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_anti_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_c*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cpow(V_tb, 2);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = m_b*mu_h*IT_0000;
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = -IT_0001*IT_0003;
    const ccomplex_t IT_0014 = IT_0009*IT_0011;
    const ccomplex_t IT_0015 = IT_0013 + IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -2);
    const ccomplex_t IT_0018 = pow(m_b, 2);
    const ccomplex_t IT_0019 = (-2)*U_sb_01*(((-0.5)*m_W*IT_0015*(IT_0005 + 1.
      /3*IT_0004*IT_0017) + IT_0000*IT_0005*IT_0009*IT_0012*IT_0018)*conj
      (U_sb_01) + -1./2*IT_0003*IT_0005*IT_0010*IT_0012*conj(U_sb_11)) + U_sb_11
      *(IT_0003*IT_0005*IT_0010*IT_0012*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0004*IT_0015*IT_0017 + (-3)*IT_0000*IT_0005*IT_0009*IT_0012*IT_0018)
      *conj(U_sb_11));
    const ccomplex_t IT_0020 = (0 + _Complex_I*-1)*(e_em*IT_0019 + (
      -1.4142135623731)*A_b*IT_0008*IT_0009*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0007*IT_0021;
    const ccomplex_t IT_0023 = pow(m_sb_2, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + (-2)*IT_0023 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*m_c*e_em*IT_0000*IT_0002
      *IT_0005*IT_0009;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0003*IT_0011;
    const ccomplex_t IT_0030 = IT_0001*IT_0009;
    const ccomplex_t IT_0031 = IT_0029 + IT_0030;
    const ccomplex_t IT_0032 = 2*U_sb_01*((IT_0000*IT_0003*IT_0005*IT_0012
      *IT_0018 + (-0.5)*m_W*(IT_0005 + 1./3*IT_0004*IT_0017)*IT_0031)*conj
      (U_sb_01) + 1./2*IT_0005*IT_0009*IT_0010*IT_0012*conj(U_sb_11)) + U_sb_11*
      (IT_0005*IT_0009*IT_0010*IT_0012*conj(U_sb_01) + 2*(IT_0000*IT_0003
      *IT_0005*IT_0012*IT_0018 + -1./3*m_W*IT_0004*IT_0017*IT_0031)*conj(U_sb_11
      ));
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(e_em*IT_0032 +
       1.4142135623731*A_b*IT_0003*IT_0008*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0028*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = -IT_0026 + -IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0016, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0005
      *IT_0016 + 1./3*IT_0004*IT_0040)*U_sb_01*conj(U_sb_01) + IT_0004*IT_0040
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0042 = (-0.333333333333333)*IT_0041;
    const ccomplex_t IT_0043 = m_c*IT_0042;
    const ccomplex_t IT_0044 = (-0.666666666666667)*IT_0041;
    const ccomplex_t IT_0045 = m_c*IT_0044;
    const ccomplex_t IT_0046 = IT_0004*IT_0040;
    const ccomplex_t IT_0047 = e_em*IT_0046;
    const ccomplex_t IT_0048 = IT_0005*IT_0016;
    const ccomplex_t IT_0049 = e_em*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0047 + (-3)*IT_0049);
    const ccomplex_t IT_0051 = (-0.166666666666667)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0040;
    const ccomplex_t IT_0056 = (-0.666666666666667)*IT_0055;
    const ccomplex_t IT_0057 = IT_0052*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = (IT_0043 + -IT_0045)*(IT_0054 + -IT_0058);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = 18*s_34;
    const ccomplex_t IT_0062 = pow(m_c, 2);
    const ccomplex_t IT_0063 = (-18)*IT_0062;
    const ccomplex_t IT_0064 = IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = V_cb*e_em*U_Wm2*conj(U_sb_01);
    const ccomplex_t IT_0066 = IT_0005*IT_0065;
    const ccomplex_t IT_0067 = IT_0005*IT_0012;
    const ccomplex_t IT_0068 = m_b*U_d2*V_cb*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + (-0.5)*IT_0070);
    const ccomplex_t IT_0072 = V_cb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0073 = IT_0005*IT_0072;
    const ccomplex_t IT_0074 = m_b*conj(U_d2)*V_cb*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0075 = IT_0067*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + (-0.5)*IT_0076);
    const ccomplex_t IT_0078 = IT_0071*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_23 + IT_0023 + IT_0062 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = V_cb*e_em*U_Wm1*conj(U_sb_01);
    const ccomplex_t IT_0083 = IT_0005*IT_0082;
    const ccomplex_t IT_0084 = m_b*U_d1*V_cb*e_em*IT_0000*conj(U_sb_11);
    const ccomplex_t IT_0085 = IT_0067*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + (-0.5)*IT_0086);
    const ccomplex_t IT_0088 = V_cb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0089 = IT_0005*IT_0088;
    const ccomplex_t IT_0090 = m_b*conj(U_d1)*V_cb*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0091 = IT_0067*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0089 + (-0.5)*IT_0092);
    const ccomplex_t IT_0094 = IT_0087*IT_0093;
    const ccomplex_t IT_0095 = cpow((-2)*s_23 + IT_0023 + IT_0062 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = -IT_0081 + -IT_0097;
    const ccomplex_t IT_0099 = IT_0044*IT_0054;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*conj
      (V_u1)*e_em*IT_0000*IT_0002*IT_0005*conj(U_sb_01);
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u1
      *e_em*IT_0000*IT_0002*IT_0005*U_sb_01;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0102*IT_0104;
    const ccomplex_t IT_0106 = IT_0095*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*conj
      (V_u2)*e_em*IT_0000*IT_0002*IT_0005*conj(U_sb_01);
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0000*IT_0002*IT_0005*U_sb_01;
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = IT_0109*IT_0111;
    const ccomplex_t IT_0113 = IT_0079*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = -IT_0107 + -IT_0114;
    const ccomplex_t IT_0116 = IT_0044*IT_0058;
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = s_34*IT_0023;
    const ccomplex_t IT_0119 = (-6)*IT_0118;
    const ccomplex_t IT_0120 = s_13*s_14;
    const ccomplex_t IT_0121 = 12*IT_0120;
    const ccomplex_t IT_0122 = IT_0119 + IT_0121;
    const ccomplex_t IT_0123 = (-18)*IT_0118;
    const ccomplex_t IT_0124 = 36*IT_0120;
    const ccomplex_t IT_0125 = IT_0123 + IT_0124;
    const ccomplex_t IT_0126 = cpow(s_12 + IT_0023 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0128 = 0.666666666666667*IT_0127;
    const ccomplex_t IT_0129 = cpow(IT_0128, 2);
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = IT_0126*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0133 = (-2)*IT_0132;
    const ccomplex_t IT_0134 = IT_0126*IT_0132;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0133*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0131 + (-0.0833333333333333)*IT_0136;
    const ccomplex_t IT_0138 = cpow(IT_0137, 2);
    const ccomplex_t IT_0139 = (-36)*IT_0118;
    const ccomplex_t IT_0140 = IT_0023*IT_0062;
    const ccomplex_t IT_0141 = (-36)*IT_0140;
    const ccomplex_t IT_0142 = 72*IT_0120;
    const ccomplex_t IT_0143 = IT_0139 + IT_0141 + IT_0142;
    const ccomplex_t IT_0144 = 0.25*IT_0136;
    const ccomplex_t IT_0145 = IT_0077*IT_0109;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0147 = IT_0079*IT_0145*IT_0146;
    const ccomplex_t IT_0148 = IT_0093*IT_0102;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0150 = IT_0095*IT_0148*IT_0149;
    const ccomplex_t IT_0151 = m_c*IT_0107;
    const ccomplex_t IT_0152 = m_c*IT_0114;
    const ccomplex_t IT_0153 = -IT_0147 + -IT_0150 + IT_0151 + IT_0152;
    const ccomplex_t IT_0154 = IT_0087*IT_0104;
    const ccomplex_t IT_0155 = IT_0095*IT_0149*IT_0154;
    const ccomplex_t IT_0156 = IT_0071*IT_0111;
    const ccomplex_t IT_0157 = IT_0079*IT_0146*IT_0156;
    const ccomplex_t IT_0158 = m_c*IT_0081;
    const ccomplex_t IT_0159 = m_c*IT_0097;
    const ccomplex_t IT_0160 = -IT_0155 + -IT_0157 + IT_0158 + IT_0159;
    const ccomplex_t IT_0161 = (-18)*IT_0160;
    const ccomplex_t IT_0162 = (-6)*IT_0060;
    const ccomplex_t IT_0163 = IT_0161 + IT_0162;
    const ccomplex_t IT_0164 = (-18)*conj(IT_0160);
    const ccomplex_t IT_0165 = (-6)*conj(IT_0060);
    const ccomplex_t IT_0166 = IT_0164 + IT_0165;
    const ccomplex_t IT_0167 = (-18)*IT_0060;
    const ccomplex_t IT_0168 = (-6)*IT_0160;
    const ccomplex_t IT_0169 = IT_0167 + IT_0168;
    const ccomplex_t IT_0170 = (-18)*conj(IT_0060);
    const ccomplex_t IT_0171 = (-6)*conj(IT_0160);
    const ccomplex_t IT_0172 = IT_0170 + IT_0171;
    const ccomplex_t IT_0173 = (-6)*IT_0100;
    const ccomplex_t IT_0174 = (-18)*IT_0098 + IT_0173;
    const ccomplex_t IT_0175 = (-18)*IT_0100;
    const ccomplex_t IT_0176 = (-6)*IT_0098;
    const ccomplex_t IT_0177 = IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = (-18)*conj(IT_0100);
    const ccomplex_t IT_0179 = (-6)*conj(IT_0098);
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = 6*s_34;
    const ccomplex_t IT_0182 = -IT_0062;
    const ccomplex_t IT_0183 = 6*IT_0039;
    const ccomplex_t IT_0184 = 6*conj(IT_0039);
    const ccomplex_t IT_0185 = m_b*e_em*mu_h*IT_0000*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0186 = IT_0005*IT_0185;
    const ccomplex_t IT_0187 = m_b*e_em*mu_h*IT_0000*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0188 = IT_0005*IT_0187;
    const ccomplex_t IT_0189 = A_b*IT_0008*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0190 = IT_0001*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = A_b*IT_0008*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0193 = IT_0001*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = (-0.5)*IT_0186 + 0.5*IT_0188 + (-0.5)*IT_0191 +
       0.5*IT_0194;
    const ccomplex_t IT_0196 = m_c*e_em*IT_0000*IT_0002*IT_0005*IT_0011;
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0195*IT_0197;
    const ccomplex_t IT_0199 = cpow((-2)*s_12 + (-2)*IT_0023 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0200 = IT_0198*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = m_c*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = IT_0001*IT_0005*IT_0012;
    const ccomplex_t IT_0205 = IT_0185*IT_0204;
    const ccomplex_t IT_0206 = IT_0187*IT_0204;
    const ccomplex_t IT_0207 = IT_0011*IT_0189;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = IT_0011*IT_0192;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = (-0.5)*IT_0205 + 0.5*IT_0206 + 0.5*IT_0208 + (
      -0.5)*IT_0210;
    const ccomplex_t IT_0212 = IT_0203*IT_0211;
    const ccomplex_t IT_0213 = pow(m_Z, 2);
    const ccomplex_t IT_0214 = cpow((-2)*s_12 + (-2)*IT_0023 + IT_0213 + 
      -reg_prop, -1);
    const ccomplex_t IT_0215 = IT_0212*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = -IT_0201 + -IT_0216;
    const ccomplex_t IT_0218 = s_13 + s_14;
    const ccomplex_t IT_0219 = m_c*IT_0218;
    const ccomplex_t IT_0220 = 18*IT_0062;
    const ccomplex_t IT_0221 = IT_0061 + IT_0220;
    const ccomplex_t IT_0222 = (-18)*s_34 + IT_0063;
    const ccomplex_t IT_0223 = s_34 + IT_0062;
    const ccomplex_t IT_0224 = 6*IT_0153;
    const ccomplex_t IT_0225 = m_c*s_14;
    const ccomplex_t IT_0226 = (-18)*IT_0153;
    const ccomplex_t IT_0227 = (-6)*IT_0153;
    const ccomplex_t IT_0228 = (-6)*IT_0140;
    const ccomplex_t IT_0229 = IT_0119 + IT_0121 + IT_0228;
    const ccomplex_t IT_0230 = conj(IT_0098) + conj(IT_0115);
    const ccomplex_t IT_0231 = conj(IT_0100) + conj(IT_0117);
    const ccomplex_t IT_0232 = IT_0100 + IT_0117;
    const ccomplex_t IT_0233 = -IT_0144;
    const ccomplex_t IT_0234 = IT_0098 + IT_0115;
    const ccomplex_t IT_0235 = -IT_0137;
    const ccomplex_t IT_0236 = (-18)*IT_0140;
    const ccomplex_t IT_0237 = IT_0123 + IT_0124 + IT_0236;
    const ccomplex_t IT_0238 = 36*s_34;
    const ccomplex_t IT_0239 = conj(IT_0098) + 3*conj(IT_0100);
    const ccomplex_t IT_0240 = 6*conj(IT_0153);
    const ccomplex_t IT_0241 = m_c*s_13;
    const ccomplex_t IT_0242 = 6*conj(IT_0060);
    const ccomplex_t IT_0243 = 18*conj(IT_0160);
    const ccomplex_t IT_0244 = 6*conj(IT_0160);
    const ccomplex_t IT_0245 = 18*conj(IT_0060);
    const ccomplex_t IT_0246 = 18*conj(IT_0153);
    const ccomplex_t IT_0247 = 18*conj(IT_0059);
    const ccomplex_t IT_0248 = 6*conj(IT_0115)*(IT_0059 + 3*IT_0153) + (18
      *conj(IT_0098) + 6*conj(IT_0100))*IT_0160 + conj(IT_0117)*(18*IT_0059 +
       IT_0224) + 6*IT_0060*IT_0239 + IT_0098*(IT_0242 + IT_0243) + IT_0100*
      (IT_0244 + IT_0245) + 6*IT_0115*(conj(IT_0059) + 0.166666666666667*IT_0246
      ) + IT_0117*(IT_0240 + IT_0247);
    const ccomplex_t IT_0249 = -s_14;
    const ccomplex_t IT_0250 = s_13 + IT_0249;
    const ccomplex_t IT_0251 = m_c*IT_0250;
    const ccomplex_t IT_0252 = 18*IT_0039;
    const ccomplex_t IT_0253 = 18*conj(IT_0039);
    const ccomplex_t IT_0254 = IT_0183*IT_0230 + IT_0184*IT_0234 + IT_0144*((
      -12)*IT_0039 + 12*conj(IT_0039) + (-6)*IT_0059 + 6*conj(IT_0059) + IT_0161
       + IT_0162 + IT_0226 + IT_0242 + IT_0243 + IT_0246) + IT_0137*((-36)
      *IT_0039 + 36*conj(IT_0039) + (-18)*IT_0059 + IT_0167 + IT_0168 + IT_0227 
      + IT_0240 + IT_0244 + IT_0245 + IT_0247) + IT_0231*IT_0252 + IT_0232
      *IT_0253;
    const ccomplex_t IT_0255 = (conj(IT_0039)*(IT_0059 + IT_0060) + IT_0039*
      (conj(IT_0059) + conj(IT_0060)))*IT_0064 + (conj(IT_0098)*IT_0100 +
       IT_0098*conj(IT_0100) + conj(IT_0115)*IT_0117 + IT_0115*conj(IT_0117))
      *IT_0122 + (IT_0098*conj(IT_0098) + IT_0100*conj(IT_0100) + IT_0115*conj
      (IT_0115) + IT_0117*conj(IT_0117))*IT_0125 + -IT_0138*IT_0143 + (-2)
      *IT_0144*(IT_0137*((-12)*IT_0118 + 24*IT_0120 + (-12)*IT_0140) + 0.5
      *IT_0143*IT_0144) + IT_0061*(IT_0059*conj(IT_0059) + IT_0060*conj(IT_0060)
       + IT_0153*conj(IT_0153) + IT_0160*conj(IT_0160)) + IT_0062*(conj(IT_0153)
      *IT_0163 + IT_0153*IT_0166 + conj(IT_0059)*IT_0169 + IT_0059*IT_0172) +
       IT_0140*(((-18)*conj(IT_0098) + (-6)*conj(IT_0100))*IT_0115 + conj
      (IT_0115)*IT_0174 + conj(IT_0117)*IT_0177 + IT_0117*IT_0180) + (conj
      (IT_0059)*IT_0153 + IT_0059*conj(IT_0153) + conj(IT_0060)*IT_0160 +
       IT_0060*conj(IT_0160))*IT_0181 + (s_34 + IT_0182)*((conj(IT_0153) + conj
      (IT_0160))*IT_0183 + (IT_0153 + IT_0160)*IT_0184) + (-6)*conj(IT_0217)*(
      (IT_0098 + 3*IT_0100 + -IT_0115 + (-3)*IT_0117)*IT_0219 + (
      -0.166666666666667)*IT_0059*IT_0221 + (-0.166666666666667)*IT_0060*IT_0222
       + (-0.166666666666667)*IT_0223*(IT_0168 + IT_0224)) + IT_0225*(conj
      (IT_0115)*IT_0163 + IT_0115*IT_0166 + conj(IT_0117)*IT_0169 + IT_0117
      *IT_0172 + conj(IT_0153)*IT_0174 + conj(IT_0059)*IT_0177 + IT_0059*IT_0180
       + conj(IT_0098)*IT_0226 + conj(IT_0100)*IT_0227) + IT_0229*(IT_0137
      *IT_0230 + IT_0144*IT_0231 + IT_0232*IT_0233 + IT_0234*IT_0235) + (IT_0144
      *IT_0230 + IT_0137*IT_0231 + IT_0233*IT_0234 + IT_0232*IT_0235)*IT_0237 +
       IT_0039*conj(IT_0039)*((-36)*IT_0062 + IT_0238) + IT_0217*(conj(IT_0059)
      *IT_0221 + conj(IT_0060)*IT_0222 + conj(IT_0217)*(36*IT_0062 + IT_0238) +
       6*IT_0219*(conj(IT_0115) + 3*conj(IT_0117) + -IT_0239) + IT_0223*(IT_0171
       + IT_0240)) + IT_0241*IT_0248 + IT_0251*IT_0254;
    return create_ccomplex_return(IT_0255);
}

