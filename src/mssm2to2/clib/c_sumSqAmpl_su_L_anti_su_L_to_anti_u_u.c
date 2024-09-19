#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_L_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_L_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = pow(m_u, 2);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0002*IT_0004*IT_0006 + (-0.5)*m_W*(IT_0002*IT_0003 + IT_0007*IT_0008)*
      (IT_0006 + -1./3*IT_0005*IT_0010));
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*m_u*e_em*IT_0000*IT_0002
      *IT_0004*IT_0006;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = pow(m_su_L, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0004*IT_0006*IT_0008 + 0.5*m_W*(IT_0002*IT_0007 + -IT_0003*IT_0008)*
      (IT_0006 + -1./3*IT_0005*IT_0010));
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_u*e_em*IT_0000*IT_0004
      *IT_0006*IT_0008;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0019 + -IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0009, -1);
    const ccomplex_t IT_0030 = IT_0005*IT_0029;
    const ccomplex_t IT_0031 = e_em*IT_0030;
    const ccomplex_t IT_0032 = IT_0006*IT_0009;
    const ccomplex_t IT_0033 = e_em*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0031 + (-3)*IT_0033);
    const ccomplex_t IT_0035 = (-0.166666666666667)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*-3)*e_em*(IT_0006*IT_0009 + -1.
      /3*IT_0005*IT_0029);
    const ccomplex_t IT_0040 = 0.166666666666667*IT_0039;
    const ccomplex_t IT_0041 = m_u*IT_0040;
    const ccomplex_t IT_0042 = 0.333333333333333*IT_0039;
    const ccomplex_t IT_0043 = m_u*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0041 + IT_0044;
    const ccomplex_t IT_0046 = IT_0038*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0029;
    const ccomplex_t IT_0048 = (-0.666666666666667)*IT_0047;
    const ccomplex_t IT_0049 = IT_0036*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0045*IT_0050;
    const ccomplex_t IT_0052 = pow(m_sG, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -IT_0052 +
       reg_prop, -1);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0055 = cpow(IT_0054, 2);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0053*IT_0056;
    const ccomplex_t IT_0058 = m_u*IT_0057;
    const ccomplex_t IT_0059 = -IT_0046 + IT_0051 + 0.5*IT_0058;
    const ccomplex_t IT_0060 = (IT_0041 + -IT_0043)*(IT_0038 + -IT_0050);
    const ccomplex_t IT_0061 = (-18)*IT_0001;
    const ccomplex_t IT_0062 = 18*s_34;
    const ccomplex_t IT_0063 = IT_0061 + IT_0062;
    const ccomplex_t IT_0064 = conj(N_B1)*e_em;
    const ccomplex_t IT_0065 = IT_0029*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W1)*e_em;
    const ccomplex_t IT_0068 = IT_0006*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + 3*IT_0069);
    const ccomplex_t IT_0071 = 0.166666666666667*IT_0070;
    const ccomplex_t IT_0072 = N_B1*e_em;
    const ccomplex_t IT_0073 = IT_0029*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = N_W1*e_em;
    const ccomplex_t IT_0076 = IT_0006*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + 3*IT_0077);
    const ccomplex_t IT_0079 = 0.166666666666667*IT_0078;
    const ccomplex_t IT_0080 = IT_0071*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = conj(N_B3)*e_em;
    const ccomplex_t IT_0085 = IT_0029*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = conj(N_W3)*e_em;
    const ccomplex_t IT_0088 = IT_0006*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + 3*IT_0089);
    const ccomplex_t IT_0091 = 0.166666666666667*IT_0090;
    const ccomplex_t IT_0092 = N_B3*e_em;
    const ccomplex_t IT_0093 = IT_0029*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = N_W3*e_em;
    const ccomplex_t IT_0096 = IT_0006*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + 3*IT_0097);
    const ccomplex_t IT_0099 = 0.166666666666667*IT_0098;
    const ccomplex_t IT_0100 = IT_0091*IT_0099;
    const ccomplex_t IT_0101 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = conj(N_B4)*e_em;
    const ccomplex_t IT_0105 = IT_0029*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = conj(N_W4)*e_em;
    const ccomplex_t IT_0108 = IT_0006*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + 3*IT_0109);
    const ccomplex_t IT_0111 = 0.166666666666667*IT_0110;
    const ccomplex_t IT_0112 = N_B4*e_em;
    const ccomplex_t IT_0113 = IT_0029*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = N_W4*e_em;
    const ccomplex_t IT_0116 = IT_0006*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0114 + 3*IT_0117);
    const ccomplex_t IT_0119 = 0.166666666666667*IT_0118;
    const ccomplex_t IT_0120 = IT_0111*IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0125 = N_B2*e_em;
    const ccomplex_t IT_0126 = IT_0029*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = N_W2*e_em;
    const ccomplex_t IT_0129 = IT_0006*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 3*IT_0130);
    const ccomplex_t IT_0132 = 0.166666666666667*IT_0131;
    const ccomplex_t IT_0133 = conj(N_B2)*e_em;
    const ccomplex_t IT_0134 = IT_0029*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = conj(N_W2)*e_em;
    const ccomplex_t IT_0137 = IT_0006*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0135 + 3*IT_0138);
    const ccomplex_t IT_0140 = 0.166666666666667*IT_0139;
    const ccomplex_t IT_0141 = IT_0132*IT_0140;
    const ccomplex_t IT_0142 = IT_0124*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = 0.166666666666667*IT_0057 + -IT_0083 + -IT_0103
       + -IT_0123 + -IT_0143;
    const ccomplex_t IT_0145 = IT_0038*IT_0042;
    const ccomplex_t IT_0146 = (-0.5)*IT_0057 + -IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0000*IT_0004*IT_0006;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0000*IT_0004*IT_0006;
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = IT_0148*IT_0150;
    const ccomplex_t IT_0152 = IT_0081*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0000*IT_0004*IT_0006;
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0000*IT_0004*IT_0006;
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = IT_0155*IT_0157;
    const ccomplex_t IT_0159 = IT_0101*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0000*IT_0004*IT_0006;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0000*IT_0004*IT_0006;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = IT_0162*IT_0164;
    const ccomplex_t IT_0166 = IT_0121*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0000*IT_0004*IT_0006;
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0000*IT_0004*IT_0006;
    const ccomplex_t IT_0171 = (-0.5)*IT_0170;
    const ccomplex_t IT_0172 = IT_0169*IT_0171;
    const ccomplex_t IT_0173 = IT_0124*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = -IT_0153 + -IT_0160 + -IT_0167 + -IT_0174;
    const ccomplex_t IT_0176 = IT_0042*IT_0050;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = s_34*IT_0016;
    const ccomplex_t IT_0179 = (-6)*IT_0178;
    const ccomplex_t IT_0180 = s_13*s_14;
    const ccomplex_t IT_0181 = 12*IT_0180;
    const ccomplex_t IT_0182 = IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = (-18)*IT_0178;
    const ccomplex_t IT_0184 = 36*IT_0180;
    const ccomplex_t IT_0185 = IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0187 = 0.666666666666667*IT_0186;
    const ccomplex_t IT_0188 = (-1.33333333333333)*IT_0186;
    const ccomplex_t IT_0189 = IT_0187*IT_0188;
    const ccomplex_t IT_0190 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0191 = IT_0189*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0194 = (-2)*IT_0193;
    const ccomplex_t IT_0195 = IT_0190*IT_0193;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = IT_0194*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0192 + (-0.0833333333333333)*IT_0197;
    const ccomplex_t IT_0199 = cpow(IT_0198, 2);
    const ccomplex_t IT_0200 = (-36)*IT_0178;
    const ccomplex_t IT_0201 = IT_0001*IT_0016;
    const ccomplex_t IT_0202 = (-36)*IT_0201;
    const ccomplex_t IT_0203 = 72*IT_0180;
    const ccomplex_t IT_0204 = IT_0200 + IT_0202 + IT_0203;
    const ccomplex_t IT_0205 = 0.25*IT_0197;
    const ccomplex_t IT_0206 = (-18)*IT_0144 + (-6)*IT_0146;
    const ccomplex_t IT_0207 = (-18)*IT_0146;
    const ccomplex_t IT_0208 = (-6)*IT_0144;
    const ccomplex_t IT_0209 = IT_0207 + IT_0208;
    const ccomplex_t IT_0210 = (-18)*conj(IT_0146);
    const ccomplex_t IT_0211 = (-6)*conj(IT_0144);
    const ccomplex_t IT_0212 = IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = (-18)*IT_0059;
    const ccomplex_t IT_0214 = IT_0079*IT_0150;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0216 = IT_0081*IT_0214*IT_0215;
    const ccomplex_t IT_0217 = IT_0132*IT_0169;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0219 = IT_0124*IT_0217*IT_0218;
    const ccomplex_t IT_0220 = IT_0099*IT_0157;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0222 = IT_0101*IT_0220*IT_0221;
    const ccomplex_t IT_0223 = IT_0119*IT_0164;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0225 = IT_0121*IT_0223*IT_0224;
    const ccomplex_t IT_0226 = m_u*IT_0083;
    const ccomplex_t IT_0227 = m_u*IT_0103;
    const ccomplex_t IT_0228 = m_u*IT_0123;
    const ccomplex_t IT_0229 = m_u*IT_0143;
    const ccomplex_t IT_0230 = (-0.166666666666667)*IT_0058 + -IT_0216 + 
      -IT_0219 + -IT_0222 + -IT_0225 + IT_0226 + IT_0227 + IT_0228 + IT_0229;
    const ccomplex_t IT_0231 = (-6)*IT_0230;
    const ccomplex_t IT_0232 = IT_0213 + IT_0231;
    const ccomplex_t IT_0233 = (-18)*conj(IT_0059);
    const ccomplex_t IT_0234 = (-6)*conj(IT_0230);
    const ccomplex_t IT_0235 = IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = IT_0071*IT_0148;
    const ccomplex_t IT_0237 = IT_0081*IT_0215*IT_0236;
    const ccomplex_t IT_0238 = IT_0091*IT_0155;
    const ccomplex_t IT_0239 = IT_0101*IT_0221*IT_0238;
    const ccomplex_t IT_0240 = IT_0111*IT_0162;
    const ccomplex_t IT_0241 = IT_0121*IT_0224*IT_0240;
    const ccomplex_t IT_0242 = m_u*IT_0153;
    const ccomplex_t IT_0243 = m_u*IT_0160;
    const ccomplex_t IT_0244 = m_u*IT_0167;
    const ccomplex_t IT_0245 = IT_0140*IT_0171;
    const ccomplex_t IT_0246 = IT_0124*IT_0218*IT_0245;
    const ccomplex_t IT_0247 = m_u*IT_0174;
    const ccomplex_t IT_0248 = -IT_0237 + -IT_0239 + -IT_0241 + IT_0242 +
       IT_0243 + IT_0244 + -IT_0246 + IT_0247;
    const ccomplex_t IT_0249 = (-18)*IT_0230;
    const ccomplex_t IT_0250 = (-6)*IT_0059;
    const ccomplex_t IT_0251 = IT_0249 + IT_0250;
    const ccomplex_t IT_0252 = (-18)*conj(IT_0230);
    const ccomplex_t IT_0253 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0254 = IT_0252 + IT_0253;
    const ccomplex_t IT_0255 = 6*s_34;
    const ccomplex_t IT_0256 = -IT_0001;
    const ccomplex_t IT_0257 = s_34 + IT_0256;
    const ccomplex_t IT_0258 = 6*IT_0028;
    const ccomplex_t IT_0259 = 6*conj(IT_0028);
    const ccomplex_t IT_0260 = m_u*s_13;
    const ccomplex_t IT_0261 = 6*conj(IT_0248);
    const ccomplex_t IT_0262 = 6*conj(IT_0230);
    const ccomplex_t IT_0263 = 18*conj(IT_0059);
    const ccomplex_t IT_0264 = 6*conj(IT_0059);
    const ccomplex_t IT_0265 = 18*conj(IT_0230);
    const ccomplex_t IT_0266 = 18*conj(IT_0248);
    const ccomplex_t IT_0267 = IT_0059*(conj(IT_0144) + 3*conj(IT_0146)) + 3*
      (conj(IT_0144) + 0.333333333333333*conj(IT_0146))*IT_0230 + 3*conj(IT_0175
      )*(0.333333333333333*IT_0060 + IT_0248) + 3*conj(IT_0177)*(IT_0060 + 1./3
      *IT_0248) + 3*IT_0177*(conj(IT_0060) + 0.0555555555555556*IT_0261) +
       0.166666666666667*IT_0146*(IT_0262 + IT_0263) + 0.166666666666667*IT_0144
      *(IT_0264 + IT_0265) + IT_0175*(conj(IT_0060) + 0.166666666666667*IT_0266);
    const ccomplex_t IT_0268 = -s_14;
    const ccomplex_t IT_0269 = s_13 + IT_0268;
    const ccomplex_t IT_0270 = m_u*IT_0269;
    const ccomplex_t IT_0271 = conj(IT_0144) + conj(IT_0175);
    const ccomplex_t IT_0272 = IT_0144 + IT_0175;
    const ccomplex_t IT_0273 = conj(IT_0146) + conj(IT_0177);
    const ccomplex_t IT_0274 = 18*IT_0028;
    const ccomplex_t IT_0275 = IT_0146 + IT_0177;
    const ccomplex_t IT_0276 = 18*conj(IT_0028);
    const ccomplex_t IT_0277 = (-18)*IT_0248;
    const ccomplex_t IT_0278 = (-6)*IT_0248;
    const ccomplex_t IT_0279 = IT_0258*IT_0271 + IT_0259*IT_0272 + IT_0273
      *IT_0274 + IT_0275*IT_0276 + IT_0205*((-12)*IT_0028 + 12*conj(IT_0028) + (
      -6)*IT_0060 + 6*conj(IT_0060) + IT_0249 + IT_0250 + IT_0264 + IT_0265 +
       IT_0266 + IT_0277) + IT_0198*((-36)*IT_0028 + 36*conj(IT_0028) + (-18)
      *IT_0060 + 18*conj(IT_0060) + IT_0213 + IT_0231 + IT_0261 + IT_0262 +
       IT_0263 + IT_0278);
    const ccomplex_t IT_0280 = m_u*s_14;
    const ccomplex_t IT_0281 = (-6)*IT_0201;
    const ccomplex_t IT_0282 = IT_0179 + IT_0181 + IT_0281;
    const ccomplex_t IT_0283 = -IT_0205;
    const ccomplex_t IT_0284 = -IT_0198;
    const ccomplex_t IT_0285 = (-18)*IT_0201;
    const ccomplex_t IT_0286 = IT_0183 + IT_0184 + IT_0285;
    const ccomplex_t IT_0287 = (36*s_34 + (-36)*IT_0001)*IT_0028*conj(IT_0028)
       + (conj(IT_0028)*(IT_0059 + IT_0060) + IT_0028*(conj(IT_0059) + conj
      (IT_0060)))*IT_0063 + (conj(IT_0144)*IT_0146 + IT_0144*conj(IT_0146) +
       conj(IT_0175)*IT_0177 + IT_0175*conj(IT_0177))*IT_0182 + (IT_0144*conj
      (IT_0144) + IT_0146*conj(IT_0146) + IT_0175*conj(IT_0175) + IT_0177*conj
      (IT_0177))*IT_0185 + -IT_0199*IT_0204 + (-2)*IT_0205*(IT_0198*((-12)
      *IT_0178 + 24*IT_0180 + (-12)*IT_0201) + 0.5*IT_0204*IT_0205) + IT_0201*((
      (-18)*conj(IT_0144) + (-6)*conj(IT_0146))*IT_0175 + conj(IT_0175)*IT_0206 
      + conj(IT_0177)*IT_0209 + IT_0177*IT_0212) + IT_0062*(IT_0059*conj(IT_0059
      ) + IT_0060*conj(IT_0060) + IT_0230*conj(IT_0230) + IT_0248*conj(IT_0248))
       + IT_0001*(conj(IT_0060)*IT_0232 + IT_0060*IT_0235 + conj(IT_0248)
      *IT_0251 + IT_0248*IT_0254) + (conj(IT_0059)*IT_0230 + IT_0059*conj
      (IT_0230) + conj(IT_0060)*IT_0248 + IT_0060*conj(IT_0248))*IT_0255 +
       IT_0257*((conj(IT_0230) + conj(IT_0248))*IT_0258 + (IT_0230 + IT_0248)
      *IT_0259) + 6*IT_0260*IT_0267 + IT_0270*IT_0279 + (conj(IT_0060)*IT_0209 +
       IT_0060*IT_0212 + conj(IT_0177)*IT_0232 + IT_0177*IT_0235 + IT_0206*conj
      (IT_0248) + conj(IT_0175)*IT_0251 + IT_0175*IT_0254 + conj(IT_0144)
      *IT_0277 + conj(IT_0146)*IT_0278)*IT_0280 + IT_0282*(IT_0198*IT_0271 +
       IT_0205*IT_0273 + IT_0275*IT_0283 + IT_0272*IT_0284) + (IT_0205*IT_0271 +
       IT_0198*IT_0273 + IT_0272*IT_0283 + IT_0275*IT_0284)*IT_0286;
    return create_ccomplex_return(IT_0287);
}

