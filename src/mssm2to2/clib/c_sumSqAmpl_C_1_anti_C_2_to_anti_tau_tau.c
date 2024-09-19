#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_tau_tau(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
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
    const ccomplex_t IT_0010 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0010;
    const ccomplex_t IT_0013 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + -IT_0012 + (-2)
      *IT_0014);
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = pow(m_C_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0018 + -IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = s_14*s_23;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0025 = IT_0016*IT_0024;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = pow(m_tau, 2);
    const ccomplex_t IT_0029 = s_12*IT_0028;
    const ccomplex_t IT_0030 = m_C_1*m_tau;
    const ccomplex_t IT_0031 = -s_24;
    const ccomplex_t IT_0032 = s_23 + IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = 4*IT_0033;
    const ccomplex_t IT_0035 = pow(m_W, -1);
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = cpow(IT_0036, -1);
    const ccomplex_t IT_0038 = sin(alpha);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em*m_tau*IT_0005*IT_0035
      *IT_0037*IT_0038;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = cos(alpha);
    const ccomplex_t IT_0042 = V_u1*U_Wm2;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = U_d2*V_Wp1;
    const ccomplex_t IT_0045 = IT_0038*IT_0044;
    const ccomplex_t IT_0046 = IT_0005*(IT_0043 + -IT_0045);
    const ccomplex_t IT_0047 = 1.4142135623731*e_em*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0040*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0018 + -IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*e_em*m_tau*IT_0005*IT_0035
      *IT_0037*IT_0041;
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = IT_0041*IT_0044;
    const ccomplex_t IT_0057 = IT_0038*IT_0042;
    const ccomplex_t IT_0058 = IT_0005*(IT_0056 + IT_0057);
    const ccomplex_t IT_0059 = 1.4142135623731*e_em*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = IT_0055*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0053 + -IT_0065;
    const ccomplex_t IT_0067 = m_C_2*m_tau;
    const ccomplex_t IT_0068 = -s_14;
    const ccomplex_t IT_0069 = s_13 + IT_0068;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = (-4)*IT_0070;
    const ccomplex_t IT_0072 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0073 = IT_0041*IT_0072;
    const ccomplex_t IT_0074 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0075 = IT_0038*IT_0074;
    const ccomplex_t IT_0076 = IT_0005*(IT_0073 + IT_0075);
    const ccomplex_t IT_0077 = 1.4142135623731*e_em*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = IT_0055*IT_0079;
    const ccomplex_t IT_0081 = IT_0063*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = IT_0041*IT_0074;
    const ccomplex_t IT_0084 = IT_0038*IT_0072;
    const ccomplex_t IT_0085 = IT_0005*(IT_0083 + -IT_0084);
    const ccomplex_t IT_0086 = 1.4142135623731*e_em*IT_0085;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0040*IT_0088;
    const ccomplex_t IT_0090 = IT_0051*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = -IT_0082 + -IT_0091;
    const ccomplex_t IT_0093 = s_34*m_C_1*m_C_2;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = cpow(s_13 + (-0.5)*IT_0018 + (-0.5)*IT_0028 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0103 = IT_0003*IT_0102;
    const ccomplex_t IT_0104 = IT_0006*IT_0102;
    const ccomplex_t IT_0105 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0106 = IT_0006*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0103 + -IT_0104 + (-2)
      *IT_0106);
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0009*IT_0108;
    const ccomplex_t IT_0110 = IT_0020*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = m_C_1*m_C_2*IT_0028;
    const ccomplex_t IT_0114 = IT_0024*IT_0108;
    const ccomplex_t IT_0115 = IT_0020*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = s_23 + s_24;
    const ccomplex_t IT_0119 = IT_0030*IT_0118;
    const ccomplex_t IT_0120 = sin(beta);
    const ccomplex_t IT_0121 = e_em*m_tau*IT_0005*IT_0035*IT_0037*IT_0120;
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = IT_0036*IT_0042;
    const ccomplex_t IT_0124 = IT_0044*IT_0120;
    const ccomplex_t IT_0125 = IT_0005*(IT_0123 + IT_0124);
    const ccomplex_t IT_0126 = 1.4142135623731*e_em*IT_0125;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0122*IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = pow(m_Z, 2);
    const ccomplex_t IT_0133 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + IT_0132 
      + -reg_prop, -1);
    const ccomplex_t IT_0134 = IT_0036*IT_0044;
    const ccomplex_t IT_0135 = IT_0042*IT_0120;
    const ccomplex_t IT_0136 = IT_0005*(IT_0134 + -IT_0135);
    const ccomplex_t IT_0137 = 1.4142135623731*e_em*IT_0136;
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = e_em*m_tau*IT_0005*IT_0035;
    const ccomplex_t IT_0140 = 0.5*IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = IT_0133*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = -IT_0131 + -IT_0143;
    const ccomplex_t IT_0145 = 16*conj(IT_0022)*IT_0023 + 8*conj(IT_0027)
      *IT_0029 + IT_0034*conj(IT_0066) + IT_0071*conj(IT_0092) + (-4)*IT_0093*
      (conj(IT_0101) + (-2)*conj(IT_0112)) + 16*IT_0113*conj(IT_0117) + (-4)
      *IT_0119*conj(IT_0144);
    const ccomplex_t IT_0146 = s_13 + s_14;
    const ccomplex_t IT_0147 = IT_0067*IT_0146;
    const ccomplex_t IT_0148 = IT_0036*IT_0074;
    const ccomplex_t IT_0149 = IT_0072*IT_0120;
    const ccomplex_t IT_0150 = IT_0005*(IT_0148 + IT_0149);
    const ccomplex_t IT_0151 = 1.4142135623731*e_em*IT_0150;
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0122*IT_0152;
    const ccomplex_t IT_0154 = IT_0129*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = IT_0036*IT_0072;
    const ccomplex_t IT_0157 = IT_0074*IT_0120;
    const ccomplex_t IT_0158 = IT_0005*(IT_0156 + -IT_0157);
    const ccomplex_t IT_0159 = 1.4142135623731*e_em*IT_0158;
    const ccomplex_t IT_0160 = 0.5*IT_0159;
    const ccomplex_t IT_0161 = IT_0140*IT_0160;
    const ccomplex_t IT_0162 = IT_0133*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = -IT_0155 + -IT_0163;
    const ccomplex_t IT_0165 = s_23*m_C_1*m_tau;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0005*IT_0035*IT_0037;
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = IT_0095*IT_0167;
    const ccomplex_t IT_0169 = IT_0098*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0005*IT_0035*IT_0037;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = IT_0167*IT_0173;
    const ccomplex_t IT_0175 = IT_0098*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = s_14*m_C_2*m_tau;
    const ccomplex_t IT_0179 = IT_0096*IT_0173;
    const ccomplex_t IT_0180 = IT_0098*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = 0.5*IT_0181;
    const ccomplex_t IT_0183 = 4*IT_0022;
    const ccomplex_t IT_0184 = s_24*m_C_1*m_tau;
    const ccomplex_t IT_0185 = s_13*s_24;
    const ccomplex_t IT_0186 = s_13*m_C_2*m_tau;
    const ccomplex_t IT_0187 = 8*IT_0027;
    const ccomplex_t IT_0188 = s_12*s_34;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = -IT_0185;
    const ccomplex_t IT_0191 = IT_0023 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = -IT_0029;
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = 2*IT_0119;
    const ccomplex_t IT_0195 = m_C_1*m_C_2;
    const ccomplex_t IT_0196 = s_34 + IT_0028;
    const ccomplex_t IT_0197 = IT_0195*IT_0196;
    const ccomplex_t IT_0198 = 2*IT_0197;
    const ccomplex_t IT_0199 = 4*IT_0147;
    const ccomplex_t IT_0200 = s_12*IT_0196;
    const ccomplex_t IT_0201 = 8*IT_0200;
    const ccomplex_t IT_0202 = (-8)*IT_0197;
    const ccomplex_t IT_0203 = (-4)*conj(IT_0022);
    const ccomplex_t IT_0204 = (-4)*IT_0147;
    const ccomplex_t IT_0205 = (-2)*IT_0147;
    const ccomplex_t IT_0206 = conj(IT_0171)*IT_0193 + (-0.5)*conj(IT_0177)
      *IT_0194 + (-0.5)*conj(IT_0182)*IT_0198 + (-0.5)*conj(IT_0112)*IT_0199 + (
      -0.5)*conj(IT_0144)*IT_0201 + (-0.5)*conj(IT_0164)*IT_0202 + (-2)*IT_0119*
      (conj(IT_0027) + 0.25*IT_0203) + (-0.5)*conj(IT_0117)*IT_0204 + (-0.5)
      *conj(IT_0101)*IT_0205;
    const ccomplex_t IT_0207 = -IT_0028;
    const ccomplex_t IT_0208 = s_34 + IT_0207;
    const ccomplex_t IT_0209 = s_12*IT_0208;
    const ccomplex_t IT_0210 = IT_0195*IT_0208;
    const ccomplex_t IT_0211 = (-8)*IT_0210;
    const ccomplex_t IT_0212 = conj(IT_0112) + conj(IT_0117);
    const ccomplex_t IT_0213 = IT_0029 + IT_0191;
    const ccomplex_t IT_0214 = 2*IT_0213;
    const ccomplex_t IT_0215 = (-2)*IT_0210;
    const ccomplex_t IT_0216 = -IT_0215;
    const ccomplex_t IT_0217 = (-2)*IT_0033;
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = (-2)*IT_0070;
    const ccomplex_t IT_0220 = -IT_0219;
    const ccomplex_t IT_0221 = 4*IT_0185;
    const ccomplex_t IT_0222 = (-8)*IT_0185;
    const ccomplex_t IT_0223 = (-4)*IT_0029;
    const ccomplex_t IT_0224 = (-4)*IT_0184;
    const ccomplex_t IT_0225 = (-4)*IT_0186;
    const ccomplex_t IT_0226 = (-4)*IT_0113;
    const ccomplex_t IT_0227 = 2*IT_0033;
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = conj(IT_0027)*IT_0113 + (-0.125)*conj(IT_0164)
      *IT_0194 + (-0.125)*IT_0093*IT_0203 + (-0.125)*conj(IT_0144)*IT_0205 + (
      -0.125)*conj(IT_0066)*IT_0220 + (-0.125)*conj(IT_0101)*IT_0221 + (-0.125)
      *conj(IT_0112)*IT_0222 + (-0.125)*conj(IT_0117)*IT_0223 + (-0.125)*conj
      (IT_0182)*IT_0224 + (-0.125)*conj(IT_0171)*IT_0225 + (-0.125)*conj(IT_0177
      )*IT_0226 + (-0.125)*conj(IT_0092)*IT_0228;
    const ccomplex_t IT_0230 = 2*IT_0070;
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = 4*IT_0178;
    const ccomplex_t IT_0233 = 4*IT_0165;
    const ccomplex_t IT_0234 = 4*IT_0093;
    const ccomplex_t IT_0235 = 4*IT_0119;
    const ccomplex_t IT_0236 = 8*conj(IT_0027);
    const ccomplex_t IT_0237 = IT_0034*conj(IT_0092);
    const ccomplex_t IT_0238 = conj(IT_0066)*IT_0071;
    const ccomplex_t IT_0239 = 8*IT_0029*conj(IT_0112) + 16*conj(IT_0022)
      *IT_0113 + 16*IT_0023*conj(IT_0117) + conj(IT_0144)*IT_0204 + conj(IT_0101
      )*IT_0223 + conj(IT_0171)*IT_0232 + conj(IT_0182)*IT_0233 + conj(IT_0177)
      *IT_0234 + conj(IT_0164)*IT_0235 + IT_0093*IT_0236 + IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = (-4)*IT_0119;
    const ccomplex_t IT_0241 = conj(IT_0022)*IT_0093 + IT_0029*conj(IT_0117) +
       0.125*IT_0113*(16*conj(IT_0027) + 8*conj(IT_0177)) + conj(IT_0182)
      *IT_0184 + 2*conj(IT_0112)*IT_0185 + conj(IT_0171)*IT_0186 + 0.125*conj
      (IT_0144)*IT_0199 + 0.125*conj(IT_0101)*IT_0222 + 0.125*IT_0237 + 0.125
      *IT_0238 + 0.125*conj(IT_0164)*IT_0240;
    const ccomplex_t IT_0242 = 4*conj(IT_0022);
    const ccomplex_t IT_0243 = 0.125*conj(IT_0066);
    const ccomplex_t IT_0244 = 0.125*conj(IT_0092);
    const ccomplex_t IT_0245 = 4*IT_0113;
    const ccomplex_t IT_0246 = 0.125*conj(IT_0171);
    const ccomplex_t IT_0247 = 4*IT_0186;
    const ccomplex_t IT_0248 = 0.125*conj(IT_0177);
    const ccomplex_t IT_0249 = 0.125*conj(IT_0101);
    const ccomplex_t IT_0250 = 0.125*conj(IT_0182);
    const ccomplex_t IT_0251 = 2*IT_0193;
    const ccomplex_t IT_0252 = 0.125*conj(IT_0164);
    const ccomplex_t IT_0253 = 0.125*conj(IT_0117);
    const ccomplex_t IT_0254 = 0.125*conj(IT_0144);
    const ccomplex_t IT_0255 = conj(IT_0112)*IT_0184 + 0.125*IT_0186*IT_0236 +
       0.125*IT_0178*IT_0242 + IT_0216*IT_0243 + IT_0214*IT_0244 + IT_0245
      *IT_0246 + IT_0247*IT_0248 + IT_0224*IT_0249 + IT_0221*IT_0250 + IT_0251
      *IT_0252 + IT_0233*IT_0253 + IT_0198*IT_0254;
    const ccomplex_t IT_0256 = 4*IT_0184;
    const ccomplex_t IT_0257 = (-2)*IT_0197;
    const ccomplex_t IT_0258 = (-0.25)*conj(IT_0144);
    const ccomplex_t IT_0259 = conj(IT_0112)*IT_0186 + 0.125*IT_0184*IT_0236 +
       0.125*IT_0165*IT_0242 + IT_0214*IT_0243 + IT_0216*IT_0244 + IT_0221
      *IT_0246 + IT_0225*IT_0249 + IT_0245*IT_0250 + IT_0232*IT_0253 + IT_0248
      *IT_0256 + IT_0252*IT_0257 + IT_0193*IT_0258;
    const ccomplex_t IT_0260 = conj(IT_0101)*IT_0194 + conj(IT_0164)*IT_0201 +
       conj(IT_0177)*IT_0205 + conj(IT_0117)*IT_0235 + conj(IT_0112)*IT_0240 + (
      -4)*IT_0147*(conj(IT_0027) + (-0.25)*IT_0242) + conj(IT_0182)*IT_0251 +
       conj(IT_0171)*IT_0257 + (-4)*IT_0202*IT_0258;
    const ccomplex_t IT_0261 = 16*IT_0027;
    const ccomplex_t IT_0262 = conj(IT_0112)*IT_0113;
    const ccomplex_t IT_0263 = IT_0185*IT_0236 + IT_0029*IT_0242 + 8*IT_0218
      *IT_0243 + 8*IT_0231*IT_0244 + 8*IT_0221*IT_0248 + 8*IT_0226*IT_0249 + 8
      *IT_0247*IT_0250 + 8*IT_0205*IT_0252 + 8*IT_0234*IT_0253 + 8*IT_0194
      *IT_0254 + 8*IT_0246*IT_0256 + 8*IT_0262;
    const ccomplex_t IT_0264 = IT_0022*IT_0145 + 8*IT_0027*(conj(IT_0022)
      *IT_0029 + 0.125*IT_0034*conj(IT_0066) + 0.125*IT_0071*conj(IT_0092) + 
      -conj(IT_0101)*IT_0113 + 0.5*IT_0119*conj(IT_0144) + (-0.5)*IT_0147*conj
      (IT_0164)) + (IT_0147*conj(IT_0164) + IT_0165*conj(IT_0171) + IT_0029*conj
      (IT_0177) + IT_0178*conj(IT_0182))*IT_0183 + (IT_0093*conj(IT_0117) + conj
      (IT_0171)*IT_0184 + conj(IT_0177)*IT_0185 + conj(IT_0182)*IT_0186)*IT_0187
       + (-2)*IT_0144*IT_0206 + IT_0066*((conj(IT_0022) + conj(IT_0027))*IT_0034
       + 8*conj(IT_0066)*IT_0209 + conj(IT_0092)*IT_0211 + IT_0071*IT_0212 +
       conj(IT_0171)*IT_0214 + conj(IT_0182)*IT_0216 + conj(IT_0177)*IT_0218 +
       conj(IT_0101)*IT_0220) + (-8)*IT_0101*IT_0229 + IT_0092*((conj(IT_0022) +
       conj(IT_0027))*IT_0071 + 8*conj(IT_0092)*IT_0209 + conj(IT_0066)*IT_0211 
      + IT_0034*IT_0212 + conj(IT_0182)*IT_0214 + conj(IT_0171)*IT_0216 + conj
      (IT_0101)*IT_0228 + conj(IT_0177)*IT_0231) + IT_0117*IT_0239 + 8*IT_0112
      *IT_0241 + 8*IT_0182*IT_0255 + 8*IT_0171*IT_0259 + IT_0164*IT_0260 +
       IT_0261*(conj(IT_0027)*IT_0185 + IT_0262) + IT_0177*IT_0263;
    return create_ccomplex_return(IT_0264);
}

