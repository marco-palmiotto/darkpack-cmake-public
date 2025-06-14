#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
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
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = e_em*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + (-3)*IT_0007);
    const ccomplex_t IT_0009 = (-0.166666666666667)*IT_0008;
    const ccomplex_t IT_0010 = pow(m_sb_2, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0005)*U_sb_01*conj(U_sb_01) + IT_0001*IT_0002
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0015 = (-0.333333333333333)*IT_0014;
    const ccomplex_t IT_0016 = m_u*IT_0015;
    const ccomplex_t IT_0017 = (-0.666666666666667)*IT_0014;
    const ccomplex_t IT_0018 = m_u*IT_0017;
    const ccomplex_t IT_0019 = IT_0013*(IT_0016 + -IT_0018);
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0022 = (-0.666666666666667)*IT_0021;
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (IT_0016 + -IT_0018)*IT_0024;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = cos(alpha);
    const ccomplex_t IT_0030 = sin(beta);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_u*e_em*IT_0005*IT_0028
      *IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = cpow(V_tb, 2);
    const ccomplex_t IT_0035 = sin(alpha);
    const ccomplex_t IT_0036 = m_b*mu_h*IT_0028;
    const ccomplex_t IT_0037 = cos(beta);
    const ccomplex_t IT_0038 = cpow(IT_0037, -1);
    const ccomplex_t IT_0039 = IT_0035*IT_0037;
    const ccomplex_t IT_0040 = IT_0029*IT_0030;
    const ccomplex_t IT_0041 = IT_0039 + IT_0040;
    const ccomplex_t IT_0042 = cpow(IT_0000, -2);
    const ccomplex_t IT_0043 = pow(m_b, 2);
    const ccomplex_t IT_0044 = 2*U_sb_01*(((-0.5)*m_W*IT_0041*(IT_0005 + 1./3
      *IT_0002*IT_0042) + IT_0005*IT_0028*IT_0035*IT_0038*IT_0043)*conj(U_sb_01)
       + 1./2*IT_0005*IT_0029*IT_0036*IT_0038*conj(U_sb_11)) + U_sb_11*(IT_0005
      *IT_0029*IT_0036*IT_0038*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0002
      *IT_0041*IT_0042 + (-3)*IT_0005*IT_0028*IT_0035*IT_0038*IT_0043)*conj
      (U_sb_11));
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(e_em*IT_0044 +
       1.4142135623731*A_b*IT_0034*IT_0035*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0033*IT_0046;
    const ccomplex_t IT_0048 = IT_0027*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + (-2)*IT_0010 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*m_u*e_em*IT_0005*IT_0028
      *IT_0031*IT_0035;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = -IT_0030*IT_0035;
    const ccomplex_t IT_0054 = IT_0029*IT_0037;
    const ccomplex_t IT_0055 = IT_0053 + IT_0054;
    const ccomplex_t IT_0056 = (-2)*U_sb_01*((IT_0005*IT_0028*IT_0029*IT_0038
      *IT_0043 + (-0.5)*m_W*(IT_0005 + 1./3*IT_0002*IT_0042)*IT_0055)*conj
      (U_sb_01) + -1./2*IT_0005*IT_0035*IT_0036*IT_0038*conj(U_sb_11)) + U_sb_11
      *(IT_0005*IT_0035*IT_0036*IT_0038*conj(U_sb_01) + (-2)*(IT_0005*IT_0028
      *IT_0029*IT_0038*IT_0043 + -1./3*m_W*IT_0002*IT_0042*IT_0055)*conj(U_sb_11
      ));
    const ccomplex_t IT_0057 = (0 + _Complex_I*-1)*(e_em*IT_0056 + (
      -1.4142135623731)*A_b*IT_0029*IT_0034*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0052*IT_0058;
    const ccomplex_t IT_0060 = IT_0050*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = -IT_0049 + -IT_0061;
    const ccomplex_t IT_0063 = 18*s_34;
    const ccomplex_t IT_0064 = pow(m_u, 2);
    const ccomplex_t IT_0065 = (-18)*IT_0064;
    const ccomplex_t IT_0066 = IT_0063 + IT_0065;
    const ccomplex_t IT_0067 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0068 = IT_0005*IT_0067;
    const ccomplex_t IT_0069 = e_em*U_Wm1*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = IT_0005*IT_0038*IT_0067;
    const ccomplex_t IT_0072 = m_b*U_d1*e_em*IT_0028*conj(U_sb_11)*V_ub_mod;
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0070 + (-0.5)*IT_0074);
    const ccomplex_t IT_0076 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0077 = IT_0005*IT_0076;
    const ccomplex_t IT_0078 = e_em*conj(U_Wm1)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = IT_0005*IT_0038*IT_0076;
    const ccomplex_t IT_0081 = m_b*conj(U_d1)*e_em*IT_0028*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0079 + (-0.5)*IT_0083);
    const ccomplex_t IT_0085 = IT_0075*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0010 + IT_0064 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = e_em*U_Wm2*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0090 = IT_0068*IT_0089;
    const ccomplex_t IT_0091 = m_b*U_d2*e_em*IT_0028*conj(U_sb_11)*V_ub_mod;
    const ccomplex_t IT_0092 = IT_0071*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + (-0.5)*IT_0093);
    const ccomplex_t IT_0095 = e_em*conj(U_Wm2)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0096 = IT_0077*IT_0095;
    const ccomplex_t IT_0097 = m_b*conj(U_d2)*e_em*IT_0028*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0098 = IT_0080*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + (-0.5)*IT_0099);
    const ccomplex_t IT_0101 = IT_0094*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_23 + IT_0010 + IT_0064 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = -IT_0088 + -IT_0104;
    const ccomplex_t IT_0106 = IT_0013*IT_0017;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *e_em*IT_0005*IT_0028*IT_0031*IT_0067*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0005*IT_0028*IT_0031*IT_0076*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0111 = 0.5*IT_0110;
    const ccomplex_t IT_0112 = IT_0109*IT_0111;
    const ccomplex_t IT_0113 = IT_0086*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0005*IT_0028*IT_0031*IT_0067*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0005*IT_0028*IT_0031*IT_0076*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = IT_0102*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = -IT_0114 + -IT_0121;
    const ccomplex_t IT_0123 = IT_0017*IT_0024;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = s_34*IT_0010;
    const ccomplex_t IT_0126 = (-6)*IT_0125;
    const ccomplex_t IT_0127 = s_13*s_14;
    const ccomplex_t IT_0128 = 12*IT_0127;
    const ccomplex_t IT_0129 = IT_0126 + IT_0128;
    const ccomplex_t IT_0130 = (-18)*IT_0125;
    const ccomplex_t IT_0131 = 36*IT_0127;
    const ccomplex_t IT_0132 = IT_0130 + IT_0131;
    const ccomplex_t IT_0133 = IT_0010*IT_0064;
    const ccomplex_t IT_0134 = IT_0105 + 3*IT_0107;
    const ccomplex_t IT_0135 = 3*conj(IT_0107);
    const ccomplex_t IT_0136 = conj(IT_0105) + IT_0135;
    const ccomplex_t IT_0137 = (-18)*IT_0105;
    const ccomplex_t IT_0138 = (-6)*IT_0107;
    const ccomplex_t IT_0139 = IT_0137 + IT_0138;
    const ccomplex_t IT_0140 = (-18)*conj(IT_0105);
    const ccomplex_t IT_0141 = (-6)*conj(IT_0107);
    const ccomplex_t IT_0142 = IT_0140 + IT_0141;
    const ccomplex_t IT_0143 = conj(IT_0019) + conj(IT_0026);
    const ccomplex_t IT_0144 = (-18)*IT_0020;
    const ccomplex_t IT_0145 = (-18)*IT_0025;
    const ccomplex_t IT_0146 = m_u*IT_0088;
    const ccomplex_t IT_0147 = m_u*IT_0104;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0149 = IT_0075*IT_0111;
    const ccomplex_t IT_0150 = IT_0086*IT_0148*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0152 = IT_0094*IT_0118;
    const ccomplex_t IT_0153 = IT_0102*IT_0151*IT_0152;
    const ccomplex_t IT_0154 = IT_0146 + IT_0147 + -IT_0150 + -IT_0153;
    const ccomplex_t IT_0155 = (-6)*IT_0154;
    const ccomplex_t IT_0156 = IT_0144 + IT_0145 + IT_0155;
    const ccomplex_t IT_0157 = IT_0019 + IT_0026;
    const ccomplex_t IT_0158 = (-18)*conj(IT_0020);
    const ccomplex_t IT_0159 = (-18)*conj(IT_0025);
    const ccomplex_t IT_0160 = (-6)*conj(IT_0154);
    const ccomplex_t IT_0161 = IT_0158 + IT_0159 + IT_0160;
    const ccomplex_t IT_0162 = IT_0084*IT_0109;
    const ccomplex_t IT_0163 = IT_0086*IT_0148*IT_0162;
    const ccomplex_t IT_0164 = IT_0100*IT_0116;
    const ccomplex_t IT_0165 = IT_0102*IT_0151*IT_0164;
    const ccomplex_t IT_0166 = -IT_0163 + -IT_0165;
    const ccomplex_t IT_0167 = IT_0114 + IT_0121;
    const ccomplex_t IT_0168 = m_u*IT_0167;
    const ccomplex_t IT_0169 = conj(IT_0166) + conj(IT_0168);
    const ccomplex_t IT_0170 = (-18)*IT_0154;
    const ccomplex_t IT_0171 = (-6)*IT_0020;
    const ccomplex_t IT_0172 = (-6)*IT_0025;
    const ccomplex_t IT_0173 = IT_0170 + IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = (-18)*conj(IT_0154);
    const ccomplex_t IT_0175 = (-6)*conj(IT_0020);
    const ccomplex_t IT_0176 = (-6)*conj(IT_0025);
    const ccomplex_t IT_0177 = IT_0174 + IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = IT_0020 + IT_0025;
    const ccomplex_t IT_0179 = conj(IT_0020) + conj(IT_0025);
    const ccomplex_t IT_0180 = 6*s_34;
    const ccomplex_t IT_0181 = conj(IT_0166)*IT_0180;
    const ccomplex_t IT_0182 = s_34 + IT_0064;
    const ccomplex_t IT_0183 = m_b*e_em*mu_h*IT_0028*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0184 = IT_0005*IT_0183;
    const ccomplex_t IT_0185 = m_b*e_em*mu_h*IT_0028*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0186 = IT_0005*IT_0185;
    const ccomplex_t IT_0187 = A_b*IT_0034*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0188 = IT_0030*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = A_b*IT_0034*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0191 = IT_0030*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = (-0.5)*IT_0184 + 0.5*IT_0186 + (-0.5)*IT_0189 +
       0.5*IT_0192;
    const ccomplex_t IT_0194 = m_u*e_em*IT_0005*IT_0028*IT_0031*IT_0037;
    const ccomplex_t IT_0195 = (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = IT_0193*IT_0195;
    const ccomplex_t IT_0197 = cpow((-2)*s_12 + (-2)*IT_0010 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0198 = IT_0196*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = pow(m_Z, 2);
    const ccomplex_t IT_0201 = cpow((-2)*s_12 + (-2)*IT_0010 + IT_0200 + 
      -reg_prop, -1);
    const ccomplex_t IT_0202 = m_u*e_em*IT_0005*IT_0028;
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = IT_0005*IT_0030*IT_0038;
    const ccomplex_t IT_0205 = IT_0183*IT_0204;
    const ccomplex_t IT_0206 = IT_0185*IT_0204;
    const ccomplex_t IT_0207 = IT_0037*IT_0187;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = IT_0037*IT_0190;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = (-0.5)*IT_0205 + 0.5*IT_0206 + 0.5*IT_0208 + (
      -0.5)*IT_0210;
    const ccomplex_t IT_0212 = IT_0203*IT_0211;
    const ccomplex_t IT_0213 = IT_0201*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = -IT_0199 + -IT_0214;
    const ccomplex_t IT_0216 = 6*conj(IT_0215);
    const ccomplex_t IT_0217 = s_13 + s_14;
    const ccomplex_t IT_0218 = m_u*IT_0217;
    const ccomplex_t IT_0219 = 6*IT_0215;
    const ccomplex_t IT_0220 = -IT_0064;
    const ccomplex_t IT_0221 = 6*IT_0062;
    const ccomplex_t IT_0222 = 6*conj(IT_0062);
    const ccomplex_t IT_0223 = IT_0063*IT_0166;
    const ccomplex_t IT_0224 = IT_0166*IT_0180;
    const ccomplex_t IT_0225 = IT_0026*IT_0063;
    const ccomplex_t IT_0226 = (-18)*s_34;
    const ccomplex_t IT_0227 = IT_0065 + IT_0226;
    const ccomplex_t IT_0228 = 18*IT_0064;
    const ccomplex_t IT_0229 = IT_0063 + IT_0228;
    const ccomplex_t IT_0230 = 6*IT_0166;
    const ccomplex_t IT_0231 = m_u*s_14;
    const ccomplex_t IT_0232 = (-18)*IT_0026;
    const ccomplex_t IT_0233 = (-18)*IT_0019;
    const ccomplex_t IT_0234 = (-6)*IT_0166;
    const ccomplex_t IT_0235 = (-18)*IT_0166;
    const ccomplex_t IT_0236 = (-6)*IT_0026;
    const ccomplex_t IT_0237 = (-6)*IT_0019;
    const ccomplex_t IT_0238 = (-0.166666666666667)*conj(IT_0124)*IT_0156 + (
      -0.166666666666667)*IT_0124*IT_0161 + 1*IT_0105*(conj(IT_0019) + conj
      (IT_0026) + 3*conj(IT_0166)) + 3*IT_0107*(conj(IT_0019) + conj(IT_0026) +
       1./3*conj(IT_0166)) + (-0.166666666666667)*IT_0142*IT_0168 + (
      -0.166666666666667)*IT_0139*conj(IT_0168) + (-0.166666666666667)*conj
      (IT_0122)*IT_0173 + (-0.166666666666667)*IT_0122*IT_0177 + (
      -0.166666666666667)*conj(IT_0107)*(IT_0232 + IT_0233 + IT_0234) + (
      -0.166666666666667)*conj(IT_0105)*(IT_0235 + IT_0236 + IT_0237);
    const ccomplex_t IT_0239 = -s_14;
    const ccomplex_t IT_0240 = s_13 + IT_0239;
    const ccomplex_t IT_0241 = m_u*IT_0240;
    const ccomplex_t IT_0242 = conj(IT_0105) + conj(IT_0122);
    const ccomplex_t IT_0243 = IT_0105 + IT_0122;
    const ccomplex_t IT_0244 = conj(IT_0107) + conj(IT_0124);
    const ccomplex_t IT_0245 = 18*IT_0062;
    const ccomplex_t IT_0246 = IT_0107 + IT_0124;
    const ccomplex_t IT_0247 = 18*conj(IT_0062);
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0249 = (-2)*IT_0248;
    const ccomplex_t IT_0250 = cpow(s_12 + IT_0010 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0251 = IT_0248*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = IT_0249*IT_0252;
    const ccomplex_t IT_0254 = (-0.0833333333333333)*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0256 = 0.666666666666667*IT_0255;
    const ccomplex_t IT_0257 = cpow(IT_0256, 2);
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = IT_0250*IT_0258;
    const ccomplex_t IT_0260 = 0.5*IT_0259;
    const ccomplex_t IT_0261 = IT_0254 + IT_0260;
    const ccomplex_t IT_0262 = 6*conj(IT_0166);
    const ccomplex_t IT_0263 = 6*conj(IT_0154);
    const ccomplex_t IT_0264 = 18*conj(IT_0020);
    const ccomplex_t IT_0265 = 18*conj(IT_0026);
    const ccomplex_t IT_0266 = 18*conj(IT_0025);
    const ccomplex_t IT_0267 = 18*conj(IT_0019);
    const ccomplex_t IT_0268 = 36*conj(IT_0062);
    const ccomplex_t IT_0269 = 0.25*IT_0253;
    const ccomplex_t IT_0270 = 6*conj(IT_0020);
    const ccomplex_t IT_0271 = 6*conj(IT_0026);
    const ccomplex_t IT_0272 = 6*conj(IT_0025);
    const ccomplex_t IT_0273 = 6*conj(IT_0019);
    const ccomplex_t IT_0274 = 18*conj(IT_0166);
    const ccomplex_t IT_0275 = 18*conj(IT_0154);
    const ccomplex_t IT_0276 = IT_0269*((-12)*IT_0062 + 12*conj(IT_0062) + (
      -18)*IT_0168 + 18*conj(IT_0168) + IT_0170 + IT_0171 + IT_0172 + IT_0235 +
       IT_0236 + IT_0237 + IT_0270 + IT_0271 + IT_0272 + IT_0273 + IT_0274 +
       IT_0275);
    const ccomplex_t IT_0277 = IT_0221*IT_0242 + IT_0222*IT_0243 + IT_0244
      *IT_0245 + IT_0246*IT_0247 + IT_0261*((-36)*IT_0062 + IT_0144 + IT_0145 +
       IT_0155 + (-6)*IT_0168 + 6*conj(IT_0168) + IT_0232 + IT_0233 + IT_0234 +
       IT_0262 + IT_0263 + IT_0264 + IT_0265 + IT_0266 + IT_0267 + IT_0268) +
       IT_0276;
    const ccomplex_t IT_0278 = (-12)*IT_0125;
    const ccomplex_t IT_0279 = (-12)*IT_0133;
    const ccomplex_t IT_0280 = 24*IT_0127;
    const ccomplex_t IT_0281 = IT_0278 + IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = (-2)*IT_0254;
    const ccomplex_t IT_0283 = (-36)*IT_0125;
    const ccomplex_t IT_0284 = (-36)*IT_0133;
    const ccomplex_t IT_0285 = 72*IT_0127;
    const ccomplex_t IT_0286 = IT_0283 + IT_0284 + IT_0285;
    const ccomplex_t IT_0287 = cpow(IT_0254, 2);
    const ccomplex_t IT_0288 = cpow(IT_0269, 2);
    const ccomplex_t IT_0289 = (-18)*IT_0133;
    const ccomplex_t IT_0290 = IT_0130 + IT_0131 + IT_0289;
    const ccomplex_t IT_0291 = -IT_0254;
    const ccomplex_t IT_0292 = -IT_0260;
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = -IT_0269;
    const ccomplex_t IT_0295 = (-6)*IT_0133;
    const ccomplex_t IT_0296 = IT_0126 + IT_0128 + IT_0295;
    const ccomplex_t IT_0297 = 36*s_34;
    const ccomplex_t IT_0298 = m_u*s_13;
    const ccomplex_t IT_0299 = conj(IT_0105)*(IT_0020 + IT_0025 + 3*IT_0154) +
       IT_0135*(IT_0020 + IT_0025 + 0.333333333333333*IT_0154) + 3*conj(IT_0122)
      *(0.333333333333333*IT_0019 + 0.333333333333333*IT_0026 + IT_0166 +
       IT_0168) + 3*conj(IT_0124)*(IT_0157 + 1./3*IT_0168 + 0.0555555555555555
      *IT_0230) + 0.166666666666667*IT_0107*(IT_0263 + IT_0264 + IT_0266) +
       IT_0124*(conj(IT_0168) + 0.166666666666667*IT_0262 + 0.166666666666667
      *IT_0265 + 0.166666666666667*IT_0267) + 3*IT_0122*(conj(IT_0168) +
       0.0555555555555556*IT_0271 + 0.0555555555555556*IT_0273 +
       0.0555555555555556*IT_0274) + 0.166666666666667*IT_0105*(IT_0270 +
       IT_0272 + IT_0275);
    const ccomplex_t IT_0300 = ((conj(IT_0019) + conj(IT_0020) + conj(IT_0025)
       + conj(IT_0026))*IT_0062 + (IT_0019 + IT_0020 + IT_0025 + IT_0026)*conj
      (IT_0062))*IT_0066 + (conj(IT_0105)*IT_0107 + IT_0105*conj(IT_0107) + conj
      (IT_0122)*IT_0124 + IT_0122*conj(IT_0124))*IT_0129 + (IT_0105*conj(IT_0105
      ) + IT_0107*conj(IT_0107) + IT_0122*conj(IT_0122) + IT_0124*conj(IT_0124))
      *IT_0132 + (-6)*IT_0133*(conj(IT_0124)*IT_0134 + IT_0124*IT_0136 + (
      -0.166666666666667)*conj(IT_0122)*IT_0139 + (-0.166666666666667)*IT_0122
      *IT_0142) + IT_0064*(IT_0143*IT_0156 + IT_0157*IT_0161 + IT_0169*IT_0173 +
       (IT_0166 + IT_0168)*IT_0177) + IT_0063*(IT_0019*IT_0143 + IT_0154*conj
      (IT_0154) + IT_0168*IT_0169 + IT_0178*IT_0179) + (IT_0143*IT_0168 +
       IT_0157*conj(IT_0168) + conj(IT_0154)*IT_0178 + IT_0154*IT_0179)*IT_0180 
      + IT_0157*IT_0181 + IT_0168*IT_0182*IT_0216 + IT_0218*(IT_0122*IT_0216 +
       conj(IT_0122)*IT_0219) + (s_34 + IT_0220)*((conj(IT_0154) + conj(IT_0166)
       + conj(IT_0168))*IT_0221 + (IT_0154 + IT_0166 + IT_0168)*IT_0222) + conj
      (IT_0166)*IT_0223 + conj(IT_0168)*(IT_0182*IT_0219 + IT_0223) + IT_0143*
      (IT_0224 + IT_0225) + conj(IT_0215)*(18*(IT_0124 + (-0.333333333333333)
      *IT_0134)*IT_0218 + IT_0178*IT_0227 + IT_0157*IT_0229 + IT_0182*(IT_0155 +
       IT_0230)) + (-6)*IT_0231*IT_0238 + IT_0241*IT_0277 + (-2)*IT_0269*IT_0281
      *(IT_0260 + (-0.5)*IT_0282) + -IT_0286*(IT_0260*(IT_0260 + -IT_0282) +
       IT_0287 + IT_0288) + IT_0290*(IT_0244*IT_0261 + IT_0242*IT_0269 + IT_0246
      *IT_0293 + IT_0243*IT_0294) + (IT_0242*IT_0261 + IT_0244*IT_0269 + IT_0243
      *IT_0293 + IT_0246*IT_0294)*IT_0296 + IT_0062*conj(IT_0062)*((-36)*IT_0064
       + IT_0297) + IT_0215*(18*(conj(IT_0124) + (-0.333333333333333)*IT_0136)
      *IT_0218 + IT_0179*IT_0227 + IT_0143*IT_0229 + IT_0182*(IT_0160 + IT_0262)
       + conj(IT_0215)*(36*IT_0064 + IT_0297)) + 6*IT_0298*IT_0299;
    return create_ccomplex_return(IT_0300);
}

