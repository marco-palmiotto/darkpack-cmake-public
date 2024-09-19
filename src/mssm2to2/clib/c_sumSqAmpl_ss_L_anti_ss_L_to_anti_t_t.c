#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_L_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_L_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = pow(m_s, 2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(alpha);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0004*IT_0006 + (-0.5)*m_W*(IT_0002*IT_0004 + IT_0007*IT_0008)*
      (IT_0006 + 1./3*IT_0005*IT_0010));
    const ccomplex_t IT_0012 = cpow(IT_0008, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*m_t*e_em*IT_0000*IT_0006
      *IT_0007*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = pow(m_ss_L, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0006*IT_0007 + (-0.5)*m_W*(IT_0002*IT_0007 + -IT_0004*IT_0008)
      *(IT_0006 + 1./3*IT_0005*IT_0010));
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_t*e_em*IT_0000*IT_0004
      *IT_0006*IT_0012;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0019 + -IT_0027;
    const ccomplex_t IT_0029 = pow(m_t, 2);
    const ccomplex_t IT_0030 = cpow(IT_0009, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*3)*e_em*(IT_0006*IT_0009 + 1./3
      *IT_0005*IT_0030);
    const ccomplex_t IT_0032 = 0.166666666666667*IT_0031;
    const ccomplex_t IT_0033 = m_t*IT_0032;
    const ccomplex_t IT_0034 = 0.333333333333333*IT_0031;
    const ccomplex_t IT_0035 = m_t*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0030;
    const ccomplex_t IT_0037 = (-0.666666666666667)*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0005*IT_0030;
    const ccomplex_t IT_0042 = e_em*IT_0041;
    const ccomplex_t IT_0043 = IT_0006*IT_0009;
    const ccomplex_t IT_0044 = e_em*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0042 + (-3)*IT_0044);
    const ccomplex_t IT_0046 = (-0.166666666666667)*IT_0045;
    const ccomplex_t IT_0047 = IT_0038*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = (IT_0033 + -IT_0035)*(IT_0040 + -IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = (-18)*IT_0029;
    const ccomplex_t IT_0053 = 18*s_34;
    const ccomplex_t IT_0054 = IT_0052 + IT_0053;
    const ccomplex_t IT_0055 = IT_0034*IT_0040;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u1
      *e_em*IT_0000*IT_0006*IT_0012;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *conj(V_u1)*e_em*IT_0000*IT_0006*IT_0012;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = IT_0058*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0016 + IT_0029 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u2
      *e_em*IT_0000*IT_0006*IT_0012;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *conj(V_u2)*e_em*IT_0000*IT_0006*IT_0012;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_23 + IT_0016 + IT_0029 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = -IT_0064 + -IT_0072;
    const ccomplex_t IT_0074 = IT_0034*IT_0048;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm1)*IT_0006;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*conj(V_ts)*e_em*U_Wm1*IT_0006;
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = IT_0062*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm2)*IT_0006;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*conj(V_ts)*e_em*U_Wm2*IT_0006;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = IT_0070*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0080 + -IT_0085;
    const ccomplex_t IT_0087 = s_34*IT_0016;
    const ccomplex_t IT_0088 = (-6)*IT_0087;
    const ccomplex_t IT_0089 = s_13*s_14;
    const ccomplex_t IT_0090 = 12*IT_0089;
    const ccomplex_t IT_0091 = IT_0088 + IT_0090;
    const ccomplex_t IT_0092 = (-18)*IT_0087;
    const ccomplex_t IT_0093 = 36*IT_0089;
    const ccomplex_t IT_0094 = IT_0092 + IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0096 = (-2)*IT_0095;
    const ccomplex_t IT_0097 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0095*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = IT_0096*IT_0099;
    const ccomplex_t IT_0101 = 0.25*IT_0100;
    const ccomplex_t IT_0102 = cpow(IT_0101, 2);
    const ccomplex_t IT_0103 = (-36)*IT_0087;
    const ccomplex_t IT_0104 = IT_0016*IT_0029;
    const ccomplex_t IT_0105 = (-36)*IT_0104;
    const ccomplex_t IT_0106 = 72*IT_0089;
    const ccomplex_t IT_0107 = IT_0103 + IT_0105 + IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0109 = 0.666666666666667*IT_0108;
    const ccomplex_t IT_0110 = cpow(IT_0109, 2);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0097*IT_0111;
    const ccomplex_t IT_0113 = (-0.0833333333333333)*IT_0100 + 0.5*IT_0112;
    const ccomplex_t IT_0114 = (-18)*IT_0051;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0116 = IT_0058*IT_0077;
    const ccomplex_t IT_0117 = IT_0062*IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0119 = IT_0066*IT_0082;
    const ccomplex_t IT_0120 = IT_0070*IT_0118*IT_0119;
    const ccomplex_t IT_0121 = m_t*IT_0080;
    const ccomplex_t IT_0122 = m_t*IT_0085;
    const ccomplex_t IT_0123 = -IT_0117 + -IT_0120 + IT_0121 + IT_0122;
    const ccomplex_t IT_0124 = (-6)*IT_0123;
    const ccomplex_t IT_0125 = IT_0114 + IT_0124;
    const ccomplex_t IT_0126 = (-18)*conj(IT_0051);
    const ccomplex_t IT_0127 = (-6)*conj(IT_0123);
    const ccomplex_t IT_0128 = IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = IT_0060*IT_0076;
    const ccomplex_t IT_0130 = IT_0062*IT_0115*IT_0129;
    const ccomplex_t IT_0131 = IT_0068*IT_0081;
    const ccomplex_t IT_0132 = IT_0070*IT_0118*IT_0131;
    const ccomplex_t IT_0133 = m_t*IT_0064;
    const ccomplex_t IT_0134 = m_t*IT_0072;
    const ccomplex_t IT_0135 = -IT_0130 + -IT_0132 + IT_0133 + IT_0134;
    const ccomplex_t IT_0136 = (-18)*IT_0123;
    const ccomplex_t IT_0137 = (-6)*IT_0051;
    const ccomplex_t IT_0138 = IT_0136 + IT_0137;
    const ccomplex_t IT_0139 = (-18)*conj(IT_0123);
    const ccomplex_t IT_0140 = (-6)*conj(IT_0051);
    const ccomplex_t IT_0141 = IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = 6*s_34;
    const ccomplex_t IT_0143 = -IT_0029;
    const ccomplex_t IT_0144 = s_34 + IT_0143;
    const ccomplex_t IT_0145 = 6*IT_0028;
    const ccomplex_t IT_0146 = 6*conj(IT_0028);
    const ccomplex_t IT_0147 = (-6)*IT_0104;
    const ccomplex_t IT_0148 = IT_0088 + IT_0090 + IT_0147;
    const ccomplex_t IT_0149 = conj(IT_0056) + conj(IT_0075);
    const ccomplex_t IT_0150 = conj(IT_0073) + conj(IT_0086);
    const ccomplex_t IT_0151 = IT_0073 + IT_0086;
    const ccomplex_t IT_0152 = -IT_0113;
    const ccomplex_t IT_0153 = IT_0056 + IT_0075;
    const ccomplex_t IT_0154 = -IT_0101;
    const ccomplex_t IT_0155 = (-18)*IT_0104;
    const ccomplex_t IT_0156 = IT_0092 + IT_0093 + IT_0155;
    const ccomplex_t IT_0157 = m_t*s_13;
    const ccomplex_t IT_0158 = 3*conj(IT_0086);
    const ccomplex_t IT_0159 = conj(IT_0075) + IT_0158;
    const ccomplex_t IT_0160 = 6*conj(IT_0123);
    const ccomplex_t IT_0161 = 18*conj(IT_0051);
    const ccomplex_t IT_0162 = 18*conj(IT_0135);
    const ccomplex_t IT_0163 = 6*conj(IT_0051);
    const ccomplex_t IT_0164 = 18*conj(IT_0123);
    const ccomplex_t IT_0165 = 6*conj(IT_0135);
    const ccomplex_t IT_0166 = 18*conj(IT_0050);
    const ccomplex_t IT_0167 = 0.166666666666667*IT_0050*(18*conj(IT_0056) + 6
      *conj(IT_0073)) + 0.166666666666667*IT_0051*(18*conj(IT_0075) + 6*conj
      (IT_0086)) + (conj(IT_0056) + 3*conj(IT_0073))*IT_0135 + IT_0123*IT_0159 +
       0.166666666666667*IT_0075*(IT_0160 + IT_0161) + IT_0073*(conj(IT_0050) +
       0.166666666666667*IT_0162) + 0.166666666666667*IT_0086*(IT_0163 + IT_0164
      ) + 0.166666666666667*IT_0056*(IT_0165 + IT_0166);
    const ccomplex_t IT_0168 = -s_14;
    const ccomplex_t IT_0169 = s_13 + IT_0168;
    const ccomplex_t IT_0170 = m_t*IT_0169;
    const ccomplex_t IT_0171 = 18*IT_0028;
    const ccomplex_t IT_0172 = 18*conj(IT_0028);
    const ccomplex_t IT_0173 = (-18)*IT_0135;
    const ccomplex_t IT_0174 = (-6)*IT_0135;
    const ccomplex_t IT_0175 = IT_0145*IT_0150 + IT_0146*IT_0151 + IT_0149
      *IT_0171 + IT_0153*IT_0172 + IT_0101*((-12)*IT_0028 + 12*conj(IT_0028) + (
      -6)*IT_0050 + 6*conj(IT_0050) + IT_0136 + IT_0137 + IT_0162 + IT_0163 +
       IT_0164 + IT_0173) + IT_0113*((-36)*IT_0028 + 36*conj(IT_0028) + (-18)
      *IT_0050 + IT_0114 + IT_0124 + IT_0160 + IT_0161 + IT_0165 + IT_0166 +
       IT_0174);
    const ccomplex_t IT_0176 = m_t*s_14;
    const ccomplex_t IT_0177 = IT_0075 + 3*IT_0086;
    const ccomplex_t IT_0178 = (-18)*IT_0075;
    const ccomplex_t IT_0179 = (-6)*IT_0086;
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = (-18)*conj(IT_0075);
    const ccomplex_t IT_0182 = (-6)*conj(IT_0086);
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = (-0.166666666666667)*IT_0056;
    const ccomplex_t IT_0185 = (-0.166666666666667)*conj(IT_0056);
    const ccomplex_t IT_0186 = (-0.166666666666667)*conj(IT_0073)*IT_0138 + (
      -0.166666666666667)*IT_0073*IT_0141 + (-0.166666666666667)*conj(IT_0086)
      *IT_0173 + (-0.166666666666667)*conj(IT_0075)*IT_0174 + conj(IT_0135)
      *IT_0177 + (-0.166666666666667)*conj(IT_0050)*IT_0180 + (
      -0.166666666666667)*IT_0050*IT_0183 + IT_0128*IT_0184 + IT_0125*IT_0185;
    const ccomplex_t IT_0187 = IT_0028*conj(IT_0028)*(36*s_34 + (-36)*IT_0029)
       + (conj(IT_0028)*(IT_0050 + IT_0051) + IT_0028*(conj(IT_0050) + conj
      (IT_0051)))*IT_0054 + (conj(IT_0056)*IT_0073 + IT_0056*conj(IT_0073) +
       conj(IT_0075)*IT_0086 + IT_0075*conj(IT_0086))*IT_0091 + (IT_0056*conj
      (IT_0056) + IT_0073*conj(IT_0073) + IT_0075*conj(IT_0075) + IT_0086*conj
      (IT_0086))*IT_0094 + -IT_0102*IT_0107 + (-2)*IT_0113*(IT_0101*((-12)
      *IT_0087 + 24*IT_0089 + (-12)*IT_0104) + 0.5*IT_0107*IT_0113) + IT_0053*
      (IT_0050*conj(IT_0050) + IT_0051*conj(IT_0051) + IT_0123*conj(IT_0123) +
       IT_0135*conj(IT_0135)) + IT_0029*(conj(IT_0050)*IT_0125 + IT_0050*IT_0128
       + conj(IT_0135)*IT_0138 + IT_0135*IT_0141) + (conj(IT_0051)*IT_0123 +
       IT_0051*conj(IT_0123) + conj(IT_0050)*IT_0135 + IT_0050*conj(IT_0135))
      *IT_0142 + IT_0144*((conj(IT_0123) + conj(IT_0135))*IT_0145 + (IT_0123 +
       IT_0135)*IT_0146) + IT_0148*(IT_0101*IT_0149 + IT_0113*IT_0150 + IT_0151
      *IT_0152 + IT_0153*IT_0154) + (IT_0113*IT_0149 + IT_0101*IT_0150 + IT_0152
      *IT_0153 + IT_0151*IT_0154)*IT_0156 + 6*IT_0157*IT_0167 + IT_0170*IT_0175 
      + (-6)*IT_0104*(IT_0073*IT_0159 + conj(IT_0073)*IT_0177 + IT_0183*IT_0184 
      + IT_0180*IT_0185) + (-6)*IT_0176*IT_0186;
    return create_ccomplex_return(IT_0187);
}

