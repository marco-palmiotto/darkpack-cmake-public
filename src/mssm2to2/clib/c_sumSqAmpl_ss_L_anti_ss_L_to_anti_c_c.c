#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_L_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_L_to_anti_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
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
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = pow(m_s, 2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0004*IT_0006 + (-0.5)*m_W*(IT_0002*IT_0004 + -IT_0007*IT_0008)
      *(IT_0006 + 1./3*IT_0005*IT_0010));
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0007, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_c*e_em*IT_0000*IT_0006
      *IT_0008*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = pow(m_ss_L, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0006*IT_0008 + (-0.5)*m_W*(IT_0004*IT_0007 + IT_0002*IT_0008)*
      (IT_0006 + 1./3*IT_0005*IT_0010));
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_c*e_em*IT_0000*IT_0004
      *IT_0006*IT_0013;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -IT_0020 + -IT_0027;
    const ccomplex_t IT_0029 = pow(m_c, 2);
    const ccomplex_t IT_0030 = cpow(IT_0009, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*3)*e_em*(IT_0006*IT_0009 + 1./3
      *IT_0005*IT_0030);
    const ccomplex_t IT_0032 = 0.166666666666667*IT_0031;
    const ccomplex_t IT_0033 = m_c*IT_0032;
    const ccomplex_t IT_0034 = 0.333333333333333*IT_0031;
    const ccomplex_t IT_0035 = m_c*IT_0034;
    const ccomplex_t IT_0036 = IT_0005*IT_0030;
    const ccomplex_t IT_0037 = e_em*IT_0036;
    const ccomplex_t IT_0038 = IT_0006*IT_0009;
    const ccomplex_t IT_0039 = e_em*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0037 + (-3)*IT_0039);
    const ccomplex_t IT_0041 = (-0.166666666666667)*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*IT_0005*IT_0030;
    const ccomplex_t IT_0046 = (-0.666666666666667)*IT_0045;
    const ccomplex_t IT_0047 = IT_0042*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = (IT_0033 + -IT_0035)*(IT_0044 + -IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = (-18)*IT_0029;
    const ccomplex_t IT_0052 = 18*s_34;
    const ccomplex_t IT_0053 = IT_0051 + IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0006;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*conj(V_cs)*e_em*U_Wm2*IT_0006;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0017 + IT_0029 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0006;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*conj(V_cs)*e_em*U_Wm1*IT_0006;
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0017 + IT_0029 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0059 + -IT_0065;
    const ccomplex_t IT_0067 = IT_0034*IT_0044;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0034*IT_0048;
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0000*IT_0006*IT_0013;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u2)*e_em*IT_0000*IT_0006*IT_0013;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = IT_0057*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0000*IT_0006*IT_0013;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u1)*e_em*IT_0000*IT_0006*IT_0013;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = IT_0063*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = -IT_0077 + -IT_0084;
    const ccomplex_t IT_0086 = s_34*IT_0017;
    const ccomplex_t IT_0087 = (-6)*IT_0086;
    const ccomplex_t IT_0088 = s_13*s_14;
    const ccomplex_t IT_0089 = 12*IT_0088;
    const ccomplex_t IT_0090 = IT_0087 + IT_0089;
    const ccomplex_t IT_0091 = (-18)*IT_0086;
    const ccomplex_t IT_0092 = 36*IT_0088;
    const ccomplex_t IT_0093 = IT_0091 + IT_0092;
    const ccomplex_t IT_0094 = cpow(s_12 + IT_0017 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0096 = 0.666666666666667*IT_0095;
    const ccomplex_t IT_0097 = cpow(IT_0096, 2);
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0094*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0101 = (-2)*IT_0100;
    const ccomplex_t IT_0102 = IT_0094*IT_0100;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0099 + (-0.0833333333333333)*IT_0104;
    const ccomplex_t IT_0106 = cpow(IT_0105, 2);
    const ccomplex_t IT_0107 = (-36)*IT_0086;
    const ccomplex_t IT_0108 = IT_0017*IT_0029;
    const ccomplex_t IT_0109 = (-36)*IT_0108;
    const ccomplex_t IT_0110 = 72*IT_0088;
    const ccomplex_t IT_0111 = IT_0107 + IT_0109 + IT_0110;
    const ccomplex_t IT_0112 = 0.25*IT_0104;
    const ccomplex_t IT_0113 = IT_0054*IT_0074;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0115 = IT_0057*IT_0113*IT_0114;
    const ccomplex_t IT_0116 = m_c*IT_0077;
    const ccomplex_t IT_0117 = IT_0060*IT_0081;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0119 = IT_0063*IT_0117*IT_0118;
    const ccomplex_t IT_0120 = m_c*IT_0084;
    const ccomplex_t IT_0121 = -IT_0115 + IT_0116 + -IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = IT_0055*IT_0072;
    const ccomplex_t IT_0123 = IT_0057*IT_0114*IT_0122;
    const ccomplex_t IT_0124 = IT_0061*IT_0079;
    const ccomplex_t IT_0125 = IT_0063*IT_0118*IT_0124;
    const ccomplex_t IT_0126 = m_c*IT_0065;
    const ccomplex_t IT_0127 = m_c*IT_0059;
    const ccomplex_t IT_0128 = -IT_0123 + -IT_0125 + IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = (-18)*IT_0128;
    const ccomplex_t IT_0130 = (-6)*IT_0050;
    const ccomplex_t IT_0131 = IT_0129 + IT_0130;
    const ccomplex_t IT_0132 = (-18)*conj(IT_0128);
    const ccomplex_t IT_0133 = (-6)*conj(IT_0050);
    const ccomplex_t IT_0134 = IT_0132 + IT_0133;
    const ccomplex_t IT_0135 = (-18)*IT_0050;
    const ccomplex_t IT_0136 = (-6)*IT_0128;
    const ccomplex_t IT_0137 = IT_0135 + IT_0136;
    const ccomplex_t IT_0138 = (-18)*conj(IT_0050);
    const ccomplex_t IT_0139 = (-6)*conj(IT_0128);
    const ccomplex_t IT_0140 = IT_0138 + IT_0139;
    const ccomplex_t IT_0141 = (-18)*IT_0066 + (-6)*IT_0068;
    const ccomplex_t IT_0142 = (-18)*IT_0068;
    const ccomplex_t IT_0143 = (-6)*IT_0066;
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = (-18)*conj(IT_0068);
    const ccomplex_t IT_0146 = (-6)*conj(IT_0066);
    const ccomplex_t IT_0147 = IT_0145 + IT_0146;
    const ccomplex_t IT_0148 = 6*s_34;
    const ccomplex_t IT_0149 = -IT_0029;
    const ccomplex_t IT_0150 = s_34 + IT_0149;
    const ccomplex_t IT_0151 = 6*IT_0028;
    const ccomplex_t IT_0152 = 6*conj(IT_0028);
    const ccomplex_t IT_0153 = m_c*s_14;
    const ccomplex_t IT_0154 = (-18)*IT_0121;
    const ccomplex_t IT_0155 = (-6)*IT_0121;
    const ccomplex_t IT_0156 = (-6)*IT_0108;
    const ccomplex_t IT_0157 = IT_0087 + IT_0089 + IT_0156;
    const ccomplex_t IT_0158 = conj(IT_0066) + conj(IT_0085);
    const ccomplex_t IT_0159 = conj(IT_0068) + conj(IT_0070);
    const ccomplex_t IT_0160 = IT_0068 + IT_0070;
    const ccomplex_t IT_0161 = -IT_0112;
    const ccomplex_t IT_0162 = IT_0066 + IT_0085;
    const ccomplex_t IT_0163 = -IT_0105;
    const ccomplex_t IT_0164 = (-18)*IT_0108;
    const ccomplex_t IT_0165 = IT_0091 + IT_0092 + IT_0164;
    const ccomplex_t IT_0166 = m_c*s_13;
    const ccomplex_t IT_0167 = 6*conj(IT_0050);
    const ccomplex_t IT_0168 = 18*conj(IT_0128);
    const ccomplex_t IT_0169 = 6*conj(IT_0128);
    const ccomplex_t IT_0170 = 18*conj(IT_0050);
    const ccomplex_t IT_0171 = 18*conj(IT_0121);
    const ccomplex_t IT_0172 = 6*conj(IT_0121);
    const ccomplex_t IT_0173 = 18*conj(IT_0049);
    const ccomplex_t IT_0174 = IT_0050*(conj(IT_0066) + 3*conj(IT_0068)) +
       0.166666666666667*IT_0049*(18*conj(IT_0070) + 6*conj(IT_0085)) + (conj
      (IT_0070) + 3*conj(IT_0085))*IT_0121 + 0.166666666666667*(18*conj(IT_0066)
       + 6*conj(IT_0068))*IT_0128 + 0.166666666666667*IT_0066*(IT_0167 + IT_0168
      ) + 0.166666666666667*IT_0068*(IT_0169 + IT_0170) + IT_0085*(conj(IT_0049)
       + 0.166666666666667*IT_0171) + 0.166666666666667*IT_0070*(IT_0172 +
       IT_0173);
    const ccomplex_t IT_0175 = -s_14;
    const ccomplex_t IT_0176 = s_13 + IT_0175;
    const ccomplex_t IT_0177 = m_c*IT_0176;
    const ccomplex_t IT_0178 = 18*IT_0028;
    const ccomplex_t IT_0179 = 18*conj(IT_0028);
    const ccomplex_t IT_0180 = IT_0151*IT_0158 + IT_0152*IT_0162 + IT_0112*((
      -12)*IT_0028 + 12*conj(IT_0028) + (-6)*IT_0049 + 6*conj(IT_0049) + IT_0129
       + IT_0130 + IT_0154 + IT_0167 + IT_0168 + IT_0171) + IT_0105*((-36)
      *IT_0028 + 36*conj(IT_0028) + (-18)*IT_0049 + IT_0135 + IT_0136 + IT_0155 
      + IT_0169 + IT_0170 + IT_0172 + IT_0173) + IT_0159*IT_0178 + IT_0160
      *IT_0179;
    const ccomplex_t IT_0181 = IT_0028*conj(IT_0028)*(36*s_34 + (-36)*IT_0029)
       + (conj(IT_0028)*(IT_0049 + IT_0050) + IT_0028*(conj(IT_0049) + conj
      (IT_0050)))*IT_0053 + (conj(IT_0066)*IT_0068 + IT_0066*conj(IT_0068) +
       conj(IT_0070)*IT_0085 + IT_0070*conj(IT_0085))*IT_0090 + (IT_0066*conj
      (IT_0066) + IT_0068*conj(IT_0068) + IT_0070*conj(IT_0070) + IT_0085*conj
      (IT_0085))*IT_0093 + -IT_0106*IT_0111 + (-2)*IT_0112*(IT_0105*((-12)
      *IT_0086 + 24*IT_0088 + (-12)*IT_0108) + 0.5*IT_0111*IT_0112) + IT_0052*
      (IT_0049*conj(IT_0049) + IT_0050*conj(IT_0050) + IT_0121*conj(IT_0121) +
       IT_0128*conj(IT_0128)) + IT_0029*(conj(IT_0121)*IT_0131 + IT_0121*IT_0134
       + conj(IT_0049)*IT_0137 + IT_0049*IT_0140) + IT_0108*(((-18)*conj(IT_0066
      ) + (-6)*conj(IT_0068))*IT_0085 + conj(IT_0085)*IT_0141 + conj(IT_0070)
      *IT_0144 + IT_0070*IT_0147) + (conj(IT_0049)*IT_0121 + IT_0049*conj
      (IT_0121) + conj(IT_0050)*IT_0128 + IT_0050*conj(IT_0128))*IT_0148 +
       IT_0150*((conj(IT_0121) + conj(IT_0128))*IT_0151 + (IT_0121 + IT_0128)
      *IT_0152) + IT_0153*(conj(IT_0085)*IT_0131 + IT_0085*IT_0134 + conj
      (IT_0070)*IT_0137 + IT_0070*IT_0140 + conj(IT_0121)*IT_0141 + conj(IT_0049
      )*IT_0144 + IT_0049*IT_0147 + conj(IT_0066)*IT_0154 + conj(IT_0068)
      *IT_0155) + IT_0157*(IT_0105*IT_0158 + IT_0112*IT_0159 + IT_0160*IT_0161 +
       IT_0162*IT_0163) + (IT_0112*IT_0158 + IT_0105*IT_0159 + IT_0161*IT_0162 +
       IT_0160*IT_0163)*IT_0165 + 6*IT_0166*IT_0174 + IT_0177*IT_0180;
    return create_ccomplex_return(IT_0181);
}

