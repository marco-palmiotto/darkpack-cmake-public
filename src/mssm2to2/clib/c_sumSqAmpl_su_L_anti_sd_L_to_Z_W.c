#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_L_to_Z_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_L_to_Z_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*V_ud*IT_0000
      *IT_0002;
    const ccomplex_t IT_0004 = (-0.166666666666667)*IT_0003;
    const ccomplex_t IT_0005 = pow(m_sd_L, 2);
    const ccomplex_t IT_0006 = pow(m_su_L, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0005 + -IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0010;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0005 + 0.5*IT_0006;
    const ccomplex_t IT_0016 = s_12 + IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = 2*IT_0013;
    const ccomplex_t IT_0019 = IT_0012*IT_0018;
    const ccomplex_t IT_0020 = s_13 + s_23;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = -IT_0011;
    const ccomplex_t IT_0023 = IT_0018*IT_0022;
    const ccomplex_t IT_0024 = s_23*IT_0023;
    const ccomplex_t IT_0025 = IT_0013*IT_0022;
    const ccomplex_t IT_0026 = s_12 + IT_0005;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (-0.5)*IT_0021 + 0.5*IT_0024 + (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0017 + IT_0028;
    const ccomplex_t IT_0030 = IT_0008*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + -IT_0005 + -IT_0006 + IT_0031 
      + -reg_prop, -1);
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = pow(m_d, 2);
    const ccomplex_t IT_0035 = pow(m_u, 2);
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = cpow(IT_0036, 2);
    const ccomplex_t IT_0038 = sin(beta);
    const ccomplex_t IT_0039 = cpow(IT_0038, 2);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0010*(IT_0033*(IT_0034 + -IT_0035) + -m_W*(IT_0037 + -IT_0039));
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*m_W*e_em*IT_0002*IT_0009;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = IT_0032*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0004 + (-2)*IT_0030 + -IT_0046;
    const ccomplex_t IT_0048 = pow(m_W, -2);
    const ccomplex_t IT_0049 = pow(m_Z, -2);
    const ccomplex_t IT_0050 = pow(s_34, 2);
    const ccomplex_t IT_0051 = IT_0048*IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (-2)*IT_0013;
    const ccomplex_t IT_0053 = IT_0012*IT_0052;
    const ccomplex_t IT_0054 = IT_0008*(IT_0014 + IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = pow(s_13, 2);
    const ccomplex_t IT_0057 = IT_0049*IT_0056;
    const ccomplex_t IT_0058 = s_13*s_14*s_34*IT_0048*IT_0049;
    const ccomplex_t IT_0059 = pow(s_14, 2);
    const ccomplex_t IT_0060 = IT_0048*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0006 + -IT_0057 + IT_0058 + IT_0061;
    const ccomplex_t IT_0063 = 3*conj(IT_0047);
    const ccomplex_t IT_0064 = 3*IT_0047;
    const ccomplex_t IT_0065 = (IT_0006 + -IT_0057)*(IT_0006 + -IT_0060);
    const ccomplex_t IT_0066 = 3*IT_0055;
    const ccomplex_t IT_0067 = IT_0022*IT_0052;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0014 + IT_0053;
    const ccomplex_t IT_0070 = IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = IT_0008*IT_0070;
    const ccomplex_t IT_0072 = pow(m_Z, 2);
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0072 + (0 + _Complex_I*1)
      *m_su_L*Gamma_ul + reg_prop, -1);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0009 + (-3)
      *IT_0001*IT_0010);
    const ccomplex_t IT_0076 = 0.333333333333333*IT_0075;
    const ccomplex_t IT_0077 = IT_0022*IT_0076;
    const ccomplex_t IT_0078 = IT_0074*IT_0077;
    const ccomplex_t IT_0079 = -IT_0071 + -IT_0078;
    const ccomplex_t IT_0080 = s_13*s_24*s_34*IT_0048*IT_0049;
    const ccomplex_t IT_0081 = s_14*s_24*IT_0048;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = s_13*s_23*IT_0049;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = s_12 + IT_0082 + IT_0084;
    const ccomplex_t IT_0086 = IT_0080 + IT_0085;
    const ccomplex_t IT_0087 = pow(s_24, 2);
    const ccomplex_t IT_0088 = IT_0048*IT_0087;
    const ccomplex_t IT_0089 = (IT_0006 + -IT_0057)*(IT_0005 + -IT_0088);
    const ccomplex_t IT_0090 = 3*IT_0079;
    const ccomplex_t IT_0091 = (IT_0006 + -IT_0057)*(s_12 + -IT_0081);
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = (-3)*IT_0055;
    const ccomplex_t IT_0094 = IT_0008*(IT_0014 + -IT_0025 + IT_0053);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = s_14*s_34*IT_0048;
    const ccomplex_t IT_0097 = (s_12 + -IT_0083)*(s_13 + -IT_0096);
    const ccomplex_t IT_0098 = -IT_0013;
    const ccomplex_t IT_0099 = IT_0012*IT_0098;
    const ccomplex_t IT_0100 = IT_0008*IT_0099;
    const ccomplex_t IT_0101 = 3*IT_0100;
    const ccomplex_t IT_0102 = s_14*s_23*s_34*IT_0048*IT_0049;
    const ccomplex_t IT_0103 = IT_0085 + IT_0102;
    const ccomplex_t IT_0104 = (s_12 + -IT_0081)*(s_12 + -IT_0083);
    const ccomplex_t IT_0105 = (IT_0006 + -IT_0060)*(s_12 + -IT_0083);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = (-3)*conj(IT_0055);
    const ccomplex_t IT_0108 = -IT_0025 + -IT_0067;
    const ccomplex_t IT_0109 = IT_0069 + IT_0108;
    const ccomplex_t IT_0110 = IT_0008*IT_0109;
    const ccomplex_t IT_0111 = cpow((-2)*s_23 + IT_0072 + (0 + _Complex_I*1)
      *m_sd_L*Gamma_dl + reg_prop, -1);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0009 + 3
      *IT_0001*IT_0010);
    const ccomplex_t IT_0114 = 0.333333333333333*IT_0113;
    const ccomplex_t IT_0115 = IT_0022*IT_0114;
    const ccomplex_t IT_0116 = IT_0112*IT_0115;
    const ccomplex_t IT_0117 = -IT_0110 + IT_0116;
    const ccomplex_t IT_0118 = s_24*s_34*IT_0048;
    const ccomplex_t IT_0119 = (s_12 + -IT_0083)*(s_23 + -IT_0118);
    const ccomplex_t IT_0120 = s_23*s_24*s_34*IT_0048*IT_0049;
    const ccomplex_t IT_0121 = -IT_0088;
    const ccomplex_t IT_0122 = pow(s_23, 2);
    const ccomplex_t IT_0123 = IT_0049*IT_0122;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0005 + IT_0120 + IT_0121 + IT_0124;
    const ccomplex_t IT_0126 = (s_12 + -IT_0083)*(IT_0005 + -IT_0088);
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = (-3)*IT_0079;
    const ccomplex_t IT_0129 = (s_12 + -IT_0081)*(IT_0005 + -IT_0123);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = (-3)*IT_0095;
    const ccomplex_t IT_0132 = (IT_0005 + -IT_0088)*(IT_0005 + -IT_0123);
    const ccomplex_t IT_0133 = IT_0022*IT_0098;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = IT_0099 + IT_0134;
    const ccomplex_t IT_0136 = IT_0008*IT_0135;
    const ccomplex_t IT_0137 = -IT_0116 + IT_0136;
    const ccomplex_t IT_0138 = (s_23 + -IT_0118)*(IT_0005 + -IT_0123);
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = 3*conj(IT_0100);
    const ccomplex_t IT_0141 = (-3)*conj(IT_0079);
    const ccomplex_t IT_0142 = (-3)*conj(IT_0095);
    const ccomplex_t IT_0143 = 3*conj(IT_0055);
    const ccomplex_t IT_0144 = 0.333333333333333*IT_0104;
    const ccomplex_t IT_0145 = 3*IT_0117*(0.333333333333333*IT_0063*IT_0125 +
       conj(IT_0117)*IT_0132 + -conj(IT_0137)*IT_0139 + 0.333333333333333
      *IT_0119*IT_0140 + 0.333333333333333*IT_0127*IT_0141 + 0.333333333333333
      *IT_0130*IT_0142 + IT_0143*IT_0144);
    const ccomplex_t IT_0146 = (IT_0006 + -IT_0060)*(IT_0005 + -IT_0123);
    const ccomplex_t IT_0147 = 3*conj(IT_0079);
    const ccomplex_t IT_0148 = s_34*IT_0048;
    const ccomplex_t IT_0149 = s_13*s_34*IT_0049;
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = s_14 + IT_0150;
    const ccomplex_t IT_0152 = IT_0148*IT_0151;
    const ccomplex_t IT_0153 = (IT_0006 + -IT_0057)*(s_13 + -IT_0096);
    const ccomplex_t IT_0154 = -IT_0153;
    const ccomplex_t IT_0155 = (IT_0006 + -IT_0057)*(s_23 + -IT_0118);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = (-0.333333333333333)*IT_0156;
    const ccomplex_t IT_0158 = IT_0048*IT_0050;
    const ccomplex_t IT_0159 = (IT_0006 + -IT_0057)*(IT_0072 + -IT_0158);
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = (-3)*conj(IT_0047);
    const ccomplex_t IT_0162 = s_23*s_34*IT_0049;
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = s_24 + IT_0163;
    const ccomplex_t IT_0165 = IT_0148*IT_0164;
    const ccomplex_t IT_0166 = (-0.333333333333333)*IT_0097;
    const ccomplex_t IT_0167 = (-0.333333333333333)*IT_0119;
    const ccomplex_t IT_0168 = (s_12 + -IT_0083)*(IT_0072 + -IT_0158);
    const ccomplex_t IT_0169 = (-0.333333333333333)*IT_0168;
    const ccomplex_t IT_0170 = (s_13 + -IT_0096)*(IT_0005 + -IT_0123);
    const ccomplex_t IT_0171 = -IT_0170;
    const ccomplex_t IT_0172 = (-0.333333333333333)*IT_0171;
    const ccomplex_t IT_0173 = (IT_0005 + -IT_0123)*(IT_0072 + -IT_0158);
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = 6*IT_0047*conj(IT_0047)*(1 + 0.5*IT_0051) +
       IT_0055*IT_0062*IT_0063 + conj(IT_0055)*(IT_0062*IT_0064 + IT_0065
      *IT_0066) + conj(IT_0079)*(IT_0064*IT_0086 + IT_0089*IT_0090 + IT_0092
      *IT_0093) + conj(IT_0095)*(IT_0097*IT_0101 + IT_0064*IT_0103 + IT_0090
      *IT_0104 + IT_0093*IT_0106) + IT_0079*(IT_0063*IT_0086 + IT_0092*IT_0107) 
      + conj(IT_0117)*(IT_0066*IT_0104 + IT_0101*IT_0119 + IT_0064*IT_0125 +
       IT_0127*IT_0128 + IT_0130*IT_0131) + IT_0145 + 3*IT_0095*
      (0.333333333333333*IT_0063*IT_0103 + 0.333333333333333*IT_0106*IT_0107 +
       0.333333333333333*IT_0097*IT_0140 + conj(IT_0095)*IT_0146 + IT_0144
      *IT_0147) + (-3)*conj(IT_0100)*(IT_0047*IT_0152 + (-0.333333333333333)
      *IT_0093*IT_0154 + IT_0128*IT_0157) + (-3)*IT_0100*((-0.333333333333333)
      *IT_0107*IT_0154 + IT_0141*IT_0157 + conj(IT_0100)*IT_0160 + (
      -0.333333333333333)*IT_0152*IT_0161) + (-3)*conj(IT_0137)*(IT_0047*IT_0165
       + IT_0066*IT_0166 + IT_0090*IT_0167 + IT_0101*IT_0169 + IT_0131*IT_0172) 
      + (-3)*IT_0137*(conj(IT_0117)*IT_0139 + (-0.333333333333333)*IT_0161
      *IT_0165 + IT_0143*IT_0166 + IT_0147*IT_0167 + IT_0140*IT_0169 + IT_0142
      *IT_0172 + conj(IT_0137)*IT_0174);
    return create_ccomplex_return(IT_0175);
}

