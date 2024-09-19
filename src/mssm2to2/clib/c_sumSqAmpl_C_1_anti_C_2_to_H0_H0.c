#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_H0_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_H0_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_H0, 2);
    const ccomplex_t IT_0001 = pow(m_C_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = U_d2*V_Wp1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = V_u1*U_Wm2;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (IT_0005 + IT_0008)*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = IT_0010*(IT_0016 + IT_0018);
    const ccomplex_t IT_0020 = 1.4142135623731*e_em*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = IT_0014*IT_0022;
    const ccomplex_t IT_0024 = IT_0002*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0027 = IT_0003*IT_0026;
    const ccomplex_t IT_0028 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0029 = IT_0006*IT_0028;
    const ccomplex_t IT_0030 = IT_0010*(IT_0027 + IT_0029);
    const ccomplex_t IT_0031 = 1.4142135623731*e_em*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = U_d2*V_Wp2;
    const ccomplex_t IT_0035 = IT_0003*IT_0034;
    const ccomplex_t IT_0036 = V_u2*U_Wm2;
    const ccomplex_t IT_0037 = IT_0006*IT_0036;
    const ccomplex_t IT_0038 = IT_0010*(IT_0035 + IT_0037);
    const ccomplex_t IT_0039 = 1.4142135623731*e_em*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0033*IT_0041;
    const ccomplex_t IT_0043 = pow(m_C_1, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0000 + IT_0043 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0042*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0023*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = -IT_0025 + IT_0046 + -IT_0049 + IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0055 = IT_0014*IT_0041;
    const ccomplex_t IT_0056 = IT_0047*IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0058 = U_d1*V_Wp1;
    const ccomplex_t IT_0059 = IT_0003*IT_0058;
    const ccomplex_t IT_0060 = V_u1*U_Wm1;
    const ccomplex_t IT_0061 = IT_0006*IT_0060;
    const ccomplex_t IT_0062 = IT_0010*(IT_0059 + IT_0061);
    const ccomplex_t IT_0063 = 1.4142135623731*e_em*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0014*IT_0065;
    const ccomplex_t IT_0067 = IT_0002*IT_0057*IT_0066;
    const ccomplex_t IT_0068 = IT_0044*IT_0054*IT_0055;
    const ccomplex_t IT_0069 = m_C_1*IT_0025;
    const ccomplex_t IT_0070 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0071 = IT_0003*IT_0070;
    const ccomplex_t IT_0072 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0073 = IT_0006*IT_0072;
    const ccomplex_t IT_0074 = IT_0010*(IT_0071 + IT_0073);
    const ccomplex_t IT_0075 = 1.4142135623731*e_em*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0014*IT_0077;
    const ccomplex_t IT_0079 = IT_0044*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = m_C_2*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0001 + -IT_0043 + -reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0003*IT_0007;
    const ccomplex_t IT_0084 = IT_0004*IT_0006;
    const ccomplex_t IT_0085 = IT_0010*(IT_0083 + -IT_0084);
    const ccomplex_t IT_0086 = 1.4142135623731*e_em*IT_0085;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cos(theta_W);
    const ccomplex_t IT_0090 = cpow(IT_0089, -2);
    const ccomplex_t IT_0091 = IT_0009*IT_0090;
    const ccomplex_t IT_0092 = IT_0010 + IT_0091;
    const ccomplex_t IT_0093 = cos(beta);
    const ccomplex_t IT_0094 = cpow(IT_0006, 3);
    const ccomplex_t IT_0095 = sin(beta);
    const ccomplex_t IT_0096 = cpow(IT_0006, 2);
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_W*e_em*IT_0092*(IT_0093
      *IT_0094 + 5*IT_0003*((-1)*IT_0003*(IT_0006*IT_0093 + 1./5*IT_0003*IT_0095
      ) + IT_0095*IT_0096));
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0088*IT_0098;
    const ccomplex_t IT_0100 = IT_0082*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = IT_0050*IT_0057*IT_0066;
    const ccomplex_t IT_0103 = IT_0033*IT_0065;
    const ccomplex_t IT_0104 = IT_0050*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = m_C_2*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_12 + -IT_0001 + -IT_0043 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0108 = cpow(IT_0003, 2);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_W*e_em*IT_0092*(IT_0094
      *IT_0095 + -IT_0003*(IT_0003*IT_0006*IT_0095 + IT_0093*(IT_0096 + -IT_0108
      )));
    const ccomplex_t IT_0110 = (-1.5)*IT_0109;
    const ccomplex_t IT_0111 = IT_0014*IT_0110;
    const ccomplex_t IT_0112 = IT_0107*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = m_C_1*IT_0049;
    const ccomplex_t IT_0115 = IT_0056 + IT_0067 + IT_0068 + IT_0069 + IT_0081
       + IT_0101 + IT_0102 + IT_0106 + IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = s_24*m_C_1;
    const ccomplex_t IT_0117 = 2*IT_0116;
    const ccomplex_t IT_0118 = IT_0022*IT_0033;
    const ccomplex_t IT_0119 = IT_0050*IT_0057*IT_0118;
    const ccomplex_t IT_0120 = IT_0002*IT_0057*IT_0118;
    const ccomplex_t IT_0121 = IT_0033*IT_0077;
    const ccomplex_t IT_0122 = IT_0044*IT_0054*IT_0121;
    const ccomplex_t IT_0123 = m_C_2*IT_0046;
    const ccomplex_t IT_0124 = IT_0003*IT_0028;
    const ccomplex_t IT_0125 = IT_0006*IT_0026;
    const ccomplex_t IT_0126 = IT_0010*(IT_0124 + -IT_0125);
    const ccomplex_t IT_0127 = 1.4142135623731*e_em*IT_0126;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0098*IT_0129;
    const ccomplex_t IT_0131 = IT_0082*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = m_C_2*IT_0052;
    const ccomplex_t IT_0134 = IT_0047*IT_0078;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = m_C_1*IT_0135;
    const ccomplex_t IT_0137 = IT_0047*IT_0054*IT_0121;
    const ccomplex_t IT_0138 = IT_0033*IT_0110;
    const ccomplex_t IT_0139 = IT_0107*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = IT_0002*IT_0103;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = m_C_1*IT_0142;
    const ccomplex_t IT_0144 = IT_0119 + IT_0120 + IT_0122 + IT_0123 + IT_0132
       + IT_0133 + IT_0136 + IT_0137 + IT_0140 + IT_0143;
    const ccomplex_t IT_0145 = s_14*m_C_2;
    const ccomplex_t IT_0146 = (-2)*IT_0145;
    const ccomplex_t IT_0147 = s_12*IT_0000;
    const ccomplex_t IT_0148 = s_14*s_24;
    const ccomplex_t IT_0149 = (-2)*IT_0148;
    const ccomplex_t IT_0150 = IT_0147 + IT_0149;
    const ccomplex_t IT_0151 = (-2)*IT_0150;
    const ccomplex_t IT_0152 = IT_0080 + IT_0105 + -IT_0135 + -IT_0142;
    const ccomplex_t IT_0153 = m_C_1*m_C_2*IT_0000;
    const ccomplex_t IT_0154 = (-2)*IT_0153;
    const ccomplex_t IT_0155 = m_C_1*m_C_2;
    const ccomplex_t IT_0156 = (-2)*IT_0155;
    const ccomplex_t IT_0157 = 0.5*conj(IT_0152);
    const ccomplex_t IT_0158 = IT_0152*(IT_0117*conj(IT_0144) + conj(IT_0115)
      *IT_0146 + IT_0151*conj(IT_0152) + conj(IT_0053)*IT_0154) + IT_0053*(conj
      (IT_0115)*IT_0117 + conj(IT_0144)*IT_0146 + conj(IT_0053)*IT_0151 + conj
      (IT_0152)*IT_0154) + 2*IT_0144*(s_12*conj(IT_0144) + 0.5*conj(IT_0053)
      *IT_0146 + 0.5*conj(IT_0115)*IT_0156 + IT_0117*IT_0157) + 2*IT_0115*(s_12
      *conj(IT_0115) + 0.5*conj(IT_0053)*IT_0117 + 0.5*conj(IT_0144)*IT_0156 +
       IT_0146*IT_0157);
    return create_ccomplex_return(IT_0158);
}

