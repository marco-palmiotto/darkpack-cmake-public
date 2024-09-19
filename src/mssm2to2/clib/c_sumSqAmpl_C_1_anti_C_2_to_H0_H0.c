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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = U_d2*V_Wp2;
    const ccomplex_t IT_0013 = IT_0000*IT_0012;
    const ccomplex_t IT_0014 = V_u2*U_Wm2;
    const ccomplex_t IT_0015 = IT_0003*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*(IT_0013 + IT_0015);
    const ccomplex_t IT_0017 = 1.4142135623731*e_em*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0019;
    const ccomplex_t IT_0021 = pow(m_H0, 2);
    const ccomplex_t IT_0022 = pow(m_C_1, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0021 + IT_0022 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = pow(m_C_2, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_23 + IT_0021 + IT_0026 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0028 = U_d2*V_Wp1;
    const ccomplex_t IT_0029 = IT_0000*IT_0028;
    const ccomplex_t IT_0030 = V_u1*U_Wm2;
    const ccomplex_t IT_0031 = IT_0003*IT_0030;
    const ccomplex_t IT_0032 = IT_0007*(IT_0029 + IT_0031);
    const ccomplex_t IT_0033 = 1.4142135623731*e_em*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0037 = IT_0000*IT_0036;
    const ccomplex_t IT_0038 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0039 = IT_0003*IT_0038;
    const ccomplex_t IT_0040 = IT_0007*(IT_0037 + IT_0039);
    const ccomplex_t IT_0041 = 1.4142135623731*e_em*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0035*IT_0043;
    const ccomplex_t IT_0045 = IT_0027*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0021 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0044*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0021 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0020*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0025 + -IT_0046 + IT_0049 + -IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0055 = U_d1*V_Wp1;
    const ccomplex_t IT_0056 = IT_0000*IT_0055;
    const ccomplex_t IT_0057 = V_u1*U_Wm1;
    const ccomplex_t IT_0058 = IT_0003*IT_0057;
    const ccomplex_t IT_0059 = IT_0007*(IT_0056 + IT_0058);
    const ccomplex_t IT_0060 = 1.4142135623731*e_em*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = IT_0035*IT_0062;
    const ccomplex_t IT_0064 = IT_0027*IT_0054*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0066 = IT_0019*IT_0035;
    const ccomplex_t IT_0067 = IT_0023*IT_0065*IT_0066;
    const ccomplex_t IT_0068 = m_C_1*IT_0046;
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0022 + -IT_0026 + -reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0000*IT_0030;
    const ccomplex_t IT_0071 = IT_0003*IT_0028;
    const ccomplex_t IT_0072 = IT_0007*(IT_0070 + -IT_0071);
    const ccomplex_t IT_0073 = 1.4142135623731*e_em*IT_0072;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = cos(theta_W);
    const ccomplex_t IT_0077 = cpow(IT_0076, -2);
    const ccomplex_t IT_0078 = IT_0006*IT_0077;
    const ccomplex_t IT_0079 = IT_0007 + IT_0078;
    const ccomplex_t IT_0080 = cos(beta);
    const ccomplex_t IT_0081 = cpow(IT_0003, 3);
    const ccomplex_t IT_0082 = sin(beta);
    const ccomplex_t IT_0083 = cpow(IT_0003, 2);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_W*e_em*IT_0079*(IT_0080
      *IT_0081 + 5*IT_0000*((-1)*IT_0000*(IT_0003*IT_0080 + 1./5*IT_0000*IT_0082
      ) + IT_0082*IT_0083));
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0075*IT_0085;
    const ccomplex_t IT_0087 = IT_0069*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + -IT_0022 + -IT_0026 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0090 = cpow(IT_0000, 2);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_W*e_em*IT_0079*(IT_0081
      *IT_0082 + -IT_0000*(IT_0000*IT_0003*IT_0082 + IT_0080*(IT_0083 + -IT_0090
      )));
    const ccomplex_t IT_0092 = (-1.5)*IT_0091;
    const ccomplex_t IT_0093 = IT_0035*IT_0092;
    const ccomplex_t IT_0094 = IT_0089*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0047*IT_0054*IT_0063;
    const ccomplex_t IT_0097 = IT_0050*IT_0065*IT_0066;
    const ccomplex_t IT_0098 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0099 = IT_0000*IT_0098;
    const ccomplex_t IT_0100 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0101 = IT_0003*IT_0100;
    const ccomplex_t IT_0102 = IT_0007*(IT_0099 + IT_0101);
    const ccomplex_t IT_0103 = 1.4142135623731*e_em*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = IT_0035*IT_0105;
    const ccomplex_t IT_0107 = IT_0023*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = m_C_2*IT_0108;
    const ccomplex_t IT_0110 = m_C_1*IT_0052;
    const ccomplex_t IT_0111 = IT_0011*IT_0062;
    const ccomplex_t IT_0112 = IT_0047*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = m_C_2*IT_0113;
    const ccomplex_t IT_0115 = IT_0064 + IT_0067 + IT_0068 + IT_0088 + IT_0095
       + IT_0096 + IT_0097 + IT_0109 + IT_0110 + IT_0114;
    const ccomplex_t IT_0116 = s_24*m_C_1;
    const ccomplex_t IT_0117 = 2*IT_0116;
    const ccomplex_t IT_0118 = IT_0000*IT_0004;
    const ccomplex_t IT_0119 = IT_0001*IT_0003;
    const ccomplex_t IT_0120 = IT_0007*(IT_0118 + -IT_0119);
    const ccomplex_t IT_0121 = 1.4142135623731*e_em*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = IT_0085*IT_0123;
    const ccomplex_t IT_0125 = IT_0069*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = IT_0011*IT_0043;
    const ccomplex_t IT_0128 = IT_0047*IT_0054*IT_0127;
    const ccomplex_t IT_0129 = IT_0027*IT_0054*IT_0127;
    const ccomplex_t IT_0130 = IT_0011*IT_0105;
    const ccomplex_t IT_0131 = IT_0023*IT_0065*IT_0130;
    const ccomplex_t IT_0132 = m_C_2*IT_0049;
    const ccomplex_t IT_0133 = IT_0027*IT_0111;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = m_C_1*IT_0134;
    const ccomplex_t IT_0136 = IT_0050*IT_0106;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = m_C_1*IT_0137;
    const ccomplex_t IT_0139 = m_C_2*IT_0025;
    const ccomplex_t IT_0140 = IT_0050*IT_0065*IT_0130;
    const ccomplex_t IT_0141 = IT_0011*IT_0092;
    const ccomplex_t IT_0142 = IT_0089*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0126 + IT_0128 + IT_0129 + IT_0131 + IT_0132
       + IT_0135 + IT_0138 + IT_0139 + IT_0140 + IT_0143;
    const ccomplex_t IT_0145 = s_14*m_C_2;
    const ccomplex_t IT_0146 = (-2)*IT_0145;
    const ccomplex_t IT_0147 = s_12*IT_0021;
    const ccomplex_t IT_0148 = s_14*s_24;
    const ccomplex_t IT_0149 = (-2)*IT_0148;
    const ccomplex_t IT_0150 = IT_0147 + IT_0149;
    const ccomplex_t IT_0151 = (-2)*IT_0150;
    const ccomplex_t IT_0152 = IT_0108 + IT_0113 + -IT_0134 + -IT_0137;
    const ccomplex_t IT_0153 = m_C_1*m_C_2*IT_0021;
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

