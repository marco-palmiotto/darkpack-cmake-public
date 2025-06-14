#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_h_h(
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
    const ccomplex_t IT_0001 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + -IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = V_u1*U_Wm2;
    const ccomplex_t IT_0013 = IT_0000*IT_0012;
    const ccomplex_t IT_0014 = U_d2*V_Wp1;
    const ccomplex_t IT_0015 = IT_0003*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*(IT_0013 + -IT_0015);
    const ccomplex_t IT_0017 = 1.4142135623731*e_em*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0019;
    const ccomplex_t IT_0021 = pow(m_h, 2);
    const ccomplex_t IT_0022 = pow(m_C_2, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0021 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0027 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0028 = IT_0000*IT_0027;
    const ccomplex_t IT_0029 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = IT_0007*(IT_0028 + -IT_0030);
    const ccomplex_t IT_0032 = 1.4142135623731*e_em*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = V_u2*U_Wm2;
    const ccomplex_t IT_0036 = IT_0000*IT_0035;
    const ccomplex_t IT_0037 = U_d2*V_Wp2;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = IT_0007*(IT_0036 + -IT_0038);
    const ccomplex_t IT_0040 = 1.4142135623731*e_em*IT_0039;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = IT_0026*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_13 + IT_0021 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0043*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0021 + IT_0022 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0020*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0025 + -IT_0045 + IT_0048 + IT_0051;
    const ccomplex_t IT_0053 = cpow(IT_0042, 2);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_C_2*IT_0053;
    const ccomplex_t IT_0055 = IT_0026*IT_0054;
    const ccomplex_t IT_0056 = IT_0046*IT_0054;
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + (-2)*IT_0022 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0000*IT_0037;
    const ccomplex_t IT_0059 = IT_0003*IT_0035;
    const ccomplex_t IT_0060 = IT_0007*(IT_0058 + IT_0059);
    const ccomplex_t IT_0061 = 1.4142135623731*e_em*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = cos(theta_W);
    const ccomplex_t IT_0065 = cpow(IT_0064, -2);
    const ccomplex_t IT_0066 = IT_0006*IT_0065;
    const ccomplex_t IT_0067 = IT_0007 + IT_0066;
    const ccomplex_t IT_0068 = sin(beta);
    const ccomplex_t IT_0069 = cpow(IT_0003, 3);
    const ccomplex_t IT_0070 = cos(beta);
    const ccomplex_t IT_0071 = cpow(IT_0000, 2);
    const ccomplex_t IT_0072 = cpow(IT_0003, 2);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_W*e_em*IT_0067*(IT_0068
      *IT_0069 + (-5)*IT_0000*(IT_0000*IT_0003*IT_0068 + (-0.2)*IT_0070*(IT_0071
       + (-5)*IT_0072)));
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0063*IT_0074;
    const ccomplex_t IT_0076 = IT_0057*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0079 = V_u2*U_Wm1;
    const ccomplex_t IT_0080 = IT_0000*IT_0079;
    const ccomplex_t IT_0081 = U_d1*V_Wp2;
    const ccomplex_t IT_0082 = IT_0003*IT_0081;
    const ccomplex_t IT_0083 = IT_0007*(IT_0080 + -IT_0082);
    const ccomplex_t IT_0084 = 1.4142135623731*e_em*IT_0083;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0019*IT_0086;
    const ccomplex_t IT_0088 = IT_0023*IT_0078*IT_0087;
    const ccomplex_t IT_0089 = IT_0049*IT_0078*IT_0087;
    const ccomplex_t IT_0090 = m_C_2*IT_0025;
    const ccomplex_t IT_0091 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0092 = IT_0000*IT_0091;
    const ccomplex_t IT_0093 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0094 = IT_0003*IT_0093;
    const ccomplex_t IT_0095 = IT_0007*(IT_0092 + -IT_0094);
    const ccomplex_t IT_0096 = 1.4142135623731*e_em*IT_0095;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0086*IT_0098;
    const ccomplex_t IT_0100 = IT_0049*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = m_C_2*IT_0101;
    const ccomplex_t IT_0103 = m_C_2*IT_0045;
    const ccomplex_t IT_0104 = m_C_2*IT_0048;
    const ccomplex_t IT_0105 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_W*e_em*IT_0067*(IT_0069
      *IT_0070 + -IT_0000*(IT_0000*(IT_0000*IT_0068 + IT_0003*IT_0070) + 
      -IT_0068*IT_0072));
    const ccomplex_t IT_0107 = 1.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0042*IT_0107;
    const ccomplex_t IT_0109 = IT_0105*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0055 + IT_0056 + IT_0077 + IT_0088 + IT_0089
       + IT_0090 + IT_0102 + IT_0103 + IT_0104 + IT_0110;
    const ccomplex_t IT_0112 = s_24*m_C_2;
    const ccomplex_t IT_0113 = 2*IT_0112;
    const ccomplex_t IT_0114 = cpow(IT_0034, 2);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*m_C_2*IT_0114;
    const ccomplex_t IT_0116 = IT_0026*IT_0115;
    const ccomplex_t IT_0117 = IT_0046*IT_0115;
    const ccomplex_t IT_0118 = IT_0011*IT_0098;
    const ccomplex_t IT_0119 = IT_0023*IT_0078*IT_0118;
    const ccomplex_t IT_0120 = IT_0049*IT_0078*IT_0118;
    const ccomplex_t IT_0121 = m_C_2*IT_0051;
    const ccomplex_t IT_0122 = IT_0023*IT_0099;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = m_C_2*IT_0123;
    const ccomplex_t IT_0125 = IT_0034*IT_0107;
    const ccomplex_t IT_0126 = IT_0105*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = IT_0000*IT_0029;
    const ccomplex_t IT_0129 = IT_0003*IT_0027;
    const ccomplex_t IT_0130 = IT_0007*(IT_0128 + IT_0129);
    const ccomplex_t IT_0131 = 1.4142135623731*e_em*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0074*IT_0133;
    const ccomplex_t IT_0135 = IT_0057*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = IT_0103 + IT_0104 + IT_0116 + IT_0117 + IT_0119
       + IT_0120 + IT_0121 + IT_0124 + IT_0127 + IT_0136;
    const ccomplex_t IT_0138 = s_14*m_C_2;
    const ccomplex_t IT_0139 = (-2)*IT_0138;
    const ccomplex_t IT_0140 = s_12*IT_0021;
    const ccomplex_t IT_0141 = s_14*s_24;
    const ccomplex_t IT_0142 = (-2)*IT_0141;
    const ccomplex_t IT_0143 = IT_0140 + IT_0142;
    const ccomplex_t IT_0144 = (-2)*IT_0143;
    const ccomplex_t IT_0145 = -IT_0045 + IT_0048 + IT_0101 + -IT_0123;
    const ccomplex_t IT_0146 = IT_0021*IT_0022;
    const ccomplex_t IT_0147 = (-2)*IT_0146;
    const ccomplex_t IT_0148 = (-2)*IT_0022;
    const ccomplex_t IT_0149 = 0.5*conj(IT_0145);
    const ccomplex_t IT_0150 = IT_0145*(IT_0113*conj(IT_0137) + conj(IT_0111)
      *IT_0139 + IT_0144*conj(IT_0145) + conj(IT_0052)*IT_0147) + IT_0052*(conj
      (IT_0111)*IT_0113 + conj(IT_0137)*IT_0139 + conj(IT_0052)*IT_0144 + conj
      (IT_0145)*IT_0147) + 2*IT_0137*(s_12*conj(IT_0137) + 0.5*conj(IT_0052)
      *IT_0139 + 0.5*conj(IT_0111)*IT_0148 + IT_0113*IT_0149) + 2*IT_0111*(s_12
      *conj(IT_0111) + 0.5*conj(IT_0052)*IT_0113 + 0.5*conj(IT_0137)*IT_0148 +
       IT_0139*IT_0149);
    return create_ccomplex_return(IT_0150);
}

