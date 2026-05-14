#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_A0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_A0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
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
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*(IT_0012 + IT_0014);
    const ccomplex_t IT_0016 = 1.4142135623731*e_em*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = pow(m_A0, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0019 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0022 = IT_0018*IT_0020*IT_0021;
    const ccomplex_t IT_0023 = pow(m_C_1, 2);
    const ccomplex_t IT_0024 = pow(m_C_2, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0023 + -IT_0024 + -reg_prop, -1);
    const ccomplex_t IT_0026 = cos(alpha);
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = sin(alpha);
    const ccomplex_t IT_0029 = IT_0004*IT_0028;
    const ccomplex_t IT_0030 = IT_0007*(IT_0027 + -IT_0029);
    const ccomplex_t IT_0031 = 1.4142135623731*e_em*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = cos(theta_W);
    const ccomplex_t IT_0035 = cpow(IT_0034, -2);
    const ccomplex_t IT_0036 = IT_0006*IT_0035;
    const ccomplex_t IT_0037 = IT_0007 + IT_0036;
    const ccomplex_t IT_0038 = cpow(IT_0003, 3);
    const ccomplex_t IT_0039 = cpow(IT_0003, 2);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_W*e_em*IT_0037*(IT_0026
      *IT_0038 + -IT_0000*(IT_0000*(IT_0003*IT_0026 + IT_0000*IT_0028) + 
      -IT_0028*IT_0039));
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0033*IT_0041;
    const ccomplex_t IT_0043 = IT_0025*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + -IT_0023 + -IT_0024 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0004*IT_0026;
    const ccomplex_t IT_0047 = IT_0001*IT_0028;
    const ccomplex_t IT_0048 = IT_0007*(IT_0046 + IT_0047);
    const ccomplex_t IT_0049 = 1.4142135623731*e_em*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = -IT_0000;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_W*e_em*IT_0037*(IT_0028
      *IT_0038 + IT_0052*(IT_0026*IT_0039 + (IT_0000*IT_0026 + -IT_0003*IT_0028)
      *IT_0052));
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0051*IT_0054;
    const ccomplex_t IT_0056 = IT_0045*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0059 = IT_0000*IT_0058;
    const ccomplex_t IT_0060 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0061 = IT_0003*IT_0060;
    const ccomplex_t IT_0062 = IT_0007*(IT_0059 + IT_0061);
    const ccomplex_t IT_0063 = 1.4142135623731*e_em*IT_0062;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = IT_0010*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0019 + IT_0024 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0065*IT_0066*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_13 + IT_0019 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0065*IT_0066*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0019 + IT_0023 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0018*IT_0021*IT_0071;
    const ccomplex_t IT_0073 = V_u1*U_Wm2;
    const ccomplex_t IT_0074 = IT_0000*IT_0073;
    const ccomplex_t IT_0075 = U_d2*V_Wp1;
    const ccomplex_t IT_0076 = IT_0003*IT_0075;
    const ccomplex_t IT_0077 = IT_0007*(IT_0074 + IT_0076);
    const ccomplex_t IT_0078 = 1.4142135623731*e_em*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = IT_0064*IT_0079;
    const ccomplex_t IT_0081 = IT_0069*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = m_C_2*IT_0082;
    const ccomplex_t IT_0084 = V_u1*U_Wm1;
    const ccomplex_t IT_0085 = IT_0000*IT_0084;
    const ccomplex_t IT_0086 = U_d1*V_Wp1;
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = IT_0007*(IT_0085 + IT_0087);
    const ccomplex_t IT_0089 = 1.4142135623731*e_em*IT_0088;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = IT_0010*IT_0090;
    const ccomplex_t IT_0092 = IT_0067*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = m_C_1*IT_0093;
    const ccomplex_t IT_0095 = V_u2*U_Wm2;
    const ccomplex_t IT_0096 = IT_0000*IT_0095;
    const ccomplex_t IT_0097 = U_d2*V_Wp2;
    const ccomplex_t IT_0098 = IT_0003*IT_0097;
    const ccomplex_t IT_0099 = IT_0007*(IT_0096 + IT_0098);
    const ccomplex_t IT_0100 = 1.4142135623731*e_em*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0010*IT_0101;
    const ccomplex_t IT_0103 = IT_0071*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = m_C_2*IT_0104;
    const ccomplex_t IT_0106 = IT_0017*IT_0079;
    const ccomplex_t IT_0107 = IT_0020*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = m_C_1*IT_0108;
    const ccomplex_t IT_0110 = IT_0022 + IT_0044 + IT_0057 + IT_0068 + IT_0070
       + IT_0072 + IT_0083 + IT_0094 + IT_0105 + IT_0109;
    const ccomplex_t IT_0111 = IT_0079*IT_0090;
    const ccomplex_t IT_0112 = IT_0066*IT_0069*IT_0111;
    const ccomplex_t IT_0113 = IT_0079*IT_0101;
    const ccomplex_t IT_0114 = IT_0020*IT_0021*IT_0113;
    const ccomplex_t IT_0115 = IT_0026*IT_0075;
    const ccomplex_t IT_0116 = IT_0028*IT_0073;
    const ccomplex_t IT_0117 = IT_0007*(IT_0115 + IT_0116);
    const ccomplex_t IT_0118 = 1.4142135623731*e_em*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = IT_0054*IT_0120;
    const ccomplex_t IT_0122 = IT_0045*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = IT_0066*IT_0067*IT_0111;
    const ccomplex_t IT_0125 = IT_0021*IT_0071*IT_0113;
    const ccomplex_t IT_0126 = IT_0067*IT_0080;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = m_C_1*IT_0127;
    const ccomplex_t IT_0129 = IT_0069*IT_0091;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = m_C_2*IT_0130;
    const ccomplex_t IT_0132 = IT_0020*IT_0102;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = m_C_1*IT_0133;
    const ccomplex_t IT_0135 = IT_0071*IT_0106;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = m_C_2*IT_0136;
    const ccomplex_t IT_0138 = IT_0026*IT_0073;
    const ccomplex_t IT_0139 = IT_0028*IT_0075;
    const ccomplex_t IT_0140 = IT_0007*(IT_0138 + -IT_0139);
    const ccomplex_t IT_0141 = 1.4142135623731*e_em*IT_0140;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0041*IT_0143;
    const ccomplex_t IT_0145 = IT_0025*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0112 + IT_0114 + IT_0123 + IT_0124 + IT_0125
       + IT_0128 + IT_0131 + IT_0134 + IT_0137 + IT_0146;
    const ccomplex_t IT_0148 = m_C_1*m_C_2;
    const ccomplex_t IT_0149 = IT_0082 + IT_0104 + -IT_0127 + -IT_0133;
    const ccomplex_t IT_0150 = s_14*m_C_2;
    const ccomplex_t IT_0151 = -IT_0093 + -IT_0108 + IT_0130 + IT_0136;
    const ccomplex_t IT_0152 = s_24*m_C_1;
    const ccomplex_t IT_0153 = 2*s_12;
    const ccomplex_t IT_0154 = s_12*IT_0019;
    const ccomplex_t IT_0155 = s_14*s_24;
    const ccomplex_t IT_0156 = (-2)*IT_0155;
    const ccomplex_t IT_0157 = IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = m_C_1*m_C_2*IT_0019;
    const ccomplex_t IT_0159 = (-2)*IT_0151;
    const ccomplex_t IT_0160 = (-2)*conj(IT_0147);
    const ccomplex_t IT_0161 = (-2)*conj(IT_0110)*(IT_0147*IT_0148 + IT_0149
      *IT_0150 + -IT_0151*IT_0152 + (-0.5)*IT_0110*IT_0153) + 2*conj(IT_0147)*
      (IT_0149*IT_0152 + 0.5*IT_0147*IT_0153) + (-2)*conj(IT_0151)*(IT_0147
      *IT_0150 + -IT_0110*IT_0152 + IT_0149*IT_0158 + (-0.5)*IT_0157*IT_0159) + 
      (-2)*conj(IT_0149)*(IT_0110*IT_0150 + -IT_0147*IT_0152 + IT_0149*IT_0157 +
       (-0.5)*IT_0158*IT_0159) + (IT_0110*IT_0148 + IT_0150*IT_0151)*IT_0160;
    return create_ccomplex_return(IT_0161);
}

