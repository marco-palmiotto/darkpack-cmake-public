#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_su_R_to_anti_d_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_su_R_to_anti_d_Z(
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
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = e_em*IT_0009;
    const ccomplex_t IT_0011 = IT_0004*IT_0007;
    const ccomplex_t IT_0012 = e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0010 + 3*IT_0012);
    const ccomplex_t IT_0014 = (-0.166666666666667)*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = pow(m_d, 2);
    const ccomplex_t IT_0017 = pow(m_C_1, 2);
    const ccomplex_t IT_0018 = pow(m_su_R, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + IT_0016 + -IT_0017 + -IT_0018 
      + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_Z, -2);
    const ccomplex_t IT_0023 = pow(m_Z, 4);
    const ccomplex_t IT_0024 = s_13*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = s_14*s_34;
    const ccomplex_t IT_0027 = 24*IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, 2);
    const ccomplex_t IT_0029 = s_13*IT_0028;
    const ccomplex_t IT_0030 = (-12)*IT_0029;
    const ccomplex_t IT_0031 = 6*IT_0025 + IT_0027 + IT_0030;
    const ccomplex_t IT_0032 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0033 = IT_0009*IT_0032;
    const ccomplex_t IT_0034 = IT_0011*IT_0032;
    const ccomplex_t IT_0035 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0036 = IT_0011*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + -IT_0034 + (-2)
      *IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0006*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0016 + IT_0018 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0016 + IT_0018 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0044 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0045 = IT_0009*IT_0044;
    const ccomplex_t IT_0046 = IT_0011*IT_0044;
    const ccomplex_t IT_0047 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0048 = IT_0011*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + -IT_0046 + (-2)
      *IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = IT_0043*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0042 + IT_0055;
    const ccomplex_t IT_0057 = (-6)*IT_0025;
    const ccomplex_t IT_0058 = IT_0027 + IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0060 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0061 = IT_0009*IT_0060;
    const ccomplex_t IT_0062 = IT_0011*IT_0060;
    const ccomplex_t IT_0063 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0064 = IT_0011*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + -IT_0062 + (-2)
      *IT_0064);
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0052*IT_0066;
    const ccomplex_t IT_0068 = IT_0043*IT_0059*IT_0067;
    const ccomplex_t IT_0069 = m_C_1*IT_0042;
    const ccomplex_t IT_0070 = m_C_1*IT_0055;
    const ccomplex_t IT_0071 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0072 = IT_0009*IT_0071;
    const ccomplex_t IT_0073 = IT_0011*IT_0071;
    const ccomplex_t IT_0074 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0075 = IT_0011*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + -IT_0073 + (-2)
      *IT_0075);
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0006*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0080 = IT_0040*IT_0078*IT_0079;
    const ccomplex_t IT_0081 = IT_0068 + IT_0069 + IT_0070 + IT_0080;
    const ccomplex_t IT_0082 = s_34*IT_0028;
    const ccomplex_t IT_0083 = IT_0022*IT_0082;
    const ccomplex_t IT_0084 = (-0.25)*IT_0083;
    const ccomplex_t IT_0085 = s_34 + IT_0084;
    const ccomplex_t IT_0086 = m_C_1*IT_0085;
    const ccomplex_t IT_0087 = (-4)*IT_0086;
    const ccomplex_t IT_0088 = 6*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0008;
    const ccomplex_t IT_0090 = 0.333333333333333*IT_0089;
    const ccomplex_t IT_0091 = IT_0006*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0093 = IT_0019*IT_0091*IT_0092;
    const ccomplex_t IT_0094 = m_d*IT_0021;
    const ccomplex_t IT_0095 = -IT_0093 + IT_0094;
    const ccomplex_t IT_0096 = s_14*IT_0028;
    const ccomplex_t IT_0097 = IT_0022*IT_0096;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = s_14 + IT_0098;
    const ccomplex_t IT_0100 = m_d*IT_0099;
    const ccomplex_t IT_0101 = 2*IT_0100;
    const ccomplex_t IT_0102 = 6*IT_0101;
    const ccomplex_t IT_0103 = 1.33333333333333*IT_0089;
    const ccomplex_t IT_0104 = cpow((-2)*s_13 + IT_0016 + IT_0017 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0006*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0103*IT_0106;
    const ccomplex_t IT_0108 = 2*IT_0021 + IT_0107;
    const ccomplex_t IT_0109 = s_14*IT_0016;
    const ccomplex_t IT_0110 = s_13*s_34;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = IT_0109 + IT_0111;
    const ccomplex_t IT_0113 = (-6)*IT_0112;
    const ccomplex_t IT_0114 = (-0.25)*IT_0097;
    const ccomplex_t IT_0115 = s_14 + IT_0114;
    const ccomplex_t IT_0116 = m_d*IT_0115;
    const ccomplex_t IT_0117 = (-4)*IT_0116;
    const ccomplex_t IT_0118 = 0.5*IT_0083;
    const ccomplex_t IT_0119 = s_34 + IT_0118;
    const ccomplex_t IT_0120 = m_C_1*IT_0119;
    const ccomplex_t IT_0121 = 2*IT_0120;
    const ccomplex_t IT_0122 = (-6)*IT_0021;
    const ccomplex_t IT_0123 = s_13*s_14;
    const ccomplex_t IT_0124 = s_34*IT_0017;
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = IT_0123 + IT_0125;
    const ccomplex_t IT_0127 = (-2)*IT_0042 + (-2)*IT_0055;
    const ccomplex_t IT_0128 = -IT_0107;
    const ccomplex_t IT_0129 = (-6)*conj(IT_0021);
    const ccomplex_t IT_0130 = pow(s_34, 2);
    const ccomplex_t IT_0131 = IT_0022*IT_0130;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0016 + IT_0132;
    const ccomplex_t IT_0134 = s_13*IT_0133;
    const ccomplex_t IT_0135 = s_14*s_34*IT_0022;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = s_13 + IT_0136;
    const ccomplex_t IT_0138 = s_13*IT_0137;
    const ccomplex_t IT_0139 = 6*IT_0138;
    const ccomplex_t IT_0140 = m_C_1*IT_0133;
    const ccomplex_t IT_0141 = 6*IT_0140;
    const ccomplex_t IT_0142 = m_d*IT_0137;
    const ccomplex_t IT_0143 = 6*IT_0142;
    const ccomplex_t IT_0144 = (-6)*IT_0138;
    const ccomplex_t IT_0145 = (-6)*IT_0108*((-0.166666666666667)*conj(IT_0056
      )*IT_0113 + (-0.166666666666667)*IT_0112*IT_0129 + conj(IT_0108)*IT_0134 +
       (-0.166666666666667)*conj(IT_0127)*IT_0139 + (-0.166666666666667)*conj
      (IT_0081)*IT_0141 + (-0.166666666666667)*conj(IT_0095)*IT_0143 + (
      -0.166666666666667)*conj(IT_0128)*IT_0144);
    const ccomplex_t IT_0146 = 12*s_13;
    const ccomplex_t IT_0147 = (-2)*IT_0026;
    const ccomplex_t IT_0148 = IT_0029 + IT_0147;
    const ccomplex_t IT_0149 = IT_0022*IT_0148;
    const ccomplex_t IT_0150 = (-6)*IT_0149;
    const ccomplex_t IT_0151 = IT_0146 + IT_0150;
    const ccomplex_t IT_0152 = m_C_1*IT_0137;
    const ccomplex_t IT_0153 = 6*IT_0152;
    const ccomplex_t IT_0154 = m_d*m_C_1;
    const ccomplex_t IT_0155 = (-18)*IT_0154;
    const ccomplex_t IT_0156 = (-6)*IT_0152;
    const ccomplex_t IT_0157 = pow(s_14, 2);
    const ccomplex_t IT_0158 = IT_0022*IT_0157;
    const ccomplex_t IT_0159 = -IT_0017;
    const ccomplex_t IT_0160 = IT_0158 + IT_0159;
    const ccomplex_t IT_0161 = s_13*IT_0160;
    const ccomplex_t IT_0162 = 6*IT_0161;
    const ccomplex_t IT_0163 = m_d*IT_0160;
    const ccomplex_t IT_0164 = 6*IT_0163;
    const ccomplex_t IT_0165 = (-6)*IT_0161;
    const ccomplex_t IT_0166 = (-6)*IT_0163;
    const ccomplex_t IT_0167 = IT_0021*(conj(IT_0021)*IT_0031 + conj(IT_0056)
      *IT_0058) + IT_0056*(IT_0031*conj(IT_0056) + conj(IT_0021)*IT_0058 + conj
      (IT_0081)*IT_0088 + conj(IT_0095)*IT_0102 + conj(IT_0108)*IT_0113) + (conj
      (IT_0108)*IT_0112 + conj(IT_0095)*IT_0117 + conj(IT_0081)*IT_0121)*IT_0122
       + 6*IT_0126*((conj(IT_0021) + conj(IT_0056))*IT_0127 + (IT_0021 + IT_0056
      )*conj(IT_0127) + -(IT_0056 + (-0.166666666666667)*IT_0122)*conj(IT_0128) 
      + -IT_0128*(conj(IT_0056) + (-0.166666666666667)*IT_0129)) + IT_0145 +
       IT_0081*(conj(IT_0056)*IT_0088 + IT_0121*IT_0129 + conj(IT_0108)*IT_0141 
      + conj(IT_0081)*IT_0151 + conj(IT_0128)*IT_0153 + conj(IT_0095)*IT_0155 +
       conj(IT_0127)*IT_0156) + IT_0127*(conj(IT_0108)*IT_0139 + conj(IT_0081)
      *IT_0156 + conj(IT_0127)*IT_0162 + conj(IT_0095)*IT_0164 + conj(IT_0128)
      *IT_0165) + IT_0128*(conj(IT_0108)*IT_0144 + conj(IT_0081)*IT_0153 + conj
      (IT_0128)*IT_0162 + conj(IT_0127)*IT_0165 + conj(IT_0095)*IT_0166) +
       IT_0095*(conj(IT_0056)*IT_0102 + IT_0117*IT_0129 + conj(IT_0108)*IT_0143 
      + conj(IT_0095)*IT_0151 + conj(IT_0081)*IT_0155 + conj(IT_0127)*IT_0164 +
       conj(IT_0128)*IT_0166);
    return create_ccomplex_return(IT_0167);
}

