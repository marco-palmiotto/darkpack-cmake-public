#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_L_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_L_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_us = param->V_us;
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
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_u, 2);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -2);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0005
      *IT_0008*IT_0010*IT_0011 + 0.5*m_W*(IT_0001*IT_0003 + -IT_0009*IT_0011)*
      (IT_0005 + -1./3*IT_0004*IT_0013));
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = pow(m_su_L, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0005*IT_0011;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0003
      *IT_0005*IT_0008*IT_0010 + (-0.5)*m_W*(IT_0003*IT_0009 + IT_0001*IT_0011)*
      (IT_0005 + -1./3*IT_0004*IT_0013));
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = pow(m_s, 2);
    const ccomplex_t IT_0031 = cpow(IT_0012, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*-3)*e_em*(IT_0005*IT_0012 + -1.
      /3*IT_0004*IT_0031);
    const ccomplex_t IT_0033 = 0.166666666666667*IT_0032;
    const ccomplex_t IT_0034 = m_s*IT_0033;
    const ccomplex_t IT_0035 = 0.333333333333333*IT_0032;
    const ccomplex_t IT_0036 = m_s*IT_0035;
    const ccomplex_t IT_0037 = IT_0004*IT_0031;
    const ccomplex_t IT_0038 = e_em*IT_0037;
    const ccomplex_t IT_0039 = IT_0005*IT_0012;
    const ccomplex_t IT_0040 = e_em*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0038 + 3*IT_0040);
    const ccomplex_t IT_0042 = (-0.166666666666667)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0031;
    const ccomplex_t IT_0047 = 0.333333333333333*IT_0046;
    const ccomplex_t IT_0048 = IT_0043*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = (IT_0034 + -IT_0036)*(IT_0045 + -IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = -(IT_0034 + -IT_0036)*(IT_0045 + -IT_0049);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = (-18)*IT_0030;
    const ccomplex_t IT_0055 = 18*s_34;
    const ccomplex_t IT_0056 = IT_0054 + IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*V_us*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*V_us*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = IT_0058*IT_0060;
    const ccomplex_t IT_0062 = IT_0057*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*V_us*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*V_us*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = IT_0064*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0063 + IT_0070;
    const ccomplex_t IT_0072 = IT_0035*IT_0045;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_us*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*V_us
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0075*IT_0077;
    const ccomplex_t IT_0079 = IT_0064*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*V_us
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_us*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = IT_0057*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0080 + IT_0087;
    const ccomplex_t IT_0089 = IT_0035*IT_0049;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = s_34*IT_0017;
    const ccomplex_t IT_0092 = (-6)*IT_0091;
    const ccomplex_t IT_0093 = s_13*s_14;
    const ccomplex_t IT_0094 = 12*IT_0093;
    const ccomplex_t IT_0095 = IT_0092 + IT_0094;
    const ccomplex_t IT_0096 = (-18)*IT_0091;
    const ccomplex_t IT_0097 = 36*IT_0093;
    const ccomplex_t IT_0098 = IT_0096 + IT_0097;
    const ccomplex_t IT_0099 = cpow(s_12 + IT_0017 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0101 = (-0.333333333333333)*IT_0100;
    const ccomplex_t IT_0102 = (-1.33333333333333)*IT_0100;
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = IT_0099*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0107 = (-2)*IT_0106;
    const ccomplex_t IT_0108 = IT_0099*IT_0106;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = 0.5*IT_0105 + (-0.0833333333333333)*IT_0110;
    const ccomplex_t IT_0112 = cpow(IT_0111, 2);
    const ccomplex_t IT_0113 = (-36)*IT_0091;
    const ccomplex_t IT_0114 = IT_0017*IT_0030;
    const ccomplex_t IT_0115 = (-36)*IT_0114;
    const ccomplex_t IT_0116 = 72*IT_0093;
    const ccomplex_t IT_0117 = IT_0113 + IT_0115 + IT_0116;
    const ccomplex_t IT_0118 = 0.25*IT_0110;
    const ccomplex_t IT_0119 = IT_0065*IT_0075;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0121 = IT_0064*IT_0119*IT_0120;
    const ccomplex_t IT_0122 = m_s*IT_0087;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0124 = IT_0058*IT_0084;
    const ccomplex_t IT_0125 = IT_0057*IT_0123*IT_0124;
    const ccomplex_t IT_0126 = m_s*IT_0080;
    const ccomplex_t IT_0127 = IT_0121 + -IT_0122 + IT_0125 + -IT_0126;
    const ccomplex_t IT_0128 = IT_0060*IT_0082;
    const ccomplex_t IT_0129 = IT_0057*IT_0123*IT_0128;
    const ccomplex_t IT_0130 = IT_0067*IT_0077;
    const ccomplex_t IT_0131 = IT_0064*IT_0120*IT_0130;
    const ccomplex_t IT_0132 = m_s*IT_0063;
    const ccomplex_t IT_0133 = m_s*IT_0070;
    const ccomplex_t IT_0134 = IT_0129 + IT_0131 + -IT_0132 + -IT_0133;
    const ccomplex_t IT_0135 = (-18)*IT_0134;
    const ccomplex_t IT_0136 = (-6)*IT_0051;
    const ccomplex_t IT_0137 = IT_0135 + IT_0136;
    const ccomplex_t IT_0138 = (-18)*conj(IT_0134);
    const ccomplex_t IT_0139 = (-6)*conj(IT_0051);
    const ccomplex_t IT_0140 = IT_0138 + IT_0139;
    const ccomplex_t IT_0141 = (-18)*IT_0051;
    const ccomplex_t IT_0142 = (-6)*IT_0134;
    const ccomplex_t IT_0143 = IT_0141 + IT_0142;
    const ccomplex_t IT_0144 = (-18)*conj(IT_0051);
    const ccomplex_t IT_0145 = (-6)*conj(IT_0134);
    const ccomplex_t IT_0146 = IT_0144 + IT_0145;
    const ccomplex_t IT_0147 = (-18)*IT_0071 + (-6)*IT_0073;
    const ccomplex_t IT_0148 = (-18)*IT_0073;
    const ccomplex_t IT_0149 = (-6)*IT_0071;
    const ccomplex_t IT_0150 = IT_0148 + IT_0149;
    const ccomplex_t IT_0151 = (-18)*conj(IT_0073);
    const ccomplex_t IT_0152 = (-6)*conj(IT_0071);
    const ccomplex_t IT_0153 = IT_0151 + IT_0152;
    const ccomplex_t IT_0154 = 6*s_34;
    const ccomplex_t IT_0155 = -IT_0030;
    const ccomplex_t IT_0156 = s_34 + IT_0155;
    const ccomplex_t IT_0157 = 6*IT_0029;
    const ccomplex_t IT_0158 = 6*conj(IT_0029);
    const ccomplex_t IT_0159 = m_s*s_14;
    const ccomplex_t IT_0160 = (-18)*IT_0127;
    const ccomplex_t IT_0161 = (-6)*IT_0127;
    const ccomplex_t IT_0162 = (-6)*IT_0114;
    const ccomplex_t IT_0163 = IT_0092 + IT_0094 + IT_0162;
    const ccomplex_t IT_0164 = conj(IT_0071) + conj(IT_0088);
    const ccomplex_t IT_0165 = conj(IT_0073) + conj(IT_0090);
    const ccomplex_t IT_0166 = IT_0073 + IT_0090;
    const ccomplex_t IT_0167 = -IT_0118;
    const ccomplex_t IT_0168 = IT_0071 + IT_0088;
    const ccomplex_t IT_0169 = -IT_0111;
    const ccomplex_t IT_0170 = (-18)*IT_0114;
    const ccomplex_t IT_0171 = IT_0096 + IT_0097 + IT_0170;
    const ccomplex_t IT_0172 = m_s*s_13;
    const ccomplex_t IT_0173 = 6*conj(IT_0051);
    const ccomplex_t IT_0174 = 18*conj(IT_0134);
    const ccomplex_t IT_0175 = 6*conj(IT_0134);
    const ccomplex_t IT_0176 = 18*conj(IT_0051);
    const ccomplex_t IT_0177 = 18*conj(IT_0127);
    const ccomplex_t IT_0178 = 6*conj(IT_0127);
    const ccomplex_t IT_0179 = 18*conj(IT_0053);
    const ccomplex_t IT_0180 = IT_0051*(conj(IT_0071) + 3*conj(IT_0073)) +
       IT_0053*(conj(IT_0088) + 3*conj(IT_0090)) + 0.166666666666667*(18*conj
      (IT_0088) + 6*conj(IT_0090))*IT_0127 + 0.166666666666667*(18*conj(IT_0071)
       + 6*conj(IT_0073))*IT_0134 + 0.166666666666667*IT_0071*(IT_0173 + IT_0174
      ) + 0.166666666666667*IT_0073*(IT_0175 + IT_0176) + IT_0088*(conj(IT_0053)
       + 0.166666666666667*IT_0177) + 0.166666666666667*IT_0090*(IT_0178 +
       IT_0179);
    const ccomplex_t IT_0181 = -s_14;
    const ccomplex_t IT_0182 = s_13 + IT_0181;
    const ccomplex_t IT_0183 = m_s*IT_0182;
    const ccomplex_t IT_0184 = 18*IT_0029;
    const ccomplex_t IT_0185 = 18*conj(IT_0029);
    const ccomplex_t IT_0186 = IT_0157*IT_0164 + IT_0158*IT_0168 + IT_0118*((
      -12)*IT_0029 + 12*conj(IT_0029) + (-6)*IT_0053 + 6*conj(IT_0053) + IT_0135
       + IT_0136 + IT_0160 + IT_0173 + IT_0174 + IT_0177) + IT_0111*((-36)
      *IT_0029 + 36*conj(IT_0029) + (-18)*IT_0053 + IT_0141 + IT_0142 + IT_0161 
      + IT_0175 + IT_0176 + IT_0178 + IT_0179) + IT_0165*IT_0184 + IT_0166
      *IT_0185;
    const ccomplex_t IT_0187 = IT_0029*conj(IT_0029)*(36*s_34 + (-36)*IT_0030)
       + (conj(IT_0029)*(IT_0051 + IT_0053) + IT_0029*(conj(IT_0051) + conj
      (IT_0053)))*IT_0056 + (conj(IT_0071)*IT_0073 + IT_0071*conj(IT_0073) +
       conj(IT_0088)*IT_0090 + IT_0088*conj(IT_0090))*IT_0095 + (IT_0071*conj
      (IT_0071) + IT_0073*conj(IT_0073) + IT_0088*conj(IT_0088) + IT_0090*conj
      (IT_0090))*IT_0098 + -IT_0112*IT_0117 + (-2)*IT_0118*(IT_0111*((-12)
      *IT_0091 + 24*IT_0093 + (-12)*IT_0114) + 0.5*IT_0117*IT_0118) + IT_0055*
      (IT_0051*conj(IT_0051) + IT_0053*conj(IT_0053) + IT_0127*conj(IT_0127) +
       IT_0134*conj(IT_0134)) + IT_0030*(conj(IT_0127)*IT_0137 + IT_0127*IT_0140
       + conj(IT_0053)*IT_0143 + IT_0053*IT_0146) + IT_0114*(((-18)*conj(IT_0071
      ) + (-6)*conj(IT_0073))*IT_0088 + conj(IT_0088)*IT_0147 + conj(IT_0090)
      *IT_0150 + IT_0090*IT_0153) + (conj(IT_0053)*IT_0127 + IT_0053*conj
      (IT_0127) + conj(IT_0051)*IT_0134 + IT_0051*conj(IT_0134))*IT_0154 +
       IT_0156*((conj(IT_0127) + conj(IT_0134))*IT_0157 + (IT_0127 + IT_0134)
      *IT_0158) + IT_0159*(conj(IT_0088)*IT_0137 + IT_0088*IT_0140 + conj
      (IT_0090)*IT_0143 + IT_0090*IT_0146 + conj(IT_0127)*IT_0147 + conj(IT_0053
      )*IT_0150 + IT_0053*IT_0153 + conj(IT_0071)*IT_0160 + conj(IT_0073)
      *IT_0161) + IT_0163*(IT_0111*IT_0164 + IT_0118*IT_0165 + IT_0166*IT_0167 +
       IT_0168*IT_0169) + (IT_0118*IT_0164 + IT_0111*IT_0165 + IT_0167*IT_0168 +
       IT_0166*IT_0169)*IT_0171 + 6*IT_0172*IT_0180 + IT_0183*IT_0186;
    return create_ccomplex_return(IT_0187);
}

