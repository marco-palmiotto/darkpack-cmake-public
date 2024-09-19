#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_L_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_L_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
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
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_su_L, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = pow(m_u, 2);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cos(alpha);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0003
      *IT_0005*IT_0006*IT_0008 + (-0.5)*m_W*(IT_0004*IT_0006 + -IT_0009*IT_0010)
      *(IT_0008 + -1./3*IT_0007*IT_0012));
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0009, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*m_b*e_em*IT_0002*IT_0008
      *IT_0010*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0003
      *IT_0005*IT_0008*IT_0010 + (-0.5)*m_W*(IT_0006*IT_0009 + IT_0004*IT_0010)*
      (IT_0008 + -1./3*IT_0007*IT_0012));
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*m_b*e_em*IT_0002*IT_0006
      *IT_0008*IT_0015;
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = IT_0021*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = pow(m_b, 2);
    const ccomplex_t IT_0031 = cpow(IT_0011, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*-3)*e_em*(IT_0008*IT_0011 + -1.
      /3*IT_0007*IT_0031);
    const ccomplex_t IT_0033 = 0.166666666666667*IT_0032;
    const ccomplex_t IT_0034 = m_b*IT_0033;
    const ccomplex_t IT_0035 = 0.333333333333333*IT_0032;
    const ccomplex_t IT_0036 = m_b*IT_0035;
    const ccomplex_t IT_0037 = IT_0007*IT_0031;
    const ccomplex_t IT_0038 = e_em*IT_0037;
    const ccomplex_t IT_0039 = IT_0008*IT_0011;
    const ccomplex_t IT_0040 = e_em*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0038 + 3*IT_0040);
    const ccomplex_t IT_0042 = (-0.166666666666667)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0007*IT_0031;
    const ccomplex_t IT_0047 = 0.333333333333333*IT_0046;
    const ccomplex_t IT_0048 = IT_0043*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = (IT_0034 + -IT_0036)*(IT_0045 + -IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = -(IT_0034 + -IT_0036)*(IT_0045 + -IT_0049);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = (-18)*IT_0030;
    const ccomplex_t IT_0057 = 18*s_34;
    const ccomplex_t IT_0058 = IT_0056 + IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0000 + IT_0030 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0060 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0008
      *IT_0060*V_ub_mod;
    const ccomplex_t IT_0062 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0008*IT_0062
      *V_ub_mod;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0061*IT_0064;
    const ccomplex_t IT_0066 = IT_0059*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_23 + IT_0000 + IT_0030 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0008
      *IT_0060*V_ub_mod;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0008*IT_0062
      *V_ub_mod;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = IT_0068*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0067 + IT_0074;
    const ccomplex_t IT_0076 = IT_0035*IT_0045;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *e_em*IT_0002*IT_0008*IT_0015*IT_0062*V_ub_mod;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*e_em
      *IT_0002*IT_0008*IT_0015*IT_0060*V_ub_mod;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = IT_0068*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *e_em*IT_0002*IT_0008*IT_0015*IT_0062*V_ub_mod;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*e_em
      *IT_0002*IT_0008*IT_0015*IT_0060*V_ub_mod;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = IT_0059*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0084 + IT_0091;
    const ccomplex_t IT_0093 = IT_0035*IT_0049;
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = s_34*IT_0000;
    const ccomplex_t IT_0096 = (-6)*IT_0095;
    const ccomplex_t IT_0097 = s_13*s_14;
    const ccomplex_t IT_0098 = 12*IT_0097;
    const ccomplex_t IT_0099 = IT_0096 + IT_0098;
    const ccomplex_t IT_0100 = (-18)*IT_0095;
    const ccomplex_t IT_0101 = 36*IT_0097;
    const ccomplex_t IT_0102 = IT_0100 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0104 = (-0.333333333333333)*IT_0103;
    const ccomplex_t IT_0105 = (-1.33333333333333)*IT_0103;
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = cpow(s_12 + IT_0000 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0111 = (-2)*IT_0110;
    const ccomplex_t IT_0112 = IT_0107*IT_0110;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0109 + (-0.0833333333333333)*IT_0114;
    const ccomplex_t IT_0116 = cpow(IT_0115, 2);
    const ccomplex_t IT_0117 = (-36)*IT_0095;
    const ccomplex_t IT_0118 = IT_0000*IT_0030;
    const ccomplex_t IT_0119 = (-36)*IT_0118;
    const ccomplex_t IT_0120 = 72*IT_0097;
    const ccomplex_t IT_0121 = IT_0117 + IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = 0.25*IT_0114;
    const ccomplex_t IT_0123 = (-18)*IT_0053;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0125 = IT_0071*IT_0081;
    const ccomplex_t IT_0126 = IT_0068*IT_0124*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0128 = IT_0064*IT_0088;
    const ccomplex_t IT_0129 = IT_0059*IT_0127*IT_0128;
    const ccomplex_t IT_0130 = m_b*IT_0074;
    const ccomplex_t IT_0131 = m_b*IT_0067;
    const ccomplex_t IT_0132 = IT_0126 + IT_0129 + -IT_0130 + -IT_0131;
    const ccomplex_t IT_0133 = (-6)*IT_0132;
    const ccomplex_t IT_0134 = IT_0123 + IT_0133;
    const ccomplex_t IT_0135 = (-18)*conj(IT_0053);
    const ccomplex_t IT_0136 = (-6)*conj(IT_0132);
    const ccomplex_t IT_0137 = IT_0135 + IT_0136;
    const ccomplex_t IT_0138 = IT_0069*IT_0079;
    const ccomplex_t IT_0139 = IT_0068*IT_0124*IT_0138;
    const ccomplex_t IT_0140 = IT_0061*IT_0086;
    const ccomplex_t IT_0141 = IT_0059*IT_0127*IT_0140;
    const ccomplex_t IT_0142 = m_b*IT_0084;
    const ccomplex_t IT_0143 = m_b*IT_0091;
    const ccomplex_t IT_0144 = IT_0139 + IT_0141 + -IT_0142 + -IT_0143;
    const ccomplex_t IT_0145 = (-18)*IT_0132;
    const ccomplex_t IT_0146 = (-6)*IT_0053;
    const ccomplex_t IT_0147 = IT_0145 + IT_0146;
    const ccomplex_t IT_0148 = (-18)*conj(IT_0132);
    const ccomplex_t IT_0149 = (-6)*conj(IT_0053);
    const ccomplex_t IT_0150 = IT_0148 + IT_0149;
    const ccomplex_t IT_0151 = (-18)*IT_0075 + (-6)*IT_0077;
    const ccomplex_t IT_0152 = (-18)*IT_0077;
    const ccomplex_t IT_0153 = (-6)*IT_0075;
    const ccomplex_t IT_0154 = IT_0152 + IT_0153;
    const ccomplex_t IT_0155 = (-18)*conj(IT_0077);
    const ccomplex_t IT_0156 = (-6)*conj(IT_0075);
    const ccomplex_t IT_0157 = IT_0155 + IT_0156;
    const ccomplex_t IT_0158 = 6*s_34;
    const ccomplex_t IT_0159 = -IT_0030;
    const ccomplex_t IT_0160 = s_34 + IT_0159;
    const ccomplex_t IT_0161 = 6*IT_0029;
    const ccomplex_t IT_0162 = 6*conj(IT_0029);
    const ccomplex_t IT_0163 = m_b*s_13;
    const ccomplex_t IT_0164 = 6*conj(IT_0144);
    const ccomplex_t IT_0165 = 6*conj(IT_0132);
    const ccomplex_t IT_0166 = 18*conj(IT_0053);
    const ccomplex_t IT_0167 = 6*conj(IT_0053);
    const ccomplex_t IT_0168 = 18*conj(IT_0132);
    const ccomplex_t IT_0169 = 18*conj(IT_0144);
    const ccomplex_t IT_0170 = IT_0053*(conj(IT_0075) + 3*conj(IT_0077)) +
       IT_0055*(conj(IT_0092) + 3*conj(IT_0094)) + 3*(conj(IT_0075) +
       0.333333333333333*conj(IT_0077))*IT_0132 + 3*(conj(IT_0092) +
       0.333333333333333*conj(IT_0094))*IT_0144 + 3*IT_0094*(conj(IT_0055) +
       0.0555555555555556*IT_0164) + 0.166666666666667*IT_0077*(IT_0165 +
       IT_0166) + 0.166666666666667*IT_0075*(IT_0167 + IT_0168) + IT_0092*(conj
      (IT_0055) + 0.166666666666667*IT_0169);
    const ccomplex_t IT_0171 = -s_14;
    const ccomplex_t IT_0172 = s_13 + IT_0171;
    const ccomplex_t IT_0173 = m_b*IT_0172;
    const ccomplex_t IT_0174 = conj(IT_0075) + conj(IT_0092);
    const ccomplex_t IT_0175 = IT_0075 + IT_0092;
    const ccomplex_t IT_0176 = conj(IT_0077) + conj(IT_0094);
    const ccomplex_t IT_0177 = 18*IT_0029;
    const ccomplex_t IT_0178 = IT_0077 + IT_0094;
    const ccomplex_t IT_0179 = 18*conj(IT_0029);
    const ccomplex_t IT_0180 = (-18)*IT_0144;
    const ccomplex_t IT_0181 = (-6)*IT_0144;
    const ccomplex_t IT_0182 = IT_0161*IT_0174 + IT_0162*IT_0175 + IT_0176
      *IT_0177 + IT_0178*IT_0179 + IT_0122*((-12)*IT_0029 + 12*conj(IT_0029) + (
      -6)*IT_0055 + 6*conj(IT_0055) + IT_0145 + IT_0146 + IT_0167 + IT_0168 +
       IT_0169 + IT_0180) + IT_0115*((-36)*IT_0029 + 36*conj(IT_0029) + (-18)
      *IT_0055 + 18*conj(IT_0055) + IT_0123 + IT_0133 + IT_0164 + IT_0165 +
       IT_0166 + IT_0181);
    const ccomplex_t IT_0183 = m_b*s_14;
    const ccomplex_t IT_0184 = (-6)*IT_0118;
    const ccomplex_t IT_0185 = IT_0096 + IT_0098 + IT_0184;
    const ccomplex_t IT_0186 = -IT_0122;
    const ccomplex_t IT_0187 = -IT_0115;
    const ccomplex_t IT_0188 = (-18)*IT_0118;
    const ccomplex_t IT_0189 = IT_0100 + IT_0101 + IT_0188;
    const ccomplex_t IT_0190 = IT_0029*conj(IT_0029)*(36*s_34 + (-36)*IT_0030)
       + (conj(IT_0029)*(IT_0053 + IT_0055) + IT_0029*(conj(IT_0053) + conj
      (IT_0055)))*IT_0058 + (conj(IT_0075)*IT_0077 + IT_0075*conj(IT_0077) +
       conj(IT_0092)*IT_0094 + IT_0092*conj(IT_0094))*IT_0099 + (IT_0075*conj
      (IT_0075) + IT_0077*conj(IT_0077) + IT_0092*conj(IT_0092) + IT_0094*conj
      (IT_0094))*IT_0102 + -IT_0116*IT_0121 + (-2)*IT_0122*(IT_0115*((-12)
      *IT_0095 + 24*IT_0097 + (-12)*IT_0118) + 0.5*IT_0121*IT_0122) + IT_0057*
      (IT_0053*conj(IT_0053) + IT_0055*conj(IT_0055) + IT_0132*conj(IT_0132) +
       IT_0144*conj(IT_0144)) + IT_0030*(conj(IT_0055)*IT_0134 + IT_0055*IT_0137
       + conj(IT_0144)*IT_0147 + IT_0144*IT_0150) + IT_0118*(((-18)*conj(IT_0075
      ) + (-6)*conj(IT_0077))*IT_0092 + conj(IT_0092)*IT_0151 + conj(IT_0094)
      *IT_0154 + IT_0094*IT_0157) + (conj(IT_0053)*IT_0132 + IT_0053*conj
      (IT_0132) + conj(IT_0055)*IT_0144 + IT_0055*conj(IT_0144))*IT_0158 +
       IT_0160*((conj(IT_0132) + conj(IT_0144))*IT_0161 + (IT_0132 + IT_0144)
      *IT_0162) + 6*IT_0163*IT_0170 + IT_0173*IT_0182 + (conj(IT_0094)*IT_0134 +
       IT_0094*IT_0137 + conj(IT_0092)*IT_0147 + IT_0092*IT_0150 + conj(IT_0144)
      *IT_0151 + conj(IT_0055)*IT_0154 + IT_0055*IT_0157 + conj(IT_0075)*IT_0180
       + conj(IT_0077)*IT_0181)*IT_0183 + IT_0185*(IT_0115*IT_0174 + IT_0122
      *IT_0176 + IT_0178*IT_0186 + IT_0175*IT_0187) + (IT_0122*IT_0174 + IT_0115
      *IT_0176 + IT_0175*IT_0186 + IT_0178*IT_0187)*IT_0189;
    return create_ccomplex_return(IT_0190);
}

