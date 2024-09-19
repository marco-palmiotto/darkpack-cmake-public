#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sd_L_to_u_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sd_L_to_u_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
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
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_sd_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + (-0.5)*IT_0014);
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0024 = IT_0007*IT_0023;
    const ccomplex_t IT_0025 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0026 = IT_0007*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + (-0.5)*IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0007;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = pow(m_N_1, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0034*IT_0040;
    const ccomplex_t IT_0042 = 2*IT_0019 + 2*IT_0032 + IT_0041;
    const ccomplex_t IT_0043 = pow(m_W, -2);
    const ccomplex_t IT_0044 = pow(s_14, 2);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = -IT_0037;
    const ccomplex_t IT_0047 = IT_0045 + IT_0046;
    const ccomplex_t IT_0048 = s_13*IT_0047;
    const ccomplex_t IT_0049 = cos(theta_W);
    const ccomplex_t IT_0050 = cpow(IT_0049, -1);
    const ccomplex_t IT_0051 = conj(N_B1)*e_em;
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W1)*e_em;
    const ccomplex_t IT_0055 = IT_0007*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 3*IT_0056);
    const ccomplex_t IT_0058 = 0.166666666666667*IT_0057;
    const ccomplex_t IT_0059 = IT_0038*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0034*IT_0060;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = m_u*m_N_1;
    const ccomplex_t IT_0064 = IT_0047*IT_0063;
    const ccomplex_t IT_0065 = -IT_0019 + -IT_0032;
    const ccomplex_t IT_0066 = s_13*s_14;
    const ccomplex_t IT_0067 = s_34*IT_0037;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = m_u*IT_0047;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0007;
    const ccomplex_t IT_0073 = IT_0016*IT_0072;
    const ccomplex_t IT_0074 = IT_0002*IT_0071*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0007;
    const ccomplex_t IT_0077 = IT_0029*IT_0076;
    const ccomplex_t IT_0078 = IT_0020*IT_0075*IT_0077;
    const ccomplex_t IT_0079 = pow(m_d, 2);
    const ccomplex_t IT_0080 = cpow((-2)*s_12 + -IT_0001 + -IT_0037 + IT_0079 
      + -reg_prop, -1);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0082 = 0.5*IT_0033;
    const ccomplex_t IT_0083 = cos(beta);
    const ccomplex_t IT_0084 = cpow(IT_0083, -1);
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0003*IT_0007*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0082*IT_0086;
    const ccomplex_t IT_0088 = IT_0080*IT_0081*IT_0087;
    const ccomplex_t IT_0089 = -IT_0074 + -IT_0078 + -IT_0088;
    const ccomplex_t IT_0090 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0091 = IT_0007*IT_0090;
    const ccomplex_t IT_0092 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0093 = IT_0007*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + 0.5*IT_0094);
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = IT_0076*IT_0096;
    const ccomplex_t IT_0098 = IT_0020*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0101 = IT_0007*IT_0100;
    const ccomplex_t IT_0102 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0103 = IT_0007*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + 0.5*IT_0104);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0072*IT_0106;
    const ccomplex_t IT_0108 = IT_0002*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0099 + IT_0109;
    const ccomplex_t IT_0111 = m_N_1*IT_0110;
    const ccomplex_t IT_0112 = IT_0089 + IT_0111;
    const ccomplex_t IT_0113 = -IT_0099 + -IT_0109;
    const ccomplex_t IT_0114 = pow(m_W, 2);
    const ccomplex_t IT_0115 = IT_0063*IT_0114;
    const ccomplex_t IT_0116 = (-3)*IT_0115;
    const ccomplex_t IT_0117 = 3*IT_0115;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0053 + (-3)*IT_0056);
    const ccomplex_t IT_0119 = 0.166666666666667*IT_0118;
    const ccomplex_t IT_0120 = IT_0082*IT_0119;
    const ccomplex_t IT_0121 = IT_0080*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = -IT_0041;
    const ccomplex_t IT_0124 = s_14*IT_0000;
    const ccomplex_t IT_0125 = s_13*s_34;
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = IT_0124 + IT_0126;
    const ccomplex_t IT_0128 = pow(m_W, 4);
    const ccomplex_t IT_0129 = s_13*IT_0128;
    const ccomplex_t IT_0130 = IT_0043*IT_0129;
    const ccomplex_t IT_0131 = s_14*s_34;
    const ccomplex_t IT_0132 = 24*IT_0131;
    const ccomplex_t IT_0133 = s_13*IT_0114;
    const ccomplex_t IT_0134 = (-12)*IT_0133;
    const ccomplex_t IT_0135 = 6*IT_0130 + IT_0132 + IT_0134;
    const ccomplex_t IT_0136 = s_14*IT_0114;
    const ccomplex_t IT_0137 = IT_0043*IT_0136;
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = s_14 + IT_0138;
    const ccomplex_t IT_0140 = m_u*IT_0139;
    const ccomplex_t IT_0141 = 2*IT_0140;
    const ccomplex_t IT_0142 = s_34*IT_0114;
    const ccomplex_t IT_0143 = IT_0043*IT_0142;
    const ccomplex_t IT_0144 = (-0.25)*IT_0143;
    const ccomplex_t IT_0145 = s_34 + IT_0144;
    const ccomplex_t IT_0146 = m_N_1*IT_0145;
    const ccomplex_t IT_0147 = (-4)*IT_0146;
    const ccomplex_t IT_0148 = m_N_1*IT_0019;
    const ccomplex_t IT_0149 = m_N_1*IT_0032;
    const ccomplex_t IT_0150 = m_u*IT_0122;
    const ccomplex_t IT_0151 = -IT_0148 + -IT_0149 + -IT_0150;
    const ccomplex_t IT_0152 = IT_0009*IT_0106;
    const ccomplex_t IT_0153 = IT_0002*IT_0071*IT_0152;
    const ccomplex_t IT_0154 = IT_0022*IT_0096;
    const ccomplex_t IT_0155 = IT_0020*IT_0075*IT_0154;
    const ccomplex_t IT_0156 = IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = IT_0151 + IT_0156;
    const ccomplex_t IT_0158 = 6*conj(IT_0065)*(IT_0042*IT_0069 + -IT_0113
      *IT_0116 + IT_0117*IT_0122 + IT_0123*IT_0127 + 0.166666666666667*IT_0065
      *IT_0135 + IT_0112*IT_0141 + IT_0147*IT_0157);
    const ccomplex_t IT_0159 = 2*IT_0099 + 2*IT_0109;
    const ccomplex_t IT_0160 = s_14*s_34*IT_0043;
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = s_13 + IT_0161;
    const ccomplex_t IT_0163 = m_N_1*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0143;
    const ccomplex_t IT_0165 = s_34 + IT_0164;
    const ccomplex_t IT_0166 = m_N_1*IT_0165;
    const ccomplex_t IT_0167 = 2*IT_0166;
    const ccomplex_t IT_0168 = (-2)*IT_0131;
    const ccomplex_t IT_0169 = IT_0133 + IT_0168;
    const ccomplex_t IT_0170 = IT_0043*IT_0169;
    const ccomplex_t IT_0171 = (-6)*IT_0170;
    const ccomplex_t IT_0172 = 12*s_13;
    const ccomplex_t IT_0173 = IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = m_u*IT_0162;
    const ccomplex_t IT_0175 = (-6)*IT_0123;
    const ccomplex_t IT_0176 = IT_0063*IT_0162;
    const ccomplex_t IT_0177 = IT_0061 + 2*IT_0122;
    const ccomplex_t IT_0178 = pow(s_34, 2);
    const ccomplex_t IT_0179 = IT_0043*IT_0178;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0000 + IT_0180;
    const ccomplex_t IT_0182 = IT_0063*IT_0181;
    const ccomplex_t IT_0183 = m_N_1*IT_0181;
    const ccomplex_t IT_0184 = (-6)*conj(IT_0123);
    const ccomplex_t IT_0185 = (-6)*conj(IT_0089);
    const ccomplex_t IT_0186 = (-18)*conj(IT_0089);
    const ccomplex_t IT_0187 = (-6)*conj(IT_0159);
    const ccomplex_t IT_0188 = s_13*IT_0162;
    const ccomplex_t IT_0189 = (-6)*conj(IT_0042);
    const ccomplex_t IT_0190 = (-6)*IT_0130;
    const ccomplex_t IT_0191 = IT_0132 + IT_0190;
    const ccomplex_t IT_0192 = (-6)*IT_0062;
    const ccomplex_t IT_0193 = 0.166666666666667*IT_0192;
    const ccomplex_t IT_0194 = IT_0159 + IT_0193;
    const ccomplex_t IT_0195 = (-0.25)*IT_0137;
    const ccomplex_t IT_0196 = s_14 + IT_0195;
    const ccomplex_t IT_0197 = m_u*IT_0196;
    const ccomplex_t IT_0198 = (-4)*IT_0197;
    const ccomplex_t IT_0199 = (-6)*IT_0177;
    const ccomplex_t IT_0200 = (-6)*IT_0089;
    const ccomplex_t IT_0201 = (-6)*conj(IT_0113)*(IT_0065*IT_0116 + (
      -0.166666666666667)*IT_0113*IT_0135 + IT_0141*IT_0157 + (
      -0.166666666666667)*IT_0122*IT_0191 + -IT_0069*IT_0194 + (
      -0.166666666666667)*IT_0127*IT_0199 + IT_0147*(IT_0111 + (
      -0.166666666666667)*IT_0200));
    const ccomplex_t IT_0202 = IT_0113 + IT_0122;
    const ccomplex_t IT_0203 = (-6)*conj(IT_0177);
    const ccomplex_t IT_0204 = (-6)*IT_0159;
    const ccomplex_t IT_0205 = (-6)*IT_0042;
    const ccomplex_t IT_0206 = 6*IT_0065;
    const ccomplex_t IT_0207 = 6*IT_0122;
    const ccomplex_t IT_0208 = (-6)*IT_0113;
    const ccomplex_t IT_0209 = IT_0063*((-18)*IT_0089 + (-18)*IT_0111) +
       IT_0157*IT_0173 + IT_0175*IT_0183 + IT_0174*IT_0199 + 6*IT_0070*(IT_0062 
      + 0.166666666666667*IT_0204) + IT_0163*IT_0205 + IT_0147*IT_0206 + IT_0198
      *IT_0207 + IT_0141*IT_0208;
    const ccomplex_t IT_0210 = 6*IT_0042*IT_0070 + IT_0063*((-18)*IT_0151 + (
      -18)*IT_0156) + IT_0112*IT_0173 + IT_0174*IT_0175 + 6*IT_0163*IT_0194 +
       IT_0183*IT_0199 + IT_0141*IT_0206 + IT_0167*IT_0207 + IT_0147*IT_0208;
    const ccomplex_t IT_0211 = 0.166666666666667*IT_0176;
    const ccomplex_t IT_0212 = s_13*IT_0181;
    const ccomplex_t IT_0213 = 0.166666666666667*IT_0212;
    const ccomplex_t IT_0214 = (-6)*conj(IT_0062);
    const ccomplex_t IT_0215 = IT_0177*IT_0188;
    const ccomplex_t IT_0216 = IT_0069*IT_0202;
    const ccomplex_t IT_0217 = 6*conj(IT_0042)*(IT_0042*IT_0048 + IT_0062
      *IT_0064 + IT_0065*IT_0069 + IT_0070*IT_0112) + IT_0158 + 6*conj(IT_0089)*
      (IT_0042*IT_0070 + IT_0065*IT_0141 + IT_0159*IT_0163 + IT_0122*IT_0167 +
       0.166666666666667*IT_0112*IT_0173 + 0.166666666666667*IT_0174*IT_0175) + 
      (IT_0112*IT_0174 + IT_0062*IT_0176 + IT_0177*IT_0182 + IT_0157*IT_0183)
      *IT_0184 + (IT_0113*IT_0147 + IT_0062*IT_0163 + IT_0177*IT_0183)*IT_0185 +
       IT_0157*(IT_0063*IT_0186 + IT_0070*IT_0187) + (IT_0064*IT_0159 + IT_0157
      *IT_0163 + IT_0176*IT_0177 + IT_0123*IT_0188)*IT_0189 + 6*conj(IT_0122)*
      (IT_0065*IT_0117 + 0.166666666666667*IT_0122*IT_0135 + IT_0112*IT_0167 +
       0.166666666666667*IT_0113*IT_0191 + IT_0069*IT_0194 + IT_0157*IT_0198 +
       0.166666666666667*IT_0127*IT_0199) + IT_0201 + (IT_0157*IT_0174 + IT_0111
      *IT_0183 + IT_0127*IT_0202)*IT_0203 + (conj(IT_0151) + conj(IT_0156))
      *IT_0209 + conj(IT_0111)*IT_0210 + 6*conj(IT_0062)*(IT_0048*IT_0062 +
       IT_0042*IT_0064 + IT_0070*IT_0157 + 0.166666666666667*IT_0163*IT_0200 +
       IT_0175*IT_0211) + 6*conj(IT_0123)*(IT_0065*IT_0127 + IT_0159*IT_0176 +
       0.166666666666667*IT_0188*IT_0205 + IT_0175*IT_0213) + 6*conj(IT_0177)*
      (0.166666666666667*IT_0175*IT_0182 + IT_0188*IT_0194 + 0.166666666666667
      *IT_0183*IT_0200 + IT_0205*IT_0211 + IT_0199*IT_0213) + IT_0214*(IT_0048
      *IT_0159 + IT_0111*IT_0163 + IT_0215 + IT_0216) + 6*conj(IT_0159)*(IT_0112
      *IT_0163 + IT_0123*IT_0176 + IT_0048*IT_0194 + 0.166666666666667*IT_0064
      *IT_0205 + IT_0215 + IT_0216);
    return create_ccomplex_return(IT_0217);
}

