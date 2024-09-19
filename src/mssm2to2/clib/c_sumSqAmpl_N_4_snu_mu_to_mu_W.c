#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_snu_mu_to_mu_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_snu_mu_to_mu_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_mu, 2);
    const ccomplex_t IT_0001 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + 0.5*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em*m_mu
      *IT_0004*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0011*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0021 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0022 = IT_0004*IT_0021;
    const ccomplex_t IT_0023 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + 0.5*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em*m_mu
      *IT_0004*IT_0012*IT_0014;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0004*IT_0012*IT_0014;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = pow(m_N_4, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0034*IT_0040;
    const ccomplex_t IT_0042 = (-2)*IT_0019 + (-2)*IT_0032 + IT_0041;
    const ccomplex_t IT_0043 = pow(m_W, -2);
    const ccomplex_t IT_0044 = pow(s_14, 2);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = -IT_0037;
    const ccomplex_t IT_0047 = IT_0045 + IT_0046;
    const ccomplex_t IT_0048 = s_13*IT_0047;
    const ccomplex_t IT_0049 = cos(theta_W);
    const ccomplex_t IT_0050 = cpow(IT_0049, -1);
    const ccomplex_t IT_0051 = conj(N_B4)*e_em;
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W4)*e_em;
    const ccomplex_t IT_0055 = IT_0004*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + IT_0056);
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0038*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0034*IT_0060;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = m_mu*m_N_4;
    const ccomplex_t IT_0064 = IT_0047*IT_0063;
    const ccomplex_t IT_0065 = IT_0019 + IT_0032;
    const ccomplex_t IT_0066 = s_13*s_14;
    const ccomplex_t IT_0067 = s_34*IT_0037;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = 0.5*IT_0033;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0053 + -IT_0056);
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = cpow(s_12 + 0.5*IT_0001 + 0.5*IT_0037 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = pow(m_W, 2);
    const ccomplex_t IT_0079 = IT_0063*IT_0078;
    const ccomplex_t IT_0080 = 3*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0004;
    const ccomplex_t IT_0082 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0083 = IT_0004*IT_0082;
    const ccomplex_t IT_0084 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0085 = IT_0004*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + (-0.5)*IT_0086);
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0081*IT_0088;
    const ccomplex_t IT_0090 = IT_0002*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0004;
    const ccomplex_t IT_0093 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0094 = IT_0004*IT_0093;
    const ccomplex_t IT_0095 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0096 = IT_0004*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + (-0.5)*IT_0097);
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0092*IT_0099;
    const ccomplex_t IT_0101 = IT_0020*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0091 + IT_0102;
    const ccomplex_t IT_0104 = (-3)*IT_0079;
    const ccomplex_t IT_0105 = -IT_0041;
    const ccomplex_t IT_0106 = s_14*IT_0000;
    const ccomplex_t IT_0107 = s_13*s_34;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = IT_0106 + IT_0108;
    const ccomplex_t IT_0110 = s_13*IT_0078;
    const ccomplex_t IT_0111 = pow(m_W, 4);
    const ccomplex_t IT_0112 = s_13*IT_0111;
    const ccomplex_t IT_0113 = IT_0043*IT_0112;
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = s_14*s_34;
    const ccomplex_t IT_0116 = (-2)*IT_0115;
    const ccomplex_t IT_0117 = IT_0110 + IT_0114 + IT_0116;
    const ccomplex_t IT_0118 = s_14*IT_0078;
    const ccomplex_t IT_0119 = IT_0043*IT_0118;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = s_14 + IT_0120;
    const ccomplex_t IT_0122 = m_mu*IT_0121;
    const ccomplex_t IT_0123 = 2*IT_0122;
    const ccomplex_t IT_0124 = IT_0011*IT_0081;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0126 = IT_0002*IT_0124*IT_0125;
    const ccomplex_t IT_0127 = IT_0027*IT_0092;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0129 = IT_0020*IT_0127*IT_0128;
    const ccomplex_t IT_0130 = IT_0126 + IT_0129;
    const ccomplex_t IT_0131 = m_N_4*IT_0103;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0130 + IT_0132;
    const ccomplex_t IT_0134 = s_34*IT_0078;
    const ccomplex_t IT_0135 = IT_0043*IT_0134;
    const ccomplex_t IT_0136 = (-0.25)*IT_0135;
    const ccomplex_t IT_0137 = s_34 + IT_0136;
    const ccomplex_t IT_0138 = m_N_4*IT_0137;
    const ccomplex_t IT_0139 = (-4)*IT_0138;
    const ccomplex_t IT_0140 = m_N_4*IT_0019;
    const ccomplex_t IT_0141 = m_N_4*IT_0032;
    const ccomplex_t IT_0142 = m_mu*IT_0076;
    const ccomplex_t IT_0143 = IT_0140 + IT_0141 + 0.5*IT_0142;
    const ccomplex_t IT_0144 = IT_0016*IT_0088;
    const ccomplex_t IT_0145 = IT_0002*IT_0125*IT_0144;
    const ccomplex_t IT_0146 = IT_0029*IT_0099;
    const ccomplex_t IT_0147 = IT_0020*IT_0128*IT_0146;
    const ccomplex_t IT_0148 = -IT_0145 + -IT_0147;
    const ccomplex_t IT_0149 = IT_0143 + IT_0148;
    const ccomplex_t IT_0150 = (-6)*conj(IT_0130);
    const ccomplex_t IT_0151 = (-2)*IT_0091 + (-2)*IT_0102;
    const ccomplex_t IT_0152 = s_14*s_34*IT_0043;
    const ccomplex_t IT_0153 = -IT_0152;
    const ccomplex_t IT_0154 = s_13 + IT_0153;
    const ccomplex_t IT_0155 = m_N_4*IT_0154;
    const ccomplex_t IT_0156 = 0.5*IT_0135;
    const ccomplex_t IT_0157 = s_34 + IT_0156;
    const ccomplex_t IT_0158 = m_N_4*IT_0157;
    const ccomplex_t IT_0159 = 2*IT_0158;
    const ccomplex_t IT_0160 = IT_0110 + IT_0116;
    const ccomplex_t IT_0161 = IT_0043*IT_0160;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = s_13 + IT_0162;
    const ccomplex_t IT_0164 = m_mu*IT_0154;
    const ccomplex_t IT_0165 = (-2)*IT_0105;
    const ccomplex_t IT_0166 = IT_0063*IT_0154;
    const ccomplex_t IT_0167 = IT_0061 + -IT_0076;
    const ccomplex_t IT_0168 = pow(s_34, 2);
    const ccomplex_t IT_0169 = IT_0043*IT_0168;
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = IT_0000 + IT_0170;
    const ccomplex_t IT_0172 = IT_0063*IT_0171;
    const ccomplex_t IT_0173 = m_N_4*IT_0171;
    const ccomplex_t IT_0174 = (-2)*conj(IT_0105);
    const ccomplex_t IT_0175 = s_13*IT_0154;
    const ccomplex_t IT_0176 = (-2)*conj(IT_0042);
    const ccomplex_t IT_0177 = m_mu*IT_0047;
    const ccomplex_t IT_0178 = (-2)*conj(IT_0151);
    const ccomplex_t IT_0179 = (-0.25)*IT_0113;
    const ccomplex_t IT_0180 = IT_0115 + IT_0179;
    const ccomplex_t IT_0181 = (-2)*IT_0062;
    const ccomplex_t IT_0182 = 0.5*IT_0181;
    const ccomplex_t IT_0183 = IT_0151 + IT_0182;
    const ccomplex_t IT_0184 = (-0.25)*IT_0119;
    const ccomplex_t IT_0185 = s_14 + IT_0184;
    const ccomplex_t IT_0186 = m_mu*IT_0185;
    const ccomplex_t IT_0187 = (-4)*IT_0186;
    const ccomplex_t IT_0188 = (-2)*IT_0167;
    const ccomplex_t IT_0189 = IT_0077 + IT_0103;
    const ccomplex_t IT_0190 = (-2)*conj(IT_0167);
    const ccomplex_t IT_0191 = (-2)*IT_0130;
    const ccomplex_t IT_0192 = (-2)*conj(IT_0103)*(IT_0065*IT_0104 + 2*IT_0103
      *IT_0117 + IT_0123*IT_0149 + (-4)*IT_0077*IT_0180 + -IT_0069*IT_0183 + (
      -0.5)*IT_0109*IT_0188 + IT_0139*(IT_0132 + (-0.5)*IT_0191));
    const ccomplex_t IT_0193 = (-2)*conj(IT_0130);
    const ccomplex_t IT_0194 = 0.5*IT_0166;
    const ccomplex_t IT_0195 = (-2)*IT_0042;
    const ccomplex_t IT_0196 = s_13*IT_0171;
    const ccomplex_t IT_0197 = 0.5*IT_0196;
    const ccomplex_t IT_0198 = 4*IT_0143;
    const ccomplex_t IT_0199 = (-2)*IT_0151;
    const ccomplex_t IT_0200 = 0.333333333333333*IT_0103;
    const ccomplex_t IT_0201 = (-0.333333333333333)*IT_0065;
    const ccomplex_t IT_0202 = (-0.333333333333333)*IT_0077;
    const ccomplex_t IT_0203 = (-0.166666666666667)*IT_0164;
    const ccomplex_t IT_0204 = (-0.166666666666667)*IT_0173;
    const ccomplex_t IT_0205 = IT_0063*IT_0133 + (-0.166666666666667)*IT_0155
      *IT_0195 + (-0.666666666666667)*IT_0163*(IT_0148 + 0.25*IT_0198) + (
      -0.333333333333333)*IT_0177*(IT_0062 + 0.5*IT_0199) + IT_0123*IT_0200 +
       IT_0139*IT_0201 + IT_0187*IT_0202 + IT_0188*IT_0203 + IT_0165*IT_0204;
    const ccomplex_t IT_0206 = (-3)*IT_0063*IT_0149 + 2*IT_0133*IT_0163 +
       IT_0042*IT_0177 + IT_0155*IT_0183 + (-3)*IT_0139*IT_0200 + (-3)*IT_0123
      *IT_0201 + (-3)*IT_0159*IT_0202 + (-3)*IT_0165*IT_0203 + (-3)*IT_0188
      *IT_0204;
    const ccomplex_t IT_0207 = (-2)*conj(IT_0062);
    const ccomplex_t IT_0208 = IT_0069*IT_0189;
    const ccomplex_t IT_0209 = 2*conj(IT_0042)*(IT_0042*IT_0048 + IT_0062
      *IT_0064 + IT_0065*IT_0069) + 2*conj(IT_0065)*(IT_0042*IT_0069 + IT_0077
      *IT_0080 + -IT_0103*IT_0104 + IT_0105*IT_0109 + (-2)*IT_0065*IT_0117 +
       IT_0123*IT_0133 + IT_0139*IT_0149) + IT_0063*IT_0149*IT_0150 + 2*conj
      (IT_0130)*(IT_0065*IT_0123 + IT_0151*IT_0155 + IT_0077*IT_0159 + 2*IT_0133
      *IT_0163 + 0.5*IT_0164*IT_0165) + (IT_0133*IT_0164 + IT_0062*IT_0166 +
       IT_0167*IT_0172 + IT_0149*IT_0173)*IT_0174 + (IT_0064*IT_0151 + IT_0149
      *IT_0155 + IT_0166*IT_0167 + IT_0105*IT_0175)*IT_0176 + 2*IT_0177*(IT_0042
      *conj(IT_0130) + conj(IT_0042)*IT_0133 + IT_0149*(conj(IT_0062) + 0.5
      *IT_0178)) + 2*conj(IT_0077)*(IT_0065*IT_0080 + (-2)*IT_0077*IT_0117 +
       IT_0133*IT_0159 + 4*IT_0103*IT_0180 + IT_0069*IT_0183 + IT_0149*IT_0187 +
       0.5*IT_0109*IT_0188) + (IT_0149*IT_0164 + IT_0132*IT_0173 + IT_0109
      *IT_0189)*IT_0190 + IT_0192 + (IT_0103*IT_0139 + IT_0062*IT_0155 + IT_0167
      *IT_0173)*IT_0193 + 2*conj(IT_0062)*(IT_0048*IT_0062 + IT_0042*IT_0064 +
       0.5*IT_0155*IT_0191 + IT_0165*IT_0194) + 2*conj(IT_0105)*(IT_0065*IT_0109
       + IT_0151*IT_0166 + 0.5*IT_0175*IT_0195 + IT_0165*IT_0197) + 2*conj
      (IT_0167)*(0.5*IT_0165*IT_0172 + IT_0175*IT_0183 + 0.5*IT_0173*IT_0191 +
       IT_0194*IT_0195 + IT_0188*IT_0197) + (-6)*(conj(IT_0143) + conj(IT_0148))
      *IT_0205 + 2*conj(IT_0132)*IT_0206 + IT_0207*(IT_0048*IT_0151 + IT_0132
      *IT_0155 + IT_0167*IT_0175 + IT_0208) + 2*conj(IT_0151)*(IT_0133*IT_0155 +
       IT_0105*IT_0166 + IT_0167*IT_0175 + IT_0048*IT_0183 + 0.5*IT_0064*IT_0195
       + IT_0208);
    return create_ccomplex_return(IT_0209);
}

