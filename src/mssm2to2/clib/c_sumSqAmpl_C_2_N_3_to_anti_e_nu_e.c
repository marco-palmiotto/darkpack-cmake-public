#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_anti_e_nu_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_anti_e_nu_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0003*IT_0005;
    const ccomplex_t IT_0009 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = (-2)*IT_0010;
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0019 = IT_0008*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = IT_0017 + IT_0020;
    const ccomplex_t IT_0022 = IT_0011 + IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0007*IT_0024;
    const ccomplex_t IT_0026 = pow(m_C_2, 2);
    const ccomplex_t IT_0027 = pow(m_N_3, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + -IT_0026 + -IT_0027 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + -IT_0026 + -IT_0027 + IT_0031 
      + -reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0000*IT_0005;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0001*IT_0005;
    const ccomplex_t IT_0036 = IT_0009*IT_0035;
    const ccomplex_t IT_0037 = (-2)*IT_0036;
    const ccomplex_t IT_0038 = IT_0001*IT_0013;
    const ccomplex_t IT_0039 = IT_0015*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = IT_0018*IT_0035;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = IT_0040 + IT_0042;
    const ccomplex_t IT_0044 = IT_0037 + IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = IT_0034*IT_0046;
    const ccomplex_t IT_0048 = IT_0032*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0030 + -IT_0049;
    const ccomplex_t IT_0051 = pow(m_e, 2);
    const ccomplex_t IT_0052 = cpow((-2)*s_13 + IT_0026 + IT_0051 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0053 = conj(N_B3)*e_em;
    const ccomplex_t IT_0054 = IT_0013*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = conj(N_W3)*e_em;
    const ccomplex_t IT_0057 = IT_0005*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + -IT_0058);
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = s_12*s_34;
    const ccomplex_t IT_0067 = s_13*s_24;
    const ccomplex_t IT_0068 = s_14*s_23;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = IT_0067 + IT_0069;
    const ccomplex_t IT_0071 = IT_0066 + IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0005;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_23 + IT_0027 + IT_0051 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_23 + IT_0027 + IT_0051 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0080 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0013;
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = IT_0062*IT_0081;
    const ccomplex_t IT_0083 = IT_0079*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0078 + IT_0084;
    const ccomplex_t IT_0086 = -IT_0066;
    const ccomplex_t IT_0087 = IT_0070 + IT_0086;
    const ccomplex_t IT_0088 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0089 = IT_0035*IT_0088;
    const ccomplex_t IT_0090 = 2*IT_0089;
    const ccomplex_t IT_0091 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0092 = IT_0038*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0095 = IT_0035*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = IT_0093 + IT_0096;
    const ccomplex_t IT_0098 = IT_0090 + IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0007*IT_0100;
    const ccomplex_t IT_0102 = IT_0028*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0008*IT_0088;
    const ccomplex_t IT_0105 = 2*IT_0104;
    const ccomplex_t IT_0106 = IT_0014*IT_0091;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = IT_0008*IT_0094;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = IT_0107 + IT_0109;
    const ccomplex_t IT_0111 = IT_0105 + IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0034*IT_0113;
    const ccomplex_t IT_0115 = IT_0032*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0103 + IT_0116;
    const ccomplex_t IT_0118 = 4*IT_0117;
    const ccomplex_t IT_0119 = 4*IT_0050;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = IT_0060*IT_0121;
    const ccomplex_t IT_0123 = IT_0052*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0129 = IT_0005*IT_0128;
    const ccomplex_t IT_0130 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0131 = IT_0005*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + (-0.5)*IT_0132);
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = IT_0127*IT_0134;
    const ccomplex_t IT_0136 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0026 + -IT_0027 + -reg_prop, -1);
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = 4*IT_0125;
    const ccomplex_t IT_0140 = 4*IT_0065;
    const ccomplex_t IT_0141 = m_e*s_14*m_N_3;
    const ccomplex_t IT_0142 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0143 = IT_0005*IT_0142;
    const ccomplex_t IT_0144 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0145 = IT_0005*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + 0.5*IT_0146);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = IT_0127*IT_0148;
    const ccomplex_t IT_0150 = IT_0136*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = 2*IT_0125;
    const ccomplex_t IT_0153 = 2*conj(IT_0125);
    const ccomplex_t IT_0154 = N_B3*e_em;
    const ccomplex_t IT_0155 = IT_0013*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = N_W3*e_em;
    const ccomplex_t IT_0158 = IT_0005*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0156 + IT_0159);
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = IT_0072*IT_0161;
    const ccomplex_t IT_0163 = IT_0076*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = IT_0062*IT_0166;
    const ccomplex_t IT_0168 = IT_0079*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = -IT_0164 + -IT_0169;
    const ccomplex_t IT_0171 = 2*IT_0170;
    const ccomplex_t IT_0172 = 2*conj(IT_0170);
    const ccomplex_t IT_0173 = 4*conj(IT_0170);
    const ccomplex_t IT_0174 = 4*IT_0170;
    const ccomplex_t IT_0175 = 0.125*conj(IT_0085);
    const ccomplex_t IT_0176 = 0.125*IT_0118*conj(IT_0138) + (0.5*conj(IT_0050
      ) + 0.5*conj(IT_0065) + conj(IT_0085))*IT_0151 + (IT_0085 + 0.125*IT_0119 
      + 0.125*IT_0140)*conj(IT_0151) + 0.5*conj(IT_0117)*(IT_0138 + 0.25*IT_0152
      ) + 0.125*IT_0117*IT_0153 + 0.125*(conj(IT_0050) + conj(IT_0065))*IT_0171 
      + 0.125*(IT_0050 + IT_0065)*IT_0172 + 0.125*IT_0085*IT_0173 + IT_0174
      *IT_0175;
    const ccomplex_t IT_0177 = 0.125*IT_0174;
    const ccomplex_t IT_0178 = IT_0151 + IT_0177;
    const ccomplex_t IT_0179 = 4*IT_0085;
    const ccomplex_t IT_0180 = m_e*s_24*m_C_2;
    const ccomplex_t IT_0181 = 4*conj(IT_0125);
    const ccomplex_t IT_0182 = (-4)*conj(IT_0117);
    const ccomplex_t IT_0183 = (-2)*IT_0117;
    const ccomplex_t IT_0184 = (-2)*conj(IT_0117);
    const ccomplex_t IT_0185 = (-0.5)*conj(IT_0065);
    const ccomplex_t IT_0186 = conj(IT_0050)*IT_0125 + IT_0050*conj(IT_0125) +
       2*(conj(IT_0050) + (-2)*conj(IT_0065) + -conj(IT_0085))*IT_0138 + 2
      *IT_0117*conj(IT_0151) + (-0.5)*conj(IT_0085)*IT_0152 + (-0.5)*IT_0085
      *IT_0153 + 2*conj(IT_0138)*(IT_0050 + (-2)*IT_0065 + (-0.25)*IT_0179) + (
      -0.5)*IT_0065*IT_0181 + (-0.5)*IT_0151*IT_0182 + (-0.5)*conj(IT_0170)
      *IT_0183 + (-0.5)*IT_0170*IT_0184 + IT_0139*IT_0185;
    const ccomplex_t IT_0187 = s_34*m_C_2*m_N_3;
    const ccomplex_t IT_0188 = IT_0065*conj(IT_0117) + 0.5*IT_0139*conj
      (IT_0151) + 0.5*IT_0153*IT_0170 + 0.5*IT_0152*conj(IT_0170) + 4*IT_0138*
      (conj(IT_0151) + 0.125*IT_0173) + 4*conj(IT_0138)*IT_0178 + 0.5*IT_0151
      *IT_0181 + 0.5*IT_0050*IT_0182 + 0.5*conj(IT_0085)*IT_0183 + 0.5*IT_0085
      *IT_0184 + (-2)*IT_0117*(conj(IT_0050) + IT_0185);
    const ccomplex_t IT_0189 = 2*(conj(IT_0065)*IT_0085 + IT_0065*conj(IT_0085
      ))*(IT_0067 + IT_0068 + IT_0086) + (-2)*IT_0050*(conj(IT_0065)*IT_0071 +
       conj(IT_0085)*IT_0087) + IT_0066*conj(IT_0117)*IT_0118 + (-2)*conj
      (IT_0050)*(IT_0065*IT_0071 + IT_0085*IT_0087 + (-0.5)*IT_0066*IT_0119) + 8
      *IT_0067*((IT_0125 + 2*IT_0138)*conj(IT_0138) + conj(IT_0125)*(IT_0138 +
       0.125*IT_0139) + 0.125*conj(IT_0065)*IT_0140) + 8*IT_0141*IT_0176 + 16
      *IT_0068*(conj(IT_0151)*(IT_0151 + 1./2*IT_0170) + 1./2*conj(IT_0170)
      *IT_0178 + 1./2*IT_0175*IT_0179) + (-2)*IT_0180*IT_0186 + 2*IT_0187*IT_0188;
    return create_ccomplex_return(IT_0189);
}

