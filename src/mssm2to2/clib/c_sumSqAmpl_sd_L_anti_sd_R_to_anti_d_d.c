#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sd_R_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sd_R_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = m_d*e_em*IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = m_d*e_em*mu_h*IT_0000;
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0009 + 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = pow(m_sd_L, 2);
    const ccomplex_t IT_0016 = pow(m_sd_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + -IT_0015 + -IT_0016 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = m_d*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0023 = IT_0008*IT_0022;
    const ccomplex_t IT_0024 = IT_0001*IT_0010;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0023 + (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0021*IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + -IT_0015 + -IT_0016 + IT_0028 
      + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0019 + -IT_0031;
    const ccomplex_t IT_0033 = pow(m_d, 2);
    const ccomplex_t IT_0034 = sin(alpha);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0034;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = cos(alpha);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0034 + m_d*e_em*mu_h*IT_0000*IT_0002*IT_0005*IT_0037);
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0015 + -IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0037;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*-1)*(m_d*e_em*mu_h*IT_0000
      *IT_0002*IT_0005*IT_0034 + (-1.4142135623731)*A_b*V_td*conj(V_td)*IT_0037);
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + -IT_0015 + -IT_0016 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = -IT_0043 + -IT_0051;
    const ccomplex_t IT_0053 = cos(theta_W);
    const ccomplex_t IT_0054 = cpow(IT_0053, -1);
    const ccomplex_t IT_0055 = conj(N_B1)*e_em;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W1)*e_em;
    const ccomplex_t IT_0059 = IT_0005*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + (-3)*IT_0060);
    const ccomplex_t IT_0062 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_23 + IT_0016 + IT_0033 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = conj(N_B2)*e_em;
    const ccomplex_t IT_0070 = IT_0054*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = conj(N_W2)*e_em;
    const ccomplex_t IT_0073 = IT_0005*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0071 + (-3)*IT_0074);
    const ccomplex_t IT_0076 = 0.166666666666667*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_23 + IT_0016 + IT_0033 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B3)*e_em;
    const ccomplex_t IT_0084 = IT_0054*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W3)*e_em;
    const ccomplex_t IT_0087 = IT_0005*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0085 + (-3)*IT_0088);
    const ccomplex_t IT_0090 = 0.166666666666667*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0016 + IT_0033 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = conj(N_B4)*e_em;
    const ccomplex_t IT_0098 = IT_0054*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = conj(N_W4)*e_em;
    const ccomplex_t IT_0101 = IT_0005*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0099 + (-3)*IT_0102);
    const ccomplex_t IT_0104 = 0.166666666666667*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0016 + IT_0033 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = -IT_0068 + -IT_0082 + -IT_0096 + -IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0054;
    const ccomplex_t IT_0113 = (-0.333333333333333)*IT_0112;
    const ccomplex_t IT_0114 = (-0.5)*IT_0091;
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = IT_0094*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0054;
    const ccomplex_t IT_0119 = (-0.333333333333333)*IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0105;
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = IT_0108*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0054;
    const ccomplex_t IT_0125 = (-0.333333333333333)*IT_0124;
    const ccomplex_t IT_0126 = (-0.5)*IT_0063;
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = IT_0066*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0054;
    const ccomplex_t IT_0131 = (-0.333333333333333)*IT_0130;
    const ccomplex_t IT_0132 = (-0.5)*IT_0077;
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = IT_0080*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = -IT_0117 + -IT_0123 + -IT_0129 + -IT_0135;
    const ccomplex_t IT_0137 = s_34*IT_0015;
    const ccomplex_t IT_0138 = (-18)*IT_0137;
    const ccomplex_t IT_0139 = s_13*s_14;
    const ccomplex_t IT_0140 = 36*IT_0139;
    const ccomplex_t IT_0141 = IT_0138 + IT_0140;
    const ccomplex_t IT_0142 = IT_0064*IT_0126;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0144 = IT_0066*IT_0142*IT_0143;
    const ccomplex_t IT_0145 = IT_0078*IT_0132;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0147 = IT_0080*IT_0145*IT_0146;
    const ccomplex_t IT_0148 = IT_0092*IT_0114;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0150 = IT_0094*IT_0148*IT_0149;
    const ccomplex_t IT_0151 = IT_0106*IT_0120;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0153 = IT_0108*IT_0151*IT_0152;
    const ccomplex_t IT_0154 = m_d*IT_0068;
    const ccomplex_t IT_0155 = m_d*IT_0082;
    const ccomplex_t IT_0156 = m_d*IT_0096;
    const ccomplex_t IT_0157 = m_d*IT_0110;
    const ccomplex_t IT_0158 = -IT_0144 + -IT_0147 + -IT_0150 + -IT_0153 +
       IT_0154 + IT_0155 + IT_0156 + IT_0157;
    const ccomplex_t IT_0159 = 18*s_34;
    const ccomplex_t IT_0160 = s_13 + s_14;
    const ccomplex_t IT_0161 = m_d*IT_0160;
    const ccomplex_t IT_0162 = 6*IT_0032;
    const ccomplex_t IT_0163 = 6*conj(IT_0032);
    const ccomplex_t IT_0164 = -s_14;
    const ccomplex_t IT_0165 = s_13 + IT_0164;
    const ccomplex_t IT_0166 = m_d*IT_0165;
    const ccomplex_t IT_0167 = 6*IT_0052;
    const ccomplex_t IT_0168 = 6*conj(IT_0052);
    const ccomplex_t IT_0169 = -IT_0033;
    const ccomplex_t IT_0170 = IT_0062*IT_0125;
    const ccomplex_t IT_0171 = IT_0066*IT_0143*IT_0170;
    const ccomplex_t IT_0172 = IT_0076*IT_0131;
    const ccomplex_t IT_0173 = IT_0080*IT_0146*IT_0172;
    const ccomplex_t IT_0174 = IT_0090*IT_0113;
    const ccomplex_t IT_0175 = IT_0094*IT_0149*IT_0174;
    const ccomplex_t IT_0176 = IT_0104*IT_0119;
    const ccomplex_t IT_0177 = IT_0108*IT_0152*IT_0176;
    const ccomplex_t IT_0178 = m_d*IT_0129;
    const ccomplex_t IT_0179 = m_d*IT_0135;
    const ccomplex_t IT_0180 = pow(m_sG, 2);
    const ccomplex_t IT_0181 = cpow((-2)*s_23 + IT_0016 + IT_0033 + -IT_0180 +
       reg_prop, -1);
    const ccomplex_t IT_0182 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0183 = cpow(IT_0182, 2);
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*m_sG*IT_0183;
    const ccomplex_t IT_0185 = IT_0181*IT_0184;
    const ccomplex_t IT_0186 = -IT_0171 + -IT_0173 + -IT_0175 + -IT_0177 +
       IT_0178 + IT_0179 + 0.166666666666667*IT_0185;
    const ccomplex_t IT_0187 = IT_0117 + IT_0123;
    const ccomplex_t IT_0188 = m_d*IT_0187;
    const ccomplex_t IT_0189 = conj(IT_0186) + conj(IT_0188);
    const ccomplex_t IT_0190 = IT_0159*IT_0186;
    const ccomplex_t IT_0191 = (-18)*IT_0158;
    const ccomplex_t IT_0192 = (-18)*conj(IT_0158);
    const ccomplex_t IT_0193 = IT_0015*IT_0033;
    const ccomplex_t IT_0194 = (-18)*IT_0111;
    const ccomplex_t IT_0195 = (-18)*conj(IT_0111);
    const ccomplex_t IT_0196 = (-0.5)*IT_0185;
    const ccomplex_t IT_0197 = (-6)*s_34*IT_0196;
    const ccomplex_t IT_0198 = m_d*s_14;
    const ccomplex_t IT_0199 = 6*IT_0196;
    const ccomplex_t IT_0200 = (-6)*IT_0196;
    const ccomplex_t IT_0201 = m_d*s_13;
    const ccomplex_t IT_0202 = 18*IT_0033;
    const ccomplex_t IT_0203 = IT_0159 + IT_0202;
    const ccomplex_t IT_0204 = (-18)*IT_0033;
    const ccomplex_t IT_0205 = IT_0159 + IT_0204;
    const ccomplex_t IT_0206 = (-6)*conj(IT_0158);
    const ccomplex_t IT_0207 = 0.166666666666667*IT_0206;
    const ccomplex_t IT_0208 = 6*s_34*(conj(IT_0186) + conj(IT_0188) + (-3)
      *IT_0196) + -(IT_0032 + -conj(IT_0032))*IT_0203 + -(IT_0052 + -conj
      (IT_0052))*IT_0205 + 6*IT_0033*(IT_0158 + IT_0207);
    const ccomplex_t IT_0209 = s_34 + IT_0033;
    const ccomplex_t IT_0210 = IT_0032*conj(IT_0032)*(36*s_34 + 36*IT_0033) + 
      (36*s_34 + (-36)*IT_0033)*IT_0052*conj(IT_0052) + (IT_0111*conj(IT_0111) +
       IT_0136*conj(IT_0136))*IT_0141 + IT_0158*conj(IT_0158)*IT_0159 + (-6)
      *IT_0161*(conj(IT_0032)*IT_0111 + IT_0032*conj(IT_0111) + (
      -0.166666666666667)*conj(IT_0136)*IT_0162 + (-0.166666666666667)*IT_0136
      *IT_0163) + IT_0166*((conj(IT_0111) + conj(IT_0136))*IT_0167 + (IT_0111 +
       IT_0136)*IT_0168) + (s_34 + IT_0169)*(IT_0168*(IT_0158 + IT_0186 +
       IT_0188) + IT_0167*(conj(IT_0158) + conj(IT_0186) + conj(IT_0188))) +
       IT_0189*IT_0190 + IT_0033*(IT_0189*IT_0191 + IT_0186*IT_0192) + IT_0193*
      (conj(IT_0136)*IT_0194 + IT_0136*IT_0195) + IT_0188*(IT_0159*IT_0189 +
       IT_0033*IT_0192 + IT_0197) + IT_0186*IT_0197 + IT_0198*(conj(IT_0136)
      *IT_0191 + IT_0136*IT_0192 + conj(IT_0188)*IT_0194 + IT_0188*IT_0195 +
       IT_0111*((-18)*conj(IT_0186) + IT_0199) + conj(IT_0111)*((-18)*IT_0186 +
       IT_0200)) + (18*conj(IT_0111)*IT_0158 + 18*IT_0111*conj(IT_0158) + conj
      (IT_0136)*(18*IT_0186 + 18*IT_0188 + IT_0199) + IT_0136*(18*conj(IT_0186) 
      + 18*conj(IT_0188) + IT_0200))*IT_0201 + IT_0196*IT_0208 + (-6)*(conj
      (IT_0032)*(IT_0158 + -IT_0186) + (-0.166666666666667)*IT_0163*IT_0188 + (
      -0.166666666666667)*IT_0162*conj(IT_0188) + -IT_0032*(conj(IT_0186) +
       IT_0207))*IT_0209;
    return create_ccomplex_return(IT_0210);
}

