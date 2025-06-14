#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sd_L_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sd_L_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_d, 2);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0003*IT_0005*IT_0008 + (-0.5)*m_W*(IT_0001*IT_0003 + -IT_0009*IT_0010)
      *(IT_0005 + 1./3*IT_0004*IT_0012));
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_sd_L, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0010;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0005*IT_0008*IT_0010 + (-0.5)*m_W*(IT_0003*IT_0009 + IT_0001*IT_0010)*
      (IT_0005 + 1./3*IT_0004*IT_0012));
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0011, -1);
    const ccomplex_t IT_0029 = IT_0004*IT_0028;
    const ccomplex_t IT_0030 = e_em*IT_0029;
    const ccomplex_t IT_0031 = IT_0005*IT_0011;
    const ccomplex_t IT_0032 = e_em*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0030 + 3*IT_0032);
    const ccomplex_t IT_0034 = (-0.166666666666667)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*3)*e_em*(IT_0005*IT_0011 + 1./3
      *IT_0004*IT_0028);
    const ccomplex_t IT_0039 = 0.166666666666667*IT_0038;
    const ccomplex_t IT_0040 = m_d*IT_0039;
    const ccomplex_t IT_0041 = 0.333333333333333*IT_0038;
    const ccomplex_t IT_0042 = m_d*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0040 + IT_0043;
    const ccomplex_t IT_0045 = IT_0037*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0028;
    const ccomplex_t IT_0047 = 0.333333333333333*IT_0046;
    const ccomplex_t IT_0048 = IT_0035*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0044*IT_0049;
    const ccomplex_t IT_0051 = pow(m_sG, 2);
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -IT_0051 +
       reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0054 = cpow(IT_0053, 2);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0052*IT_0055;
    const ccomplex_t IT_0057 = m_d*IT_0056;
    const ccomplex_t IT_0058 = -IT_0045 + IT_0050 + 0.5*IT_0057;
    const ccomplex_t IT_0059 = (IT_0040 + -IT_0042)*(IT_0037 + -IT_0049);
    const ccomplex_t IT_0060 = (-18)*IT_0008;
    const ccomplex_t IT_0061 = 18*s_34;
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = conj(N_B1)*e_em;
    const ccomplex_t IT_0064 = IT_0028*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W1)*e_em;
    const ccomplex_t IT_0067 = IT_0005*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + (-3)*IT_0068);
    const ccomplex_t IT_0070 = 0.166666666666667*IT_0069;
    const ccomplex_t IT_0071 = N_B1*e_em;
    const ccomplex_t IT_0072 = IT_0028*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = N_W1*e_em;
    const ccomplex_t IT_0075 = IT_0005*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + (-3)*IT_0076);
    const ccomplex_t IT_0078 = 0.166666666666667*IT_0077;
    const ccomplex_t IT_0079 = IT_0070*IT_0078;
    const ccomplex_t IT_0080 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B2)*e_em;
    const ccomplex_t IT_0084 = IT_0028*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W2)*e_em;
    const ccomplex_t IT_0087 = IT_0005*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0085 + (-3)*IT_0088);
    const ccomplex_t IT_0090 = 0.166666666666667*IT_0089;
    const ccomplex_t IT_0091 = N_B2*e_em;
    const ccomplex_t IT_0092 = IT_0028*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = N_W2*e_em;
    const ccomplex_t IT_0095 = IT_0005*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + (-3)*IT_0096);
    const ccomplex_t IT_0098 = 0.166666666666667*IT_0097;
    const ccomplex_t IT_0099 = IT_0090*IT_0098;
    const ccomplex_t IT_0100 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0101 = IT_0099*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = conj(N_B3)*e_em;
    const ccomplex_t IT_0104 = IT_0028*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = conj(N_W3)*e_em;
    const ccomplex_t IT_0107 = IT_0005*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + (-3)*IT_0108);
    const ccomplex_t IT_0110 = 0.166666666666667*IT_0109;
    const ccomplex_t IT_0111 = N_B3*e_em;
    const ccomplex_t IT_0112 = IT_0028*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = N_W3*e_em;
    const ccomplex_t IT_0115 = IT_0005*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*(IT_0113 + (-3)*IT_0116);
    const ccomplex_t IT_0118 = 0.166666666666667*IT_0117;
    const ccomplex_t IT_0119 = IT_0110*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B4)*e_em;
    const ccomplex_t IT_0124 = IT_0028*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = conj(N_W4)*e_em;
    const ccomplex_t IT_0127 = IT_0005*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + (-3)*IT_0128);
    const ccomplex_t IT_0130 = 0.166666666666667*IT_0129;
    const ccomplex_t IT_0131 = N_B4*e_em;
    const ccomplex_t IT_0132 = IT_0028*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = N_W4*e_em;
    const ccomplex_t IT_0135 = IT_0005*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + (-3)*IT_0136);
    const ccomplex_t IT_0138 = 0.166666666666667*IT_0137;
    const ccomplex_t IT_0139 = IT_0130*IT_0138;
    const ccomplex_t IT_0140 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = 0.166666666666667*IT_0056 + -IT_0082 + -IT_0102
       + -IT_0122 + -IT_0142;
    const ccomplex_t IT_0144 = IT_0037*IT_0041;
    const ccomplex_t IT_0145 = (-0.5)*IT_0056 + -IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0147 = (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = IT_0147*IT_0149;
    const ccomplex_t IT_0151 = IT_0080*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = IT_0154*IT_0156;
    const ccomplex_t IT_0158 = IT_0100*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = IT_0161*IT_0163;
    const ccomplex_t IT_0165 = IT_0120*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0168 = (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = IT_0168*IT_0170;
    const ccomplex_t IT_0172 = IT_0140*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = -IT_0152 + -IT_0159 + -IT_0166 + -IT_0173;
    const ccomplex_t IT_0175 = IT_0041*IT_0049;
    const ccomplex_t IT_0176 = -IT_0175;
    const ccomplex_t IT_0177 = s_34*IT_0016;
    const ccomplex_t IT_0178 = (-6)*IT_0177;
    const ccomplex_t IT_0179 = s_13*s_14;
    const ccomplex_t IT_0180 = 12*IT_0179;
    const ccomplex_t IT_0181 = IT_0178 + IT_0180;
    const ccomplex_t IT_0182 = (-18)*IT_0177;
    const ccomplex_t IT_0183 = 36*IT_0179;
    const ccomplex_t IT_0184 = IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0186 = 0.666666666666667*IT_0185;
    const ccomplex_t IT_0187 = (-0.333333333333333)*IT_0185;
    const ccomplex_t IT_0188 = IT_0186*IT_0187;
    const ccomplex_t IT_0189 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0190 = IT_0188*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0193 = (-2)*IT_0192;
    const ccomplex_t IT_0194 = IT_0189*IT_0192;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = IT_0193*IT_0195;
    const ccomplex_t IT_0197 = 0.5*IT_0191 + (-0.0833333333333333)*IT_0196;
    const ccomplex_t IT_0198 = cpow(IT_0197, 2);
    const ccomplex_t IT_0199 = (-36)*IT_0177;
    const ccomplex_t IT_0200 = IT_0008*IT_0016;
    const ccomplex_t IT_0201 = (-36)*IT_0200;
    const ccomplex_t IT_0202 = 72*IT_0179;
    const ccomplex_t IT_0203 = IT_0199 + IT_0201 + IT_0202;
    const ccomplex_t IT_0204 = 0.25*IT_0196;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0206 = IT_0070*IT_0147;
    const ccomplex_t IT_0207 = IT_0080*IT_0205*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0209 = IT_0090*IT_0154;
    const ccomplex_t IT_0210 = IT_0100*IT_0208*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0212 = IT_0110*IT_0161;
    const ccomplex_t IT_0213 = IT_0120*IT_0211*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0215 = IT_0130*IT_0168;
    const ccomplex_t IT_0216 = IT_0140*IT_0214*IT_0215;
    const ccomplex_t IT_0217 = m_d*IT_0152;
    const ccomplex_t IT_0218 = m_d*IT_0159;
    const ccomplex_t IT_0219 = m_d*IT_0166;
    const ccomplex_t IT_0220 = m_d*IT_0173;
    const ccomplex_t IT_0221 = -IT_0207 + -IT_0210 + -IT_0213 + -IT_0216 +
       IT_0217 + IT_0218 + IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = IT_0078*IT_0149;
    const ccomplex_t IT_0223 = IT_0080*IT_0205*IT_0222;
    const ccomplex_t IT_0224 = IT_0098*IT_0156;
    const ccomplex_t IT_0225 = IT_0100*IT_0208*IT_0224;
    const ccomplex_t IT_0226 = IT_0118*IT_0163;
    const ccomplex_t IT_0227 = IT_0120*IT_0211*IT_0226;
    const ccomplex_t IT_0228 = IT_0138*IT_0170;
    const ccomplex_t IT_0229 = IT_0140*IT_0214*IT_0228;
    const ccomplex_t IT_0230 = m_d*IT_0082;
    const ccomplex_t IT_0231 = m_d*IT_0102;
    const ccomplex_t IT_0232 = m_d*IT_0122;
    const ccomplex_t IT_0233 = m_d*IT_0142;
    const ccomplex_t IT_0234 = (-0.166666666666667)*IT_0057 + -IT_0223 + 
      -IT_0225 + -IT_0227 + -IT_0229 + IT_0230 + IT_0231 + IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = (-18)*IT_0234;
    const ccomplex_t IT_0236 = (-6)*IT_0058;
    const ccomplex_t IT_0237 = IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = (-18)*conj(IT_0234);
    const ccomplex_t IT_0239 = (-6)*conj(IT_0058);
    const ccomplex_t IT_0240 = IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = (-18)*IT_0058;
    const ccomplex_t IT_0242 = (-6)*IT_0234;
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = (-18)*conj(IT_0058);
    const ccomplex_t IT_0245 = (-6)*conj(IT_0234);
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = (-18)*IT_0143 + (-6)*IT_0145;
    const ccomplex_t IT_0248 = (-18)*IT_0145;
    const ccomplex_t IT_0249 = (-6)*IT_0143;
    const ccomplex_t IT_0250 = IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = (-18)*conj(IT_0145);
    const ccomplex_t IT_0252 = (-6)*conj(IT_0143);
    const ccomplex_t IT_0253 = IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 6*s_34;
    const ccomplex_t IT_0255 = -IT_0008;
    const ccomplex_t IT_0256 = s_34 + IT_0255;
    const ccomplex_t IT_0257 = 6*IT_0027;
    const ccomplex_t IT_0258 = 6*conj(IT_0027);
    const ccomplex_t IT_0259 = m_d*s_14;
    const ccomplex_t IT_0260 = (-18)*IT_0221;
    const ccomplex_t IT_0261 = (-6)*IT_0221;
    const ccomplex_t IT_0262 = (-6)*IT_0200;
    const ccomplex_t IT_0263 = IT_0178 + IT_0180 + IT_0262;
    const ccomplex_t IT_0264 = conj(IT_0143) + conj(IT_0174);
    const ccomplex_t IT_0265 = conj(IT_0145) + conj(IT_0176);
    const ccomplex_t IT_0266 = IT_0145 + IT_0176;
    const ccomplex_t IT_0267 = -IT_0204;
    const ccomplex_t IT_0268 = IT_0143 + IT_0174;
    const ccomplex_t IT_0269 = -IT_0197;
    const ccomplex_t IT_0270 = (-18)*IT_0200;
    const ccomplex_t IT_0271 = IT_0182 + IT_0183 + IT_0270;
    const ccomplex_t IT_0272 = m_d*s_13;
    const ccomplex_t IT_0273 = 6*conj(IT_0058);
    const ccomplex_t IT_0274 = 18*conj(IT_0234);
    const ccomplex_t IT_0275 = 6*conj(IT_0234);
    const ccomplex_t IT_0276 = 18*conj(IT_0058);
    const ccomplex_t IT_0277 = 18*conj(IT_0221);
    const ccomplex_t IT_0278 = 6*conj(IT_0221);
    const ccomplex_t IT_0279 = 18*conj(IT_0059);
    const ccomplex_t IT_0280 = IT_0058*(conj(IT_0143) + 3*conj(IT_0145)) +
       IT_0059*(conj(IT_0174) + 3*conj(IT_0176)) + 0.166666666666667*(18*conj
      (IT_0174) + 6*conj(IT_0176))*IT_0221 + 0.166666666666667*(18*conj(IT_0143)
       + 6*conj(IT_0145))*IT_0234 + 0.166666666666667*IT_0143*(IT_0273 + IT_0274
      ) + 0.166666666666667*IT_0145*(IT_0275 + IT_0276) + IT_0174*(conj(IT_0059)
       + 0.166666666666667*IT_0277) + 0.166666666666667*IT_0176*(IT_0278 +
       IT_0279);
    const ccomplex_t IT_0281 = -s_14;
    const ccomplex_t IT_0282 = s_13 + IT_0281;
    const ccomplex_t IT_0283 = m_d*IT_0282;
    const ccomplex_t IT_0284 = 18*IT_0027;
    const ccomplex_t IT_0285 = 18*conj(IT_0027);
    const ccomplex_t IT_0286 = IT_0257*IT_0264 + IT_0258*IT_0268 + IT_0204*((
      -12)*IT_0027 + 12*conj(IT_0027) + (-6)*IT_0059 + 6*conj(IT_0059) + IT_0235
       + IT_0236 + IT_0260 + IT_0273 + IT_0274 + IT_0277) + IT_0197*((-36)
      *IT_0027 + 36*conj(IT_0027) + (-18)*IT_0059 + IT_0241 + IT_0242 + IT_0261 
      + IT_0275 + IT_0276 + IT_0278 + IT_0279) + IT_0265*IT_0284 + IT_0266
      *IT_0285;
    const ccomplex_t IT_0287 = (36*s_34 + (-36)*IT_0008)*IT_0027*conj(IT_0027)
       + (conj(IT_0027)*(IT_0058 + IT_0059) + IT_0027*(conj(IT_0058) + conj
      (IT_0059)))*IT_0062 + (conj(IT_0143)*IT_0145 + IT_0143*conj(IT_0145) +
       conj(IT_0174)*IT_0176 + IT_0174*conj(IT_0176))*IT_0181 + (IT_0143*conj
      (IT_0143) + IT_0145*conj(IT_0145) + IT_0174*conj(IT_0174) + IT_0176*conj
      (IT_0176))*IT_0184 + -IT_0198*IT_0203 + (-2)*IT_0204*(IT_0197*((-12)
      *IT_0177 + 24*IT_0179 + (-12)*IT_0200) + 0.5*IT_0203*IT_0204) + IT_0061*
      (IT_0058*conj(IT_0058) + IT_0059*conj(IT_0059) + IT_0221*conj(IT_0221) +
       IT_0234*conj(IT_0234)) + IT_0008*(conj(IT_0221)*IT_0237 + IT_0221*IT_0240
       + conj(IT_0059)*IT_0243 + IT_0059*IT_0246) + IT_0200*(((-18)*conj(IT_0143
      ) + (-6)*conj(IT_0145))*IT_0174 + conj(IT_0174)*IT_0247 + conj(IT_0176)
      *IT_0250 + IT_0176*IT_0253) + (conj(IT_0059)*IT_0221 + IT_0059*conj
      (IT_0221) + conj(IT_0058)*IT_0234 + IT_0058*conj(IT_0234))*IT_0254 +
       IT_0256*((conj(IT_0221) + conj(IT_0234))*IT_0257 + (IT_0221 + IT_0234)
      *IT_0258) + IT_0259*(conj(IT_0174)*IT_0237 + IT_0174*IT_0240 + conj
      (IT_0176)*IT_0243 + IT_0176*IT_0246 + conj(IT_0221)*IT_0247 + conj(IT_0059
      )*IT_0250 + IT_0059*IT_0253 + conj(IT_0143)*IT_0260 + conj(IT_0145)
      *IT_0261) + IT_0263*(IT_0197*IT_0264 + IT_0204*IT_0265 + IT_0266*IT_0267 +
       IT_0268*IT_0269) + (IT_0204*IT_0264 + IT_0197*IT_0265 + IT_0267*IT_0268 +
       IT_0266*IT_0269)*IT_0271 + 6*IT_0272*IT_0280 + IT_0283*IT_0286;
    return create_ccomplex_return(IT_0287);
}

