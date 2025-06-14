#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_se_L_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_se_L_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
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
    const creal_t m_se_L = param->m_se_L;
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
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_e, 2);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0003*IT_0005*IT_0008 + (-0.5)*m_W*(IT_0001*IT_0003 + -IT_0009*IT_0010)
      *(IT_0005 + -IT_0004*IT_0012));
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_se_L, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0010;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0005*IT_0008*IT_0010 + (-0.5)*m_W*(IT_0003*IT_0009 + IT_0001*IT_0010)*
      (IT_0005 + -IT_0004*IT_0012));
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0011, -1);
    const ccomplex_t IT_0029 = N_B1*e_em;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = N_W1*e_em;
    const ccomplex_t IT_0033 = IT_0005*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + IT_0034);
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*m_e*N_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0040*IT_0041;
    const ccomplex_t IT_0043 = N_B2*e_em;
    const ccomplex_t IT_0044 = IT_0028*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = N_W2*e_em;
    const ccomplex_t IT_0047 = IT_0005*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*m_e*N_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0053*IT_0054*IT_0055;
    const ccomplex_t IT_0057 = N_B3*e_em;
    const ccomplex_t IT_0058 = IT_0028*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = N_W3*e_em;
    const ccomplex_t IT_0061 = IT_0005*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + IT_0062);
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0067*IT_0068*IT_0069;
    const ccomplex_t IT_0071 = N_B4*e_em;
    const ccomplex_t IT_0072 = IT_0028*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = N_W4*e_em;
    const ccomplex_t IT_0075 = IT_0005*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + IT_0076);
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0083 = cpow((-2)*s_23 + IT_0008 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0081*IT_0082*IT_0083;
    const ccomplex_t IT_0085 = IT_0004*IT_0028;
    const ccomplex_t IT_0086 = e_em*IT_0085;
    const ccomplex_t IT_0087 = IT_0005*IT_0011;
    const ccomplex_t IT_0088 = e_em*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0086 + -IT_0088);
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*-1)*e_em*(IT_0005*IT_0011 + 
      -IT_0004*IT_0028);
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = m_e*IT_0095;
    const ccomplex_t IT_0097 = -IT_0094;
    const ccomplex_t IT_0098 = m_e*IT_0097;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0096 + IT_0099;
    const ccomplex_t IT_0101 = IT_0093*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0028;
    const ccomplex_t IT_0103 = IT_0091*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0100*IT_0104;
    const ccomplex_t IT_0106 = conj(N_B1)*e_em;
    const ccomplex_t IT_0107 = IT_0028*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = conj(N_W1)*e_em;
    const ccomplex_t IT_0110 = IT_0005*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + IT_0111);
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0036*IT_0113;
    const ccomplex_t IT_0115 = IT_0041*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = m_e*IT_0116;
    const ccomplex_t IT_0118 = conj(N_B2)*e_em;
    const ccomplex_t IT_0119 = IT_0028*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = conj(N_W2)*e_em;
    const ccomplex_t IT_0122 = IT_0005*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + IT_0123);
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0050*IT_0125;
    const ccomplex_t IT_0127 = IT_0055*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = m_e*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B3)*e_em;
    const ccomplex_t IT_0131 = IT_0028*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = conj(N_W3)*e_em;
    const ccomplex_t IT_0134 = IT_0005*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + IT_0135);
    const ccomplex_t IT_0137 = (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = IT_0064*IT_0137;
    const ccomplex_t IT_0139 = IT_0069*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_e*IT_0140;
    const ccomplex_t IT_0142 = conj(N_B4)*e_em;
    const ccomplex_t IT_0143 = IT_0028*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = conj(N_W4)*e_em;
    const ccomplex_t IT_0146 = IT_0005*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0144 + IT_0147);
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = IT_0078*IT_0149;
    const ccomplex_t IT_0151 = IT_0083*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = m_e*IT_0152;
    const ccomplex_t IT_0154 = -IT_0042 + -IT_0056 + -IT_0070 + -IT_0084 + 
      -IT_0101 + IT_0105 + IT_0117 + IT_0129 + IT_0141 + IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = IT_0113*IT_0156;
    const ccomplex_t IT_0158 = IT_0040*IT_0041*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = IT_0125*IT_0160;
    const ccomplex_t IT_0162 = IT_0054*IT_0055*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = IT_0137*IT_0164;
    const ccomplex_t IT_0166 = IT_0068*IT_0069*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0168 = (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = IT_0149*IT_0168;
    const ccomplex_t IT_0170 = IT_0082*IT_0083*IT_0169;
    const ccomplex_t IT_0171 = IT_0038*IT_0156;
    const ccomplex_t IT_0172 = IT_0041*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = m_e*IT_0173;
    const ccomplex_t IT_0175 = IT_0052*IT_0160;
    const ccomplex_t IT_0176 = IT_0055*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = m_e*IT_0177;
    const ccomplex_t IT_0179 = IT_0066*IT_0164;
    const ccomplex_t IT_0180 = IT_0069*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = m_e*IT_0181;
    const ccomplex_t IT_0183 = IT_0080*IT_0168;
    const ccomplex_t IT_0184 = IT_0083*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = m_e*IT_0185;
    const ccomplex_t IT_0187 = IT_0101 + -IT_0105 + -IT_0158 + -IT_0162 + 
      -IT_0166 + -IT_0170 + IT_0174 + IT_0178 + IT_0182 + IT_0186;
    const ccomplex_t IT_0188 = -IT_0008;
    const ccomplex_t IT_0189 = s_34 + IT_0188;
    const ccomplex_t IT_0190 = IT_0093*IT_0097;
    const ccomplex_t IT_0191 = -IT_0116 + -IT_0128 + -IT_0140 + -IT_0152 + 
      -IT_0190;
    const ccomplex_t IT_0192 = m_e*s_13;
    const ccomplex_t IT_0193 = m_e*s_14;
    const ccomplex_t IT_0194 = IT_0097*IT_0104;
    const ccomplex_t IT_0195 = -IT_0173 + -IT_0177 + -IT_0181 + -IT_0185 + 
      -IT_0194;
    const ccomplex_t IT_0196 = IT_0008*IT_0016;
    const ccomplex_t IT_0197 = s_13*s_14;
    const ccomplex_t IT_0198 = s_34*IT_0016;
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0197 + IT_0199;
    const ccomplex_t IT_0201 = 2*s_34;
    const ccomplex_t IT_0202 = (-0.5)*IT_0196;
    const ccomplex_t IT_0203 = IT_0197 + IT_0199 + IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = 2*IT_0204;
    const ccomplex_t IT_0207 = IT_0205*IT_0206;
    const ccomplex_t IT_0208 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0209 = IT_0207*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = 0.5*IT_0210;
    const ccomplex_t IT_0212 = cpow(IT_0211, 2);
    const ccomplex_t IT_0213 = conj(IT_0191) + conj(IT_0195);
    const ccomplex_t IT_0214 = 4*IT_0211;
    const ccomplex_t IT_0215 = IT_0191 + IT_0195;
    const ccomplex_t IT_0216 = (-4)*IT_0211;
    const ccomplex_t IT_0217 = (-2)*IT_0154;
    const ccomplex_t IT_0218 = (-2)*conj(IT_0154);
    const ccomplex_t IT_0219 = (-2)*IT_0187;
    const ccomplex_t IT_0220 = -s_14;
    const ccomplex_t IT_0221 = s_13 + IT_0220;
    const ccomplex_t IT_0222 = m_e*IT_0221;
    const ccomplex_t IT_0223 = 2*(conj(IT_0027)*(IT_0154 + IT_0187) + 2
      *IT_0027*(conj(IT_0027) + 0.5*conj(IT_0154) + 0.5*conj(IT_0187)))*IT_0189 
      + 2*IT_0191*(conj(IT_0154)*IT_0192 + -conj(IT_0187)*IT_0193 + -conj
      (IT_0195)*IT_0196 + 2*conj(IT_0191)*IT_0200) + 2*IT_0154*(conj(IT_0191)
      *IT_0192 + 0.5*conj(IT_0154)*IT_0201) + (-8)*IT_0203*(IT_0212 + (-0.125)
      *IT_0213*IT_0214 + (-0.125)*IT_0215*IT_0216) + (IT_0008*conj(IT_0187) +
       IT_0193*conj(IT_0195))*IT_0217 + 2*IT_0187*(IT_0192*conj(IT_0195) + 0.5
      *conj(IT_0187)*IT_0201 + 0.5*IT_0008*IT_0218) + 2*IT_0195*(conj(IT_0187)
      *IT_0192 + -conj(IT_0191)*IT_0196 + 2*conj(IT_0195)*IT_0200 + 0.5*IT_0193
      *IT_0218) + conj(IT_0191)*IT_0193*IT_0219 + 2*(IT_0027*IT_0213 + conj
      (IT_0027)*IT_0215 + (-2)*IT_0211*(IT_0027 + -conj(IT_0027) + (-0.5)*conj
      (IT_0154) + (-0.5)*conj(IT_0187) + (-0.25)*IT_0217 + (-0.25)*IT_0219))
      *IT_0222;
    return create_ccomplex_return(IT_0223);
}

