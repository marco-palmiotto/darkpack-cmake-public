#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_1_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_1_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*IT_0000*IT_0002*U_st_00*conj
      (U_st_00);
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = pow(m_st_1, 2);
    const ccomplex_t IT_0006 = cpow(s_12 + IT_0005 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0009 = (-2)*IT_0008;
    const ccomplex_t IT_0010 = (-0.666666666666667)*IT_0008;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = s_13 + s_23;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = -IT_0008;
    const ccomplex_t IT_0015 = (-1.33333333333333)*IT_0008;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = s_12 + IT_0005;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = IT_0009*IT_0015;
    const ccomplex_t IT_0020 = s_23*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0013 + (-0.5)*IT_0018 + 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0010*IT_0014;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = IT_0021 + IT_0023;
    const ccomplex_t IT_0025 = IT_0007*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0005 + -reg_prop, -1);
    const ccomplex_t IT_0027 = cos(alpha);
    const ccomplex_t IT_0028 = sin(beta);
    const ccomplex_t IT_0029 = cos(beta);
    const ccomplex_t IT_0030 = sin(alpha);
    const ccomplex_t IT_0031 = cpow(IT_0001, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_W*e_em*(IT_0027*IT_0028 + 
      -IT_0029*IT_0030)*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = m_t*mu_h*IT_0033;
    const ccomplex_t IT_0035 = cpow(IT_0028, -1);
    const ccomplex_t IT_0036 = IT_0029*IT_0030;
    const ccomplex_t IT_0037 = IT_0027*IT_0028;
    const ccomplex_t IT_0038 = IT_0036 + IT_0037;
    const ccomplex_t IT_0039 = cos(theta_W);
    const ccomplex_t IT_0040 = cpow(IT_0039, -2);
    const ccomplex_t IT_0041 = pow(m_t, 2);
    const ccomplex_t IT_0042 = 2*U_st_00*(((-0.5)*m_W*IT_0038*(IT_0031 + -1./3
      *IT_0001*IT_0040) + IT_0027*IT_0031*IT_0033*IT_0035*IT_0041)*conj(U_st_00)
       + 1./2*IT_0030*IT_0031*IT_0034*IT_0035*conj(U_st_10)) + U_st_10*(IT_0030
      *IT_0031*IT_0034*IT_0035*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0001
      *IT_0038*IT_0040 + (-1.5)*IT_0027*IT_0031*IT_0033*IT_0035*IT_0041)*conj
      (U_st_10));
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(e_em*IT_0042 +
       1.4142135623731*A_t*IT_0027*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0032*IT_0044;
    const ccomplex_t IT_0046 = IT_0026*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + (-2)*IT_0005 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_W*e_em*(IT_0027*IT_0029 +
       IT_0028*IT_0030)*IT_0031;
    const ccomplex_t IT_0050 = IT_0027*IT_0029;
    const ccomplex_t IT_0051 = -IT_0028*IT_0030 + IT_0050;
    const ccomplex_t IT_0052 = (-2)*U_st_00*((IT_0030*IT_0031*IT_0033*IT_0035
      *IT_0041 + 0.5*m_W*(IT_0031 + -1./3*IT_0001*IT_0040)*IT_0051)*conj(U_st_00
      ) + -1./2*IT_0027*IT_0031*IT_0034*IT_0035*conj(U_st_10)) + U_st_10*
      (IT_0027*IT_0031*IT_0034*IT_0035*conj(U_st_00) + (-1.33333333333333)*(1.5
      *IT_0030*IT_0031*IT_0033*IT_0035*IT_0041 + m_W*IT_0001*IT_0040*IT_0051)
      *conj(U_st_10));
    const ccomplex_t IT_0053 = (0 + _Complex_I*-1)*(e_em*IT_0052 + (
      -1.4142135623731)*A_t*IT_0030*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0049*IT_0054;
    const ccomplex_t IT_0056 = IT_0048*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*e_em*IT_0031*IT_0039;
    const ccomplex_t IT_0059 = (-2)*IT_0058;
    const ccomplex_t IT_0060 = cpow(IT_0039, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0031
      *IT_0039 + -1./3*IT_0001*IT_0060)*U_st_00*conj(U_st_00) + IT_0001*IT_0060
      *U_st_10*conj(U_st_10));
    const ccomplex_t IT_0062 = 0.666666666666667*IT_0061;
    const ccomplex_t IT_0063 = IT_0059*IT_0062;
    const ccomplex_t IT_0064 = IT_0012*IT_0063;
    const ccomplex_t IT_0065 = -IT_0058;
    const ccomplex_t IT_0066 = 1.33333333333333*IT_0061;
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = IT_0017*IT_0067;
    const ccomplex_t IT_0069 = IT_0064 + IT_0068;
    const ccomplex_t IT_0070 = IT_0062*IT_0065;
    const ccomplex_t IT_0071 = IT_0017*IT_0070;
    const ccomplex_t IT_0072 = IT_0059*IT_0066;
    const ccomplex_t IT_0073 = s_23*IT_0072;
    const ccomplex_t IT_0074 = (-2)*IT_0071 + -IT_0073;
    const ccomplex_t IT_0075 = IT_0069 + IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0005 + -reg_prop, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = IT_0075*IT_0077;
    const ccomplex_t IT_0079 = IT_0004 + IT_0025 + -IT_0047 + -IT_0057 +
       IT_0078;
    const ccomplex_t IT_0080 = pow(m_W, -4);
    const ccomplex_t IT_0081 = pow(s_34, 2);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = 2*IT_0008;
    const ccomplex_t IT_0084 = IT_0010*IT_0083;
    const ccomplex_t IT_0085 = IT_0022 + IT_0084;
    const ccomplex_t IT_0086 = IT_0007*IT_0085;
    const ccomplex_t IT_0087 = 2*IT_0058;
    const ccomplex_t IT_0088 = IT_0062*IT_0087;
    const ccomplex_t IT_0089 = IT_0070 + IT_0088;
    const ccomplex_t IT_0090 = IT_0077*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0086 + -IT_0090;
    const ccomplex_t IT_0092 = pow(m_W, -2);
    const ccomplex_t IT_0093 = pow(s_13, 2);
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = -IT_0005;
    const ccomplex_t IT_0096 = pow(s_14, 2);
    const ccomplex_t IT_0097 = IT_0092*IT_0096;
    const ccomplex_t IT_0098 = (IT_0094 + IT_0095)*(IT_0095 + IT_0097);
    const ccomplex_t IT_0099 = 3*IT_0091;
    const ccomplex_t IT_0100 = s_13*s_14*s_34*IT_0080;
    const ccomplex_t IT_0101 = IT_0094 + IT_0095 + IT_0097 + -IT_0100;
    const ccomplex_t IT_0102 = (-3)*IT_0079;
    const ccomplex_t IT_0103 = (-3)*conj(IT_0079);
    const ccomplex_t IT_0104 = IT_0015*IT_0083;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = IT_0085 + IT_0105;
    const ccomplex_t IT_0107 = IT_0007*IT_0106;
    const ccomplex_t IT_0108 = IT_0066*IT_0087;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = IT_0089 + IT_0109;
    const ccomplex_t IT_0111 = IT_0077*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0107 + -IT_0111;
    const ccomplex_t IT_0113 = s_13*s_24*s_34*IT_0080;
    const ccomplex_t IT_0114 = s_14*s_24*IT_0092;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = s_13*s_23*IT_0092;
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = s_12 + IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = IT_0113 + IT_0118;
    const ccomplex_t IT_0120 = 3*IT_0079;
    const ccomplex_t IT_0121 = pow(s_24, 2);
    const ccomplex_t IT_0122 = IT_0092*IT_0121;
    const ccomplex_t IT_0123 = -(IT_0005 + -IT_0094)*(IT_0005 + -IT_0122);
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = 3*IT_0112;
    const ccomplex_t IT_0126 = -(IT_0005 + -IT_0094)*(s_12 + -IT_0114);
    const ccomplex_t IT_0127 = (-3)*IT_0091;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0031*conj(U_st_00);
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0031*U_st_00;
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = pow(m_W, 2);
    const ccomplex_t IT_0134 = cpow((-2)*s_23 + IT_0005 + IT_0133 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0132*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0031*conj(U_st_00);
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0031*U_st_00;
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = cpow((-2)*s_23 + IT_0005 + IT_0133 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0141*IT_0143;
    const ccomplex_t IT_0145 = -IT_0016;
    const ccomplex_t IT_0146 = IT_0085 + IT_0145;
    const ccomplex_t IT_0147 = IT_0007*IT_0146;
    const ccomplex_t IT_0148 = -IT_0067;
    const ccomplex_t IT_0149 = IT_0089 + IT_0148;
    const ccomplex_t IT_0150 = IT_0077*IT_0149;
    const ccomplex_t IT_0151 = cpow((-2)*s_23 + IT_0005 + IT_0133 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0031*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0031*U_sb_00*conj(U_st_00);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0154*IT_0156;
    const ccomplex_t IT_0158 = IT_0152*IT_0157;
    const ccomplex_t IT_0159 = cpow((-2)*s_23 + IT_0005 + IT_0133 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0031*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0031*U_sb_01*conj(U_st_00);
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = IT_0162*IT_0164;
    const ccomplex_t IT_0166 = IT_0160*IT_0165;
    const ccomplex_t IT_0167 = -IT_0136 + -IT_0144 + 0.5*IT_0147 + -IT_0150 + 
      -IT_0158 + -IT_0166;
    const ccomplex_t IT_0168 = s_14*s_23*s_34*IT_0080;
    const ccomplex_t IT_0169 = IT_0118 + IT_0168;
    const ccomplex_t IT_0170 = (s_12 + -IT_0114)*(s_12 + -IT_0116);
    const ccomplex_t IT_0171 = s_14*s_34*IT_0092;
    const ccomplex_t IT_0172 = (s_12 + -IT_0116)*(s_13 + -IT_0171);
    const ccomplex_t IT_0173 = IT_0008*IT_0010;
    const ccomplex_t IT_0174 = IT_0007*IT_0173;
    const ccomplex_t IT_0175 = IT_0058*IT_0062;
    const ccomplex_t IT_0176 = IT_0077*IT_0175;
    const ccomplex_t IT_0177 = (-0.5)*IT_0174 + IT_0176;
    const ccomplex_t IT_0178 = 3*IT_0177;
    const ccomplex_t IT_0179 = -(IT_0005 + -IT_0097)*(s_12 + -IT_0116);
    const ccomplex_t IT_0180 = 3*conj(IT_0079);
    const ccomplex_t IT_0181 = (-3)*conj(IT_0091);
    const ccomplex_t IT_0182 = s_34*IT_0092;
    const ccomplex_t IT_0183 = s_13*s_34*IT_0092;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = s_14 + IT_0184;
    const ccomplex_t IT_0186 = IT_0182*IT_0185;
    const ccomplex_t IT_0187 = -(IT_0005 + -IT_0094)*(s_13 + -IT_0171);
    const ccomplex_t IT_0188 = s_24*s_34*IT_0092;
    const ccomplex_t IT_0189 = -(IT_0005 + -IT_0094)*(s_23 + -IT_0188);
    const ccomplex_t IT_0190 = (-3)*IT_0112;
    const ccomplex_t IT_0191 = IT_0081*IT_0092;
    const ccomplex_t IT_0192 = -(IT_0005 + -IT_0094)*(IT_0133 + -IT_0191);
    const ccomplex_t IT_0193 = (-3)*conj(IT_0112);
    const ccomplex_t IT_0194 = -IT_0130;
    const ccomplex_t IT_0195 = IT_0129*IT_0194;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0132 + IT_0196;
    const ccomplex_t IT_0198 = IT_0135*IT_0197;
    const ccomplex_t IT_0199 = -IT_0139;
    const ccomplex_t IT_0200 = IT_0138*IT_0199;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0141 + IT_0201;
    const ccomplex_t IT_0203 = IT_0143*IT_0202;
    const ccomplex_t IT_0204 = -IT_0016 + -IT_0104;
    const ccomplex_t IT_0205 = IT_0085 + IT_0204;
    const ccomplex_t IT_0206 = IT_0007*IT_0205;
    const ccomplex_t IT_0207 = -IT_0153;
    const ccomplex_t IT_0208 = IT_0156*IT_0207;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = IT_0157 + IT_0209;
    const ccomplex_t IT_0211 = IT_0152*IT_0210;
    const ccomplex_t IT_0212 = -IT_0161;
    const ccomplex_t IT_0213 = IT_0164*IT_0212;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = IT_0165 + IT_0214;
    const ccomplex_t IT_0216 = IT_0160*IT_0215;
    const ccomplex_t IT_0217 = -IT_0067 + -IT_0108;
    const ccomplex_t IT_0218 = IT_0089 + IT_0217;
    const ccomplex_t IT_0219 = IT_0077*IT_0218;
    const ccomplex_t IT_0220 = -IT_0198 + -IT_0203 + 0.5*IT_0206 + -IT_0211 + 
      -IT_0216 + -IT_0219;
    const ccomplex_t IT_0221 = s_23*s_24*s_34*IT_0080;
    const ccomplex_t IT_0222 = -IT_0122;
    const ccomplex_t IT_0223 = pow(s_23, 2);
    const ccomplex_t IT_0224 = IT_0092*IT_0223;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = IT_0005 + IT_0221 + IT_0222 + IT_0225;
    const ccomplex_t IT_0227 = (s_12 + -IT_0116)*(s_23 + -IT_0188);
    const ccomplex_t IT_0228 = (s_12 + -IT_0116)*(IT_0005 + -IT_0122);
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = (s_12 + -IT_0114)*(IT_0005 + -IT_0224);
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = (-3)*IT_0167;
    const ccomplex_t IT_0233 = IT_0008*IT_0015;
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0173 + IT_0234;
    const ccomplex_t IT_0236 = IT_0007*IT_0235;
    const ccomplex_t IT_0237 = IT_0058*IT_0066;
    const ccomplex_t IT_0238 = -IT_0237;
    const ccomplex_t IT_0239 = IT_0175 + IT_0238;
    const ccomplex_t IT_0240 = IT_0077*IT_0239;
    const ccomplex_t IT_0241 = IT_0198 + IT_0203 + IT_0211 + IT_0216 + (-0.5)
      *IT_0236 + IT_0240;
    const ccomplex_t IT_0242 = (s_12 + -IT_0116)*(IT_0133 + -IT_0191);
    const ccomplex_t IT_0243 = s_23*s_34*IT_0092;
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = s_24 + IT_0244;
    const ccomplex_t IT_0246 = IT_0182*IT_0245;
    const ccomplex_t IT_0247 = (s_13 + -IT_0171)*(IT_0005 + -IT_0224);
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = (s_23 + -IT_0188)*(IT_0005 + -IT_0224);
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = (IT_0133 + -IT_0191)*(IT_0005 + -IT_0224);
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = 3*conj(IT_0091);
    const ccomplex_t IT_0254 = 3*conj(IT_0112);
    const ccomplex_t IT_0255 = 3*conj(IT_0177);
    const ccomplex_t IT_0256 = (-3)*conj(IT_0167);
    const ccomplex_t IT_0257 = (-3)*IT_0241*((-0.333333333333333)*IT_0103
      *IT_0246 + conj(IT_0220)*IT_0250 + conj(IT_0241)*IT_0252 + (
      -0.333333333333333)*IT_0172*IT_0253 + (-0.333333333333333)*IT_0227*IT_0254
       + (-0.333333333333333)*IT_0242*IT_0255 + (-0.333333333333333)*IT_0248
      *IT_0256);
    const ccomplex_t IT_0258 = (IT_0005 + -IT_0122)*(IT_0005 + -IT_0224);
    const ccomplex_t IT_0259 = 0.333333333333333*IT_0170;
    const ccomplex_t IT_0260 = 3*IT_0220*(0.333333333333333*IT_0180*IT_0226 +
       0.333333333333333*IT_0193*IT_0229 + -conj(IT_0241)*IT_0250 +
       0.333333333333333*IT_0227*IT_0255 + 0.333333333333333*IT_0231*IT_0256 +
       conj(IT_0220)*IT_0258 + IT_0253*IT_0259);
    const ccomplex_t IT_0261 = -(IT_0005 + -IT_0097)*(IT_0005 + -IT_0224);
    const ccomplex_t IT_0262 = -IT_0261;
    const ccomplex_t IT_0263 = 6*IT_0079*conj(IT_0079)*(1 + 0.5*IT_0082) +
       conj(IT_0091)*(IT_0098*IT_0099 + IT_0101*IT_0102) + IT_0091*IT_0101
      *IT_0103 + conj(IT_0112)*(IT_0119*IT_0120 + IT_0124*IT_0125 + IT_0126
      *IT_0127) + conj(IT_0167)*(IT_0120*IT_0169 + IT_0125*IT_0170 + IT_0172
      *IT_0178 + IT_0127*IT_0179) + IT_0112*(IT_0119*IT_0180 + IT_0126*IT_0181) 
      + conj(IT_0177)*(IT_0102*IT_0186 + IT_0127*IT_0187 + IT_0189*IT_0190) + (
      -3)*IT_0177*((-0.333333333333333)*IT_0103*IT_0186 + (-0.333333333333333)
      *IT_0181*IT_0187 + conj(IT_0177)*IT_0192 + (-0.333333333333333)*IT_0189
      *IT_0193) + conj(IT_0220)*(IT_0099*IT_0170 + IT_0120*IT_0226 + IT_0178
      *IT_0227 + IT_0190*IT_0229 + IT_0231*IT_0232) + conj(IT_0241)*(IT_0099
      *IT_0172 + IT_0125*IT_0227 + IT_0178*IT_0242 + IT_0102*IT_0246 + IT_0232
      *IT_0248) + IT_0257 + IT_0260 + 3*IT_0167*(0.333333333333333*IT_0169
      *IT_0180 + 0.333333333333333*IT_0179*IT_0181 + 0.333333333333333*IT_0172
      *IT_0255 + IT_0254*IT_0259 + conj(IT_0167)*IT_0262);
    return create_ccomplex_return(IT_0263);
}

