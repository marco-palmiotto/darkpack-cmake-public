#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_anti_W_W(
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
    const creal_t m_st_2 = param->m_st_2;
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
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*IT_0000*IT_0002*U_st_01*conj
      (U_st_01);
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = cpow(IT_0001, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*m_W*e_em*(IT_0005*IT_0006 + 
      -IT_0007*IT_0008)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = m_t*mu_h*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0006, -1);
    const ccomplex_t IT_0014 = IT_0007*IT_0008;
    const ccomplex_t IT_0015 = IT_0005*IT_0006;
    const ccomplex_t IT_0016 = IT_0014 + IT_0015;
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -2);
    const ccomplex_t IT_0019 = pow(m_t, 2);
    const ccomplex_t IT_0020 = 2*U_st_01*(((-0.5)*m_W*IT_0016*(IT_0009 + -1./3
      *IT_0001*IT_0018) + IT_0005*IT_0009*IT_0011*IT_0013*IT_0019)*conj(U_st_01)
       + 1./2*IT_0008*IT_0009*IT_0012*IT_0013*conj(U_st_11)) + U_st_11*(IT_0008
      *IT_0009*IT_0012*IT_0013*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0001
      *IT_0016*IT_0018 + (-1.5)*IT_0005*IT_0009*IT_0011*IT_0013*IT_0019)*conj
      (U_st_11));
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(e_em*IT_0020 +
       1.4142135623731*A_t*IT_0005*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0010*IT_0022;
    const ccomplex_t IT_0024 = pow(m_st_2, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0024 + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_W*e_em*(IT_0005*IT_0007 +
       IT_0006*IT_0008)*IT_0009;
    const ccomplex_t IT_0029 = IT_0005*IT_0007;
    const ccomplex_t IT_0030 = -IT_0006*IT_0008 + IT_0029;
    const ccomplex_t IT_0031 = (-2)*U_st_01*((IT_0008*IT_0009*IT_0011*IT_0013
      *IT_0019 + 0.5*m_W*(IT_0009 + -1./3*IT_0001*IT_0018)*IT_0030)*conj(U_st_01
      ) + -1./2*IT_0005*IT_0009*IT_0012*IT_0013*conj(U_st_11)) + U_st_11*
      (IT_0005*IT_0009*IT_0012*IT_0013*conj(U_st_01) + (-1.33333333333333)*(1.5
      *IT_0008*IT_0009*IT_0011*IT_0013*IT_0019 + m_W*IT_0001*IT_0018*IT_0030)
      *conj(U_st_11));
    const ccomplex_t IT_0032 = (0 + _Complex_I*-1)*(e_em*IT_0031 + (
      -1.4142135623731)*A_t*IT_0008*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + (-2)*IT_0024 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow(s_12 + IT_0024 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = (-0.666666666666667)*IT_0040;
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = s_12 + IT_0024;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (-2)*IT_0040;
    const ccomplex_t IT_0047 = IT_0042*IT_0046;
    const ccomplex_t IT_0048 = s_13 + s_23;
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (-1.33333333333333)*IT_0040;
    const ccomplex_t IT_0051 = IT_0041*IT_0050;
    const ccomplex_t IT_0052 = IT_0044*IT_0051;
    const ccomplex_t IT_0053 = IT_0046*IT_0050;
    const ccomplex_t IT_0054 = s_23*IT_0053;
    const ccomplex_t IT_0055 = (-0.5)*IT_0049 + (-0.5)*IT_0052 + 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0045 + IT_0055;
    const ccomplex_t IT_0057 = IT_0039*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0024 + -reg_prop, -1);
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0017;
    const ccomplex_t IT_0061 = (-2)*IT_0060;
    const ccomplex_t IT_0062 = cpow(IT_0017, -1);
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0009
      *IT_0017 + -1./3*IT_0001*IT_0062)*U_st_01*conj(U_st_01) + IT_0001*IT_0062
      *U_st_11*conj(U_st_11));
    const ccomplex_t IT_0064 = 0.666666666666667*IT_0063;
    const ccomplex_t IT_0065 = IT_0061*IT_0064;
    const ccomplex_t IT_0066 = IT_0048*IT_0065;
    const ccomplex_t IT_0067 = -IT_0060;
    const ccomplex_t IT_0068 = 1.33333333333333*IT_0063;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = IT_0044*IT_0069;
    const ccomplex_t IT_0071 = IT_0066 + IT_0070;
    const ccomplex_t IT_0072 = IT_0064*IT_0067;
    const ccomplex_t IT_0073 = IT_0044*IT_0072;
    const ccomplex_t IT_0074 = IT_0061*IT_0068;
    const ccomplex_t IT_0075 = s_23*IT_0074;
    const ccomplex_t IT_0076 = (-2)*IT_0073 + -IT_0075;
    const ccomplex_t IT_0077 = IT_0071 + IT_0076;
    const ccomplex_t IT_0078 = IT_0059*IT_0077;
    const ccomplex_t IT_0079 = IT_0004 + -IT_0027 + -IT_0037 + IT_0057 +
       IT_0078;
    const ccomplex_t IT_0080 = pow(m_W, -4);
    const ccomplex_t IT_0081 = pow(s_34, 2);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = 2*IT_0040;
    const ccomplex_t IT_0084 = IT_0042*IT_0083;
    const ccomplex_t IT_0085 = IT_0043 + IT_0084;
    const ccomplex_t IT_0086 = IT_0039*IT_0085;
    const ccomplex_t IT_0087 = 2*IT_0060;
    const ccomplex_t IT_0088 = IT_0064*IT_0087;
    const ccomplex_t IT_0089 = IT_0072 + IT_0088;
    const ccomplex_t IT_0090 = IT_0059*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0086 + -IT_0090;
    const ccomplex_t IT_0092 = pow(m_W, -2);
    const ccomplex_t IT_0093 = pow(s_13, 2);
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = -IT_0024;
    const ccomplex_t IT_0096 = pow(s_14, 2);
    const ccomplex_t IT_0097 = IT_0092*IT_0096;
    const ccomplex_t IT_0098 = (IT_0094 + IT_0095)*(IT_0095 + IT_0097);
    const ccomplex_t IT_0099 = 3*IT_0091;
    const ccomplex_t IT_0100 = s_13*s_14*s_34*IT_0080;
    const ccomplex_t IT_0101 = IT_0094 + IT_0095 + IT_0097 + -IT_0100;
    const ccomplex_t IT_0102 = (-3)*IT_0079;
    const ccomplex_t IT_0103 = (-3)*conj(IT_0079);
    const ccomplex_t IT_0104 = -IT_0051;
    const ccomplex_t IT_0105 = IT_0085 + IT_0104;
    const ccomplex_t IT_0106 = IT_0039*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0009*conj(U_st_01);
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0009*U_st_01;
    const ccomplex_t IT_0110 = (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = IT_0108*IT_0110;
    const ccomplex_t IT_0112 = pow(m_W, 2);
    const ccomplex_t IT_0113 = cpow((-2)*s_23 + IT_0024 + IT_0112 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0111*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0009*conj(U_st_01);
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0009*U_st_01;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = IT_0117*IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_23 + IT_0024 + IT_0112 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0009*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0009*U_sb_00*conj(U_st_01);
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = IT_0125*IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_23 + IT_0024 + IT_0112 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0009*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0009*U_sb_01*conj(U_st_01);
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = IT_0133*IT_0135;
    const ccomplex_t IT_0137 = cpow((-2)*s_23 + IT_0024 + IT_0112 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = IT_0136*IT_0138;
    const ccomplex_t IT_0140 = -IT_0069;
    const ccomplex_t IT_0141 = IT_0089 + IT_0140;
    const ccomplex_t IT_0142 = IT_0059*IT_0141;
    const ccomplex_t IT_0143 = 0.5*IT_0106 + -IT_0115 + -IT_0123 + -IT_0131 + 
      -IT_0139 + -IT_0142;
    const ccomplex_t IT_0144 = s_14*s_23*s_34*IT_0080;
    const ccomplex_t IT_0145 = s_14*s_24*IT_0092;
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = s_13*s_23*IT_0092;
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = s_12 + IT_0146 + IT_0148;
    const ccomplex_t IT_0150 = IT_0144 + IT_0149;
    const ccomplex_t IT_0151 = 3*IT_0079;
    const ccomplex_t IT_0152 = pow(s_23, 2);
    const ccomplex_t IT_0153 = IT_0092*IT_0152;
    const ccomplex_t IT_0154 = -(IT_0024 + -IT_0097)*(IT_0024 + -IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = 3*IT_0143;
    const ccomplex_t IT_0157 = -(IT_0024 + -IT_0097)*(s_12 + -IT_0147);
    const ccomplex_t IT_0158 = (-3)*IT_0091;
    const ccomplex_t IT_0159 = 3*conj(IT_0079);
    const ccomplex_t IT_0160 = (-3)*conj(IT_0091);
    const ccomplex_t IT_0161 = IT_0050*IT_0083;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0085 + IT_0162;
    const ccomplex_t IT_0164 = IT_0039*IT_0163;
    const ccomplex_t IT_0165 = IT_0068*IT_0087;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = IT_0089 + IT_0166;
    const ccomplex_t IT_0168 = IT_0059*IT_0167;
    const ccomplex_t IT_0169 = 0.5*IT_0164 + -IT_0168;
    const ccomplex_t IT_0170 = s_13*s_24*s_34*IT_0080;
    const ccomplex_t IT_0171 = IT_0149 + IT_0170;
    const ccomplex_t IT_0172 = (s_12 + -IT_0145)*(s_12 + -IT_0147);
    const ccomplex_t IT_0173 = 3*conj(IT_0143);
    const ccomplex_t IT_0174 = -(IT_0024 + -IT_0094)*(s_12 + -IT_0145);
    const ccomplex_t IT_0175 = -IT_0109;
    const ccomplex_t IT_0176 = IT_0108*IT_0175;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = IT_0111 + IT_0177;
    const ccomplex_t IT_0179 = IT_0114*IT_0178;
    const ccomplex_t IT_0180 = -IT_0118;
    const ccomplex_t IT_0181 = IT_0117*IT_0180;
    const ccomplex_t IT_0182 = -IT_0181;
    const ccomplex_t IT_0183 = IT_0120 + IT_0182;
    const ccomplex_t IT_0184 = IT_0122*IT_0183;
    const ccomplex_t IT_0185 = -IT_0124;
    const ccomplex_t IT_0186 = IT_0127*IT_0185;
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = IT_0128 + IT_0187;
    const ccomplex_t IT_0189 = IT_0130*IT_0188;
    const ccomplex_t IT_0190 = -IT_0132;
    const ccomplex_t IT_0191 = IT_0135*IT_0190;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = IT_0136 + IT_0192;
    const ccomplex_t IT_0194 = IT_0138*IT_0193;
    const ccomplex_t IT_0195 = IT_0040*IT_0042;
    const ccomplex_t IT_0196 = IT_0040*IT_0050;
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = IT_0195 + IT_0197;
    const ccomplex_t IT_0199 = IT_0039*IT_0198;
    const ccomplex_t IT_0200 = IT_0060*IT_0064;
    const ccomplex_t IT_0201 = IT_0060*IT_0068;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = IT_0059*IT_0203;
    const ccomplex_t IT_0205 = IT_0179 + IT_0184 + IT_0189 + IT_0194 + (-0.5)
      *IT_0199 + IT_0204;
    const ccomplex_t IT_0206 = s_14*s_34*IT_0092;
    const ccomplex_t IT_0207 = (s_12 + -IT_0147)*(s_13 + -IT_0206);
    const ccomplex_t IT_0208 = 3*conj(IT_0091);
    const ccomplex_t IT_0209 = s_34*IT_0092;
    const ccomplex_t IT_0210 = s_23*s_34*IT_0092;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = s_24 + IT_0211;
    const ccomplex_t IT_0213 = IT_0209*IT_0212;
    const ccomplex_t IT_0214 = (IT_0024 + -IT_0153)*(s_13 + -IT_0206);
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = (-3)*conj(IT_0143);
    const ccomplex_t IT_0217 = -IT_0051 + -IT_0161;
    const ccomplex_t IT_0218 = IT_0085 + IT_0217;
    const ccomplex_t IT_0219 = IT_0039*IT_0218;
    const ccomplex_t IT_0220 = -IT_0069 + -IT_0165;
    const ccomplex_t IT_0221 = IT_0089 + IT_0220;
    const ccomplex_t IT_0222 = IT_0059*IT_0221;
    const ccomplex_t IT_0223 = -IT_0179 + -IT_0184 + -IT_0189 + -IT_0194 + 0.5
      *IT_0219 + -IT_0222;
    const ccomplex_t IT_0224 = s_23*s_24*s_34*IT_0080;
    const ccomplex_t IT_0225 = pow(s_24, 2);
    const ccomplex_t IT_0226 = IT_0092*IT_0225;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = -IT_0153;
    const ccomplex_t IT_0229 = IT_0024 + IT_0224 + IT_0227 + IT_0228;
    const ccomplex_t IT_0230 = (s_12 + -IT_0145)*(IT_0024 + -IT_0153);
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = (s_12 + -IT_0147)*(IT_0024 + -IT_0226);
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = (-3)*conj(IT_0169);
    const ccomplex_t IT_0235 = IT_0039*IT_0195;
    const ccomplex_t IT_0236 = IT_0059*IT_0200;
    const ccomplex_t IT_0237 = (-0.5)*IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = s_13*s_34*IT_0092;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = s_14 + IT_0239;
    const ccomplex_t IT_0241 = IT_0209*IT_0240;
    const ccomplex_t IT_0242 = -(IT_0024 + -IT_0094)*(s_13 + -IT_0206);
    const ccomplex_t IT_0243 = s_24*s_34*IT_0092;
    const ccomplex_t IT_0244 = -(IT_0024 + -IT_0094)*(s_23 + -IT_0243);
    const ccomplex_t IT_0245 = (s_12 + -IT_0147)*(s_23 + -IT_0243);
    const ccomplex_t IT_0246 = (IT_0024 + -IT_0153)*(IT_0024 + -IT_0226);
    const ccomplex_t IT_0247 = 3*IT_0223;
    const ccomplex_t IT_0248 = (-3)*IT_0143;
    const ccomplex_t IT_0249 = (-3)*IT_0169;
    const ccomplex_t IT_0250 = (IT_0024 + -IT_0153)*(s_23 + -IT_0243);
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = -IT_0205;
    const ccomplex_t IT_0253 = 0.333333333333333*IT_0172;
    const ccomplex_t IT_0254 = -(IT_0024 + -IT_0094)*(IT_0024 + -IT_0226);
    const ccomplex_t IT_0255 = -IT_0254;
    const ccomplex_t IT_0256 = 3*IT_0169;
    const ccomplex_t IT_0257 = IT_0081*IT_0092;
    const ccomplex_t IT_0258 = (s_12 + -IT_0147)*(IT_0112 + -IT_0257);
    const ccomplex_t IT_0259 = (IT_0024 + -IT_0153)*(IT_0112 + -IT_0257);
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = (-0.333333333333333)*IT_0245;
    const ccomplex_t IT_0262 = -(IT_0024 + -IT_0094)*(IT_0112 + -IT_0257);
    const ccomplex_t IT_0263 = 6*IT_0079*conj(IT_0079)*(1 + 0.5*IT_0082) +
       conj(IT_0091)*(IT_0098*IT_0099 + IT_0101*IT_0102) + IT_0091*IT_0101
      *IT_0103 + conj(IT_0143)*(IT_0150*IT_0151 + IT_0155*IT_0156 + IT_0157
      *IT_0158) + IT_0143*(IT_0150*IT_0159 + IT_0157*IT_0160) + IT_0169*(IT_0159
      *IT_0171 + IT_0172*IT_0173 + IT_0160*IT_0174) + IT_0205*(IT_0207*IT_0208 +
       IT_0103*IT_0213 + IT_0215*IT_0216) + IT_0223*(IT_0172*IT_0208 + IT_0159
      *IT_0229 + IT_0216*IT_0231 + IT_0233*IT_0234) + IT_0237*(IT_0173*IT_0207 +
       IT_0103*IT_0241 + IT_0160*IT_0242 + IT_0234*IT_0244) + 3*conj(IT_0223)*
      (0.333333333333333*IT_0151*IT_0229 + IT_0237*IT_0245 + 0.333333333333333
      *IT_0246*IT_0247 + 0.333333333333333*IT_0231*IT_0248 + 0.333333333333333
      *IT_0233*IT_0249 + IT_0251*IT_0252 + IT_0099*IT_0253) + 3*conj(IT_0169)*
      (0.333333333333333*IT_0151*IT_0171 + 0.333333333333333*IT_0158*IT_0174 +
       IT_0205*IT_0245 + IT_0156*IT_0253 + 0.333333333333333*IT_0255*IT_0256) + 
      (-3)*conj(IT_0205)*((-0.333333333333333)*IT_0099*IT_0207 + (
      -0.333333333333333)*IT_0102*IT_0213 + (-0.333333333333333)*IT_0215*IT_0248
       + IT_0223*IT_0251 + -IT_0237*IT_0258 + IT_0205*IT_0260 + IT_0256*IT_0261)
       + (-3)*conj(IT_0237)*((-0.333333333333333)*IT_0156*IT_0207 + (
      -0.333333333333333)*IT_0102*IT_0241 + (-0.333333333333333)*IT_0158*IT_0242
       + (-0.333333333333333)*IT_0244*IT_0249 + IT_0252*IT_0258 + IT_0247
      *IT_0261 + IT_0237*IT_0262);
    return create_ccomplex_return(IT_0263);
}

