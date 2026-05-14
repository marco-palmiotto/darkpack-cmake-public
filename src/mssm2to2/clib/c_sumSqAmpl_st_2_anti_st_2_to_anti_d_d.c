#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_01*conj(U_st_01) + IT_0001*IT_0002*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0005 = 1.33333333333333*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0002;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0003;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0007 + 3*IT_0009);
    const ccomplex_t IT_0011 = (-0.166666666666667)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_st_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0005*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = cos(beta);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_td*e_em*IT_0003*IT_0018*IT_0020*conj(U_st_01);
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*conj
      (V_td)*e_em*IT_0003*IT_0018*IT_0020*U_st_01;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = pow(m_d, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_23 + IT_0012 + IT_0026 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0003*IT_0018*IT_0020*conj(U_st_01);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0003*IT_0018*IT_0020*U_st_01;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0012 + IT_0026 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0029 + IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0040 = 0.333333333333333*IT_0039;
    const ccomplex_t IT_0041 = IT_0013*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0005*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = sin(beta);
    const ccomplex_t IT_0046 = cpow(IT_0045, -1);
    const ccomplex_t IT_0047 = IT_0003*IT_0046;
    const ccomplex_t IT_0048 = m_t*V_td*V_u2*e_em*IT_0018*conj(U_st_11);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = V_td*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0054 = IT_0003*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0052 + IT_0057;
    const ccomplex_t IT_0059 = conj(V_td)*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0060 = IT_0003*IT_0059;
    const ccomplex_t IT_0061 = m_t*conj(V_td)*conj(V_u2)*e_em*IT_0018*U_st_11;
    const ccomplex_t IT_0062 = IT_0047*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0060 + (-0.5)*IT_0063);
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0035*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0058*IT_0067;
    const ccomplex_t IT_0069 = conj(V_td)*e_em*conj(V_Wp1)*U_st_01;
    const ccomplex_t IT_0070 = IT_0003*IT_0069;
    const ccomplex_t IT_0071 = m_t*conj(V_td)*conj(V_u1)*e_em*IT_0018*U_st_11;
    const ccomplex_t IT_0072 = IT_0047*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + (-0.5)*IT_0073);
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0027*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = m_t*V_td*V_u1*e_em*IT_0018*conj(U_st_11);
    const ccomplex_t IT_0079 = IT_0047*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = V_td*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0084 = IT_0003*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0082 + IT_0087;
    const ccomplex_t IT_0089 = IT_0077*IT_0088;
    const ccomplex_t IT_0090 = IT_0068 + IT_0089;
    const ccomplex_t IT_0091 = IT_0012*IT_0026;
    const ccomplex_t IT_0092 = cos(alpha);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_d*e_em*IT_0003*IT_0018
      *IT_0020*IT_0092;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = sin(alpha);
    const ccomplex_t IT_0096 = m_t*mu_h*IT_0018;
    const ccomplex_t IT_0097 = IT_0019*IT_0092;
    const ccomplex_t IT_0098 = -IT_0045*IT_0095 + IT_0097;
    const ccomplex_t IT_0099 = cpow(IT_0000, -2);
    const ccomplex_t IT_0100 = pow(m_t, 2);
    const ccomplex_t IT_0101 = (-2)*U_st_01*((0.5*m_W*IT_0098*(IT_0003 + -1./3
      *IT_0002*IT_0099) + IT_0003*IT_0018*IT_0046*IT_0095*IT_0100)*conj(U_st_01)
       + -1./2*IT_0003*IT_0046*IT_0092*IT_0096*conj(U_st_11)) + U_st_11*(IT_0003
      *IT_0046*IT_0092*IT_0096*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0098*IT_0099 + 1.5*IT_0003*IT_0018*IT_0046*IT_0095*IT_0100)*conj
      (U_st_11));
    const ccomplex_t IT_0102 = (0 + _Complex_I*-1)*(e_em*IT_0101 + (
      -1.4142135623731)*A_t*IT_0095*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0094*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_12 + (-2)*IT_0012 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_d*e_em*IT_0003*IT_0018
      *IT_0020*IT_0095;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = IT_0019*IT_0095;
    const ccomplex_t IT_0111 = IT_0045*IT_0092;
    const ccomplex_t IT_0112 = IT_0110 + IT_0111;
    const ccomplex_t IT_0113 = 2*U_st_01*((IT_0003*IT_0018*IT_0046*IT_0092
      *IT_0100 + (-0.5)*m_W*(IT_0003 + -1./3*IT_0002*IT_0099)*IT_0112)*conj
      (U_st_01) + 1./2*IT_0003*IT_0046*IT_0095*IT_0096*conj(U_st_11)) + U_st_11*
      (IT_0003*IT_0046*IT_0095*IT_0096*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0003*IT_0018*IT_0046*IT_0092*IT_0100 + m_W*IT_0002*IT_0099
      *IT_0112)*conj(U_st_11));
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(e_em*IT_0113 +
       1.4142135623731*A_t*IT_0092*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = IT_0109*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = -IT_0107 + -IT_0119;
    const ccomplex_t IT_0121 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0122 = m_d*IT_0121;
    const ccomplex_t IT_0123 = m_d*IT_0005;
    const ccomplex_t IT_0124 = IT_0015*(IT_0122 + -IT_0123);
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = IT_0042*(IT_0122 + -IT_0123);
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = 18*s_34;
    const ccomplex_t IT_0129 = (-18)*IT_0026;
    const ccomplex_t IT_0130 = IT_0128 + IT_0129;
    const ccomplex_t IT_0131 = s_34*IT_0012;
    const ccomplex_t IT_0132 = (-6)*IT_0131;
    const ccomplex_t IT_0133 = s_13*s_14;
    const ccomplex_t IT_0134 = 12*IT_0133;
    const ccomplex_t IT_0135 = IT_0132 + IT_0134;
    const ccomplex_t IT_0136 = (-18)*IT_0131;
    const ccomplex_t IT_0137 = 36*IT_0133;
    const ccomplex_t IT_0138 = IT_0136 + IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0140 = IT_0024*IT_0088;
    const ccomplex_t IT_0141 = IT_0027*IT_0139*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0143 = IT_0033*IT_0058;
    const ccomplex_t IT_0144 = IT_0035*IT_0142*IT_0143;
    const ccomplex_t IT_0145 = m_d*IT_0082;
    const ccomplex_t IT_0146 = m_d*IT_0086;
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = IT_0145 + IT_0147;
    const ccomplex_t IT_0149 = IT_0077*IT_0148;
    const ccomplex_t IT_0150 = m_d*IT_0056;
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = m_d*IT_0052;
    const ccomplex_t IT_0153 = IT_0151 + IT_0152;
    const ccomplex_t IT_0154 = IT_0067*IT_0153;
    const ccomplex_t IT_0155 = -IT_0141 + -IT_0144 + -IT_0149 + -IT_0154;
    const ccomplex_t IT_0156 = IT_0124 + IT_0127;
    const ccomplex_t IT_0157 = conj(IT_0124) + conj(IT_0127);
    const ccomplex_t IT_0158 = conj(IT_0125) + conj(IT_0126);
    const ccomplex_t IT_0159 = m_d*IT_0038;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0022*IT_0075;
    const ccomplex_t IT_0162 = IT_0027*IT_0139*IT_0161;
    const ccomplex_t IT_0163 = IT_0031*IT_0065;
    const ccomplex_t IT_0164 = IT_0035*IT_0142*IT_0163;
    const ccomplex_t IT_0165 = -IT_0162 + -IT_0164;
    const ccomplex_t IT_0166 = conj(IT_0160) + conj(IT_0165);
    const ccomplex_t IT_0167 = IT_0125 + IT_0126;
    const ccomplex_t IT_0168 = IT_0160 + IT_0165;
    const ccomplex_t IT_0169 = 6*s_34;
    const ccomplex_t IT_0170 = conj(IT_0124)*IT_0169;
    const ccomplex_t IT_0171 = -IT_0026;
    const ccomplex_t IT_0172 = 6*IT_0120;
    const ccomplex_t IT_0173 = 6*conj(IT_0120);
    const ccomplex_t IT_0174 = IT_0125*IT_0128;
    const ccomplex_t IT_0175 = IT_0128*IT_0165;
    const ccomplex_t IT_0176 = m_d*s_13;
    const ccomplex_t IT_0177 = 6*IT_0165;
    const ccomplex_t IT_0178 = 6*IT_0160;
    const ccomplex_t IT_0179 = 6*conj(IT_0125);
    const ccomplex_t IT_0180 = 6*conj(IT_0126);
    const ccomplex_t IT_0181 = 18*conj(IT_0155);
    const ccomplex_t IT_0182 = 6*conj(IT_0155);
    const ccomplex_t IT_0183 = 18*conj(IT_0125);
    const ccomplex_t IT_0184 = 18*conj(IT_0126);
    const ccomplex_t IT_0185 = 6*conj(IT_0165);
    const ccomplex_t IT_0186 = 6*conj(IT_0160);
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = 18*conj(IT_0124);
    const ccomplex_t IT_0189 = 18*conj(IT_0127);
    const ccomplex_t IT_0190 = 6*conj(IT_0124);
    const ccomplex_t IT_0191 = 6*conj(IT_0127);
    const ccomplex_t IT_0192 = 18*conj(IT_0165);
    const ccomplex_t IT_0193 = 18*conj(IT_0160);
    const ccomplex_t IT_0194 = 6*conj(IT_0090)*(IT_0125 + IT_0126 + 3*IT_0155)
       + conj(IT_0017)*(18*IT_0125 + 18*IT_0126 + 6*IT_0155) + 6*conj(IT_0038)*
      (IT_0124 + IT_0127 + 3*IT_0160 + 3*IT_0165) + conj(IT_0044)*(18*IT_0124 +
       18*IT_0127 + IT_0177 + IT_0178) + IT_0090*(IT_0179 + IT_0180 + IT_0181) +
       IT_0017*(IT_0182 + IT_0183 + IT_0184) + IT_0044*(IT_0187 + IT_0188 +
       IT_0189) + IT_0038*(IT_0190 + IT_0191 + IT_0192 + IT_0193);
    const ccomplex_t IT_0195 = (-18)*IT_0125;
    const ccomplex_t IT_0196 = (-18)*IT_0126;
    const ccomplex_t IT_0197 = (-6)*IT_0155;
    const ccomplex_t IT_0198 = IT_0195 + IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = (-18)*conj(IT_0125);
    const ccomplex_t IT_0200 = (-18)*conj(IT_0126);
    const ccomplex_t IT_0201 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0202 = IT_0199 + IT_0200 + IT_0201;
    const ccomplex_t IT_0203 = (-18)*IT_0155;
    const ccomplex_t IT_0204 = (-6)*IT_0125;
    const ccomplex_t IT_0205 = (-6)*IT_0126;
    const ccomplex_t IT_0206 = IT_0203 + IT_0204 + IT_0205;
    const ccomplex_t IT_0207 = (-18)*conj(IT_0155);
    const ccomplex_t IT_0208 = (-6)*conj(IT_0125);
    const ccomplex_t IT_0209 = (-6)*conj(IT_0126);
    const ccomplex_t IT_0210 = IT_0207 + IT_0208 + IT_0209;
    const ccomplex_t IT_0211 = (-18)*IT_0124;
    const ccomplex_t IT_0212 = (-18)*conj(IT_0124);
    const ccomplex_t IT_0213 = (-6)*IT_0124;
    const ccomplex_t IT_0214 = (-6)*conj(IT_0124);
    const ccomplex_t IT_0215 = m_d*e_em*IT_0003*IT_0018*IT_0020*IT_0045;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = m_t*e_em*mu_h*IT_0018*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0218 = IT_0003*IT_0217;
    const ccomplex_t IT_0219 = m_t*e_em*mu_h*IT_0018*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0220 = IT_0003*IT_0219;
    const ccomplex_t IT_0221 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0222 = IT_0019*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0225 = IT_0019*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (-0.5)*IT_0218 + 0.5*IT_0220 + (-0.5)*IT_0223 +
       0.5*IT_0226;
    const ccomplex_t IT_0228 = IT_0216*IT_0227;
    const ccomplex_t IT_0229 = cpow((-2)*s_12 + (-2)*IT_0012 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0230 = IT_0228*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = m_d*e_em*IT_0003*IT_0018;
    const ccomplex_t IT_0233 = 0.5*IT_0232;
    const ccomplex_t IT_0234 = IT_0003*IT_0019*IT_0046;
    const ccomplex_t IT_0235 = IT_0217*IT_0234;
    const ccomplex_t IT_0236 = IT_0219*IT_0234;
    const ccomplex_t IT_0237 = IT_0045*IT_0221;
    const ccomplex_t IT_0238 = 1.4142135623731*IT_0237;
    const ccomplex_t IT_0239 = IT_0045*IT_0224;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = 0.5*IT_0235 + (-0.5)*IT_0236 + (-0.5)*IT_0238 +
       0.5*IT_0240;
    const ccomplex_t IT_0242 = IT_0233*IT_0241;
    const ccomplex_t IT_0243 = pow(m_Z, 2);
    const ccomplex_t IT_0244 = cpow((-2)*s_12 + (-2)*IT_0012 + IT_0243 + 
      -reg_prop, -1);
    const ccomplex_t IT_0245 = IT_0242*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*IT_0245;
    const ccomplex_t IT_0247 = -IT_0231 + -IT_0246;
    const ccomplex_t IT_0248 = s_13 + s_14;
    const ccomplex_t IT_0249 = m_d*IT_0248;
    const ccomplex_t IT_0250 = 18*IT_0026;
    const ccomplex_t IT_0251 = IT_0128 + IT_0250;
    const ccomplex_t IT_0252 = (-18)*s_34;
    const ccomplex_t IT_0253 = IT_0129 + IT_0252;
    const ccomplex_t IT_0254 = s_34 + IT_0026;
    const ccomplex_t IT_0255 = m_d*s_14;
    const ccomplex_t IT_0256 = (-18)*IT_0165;
    const ccomplex_t IT_0257 = (-18)*IT_0160;
    const ccomplex_t IT_0258 = (-6)*IT_0127;
    const ccomplex_t IT_0259 = (-18)*IT_0127;
    const ccomplex_t IT_0260 = (-6)*IT_0165;
    const ccomplex_t IT_0261 = (-6)*IT_0160;
    const ccomplex_t IT_0262 = conj(IT_0044)*IT_0198 + IT_0044*IT_0202 + conj
      (IT_0038)*IT_0206 + IT_0038*IT_0210 + (-18)*IT_0017*(conj(IT_0127) + 1./3
      *conj(IT_0160) + 1./3*conj(IT_0165) + (-0.0555555555555555)*IT_0212) + (
      -18)*IT_0090*(0.333333333333333*conj(IT_0127) + conj(IT_0160) + conj
      (IT_0165) + (-0.0555555555555556)*IT_0214) + conj(IT_0090)*(IT_0213 +
       IT_0256 + IT_0257 + IT_0258) + conj(IT_0017)*(IT_0211 + IT_0259 + IT_0260
       + IT_0261);
    const ccomplex_t IT_0263 = -s_14;
    const ccomplex_t IT_0264 = s_13 + IT_0263;
    const ccomplex_t IT_0265 = m_d*IT_0264;
    const ccomplex_t IT_0266 = conj(IT_0038) + conj(IT_0090);
    const ccomplex_t IT_0267 = IT_0038 + IT_0090;
    const ccomplex_t IT_0268 = conj(IT_0017) + conj(IT_0044);
    const ccomplex_t IT_0269 = 18*IT_0120;
    const ccomplex_t IT_0270 = IT_0017 + IT_0044;
    const ccomplex_t IT_0271 = 18*conj(IT_0120);
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0273 = (-2)*IT_0272;
    const ccomplex_t IT_0274 = cpow(s_12 + IT_0012 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0275 = IT_0272*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*IT_0275;
    const ccomplex_t IT_0277 = IT_0273*IT_0276;
    const ccomplex_t IT_0278 = 0.25*IT_0277;
    const ccomplex_t IT_0279 = (-0.0833333333333333)*IT_0277;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0281 = (-0.333333333333333)*IT_0280;
    const ccomplex_t IT_0282 = (-1.33333333333333)*IT_0280;
    const ccomplex_t IT_0283 = IT_0281*IT_0282;
    const ccomplex_t IT_0284 = IT_0274*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = 0.5*IT_0285;
    const ccomplex_t IT_0287 = IT_0279 + IT_0286;
    const ccomplex_t IT_0288 = 36*conj(IT_0120);
    const ccomplex_t IT_0289 = IT_0172*IT_0266 + IT_0173*IT_0267 + IT_0268
      *IT_0269 + IT_0270*IT_0271 + ((-12)*IT_0120 + 12*conj(IT_0120) + IT_0179 +
       IT_0180 + IT_0181 + IT_0190 + IT_0191 + IT_0192 + IT_0193 + IT_0203 +
       IT_0204 + IT_0205 + IT_0213 + IT_0256 + IT_0257 + IT_0258)*IT_0278 +
       IT_0287*((-36)*IT_0120 + IT_0182 + IT_0183 + IT_0184 + IT_0185 + IT_0186 
      + IT_0188 + IT_0189 + IT_0195 + IT_0196 + IT_0197 + IT_0211 + IT_0259 +
       IT_0260 + IT_0261 + IT_0288);
    const ccomplex_t IT_0290 = (-12)*IT_0131;
    const ccomplex_t IT_0291 = (-12)*IT_0091;
    const ccomplex_t IT_0292 = 24*IT_0133;
    const ccomplex_t IT_0293 = IT_0290 + IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = (-2)*IT_0279;
    const ccomplex_t IT_0295 = (-36)*IT_0131;
    const ccomplex_t IT_0296 = (-36)*IT_0091;
    const ccomplex_t IT_0297 = 72*IT_0133;
    const ccomplex_t IT_0298 = IT_0295 + IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = cpow(IT_0279, 2);
    const ccomplex_t IT_0300 = cpow(IT_0278, 2);
    const ccomplex_t IT_0301 = (-6)*IT_0091;
    const ccomplex_t IT_0302 = IT_0132 + IT_0134 + IT_0301;
    const ccomplex_t IT_0303 = -IT_0279;
    const ccomplex_t IT_0304 = -IT_0286;
    const ccomplex_t IT_0305 = IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = -IT_0278;
    const ccomplex_t IT_0307 = (-18)*IT_0091;
    const ccomplex_t IT_0308 = IT_0136 + IT_0137 + IT_0307;
    const ccomplex_t IT_0309 = 36*s_34;
    const ccomplex_t IT_0310 = IT_0247*(((-18)*conj(IT_0017) + 6*conj(IT_0038)
       + 18*conj(IT_0044) + (-6)*conj(IT_0090))*IT_0249 + IT_0157*IT_0251 +
       IT_0158*IT_0253 + (IT_0187 + IT_0201)*IT_0254 + conj(IT_0247)*(36*IT_0026
       + IT_0309));
    const ccomplex_t IT_0311 = (-6)*(conj(IT_0017)*(IT_0038 + 3*IT_0044) +
       IT_0017*(conj(IT_0038) + 3*conj(IT_0044)) + 3*(conj(IT_0038) +
       0.333333333333333*conj(IT_0044))*IT_0090 + 3*(IT_0038 + 0.333333333333333
      *IT_0044)*conj(IT_0090))*IT_0091 + (conj(IT_0120)*(IT_0124 + IT_0125 +
       IT_0126 + IT_0127) + IT_0120*(conj(IT_0124) + conj(IT_0125) + conj
      (IT_0126) + conj(IT_0127)))*IT_0130 + (conj(IT_0038)*IT_0044 + IT_0038
      *conj(IT_0044) + conj(IT_0017)*IT_0090 + IT_0017*conj(IT_0090))*IT_0135 + 
      (IT_0017*conj(IT_0017) + IT_0038*conj(IT_0038) + IT_0044*conj(IT_0044) +
       IT_0090*conj(IT_0090))*IT_0138 + IT_0128*(IT_0155*conj(IT_0155) + IT_0156
      *IT_0157 + IT_0126*IT_0158 + IT_0160*IT_0166) + (IT_0155*IT_0158 + IT_0156
      *IT_0166 + conj(IT_0155)*IT_0167 + conj(IT_0127)*IT_0168)*IT_0169 +
       IT_0168*IT_0170 + (s_34 + IT_0171)*((conj(IT_0155) + conj(IT_0160) + conj
      (IT_0165))*IT_0172 + (IT_0155 + IT_0160 + IT_0165)*IT_0173) + IT_0158
      *IT_0174 + IT_0166*IT_0175 + IT_0176*IT_0194 + IT_0026*(conj(IT_0127)
      *IT_0198 + IT_0127*IT_0202 + IT_0166*IT_0206 + IT_0168*IT_0210 + IT_0158
      *IT_0211 + IT_0167*IT_0212 + conj(IT_0155)*IT_0213 + IT_0155*IT_0214) +
       conj(IT_0247)*(((-18)*IT_0017 + 6*IT_0038 + 18*IT_0044 + (-6)*IT_0090)
      *IT_0249 + IT_0156*IT_0251 + IT_0167*IT_0253 + (IT_0177 + IT_0178 +
       IT_0197)*IT_0254) + IT_0255*IT_0262 + IT_0265*IT_0289 + (-2)*IT_0278
      *IT_0293*(IT_0286 + (-0.5)*IT_0294) + -IT_0298*(IT_0286*(IT_0286 + 
      -IT_0294) + IT_0299 + IT_0300) + IT_0302*(IT_0268*IT_0278 + IT_0266
      *IT_0287 + IT_0267*IT_0305 + IT_0270*IT_0306) + (IT_0266*IT_0278 + IT_0268
      *IT_0287 + IT_0270*IT_0305 + IT_0267*IT_0306)*IT_0308 + IT_0120*conj
      (IT_0120)*((-36)*IT_0026 + IT_0309) + IT_0310;
    return create_ccomplex_return(IT_0311);
}

