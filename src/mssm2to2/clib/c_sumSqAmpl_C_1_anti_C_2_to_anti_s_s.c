#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_cs)*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0004;
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_s, 2);
    const ccomplex_t IT_0011 = pow(m_C_1, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0016;
    const ccomplex_t IT_0018 = 0.333333333333333*IT_0017;
    const ccomplex_t IT_0019 = IT_0003*IT_0016;
    const ccomplex_t IT_0020 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0004*IT_0015;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0021 + -IT_0023 + (-2)
      *IT_0025);
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0018*IT_0027;
    const ccomplex_t IT_0029 = pow(m_C_2, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0011 + -IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = m_s*s_24*m_C_1;
    const ccomplex_t IT_0034 = m_C_1*m_C_2;
    const ccomplex_t IT_0035 = -IT_0010;
    const ccomplex_t IT_0036 = s_34 + IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = 2*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + -IT_0011 + -IT_0029 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0040 = cos(alpha);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0004*IT_0040;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0044 = IT_0040*IT_0043;
    const ccomplex_t IT_0045 = sin(alpha);
    const ccomplex_t IT_0046 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = IT_0004*(IT_0044 + IT_0047);
    const ccomplex_t IT_0049 = 1.4142135623731*e_em*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = IT_0042*IT_0051;
    const ccomplex_t IT_0053 = IT_0039*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0004*IT_0045;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0040*IT_0046;
    const ccomplex_t IT_0058 = IT_0043*IT_0045;
    const ccomplex_t IT_0059 = IT_0004*(IT_0057 + -IT_0058);
    const ccomplex_t IT_0060 = 1.4142135623731*e_em*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0056*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0011 + -IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = -IT_0054 + -IT_0066;
    const ccomplex_t IT_0068 = s_34 + IT_0010;
    const ccomplex_t IT_0069 = IT_0034*IT_0068;
    const ccomplex_t IT_0070 = sin(beta);
    const ccomplex_t IT_0071 = m_s*e_em*IT_0000*IT_0002*IT_0004*IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = IT_0001*IT_0046;
    const ccomplex_t IT_0074 = IT_0043*IT_0070;
    const ccomplex_t IT_0075 = IT_0004*(IT_0073 + IT_0074);
    const ccomplex_t IT_0076 = 1.4142135623731*e_em*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = IT_0072*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + -IT_0011 + -IT_0029 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = IT_0001*IT_0043;
    const ccomplex_t IT_0083 = IT_0046*IT_0070;
    const ccomplex_t IT_0084 = IT_0004*(IT_0082 + -IT_0083);
    const ccomplex_t IT_0085 = 1.4142135623731*e_em*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = m_s*e_em*IT_0000*IT_0004;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = pow(m_Z, 2);
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + -IT_0011 + -IT_0029 + IT_0090 
      + -reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0089*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = -IT_0081 + -IT_0093;
    const ccomplex_t IT_0095 = m_s*s_13*m_C_2;
    const ccomplex_t IT_0096 = e_em*IT_0019;
    const ccomplex_t IT_0097 = e_em*IT_0022;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0096 + 3*IT_0097);
    const ccomplex_t IT_0099 = (-0.166666666666667)*IT_0098;
    const ccomplex_t IT_0100 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0101 = IT_0019*IT_0100;
    const ccomplex_t IT_0102 = IT_0022*IT_0100;
    const ccomplex_t IT_0103 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0104 = IT_0022*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + -IT_0102 + (-2)
      *IT_0104);
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0099*IT_0106;
    const ccomplex_t IT_0108 = IT_0030*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = V_u1*U_Wm2;
    const ccomplex_t IT_0112 = IT_0040*IT_0111;
    const ccomplex_t IT_0113 = U_d2*V_Wp1;
    const ccomplex_t IT_0114 = IT_0045*IT_0113;
    const ccomplex_t IT_0115 = IT_0004*(IT_0112 + -IT_0114);
    const ccomplex_t IT_0116 = 1.4142135623731*e_em*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0056*IT_0118;
    const ccomplex_t IT_0120 = IT_0064*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = IT_0040*IT_0113;
    const ccomplex_t IT_0123 = IT_0045*IT_0111;
    const ccomplex_t IT_0124 = IT_0004*(IT_0122 + IT_0123);
    const ccomplex_t IT_0125 = 1.4142135623731*e_em*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0042*IT_0127;
    const ccomplex_t IT_0129 = IT_0039*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = -IT_0121 + -IT_0130;
    const ccomplex_t IT_0132 = s_12*IT_0010;
    const ccomplex_t IT_0133 = s_14*s_23;
    const ccomplex_t IT_0134 = s_12*s_34;
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = s_13*s_24;
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = IT_0132 + IT_0133 + IT_0135 + IT_0137;
    const ccomplex_t IT_0139 = IT_0001*IT_0113;
    const ccomplex_t IT_0140 = IT_0070*IT_0111;
    const ccomplex_t IT_0141 = IT_0004*(IT_0139 + -IT_0140);
    const ccomplex_t IT_0142 = 1.4142135623731*e_em*IT_0141;
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = IT_0088*IT_0143;
    const ccomplex_t IT_0145 = IT_0091*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0001*IT_0111;
    const ccomplex_t IT_0148 = IT_0070*IT_0113;
    const ccomplex_t IT_0149 = IT_0004*(IT_0147 + IT_0148);
    const ccomplex_t IT_0150 = 1.4142135623731*e_em*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0072*IT_0151;
    const ccomplex_t IT_0153 = IT_0079*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = -IT_0146 + -IT_0154;
    const ccomplex_t IT_0156 = -IT_0133;
    const ccomplex_t IT_0157 = IT_0132 + IT_0134 + IT_0136 + IT_0156;
    const ccomplex_t IT_0158 = m_s*m_C_1;
    const ccomplex_t IT_0159 = -s_24;
    const ccomplex_t IT_0160 = s_23 + IT_0159;
    const ccomplex_t IT_0161 = IT_0158*IT_0160;
    const ccomplex_t IT_0162 = 4*IT_0161;
    const ccomplex_t IT_0163 = m_s*m_C_2;
    const ccomplex_t IT_0164 = -s_14;
    const ccomplex_t IT_0165 = s_13 + IT_0164;
    const ccomplex_t IT_0166 = IT_0163*IT_0165;
    const ccomplex_t IT_0167 = (-4)*IT_0166;
    const ccomplex_t IT_0168 = 3*IT_0032;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_us*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*V_us
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0172 = 0.5*IT_0171;
    const ccomplex_t IT_0173 = IT_0170*IT_0172;
    const ccomplex_t IT_0174 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0175 = IT_0173*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = m_C_1*m_C_2*IT_0010;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0000*IT_0002*IT_0004*conj(U_st_00);
    const ccomplex_t IT_0180 = (-0.5)*IT_0179;
    const ccomplex_t IT_0181 = conj(V_ts)*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0182 = IT_0004*IT_0181;
    const ccomplex_t IT_0183 = cpow(IT_0070, -1);
    const ccomplex_t IT_0184 = IT_0004*IT_0183;
    const ccomplex_t IT_0185 = m_t*conj(V_ts)*conj(V_u2)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0186 = IT_0184*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0182 + (-0.5)*IT_0187);
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0180*IT_0189;
    const ccomplex_t IT_0191 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0192 = IT_0190*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0000*IT_0004*IT_0183;
    const ccomplex_t IT_0196 = 0.5*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u2)*e_em*IT_0000*IT_0004*IT_0183;
    const ccomplex_t IT_0198 = (-0.5)*IT_0197;
    const ccomplex_t IT_0199 = IT_0196*IT_0198;
    const ccomplex_t IT_0200 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0201 = IT_0199*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp2)
      *IT_0004;
    const ccomplex_t IT_0204 = IT_0008*IT_0203;
    const ccomplex_t IT_0205 = IT_0012*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = -IT_0202 + -IT_0206;
    const ccomplex_t IT_0208 = 12*IT_0014;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0006*IT_0210;
    const ccomplex_t IT_0212 = IT_0012*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0000*IT_0002*IT_0004*U_st_00;
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = IT_0180*IT_0215;
    const ccomplex_t IT_0217 = IT_0191*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0000*IT_0002*IT_0004*conj(U_st_01);
    const ccomplex_t IT_0220 = (-0.5)*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0000*IT_0002*IT_0004*U_st_01;
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = IT_0220*IT_0222;
    const ccomplex_t IT_0224 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0225 = IT_0223*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = IT_0213 + IT_0218 + IT_0226;
    const ccomplex_t IT_0228 = m_s*s_23*m_C_1;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*V_us*e_em*V_Wp1*IT_0004;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0172*IT_0230;
    const ccomplex_t IT_0232 = IT_0174*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = m_t*V_ts*V_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0235 = IT_0184*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = (-0.5)*IT_0237;
    const ccomplex_t IT_0239 = V_ts*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0240 = IT_0004*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = -IT_0241;
    const ccomplex_t IT_0243 = -IT_0242;
    const ccomplex_t IT_0244 = IT_0238 + IT_0243;
    const ccomplex_t IT_0245 = IT_0215*IT_0244;
    const ccomplex_t IT_0246 = IT_0191*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = m_t*V_ts*V_u1*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0249 = IT_0184*IT_0248;
    const ccomplex_t IT_0250 = 1.4142135623731*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = (-0.5)*IT_0251;
    const ccomplex_t IT_0253 = V_ts*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0254 = IT_0004*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = -IT_0255;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0252 + IT_0257;
    const ccomplex_t IT_0259 = IT_0222*IT_0258;
    const ccomplex_t IT_0260 = IT_0224*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = -IT_0233 + IT_0247 + IT_0261;
    const ccomplex_t IT_0263 = s_13 + s_14;
    const ccomplex_t IT_0264 = IT_0163*IT_0263;
    const ccomplex_t IT_0265 = s_34*m_C_1*m_C_2;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*V_us*e_em*conj(V_Wp2)*IT_0004;
    const ccomplex_t IT_0267 = IT_0230*IT_0266;
    const ccomplex_t IT_0268 = IT_0174*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = IT_0189*IT_0244;
    const ccomplex_t IT_0271 = IT_0191*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = conj(V_ts)*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0274 = IT_0004*IT_0273;
    const ccomplex_t IT_0275 = m_t*conj(V_ts)*conj(V_u2)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0276 = IT_0184*IT_0275;
    const ccomplex_t IT_0277 = 1.4142135623731*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*(IT_0274 + (-0.5)*IT_0277);
    const ccomplex_t IT_0279 = -IT_0278;
    const ccomplex_t IT_0280 = IT_0258*IT_0279;
    const ccomplex_t IT_0281 = IT_0224*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_us
      *e_em*IT_0000*IT_0004*IT_0183;
    const ccomplex_t IT_0284 = 0.5*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *V_us*e_em*IT_0000*IT_0004*IT_0183;
    const ccomplex_t IT_0286 = (-0.5)*IT_0285;
    const ccomplex_t IT_0287 = IT_0284*IT_0286;
    const ccomplex_t IT_0288 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0289 = IT_0287*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*IT_0289;
    const ccomplex_t IT_0291 = IT_0269 + IT_0272 + IT_0282 + IT_0290;
    const ccomplex_t IT_0292 = m_s*s_14*m_C_2;
    const ccomplex_t IT_0293 = IT_0203*IT_0210;
    const ccomplex_t IT_0294 = IT_0012*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*IT_0294;
    const ccomplex_t IT_0296 = IT_0170*IT_0266;
    const ccomplex_t IT_0297 = IT_0174*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = IT_0220*IT_0279;
    const ccomplex_t IT_0300 = IT_0224*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*IT_0300;
    const ccomplex_t IT_0302 = -IT_0295 + -IT_0298 + IT_0301;
    const ccomplex_t IT_0303 = IT_0027*IT_0099;
    const ccomplex_t IT_0304 = IT_0030*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*IT_0304;
    const ccomplex_t IT_0306 = 12*IT_0305;
    const ccomplex_t IT_0307 = IT_0018*IT_0106;
    const ccomplex_t IT_0308 = IT_0030*IT_0307;
    const ccomplex_t IT_0309 = (0 + _Complex_I*1)*IT_0308;
    const ccomplex_t IT_0310 = -IT_0309;
    const ccomplex_t IT_0311 = 24*IT_0032;
    const ccomplex_t IT_0312 = 48*IT_0032;
    const ccomplex_t IT_0313 = s_23 + s_24;
    const ccomplex_t IT_0314 = IT_0158*IT_0313;
    const ccomplex_t IT_0315 = 24*conj(IT_0305);
    const ccomplex_t IT_0316 = (-24)*conj(IT_0014);
    const ccomplex_t IT_0317 = (-24)*IT_0032;
    const ccomplex_t IT_0318 = (-12)*IT_0014;
    const ccomplex_t IT_0319 = (-12)*conj(IT_0014);
    const ccomplex_t IT_0320 = IT_0305*(24*conj(IT_0032)*IT_0132 + 3*conj
      (IT_0131)*IT_0162 + 3*conj(IT_0067)*IT_0167 + 24*conj(IT_0110)*IT_0265 +
       48*IT_0133*conj(IT_0305) + 48*IT_0178*conj(IT_0310) + IT_0228*IT_0319);
    const ccomplex_t IT_0321 = (-12)*IT_0305;
    const ccomplex_t IT_0322 = 3*IT_0162;
    const ccomplex_t IT_0323 = conj(IT_0067)*IT_0322;
    const ccomplex_t IT_0324 = 3*IT_0167;
    const ccomplex_t IT_0325 = conj(IT_0131)*IT_0324;
    const ccomplex_t IT_0326 = 12*IT_0264;
    const ccomplex_t IT_0327 = 24*IT_0132;
    const ccomplex_t IT_0328 = 24*IT_0136;
    const ccomplex_t IT_0329 = 24*IT_0033;
    const ccomplex_t IT_0330 = 24*IT_0095;
    const ccomplex_t IT_0331 = 24*IT_0178;
    const ccomplex_t IT_0332 = (-24)*IT_0136;
    const ccomplex_t IT_0333 = (-24)*IT_0033;
    const ccomplex_t IT_0334 = (-24)*IT_0178;
    const ccomplex_t IT_0335 = (-12)*IT_0314;
    const ccomplex_t IT_0336 = 48*conj(IT_0110)*IT_0136 + 48*conj(IT_0032)
      *IT_0178 + IT_0265*IT_0315 + IT_0095*IT_0316 + IT_0323 + IT_0325 + conj
      (IT_0155)*IT_0326 + conj(IT_0310)*IT_0327 + conj(IT_0291)*IT_0328 + conj
      (IT_0302)*IT_0329 + conj(IT_0262)*IT_0330 + conj(IT_0227)*IT_0331 + conj
      (IT_0207)*IT_0332 + conj(IT_0194)*IT_0333 + conj(IT_0177)*IT_0334 + conj
      (IT_0094)*IT_0335;
    const ccomplex_t IT_0337 = 12*IT_0132;
    const ccomplex_t IT_0338 = 12*IT_0314;
    const ccomplex_t IT_0339 = 12*IT_0292;
    const ccomplex_t IT_0340 = 12*IT_0228;
    const ccomplex_t IT_0341 = 12*IT_0265;
    const ccomplex_t IT_0342 = 24*conj(IT_0032);
    const ccomplex_t IT_0343 = (-12)*IT_0132;
    const ccomplex_t IT_0344 = (-12)*IT_0228;
    const ccomplex_t IT_0345 = (-12)*IT_0265;
    const ccomplex_t IT_0346 = (-12)*IT_0264;
    const ccomplex_t IT_0347 = 48*IT_0178*conj(IT_0305) + 48*IT_0133*conj
      (IT_0310) + IT_0292*IT_0319 + IT_0323 + IT_0325 + conj(IT_0110)*IT_0327 +
       conj(IT_0291)*IT_0337 + conj(IT_0094)*IT_0338 + conj(IT_0262)*IT_0339 +
       conj(IT_0302)*IT_0340 + conj(IT_0227)*IT_0341 + IT_0265*IT_0342 + conj
      (IT_0207)*IT_0343 + conj(IT_0194)*IT_0344 + conj(IT_0177)*IT_0345 + conj
      (IT_0155)*IT_0346;
    const ccomplex_t IT_0348 = 6*IT_0157;
    const ccomplex_t IT_0349 = 6*IT_0314;
    const ccomplex_t IT_0350 = 6*IT_0264;
    const ccomplex_t IT_0351 = 12*conj(IT_0305);
    const ccomplex_t IT_0352 = s_12*IT_0068;
    const ccomplex_t IT_0353 = 24*IT_0352;
    const ccomplex_t IT_0354 = (-24)*IT_0069;
    const ccomplex_t IT_0355 = (-6)*IT_0157;
    const ccomplex_t IT_0356 = (-6)*IT_0314;
    const ccomplex_t IT_0357 = (-6)*IT_0069;
    const ccomplex_t IT_0358 = (-6)*IT_0264;
    const ccomplex_t IT_0359 = conj(IT_0014)*IT_0069 + 0.166666666666667*conj
      (IT_0110)*IT_0335 + 0.166666666666667*conj(IT_0310)*IT_0338 +
       0.166666666666667*conj(IT_0194)*IT_0348 + 0.166666666666667*conj(IT_0207)
      *IT_0349 + 0.166666666666667*conj(IT_0177)*IT_0350 + (-2)*IT_0264*(conj
      (IT_0032) + (-0.0833333333333333)*IT_0351) + 0.166666666666667*conj
      (IT_0094)*IT_0353 + 0.166666666666667*conj(IT_0155)*IT_0354 +
       0.166666666666667*conj(IT_0302)*IT_0355 + 0.166666666666667*conj(IT_0291)
      *IT_0356 + 0.166666666666667*conj(IT_0262)*IT_0357 + 0.166666666666667
      *conj(IT_0227)*IT_0358;
    const ccomplex_t IT_0360 = 6*IT_0069;
    const ccomplex_t IT_0361 = (-12)*conj(IT_0305);
    const ccomplex_t IT_0362 = (-6)*conj(IT_0014);
    const ccomplex_t IT_0363 = conj(IT_0110)*IT_0326 + conj(IT_0310)*IT_0346 +
       conj(IT_0262)*IT_0348 + conj(IT_0227)*IT_0349 + conj(IT_0291)*IT_0350 +
       conj(IT_0155)*IT_0353 + conj(IT_0094)*IT_0354 + conj(IT_0177)*IT_0356 +
       conj(IT_0194)*IT_0357 + conj(IT_0207)*IT_0358 + conj(IT_0302)*IT_0360 +
       IT_0314*(12*conj(IT_0032) + IT_0361) + IT_0157*IT_0362;
    const ccomplex_t IT_0364 = 12*conj(IT_0014);
    const ccomplex_t IT_0365 = 12*IT_0136;
    const ccomplex_t IT_0366 = 12*IT_0033;
    const ccomplex_t IT_0367 = 12*IT_0095;
    const ccomplex_t IT_0368 = (-24)*conj(IT_0032);
    const ccomplex_t IT_0369 = (-12)*IT_0136;
    const ccomplex_t IT_0370 = (-12)*IT_0033;
    const ccomplex_t IT_0371 = (-12)*IT_0095;
    const ccomplex_t IT_0372 = (-12)*IT_0178;
    const ccomplex_t IT_0373 = (-6)*IT_0138;
    const ccomplex_t IT_0374 = (-3)*IT_0038;
    const ccomplex_t IT_0375 = conj(IT_0110)*IT_0333 + conj(IT_0310)*IT_0344 +
       conj(IT_0094)*IT_0348 + conj(IT_0155)*IT_0357 + IT_0292*IT_0361 + IT_0178
      *IT_0364 + conj(IT_0194)*IT_0365 + conj(IT_0207)*IT_0366 + conj(IT_0177)
      *IT_0367 + IT_0095*IT_0368 + conj(IT_0302)*IT_0369 + conj(IT_0291)*IT_0370
       + conj(IT_0227)*IT_0371 + conj(IT_0262)*IT_0372 + conj(IT_0067)*IT_0373 +
       conj(IT_0131)*IT_0374;
    const ccomplex_t IT_0376 = 6*IT_0138;
    const ccomplex_t IT_0377 = 12*IT_0178;
    const ccomplex_t IT_0378 = (-2)*IT_0037;
    const ccomplex_t IT_0379 = (-3)*IT_0378;
    const ccomplex_t IT_0380 = IT_0178*IT_0319 + conj(IT_0110)*IT_0329 + conj
      (IT_0310)*IT_0340 + IT_0095*IT_0342 + IT_0292*IT_0351 + conj(IT_0094)
      *IT_0355 + conj(IT_0155)*IT_0360 + conj(IT_0302)*IT_0365 + conj(IT_0291)
      *IT_0366 + conj(IT_0227)*IT_0367 + conj(IT_0194)*IT_0369 + conj(IT_0207)
      *IT_0370 + conj(IT_0177)*IT_0371 + conj(IT_0067)*IT_0376 + conj(IT_0262)
      *IT_0377 + conj(IT_0131)*IT_0379;
    const ccomplex_t IT_0381 = IT_0136*IT_0319 + conj(IT_0110)*IT_0330 + conj
      (IT_0310)*IT_0339 + IT_0033*IT_0342 + conj(IT_0155)*IT_0348 + IT_0228
      *IT_0351 + conj(IT_0094)*IT_0357 + conj(IT_0262)*IT_0365 + conj(IT_0227)
      *IT_0366 + conj(IT_0291)*IT_0367 + conj(IT_0177)*IT_0370 + conj(IT_0207)
      *IT_0371 + conj(IT_0194)*IT_0372 + conj(IT_0131)*IT_0376 + conj(IT_0302)
      *IT_0377 + conj(IT_0067)*IT_0379;
    const ccomplex_t IT_0382 = (-2)*IT_0161;
    const ccomplex_t IT_0383 = (-3)*IT_0382;
    const ccomplex_t IT_0384 = 2*IT_0166;
    const ccomplex_t IT_0385 = (-3)*IT_0384;
    const ccomplex_t IT_0386 = IT_0095*IT_0319 + conj(IT_0110)*IT_0328 + conj
      (IT_0310)*IT_0337 + IT_0178*IT_0342 + conj(IT_0155)*IT_0350 + IT_0265
      *IT_0351 + conj(IT_0094)*IT_0356 + conj(IT_0291)*IT_0365 + conj(IT_0302)
      *IT_0366 + conj(IT_0262)*IT_0367 + conj(IT_0207)*IT_0369 + conj(IT_0194)
      *IT_0370 + conj(IT_0177)*IT_0372 + conj(IT_0227)*IT_0377 + conj(IT_0067)
      *IT_0383 + conj(IT_0131)*IT_0385;
    const ccomplex_t IT_0387 = IT_0033*IT_0319 + conj(IT_0110)*IT_0331 + conj
      (IT_0310)*IT_0341 + IT_0136*IT_0342 + conj(IT_0155)*IT_0349 + IT_0132
      *IT_0351 + conj(IT_0094)*IT_0358 + conj(IT_0227)*IT_0365 + conj(IT_0262)
      *IT_0366 + conj(IT_0302)*IT_0367 + conj(IT_0177)*IT_0369 + conj(IT_0194)
      *IT_0371 + conj(IT_0207)*IT_0372 + conj(IT_0291)*IT_0377 + conj(IT_0131)
      *IT_0383 + conj(IT_0067)*IT_0385;
    const ccomplex_t IT_0388 = 2*IT_0161;
    const ccomplex_t IT_0389 = (-3)*IT_0388;
    const ccomplex_t IT_0390 = (-2)*IT_0166;
    const ccomplex_t IT_0391 = (-3)*IT_0390;
    const ccomplex_t IT_0392 = conj(IT_0110)*IT_0332 + conj(IT_0310)*IT_0343 +
       conj(IT_0094)*IT_0349 + conj(IT_0155)*IT_0358 + IT_0265*IT_0361 + IT_0095
      *IT_0364 + conj(IT_0207)*IT_0365 + conj(IT_0194)*IT_0366 + IT_0178*IT_0368
       + conj(IT_0291)*IT_0369 + conj(IT_0302)*IT_0370 + conj(IT_0262)*IT_0371 +
       conj(IT_0227)*IT_0372 + conj(IT_0177)*IT_0377 + conj(IT_0067)*IT_0389 +
       conj(IT_0131)*IT_0391;
    const ccomplex_t IT_0393 = conj(IT_0110)*IT_0334 + conj(IT_0310)*IT_0345 +
       conj(IT_0094)*IT_0350 + conj(IT_0155)*IT_0356 + IT_0132*IT_0361 + IT_0033
      *IT_0364 + conj(IT_0177)*IT_0365 + conj(IT_0194)*IT_0367 + IT_0136*IT_0368
       + conj(IT_0227)*IT_0369 + conj(IT_0262)*IT_0370 + conj(IT_0302)*IT_0371 +
       conj(IT_0291)*IT_0372 + conj(IT_0207)*IT_0377 + conj(IT_0131)*IT_0389 +
       conj(IT_0067)*IT_0391;
    const ccomplex_t IT_0394 = conj(IT_0110) + conj(IT_0310);
    const ccomplex_t IT_0395 = (-8)*IT_0037;
    const ccomplex_t IT_0396 = 3*IT_0395;
    const ccomplex_t IT_0397 = -s_34;
    const ccomplex_t IT_0398 = IT_0010 + IT_0397;
    const ccomplex_t IT_0399 = s_12*IT_0398;
    const ccomplex_t IT_0400 = (-24)*IT_0399;
    const ccomplex_t IT_0401 = 3*conj(IT_0032);
    const ccomplex_t IT_0402 = 0.333333333333333*IT_0401;
    const ccomplex_t IT_0403 = IT_0138*IT_0362 + conj(IT_0194)*IT_0374 + conj
      (IT_0262)*IT_0376 + conj(IT_0302)*IT_0379 + conj(IT_0227)*IT_0383 + conj
      (IT_0291)*IT_0385 + conj(IT_0177)*IT_0389 + conj(IT_0207)*IT_0391 +
       IT_0324*IT_0394 + conj(IT_0067)*IT_0396 + conj(IT_0131)*IT_0400 + 3
      *IT_0162*(conj(IT_0305) + IT_0402);
    const ccomplex_t IT_0404 = conj(IT_0014)*IT_0038 + (-0.333333333333333)
      *conj(IT_0194)*IT_0373 + (-0.333333333333333)*conj(IT_0302)*IT_0376 + (
      -0.333333333333333)*conj(IT_0262)*IT_0379 + (-0.333333333333333)*conj
      (IT_0291)*IT_0383 + (-0.333333333333333)*conj(IT_0227)*IT_0385 + (
      -0.333333333333333)*conj(IT_0207)*IT_0389 + (-0.333333333333333)*conj
      (IT_0177)*IT_0391 + (-0.333333333333333)*IT_0322*IT_0394 + (
      -0.333333333333333)*conj(IT_0131)*IT_0396 + (-0.333333333333333)*conj
      (IT_0067)*IT_0400 + -IT_0167*(conj(IT_0305) + IT_0402);
    const ccomplex_t IT_0405 = IT_0014*((-24)*conj(IT_0032)*IT_0033 + (-3)
      *IT_0038*conj(IT_0067) + 6*IT_0069*conj(IT_0094) + (-24)*IT_0095*conj
      (IT_0110) + (-6)*conj(IT_0131)*IT_0138 + (-6)*conj(IT_0155)*IT_0157) + 
      (conj(IT_0131)*IT_0162 + conj(IT_0067)*IT_0167)*IT_0168 + (conj(IT_0014)
      *IT_0136 + IT_0033*conj(IT_0177) + IT_0178*conj(IT_0194) + IT_0095*conj
      (IT_0207))*IT_0208 + (IT_0132*conj(IT_0227) + IT_0228*conj(IT_0262) + conj
      (IT_0094)*IT_0264 + IT_0265*conj(IT_0291) + IT_0292*conj(IT_0302))*IT_0306
       + (IT_0136*conj(IT_0227) + IT_0033*conj(IT_0262) + IT_0178*conj(IT_0291) 
      + IT_0095*conj(IT_0302) + IT_0265*conj(IT_0310))*IT_0311 + (conj(IT_0032)
      *IT_0136 + conj(IT_0110)*IT_0178)*IT_0312 + IT_0032*((-12)*conj(IT_0094)
      *IT_0264 + 12*conj(IT_0155)*IT_0314 + IT_0132*IT_0315 + IT_0033*IT_0316) +
       (IT_0136*conj(IT_0177) + IT_0095*conj(IT_0194) + IT_0178*conj(IT_0207))
      *IT_0317 + (IT_0033*conj(IT_0227) + IT_0136*conj(IT_0262) + IT_0095*conj
      (IT_0291) + IT_0178*conj(IT_0302) + IT_0228*conj(IT_0305) + IT_0292*conj
      (IT_0310))*IT_0318 + IT_0320 + (IT_0132*conj(IT_0177) + conj(IT_0207)
      *IT_0265 + conj(IT_0194)*IT_0292 + conj(IT_0155)*IT_0314)*IT_0321 +
       IT_0110*IT_0336 + IT_0310*IT_0347 + 6*IT_0094*IT_0359 + IT_0155*IT_0363 +
       IT_0194*IT_0375 + IT_0302*IT_0380 + IT_0262*IT_0381 + IT_0291*IT_0386 +
       IT_0227*IT_0387 + IT_0207*IT_0392 + IT_0177*IT_0393 + IT_0131*IT_0403 + (
      -3)*IT_0067*IT_0404;
    return create_ccomplex_return(IT_0405);
}

