#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = e_em*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + 3*IT_0007);
    const ccomplex_t IT_0009 = (-0.166666666666667)*IT_0008;
    const ccomplex_t IT_0010 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0010;
    const ccomplex_t IT_0013 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + -IT_0012 + (-2)
      *IT_0014);
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_2, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0023 = 0.333333333333333*IT_0022;
    const ccomplex_t IT_0024 = IT_0016*IT_0023;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = pow(m_d, 2);
    const ccomplex_t IT_0028 = s_12*IT_0027;
    const ccomplex_t IT_0029 = s_14*s_23;
    const ccomplex_t IT_0030 = m_d*m_C_2;
    const ccomplex_t IT_0031 = -s_24;
    const ccomplex_t IT_0032 = s_23 + IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = 4*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + (-2)*IT_0018 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0036 = pow(m_W, -1);
    const ccomplex_t IT_0037 = cos(beta);
    const ccomplex_t IT_0038 = cpow(IT_0037, -1);
    const ccomplex_t IT_0039 = cos(alpha);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_d*e_em*IT_0005*IT_0036
      *IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = U_d2*V_Wp2;
    const ccomplex_t IT_0043 = IT_0039*IT_0042;
    const ccomplex_t IT_0044 = sin(alpha);
    const ccomplex_t IT_0045 = V_u2*U_Wm2;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = IT_0005*(IT_0043 + IT_0046);
    const ccomplex_t IT_0048 = 1.4142135623731*e_em*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0041*IT_0050;
    const ccomplex_t IT_0052 = IT_0035*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_d*e_em*IT_0005*IT_0036
      *IT_0038*IT_0044;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0039*IT_0045;
    const ccomplex_t IT_0058 = IT_0042*IT_0044;
    const ccomplex_t IT_0059 = IT_0005*(IT_0057 + -IT_0058);
    const ccomplex_t IT_0060 = 1.4142135623731*e_em*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0056*IT_0062;
    const ccomplex_t IT_0064 = IT_0054*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0053 + -IT_0065;
    const ccomplex_t IT_0067 = s_34*IT_0018;
    const ccomplex_t IT_0068 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0069 = IT_0003*IT_0068;
    const ccomplex_t IT_0070 = IT_0006*IT_0068;
    const ccomplex_t IT_0071 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0072 = IT_0006*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + -IT_0070 + (-2)
      *IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0009*IT_0074;
    const ccomplex_t IT_0076 = IT_0019*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = IT_0018*IT_0027;
    const ccomplex_t IT_0080 = IT_0023*IT_0074;
    const ccomplex_t IT_0081 = IT_0019*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = -s_14;
    const ccomplex_t IT_0085 = s_13 + IT_0084;
    const ccomplex_t IT_0086 = IT_0030*IT_0085;
    const ccomplex_t IT_0087 = (-4)*IT_0086;
    const ccomplex_t IT_0088 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0089 = IT_0039*IT_0088;
    const ccomplex_t IT_0090 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0091 = IT_0044*IT_0090;
    const ccomplex_t IT_0092 = IT_0005*(IT_0089 + IT_0091);
    const ccomplex_t IT_0093 = 1.4142135623731*e_em*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = IT_0041*IT_0095;
    const ccomplex_t IT_0097 = IT_0035*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0039*IT_0090;
    const ccomplex_t IT_0100 = IT_0044*IT_0088;
    const ccomplex_t IT_0101 = IT_0005*(IT_0099 + -IT_0100);
    const ccomplex_t IT_0102 = 1.4142135623731*e_em*IT_0101;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0056*IT_0104;
    const ccomplex_t IT_0106 = IT_0054*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = -IT_0098 + -IT_0107;
    const ccomplex_t IT_0109 = s_23 + s_24;
    const ccomplex_t IT_0110 = IT_0030*IT_0109;
    const ccomplex_t IT_0111 = cpow((-2)*s_12 + (-2)*IT_0018 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0112 = sin(beta);
    const ccomplex_t IT_0113 = m_d*e_em*IT_0005*IT_0036*IT_0038*IT_0112;
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = IT_0037*IT_0045;
    const ccomplex_t IT_0116 = IT_0042*IT_0112;
    const ccomplex_t IT_0117 = IT_0005*(IT_0115 + IT_0116);
    const ccomplex_t IT_0118 = 1.4142135623731*e_em*IT_0117;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0114*IT_0119;
    const ccomplex_t IT_0121 = IT_0111*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = pow(m_Z, 2);
    const ccomplex_t IT_0124 = cpow((-2)*s_12 + (-2)*IT_0018 + IT_0123 + 
      -reg_prop, -1);
    const ccomplex_t IT_0125 = IT_0037*IT_0042;
    const ccomplex_t IT_0126 = IT_0045*IT_0112;
    const ccomplex_t IT_0127 = IT_0005*(IT_0125 + -IT_0126);
    const ccomplex_t IT_0128 = 1.4142135623731*e_em*IT_0127;
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = m_d*e_em*IT_0005*IT_0036;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = IT_0124*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = -IT_0122 + -IT_0134;
    const ccomplex_t IT_0136 = 24*conj(IT_0026)*IT_0028 + 48*conj(IT_0021)
      *IT_0029 + 3*IT_0034*conj(IT_0066) + 24*IT_0067*conj(IT_0078) + 48*IT_0079
      *conj(IT_0083) + 3*IT_0087*conj(IT_0108) + (-12)*IT_0110*conj(IT_0135);
    const ccomplex_t IT_0137 = 3*IT_0026;
    const ccomplex_t IT_0138 = m_d*s_23*m_C_2;
    const ccomplex_t IT_0139 = cpow((-2)*s_13 + IT_0018 + IT_0027 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*V_ud
      *e_em*IT_0005*IT_0036*IT_0038;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = IT_0141*IT_0143;
    const ccomplex_t IT_0145 = IT_0139*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = cpow((-2)*s_13 + IT_0018 + IT_0027 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*V_cd*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_cd)*e_em*IT_0005*IT_0036*IT_0038;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0149*IT_0151;
    const ccomplex_t IT_0153 = IT_0147*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = cpow((-2)*s_13 + IT_0018 + IT_0027 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0005*IT_0036*IT_0038*U_st_01;
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = cpow(IT_0112, -1);
    const ccomplex_t IT_0159 = IT_0005*IT_0158;
    const ccomplex_t IT_0160 = m_t*V_td*V_u2*e_em*IT_0036*conj(U_st_11);
    const ccomplex_t IT_0161 = IT_0159*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = V_td*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0166 = IT_0005*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = -IT_0168;
    const ccomplex_t IT_0170 = IT_0164 + IT_0169;
    const ccomplex_t IT_0171 = IT_0157*IT_0170;
    const ccomplex_t IT_0172 = IT_0155*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = cpow((-2)*s_13 + IT_0018 + IT_0027 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0175 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0005*IT_0036*IT_0038*U_st_00;
    const ccomplex_t IT_0176 = 0.5*IT_0175;
    const ccomplex_t IT_0177 = m_t*V_td*V_u2*e_em*IT_0036*conj(U_st_10);
    const ccomplex_t IT_0178 = IT_0159*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = V_td*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0183 = IT_0005*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = IT_0181 + IT_0186;
    const ccomplex_t IT_0188 = IT_0176*IT_0187;
    const ccomplex_t IT_0189 = IT_0174*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = -IT_0146 + -IT_0154 + IT_0173 + IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0005*IT_0036*IT_0038;
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = IT_0141*IT_0193;
    const ccomplex_t IT_0195 = IT_0139*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0005*IT_0036*IT_0038*conj(U_st_00);
    const ccomplex_t IT_0198 = (-0.5)*IT_0197;
    const ccomplex_t IT_0199 = IT_0176*IT_0198;
    const ccomplex_t IT_0200 = IT_0174*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0005*IT_0036*IT_0038*conj(U_st_01);
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = IT_0157*IT_0203;
    const ccomplex_t IT_0205 = IT_0155*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_cd*e_em*IT_0005*IT_0036*IT_0038;
    const ccomplex_t IT_0208 = (-0.5)*IT_0207;
    const ccomplex_t IT_0209 = IT_0151*IT_0208;
    const ccomplex_t IT_0210 = IT_0147*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = IT_0196 + IT_0201 + IT_0206 + IT_0211;
    const ccomplex_t IT_0213 = s_13 + s_14;
    const ccomplex_t IT_0214 = IT_0030*IT_0213;
    const ccomplex_t IT_0215 = IT_0037*IT_0088;
    const ccomplex_t IT_0216 = IT_0090*IT_0112;
    const ccomplex_t IT_0217 = IT_0005*(IT_0215 + -IT_0216);
    const ccomplex_t IT_0218 = 1.4142135623731*e_em*IT_0217;
    const ccomplex_t IT_0219 = 0.5*IT_0218;
    const ccomplex_t IT_0220 = IT_0131*IT_0219;
    const ccomplex_t IT_0221 = IT_0124*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = IT_0037*IT_0090;
    const ccomplex_t IT_0224 = IT_0088*IT_0112;
    const ccomplex_t IT_0225 = IT_0005*(IT_0223 + IT_0224);
    const ccomplex_t IT_0226 = 1.4142135623731*e_em*IT_0225;
    const ccomplex_t IT_0227 = (-0.5)*IT_0226;
    const ccomplex_t IT_0228 = IT_0114*IT_0227;
    const ccomplex_t IT_0229 = IT_0111*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = -IT_0222 + -IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u2
      *e_em*IT_0005*IT_0036*IT_0158;
    const ccomplex_t IT_0233 = 0.5*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u2)*e_em*IT_0005*IT_0036*IT_0158;
    const ccomplex_t IT_0235 = (-0.5)*IT_0234;
    const ccomplex_t IT_0236 = IT_0233*IT_0235;
    const ccomplex_t IT_0237 = cpow((-2)*s_13 + IT_0018 + IT_0027 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0238 = IT_0236*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0005*IT_0036*IT_0158;
    const ccomplex_t IT_0241 = 0.5*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *V_ud*e_em*IT_0005*IT_0036*IT_0158;
    const ccomplex_t IT_0243 = (-0.5)*IT_0242;
    const ccomplex_t IT_0244 = IT_0241*IT_0243;
    const ccomplex_t IT_0245 = cpow((-2)*s_13 + IT_0018 + IT_0027 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0246 = IT_0244*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*conj(V_cd)*e_em*conj(V_Wp2)
      *IT_0005;
    const ccomplex_t IT_0249 = IT_0149*IT_0248;
    const ccomplex_t IT_0250 = IT_0147*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0253 = IT_0143*IT_0252;
    const ccomplex_t IT_0254 = IT_0139*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = conj(V_td)*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0257 = IT_0005*IT_0256;
    const ccomplex_t IT_0258 = m_t*conj(V_td)*conj(V_u2)*e_em*IT_0036*U_st_10;
    const ccomplex_t IT_0259 = IT_0159*IT_0258;
    const ccomplex_t IT_0260 = 1.4142135623731*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*(IT_0257 + (-0.5)*IT_0260);
    const ccomplex_t IT_0262 = -IT_0261;
    const ccomplex_t IT_0263 = IT_0187*IT_0262;
    const ccomplex_t IT_0264 = IT_0174*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = conj(V_td)*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0267 = IT_0005*IT_0266;
    const ccomplex_t IT_0268 = m_t*conj(V_td)*conj(V_u2)*e_em*IT_0036*U_st_11;
    const ccomplex_t IT_0269 = IT_0159*IT_0268;
    const ccomplex_t IT_0270 = 1.4142135623731*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*(IT_0267 + (-0.5)*IT_0270);
    const ccomplex_t IT_0272 = -IT_0271;
    const ccomplex_t IT_0273 = IT_0170*IT_0272;
    const ccomplex_t IT_0274 = IT_0155*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = IT_0239 + IT_0247 + IT_0251 + IT_0255 + IT_0265
       + IT_0275;
    const ccomplex_t IT_0277 = m_d*s_14*m_C_2;
    const ccomplex_t IT_0278 = IT_0208*IT_0248;
    const ccomplex_t IT_0279 = IT_0147*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*IT_0279;
    const ccomplex_t IT_0281 = IT_0193*IT_0252;
    const ccomplex_t IT_0282 = IT_0139*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = IT_0203*IT_0272;
    const ccomplex_t IT_0285 = IT_0155*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = IT_0198*IT_0262;
    const ccomplex_t IT_0288 = IT_0174*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = -IT_0280 + -IT_0283 + IT_0286 + IT_0289;
    const ccomplex_t IT_0291 = 12*IT_0021;
    const ccomplex_t IT_0292 = 24*conj(IT_0021);
    const ccomplex_t IT_0293 = m_d*s_24*m_C_2;
    const ccomplex_t IT_0294 = s_13*s_24;
    const ccomplex_t IT_0295 = m_d*s_13*m_C_2;
    const ccomplex_t IT_0296 = 24*IT_0026;
    const ccomplex_t IT_0297 = 48*IT_0026;
    const ccomplex_t IT_0298 = 24*IT_0079;
    const ccomplex_t IT_0299 = 24*IT_0294;
    const ccomplex_t IT_0300 = 24*IT_0028;
    const ccomplex_t IT_0301 = 48*IT_0294;
    const ccomplex_t IT_0302 = 3*IT_0034;
    const ccomplex_t IT_0303 = conj(IT_0108)*IT_0302;
    const ccomplex_t IT_0304 = 3*IT_0087;
    const ccomplex_t IT_0305 = conj(IT_0066)*IT_0304;
    const ccomplex_t IT_0306 = 12*IT_0214;
    const ccomplex_t IT_0307 = 24*IT_0293;
    const ccomplex_t IT_0308 = 24*IT_0295;
    const ccomplex_t IT_0309 = (-12)*IT_0110;
    const ccomplex_t IT_0310 = 48*conj(IT_0026)*IT_0079 + IT_0067*IT_0292 +
       conj(IT_0212)*IT_0298 + conj(IT_0276)*IT_0299 + conj(IT_0083)*IT_0300 +
       conj(IT_0078)*IT_0301 + IT_0303 + IT_0305 + conj(IT_0135)*IT_0306 + conj
      (IT_0290)*IT_0307 + conj(IT_0191)*IT_0308 + conj(IT_0231)*IT_0309;
    const ccomplex_t IT_0311 = s_12*s_34;
    const ccomplex_t IT_0312 = -IT_0029;
    const ccomplex_t IT_0313 = IT_0028 + IT_0312;
    const ccomplex_t IT_0314 = IT_0294 + IT_0311 + IT_0313;
    const ccomplex_t IT_0315 = 6*IT_0110;
    const ccomplex_t IT_0316 = 6*IT_0214;
    const ccomplex_t IT_0317 = s_34 + IT_0027;
    const ccomplex_t IT_0318 = IT_0018*IT_0317;
    const ccomplex_t IT_0319 = 6*IT_0318;
    const ccomplex_t IT_0320 = s_12*IT_0317;
    const ccomplex_t IT_0321 = 24*IT_0320;
    const ccomplex_t IT_0322 = (-24)*IT_0318;
    const ccomplex_t IT_0323 = (-12)*IT_0214;
    const ccomplex_t IT_0324 = ((-12)*conj(IT_0021) + 12*conj(IT_0026))
      *IT_0110 + conj(IT_0078)*IT_0306 + 6*conj(IT_0191)*IT_0314 + conj(IT_0212)
      *IT_0315 + conj(IT_0276)*IT_0316 + conj(IT_0290)*IT_0319 + conj(IT_0135)
      *IT_0321 + conj(IT_0231)*IT_0322 + conj(IT_0083)*IT_0323;
    const ccomplex_t IT_0325 = 12*IT_0067;
    const ccomplex_t IT_0326 = 12*IT_0028;
    const ccomplex_t IT_0327 = 48*IT_0029;
    const ccomplex_t IT_0328 = 12*IT_0277;
    const ccomplex_t IT_0329 = 12*IT_0138;
    const ccomplex_t IT_0330 = 12*IT_0110;
    const ccomplex_t IT_0331 = 24*conj(IT_0026);
    const ccomplex_t IT_0332 = 48*conj(IT_0021)*IT_0079 + conj(IT_0078)
      *IT_0300 + IT_0303 + IT_0305 + conj(IT_0135)*IT_0323 + conj(IT_0212)
      *IT_0325 + conj(IT_0276)*IT_0326 + conj(IT_0083)*IT_0327 + conj(IT_0191)
      *IT_0328 + conj(IT_0290)*IT_0329 + conj(IT_0231)*IT_0330 + IT_0067*IT_0331;
    const ccomplex_t IT_0333 = -IT_0311;
    const ccomplex_t IT_0334 = -IT_0294;
    const ccomplex_t IT_0335 = IT_0028 + IT_0029 + IT_0333 + IT_0334;
    const ccomplex_t IT_0336 = 6*IT_0335;
    const ccomplex_t IT_0337 = 12*conj(IT_0021);
    const ccomplex_t IT_0338 = 12*IT_0293;
    const ccomplex_t IT_0339 = 12*IT_0295;
    const ccomplex_t IT_0340 = 12*IT_0079;
    const ccomplex_t IT_0341 = 12*IT_0294;
    const ccomplex_t IT_0342 = -s_34;
    const ccomplex_t IT_0343 = IT_0027 + IT_0342;
    const ccomplex_t IT_0344 = IT_0018*IT_0343;
    const ccomplex_t IT_0345 = (-6)*IT_0344;
    const ccomplex_t IT_0346 = (-0.166666666666667)*conj(IT_0078)*IT_0307 +
       conj(IT_0231)*IT_0314 + (-0.166666666666667)*conj(IT_0135)*IT_0319 + (
      -0.166666666666667)*conj(IT_0083)*IT_0329 + (-0.166666666666667)*IT_0295
      *IT_0331 + (-0.166666666666667)*conj(IT_0108)*IT_0336 + (
      -0.166666666666667)*IT_0277*IT_0337 + (-0.166666666666667)*conj(IT_0276)
      *IT_0338 + (-0.166666666666667)*conj(IT_0212)*IT_0339 + (
      -0.166666666666667)*conj(IT_0191)*IT_0340 + (-0.166666666666667)*conj
      (IT_0290)*IT_0341 + (-0.166666666666667)*conj(IT_0066)*IT_0345;
    const ccomplex_t IT_0347 = (-6)*IT_0318;
    const ccomplex_t IT_0348 = 0.166666666666667*conj(IT_0078)*IT_0308 + conj
      (IT_0135)*IT_0314 + 0.166666666666667*conj(IT_0083)*IT_0328 +
       0.166666666666667*IT_0293*IT_0331 + 0.166666666666667*conj(IT_0066)
      *IT_0336 + 0.166666666666667*IT_0138*IT_0337 + 0.166666666666667*conj
      (IT_0212)*IT_0338 + 0.166666666666667*conj(IT_0276)*IT_0339 +
       0.166666666666667*conj(IT_0290)*IT_0340 + 0.166666666666667*conj(IT_0191)
      *IT_0341 + 0.166666666666667*conj(IT_0108)*IT_0345 + 0.166666666666667
      *conj(IT_0231)*IT_0347;
    const ccomplex_t IT_0349 = (-6)*IT_0214;
    const ccomplex_t IT_0350 = (-6)*IT_0110;
    const ccomplex_t IT_0351 = (-0.166666666666667)*conj(IT_0078)*IT_0309 +
       conj(IT_0290)*IT_0314 + (-0.166666666666667)*conj(IT_0231)*IT_0321 + (
      -0.166666666666667)*conj(IT_0135)*IT_0322 + (-0.166666666666667)*conj
      (IT_0083)*IT_0330 + (-0.166666666666667)*IT_0214*((-12)*conj(IT_0026) +
       IT_0337) + (-0.166666666666667)*conj(IT_0191)*IT_0347 + (
      -0.166666666666667)*conj(IT_0212)*IT_0349 + (-0.166666666666667)*conj
      (IT_0276)*IT_0350;
    const ccomplex_t IT_0352 = (-2)*IT_0033;
    const ccomplex_t IT_0353 = (-3)*IT_0352;
    const ccomplex_t IT_0354 = 2*IT_0086;
    const ccomplex_t IT_0355 = (-3)*IT_0354;
    const ccomplex_t IT_0356 = conj(IT_0078)*IT_0299 + conj(IT_0135)*IT_0316 +
       conj(IT_0083)*IT_0326 + IT_0079*IT_0331 + IT_0067*IT_0337 + conj(IT_0290)
      *IT_0338 + conj(IT_0191)*IT_0339 + conj(IT_0212)*IT_0340 + conj(IT_0276)
      *IT_0341 + conj(IT_0231)*IT_0350 + conj(IT_0108)*IT_0353 + conj(IT_0066)
      *IT_0355;
    const ccomplex_t IT_0357 = conj(IT_0078)*IT_0298 + conj(IT_0135)*IT_0315 +
       conj(IT_0083)*IT_0325 + IT_0294*IT_0331 + IT_0028*IT_0337 + conj(IT_0191)
      *IT_0338 + conj(IT_0290)*IT_0339 + conj(IT_0276)*IT_0340 + conj(IT_0212)
      *IT_0341 + conj(IT_0231)*IT_0349 + conj(IT_0066)*IT_0353 + conj(IT_0108)
      *IT_0355;
    const ccomplex_t IT_0358 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0359 = (-0.333333333333333)*IT_0358;
    const ccomplex_t IT_0360 = IT_0358*IT_0359;
    const ccomplex_t IT_0361 = cpow(s_12 + IT_0018 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0362 = IT_0360*IT_0361;
    const ccomplex_t IT_0363 = (0 + _Complex_I*1)*IT_0362;
    const ccomplex_t IT_0364 = 0.5*IT_0363;
    const ccomplex_t IT_0365 = IT_0298 + IT_0299 + IT_0325 + IT_0326;
    const ccomplex_t IT_0366 = 24*IT_0067;
    const ccomplex_t IT_0367 = 48*IT_0079;
    const ccomplex_t IT_0368 = IT_0300 + IT_0327 + IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = IT_0300 + IT_0301 + IT_0366 + IT_0367;
    const ccomplex_t IT_0370 = -s_14 + -s_23;
    const ccomplex_t IT_0371 = s_13 + s_24;
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = IT_0030*IT_0372;
    const ccomplex_t IT_0374 = (-8)*IT_0373;
    const ccomplex_t IT_0375 = 3*IT_0374;
    const ccomplex_t IT_0376 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0377 = IT_0371 + IT_0376;
    const ccomplex_t IT_0378 = IT_0030*IT_0377;
    const ccomplex_t IT_0379 = 24*IT_0378;
    const ccomplex_t IT_0380 = (-24)*IT_0378;
    const ccomplex_t IT_0381 = (-3)*IT_0374;
    const ccomplex_t IT_0382 = (96*IT_0028 + 96*IT_0029 + 96*IT_0067 + 192
      *IT_0079 + 96*IT_0294)*IT_0364 + (IT_0212 + -conj(IT_0212))*IT_0365 + 
      (IT_0276 + -conj(IT_0276))*IT_0365 + (IT_0021 + -conj(IT_0021) + IT_0083 +
       -conj(IT_0083))*IT_0368 + (IT_0026 + -conj(IT_0026) + IT_0078 + -conj
      (IT_0078))*IT_0369 + -(conj(IT_0066) + conj(IT_0108))*IT_0375 + -(conj
      (IT_0191) + conj(IT_0290))*IT_0379 + -(IT_0191 + IT_0290)*IT_0380 + -
      (IT_0066 + IT_0108)*IT_0381;
    const ccomplex_t IT_0383 = s_12*IT_0343;
    const ccomplex_t IT_0384 = (-24)*IT_0383;
    const ccomplex_t IT_0385 = conj(IT_0078) + conj(IT_0083);
    const ccomplex_t IT_0386 = 0.333333333333333*IT_0385;
    const ccomplex_t IT_0387 = 24*IT_0344;
    const ccomplex_t IT_0388 = 0.333333333333333*conj(IT_0108);
    const ccomplex_t IT_0389 = 0.333333333333333*conj(IT_0191);
    const ccomplex_t IT_0390 = 0.333333333333333*conj(IT_0212);
    const ccomplex_t IT_0391 = 0.333333333333333*conj(IT_0276);
    const ccomplex_t IT_0392 = 0.333333333333333*conj(IT_0290);
    const ccomplex_t IT_0393 = 3*conj(IT_0026);
    const ccomplex_t IT_0394 = 0.333333333333333*IT_0393;
    const ccomplex_t IT_0395 = conj(IT_0066)*IT_0384 + 3*IT_0304*IT_0386 + 3
      *IT_0387*IT_0388 + 3*IT_0336*IT_0389 + 3*IT_0353*IT_0390 + 3*IT_0355
      *IT_0391 + 3*IT_0345*IT_0392 + 3*IT_0034*(conj(IT_0021) + IT_0394);
    const ccomplex_t IT_0396 = 3*IT_0302*IT_0386 + conj(IT_0066)*IT_0387 + 3
      *IT_0384*IT_0388 + 3*IT_0345*IT_0389 + 3*IT_0355*IT_0390 + 3*IT_0353
      *IT_0391 + 3*IT_0336*IT_0392 + 3*IT_0087*(conj(IT_0021) + IT_0394);
    const ccomplex_t IT_0397 = IT_0021*IT_0136 + (IT_0034*conj(IT_0066) +
       IT_0087*conj(IT_0108))*IT_0137 + (IT_0138*conj(IT_0191) + IT_0028*conj
      (IT_0212) + IT_0214*conj(IT_0231) + IT_0067*conj(IT_0276) + IT_0277*conj
      (IT_0290))*IT_0291 + IT_0026*(12*IT_0110*conj(IT_0135) + (-12)*IT_0214
      *conj(IT_0231) + IT_0028*IT_0292) + (IT_0067*conj(IT_0083) + IT_0079*conj
      (IT_0276) + conj(IT_0191)*IT_0293 + conj(IT_0212)*IT_0294 + conj(IT_0290)
      *IT_0295)*IT_0296 + (conj(IT_0078)*IT_0079 + conj(IT_0026)*IT_0294)
      *IT_0297 + IT_0078*IT_0310 + IT_0135*IT_0324 + IT_0083*IT_0332 + (-6)
      *IT_0290*IT_0346 + 6*IT_0191*IT_0348 + (-6)*IT_0231*IT_0351 + IT_0276
      *IT_0356 + IT_0212*IT_0357 + -IT_0364*IT_0382 + IT_0066*IT_0395 + IT_0108
      *IT_0396;
    return create_ccomplex_return(IT_0397);
}

