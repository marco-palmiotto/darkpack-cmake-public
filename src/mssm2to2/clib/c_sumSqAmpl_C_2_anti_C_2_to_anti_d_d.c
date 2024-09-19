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
    const ccomplex_t IT_0000 = pow(m_d, 2);
    const ccomplex_t IT_0001 = pow(m_C_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*V_cd*e_em*V_Wp2*IT_0004;
    const ccomplex_t IT_0006 = -IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_cd)*e_em*IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0011;
    const ccomplex_t IT_0013 = IT_0002*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0016;
    const ccomplex_t IT_0018 = 0.333333333333333*IT_0017;
    const ccomplex_t IT_0019 = IT_0003*IT_0016;
    const ccomplex_t IT_0020 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0004*IT_0015;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0021 + -IT_0023 + (-2)
      *IT_0025);
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0018*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0001 + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = m_d*s_24*m_C_2;
    const ccomplex_t IT_0033 = -s_34;
    const ccomplex_t IT_0034 = IT_0000 + IT_0033;
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + (-2)*IT_0001 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0037 = cos(alpha);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_d*e_em*IT_0004*IT_0007
      *IT_0009*IT_0037;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = sin(alpha);
    const ccomplex_t IT_0043 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = IT_0004*(IT_0041 + IT_0044);
    const ccomplex_t IT_0046 = 1.4142135623731*e_em*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0039*IT_0048;
    const ccomplex_t IT_0050 = IT_0036*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_d*e_em*IT_0004*IT_0007
      *IT_0009*IT_0042;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0037*IT_0043;
    const ccomplex_t IT_0055 = IT_0040*IT_0042;
    const ccomplex_t IT_0056 = IT_0004*(IT_0054 + -IT_0055);
    const ccomplex_t IT_0057 = 1.4142135623731*e_em*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0053*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0001 + -reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = -IT_0051 + -IT_0063;
    const ccomplex_t IT_0065 = s_34 + IT_0000;
    const ccomplex_t IT_0066 = IT_0001*IT_0065;
    const ccomplex_t IT_0067 = IT_0008*IT_0040;
    const ccomplex_t IT_0068 = sin(beta);
    const ccomplex_t IT_0069 = IT_0043*IT_0068;
    const ccomplex_t IT_0070 = IT_0004*(IT_0067 + -IT_0069);
    const ccomplex_t IT_0071 = 1.4142135623731*e_em*IT_0070;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = m_d*e_em*IT_0004*IT_0007;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = pow(m_Z, 2);
    const ccomplex_t IT_0077 = cpow((-2)*s_12 + (-2)*IT_0001 + IT_0076 + 
      -reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0075*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = m_d*e_em*IT_0004*IT_0007*IT_0009*IT_0068;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = IT_0008*IT_0043;
    const ccomplex_t IT_0083 = IT_0040*IT_0068;
    const ccomplex_t IT_0084 = IT_0004*(IT_0082 + IT_0083);
    const ccomplex_t IT_0085 = 1.4142135623731*e_em*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0081*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_12 + (-2)*IT_0001 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = -IT_0079 + -IT_0090;
    const ccomplex_t IT_0092 = m_d*s_13*m_C_2;
    const ccomplex_t IT_0093 = e_em*IT_0019;
    const ccomplex_t IT_0094 = e_em*IT_0022;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0093 + 3*IT_0094);
    const ccomplex_t IT_0096 = (-0.166666666666667)*IT_0095;
    const ccomplex_t IT_0097 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0098 = IT_0019*IT_0097;
    const ccomplex_t IT_0099 = IT_0022*IT_0097;
    const ccomplex_t IT_0100 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0101 = IT_0022*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + -IT_0099 + (-2)
      *IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0096*IT_0103;
    const ccomplex_t IT_0105 = IT_0029*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = V_u2*U_Wm2;
    const ccomplex_t IT_0109 = IT_0037*IT_0108;
    const ccomplex_t IT_0110 = U_d2*V_Wp2;
    const ccomplex_t IT_0111 = IT_0042*IT_0110;
    const ccomplex_t IT_0112 = IT_0004*(IT_0109 + -IT_0111);
    const ccomplex_t IT_0113 = 1.4142135623731*e_em*IT_0112;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = IT_0053*IT_0115;
    const ccomplex_t IT_0117 = IT_0061*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0037*IT_0110;
    const ccomplex_t IT_0120 = IT_0042*IT_0108;
    const ccomplex_t IT_0121 = IT_0004*(IT_0119 + IT_0120);
    const ccomplex_t IT_0122 = 1.4142135623731*e_em*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = IT_0039*IT_0124;
    const ccomplex_t IT_0126 = IT_0036*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = -IT_0118 + -IT_0127;
    const ccomplex_t IT_0129 = s_12*IT_0000;
    const ccomplex_t IT_0130 = s_14*s_23;
    const ccomplex_t IT_0131 = s_12*s_34;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = s_13*s_24;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = IT_0129 + IT_0130 + IT_0132 + IT_0134;
    const ccomplex_t IT_0136 = IT_0008*IT_0108;
    const ccomplex_t IT_0137 = IT_0068*IT_0110;
    const ccomplex_t IT_0138 = IT_0004*(IT_0136 + IT_0137);
    const ccomplex_t IT_0139 = 1.4142135623731*e_em*IT_0138;
    const ccomplex_t IT_0140 = 0.5*IT_0139;
    const ccomplex_t IT_0141 = IT_0081*IT_0140;
    const ccomplex_t IT_0142 = IT_0088*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0008*IT_0110;
    const ccomplex_t IT_0145 = IT_0068*IT_0108;
    const ccomplex_t IT_0146 = IT_0004*(IT_0144 + -IT_0145);
    const ccomplex_t IT_0147 = 1.4142135623731*e_em*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0074*IT_0148;
    const ccomplex_t IT_0150 = IT_0077*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = -IT_0143 + -IT_0151;
    const ccomplex_t IT_0153 = -IT_0130;
    const ccomplex_t IT_0154 = IT_0129 + IT_0131 + IT_0133 + IT_0153;
    const ccomplex_t IT_0155 = m_d*m_C_2;
    const ccomplex_t IT_0156 = -s_14;
    const ccomplex_t IT_0157 = s_13 + IT_0156;
    const ccomplex_t IT_0158 = IT_0155*IT_0157;
    const ccomplex_t IT_0159 = (-4)*IT_0158;
    const ccomplex_t IT_0160 = -s_24;
    const ccomplex_t IT_0161 = s_23 + IT_0160;
    const ccomplex_t IT_0162 = IT_0155*IT_0161;
    const ccomplex_t IT_0163 = 4*IT_0162;
    const ccomplex_t IT_0164 = 3*IT_0031;
    const ccomplex_t IT_0165 = m_d*s_23*m_C_2;
    const ccomplex_t IT_0166 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0167 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0004*IT_0007*IT_0009*U_st_01;
    const ccomplex_t IT_0168 = 0.5*IT_0167;
    const ccomplex_t IT_0169 = cpow(IT_0068, -1);
    const ccomplex_t IT_0170 = IT_0004*IT_0169;
    const ccomplex_t IT_0171 = m_t*V_td*V_u2*e_em*IT_0007*conj(U_st_11);
    const ccomplex_t IT_0172 = IT_0170*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = (-0.5)*IT_0174;
    const ccomplex_t IT_0176 = V_td*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0177 = IT_0004*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0175 + IT_0180;
    const ccomplex_t IT_0182 = IT_0168*IT_0181;
    const ccomplex_t IT_0183 = IT_0166*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0004;
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*V_ud
      *e_em*IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0189 = 0.5*IT_0188;
    const ccomplex_t IT_0190 = IT_0187*IT_0189;
    const ccomplex_t IT_0191 = IT_0185*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0194 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0004*IT_0007*IT_0009*U_st_00;
    const ccomplex_t IT_0195 = 0.5*IT_0194;
    const ccomplex_t IT_0196 = m_t*V_td*V_u2*e_em*IT_0007*conj(U_st_10);
    const ccomplex_t IT_0197 = IT_0170*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = V_td*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0202 = IT_0004*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0200 + IT_0205;
    const ccomplex_t IT_0207 = IT_0195*IT_0206;
    const ccomplex_t IT_0208 = IT_0193*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = IT_0184 + -IT_0192 + IT_0209;
    const ccomplex_t IT_0211 = s_34*IT_0001;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0004*IT_0007*IT_0169;
    const ccomplex_t IT_0213 = 0.5*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *V_ud*e_em*IT_0004*IT_0007*IT_0169;
    const ccomplex_t IT_0215 = (-0.5)*IT_0214;
    const ccomplex_t IT_0216 = IT_0213*IT_0215;
    const ccomplex_t IT_0217 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0218 = IT_0216*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp2)*IT_0004;
    const ccomplex_t IT_0221 = IT_0187*IT_0220;
    const ccomplex_t IT_0222 = IT_0185*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = conj(V_td)*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0225 = IT_0004*IT_0224;
    const ccomplex_t IT_0226 = m_t*conj(V_td)*conj(V_u2)*e_em*IT_0007*U_st_10;
    const ccomplex_t IT_0227 = IT_0170*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*(IT_0225 + (-0.5)*IT_0228);
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0206*IT_0230;
    const ccomplex_t IT_0232 = IT_0193*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = conj(V_td)*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0235 = IT_0004*IT_0234;
    const ccomplex_t IT_0236 = m_t*conj(V_td)*conj(V_u2)*e_em*IT_0007*U_st_11;
    const ccomplex_t IT_0237 = IT_0170*IT_0236;
    const ccomplex_t IT_0238 = 1.4142135623731*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*(IT_0235 + (-0.5)*IT_0238);
    const ccomplex_t IT_0240 = -IT_0239;
    const ccomplex_t IT_0241 = IT_0181*IT_0240;
    const ccomplex_t IT_0242 = IT_0166*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = IT_0219 + IT_0223 + IT_0233 + IT_0243;
    const ccomplex_t IT_0245 = s_13 + s_14;
    const ccomplex_t IT_0246 = IT_0155*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_cd*e_em*IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0248 = (-0.5)*IT_0247;
    const ccomplex_t IT_0249 = IT_0011*IT_0248;
    const ccomplex_t IT_0250 = IT_0002*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0004*IT_0007*IT_0009*conj(U_st_01);
    const ccomplex_t IT_0253 = (-0.5)*IT_0252;
    const ccomplex_t IT_0254 = IT_0168*IT_0253;
    const ccomplex_t IT_0255 = IT_0166*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0004*IT_0007*IT_0009*conj(U_st_00);
    const ccomplex_t IT_0258 = (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = IT_0195*IT_0258;
    const ccomplex_t IT_0260 = IT_0193*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = IT_0251 + IT_0256 + IT_0261;
    const ccomplex_t IT_0263 = m_d*s_14*m_C_2;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0265 = (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = IT_0220*IT_0265;
    const ccomplex_t IT_0267 = IT_0185*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*conj(V_cd)*e_em*conj(V_Wp2)
      *IT_0004;
    const ccomplex_t IT_0270 = IT_0248*IT_0269;
    const ccomplex_t IT_0271 = IT_0002*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = IT_0240*IT_0253;
    const ccomplex_t IT_0274 = IT_0166*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = -IT_0268 + -IT_0272 + IT_0275;
    const ccomplex_t IT_0277 = IT_0027*IT_0096;
    const ccomplex_t IT_0278 = IT_0029*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = 12*IT_0279;
    const ccomplex_t IT_0281 = 24*IT_0279;
    const ccomplex_t IT_0282 = s_23 + s_24;
    const ccomplex_t IT_0283 = IT_0155*IT_0282;
    const ccomplex_t IT_0284 = 24*conj(IT_0279);
    const ccomplex_t IT_0285 = (-24)*conj(IT_0014);
    const ccomplex_t IT_0286 = (-12)*IT_0014;
    const ccomplex_t IT_0287 = IT_0018*IT_0103;
    const ccomplex_t IT_0288 = IT_0029*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = -IT_0289;
    const ccomplex_t IT_0291 = 24*IT_0031;
    const ccomplex_t IT_0292 = 12*IT_0014;
    const ccomplex_t IT_0293 = 48*IT_0031;
    const ccomplex_t IT_0294 = IT_0189*IT_0265;
    const ccomplex_t IT_0295 = IT_0185*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = -IT_0296;
    const ccomplex_t IT_0298 = (-24)*IT_0031;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u2
      *e_em*IT_0004*IT_0007*IT_0169;
    const ccomplex_t IT_0300 = 0.5*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u2)*e_em*IT_0004*IT_0007*IT_0169;
    const ccomplex_t IT_0302 = (-0.5)*IT_0301;
    const ccomplex_t IT_0303 = IT_0300*IT_0302;
    const ccomplex_t IT_0304 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0305 = IT_0303*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = IT_0006*IT_0269;
    const ccomplex_t IT_0308 = IT_0002*IT_0307;
    const ccomplex_t IT_0309 = (0 + _Complex_I*1)*IT_0308;
    const ccomplex_t IT_0310 = -IT_0306 + -IT_0309;
    const ccomplex_t IT_0311 = IT_0230*IT_0258;
    const ccomplex_t IT_0312 = IT_0193*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*IT_0312;
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = IT_0000*IT_0001;
    const ccomplex_t IT_0316 = (-12)*IT_0279;
    const ccomplex_t IT_0317 = 24*IT_0315;
    const ccomplex_t IT_0318 = 24*IT_0133;
    const ccomplex_t IT_0319 = 24*IT_0129;
    const ccomplex_t IT_0320 = 48*IT_0133;
    const ccomplex_t IT_0321 = (-24)*IT_0315;
    const ccomplex_t IT_0322 = (-24)*IT_0133;
    const ccomplex_t IT_0323 = 3*IT_0163;
    const ccomplex_t IT_0324 = conj(IT_0064)*IT_0323;
    const ccomplex_t IT_0325 = 3*IT_0159;
    const ccomplex_t IT_0326 = conj(IT_0128)*IT_0325;
    const ccomplex_t IT_0327 = 12*IT_0246;
    const ccomplex_t IT_0328 = 24*IT_0092;
    const ccomplex_t IT_0329 = 24*IT_0032;
    const ccomplex_t IT_0330 = (-24)*IT_0032;
    const ccomplex_t IT_0331 = (-12)*IT_0283;
    const ccomplex_t IT_0332 = IT_0211*IT_0284 + IT_0092*IT_0285 + 48*conj
      (IT_0031)*IT_0315 + conj(IT_0262)*IT_0317 + conj(IT_0244)*IT_0318 + conj
      (IT_0290)*IT_0319 + conj(IT_0107)*IT_0320 + conj(IT_0297)*IT_0321 + conj
      (IT_0310)*IT_0322 + IT_0324 + IT_0326 + conj(IT_0152)*IT_0327 + conj
      (IT_0210)*IT_0328 + conj(IT_0276)*IT_0329 + conj(IT_0314)*IT_0330 + conj
      (IT_0091)*IT_0331;
    const ccomplex_t IT_0333 = 6*IT_0135;
    const ccomplex_t IT_0334 = 6*IT_0066;
    const ccomplex_t IT_0335 = 12*conj(IT_0279);
    const ccomplex_t IT_0336 = 12*IT_0133;
    const ccomplex_t IT_0337 = 12*IT_0092;
    const ccomplex_t IT_0338 = 12*IT_0032;
    const ccomplex_t IT_0339 = 12*IT_0315;
    const ccomplex_t IT_0340 = 12*IT_0165;
    const ccomplex_t IT_0341 = 24*conj(IT_0031);
    const ccomplex_t IT_0342 = (-12)*conj(IT_0014);
    const ccomplex_t IT_0343 = (-12)*IT_0133;
    const ccomplex_t IT_0344 = (-12)*IT_0092;
    const ccomplex_t IT_0345 = (-12)*IT_0032;
    const ccomplex_t IT_0346 = (-6)*IT_0154;
    const ccomplex_t IT_0347 = (-6)*IT_0035;
    const ccomplex_t IT_0348 = conj(IT_0107)*IT_0329 + conj(IT_0064)*IT_0333 +
       conj(IT_0152)*IT_0334 + IT_0263*IT_0335 + conj(IT_0276)*IT_0336 + conj
      (IT_0262)*IT_0337 + conj(IT_0244)*IT_0338 + conj(IT_0210)*IT_0339 + conj
      (IT_0290)*IT_0340 + IT_0092*IT_0341 + IT_0315*IT_0342 + conj(IT_0314)
      *IT_0343 + conj(IT_0297)*IT_0344 + conj(IT_0310)*IT_0345 + conj(IT_0091)
      *IT_0346 + conj(IT_0128)*IT_0347;
    const ccomplex_t IT_0349 = 6*IT_0154;
    const ccomplex_t IT_0350 = 6*IT_0035;
    const ccomplex_t IT_0351 = 12*conj(IT_0014);
    const ccomplex_t IT_0352 = (-24)*conj(IT_0031);
    const ccomplex_t IT_0353 = (-12)*conj(IT_0279);
    const ccomplex_t IT_0354 = (-12)*IT_0315;
    const ccomplex_t IT_0355 = (-12)*IT_0165;
    const ccomplex_t IT_0356 = (-6)*IT_0135;
    const ccomplex_t IT_0357 = (-6)*IT_0066;
    const ccomplex_t IT_0358 = conj(IT_0107)*IT_0330 + conj(IT_0314)*IT_0336 +
       conj(IT_0297)*IT_0337 + conj(IT_0310)*IT_0338 + conj(IT_0276)*IT_0343 +
       conj(IT_0262)*IT_0344 + conj(IT_0244)*IT_0345 + conj(IT_0091)*IT_0349 +
       conj(IT_0128)*IT_0350 + IT_0315*IT_0351 + IT_0092*IT_0352 + IT_0263
      *IT_0353 + conj(IT_0210)*IT_0354 + conj(IT_0290)*IT_0355 + conj(IT_0064)
      *IT_0356 + conj(IT_0152)*IT_0357;
    const ccomplex_t IT_0359 = 12*IT_0263;
    const ccomplex_t IT_0360 = 12*IT_0133*conj(IT_0210) + conj(IT_0107)
      *IT_0328 + conj(IT_0128)*IT_0333 + IT_0165*IT_0335 + conj(IT_0244)*IT_0337
       + conj(IT_0262)*IT_0338 + conj(IT_0276)*IT_0339 + IT_0032*IT_0341 +
       IT_0133*IT_0342 + conj(IT_0310)*IT_0344 + conj(IT_0297)*IT_0345 + conj
      (IT_0064)*IT_0347 + conj(IT_0152)*IT_0349 + conj(IT_0314)*IT_0354 + conj
      (IT_0091)*IT_0357 + conj(IT_0290)*IT_0359;
    const ccomplex_t IT_0361 = 6*IT_0246;
    const ccomplex_t IT_0362 = 6*IT_0283;
    const ccomplex_t IT_0363 = s_12*IT_0065;
    const ccomplex_t IT_0364 = 24*IT_0363;
    const ccomplex_t IT_0365 = (-24)*IT_0066;
    const ccomplex_t IT_0366 = (-12)*IT_0246;
    const ccomplex_t IT_0367 = (-6)*conj(IT_0014);
    const ccomplex_t IT_0368 = (-6)*IT_0246;
    const ccomplex_t IT_0369 = (-6)*IT_0283;
    const ccomplex_t IT_0370 = conj(IT_0107)*IT_0327 + conj(IT_0276)*IT_0334 +
       conj(IT_0210)*IT_0349 + IT_0283*(12*conj(IT_0031) + IT_0353) + conj
      (IT_0314)*IT_0357 + conj(IT_0244)*IT_0361 + conj(IT_0262)*IT_0362 + conj
      (IT_0152)*IT_0364 + conj(IT_0091)*IT_0365 + conj(IT_0290)*IT_0366 +
       IT_0154*IT_0367 + conj(IT_0310)*IT_0368 + conj(IT_0297)*IT_0369;
    const ccomplex_t IT_0371 = 6*conj(IT_0014);
    const ccomplex_t IT_0372 = 12*IT_0283;
    const ccomplex_t IT_0373 = conj(IT_0107)*IT_0331 + IT_0246*((-12)*conj
      (IT_0031) + IT_0335) + conj(IT_0276)*IT_0346 + conj(IT_0314)*IT_0349 +
       conj(IT_0210)*IT_0357 + conj(IT_0297)*IT_0361 + conj(IT_0310)*IT_0362 +
       conj(IT_0091)*IT_0364 + conj(IT_0152)*IT_0365 + conj(IT_0262)*IT_0368 +
       conj(IT_0244)*IT_0369 + IT_0066*IT_0371 + conj(IT_0290)*IT_0372;
    const ccomplex_t IT_0374 = 12*IT_0129;
    const ccomplex_t IT_0375 = 2*IT_0158;
    const ccomplex_t IT_0376 = (-3)*IT_0375;
    const ccomplex_t IT_0377 = (-2)*IT_0162;
    const ccomplex_t IT_0378 = (-3)*IT_0377;
    const ccomplex_t IT_0379 = conj(IT_0107)*IT_0318 + IT_0211*IT_0335 + conj
      (IT_0244)*IT_0336 + conj(IT_0210)*IT_0337 + conj(IT_0276)*IT_0338 + conj
      (IT_0262)*IT_0339 + IT_0315*IT_0341 + IT_0092*IT_0342 + conj(IT_0310)
      *IT_0343 + conj(IT_0314)*IT_0345 + conj(IT_0297)*IT_0354 + conj(IT_0152)
      *IT_0361 + conj(IT_0091)*IT_0369 + conj(IT_0290)*IT_0374 + conj(IT_0128)
      *IT_0376 + conj(IT_0064)*IT_0378;
    const ccomplex_t IT_0380 = 12*IT_0211;
    const ccomplex_t IT_0381 = conj(IT_0107)*IT_0317 + IT_0129*IT_0335 + conj
      (IT_0262)*IT_0336 + conj(IT_0276)*IT_0337 + conj(IT_0210)*IT_0338 + conj
      (IT_0244)*IT_0339 + IT_0133*IT_0341 + IT_0032*IT_0342 + conj(IT_0297)
      *IT_0343 + conj(IT_0314)*IT_0344 + conj(IT_0310)*IT_0354 + conj(IT_0152)
      *IT_0362 + conj(IT_0091)*IT_0368 + conj(IT_0064)*IT_0376 + conj(IT_0128)
      *IT_0378 + conj(IT_0290)*IT_0380;
    const ccomplex_t IT_0382 = (-12)*IT_0129;
    const ccomplex_t IT_0383 = (-2)*IT_0158;
    const ccomplex_t IT_0384 = (-3)*IT_0383;
    const ccomplex_t IT_0385 = 2*IT_0162;
    const ccomplex_t IT_0386 = (-3)*IT_0385;
    const ccomplex_t IT_0387 = conj(IT_0107)*IT_0322 + conj(IT_0310)*IT_0336 +
       conj(IT_0314)*IT_0338 + conj(IT_0297)*IT_0339 + conj(IT_0244)*IT_0343 +
       conj(IT_0210)*IT_0344 + conj(IT_0276)*IT_0345 + IT_0092*IT_0351 + IT_0315
      *IT_0352 + IT_0211*IT_0353 + conj(IT_0262)*IT_0354 + conj(IT_0091)*IT_0362
       + conj(IT_0152)*IT_0368 + conj(IT_0290)*IT_0382 + conj(IT_0128)*IT_0384 +
       conj(IT_0064)*IT_0386;
    const ccomplex_t IT_0388 = (-12)*IT_0211;
    const ccomplex_t IT_0389 = conj(IT_0107)*IT_0321 + conj(IT_0297)*IT_0336 +
       conj(IT_0314)*IT_0337 + conj(IT_0310)*IT_0339 + conj(IT_0262)*IT_0343 +
       conj(IT_0276)*IT_0344 + conj(IT_0210)*IT_0345 + IT_0032*IT_0351 + IT_0133
      *IT_0352 + IT_0129*IT_0353 + conj(IT_0244)*IT_0354 + conj(IT_0091)*IT_0361
       + conj(IT_0152)*IT_0369 + conj(IT_0064)*IT_0384 + conj(IT_0128)*IT_0386 +
       conj(IT_0290)*IT_0388;
    const ccomplex_t IT_0390 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0391 = (-0.333333333333333)*IT_0390;
    const ccomplex_t IT_0392 = IT_0390*IT_0391;
    const ccomplex_t IT_0393 = cpow(s_12 + IT_0001 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0394 = IT_0392*IT_0393;
    const ccomplex_t IT_0395 = (0 + _Complex_I*1)*IT_0394;
    const ccomplex_t IT_0396 = 0.5*IT_0395;
    const ccomplex_t IT_0397 = IT_0317 + IT_0318 + IT_0374 + IT_0380;
    const ccomplex_t IT_0398 = 24*IT_0211;
    const ccomplex_t IT_0399 = 48*IT_0315;
    const ccomplex_t IT_0400 = 48*IT_0130;
    const ccomplex_t IT_0401 = IT_0319 + IT_0398 + IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = IT_0319 + IT_0320 + IT_0398 + IT_0399;
    const ccomplex_t IT_0403 = IT_0321 + IT_0322 + IT_0382 + IT_0388;
    const ccomplex_t IT_0404 = s_13 + s_24;
    const ccomplex_t IT_0405 = -s_14 + -s_23;
    const ccomplex_t IT_0406 = IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = IT_0155*IT_0406;
    const ccomplex_t IT_0408 = (-8)*IT_0407;
    const ccomplex_t IT_0409 = 3*IT_0408;
    const ccomplex_t IT_0410 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0411 = IT_0404 + IT_0410;
    const ccomplex_t IT_0412 = IT_0155*IT_0411;
    const ccomplex_t IT_0413 = 24*IT_0412;
    const ccomplex_t IT_0414 = (-24)*IT_0412;
    const ccomplex_t IT_0415 = (-3)*IT_0408;
    const ccomplex_t IT_0416 = -IT_0397;
    const ccomplex_t IT_0417 = -IT_0403;
    const ccomplex_t IT_0418 = (96*IT_0129 + 96*IT_0130 + 96*IT_0133 + 96
      *IT_0211 + 192*IT_0315)*IT_0396 + -(conj(IT_0244) + conj(IT_0262))*IT_0397
       + (IT_0279 + -conj(IT_0279) + IT_0290 + -conj(IT_0290))*IT_0401 + 
      (IT_0031 + -conj(IT_0031) + IT_0107 + -conj(IT_0107))*IT_0402 + -(conj
      (IT_0297) + conj(IT_0310))*IT_0403 + -(conj(IT_0064) + conj(IT_0128))
      *IT_0409 + -(24*IT_0014 + IT_0285)*IT_0412 + -(conj(IT_0210) + conj
      (IT_0276) + IT_0314)*IT_0413 + -(IT_0210 + IT_0276 + conj(IT_0314))
      *IT_0414 + -(IT_0064 + IT_0128)*IT_0415 + -(IT_0244 + IT_0262)*IT_0416 + -
      (IT_0297 + IT_0310)*IT_0417;
    const ccomplex_t IT_0419 = 48*conj(IT_0279);
    const ccomplex_t IT_0420 = conj(IT_0107)*IT_0319 + IT_0324 + IT_0326 +
       conj(IT_0276)*IT_0340 + IT_0211*IT_0341 + IT_0263*IT_0342 + conj(IT_0314)
      *IT_0355 + conj(IT_0210)*IT_0359 + conj(IT_0152)*IT_0366 + conj(IT_0091)
      *IT_0372 + conj(IT_0244)*IT_0374 + conj(IT_0262)*IT_0380 + conj(IT_0310)
      *IT_0382 + conj(IT_0297)*IT_0388 + conj(IT_0290)*IT_0400 + IT_0315*IT_0419;
    const ccomplex_t IT_0421 = 24*IT_0035;
    const ccomplex_t IT_0422 = 0.333333333333333*conj(IT_0064);
    const ccomplex_t IT_0423 = 0.333333333333333*conj(IT_0297);
    const ccomplex_t IT_0424 = s_12*IT_0034;
    const ccomplex_t IT_0425 = (-24)*IT_0424;
    const ccomplex_t IT_0426 = 0.333333333333333*conj(IT_0128);
    const ccomplex_t IT_0427 = 0.333333333333333*conj(IT_0210);
    const ccomplex_t IT_0428 = 0.333333333333333*conj(IT_0244);
    const ccomplex_t IT_0429 = 0.333333333333333*conj(IT_0310);
    const ccomplex_t IT_0430 = 0.333333333333333*conj(IT_0262);
    const ccomplex_t IT_0431 = 0.333333333333333*conj(IT_0276);
    const ccomplex_t IT_0432 = 0.333333333333333*conj(IT_0314);
    const ccomplex_t IT_0433 = conj(IT_0107) + conj(IT_0290);
    const ccomplex_t IT_0434 = 0.333333333333333*IT_0433;
    const ccomplex_t IT_0435 = 3*conj(IT_0031);
    const ccomplex_t IT_0436 = 0.333333333333333*IT_0435;
    const ccomplex_t IT_0437 = IT_0135*IT_0367 + 3*IT_0421*IT_0422 + 3*IT_0386
      *IT_0423 + 3*IT_0425*IT_0426 + 3*IT_0333*IT_0427 + 3*IT_0376*IT_0428 + 3
      *IT_0384*IT_0429 + 3*IT_0378*IT_0430 + 3*IT_0347*IT_0431 + 3*IT_0350
      *IT_0432 + 3*IT_0325*IT_0434 + 3*IT_0163*(conj(IT_0279) + IT_0436);
    const ccomplex_t IT_0438 = IT_0035*IT_0371 + 3*IT_0384*IT_0423 + 3*IT_0422
      *IT_0425 + 3*IT_0421*IT_0426 + 3*IT_0347*IT_0427 + 3*IT_0378*IT_0428 + 3
      *IT_0386*IT_0429 + 3*IT_0376*IT_0430 + 3*IT_0333*IT_0431 + 3*IT_0356
      *IT_0432 + 3*IT_0323*IT_0434 + 3*IT_0159*(conj(IT_0279) + IT_0436);
    const ccomplex_t IT_0439 = IT_0014*((-24)*conj(IT_0031)*IT_0032 + 6
      *IT_0035*conj(IT_0064) + 6*IT_0066*conj(IT_0091) + (-24)*IT_0092*conj
      (IT_0107) + (-6)*conj(IT_0128)*IT_0135 + (-6)*conj(IT_0152)*IT_0154) + 
      (conj(IT_0064)*IT_0159 + conj(IT_0128)*IT_0163)*IT_0164 + (IT_0165*conj
      (IT_0210) + IT_0211*conj(IT_0244) + conj(IT_0091)*IT_0246 + IT_0129*conj
      (IT_0262) + IT_0263*conj(IT_0276))*IT_0280 + (conj(IT_0031)*IT_0129 + conj
      (IT_0107)*IT_0211)*IT_0281 + IT_0031*((-12)*conj(IT_0091)*IT_0246 + 12
      *conj(IT_0152)*IT_0283 + IT_0129*IT_0284 + IT_0032*IT_0285) + IT_0165*conj
      (IT_0279)*IT_0286 + conj(IT_0290)*(IT_0263*IT_0286 + IT_0211*IT_0291) +
       IT_0032*(conj(IT_0262)*IT_0286 + conj(IT_0210)*IT_0291 + IT_0292*conj
      (IT_0297)) + IT_0133*(conj(IT_0210)*IT_0286 + conj(IT_0262)*IT_0291 + conj
      (IT_0014)*IT_0292 + conj(IT_0031)*IT_0293 + conj(IT_0297)*IT_0298) +
       IT_0092*(conj(IT_0244)*IT_0286 + conj(IT_0276)*IT_0291 + IT_0292*conj
      (IT_0310) + IT_0298*conj(IT_0314)) + (conj(IT_0276)*IT_0286 + 48*IT_0279
      *conj(IT_0290) + conj(IT_0244)*IT_0291 + conj(IT_0107)*IT_0293 + IT_0298
      *conj(IT_0310) + IT_0292*conj(IT_0314))*IT_0315 + (conj(IT_0152)*IT_0283 +
       IT_0129*conj(IT_0297) + IT_0211*conj(IT_0310) + IT_0263*conj(IT_0314))
      *IT_0316 + IT_0107*IT_0332 + IT_0276*IT_0348 + IT_0314*IT_0358 + IT_0210
      *IT_0360 + IT_0152*IT_0370 + IT_0091*IT_0373 + IT_0244*IT_0379 + IT_0262
      *IT_0381 + IT_0310*IT_0387 + IT_0297*IT_0389 + -IT_0396*IT_0418 + 3
      *IT_0279*(conj(IT_0064)*IT_0159 + conj(IT_0128)*IT_0163 +
       0.333333333333333*IT_0165*IT_0342 + 0.333333333333333*IT_0130*IT_0419) +
       IT_0290*IT_0420 + IT_0128*IT_0437 + IT_0064*IT_0438;
    return create_ccomplex_return(IT_0439);
}

