#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_anti_s_s(
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
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = m_s*e_em*IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = V_u1*U_Wm1;
    const ccomplex_t IT_0009 = IT_0001*IT_0008;
    const ccomplex_t IT_0010 = U_d1*V_Wp1;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*(IT_0009 + IT_0011);
    const ccomplex_t IT_0013 = 1.4142135623731*e_em*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = pow(m_C_1, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + (-2)*IT_0016 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_Z, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + (-2)*IT_0016 + IT_0020 + 
      -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0001*IT_0010;
    const ccomplex_t IT_0023 = IT_0003*IT_0008;
    const ccomplex_t IT_0024 = IT_0005*(IT_0022 + -IT_0023);
    const ccomplex_t IT_0025 = 1.4142135623731*e_em*IT_0024;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = m_s*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0019 + -IT_0031;
    const ccomplex_t IT_0033 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0036 = IT_0003*IT_0035;
    const ccomplex_t IT_0037 = IT_0005*(IT_0034 + IT_0036);
    const ccomplex_t IT_0038 = 1.4142135623731*e_em*IT_0037;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0007*IT_0039;
    const ccomplex_t IT_0041 = IT_0017*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0001*IT_0035;
    const ccomplex_t IT_0044 = IT_0003*IT_0033;
    const ccomplex_t IT_0045 = IT_0005*(IT_0043 + -IT_0044);
    const ccomplex_t IT_0046 = 1.4142135623731*e_em*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0028*IT_0047;
    const ccomplex_t IT_0049 = IT_0021*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = -IT_0042 + -IT_0050;
    const ccomplex_t IT_0052 = pow(m_s, 2);
    const ccomplex_t IT_0053 = s_34 + IT_0052;
    const ccomplex_t IT_0054 = s_12*IT_0053;
    const ccomplex_t IT_0055 = cos(alpha);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0005*IT_0055;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0035*IT_0055;
    const ccomplex_t IT_0059 = sin(alpha);
    const ccomplex_t IT_0060 = IT_0033*IT_0059;
    const ccomplex_t IT_0061 = IT_0005*(IT_0058 + IT_0060);
    const ccomplex_t IT_0062 = 1.4142135623731*e_em*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0057*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0005*IT_0059;
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0033*IT_0055;
    const ccomplex_t IT_0072 = IT_0035*IT_0059;
    const ccomplex_t IT_0073 = IT_0005*(IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = 1.4142135623731*e_em*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0070*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = -IT_0068 + -IT_0080;
    const ccomplex_t IT_0082 = IT_0008*IT_0055;
    const ccomplex_t IT_0083 = IT_0010*IT_0059;
    const ccomplex_t IT_0084 = IT_0005*(IT_0082 + -IT_0083);
    const ccomplex_t IT_0085 = 1.4142135623731*e_em*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0070*IT_0087;
    const ccomplex_t IT_0089 = IT_0078*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0010*IT_0055;
    const ccomplex_t IT_0092 = IT_0008*IT_0059;
    const ccomplex_t IT_0093 = IT_0005*(IT_0091 + IT_0092);
    const ccomplex_t IT_0094 = 1.4142135623731*e_em*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0057*IT_0096;
    const ccomplex_t IT_0098 = IT_0066*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0090 + -IT_0099;
    const ccomplex_t IT_0101 = -s_34;
    const ccomplex_t IT_0102 = IT_0052 + IT_0101;
    const ccomplex_t IT_0103 = s_12*IT_0102;
    const ccomplex_t IT_0104 = cpow(IT_0003, -1);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_us
      *e_em*IT_0000*IT_0005*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *V_us*e_em*IT_0000*IT_0005*IT_0104;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0106*IT_0108;
    const ccomplex_t IT_0110 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*V_us*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*V_us*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0121 = conj(V_ts)*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0122 = IT_0005*IT_0121;
    const ccomplex_t IT_0123 = IT_0005*IT_0104;
    const ccomplex_t IT_0124 = m_t*conj(V_ts)*conj(V_u1)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0122 + (-0.5)*IT_0126);
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = m_t*V_ts*V_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0130 = IT_0123*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = V_ts*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0135 = IT_0005*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0133 + IT_0138;
    const ccomplex_t IT_0140 = IT_0128*IT_0139;
    const ccomplex_t IT_0141 = IT_0120*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0144 = conj(V_ts)*e_em*conj(V_Wp1)*U_st_01;
    const ccomplex_t IT_0145 = IT_0005*IT_0144;
    const ccomplex_t IT_0146 = m_t*conj(V_ts)*conj(V_u1)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0147 = IT_0123*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + (-0.5)*IT_0148);
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = m_t*V_ts*V_u1*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0152 = IT_0123*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = V_ts*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0157 = IT_0005*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0155 + IT_0160;
    const ccomplex_t IT_0162 = IT_0150*IT_0161;
    const ccomplex_t IT_0163 = IT_0143*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0112 + IT_0119 + IT_0142 + IT_0164;
    const ccomplex_t IT_0166 = m_s*m_C_1;
    const ccomplex_t IT_0167 = s_23 + s_24;
    const ccomplex_t IT_0168 = IT_0166*IT_0167;
    const ccomplex_t IT_0169 = -s_14;
    const ccomplex_t IT_0170 = s_13 + IT_0169;
    const ccomplex_t IT_0171 = IT_0166*IT_0170;
    const ccomplex_t IT_0172 = (-4)*IT_0171;
    const ccomplex_t IT_0173 = cos(theta_W);
    const ccomplex_t IT_0174 = cpow(IT_0173, -1);
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0174;
    const ccomplex_t IT_0176 = 0.333333333333333*IT_0175;
    const ccomplex_t IT_0177 = IT_0004*IT_0174;
    const ccomplex_t IT_0178 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0179 = IT_0177*IT_0178;
    const ccomplex_t IT_0180 = IT_0005*IT_0173;
    const ccomplex_t IT_0181 = IT_0178*IT_0180;
    const ccomplex_t IT_0182 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0183 = IT_0180*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*(IT_0179 + -IT_0181 + (-2)
      *IT_0183);
    const ccomplex_t IT_0185 = 0.5*IT_0184;
    const ccomplex_t IT_0186 = IT_0176*IT_0185;
    const ccomplex_t IT_0187 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0188 = IT_0186*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = 3*IT_0081;
    const ccomplex_t IT_0191 = 3*conj(IT_0081);
    const ccomplex_t IT_0192 = e_em*IT_0177;
    const ccomplex_t IT_0193 = e_em*IT_0180;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*(IT_0192 + 3*IT_0193);
    const ccomplex_t IT_0195 = (-0.166666666666667)*IT_0194;
    const ccomplex_t IT_0196 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0197 = IT_0177*IT_0196;
    const ccomplex_t IT_0198 = IT_0180*IT_0196;
    const ccomplex_t IT_0199 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0200 = IT_0180*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0197 + -IT_0198 + (-2)
      *IT_0200);
    const ccomplex_t IT_0202 = (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = IT_0195*IT_0202;
    const ccomplex_t IT_0204 = IT_0187*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = -IT_0205;
    const ccomplex_t IT_0207 = 3*IT_0100;
    const ccomplex_t IT_0208 = 3*conj(IT_0100);
    const ccomplex_t IT_0209 = -s_24;
    const ccomplex_t IT_0210 = s_23 + IT_0209;
    const ccomplex_t IT_0211 = IT_0166*IT_0210;
    const ccomplex_t IT_0212 = 4*IT_0211;
    const ccomplex_t IT_0213 = s_13*s_24;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp1)
      *IT_0005;
    const ccomplex_t IT_0217 = IT_0215*IT_0216;
    const ccomplex_t IT_0218 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0219 = IT_0217*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0000*IT_0005*IT_0104;
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u1)*e_em*IT_0000*IT_0005*IT_0104;
    const ccomplex_t IT_0224 = (-0.5)*IT_0223;
    const ccomplex_t IT_0225 = IT_0222*IT_0224;
    const ccomplex_t IT_0226 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0227 = IT_0225*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = -IT_0220 + -IT_0228;
    const ccomplex_t IT_0230 = (-24)*conj(IT_0229);
    const ccomplex_t IT_0231 = 48*conj(IT_0206);
    const ccomplex_t IT_0232 = 24*conj(IT_0165) + IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_00);
    const ccomplex_t IT_0234 = (-0.5)*IT_0233;
    const ccomplex_t IT_0235 = IT_0128*IT_0234;
    const ccomplex_t IT_0236 = IT_0120*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = -IT_0237;
    const ccomplex_t IT_0239 = 12*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0241 = (-0.5)*IT_0240;
    const ccomplex_t IT_0242 = IT_0216*IT_0241;
    const ccomplex_t IT_0243 = IT_0218*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_us*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0246 = (-0.5)*IT_0245;
    const ccomplex_t IT_0247 = IT_0115*IT_0246;
    const ccomplex_t IT_0248 = IT_0117*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_01);
    const ccomplex_t IT_0251 = (-0.5)*IT_0250;
    const ccomplex_t IT_0252 = IT_0150*IT_0251;
    const ccomplex_t IT_0253 = IT_0143*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = -IT_0244 + -IT_0249 + IT_0254;
    const ccomplex_t IT_0256 = (-12)*IT_0255;
    const ccomplex_t IT_0257 = IT_0239 + IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*V_us
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0259 = 0.5*IT_0258;
    const ccomplex_t IT_0260 = IT_0246*IT_0259;
    const ccomplex_t IT_0261 = IT_0117*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = 12*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_cs)*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0266 = 0.5*IT_0265;
    const ccomplex_t IT_0267 = IT_0241*IT_0266;
    const ccomplex_t IT_0268 = IT_0218*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0000*IT_0002*IT_0005*U_st_00;
    const ccomplex_t IT_0271 = 0.5*IT_0270;
    const ccomplex_t IT_0272 = IT_0234*IT_0271;
    const ccomplex_t IT_0273 = IT_0120*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0000*IT_0002*IT_0005*U_st_01;
    const ccomplex_t IT_0276 = 0.5*IT_0275;
    const ccomplex_t IT_0277 = IT_0251*IT_0276;
    const ccomplex_t IT_0278 = IT_0143*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = IT_0269 + IT_0274 + IT_0279;
    const ccomplex_t IT_0281 = (-12)*IT_0280;
    const ccomplex_t IT_0282 = IT_0264 + IT_0281;
    const ccomplex_t IT_0283 = 12*IT_0280;
    const ccomplex_t IT_0284 = IT_0114*IT_0259;
    const ccomplex_t IT_0285 = IT_0117*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = IT_0139*IT_0271;
    const ccomplex_t IT_0288 = IT_0120*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0161*IT_0276;
    const ccomplex_t IT_0291 = IT_0143*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = -IT_0286 + IT_0289 + IT_0292;
    const ccomplex_t IT_0294 = 12*IT_0293;
    const ccomplex_t IT_0295 = 12*IT_0255;
    const ccomplex_t IT_0296 = 12*IT_0229;
    const ccomplex_t IT_0297 = (-24)*IT_0229;
    const ccomplex_t IT_0298 = (-24)*IT_0263;
    const ccomplex_t IT_0299 = (-24)*conj(IT_0263);
    const ccomplex_t IT_0300 = (-12)*conj(IT_0280);
    const ccomplex_t IT_0301 = (-12)*conj(IT_0255);
    const ccomplex_t IT_0302 = IT_0206*IT_0232 + conj(IT_0238)*IT_0257 + conj
      (IT_0263)*IT_0282 + conj(IT_0280)*IT_0283 + conj(IT_0293)*IT_0294 + conj
      (IT_0255)*IT_0295 + conj(IT_0229)*IT_0296 + conj(IT_0206)*(24*IT_0165 +
       IT_0297) + conj(IT_0189)*(24*IT_0280 + IT_0298) + IT_0189*(48*conj
      (IT_0189) + 24*conj(IT_0280) + IT_0299) + IT_0263*IT_0300 + IT_0238*IT_0301;
    const ccomplex_t IT_0303 = 12*IT_0165;
    const ccomplex_t IT_0304 = IT_0215*IT_0266;
    const ccomplex_t IT_0305 = IT_0218*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = (-12)*IT_0293;
    const ccomplex_t IT_0308 = (-12)*conj(IT_0293);
    const ccomplex_t IT_0309 = (-12)*IT_0165;
    const ccomplex_t IT_0310 = (-12)*conj(IT_0165);
    const ccomplex_t IT_0311 = IT_0176*IT_0202;
    const ccomplex_t IT_0312 = IT_0187*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*IT_0312;
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = s_12*IT_0052;
    const ccomplex_t IT_0316 = (-12)*IT_0229;
    const ccomplex_t IT_0317 = 24*IT_0206;
    const ccomplex_t IT_0318 = IT_0303 + IT_0316 + IT_0317;
    const ccomplex_t IT_0319 = m_s*s_14*m_C_1;
    const ccomplex_t IT_0320 = (-12)*IT_0306;
    const ccomplex_t IT_0321 = IT_0294 + IT_0320;
    const ccomplex_t IT_0322 = s_34*IT_0016;
    const ccomplex_t IT_0323 = (-12)*IT_0263;
    const ccomplex_t IT_0324 = 24*IT_0189;
    const ccomplex_t IT_0325 = IT_0283 + IT_0323 + IT_0324;
    const ccomplex_t IT_0326 = m_s*s_23*m_C_1;
    const ccomplex_t IT_0327 = (-12)*IT_0238;
    const ccomplex_t IT_0328 = IT_0295 + IT_0327;
    const ccomplex_t IT_0329 = 12*IT_0051;
    const ccomplex_t IT_0330 = s_13 + s_14;
    const ccomplex_t IT_0331 = IT_0166*IT_0330;
    const ccomplex_t IT_0332 = (-12)*IT_0032;
    const ccomplex_t IT_0333 = IT_0185*IT_0195;
    const ccomplex_t IT_0334 = IT_0187*IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = IT_0016*IT_0102;
    const ccomplex_t IT_0337 = 6*conj(IT_0081);
    const ccomplex_t IT_0338 = 6*IT_0100;
    const ccomplex_t IT_0339 = (-6)*conj(IT_0081);
    const ccomplex_t IT_0340 = (-6)*IT_0100;
    const ccomplex_t IT_0341 = (-6)*conj(IT_0100);
    const ccomplex_t IT_0342 = 24*conj(IT_0081)*IT_0100 + 6*conj(IT_0100)
      *IT_0238 + IT_0081*(24*conj(IT_0100) + (-6)*conj(IT_0293) + 6*conj(IT_0306
      )) + IT_0306*IT_0337 + conj(IT_0238)*IT_0338 + IT_0293*IT_0339 + conj
      (IT_0255)*IT_0340 + IT_0255*IT_0341;
    const ccomplex_t IT_0343 = s_14*s_23;
    const ccomplex_t IT_0344 = s_12*s_34;
    const ccomplex_t IT_0345 = -IT_0344;
    const ccomplex_t IT_0346 = -IT_0213;
    const ccomplex_t IT_0347 = IT_0315 + IT_0343 + IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = -IT_0081*(conj(IT_0238) + -conj(IT_0255)) +
       conj(IT_0100)*IT_0293 + 0.166666666666667*IT_0255*IT_0337 +
       0.166666666666667*conj(IT_0293)*IT_0338 + 0.166666666666667*IT_0238
      *IT_0339 + 0.166666666666667*conj(IT_0306)*IT_0340 + 0.166666666666667
      *IT_0306*IT_0341;
    const ccomplex_t IT_0349 = IT_0016*IT_0053;
    const ccomplex_t IT_0350 = 6*IT_0032;
    const ccomplex_t IT_0351 = 6*conj(IT_0032);
    const ccomplex_t IT_0352 = 6*IT_0051;
    const ccomplex_t IT_0353 = 6*conj(IT_0051);
    const ccomplex_t IT_0354 = (-6)*IT_0032;
    const ccomplex_t IT_0355 = (-6)*conj(IT_0032);
    const ccomplex_t IT_0356 = 4*conj(IT_0051)*(IT_0032 + 1./4*IT_0293) + 4
      *IT_0051*(conj(IT_0032) + 1./4*conj(IT_0293)) + (-0.166666666666667)*conj
      (IT_0255)*IT_0350 + (-0.166666666666667)*IT_0255*IT_0351 + (
      -0.166666666666667)*conj(IT_0306)*IT_0352 + (-0.166666666666667)*IT_0306
      *IT_0353 + (-0.166666666666667)*conj(IT_0238)*IT_0354 + (
      -0.166666666666667)*IT_0238*IT_0355;
    const ccomplex_t IT_0357 = 12*IT_0032;
    const ccomplex_t IT_0358 = 12*conj(IT_0032);
    const ccomplex_t IT_0359 = (-12)*IT_0051;
    const ccomplex_t IT_0360 = (-12)*conj(IT_0051);
    const ccomplex_t IT_0361 = conj(IT_0051)*IT_0229 + IT_0051*conj(IT_0229) +
       0.166666666666667*conj(IT_0280)*IT_0350 + 0.166666666666667*IT_0280
      *IT_0351 + 0.166666666666667*conj(IT_0263)*IT_0354 + 0.166666666666667
      *IT_0263*IT_0355 + 0.166666666666667*conj(IT_0189)*IT_0357 +
       0.166666666666667*IT_0189*IT_0358 + 0.166666666666667*conj(IT_0206)
      *IT_0359 + 0.166666666666667*IT_0206*IT_0360;
    const ccomplex_t IT_0362 = -IT_0343;
    const ccomplex_t IT_0363 = IT_0315 + IT_0362;
    const ccomplex_t IT_0364 = IT_0213 + IT_0344 + IT_0363;
    const ccomplex_t IT_0365 = (-6)*IT_0051;
    const ccomplex_t IT_0366 = (-6)*conj(IT_0051);
    const ccomplex_t IT_0367 = conj(IT_0165)*IT_0350 + IT_0165*IT_0351 + conj
      (IT_0263)*IT_0352 + IT_0263*IT_0353 + conj(IT_0229)*IT_0354 + IT_0229
      *IT_0355 + conj(IT_0206)*IT_0357 + IT_0206*IT_0358 + conj(IT_0189)*IT_0359
       + IT_0189*IT_0360 + conj(IT_0280)*IT_0365 + IT_0280*IT_0366;
    const ccomplex_t IT_0368 = (-2)*IT_0211;
    const ccomplex_t IT_0369 = (-3)*IT_0081;
    const ccomplex_t IT_0370 = (-3)*conj(IT_0081);
    const ccomplex_t IT_0371 = (-3)*IT_0100;
    const ccomplex_t IT_0372 = (-3)*conj(IT_0100);
    const ccomplex_t IT_0373 = 2*IT_0171;
    const ccomplex_t IT_0374 = (-2)*IT_0171;
    const ccomplex_t IT_0375 = 2*IT_0211;
    const ccomplex_t IT_0376 = 24*IT_0322;
    const ccomplex_t IT_0377 = IT_0016*IT_0052;
    const ccomplex_t IT_0378 = 48*IT_0377;
    const ccomplex_t IT_0379 = 24*IT_0315;
    const ccomplex_t IT_0380 = 48*IT_0213;
    const ccomplex_t IT_0381 = IT_0376 + IT_0378 + IT_0379 + IT_0380;
    const ccomplex_t IT_0382 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0383 = (-0.333333333333333)*IT_0382;
    const ccomplex_t IT_0384 = IT_0382*IT_0383;
    const ccomplex_t IT_0385 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0386 = IT_0384*IT_0385;
    const ccomplex_t IT_0387 = (0 + _Complex_I*1)*IT_0386;
    const ccomplex_t IT_0388 = 0.5*IT_0387;
    const ccomplex_t IT_0389 = -IT_0388;
    const ccomplex_t IT_0390 = m_s*s_13*m_C_1;
    const ccomplex_t IT_0391 = 12*conj(IT_0280);
    const ccomplex_t IT_0392 = (-12)*conj(IT_0263);
    const ccomplex_t IT_0393 = 24*conj(IT_0189);
    const ccomplex_t IT_0394 = IT_0391 + IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = 12*conj(IT_0293);
    const ccomplex_t IT_0396 = (-24)*IT_0238;
    const ccomplex_t IT_0397 = (-24)*conj(IT_0238);
    const ccomplex_t IT_0398 = 12*conj(IT_0263);
    const ccomplex_t IT_0399 = 12*conj(IT_0229);
    const ccomplex_t IT_0400 = 24*conj(IT_0206);
    const ccomplex_t IT_0401 = conj(IT_0238)*IT_0282 + conj(IT_0165)*IT_0294 +
       conj(IT_0229)*IT_0307 + IT_0229*IT_0308 + conj(IT_0306)*((-24)*IT_0206 +
       IT_0296 + IT_0309) + conj(IT_0293)*IT_0317 + conj(IT_0255)*IT_0325 +
       IT_0255*IT_0394 + IT_0165*IT_0395 + conj(IT_0189)*IT_0396 + IT_0189
      *IT_0397 + IT_0238*(IT_0300 + IT_0398) + IT_0306*((-24)*conj(IT_0206) +
       IT_0310 + IT_0399) + IT_0293*IT_0400;
    const ccomplex_t IT_0402 = 12*IT_0322;
    const ccomplex_t IT_0403 = 24*IT_0213;
    const ccomplex_t IT_0404 = 24*IT_0377;
    const ccomplex_t IT_0405 = 12*IT_0315;
    const ccomplex_t IT_0406 = IT_0402 + IT_0403 + IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = 48*IT_0343;
    const ccomplex_t IT_0408 = IT_0376 + IT_0378 + IT_0379 + IT_0407;
    const ccomplex_t IT_0409 = (-12)*IT_0322;
    const ccomplex_t IT_0410 = (-24)*IT_0213;
    const ccomplex_t IT_0411 = (-24)*IT_0377;
    const ccomplex_t IT_0412 = (-12)*IT_0315;
    const ccomplex_t IT_0413 = IT_0409 + IT_0410 + IT_0411 + IT_0412;
    const ccomplex_t IT_0414 = -s_14 + -s_23;
    const ccomplex_t IT_0415 = s_13 + s_24;
    const ccomplex_t IT_0416 = IT_0414 + IT_0415;
    const ccomplex_t IT_0417 = IT_0166*IT_0416;
    const ccomplex_t IT_0418 = (-8)*IT_0417;
    const ccomplex_t IT_0419 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0420 = IT_0415 + IT_0419;
    const ccomplex_t IT_0421 = IT_0166*IT_0420;
    const ccomplex_t IT_0422 = 24*conj(IT_0293);
    const ccomplex_t IT_0423 = -(conj(IT_0189) + conj(IT_0206))*IT_0381 + (96
      *IT_0213 + 96*IT_0315 + 96*IT_0322 + 96*IT_0343 + 192*IT_0377)*IT_0388 + 
      (IT_0165 + -conj(IT_0165) + IT_0280 + -conj(IT_0280))*IT_0406 + -(conj
      (IT_0314) + conj(IT_0335))*IT_0408 + (IT_0229 + -conj(IT_0229) + IT_0263 +
       -conj(IT_0263))*IT_0413 + -(IT_0191 + IT_0208 + IT_0369 + IT_0371)
      *IT_0418 + -IT_0421*(24*IT_0238 + (-24)*IT_0255 + 24*conj(IT_0255) + (-24)
      *IT_0293 + 24*IT_0306 + (-24)*conj(IT_0306) + IT_0397 + IT_0422);
    const ccomplex_t IT_0424 = m_s*s_24*m_C_1;
    const ccomplex_t IT_0425 = 12*conj(IT_0165);
    const ccomplex_t IT_0426 = (-12)*conj(IT_0229);
    const ccomplex_t IT_0427 = IT_0400 + IT_0425 + IT_0426;
    const ccomplex_t IT_0428 = IT_0296 + IT_0309;
    const ccomplex_t IT_0429 = IT_0310 + IT_0399;
    const ccomplex_t IT_0430 = 24*conj(IT_0189)*IT_0293 + conj(IT_0280)
      *IT_0294 + ((-24)*IT_0189 + IT_0264 + IT_0281)*conj(IT_0306) + conj
      (IT_0263)*IT_0307 + IT_0263*IT_0308 + conj(IT_0255)*IT_0318 + IT_0280
      *IT_0395 + conj(IT_0206)*IT_0396 + IT_0206*IT_0397 + IT_0306*((-24)*conj
      (IT_0189) + IT_0300 + IT_0398) + IT_0189*IT_0422 + IT_0255*IT_0427 + conj
      (IT_0238)*IT_0428 + IT_0238*IT_0429;
    const ccomplex_t IT_0431 = 48*conj(IT_0314);
    const ccomplex_t IT_0432 = IT_0189*IT_0232 + conj(IT_0255)*IT_0294 + conj
      (IT_0189)*(24*IT_0165 + 48*IT_0206 + IT_0297) + conj(IT_0206)*IT_0298 +
       IT_0206*IT_0299 + IT_0301*IT_0306 + IT_0257*conj(IT_0306) + conj(IT_0238)
      *(12*IT_0306 + IT_0307) + IT_0238*IT_0308 + conj(IT_0280)*IT_0318 + 48
      *IT_0314*conj(IT_0335) + IT_0255*IT_0395 + IT_0280*IT_0427 + conj(IT_0263)
      *IT_0428 + IT_0263*IT_0429 + IT_0335*IT_0431;
    const ccomplex_t IT_0433 = (-12)*conj(IT_0306);
    const ccomplex_t IT_0434 = IT_0395 + IT_0433;
    const ccomplex_t IT_0435 = 12*conj(IT_0255);
    const ccomplex_t IT_0436 = (-12)*conj(IT_0238);
    const ccomplex_t IT_0437 = IT_0435 + IT_0436;
    const ccomplex_t IT_0438 = 12*conj(IT_0051);
    const ccomplex_t IT_0439 = (-12)*conj(IT_0032);
    const ccomplex_t IT_0440 = IT_0389*IT_0408;
    const ccomplex_t IT_0441 = IT_0335*(IT_0172*IT_0191 + IT_0208*IT_0212 + 48
      *conj(IT_0335)*IT_0343 + IT_0315*IT_0394 + IT_0322*IT_0427 + IT_0326
      *IT_0434 + IT_0319*IT_0437 + IT_0331*IT_0438 + IT_0168*IT_0439 + IT_0440);
    const ccomplex_t IT_0442 = (24*IT_0032*conj(IT_0032) + 24*IT_0051*conj
      (IT_0051))*IT_0054 + ((-24)*IT_0081*conj(IT_0081) + (-24)*IT_0100*conj
      (IT_0100))*IT_0103 + (-6)*(conj(IT_0051)*IT_0165 + IT_0051*conj(IT_0165))
      *IT_0168 + IT_0172*(conj(IT_0189)*IT_0190 + IT_0189*IT_0191 + conj(IT_0206
      )*IT_0207 + IT_0206*IT_0208) + (IT_0191*IT_0206 + IT_0190*conj(IT_0206) +
       conj(IT_0189)*IT_0207 + IT_0189*IT_0208)*IT_0212 + IT_0213*IT_0302 +
       IT_0213*(conj(IT_0165)*IT_0303 + conj(IT_0306)*IT_0307 + IT_0306*(12*conj
      (IT_0306) + IT_0308) + conj(IT_0229)*IT_0309 + IT_0229*IT_0310) + conj
      (IT_0314)*(IT_0172*IT_0207 + IT_0190*IT_0212 + IT_0315*IT_0318 + IT_0319
      *IT_0321 + IT_0322*IT_0325 + IT_0326*IT_0328 + IT_0168*IT_0329 + IT_0331
      *IT_0332) + (IT_0172*IT_0190 + IT_0207*IT_0212 + IT_0318*IT_0322 + IT_0315
      *IT_0325 + IT_0321*IT_0326 + IT_0319*IT_0328 + IT_0329*IT_0331 + IT_0168
      *IT_0332)*conj(IT_0335) + IT_0336*IT_0342 + 6*IT_0347*IT_0348 + (-6)
      *IT_0349*IT_0356 + 6*IT_0168*IT_0361 + (conj(IT_0293)*IT_0350 + IT_0293
      *IT_0351 + conj(IT_0306)*IT_0354 + IT_0306*IT_0355)*IT_0364 + IT_0364*
      (conj(IT_0238)*IT_0352 + IT_0238*IT_0353 + conj(IT_0255)*IT_0365 + IT_0255
      *IT_0366) + IT_0331*IT_0367 + IT_0368*(conj(IT_0165)*IT_0369 + IT_0165
      *IT_0370 + conj(IT_0280)*IT_0371 + IT_0280*IT_0372) + (conj(IT_0280)
      *IT_0369 + IT_0280*IT_0370 + conj(IT_0165)*IT_0371 + IT_0165*IT_0372)
      *IT_0373 + (conj(IT_0263)*IT_0369 + IT_0263*IT_0370 + conj(IT_0229)
      *IT_0371 + IT_0229*IT_0372)*IT_0374 + (conj(IT_0229)*IT_0369 + IT_0229
      *IT_0370 + conj(IT_0263)*IT_0371 + IT_0263*IT_0372)*IT_0375 + (IT_0189 +
       IT_0206)*IT_0381*IT_0389 + IT_0390*IT_0401 + -IT_0388*IT_0423 + IT_0424
      *IT_0430 + IT_0377*IT_0432 + IT_0314*(IT_0172*IT_0208 + IT_0191*IT_0212 +
       IT_0322*IT_0394 + IT_0315*IT_0427 + IT_0343*IT_0431 + IT_0319*IT_0434 +
       IT_0326*IT_0437 + IT_0168*IT_0438 + IT_0331*IT_0439 + IT_0440) + IT_0441;
    return create_ccomplex_return(IT_0442);
}

