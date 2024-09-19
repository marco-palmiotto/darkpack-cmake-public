#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_anti_d_d(
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
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_td = param->V_td;
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
    const ccomplex_t IT_0006 = m_d*e_em*IT_0000*IT_0002*IT_0003*IT_0005;
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
    const ccomplex_t IT_0027 = m_d*e_em*IT_0000*IT_0005;
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
    const ccomplex_t IT_0052 = pow(m_d, 2);
    const ccomplex_t IT_0053 = s_34 + IT_0052;
    const ccomplex_t IT_0054 = s_12*IT_0053;
    const ccomplex_t IT_0055 = sin(alpha);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0055;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = cos(alpha);
    const ccomplex_t IT_0059 = IT_0008*IT_0058;
    const ccomplex_t IT_0060 = IT_0010*IT_0055;
    const ccomplex_t IT_0061 = IT_0005*(IT_0059 + -IT_0060);
    const ccomplex_t IT_0062 = 1.4142135623731*e_em*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0057*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0058;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = IT_0010*IT_0058;
    const ccomplex_t IT_0072 = IT_0008*IT_0055;
    const ccomplex_t IT_0073 = IT_0005*(IT_0071 + IT_0072);
    const ccomplex_t IT_0074 = 1.4142135623731*e_em*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0070*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_12 + (-2)*IT_0016 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = -IT_0068 + -IT_0080;
    const ccomplex_t IT_0082 = IT_0033*IT_0058;
    const ccomplex_t IT_0083 = IT_0035*IT_0055;
    const ccomplex_t IT_0084 = IT_0005*(IT_0082 + -IT_0083);
    const ccomplex_t IT_0085 = 1.4142135623731*e_em*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0057*IT_0087;
    const ccomplex_t IT_0089 = IT_0066*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0035*IT_0058;
    const ccomplex_t IT_0092 = IT_0033*IT_0055;
    const ccomplex_t IT_0093 = IT_0005*(IT_0091 + IT_0092);
    const ccomplex_t IT_0094 = 1.4142135623731*e_em*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0070*IT_0096;
    const ccomplex_t IT_0098 = IT_0078*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0090 + -IT_0099;
    const ccomplex_t IT_0101 = -s_34;
    const ccomplex_t IT_0102 = IT_0052 + IT_0101;
    const ccomplex_t IT_0103 = s_12*IT_0102;
    const ccomplex_t IT_0104 = m_d*m_C_1;
    const ccomplex_t IT_0105 = -s_24;
    const ccomplex_t IT_0106 = s_23 + IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = 4*IT_0107;
    const ccomplex_t IT_0109 = cos(theta_W);
    const ccomplex_t IT_0110 = cpow(IT_0109, -1);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0110;
    const ccomplex_t IT_0112 = 0.333333333333333*IT_0111;
    const ccomplex_t IT_0113 = IT_0004*IT_0110;
    const ccomplex_t IT_0114 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = IT_0005*IT_0109;
    const ccomplex_t IT_0117 = IT_0114*IT_0116;
    const ccomplex_t IT_0118 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0115 + -IT_0117 + (-2)
      *IT_0119);
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = IT_0112*IT_0121;
    const ccomplex_t IT_0123 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0016 + -reg_prop, -1);
    const ccomplex_t IT_0124 = IT_0122*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = 3*IT_0081;
    const ccomplex_t IT_0127 = 3*conj(IT_0081);
    const ccomplex_t IT_0128 = e_em*IT_0113;
    const ccomplex_t IT_0129 = e_em*IT_0116;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0128 + 3*IT_0129);
    const ccomplex_t IT_0131 = (-0.166666666666667)*IT_0130;
    const ccomplex_t IT_0132 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0133 = IT_0113*IT_0132;
    const ccomplex_t IT_0134 = IT_0116*IT_0132;
    const ccomplex_t IT_0135 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0136 = IT_0116*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + -IT_0134 + (-2)
      *IT_0136);
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = IT_0131*IT_0138;
    const ccomplex_t IT_0140 = IT_0123*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = 3*IT_0100;
    const ccomplex_t IT_0144 = 3*conj(IT_0100);
    const ccomplex_t IT_0145 = -s_14;
    const ccomplex_t IT_0146 = s_13 + IT_0145;
    const ccomplex_t IT_0147 = IT_0104*IT_0146;
    const ccomplex_t IT_0148 = (-4)*IT_0147;
    const ccomplex_t IT_0149 = IT_0121*IT_0131;
    const ccomplex_t IT_0150 = IT_0123*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = s_23 + s_24;
    const ccomplex_t IT_0153 = IT_0104*IT_0152;
    const ccomplex_t IT_0154 = s_34*IT_0016;
    const ccomplex_t IT_0155 = cpow(IT_0003, -1);
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u1
      *e_em*IT_0000*IT_0005*IT_0155;
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u1)*e_em*IT_0000*IT_0005*IT_0155;
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0157*IT_0159;
    const ccomplex_t IT_0161 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0162 = IT_0160*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*V_cd*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*conj(V_cd)*e_em*conj(V_Wp1)
      *IT_0005;
    const ccomplex_t IT_0167 = IT_0165*IT_0166;
    const ccomplex_t IT_0168 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0169 = IT_0167*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0174 = IT_0172*IT_0173;
    const ccomplex_t IT_0175 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0176 = IT_0174*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0179 = conj(V_td)*e_em*conj(V_Wp1)*U_st_01;
    const ccomplex_t IT_0180 = IT_0005*IT_0179;
    const ccomplex_t IT_0181 = IT_0005*IT_0155;
    const ccomplex_t IT_0182 = m_t*conj(V_td)*conj(V_u1)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0183 = IT_0181*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0180 + (-0.5)*IT_0184);
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = m_t*V_td*V_u1*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0188 = IT_0181*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = (-0.5)*IT_0190;
    const ccomplex_t IT_0192 = V_td*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0193 = IT_0005*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0191 + IT_0196;
    const ccomplex_t IT_0198 = IT_0186*IT_0197;
    const ccomplex_t IT_0199 = IT_0178*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = IT_0163 + IT_0170 + IT_0177 + IT_0200;
    const ccomplex_t IT_0202 = 12*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0000*IT_0005*IT_0155;
    const ccomplex_t IT_0204 = 0.5*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *V_ud*e_em*IT_0000*IT_0005*IT_0155;
    const ccomplex_t IT_0206 = (-0.5)*IT_0205;
    const ccomplex_t IT_0207 = IT_0204*IT_0206;
    const ccomplex_t IT_0208 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0209 = IT_0207*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = cpow((-2)*s_13 + IT_0016 + IT_0052 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0212 = conj(V_td)*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0213 = IT_0005*IT_0212;
    const ccomplex_t IT_0214 = m_t*conj(V_td)*conj(V_u1)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0215 = IT_0181*IT_0214;
    const ccomplex_t IT_0216 = 1.4142135623731*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*(IT_0213 + (-0.5)*IT_0216);
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = m_t*V_td*V_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0220 = IT_0181*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = (-0.5)*IT_0222;
    const ccomplex_t IT_0224 = V_td*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0225 = IT_0005*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = IT_0223 + IT_0228;
    const ccomplex_t IT_0230 = IT_0218*IT_0229;
    const ccomplex_t IT_0231 = IT_0211*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = -IT_0210 + -IT_0232;
    const ccomplex_t IT_0234 = (-12)*IT_0233;
    const ccomplex_t IT_0235 = 24*IT_0142;
    const ccomplex_t IT_0236 = IT_0202 + IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = m_d*s_23*m_C_1;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*conj
      (V_cd)*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0239 = 0.5*IT_0238;
    const ccomplex_t IT_0240 = IT_0165*IT_0239;
    const ccomplex_t IT_0241 = IT_0168*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*conj
      (V_td)*e_em*IT_0000*IT_0002*IT_0005*U_st_00;
    const ccomplex_t IT_0244 = 0.5*IT_0243;
    const ccomplex_t IT_0245 = IT_0229*IT_0244;
    const ccomplex_t IT_0246 = IT_0211*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*conj
      (V_td)*e_em*IT_0000*IT_0002*IT_0005*U_st_01;
    const ccomplex_t IT_0249 = 0.5*IT_0248;
    const ccomplex_t IT_0250 = IT_0197*IT_0249;
    const ccomplex_t IT_0251 = IT_0178*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = -IT_0242 + IT_0247 + IT_0252;
    const ccomplex_t IT_0254 = 12*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*V_ud
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0256 = 0.5*IT_0255;
    const ccomplex_t IT_0257 = IT_0172*IT_0256;
    const ccomplex_t IT_0258 = IT_0175*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = (-12)*IT_0259;
    const ccomplex_t IT_0261 = IT_0254 + IT_0260;
    const ccomplex_t IT_0262 = s_12*IT_0052;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_cd*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0264 = (-0.5)*IT_0263;
    const ccomplex_t IT_0265 = IT_0239*IT_0264;
    const ccomplex_t IT_0266 = IT_0168*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = IT_0256*IT_0269;
    const ccomplex_t IT_0271 = IT_0175*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_td*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_01);
    const ccomplex_t IT_0274 = (-0.5)*IT_0273;
    const ccomplex_t IT_0275 = IT_0249*IT_0274;
    const ccomplex_t IT_0276 = IT_0178*IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*IT_0276;
    const ccomplex_t IT_0278 = IT_0267 + IT_0272 + IT_0277;
    const ccomplex_t IT_0279 = 12*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_td*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_00);
    const ccomplex_t IT_0281 = (-0.5)*IT_0280;
    const ccomplex_t IT_0282 = IT_0244*IT_0281;
    const ccomplex_t IT_0283 = IT_0211*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = -IT_0284;
    const ccomplex_t IT_0286 = (-12)*IT_0285;
    const ccomplex_t IT_0287 = 24*IT_0125;
    const ccomplex_t IT_0288 = IT_0279 + IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = m_d*s_14*m_C_1;
    const ccomplex_t IT_0290 = IT_0218*IT_0281;
    const ccomplex_t IT_0291 = IT_0211*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = IT_0166*IT_0264;
    const ccomplex_t IT_0294 = IT_0168*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*IT_0294;
    const ccomplex_t IT_0296 = IT_0292 + -IT_0295;
    const ccomplex_t IT_0297 = 12*IT_0296;
    const ccomplex_t IT_0298 = IT_0173*IT_0269;
    const ccomplex_t IT_0299 = IT_0175*IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = IT_0186*IT_0274;
    const ccomplex_t IT_0302 = IT_0178*IT_0301;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*IT_0302;
    const ccomplex_t IT_0304 = IT_0300 + -IT_0303;
    const ccomplex_t IT_0305 = (-12)*IT_0304;
    const ccomplex_t IT_0306 = IT_0297 + IT_0305;
    const ccomplex_t IT_0307 = s_13 + s_14;
    const ccomplex_t IT_0308 = IT_0104*IT_0307;
    const ccomplex_t IT_0309 = 12*IT_0051;
    const ccomplex_t IT_0310 = IT_0112*IT_0138;
    const ccomplex_t IT_0311 = IT_0123*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*IT_0311;
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = s_14*s_23;
    const ccomplex_t IT_0315 = 48*IT_0313;
    const ccomplex_t IT_0316 = IT_0016*IT_0052;
    const ccomplex_t IT_0317 = (-24)*conj(IT_0285);
    const ccomplex_t IT_0318 = 48*conj(IT_0125);
    const ccomplex_t IT_0319 = 24*conj(IT_0278) + IT_0317 + IT_0318;
    const ccomplex_t IT_0320 = 12*conj(IT_0278);
    const ccomplex_t IT_0321 = (-12)*conj(IT_0285);
    const ccomplex_t IT_0322 = 24*conj(IT_0125);
    const ccomplex_t IT_0323 = IT_0320 + IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = 12*IT_0285;
    const ccomplex_t IT_0325 = (-12)*IT_0278;
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = 12*conj(IT_0285);
    const ccomplex_t IT_0328 = (-12)*conj(IT_0278);
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = 12*IT_0304;
    const ccomplex_t IT_0331 = (-12)*IT_0296;
    const ccomplex_t IT_0332 = IT_0330 + IT_0331;
    const ccomplex_t IT_0333 = 12*conj(IT_0253);
    const ccomplex_t IT_0334 = (-24)*IT_0285;
    const ccomplex_t IT_0335 = (-24)*IT_0233;
    const ccomplex_t IT_0336 = (-24)*conj(IT_0233);
    const ccomplex_t IT_0337 = (-12)*IT_0253;
    const ccomplex_t IT_0338 = (-12)*conj(IT_0253);
    const ccomplex_t IT_0339 = (-12)*conj(IT_0296);
    const ccomplex_t IT_0340 = conj(IT_0201)*IT_0288 + IT_0254*conj(IT_0296) +
       48*IT_0151*conj(IT_0313) + conj(IT_0151)*IT_0315 + IT_0142*IT_0319 +
       IT_0201*IT_0323 + conj(IT_0233)*IT_0326 + IT_0233*IT_0329 + conj(IT_0259)
      *IT_0332 + IT_0296*IT_0333 + conj(IT_0142)*(48*IT_0125 + 24*IT_0278 +
       IT_0334) + conj(IT_0125)*IT_0335 + IT_0125*IT_0336 + conj(IT_0304)*(12
      *IT_0259 + IT_0337) + IT_0304*IT_0338 + IT_0259*IT_0339;
    const ccomplex_t IT_0341 = s_13*s_24;
    const ccomplex_t IT_0342 = 12*IT_0233;
    const ccomplex_t IT_0343 = (-12)*IT_0201;
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = (-12)*conj(IT_0201);
    const ccomplex_t IT_0346 = conj(IT_0201)*IT_0202 + conj(IT_0253)*IT_0254 +
       conj(IT_0278)*IT_0279 + conj(IT_0296)*IT_0297 + IT_0125*IT_0319 + conj
      (IT_0285)*IT_0326 + IT_0285*IT_0328 + conj(IT_0304)*IT_0332 + conj(IT_0125
      )*(24*IT_0278 + IT_0334) + conj(IT_0142)*(24*IT_0201 + IT_0335) + IT_0142*
      (48*conj(IT_0142) + 24*conj(IT_0201) + IT_0336) + conj(IT_0259)*IT_0337 +
       IT_0259*(12*conj(IT_0259) + IT_0338) + IT_0304*IT_0339 + conj(IT_0233)
      *IT_0344 + IT_0233*IT_0345;
    const ccomplex_t IT_0347 = 12*conj(IT_0032);
    const ccomplex_t IT_0348 = (-6)*IT_0032;
    const ccomplex_t IT_0349 = (-6)*conj(IT_0032);
    const ccomplex_t IT_0350 = IT_0016*IT_0053;
    const ccomplex_t IT_0351 = 0.166666666666667*conj(IT_0051)*((-24)*IT_0032 
      + (-6)*IT_0253 + 6*IT_0259) + 0.166666666666667*IT_0051*((-24)*conj
      (IT_0032) + (-6)*conj(IT_0253) + 6*conj(IT_0259)) + conj(IT_0032)*IT_0296 
      + IT_0032*conj(IT_0296) + 0.166666666666667*conj(IT_0304)*IT_0348 +
       0.166666666666667*IT_0304*IT_0349;
    const ccomplex_t IT_0352 = (-12)*IT_0051;
    const ccomplex_t IT_0353 = (-12)*conj(IT_0051);
    const ccomplex_t IT_0354 = (-6)*IT_0051;
    const ccomplex_t IT_0355 = (-6)*conj(IT_0051);
    const ccomplex_t IT_0356 = IT_0016*IT_0102;
    const ccomplex_t IT_0357 = (-6)*IT_0081;
    const ccomplex_t IT_0358 = (-6)*conj(IT_0081);
    const ccomplex_t IT_0359 = (-6)*IT_0100;
    const ccomplex_t IT_0360 = (-6)*conj(IT_0100);
    const ccomplex_t IT_0361 = 6*conj(IT_0100)*IT_0259 + 6*IT_0100*conj
      (IT_0259) + conj(IT_0081)*(24*IT_0100 + 6*IT_0304) + IT_0081*(24*conj
      (IT_0100) + 6*conj(IT_0304)) + conj(IT_0296)*IT_0357 + IT_0296*IT_0358 +
       conj(IT_0253)*IT_0359 + IT_0253*IT_0360;
    const ccomplex_t IT_0362 = 2*IT_0147;
    const ccomplex_t IT_0363 = (-3)*IT_0081;
    const ccomplex_t IT_0364 = (-3)*conj(IT_0081);
    const ccomplex_t IT_0365 = (-3)*IT_0100;
    const ccomplex_t IT_0366 = (-3)*conj(IT_0100);
    const ccomplex_t IT_0367 = (-2)*IT_0147;
    const ccomplex_t IT_0368 = (-2)*IT_0107;
    const ccomplex_t IT_0369 = 2*IT_0107;
    const ccomplex_t IT_0370 = 24*IT_0154;
    const ccomplex_t IT_0371 = 48*IT_0316;
    const ccomplex_t IT_0372 = 24*IT_0262;
    const ccomplex_t IT_0373 = 48*IT_0341;
    const ccomplex_t IT_0374 = IT_0370 + IT_0371 + IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0376 = (-0.333333333333333)*IT_0375;
    const ccomplex_t IT_0377 = IT_0375*IT_0376;
    const ccomplex_t IT_0378 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0379 = IT_0377*IT_0378;
    const ccomplex_t IT_0380 = (0 + _Complex_I*1)*IT_0379;
    const ccomplex_t IT_0381 = 0.5*IT_0380;
    const ccomplex_t IT_0382 = -IT_0381;
    const ccomplex_t IT_0383 = m_d*s_13*m_C_1;
    const ccomplex_t IT_0384 = (-24)*IT_0304;
    const ccomplex_t IT_0385 = (-24)*conj(IT_0304);
    const ccomplex_t IT_0386 = 12*conj(IT_0233);
    const ccomplex_t IT_0387 = 24*conj(IT_0142);
    const ccomplex_t IT_0388 = IT_0235*conj(IT_0253) + conj(IT_0201)*IT_0254 +
       IT_0288*conj(IT_0296) + IT_0296*IT_0323 + conj(IT_0304)*IT_0326 + IT_0304
      *IT_0329 + IT_0201*IT_0333 + conj(IT_0233)*IT_0337 + IT_0233*IT_0338 +
       conj(IT_0259)*((-24)*IT_0142 + IT_0342 + IT_0343) + conj(IT_0125)*IT_0384
       + IT_0125*IT_0385 + IT_0259*((-24)*conj(IT_0142) + IT_0345 + IT_0386) +
       IT_0253*IT_0387;
    const ccomplex_t IT_0389 = m_d*s_24*m_C_1;
    const ccomplex_t IT_0390 = 12*conj(IT_0201);
    const ccomplex_t IT_0391 = (-12)*conj(IT_0233);
    const ccomplex_t IT_0392 = IT_0387 + IT_0390 + IT_0391;
    const ccomplex_t IT_0393 = 24*conj(IT_0253);
    const ccomplex_t IT_0394 = 24*conj(IT_0125)*IT_0253 + IT_0254*conj(IT_0278
      ) + IT_0236*conj(IT_0296) + conj(IT_0259)*((-24)*IT_0125 + IT_0324 +
       IT_0325) + IT_0259*((-24)*conj(IT_0125) + IT_0327 + IT_0328) + IT_0278
      *IT_0333 + conj(IT_0285)*IT_0337 + IT_0285*IT_0338 + conj(IT_0304)*IT_0344
       + conj(IT_0142)*IT_0384 + IT_0142*IT_0385 + IT_0304*(IT_0345 + IT_0386) +
       IT_0296*IT_0392 + IT_0125*IT_0393;
    const ccomplex_t IT_0395 = -conj(IT_0285);
    const ccomplex_t IT_0396 = -conj(IT_0032)*IT_0201 + (-2)*IT_0032*(conj
      (IT_0142) + 0.5*conj(IT_0201)) + conj(IT_0051)*(IT_0278 + -IT_0285) + (
      -0.166666666666667)*IT_0142*IT_0347 + (-0.166666666666667)*conj(IT_0233)
      *IT_0348 + (-0.166666666666667)*IT_0233*IT_0349 + (-0.166666666666667)
      *conj(IT_0125)*IT_0352 + (-0.166666666666667)*IT_0125*IT_0353 + IT_0051*
      (conj(IT_0278) + IT_0395);
    const ccomplex_t IT_0397 = 12*IT_0154;
    const ccomplex_t IT_0398 = 24*IT_0316;
    const ccomplex_t IT_0399 = 24*IT_0341;
    const ccomplex_t IT_0400 = 12*IT_0262;
    const ccomplex_t IT_0401 = IT_0397 + IT_0398 + IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = 48*IT_0314;
    const ccomplex_t IT_0403 = IT_0370 + IT_0371 + IT_0372 + IT_0402;
    const ccomplex_t IT_0404 = -s_14 + -s_23;
    const ccomplex_t IT_0405 = s_13 + s_24;
    const ccomplex_t IT_0406 = IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = IT_0104*IT_0406;
    const ccomplex_t IT_0408 = (-8)*IT_0407;
    const ccomplex_t IT_0409 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0410 = IT_0405 + IT_0409;
    const ccomplex_t IT_0411 = IT_0104*IT_0410;
    const ccomplex_t IT_0412 = (-12)*IT_0154;
    const ccomplex_t IT_0413 = (-24)*IT_0316;
    const ccomplex_t IT_0414 = (-24)*IT_0341;
    const ccomplex_t IT_0415 = (-12)*IT_0262;
    const ccomplex_t IT_0416 = IT_0412 + IT_0413 + IT_0414 + IT_0415;
    const ccomplex_t IT_0417 = -(conj(IT_0125) + conj(IT_0142))*IT_0374 + (96
      *IT_0154 + 96*IT_0262 + 96*IT_0314 + 192*IT_0316 + 96*IT_0341)*IT_0381 + 
      (IT_0201 + -conj(IT_0201) + IT_0278 + -conj(IT_0278))*IT_0401 + -(conj
      (IT_0151) + conj(IT_0313))*IT_0403 + -(IT_0127 + IT_0144 + IT_0363 +
       IT_0365)*IT_0408 + -((-24)*IT_0253 + 24*IT_0259 + (-24)*conj(IT_0259) + (
      -24)*IT_0296 + 24*conj(IT_0296) + 24*IT_0304 + IT_0385 + IT_0393)*IT_0411 
      + (IT_0233 + -conj(IT_0233) + IT_0285 + IT_0395)*IT_0416;
    const ccomplex_t IT_0418 = s_12*s_34;
    const ccomplex_t IT_0419 = -IT_0314;
    const ccomplex_t IT_0420 = IT_0262 + IT_0341 + IT_0418 + IT_0419;
    const ccomplex_t IT_0421 = 0.166666666666667*IT_0259;
    const ccomplex_t IT_0422 = 0.166666666666667*conj(IT_0259);
    const ccomplex_t IT_0423 = -IT_0418;
    const ccomplex_t IT_0424 = -IT_0341;
    const ccomplex_t IT_0425 = IT_0262 + IT_0314 + IT_0423 + IT_0424;
    const ccomplex_t IT_0426 = (-12)*conj(IT_0259);
    const ccomplex_t IT_0427 = IT_0333 + IT_0426;
    const ccomplex_t IT_0428 = 12*conj(IT_0296);
    const ccomplex_t IT_0429 = (-12)*conj(IT_0304);
    const ccomplex_t IT_0430 = IT_0428 + IT_0429;
    const ccomplex_t IT_0431 = 12*conj(IT_0051);
    const ccomplex_t IT_0432 = (-12)*conj(IT_0032);
    const ccomplex_t IT_0433 = IT_0382*IT_0403;
    const ccomplex_t IT_0434 = IT_0151*(IT_0108*IT_0127 + IT_0144*IT_0148 + 48
      *conj(IT_0151)*IT_0314 + IT_0262*IT_0323 + IT_0154*IT_0392 + IT_0237
      *IT_0427 + IT_0289*IT_0430 + IT_0308*IT_0431 + IT_0153*IT_0432 + IT_0433);
    const ccomplex_t IT_0435 = (24*IT_0032*conj(IT_0032) + 24*IT_0051*conj
      (IT_0051))*IT_0054 + ((-24)*IT_0081*conj(IT_0081) + (-24)*IT_0100*conj
      (IT_0100))*IT_0103 + IT_0108*(conj(IT_0125)*IT_0126 + IT_0125*IT_0127 +
       conj(IT_0142)*IT_0143 + IT_0142*IT_0144) + (IT_0127*IT_0142 + IT_0126
      *conj(IT_0142) + conj(IT_0125)*IT_0143 + IT_0125*IT_0144)*IT_0148 + conj
      (IT_0151)*(IT_0108*IT_0126 + IT_0143*IT_0148 + (-12)*IT_0032*IT_0153 +
       IT_0154*IT_0236 + IT_0237*IT_0261 + IT_0262*IT_0288 + IT_0289*IT_0306 +
       IT_0308*IT_0309) + conj(IT_0313)*(IT_0108*IT_0143 + IT_0126*IT_0148 +
       IT_0236*IT_0262 + IT_0154*IT_0288 + IT_0261*IT_0289 + IT_0237*IT_0306 + (
      -12)*IT_0032*IT_0308 + IT_0153*IT_0309 + IT_0314*IT_0315) + IT_0316
      *IT_0340 + IT_0341*IT_0346 + 6*IT_0153*(conj(IT_0032)*IT_0278 + 2*IT_0032*
      (conj(IT_0125) + 1./2*conj(IT_0278)) + 0.166666666666667*IT_0125*IT_0347 +
       0.166666666666667*conj(IT_0285)*IT_0348 + 0.166666666666667*IT_0285
      *IT_0349) + 6*IT_0350*IT_0351 + 6*IT_0153*(conj(IT_0051)*IT_0233 + IT_0051
      *conj(IT_0233) + 0.166666666666667*conj(IT_0142)*IT_0352 +
       0.166666666666667*IT_0142*IT_0353 + 0.166666666666667*conj(IT_0201)
      *IT_0354 + 0.166666666666667*IT_0201*IT_0355) + IT_0356*IT_0361 + IT_0362*
      (conj(IT_0201)*IT_0363 + IT_0201*IT_0364 + conj(IT_0278)*IT_0365 + IT_0278
      *IT_0366) + (conj(IT_0233)*IT_0363 + IT_0233*IT_0364 + conj(IT_0285)
      *IT_0365 + IT_0285*IT_0366)*IT_0367 + (conj(IT_0278)*IT_0363 + IT_0278
      *IT_0364 + conj(IT_0201)*IT_0365 + IT_0201*IT_0366)*IT_0368 + (conj
      (IT_0285)*IT_0363 + IT_0285*IT_0364 + conj(IT_0233)*IT_0365 + IT_0233
      *IT_0366)*IT_0369 + (IT_0125 + IT_0142)*IT_0374*IT_0382 + IT_0383*IT_0388 
      + IT_0389*IT_0394 + (-6)*IT_0308*IT_0396 + -IT_0381*IT_0417 + 6*IT_0420*
      (conj(IT_0032)*IT_0253 + IT_0032*conj(IT_0253) + conj(IT_0051)*IT_0304 +
       IT_0051*conj(IT_0304) + 0.166666666666667*conj(IT_0296)*IT_0354 +
       0.166666666666667*IT_0296*IT_0355 + IT_0349*IT_0421 + IT_0348*IT_0422) +
       6*(conj(IT_0081)*IT_0253 + IT_0081*conj(IT_0253) + conj(IT_0100)*IT_0296 
      + IT_0100*conj(IT_0296) + 0.166666666666667*conj(IT_0304)*IT_0359 +
       0.166666666666667*IT_0304*IT_0360 + IT_0358*IT_0421 + IT_0357*IT_0422)
      *IT_0425 + IT_0313*(IT_0108*IT_0144 + IT_0127*IT_0148 + IT_0154*IT_0323 +
       IT_0262*IT_0392 + IT_0289*IT_0427 + IT_0237*IT_0430 + IT_0153*IT_0431 +
       IT_0308*IT_0432 + IT_0433) + IT_0434;
    return create_ccomplex_return(IT_0435);
}

