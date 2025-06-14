#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_h_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_h_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(alpha);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0000*conj(U_sb_00);
    const ccomplex_t IT_0002 = 0.5*IT_0001;
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0003
      *U_sb_00;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = pow(m_h, 2);
    const ccomplex_t IT_0008 = pow(m_sb_1, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0007 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = pow(m_W, -1);
    const ccomplex_t IT_0016 = m_b*e_em*mu_h*IT_0015*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = m_b*e_em*mu_h*IT_0015*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = cpow(V_tb, 2);
    const ccomplex_t IT_0021 = A_b*IT_0020*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0022 = IT_0003*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = A_b*IT_0020*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (-0.5)*IT_0017 + 0.5*IT_0019 + (-0.5)*IT_0023 +
       0.5*IT_0026;
    const ccomplex_t IT_0028 = m_b*mu_h*IT_0015;
    const ccomplex_t IT_0029 = cos(beta);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = cos(alpha);
    const ccomplex_t IT_0032 = IT_0000*IT_0029;
    const ccomplex_t IT_0033 = IT_0003*IT_0031;
    const ccomplex_t IT_0034 = IT_0032 + IT_0033;
    const ccomplex_t IT_0035 = cos(theta_W);
    const ccomplex_t IT_0036 = cpow(IT_0035, -2);
    const ccomplex_t IT_0037 = pow(m_b, 2);
    const ccomplex_t IT_0038 = 2*U_sb_00*(((-0.5)*m_W*IT_0034*(IT_0014 + 1./3
      *IT_0013*IT_0036) + IT_0000*IT_0014*IT_0015*IT_0030*IT_0037)*conj(U_sb_00)
       + 1./2*IT_0014*IT_0028*IT_0030*IT_0031*conj(U_sb_10)) + U_sb_10*(IT_0014
      *IT_0028*IT_0030*IT_0031*conj(U_sb_00) + (-0.666666666666667)*(m_W*IT_0013
      *IT_0034*IT_0036 + (-3)*IT_0000*IT_0014*IT_0015*IT_0030*IT_0037)*conj
      (U_sb_10));
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(e_em*IT_0038 +
       1.4142135623731*A_b*IT_0000*IT_0020*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0027*IT_0040;
    const ccomplex_t IT_0042 = IT_0012*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0007 + (0 + _Complex_I*1)
      *m_sb_1*Gamma_b1 + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0041*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*A_b*V_tb*V_td*IT_0003*conj
      (U_sb_10);
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0000*U_sb_10;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*A_b*V_tb*V_ts*IT_0003*conj
      (U_sb_00);
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0000*U_sb_00;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0056*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*A_b*V_tb*V_td*IT_0003*conj
      (U_sb_00);
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0000*U_sb_00;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*A_b*V_tb*V_ts*IT_0003*conj
      (U_sb_10);
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0000*U_sb_10;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = m_b*e_em*mu_h*IT_0015*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0080 = IT_0014*IT_0079;
    const ccomplex_t IT_0081 = m_b*e_em*mu_h*IT_0015*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0082 = IT_0014*IT_0081;
    const ccomplex_t IT_0083 = A_b*IT_0020*conj(U_sb_01)*U_sb_10;
    const ccomplex_t IT_0084 = IT_0003*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = A_b*IT_0020*U_sb_00*conj(U_sb_11);
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0080 + 0.5*IT_0082 + (-0.5)*IT_0085 +
       0.5*IT_0088;
    const ccomplex_t IT_0090 = 2*conj(U_sb_00)*(((-0.5)*m_W*IT_0034*(IT_0014 +
       1./3*IT_0013*IT_0036) + IT_0000*IT_0014*IT_0015*IT_0030*IT_0037)*U_sb_01 
      + 1./2*IT_0014*IT_0028*IT_0030*IT_0031*U_sb_11) + conj(U_sb_10)*(IT_0014
      *IT_0028*IT_0030*IT_0031*U_sb_01 + (-0.666666666666667)*(m_W*IT_0013
      *IT_0034*IT_0036 + (-3)*IT_0000*IT_0014*IT_0015*IT_0030*IT_0037)*U_sb_11);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(e_em*IT_0090 +
       1.4142135623731*A_b*IT_0000*IT_0020*(U_sb_01*conj(U_sb_10) + conj(U_sb_00
      )*U_sb_11));
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0089*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = m_b*e_em*mu_h*IT_0015*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0098 = IT_0014*IT_0097;
    const ccomplex_t IT_0099 = m_b*e_em*mu_h*IT_0015*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0100 = IT_0014*IT_0099;
    const ccomplex_t IT_0101 = A_b*IT_0020*U_sb_01*conj(U_sb_10);
    const ccomplex_t IT_0102 = IT_0003*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = A_b*IT_0020*conj(U_sb_00)*U_sb_11;
    const ccomplex_t IT_0105 = IT_0003*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0098 + (-0.5)*IT_0100 + 0.5*IT_0103 + (
      -0.5)*IT_0106;
    const ccomplex_t IT_0108 = 2*U_sb_00*(((-0.5)*m_W*IT_0034*(IT_0014 + 1./3
      *IT_0013*IT_0036) + IT_0000*IT_0014*IT_0015*IT_0030*IT_0037)*conj(U_sb_01)
       + 1./2*IT_0014*IT_0028*IT_0030*IT_0031*conj(U_sb_11)) + U_sb_10*(IT_0014
      *IT_0028*IT_0030*IT_0031*conj(U_sb_01) + (-0.666666666666667)*(m_W*IT_0013
      *IT_0034*IT_0036 + (-3)*IT_0000*IT_0014*IT_0015*IT_0030*IT_0037)*conj
      (U_sb_11));
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(e_em*IT_0108 +
       1.4142135623731*A_b*IT_0000*IT_0020*(conj(U_sb_01)*U_sb_10 + U_sb_00*conj
      (U_sb_11)));
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0107*IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0113 = IT_0111*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0013*IT_0036;
    const ccomplex_t IT_0116 = IT_0014 + IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*m_W*e_em*IT_0003*IT_0029*
      (IT_0000*IT_0029 + IT_0003*IT_0031)*IT_0116;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = IT_0003*IT_0014*IT_0030;
    const ccomplex_t IT_0120 = IT_0016*IT_0119;
    const ccomplex_t IT_0121 = IT_0018*IT_0119;
    const ccomplex_t IT_0122 = IT_0021*IT_0029;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = IT_0024*IT_0029;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (-0.5)*IT_0120 + 0.5*IT_0121 + 0.5*IT_0123 + (
      -0.5)*IT_0125;
    const ccomplex_t IT_0127 = IT_0118*IT_0126;
    const ccomplex_t IT_0128 = pow(m_Z, 2);
    const ccomplex_t IT_0129 = cpow((-2)*s_12 + (-2)*IT_0008 + IT_0128 + 
      -reg_prop, -1);
    const ccomplex_t IT_0130 = IT_0127*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0000*conj(U_sb_00);
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0003
      *U_sb_00;
    const ccomplex_t IT_0135 = 0.5*IT_0134;
    const ccomplex_t IT_0136 = IT_0133*IT_0135;
    const ccomplex_t IT_0137 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0138 = IT_0136*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = cpow(IT_0035, -1);
    const ccomplex_t IT_0141 = e_em*IT_0014*(IT_0000*IT_0003 + IT_0029*IT_0031
      )*IT_0140;
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0014
      *IT_0035 + 1./3*IT_0013*IT_0140)*U_sb_00*conj(U_sb_00) + IT_0013*IT_0140
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0144 = (-0.333333333333333)*IT_0143;
    const ccomplex_t IT_0145 = IT_0142*IT_0144;
    const ccomplex_t IT_0146 = s_12 + IT_0008;
    const ccomplex_t IT_0147 = IT_0145*IT_0146;
    const ccomplex_t IT_0148 = -IT_0141;
    const ccomplex_t IT_0149 = IT_0144*IT_0148;
    const ccomplex_t IT_0150 = s_13 + s_23;
    const ccomplex_t IT_0151 = IT_0149*IT_0150;
    const ccomplex_t IT_0152 = (-0.666666666666667)*IT_0143;
    const ccomplex_t IT_0153 = IT_0142*IT_0152;
    const ccomplex_t IT_0154 = IT_0146*IT_0153;
    const ccomplex_t IT_0155 = IT_0148*IT_0152;
    const ccomplex_t IT_0156 = s_23*IT_0155;
    const ccomplex_t IT_0157 = (-0.5)*IT_0151 + (-0.5)*IT_0154 + 0.5*IT_0156;
    const ccomplex_t IT_0158 = IT_0147 + IT_0157;
    const ccomplex_t IT_0159 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0160 = IT_0158*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = cpow(IT_0003, 3);
    const ccomplex_t IT_0163 = cpow(IT_0003, 2);
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*m_W*e_em*IT_0116*(IT_0031
      *IT_0162 + -IT_0029*(IT_0029*(IT_0000*IT_0029 + IT_0003*IT_0031) + 
      -IT_0000*IT_0163));
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = IT_0027*IT_0165;
    const ccomplex_t IT_0167 = cpow((-2)*s_12 + (-2)*IT_0008 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0168 = IT_0166*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*A_b*V_tb*conj(V_td)*IT_0003
      *U_sb_10;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = IT_0171*IT_0173;
    const ccomplex_t IT_0175 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0176 = IT_0174*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_ts
      *IT_0000*conj(U_sb_10);
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*A_b*V_tb*conj(V_ts)*IT_0003
      *U_sb_10;
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = IT_0179*IT_0181;
    const ccomplex_t IT_0183 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0184 = IT_0182*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0011 + IT_0043 + IT_0046 + IT_0054 + IT_0062
       + IT_0070 + IT_0078 + IT_0096 + IT_0114 + -IT_0131 + IT_0139 + (-2)
      *IT_0161 + -IT_0169 + IT_0177 + IT_0185;
    return create_ccomplex_return(3*IT_0186*conj(IT_0186));
}

