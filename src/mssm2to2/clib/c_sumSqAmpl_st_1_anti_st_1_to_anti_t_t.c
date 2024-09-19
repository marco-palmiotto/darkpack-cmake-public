#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_1_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_1_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
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
    const creal_t m_st_1 = param->m_st_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_t, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_st_1, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*m_sG;
    const ccomplex_t IT_0008 = IT_0002*IT_0006*IT_0007;
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = e_em*IT_0012;
    const ccomplex_t IT_0014 = cpow(IT_0011, -1);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = e_em*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + (-3)*IT_0016);
    const ccomplex_t IT_0018 = (-0.166666666666667)*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0005 + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0010*IT_0011 
      + (-3)*IT_0009*IT_0014)*U_st_00*conj(U_st_00) + IT_0010*IT_0011*U_st_10
      *conj(U_st_10));
    const ccomplex_t IT_0023 = 0.666666666666667*IT_0022;
    const ccomplex_t IT_0024 = m_t*IT_0023;
    const ccomplex_t IT_0025 = 1.33333333333333*IT_0022;
    const ccomplex_t IT_0026 = m_t*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0024 + IT_0027;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*IT_0010*IT_0011;
    const ccomplex_t IT_0031 = (-0.666666666666667)*IT_0030;
    const ccomplex_t IT_0032 = IT_0019*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0036 = IT_0000*IT_0035;
    const ccomplex_t IT_0037 = IT_0006*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = m_t*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0008 + -IT_0029 + IT_0034 + 0.5
      *IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0042 = IT_0035*IT_0041;
    const ccomplex_t IT_0043 = IT_0006*IT_0007*IT_0042;
    const ccomplex_t IT_0044 = IT_0001*IT_0041;
    const ccomplex_t IT_0045 = IT_0006*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = m_t*IT_0046;
    const ccomplex_t IT_0048 = IT_0029 + -IT_0034 + (-0.5)*IT_0043 + 0.5
      *IT_0047;
    const ccomplex_t IT_0049 = pow(m_W, -1);
    const ccomplex_t IT_0050 = sin(beta);
    const ccomplex_t IT_0051 = cpow(IT_0050, -1);
    const ccomplex_t IT_0052 = sin(alpha);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_t*e_em*IT_0014*IT_0049
      *IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = m_t*mu_h*IT_0049;
    const ccomplex_t IT_0056 = cos(alpha);
    const ccomplex_t IT_0057 = cos(beta);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = -IT_0050*IT_0052 + IT_0058;
    const ccomplex_t IT_0060 = cpow(IT_0009, -2);
    const ccomplex_t IT_0061 = (-2)*U_st_00*((IT_0003*IT_0014*IT_0049*IT_0051
      *IT_0052 + 0.5*m_W*IT_0059*(IT_0014 + -1./3*IT_0011*IT_0060))*conj(U_st_00
      ) + -1./2*IT_0014*IT_0051*IT_0055*IT_0056*conj(U_st_10)) + U_st_10*
      (IT_0014*IT_0051*IT_0055*IT_0056*conj(U_st_00) + (-1.33333333333333)*(1.5
      *IT_0003*IT_0014*IT_0049*IT_0051*IT_0052 + m_W*IT_0011*IT_0059*IT_0060)
      *conj(U_st_10));
    const ccomplex_t IT_0062 = (0 + _Complex_I*-1)*(e_em*IT_0061 + (
      -1.4142135623731)*A_t*IT_0052*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0054*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + (-2)*IT_0005 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*m_t*e_em*IT_0014*IT_0049
      *IT_0051*IT_0056;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0052*IT_0057;
    const ccomplex_t IT_0071 = IT_0050*IT_0056;
    const ccomplex_t IT_0072 = IT_0070 + IT_0071;
    const ccomplex_t IT_0073 = 2*U_st_00*((IT_0003*IT_0014*IT_0049*IT_0051
      *IT_0056 + (-0.5)*m_W*(IT_0014 + -1./3*IT_0011*IT_0060)*IT_0072)*conj
      (U_st_00) + 1./2*IT_0014*IT_0051*IT_0052*IT_0055*conj(U_st_10)) + U_st_10*
      (IT_0014*IT_0051*IT_0052*IT_0055*conj(U_st_00) + (-1.33333333333333)*((
      -1.5)*IT_0003*IT_0014*IT_0049*IT_0051*IT_0056 + m_W*IT_0011*IT_0060
      *IT_0072)*conj(U_st_10));
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(e_em*IT_0073 +
       1.4142135623731*A_t*IT_0056*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = IT_0069*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0005 + -reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = -IT_0067 + -IT_0079;
    const ccomplex_t IT_0081 = 18*s_34;
    const ccomplex_t IT_0082 = (-18)*IT_0003;
    const ccomplex_t IT_0083 = IT_0081 + IT_0082;
    const ccomplex_t IT_0084 = IT_0021*IT_0025;
    const ccomplex_t IT_0085 = (-0.5)*IT_0038 + -IT_0084;
    const ccomplex_t IT_0086 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0087 = IT_0010*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0090 = IT_0014*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = IT_0014*IT_0051;
    const ccomplex_t IT_0093 = m_t*N_u1*e_em*IT_0049*conj(U_st_10);
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0088 + 3*IT_0091 + 3
      *IT_0095);
    const ccomplex_t IT_0097 = 0.166666666666667*IT_0096;
    const ccomplex_t IT_0098 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0099 = IT_0010*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0102 = IT_0014*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = m_t*conj(N_u1)*e_em*IT_0049*U_st_10;
    const ccomplex_t IT_0105 = IT_0092*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0100 + 3*IT_0103 + 3
      *IT_0106);
    const ccomplex_t IT_0108 = 0.166666666666667*IT_0107;
    const ccomplex_t IT_0109 = IT_0097*IT_0108;
    const ccomplex_t IT_0110 = cpow((-2)*s_23 + IT_0003 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0114 = IT_0010*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0117 = IT_0014*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = m_t*N_u2*e_em*IT_0049*conj(U_st_10);
    const ccomplex_t IT_0120 = IT_0092*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0115 + 3*IT_0118 + 3
      *IT_0121);
    const ccomplex_t IT_0123 = 0.166666666666667*IT_0122;
    const ccomplex_t IT_0124 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0125 = IT_0010*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0128 = IT_0014*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = m_t*conj(N_u2)*e_em*IT_0049*U_st_10;
    const ccomplex_t IT_0131 = IT_0092*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0126 + 3*IT_0129 + 3
      *IT_0132);
    const ccomplex_t IT_0134 = 0.166666666666667*IT_0133;
    const ccomplex_t IT_0135 = IT_0123*IT_0134;
    const ccomplex_t IT_0136 = cpow((-2)*s_23 + IT_0003 + IT_0005 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0140 = IT_0010*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0143 = IT_0014*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = m_t*N_u3*e_em*IT_0049*conj(U_st_10);
    const ccomplex_t IT_0146 = IT_0092*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0141 + 3*IT_0144 + 3
      *IT_0147);
    const ccomplex_t IT_0149 = 0.166666666666667*IT_0148;
    const ccomplex_t IT_0150 = conj(N_B3)*e_em*U_st_00;
    const ccomplex_t IT_0151 = IT_0010*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = conj(N_W3)*e_em*U_st_00;
    const ccomplex_t IT_0154 = IT_0014*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = m_t*conj(N_u3)*e_em*IT_0049*U_st_10;
    const ccomplex_t IT_0157 = IT_0092*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0152 + 3*IT_0155 + 3
      *IT_0158);
    const ccomplex_t IT_0160 = 0.166666666666667*IT_0159;
    const ccomplex_t IT_0161 = IT_0149*IT_0160;
    const ccomplex_t IT_0162 = cpow((-2)*s_23 + IT_0003 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0166 = IT_0010*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0169 = IT_0014*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = m_t*N_u4*e_em*IT_0049*conj(U_st_10);
    const ccomplex_t IT_0172 = IT_0092*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0167 + 3*IT_0170 + 3
      *IT_0173);
    const ccomplex_t IT_0175 = 0.166666666666667*IT_0174;
    const ccomplex_t IT_0176 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0177 = IT_0010*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0180 = IT_0014*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = m_t*conj(N_u4)*e_em*IT_0049*U_st_10;
    const ccomplex_t IT_0183 = IT_0092*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0178 + 3*IT_0181 + 3
      *IT_0184);
    const ccomplex_t IT_0186 = 0.166666666666667*IT_0185;
    const ccomplex_t IT_0187 = IT_0175*IT_0186;
    const ccomplex_t IT_0188 = cpow((-2)*s_23 + IT_0003 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0189 = IT_0187*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = 0.166666666666667*IT_0038 + -IT_0112 + -IT_0138
       + -IT_0164 + -IT_0190;
    const ccomplex_t IT_0192 = IT_0025*IT_0033;
    const ccomplex_t IT_0193 = (-0.5)*IT_0046 + -IT_0192;
    const ccomplex_t IT_0194 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0195 = IT_0010*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = m_t*conj(N_u1)*e_em*IT_0049*conj(U_st_00);
    const ccomplex_t IT_0198 = IT_0092*IT_0197;
    const ccomplex_t IT_0199 = 1.4142135623731*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*(IT_0196 + (-0.75)*IT_0199);
    const ccomplex_t IT_0201 = 0.666666666666667*IT_0200;
    const ccomplex_t IT_0202 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0203 = IT_0010*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = m_t*N_u1*e_em*IT_0049*U_st_00;
    const ccomplex_t IT_0206 = IT_0092*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0204 + (-0.75)*IT_0207);
    const ccomplex_t IT_0209 = 0.666666666666667*IT_0208;
    const ccomplex_t IT_0210 = IT_0201*IT_0209;
    const ccomplex_t IT_0211 = IT_0110*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0214 = IT_0010*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = m_t*conj(N_u2)*e_em*IT_0049*conj(U_st_00);
    const ccomplex_t IT_0217 = IT_0092*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*(IT_0215 + (-0.75)*IT_0218);
    const ccomplex_t IT_0220 = 0.666666666666667*IT_0219;
    const ccomplex_t IT_0221 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0222 = IT_0010*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = m_t*N_u2*e_em*IT_0049*U_st_00;
    const ccomplex_t IT_0225 = IT_0092*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + (-0.75)*IT_0226);
    const ccomplex_t IT_0228 = 0.666666666666667*IT_0227;
    const ccomplex_t IT_0229 = IT_0220*IT_0228;
    const ccomplex_t IT_0230 = IT_0136*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0233 = IT_0010*IT_0232;
    const ccomplex_t IT_0234 = 1.4142135623731*IT_0233;
    const ccomplex_t IT_0235 = m_t*conj(N_u3)*e_em*IT_0049*conj(U_st_00);
    const ccomplex_t IT_0236 = IT_0092*IT_0235;
    const ccomplex_t IT_0237 = 1.4142135623731*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*(IT_0234 + (-0.75)*IT_0237);
    const ccomplex_t IT_0239 = 0.666666666666667*IT_0238;
    const ccomplex_t IT_0240 = N_B3*e_em*U_st_10;
    const ccomplex_t IT_0241 = IT_0010*IT_0240;
    const ccomplex_t IT_0242 = 1.4142135623731*IT_0241;
    const ccomplex_t IT_0243 = m_t*N_u3*e_em*IT_0049*U_st_00;
    const ccomplex_t IT_0244 = IT_0092*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*(IT_0242 + (-0.75)*IT_0245);
    const ccomplex_t IT_0247 = 0.666666666666667*IT_0246;
    const ccomplex_t IT_0248 = IT_0239*IT_0247;
    const ccomplex_t IT_0249 = IT_0162*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0252 = IT_0010*IT_0251;
    const ccomplex_t IT_0253 = 1.4142135623731*IT_0252;
    const ccomplex_t IT_0254 = m_t*conj(N_u4)*e_em*IT_0049*conj(U_st_00);
    const ccomplex_t IT_0255 = IT_0092*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*(IT_0253 + (-0.75)*IT_0256);
    const ccomplex_t IT_0258 = 0.666666666666667*IT_0257;
    const ccomplex_t IT_0259 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0260 = IT_0010*IT_0259;
    const ccomplex_t IT_0261 = 1.4142135623731*IT_0260;
    const ccomplex_t IT_0262 = m_t*N_u4*e_em*IT_0049*U_st_00;
    const ccomplex_t IT_0263 = IT_0092*IT_0262;
    const ccomplex_t IT_0264 = 1.4142135623731*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*(IT_0261 + (-0.75)*IT_0264);
    const ccomplex_t IT_0266 = 0.666666666666667*IT_0265;
    const ccomplex_t IT_0267 = IT_0258*IT_0266;
    const ccomplex_t IT_0268 = IT_0188*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = 0.166666666666667*IT_0046 + -IT_0212 + -IT_0231
       + -IT_0250 + -IT_0269;
    const ccomplex_t IT_0271 = s_34*IT_0005;
    const ccomplex_t IT_0272 = (-6)*IT_0271;
    const ccomplex_t IT_0273 = s_13*s_14;
    const ccomplex_t IT_0274 = 12*IT_0273;
    const ccomplex_t IT_0275 = IT_0272 + IT_0274;
    const ccomplex_t IT_0276 = (-18)*IT_0271;
    const ccomplex_t IT_0277 = 36*IT_0273;
    const ccomplex_t IT_0278 = IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0280 = 0.666666666666667*IT_0279;
    const ccomplex_t IT_0281 = (-1.33333333333333)*IT_0279;
    const ccomplex_t IT_0282 = IT_0280*IT_0281;
    const ccomplex_t IT_0283 = cpow(s_12 + IT_0005 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0284 = IT_0282*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0287 = (-2)*IT_0286;
    const ccomplex_t IT_0288 = IT_0283*IT_0286;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0287*IT_0289;
    const ccomplex_t IT_0291 = 0.5*IT_0285 + (-0.0833333333333333)*IT_0290;
    const ccomplex_t IT_0292 = cpow(IT_0291, 2);
    const ccomplex_t IT_0293 = (-36)*IT_0271;
    const ccomplex_t IT_0294 = IT_0003*IT_0005;
    const ccomplex_t IT_0295 = (-36)*IT_0294;
    const ccomplex_t IT_0296 = 72*IT_0273;
    const ccomplex_t IT_0297 = IT_0293 + IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = 0.25*IT_0290;
    const ccomplex_t IT_0299 = (-18)*IT_0040;
    const ccomplex_t IT_0300 = IT_0097*IT_0209;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0302 = IT_0110*IT_0300*IT_0301;
    const ccomplex_t IT_0303 = IT_0123*IT_0228;
    const ccomplex_t IT_0304 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0305 = IT_0136*IT_0303*IT_0304;
    const ccomplex_t IT_0306 = IT_0149*IT_0247;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0308 = IT_0162*IT_0306*IT_0307;
    const ccomplex_t IT_0309 = IT_0175*IT_0266;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0311 = IT_0188*IT_0309*IT_0310;
    const ccomplex_t IT_0312 = m_t*IT_0112;
    const ccomplex_t IT_0313 = m_t*IT_0138;
    const ccomplex_t IT_0314 = m_t*IT_0164;
    const ccomplex_t IT_0315 = m_t*IT_0190;
    const ccomplex_t IT_0316 = 0.166666666666667*IT_0008 + (-0.166666666666667
      )*IT_0039 + -IT_0302 + -IT_0305 + -IT_0308 + -IT_0311 + IT_0312 + IT_0313 
      + IT_0314 + IT_0315;
    const ccomplex_t IT_0317 = (-6)*IT_0316;
    const ccomplex_t IT_0318 = IT_0299 + IT_0317;
    const ccomplex_t IT_0319 = (-18)*conj(IT_0040);
    const ccomplex_t IT_0320 = (-6)*conj(IT_0316);
    const ccomplex_t IT_0321 = IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = IT_0108*IT_0201;
    const ccomplex_t IT_0323 = IT_0110*IT_0301*IT_0322;
    const ccomplex_t IT_0324 = IT_0134*IT_0220;
    const ccomplex_t IT_0325 = IT_0136*IT_0304*IT_0324;
    const ccomplex_t IT_0326 = IT_0160*IT_0239;
    const ccomplex_t IT_0327 = IT_0162*IT_0307*IT_0326;
    const ccomplex_t IT_0328 = m_t*IT_0212;
    const ccomplex_t IT_0329 = m_t*IT_0231;
    const ccomplex_t IT_0330 = m_t*IT_0250;
    const ccomplex_t IT_0331 = IT_0186*IT_0258;
    const ccomplex_t IT_0332 = IT_0188*IT_0310*IT_0331;
    const ccomplex_t IT_0333 = m_t*IT_0269;
    const ccomplex_t IT_0334 = 0.166666666666667*IT_0043 + (-0.166666666666667
      )*IT_0047 + -IT_0323 + -IT_0325 + -IT_0327 + IT_0328 + IT_0329 + IT_0330 +
       -IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = (-18)*IT_0316;
    const ccomplex_t IT_0336 = (-6)*IT_0040;
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = (-18)*conj(IT_0316);
    const ccomplex_t IT_0339 = (-6)*conj(IT_0040);
    const ccomplex_t IT_0340 = IT_0338 + IT_0339;
    const ccomplex_t IT_0341 = 6*s_34;
    const ccomplex_t IT_0342 = -IT_0003;
    const ccomplex_t IT_0343 = 6*IT_0080;
    const ccomplex_t IT_0344 = 6*conj(IT_0080);
    const ccomplex_t IT_0345 = -s_14;
    const ccomplex_t IT_0346 = s_13 + IT_0345;
    const ccomplex_t IT_0347 = m_t*IT_0346;
    const ccomplex_t IT_0348 = conj(IT_0191) + conj(IT_0270);
    const ccomplex_t IT_0349 = IT_0191 + IT_0270;
    const ccomplex_t IT_0350 = conj(IT_0085) + conj(IT_0193);
    const ccomplex_t IT_0351 = 18*IT_0080;
    const ccomplex_t IT_0352 = IT_0085 + IT_0193;
    const ccomplex_t IT_0353 = 18*conj(IT_0080);
    const ccomplex_t IT_0354 = 6*conj(IT_0040);
    const ccomplex_t IT_0355 = 6*conj(IT_0048);
    const ccomplex_t IT_0356 = 18*conj(IT_0316);
    const ccomplex_t IT_0357 = 18*conj(IT_0334);
    const ccomplex_t IT_0358 = (-18)*IT_0334;
    const ccomplex_t IT_0359 = (-6)*IT_0048;
    const ccomplex_t IT_0360 = 6*conj(IT_0316);
    const ccomplex_t IT_0361 = 6*conj(IT_0334);
    const ccomplex_t IT_0362 = 18*conj(IT_0040);
    const ccomplex_t IT_0363 = 18*conj(IT_0048);
    const ccomplex_t IT_0364 = (-18)*IT_0048;
    const ccomplex_t IT_0365 = (-6)*IT_0334;
    const ccomplex_t IT_0366 = IT_0343*IT_0348 + IT_0344*IT_0349 + IT_0350
      *IT_0351 + IT_0352*IT_0353 + IT_0298*((-12)*IT_0080 + 12*conj(IT_0080) +
       IT_0335 + IT_0336 + IT_0354 + IT_0355 + IT_0356 + IT_0357 + IT_0358 +
       IT_0359) + IT_0291*((-36)*IT_0080 + 36*conj(IT_0080) + IT_0299 + IT_0317 
      + IT_0360 + IT_0361 + IT_0362 + IT_0363 + IT_0364 + IT_0365);
    const ccomplex_t IT_0367 = (-18)*IT_0294;
    const ccomplex_t IT_0368 = IT_0276 + IT_0277 + IT_0367;
    const ccomplex_t IT_0369 = -IT_0298;
    const ccomplex_t IT_0370 = -IT_0291;
    const ccomplex_t IT_0371 = (-6)*IT_0294;
    const ccomplex_t IT_0372 = IT_0272 + IT_0274 + IT_0371;
    const ccomplex_t IT_0373 = 36*s_34;
    const ccomplex_t IT_0374 = m_t*s_13;
    const ccomplex_t IT_0375 = 6*IT_0334;
    const ccomplex_t IT_0376 = 6*conj(IT_0270);
    const ccomplex_t IT_0377 = conj(IT_0191)*(IT_0040 + 3*IT_0316) +
       0.166666666666667*conj(IT_0085)*(18*IT_0040 + 6*IT_0316) +
       0.166666666666667*IT_0191*(IT_0354 + IT_0356) + 0.166666666666667*IT_0270
      *(IT_0355 + IT_0357) + 0.166666666666667*IT_0085*(IT_0360 + IT_0362) +
       0.166666666666667*IT_0193*(IT_0361 + IT_0363) + 0.166666666666667*conj
      (IT_0193)*(18*IT_0048 + IT_0375) + 0.166666666666667*(IT_0048 + 3*IT_0334)
      *IT_0376;
    const ccomplex_t IT_0378 = m_t*e_em*IT_0014*IT_0049;
    const ccomplex_t IT_0379 = (-0.5)*IT_0378;
    const ccomplex_t IT_0380 = m_t*e_em*mu_h*IT_0049*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0381 = IT_0014*IT_0051*IT_0057;
    const ccomplex_t IT_0382 = IT_0380*IT_0381;
    const ccomplex_t IT_0383 = m_t*e_em*mu_h*IT_0049*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0384 = IT_0381*IT_0383;
    const ccomplex_t IT_0385 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0386 = IT_0050*IT_0385;
    const ccomplex_t IT_0387 = 1.4142135623731*IT_0386;
    const ccomplex_t IT_0388 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0389 = IT_0050*IT_0388;
    const ccomplex_t IT_0390 = 1.4142135623731*IT_0389;
    const ccomplex_t IT_0391 = 0.5*IT_0382 + (-0.5)*IT_0384 + (-0.5)*IT_0387 +
       0.5*IT_0390;
    const ccomplex_t IT_0392 = IT_0379*IT_0391;
    const ccomplex_t IT_0393 = pow(m_Z, 2);
    const ccomplex_t IT_0394 = cpow((-2)*s_12 + (-2)*IT_0005 + IT_0393 + 
      -reg_prop, -1);
    const ccomplex_t IT_0395 = IT_0392*IT_0394;
    const ccomplex_t IT_0396 = (0 + _Complex_I*1)*IT_0395;
    const ccomplex_t IT_0397 = IT_0014*IT_0380;
    const ccomplex_t IT_0398 = IT_0014*IT_0383;
    const ccomplex_t IT_0399 = IT_0057*IT_0385;
    const ccomplex_t IT_0400 = 1.4142135623731*IT_0399;
    const ccomplex_t IT_0401 = IT_0057*IT_0388;
    const ccomplex_t IT_0402 = 1.4142135623731*IT_0401;
    const ccomplex_t IT_0403 = (-0.5)*IT_0397 + 0.5*IT_0398 + (-0.5)*IT_0400 +
       0.5*IT_0402;
    const ccomplex_t IT_0404 = m_t*e_em*IT_0014*IT_0049*IT_0051*IT_0057;
    const ccomplex_t IT_0405 = (-0.5)*IT_0404;
    const ccomplex_t IT_0406 = IT_0403*IT_0405;
    const ccomplex_t IT_0407 = cpow((-2)*s_12 + (-2)*IT_0005 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0408 = IT_0406*IT_0407;
    const ccomplex_t IT_0409 = (0 + _Complex_I*1)*IT_0408;
    const ccomplex_t IT_0410 = -IT_0396 + -IT_0409;
    const ccomplex_t IT_0411 = 18*IT_0003;
    const ccomplex_t IT_0412 = IT_0081 + IT_0411;
    const ccomplex_t IT_0413 = (-18)*s_34 + IT_0082;
    const ccomplex_t IT_0414 = s_34 + IT_0003;
    const ccomplex_t IT_0415 = s_13 + s_14;
    const ccomplex_t IT_0416 = m_t*IT_0415;
    const ccomplex_t IT_0417 = (-6)*IT_0191;
    const ccomplex_t IT_0418 = m_t*s_14;
    const ccomplex_t IT_0419 = conj(IT_0193)*IT_0318 + IT_0193*IT_0321 +
       IT_0085*((-18)*conj(IT_0048) + (-6)*conj(IT_0334)) + conj(IT_0270)
      *IT_0337 + IT_0270*IT_0340 + conj(IT_0191)*(IT_0358 + IT_0359) + conj
      (IT_0085)*(IT_0364 + IT_0365) + (conj(IT_0048) + 3*conj(IT_0334))*IT_0417;
    const ccomplex_t IT_0420 = (-6)*conj(IT_0191);
    const ccomplex_t IT_0421 = ((conj(IT_0040) + conj(IT_0048))*IT_0080 + 
      (IT_0040 + IT_0048)*conj(IT_0080))*IT_0083 + (conj(IT_0085)*IT_0191 +
       IT_0085*conj(IT_0191) + conj(IT_0193)*IT_0270 + IT_0193*conj(IT_0270))
      *IT_0275 + (IT_0085*conj(IT_0085) + IT_0191*conj(IT_0191) + IT_0193*conj
      (IT_0193) + IT_0270*conj(IT_0270))*IT_0278 + -IT_0292*IT_0297 + (-2)
      *IT_0298*(IT_0291*((-12)*IT_0271 + 24*IT_0273 + (-12)*IT_0294) + 0.5
      *IT_0297*IT_0298) + IT_0081*(IT_0040*conj(IT_0040) + IT_0048*conj(IT_0048)
       + IT_0316*conj(IT_0316) + IT_0334*conj(IT_0334)) + IT_0003*(conj(IT_0048)
      *IT_0318 + IT_0048*IT_0321 + conj(IT_0334)*IT_0337 + IT_0334*IT_0340) + 
      (conj(IT_0040)*IT_0316 + IT_0040*conj(IT_0316) + conj(IT_0048)*IT_0334 +
       IT_0048*conj(IT_0334))*IT_0341 + (s_34 + IT_0342)*((conj(IT_0316) + conj
      (IT_0334))*IT_0343 + (IT_0316 + IT_0334)*IT_0344) + IT_0347*IT_0366 +
       IT_0368*(IT_0298*IT_0348 + IT_0291*IT_0350 + IT_0349*IT_0369 + IT_0352
      *IT_0370) + (IT_0291*IT_0348 + IT_0298*IT_0350 + IT_0352*IT_0369 + IT_0349
      *IT_0370)*IT_0372 + IT_0080*conj(IT_0080)*((-36)*IT_0003 + IT_0373) + 6
      *IT_0374*IT_0377 + conj(IT_0410)*(IT_0048*IT_0412 + IT_0040*IT_0413 + 
      (IT_0317 + IT_0375)*IT_0414 + IT_0416*((-18)*IT_0085 + 18*IT_0193 + 6
      *IT_0270 + IT_0417)) + IT_0418*IT_0419 + IT_0410*((36*IT_0003 + IT_0373)
      *conj(IT_0410) + conj(IT_0048)*IT_0412 + conj(IT_0040)*IT_0413 + (IT_0320 
      + IT_0361)*IT_0414 + IT_0416*((-18)*conj(IT_0085) + 18*conj(IT_0193) +
       IT_0376 + IT_0420)) + IT_0294*(conj(IT_0085)*((-18)*IT_0193 + (-6)
      *IT_0270) + IT_0085*((-18)*conj(IT_0193) + (-6)*conj(IT_0270)) + (conj
      (IT_0193) + 3*conj(IT_0270))*IT_0417 + (IT_0193 + 3*IT_0270)*IT_0420);
    return create_ccomplex_return(IT_0421);
}

