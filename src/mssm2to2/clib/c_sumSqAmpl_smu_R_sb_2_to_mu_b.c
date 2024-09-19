#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_sb_2_to_mu_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_sb_2_to_mu_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_mu, 2);
    const ccomplex_t IT_0001 = pow(m_smu_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = IT_0005*IT_0007;
    const ccomplex_t IT_0016 = m_b*N_d3*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0014 + 1.5*IT_0018);
    const ccomplex_t IT_0020 = (-0.333333333333333)*IT_0019;
    const ccomplex_t IT_0021 = IT_0009*IT_0020;
    const ccomplex_t IT_0022 = IT_0002*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0025 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0026 = IT_0011*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = m_b*N_d2*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0029 = IT_0015*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0027 + 1.5*IT_0030);
    const ccomplex_t IT_0032 = (-0.333333333333333)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = IT_0024*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0042 = IT_0011*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = m_b*N_d1*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0045 = IT_0015*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + 1.5*IT_0046);
    const ccomplex_t IT_0048 = (-0.333333333333333)*IT_0047;
    const ccomplex_t IT_0049 = IT_0040*IT_0048;
    const ccomplex_t IT_0050 = IT_0038*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0056 = IT_0011*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = m_b*N_d4*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0059 = IT_0015*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + 1.5*IT_0060);
    const ccomplex_t IT_0062 = (-0.333333333333333)*IT_0061;
    const ccomplex_t IT_0063 = IT_0054*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = -IT_0023 + -IT_0037 + -IT_0051 + -IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0070 = IT_0011*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0073 = IT_0007*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = m_b*conj(N_d4)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0076 = IT_0015*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0071 + (-3)*IT_0074 + 3
      *IT_0077);
    const ccomplex_t IT_0079 = 0.166666666666667*IT_0078;
    const ccomplex_t IT_0080 = IT_0068*IT_0079;
    const ccomplex_t IT_0081 = IT_0052*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0011;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = m_b*conj(N_d1)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0086 = IT_0015*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0089 = IT_0011*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0092 = IT_0007*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0087 + 0.333333333333333
      *IT_0090 + -IT_0093);
    const ccomplex_t IT_0095 = 3*IT_0094;
    const ccomplex_t IT_0096 = 0.166666666666667*IT_0095;
    const ccomplex_t IT_0097 = IT_0084*IT_0096;
    const ccomplex_t IT_0098 = IT_0038*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0101 = IT_0007*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = m_b*conj(N_d2)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0104 = IT_0015*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0107 = IT_0011*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0102 + -IT_0105 + (
      -0.333333333333333)*IT_0108);
    const ccomplex_t IT_0110 = (-3)*IT_0109;
    const ccomplex_t IT_0111 = 0.166666666666667*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0011;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = IT_0024*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0011;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = m_b*conj(N_d3)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0120 = IT_0015*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0123 = IT_0007*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0126 = IT_0011*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0121 + -IT_0124 +
       0.333333333333333*IT_0127);
    const ccomplex_t IT_0129 = 3*IT_0128;
    const ccomplex_t IT_0130 = 0.166666666666667*IT_0129;
    const ccomplex_t IT_0131 = IT_0118*IT_0130;
    const ccomplex_t IT_0132 = IT_0002*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = -IT_0082 + -IT_0099 + -IT_0116 + -IT_0133;
    const ccomplex_t IT_0135 = m_b*m_mu*IT_0001;
    const ccomplex_t IT_0136 = s_34*IT_0001;
    const ccomplex_t IT_0137 = (-6)*IT_0136;
    const ccomplex_t IT_0138 = s_13*s_14;
    const ccomplex_t IT_0139 = 12*IT_0138;
    const ccomplex_t IT_0140 = IT_0137 + IT_0139;
    const ccomplex_t IT_0141 = m_b*s_13;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0143 = IT_0009*IT_0130;
    const ccomplex_t IT_0144 = IT_0002*IT_0142*IT_0143;
    const ccomplex_t IT_0145 = m_mu*IT_0099;
    const ccomplex_t IT_0146 = m_mu*IT_0116;
    const ccomplex_t IT_0147 = m_mu*IT_0133;
    const ccomplex_t IT_0148 = m_mu*IT_0082;
    const ccomplex_t IT_0149 = IT_0144 + -IT_0145 + -IT_0146 + -IT_0147 + 
      -IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0151 = IT_0034*IT_0111;
    const ccomplex_t IT_0152 = IT_0024*IT_0150*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0154 = IT_0040*IT_0096;
    const ccomplex_t IT_0155 = IT_0038*IT_0153*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0157 = IT_0054*IT_0079;
    const ccomplex_t IT_0158 = IT_0052*IT_0156*IT_0157;
    const ccomplex_t IT_0159 = -IT_0152 + -IT_0155 + -IT_0158;
    const ccomplex_t IT_0160 = IT_0149 + -IT_0159;
    const ccomplex_t IT_0161 = -conj(IT_0159);
    const ccomplex_t IT_0162 = conj(IT_0149) + IT_0161;
    const ccomplex_t IT_0163 = IT_0048*IT_0084;
    const ccomplex_t IT_0164 = IT_0038*IT_0153*IT_0163;
    const ccomplex_t IT_0165 = m_mu*IT_0023;
    const ccomplex_t IT_0166 = IT_0062*IT_0068;
    const ccomplex_t IT_0167 = IT_0052*IT_0156*IT_0166;
    const ccomplex_t IT_0168 = m_mu*IT_0051;
    const ccomplex_t IT_0169 = m_mu*IT_0037;
    const ccomplex_t IT_0170 = m_mu*IT_0065;
    const ccomplex_t IT_0171 = IT_0164 + -IT_0165 + IT_0167 + -IT_0168 + 
      -IT_0169 + -IT_0170;
    const ccomplex_t IT_0172 = 6*IT_0171;
    const ccomplex_t IT_0173 = IT_0020*IT_0118;
    const ccomplex_t IT_0174 = IT_0002*IT_0142*IT_0173;
    const ccomplex_t IT_0175 = IT_0032*IT_0113;
    const ccomplex_t IT_0176 = IT_0024*IT_0150*IT_0175;
    const ccomplex_t IT_0177 = -IT_0174 + -IT_0176;
    const ccomplex_t IT_0178 = (-6)*IT_0177;
    const ccomplex_t IT_0179 = IT_0172 + IT_0178;
    const ccomplex_t IT_0180 = 6*conj(IT_0171);
    const ccomplex_t IT_0181 = (-6)*conj(IT_0177);
    const ccomplex_t IT_0182 = IT_0180 + IT_0181;
    const ccomplex_t IT_0183 = m_mu*s_14;
    const ccomplex_t IT_0184 = 6*IT_0177;
    const ccomplex_t IT_0185 = (-6)*IT_0171;
    const ccomplex_t IT_0186 = IT_0184 + IT_0185;
    const ccomplex_t IT_0187 = 6*conj(IT_0177);
    const ccomplex_t IT_0188 = (-6)*conj(IT_0171);
    const ccomplex_t IT_0189 = IT_0187 + IT_0188;
    const ccomplex_t IT_0190 = m_b*m_mu;
    const ccomplex_t IT_0191 = 6*s_34;
    const ccomplex_t IT_0192 = (-6)*s_34;
    const ccomplex_t IT_0193 = (-6)*conj(IT_0066)*(IT_0134*IT_0135 + (
      -0.166666666666667)*IT_0066*IT_0140) + (-6)*conj(IT_0134)*(IT_0066*IT_0135
       + (-0.166666666666667)*IT_0134*IT_0140) + 6*IT_0141*(conj(IT_0066)
      *IT_0160 + IT_0066*IT_0162 + 0.166666666666667*conj(IT_0134)*IT_0179 +
       0.166666666666667*IT_0134*IT_0182) + (-6)*IT_0183*(conj(IT_0134)*IT_0160 
      + IT_0134*IT_0162 + (-0.166666666666667)*conj(IT_0066)*IT_0186 + (
      -0.166666666666667)*IT_0066*IT_0189) + (conj(IT_0159)*IT_0179 + IT_0159
      *IT_0182 + conj(IT_0149)*IT_0186 + IT_0149*IT_0189)*IT_0190 + (IT_0149
      *conj(IT_0149) + IT_0159*conj(IT_0159) + IT_0171*conj(IT_0171) + IT_0177
      *conj(IT_0177))*IT_0191 + (conj(IT_0149)*IT_0159 + IT_0149*conj(IT_0159) +
       conj(IT_0171)*IT_0177 + IT_0171*conj(IT_0177))*IT_0192;
    return create_ccomplex_return(IT_0193);
}

