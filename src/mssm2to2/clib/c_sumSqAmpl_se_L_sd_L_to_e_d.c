#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_sd_L_to_e_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_sd_L_to_e_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_L = param->m_se_L;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B2)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W2)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + (-3)*IT_0009);
    const ccomplex_t IT_0011 = 0.166666666666667*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_e*N_d2*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0011*IT_0016;
    const ccomplex_t IT_0018 = pow(m_e, 2);
    const ccomplex_t IT_0019 = pow(m_se_L, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B3)*e_em;
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W3)*e_em;
    const ccomplex_t IT_0027 = IT_0006*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + (-3)*IT_0028);
    const ccomplex_t IT_0030 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B4)*e_em;
    const ccomplex_t IT_0038 = IT_0001*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W4)*e_em;
    const ccomplex_t IT_0041 = IT_0006*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + (-3)*IT_0042);
    const ccomplex_t IT_0044 = 0.166666666666667*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_13 + IT_0018 + IT_0019 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0052 = conj(N_B1)*e_em;
    const ccomplex_t IT_0053 = IT_0001*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = conj(N_W1)*e_em;
    const ccomplex_t IT_0056 = IT_0006*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + (-3)*IT_0057);
    const ccomplex_t IT_0059 = 0.166666666666667*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*m_e*N_d1*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = IT_0051*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = -IT_0022 + -IT_0036 + -IT_0050 + -IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0004 + IT_0009);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = IT_0020*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0025 + IT_0028);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = IT_0034*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0039 + IT_0042);
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = IT_0048*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0054 + IT_0057);
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = IT_0051*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = -IT_0072 + -IT_0079 + -IT_0086 + -IT_0093;
    const ccomplex_t IT_0095 = m_d*m_e*IT_0019;
    const ccomplex_t IT_0096 = s_34*IT_0019;
    const ccomplex_t IT_0097 = (-6)*IT_0096;
    const ccomplex_t IT_0098 = s_13*s_14;
    const ccomplex_t IT_0099 = 12*IT_0098;
    const ccomplex_t IT_0100 = IT_0097 + IT_0099;
    const ccomplex_t IT_0101 = m_d*m_e;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0103 = IT_0011*IT_0067;
    const ccomplex_t IT_0104 = IT_0020*IT_0102*IT_0103;
    const ccomplex_t IT_0105 = m_e*IT_0022;
    const ccomplex_t IT_0106 = m_e*IT_0036;
    const ccomplex_t IT_0107 = m_e*IT_0050;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0109 = IT_0030*IT_0074;
    const ccomplex_t IT_0110 = IT_0034*IT_0108*IT_0109;
    const ccomplex_t IT_0111 = m_e*IT_0064;
    const ccomplex_t IT_0112 = IT_0104 + -IT_0105 + -IT_0106 + -IT_0107 +
       IT_0110 + -IT_0111;
    const ccomplex_t IT_0113 = IT_0016*IT_0069;
    const ccomplex_t IT_0114 = IT_0020*IT_0102*IT_0113;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = 6*IT_0115;
    const ccomplex_t IT_0117 = IT_0032*IT_0076;
    const ccomplex_t IT_0118 = IT_0034*IT_0108*IT_0117;
    const ccomplex_t IT_0119 = IT_0046*IT_0083;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0121 = IT_0048*IT_0119*IT_0120;
    const ccomplex_t IT_0122 = m_e*IT_0072;
    const ccomplex_t IT_0123 = m_e*IT_0079;
    const ccomplex_t IT_0124 = m_e*IT_0086;
    const ccomplex_t IT_0125 = IT_0061*IT_0090;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0127 = IT_0051*IT_0125*IT_0126;
    const ccomplex_t IT_0128 = m_e*IT_0093;
    const ccomplex_t IT_0129 = IT_0118 + IT_0121 + -IT_0122 + -IT_0123 + 
      -IT_0124 + IT_0127 + -IT_0128;
    const ccomplex_t IT_0130 = (-6)*IT_0129;
    const ccomplex_t IT_0131 = IT_0116 + IT_0130;
    const ccomplex_t IT_0132 = 6*conj(IT_0115);
    const ccomplex_t IT_0133 = (-6)*conj(IT_0129);
    const ccomplex_t IT_0134 = IT_0132 + IT_0133;
    const ccomplex_t IT_0135 = IT_0059*IT_0088;
    const ccomplex_t IT_0136 = IT_0051*IT_0126*IT_0135;
    const ccomplex_t IT_0137 = IT_0044*IT_0081;
    const ccomplex_t IT_0138 = IT_0048*IT_0120*IT_0137;
    const ccomplex_t IT_0139 = -IT_0136 + -IT_0138;
    const ccomplex_t IT_0140 = 6*IT_0129;
    const ccomplex_t IT_0141 = (-6)*IT_0115;
    const ccomplex_t IT_0142 = IT_0140 + IT_0141;
    const ccomplex_t IT_0143 = 6*conj(IT_0129);
    const ccomplex_t IT_0144 = (-6)*conj(IT_0115);
    const ccomplex_t IT_0145 = IT_0143 + IT_0144;
    const ccomplex_t IT_0146 = m_d*s_13;
    const ccomplex_t IT_0147 = IT_0112 + -IT_0139;
    const ccomplex_t IT_0148 = -conj(IT_0139);
    const ccomplex_t IT_0149 = conj(IT_0112) + IT_0148;
    const ccomplex_t IT_0150 = 6*s_34;
    const ccomplex_t IT_0151 = (-6)*s_34;
    const ccomplex_t IT_0152 = (-0.166666666666667)*IT_0094;
    const ccomplex_t IT_0153 = m_e*s_14;
    const ccomplex_t IT_0154 = (-6)*conj(IT_0065)*(IT_0094*IT_0095 + (
      -0.166666666666667)*IT_0065*IT_0100) + IT_0101*(conj(IT_0112)*IT_0131 +
       IT_0112*IT_0134 + conj(IT_0139)*IT_0142 + IT_0139*IT_0145) + 6*IT_0146*
      (0.166666666666667*conj(IT_0065)*IT_0142 + 0.166666666666667*IT_0065
      *IT_0145 + conj(IT_0094)*IT_0147 + IT_0094*IT_0149) + (IT_0112*conj
      (IT_0112) + IT_0115*conj(IT_0115) + IT_0129*conj(IT_0129) + IT_0139*conj
      (IT_0139))*IT_0150 + (conj(IT_0115)*IT_0129 + IT_0115*conj(IT_0129) + conj
      (IT_0112)*IT_0139 + IT_0112*conj(IT_0139))*IT_0151 + (-6)*conj(IT_0094)*
      (IT_0065*IT_0095 + IT_0100*IT_0152) + (-6)*((-0.166666666666667)*conj
      (IT_0094)*IT_0131 + conj(IT_0065)*IT_0147 + IT_0065*IT_0149 + IT_0134
      *IT_0152)*IT_0153;
    return create_ccomplex_return(IT_0154);
}

