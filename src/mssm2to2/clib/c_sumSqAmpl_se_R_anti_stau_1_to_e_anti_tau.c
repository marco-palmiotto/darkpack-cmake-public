#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_stau_1_to_e_anti_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_stau_1_to_e_anti_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
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
    const creal_t m_tau = param->m_tau;
    const creal_t m_se_R = param->m_se_R;
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
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = N_B1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = N_W1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = IT_0002*IT_0004;
    const ccomplex_t IT_0016 = N_d1*e_em*m_tau*IT_0000*conj(U_stau_10);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0011 + IT_0014 + 
      -IT_0018);
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0006*IT_0020;
    const ccomplex_t IT_0022 = pow(m_e, 2);
    const ccomplex_t IT_0023 = pow(m_se_R, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = N_B2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0030 = IT_0008*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = N_W2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0033 = IT_0004*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = N_d2*e_em*m_tau*IT_0000*conj(U_stau_10);
    const ccomplex_t IT_0036 = IT_0015*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0031 + IT_0034 + 
      -IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = IT_0028*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = N_B3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0047 = IT_0008*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = N_W3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0050 = IT_0004*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = N_d3*e_em*m_tau*IT_0000*conj(U_stau_10);
    const ccomplex_t IT_0053 = IT_0015*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0048 + IT_0051 + 
      -IT_0054);
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0045*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = N_B4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0064 = IT_0008*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = N_W4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0067 = IT_0004*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = N_d4*e_em*m_tau*IT_0000*conj(U_stau_10);
    const ccomplex_t IT_0070 = IT_0015*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0065 + IT_0068 + 
      -IT_0071);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = IT_0062*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0026 + -IT_0043 + -IT_0060 + -IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0008;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = IT_0020*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0083 = IT_0024*IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0008;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0039*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0088 = IT_0041*IT_0086*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0008;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = IT_0056*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0093 = IT_0058*IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0008;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0073*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0098 = IT_0075*IT_0096*IT_0097;
    const ccomplex_t IT_0099 = m_e*IT_0026;
    const ccomplex_t IT_0100 = m_e*IT_0043;
    const ccomplex_t IT_0101 = m_e*IT_0060;
    const ccomplex_t IT_0102 = m_e*IT_0077;
    const ccomplex_t IT_0103 = -IT_0083 + -IT_0088 + -IT_0093 + -IT_0098 +
       IT_0099 + IT_0100 + IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = m_e*s_14;
    const ccomplex_t IT_0105 = conj(N_B1)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0106 = IT_0008*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = conj(N_d1)*e_em*m_tau*IT_0000*conj(U_stau_00);
    const ccomplex_t IT_0109 = IT_0015*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + 0.5*IT_0110);
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0006*IT_0112;
    const ccomplex_t IT_0114 = IT_0024*IT_0082*IT_0113;
    const ccomplex_t IT_0115 = conj(N_B2)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0116 = IT_0008*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = conj(N_d2)*e_em*m_tau*IT_0000*conj(U_stau_00);
    const ccomplex_t IT_0119 = IT_0015*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0117 + 0.5*IT_0120);
    const ccomplex_t IT_0122 = -IT_0121;
    const ccomplex_t IT_0123 = IT_0028*IT_0122;
    const ccomplex_t IT_0124 = IT_0041*IT_0087*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B3)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0126 = IT_0008*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = conj(N_d3)*e_em*m_tau*IT_0000*conj(U_stau_00);
    const ccomplex_t IT_0129 = IT_0015*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 0.5*IT_0130);
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0045*IT_0132;
    const ccomplex_t IT_0134 = IT_0058*IT_0092*IT_0133;
    const ccomplex_t IT_0135 = conj(N_B4)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0136 = IT_0008*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = conj(N_d4)*e_em*m_tau*IT_0000*conj(U_stau_00);
    const ccomplex_t IT_0139 = IT_0015*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0137 + 0.5*IT_0140);
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = IT_0062*IT_0142;
    const ccomplex_t IT_0144 = IT_0075*IT_0097*IT_0143;
    const ccomplex_t IT_0145 = IT_0080*IT_0112;
    const ccomplex_t IT_0146 = IT_0024*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = m_e*IT_0147;
    const ccomplex_t IT_0149 = IT_0085*IT_0122;
    const ccomplex_t IT_0150 = IT_0041*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = m_e*IT_0151;
    const ccomplex_t IT_0153 = IT_0090*IT_0132;
    const ccomplex_t IT_0154 = IT_0058*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = m_e*IT_0155;
    const ccomplex_t IT_0157 = IT_0095*IT_0142;
    const ccomplex_t IT_0158 = IT_0075*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = m_e*IT_0159;
    const ccomplex_t IT_0161 = -IT_0114 + -IT_0124 + -IT_0134 + -IT_0144 +
       IT_0148 + IT_0152 + IT_0156 + IT_0160;
    const ccomplex_t IT_0162 = s_13*m_tau;
    const ccomplex_t IT_0163 = -IT_0147 + -IT_0151 + -IT_0155 + -IT_0159;
    const ccomplex_t IT_0164 = m_e*m_tau*IT_0023;
    const ccomplex_t IT_0165 = s_13*s_14;
    const ccomplex_t IT_0166 = s_34*IT_0023;
    const ccomplex_t IT_0167 = (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = m_e*m_tau;
    const ccomplex_t IT_0170 = 2*s_34;
    const ccomplex_t IT_0171 = (-2)*conj(IT_0103);
    const ccomplex_t IT_0172 = 2*IT_0078*(conj(IT_0103)*IT_0104 + -conj
      (IT_0161)*IT_0162 + -conj(IT_0163)*IT_0164 + 2*conj(IT_0078)*IT_0168) + 2
      *IT_0103*(conj(IT_0078)*IT_0104 + -IT_0162*conj(IT_0163) + -conj(IT_0161)
      *IT_0169 + 0.5*conj(IT_0103)*IT_0170) + 2*IT_0163*(IT_0104*conj(IT_0161) +
       -conj(IT_0078)*IT_0164 + 2*conj(IT_0163)*IT_0168 + 0.5*IT_0162*IT_0171) +
       (-2)*IT_0161*(conj(IT_0078)*IT_0162 + -IT_0104*conj(IT_0163) + (-0.5)
      *conj(IT_0161)*IT_0170 + (-0.5)*IT_0169*IT_0171);
    return create_ccomplex_return(IT_0172);
}

