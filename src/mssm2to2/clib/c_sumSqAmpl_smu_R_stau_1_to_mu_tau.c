#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_stau_1_to_mu_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_stau_1_to_mu_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_1 = param->m_stau_1;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = pow(m_W, -1);
    const ccomplex_t IT_0016 = conj(N_d4)*e_em*m_tau*IT_0015*U_stau_10;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0006 + IT_0011 + 
      -IT_0018);
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0003*IT_0020;
    const ccomplex_t IT_0022 = pow(m_mu, 2);
    const ccomplex_t IT_0023 = pow(m_smu_R, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0031 = IT_0001*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0034 = IT_0008*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = conj(N_d3)*e_em*m_tau*IT_0015*U_stau_10;
    const ccomplex_t IT_0037 = IT_0014*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0032 + IT_0035 + 
      -IT_0038);
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0029*IT_0040;
    const ccomplex_t IT_0042 = IT_0027*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0048 = IT_0001*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0051 = IT_0008*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = conj(N_d1)*e_em*m_tau*IT_0015*U_stau_10;
    const ccomplex_t IT_0054 = IT_0014*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0049 + IT_0052 + 
      -IT_0055);
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0046*IT_0057;
    const ccomplex_t IT_0059 = IT_0044*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0064 = IT_0001*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0067 = IT_0008*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = conj(N_d2)*e_em*m_tau*IT_0015*U_stau_10;
    const ccomplex_t IT_0070 = IT_0014*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0065 + IT_0068 + 
      -IT_0071);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = IT_0062*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_13 + IT_0022 + IT_0023 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0026 + -IT_0043 + -IT_0060 + -IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0082 = IT_0001*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = N_d4*e_em*m_tau*IT_0015*U_stau_00;
    const ccomplex_t IT_0085 = IT_0014*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + 0.5*IT_0086);
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0080*IT_0088;
    const ccomplex_t IT_0090 = IT_0024*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0093 = IT_0001*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = N_d3*e_em*m_tau*IT_0015*U_stau_00;
    const ccomplex_t IT_0096 = IT_0014*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + 0.5*IT_0097);
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0099*IT_0101;
    const ccomplex_t IT_0103 = IT_0027*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0108 = IT_0001*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = N_d2*e_em*m_tau*IT_0015*U_stau_00;
    const ccomplex_t IT_0111 = IT_0014*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + 0.5*IT_0112);
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0106*IT_0114;
    const ccomplex_t IT_0116 = IT_0075*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0008*IT_0013*IT_0015;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0121 = IT_0001*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = N_d1*e_em*m_tau*IT_0015*U_stau_00;
    const ccomplex_t IT_0124 = IT_0014*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0122 + 0.5*IT_0125);
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = IT_0119*IT_0127;
    const ccomplex_t IT_0129 = IT_0044*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = -IT_0091 + -IT_0104 + -IT_0117 + -IT_0130;
    const ccomplex_t IT_0132 = pow(m_stau_1, 2);
    const ccomplex_t IT_0133 = m_mu*m_tau*IT_0132;
    const ccomplex_t IT_0134 = s_23*s_24;
    const ccomplex_t IT_0135 = s_34*IT_0132;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = IT_0134 + IT_0136;
    const ccomplex_t IT_0138 = s_23*m_tau;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0140 = IT_0040*IT_0101;
    const ccomplex_t IT_0141 = IT_0027*IT_0139*IT_0140;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0144 = IT_0020*IT_0080;
    const ccomplex_t IT_0145 = IT_0024*IT_0143*IT_0144;
    const ccomplex_t IT_0146 = m_tau*IT_0117;
    const ccomplex_t IT_0147 = m_tau*IT_0091;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0149 = IT_0057*IT_0119;
    const ccomplex_t IT_0150 = IT_0044*IT_0148*IT_0149;
    const ccomplex_t IT_0151 = m_tau*IT_0104;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0153 = IT_0073*IT_0106;
    const ccomplex_t IT_0154 = IT_0075*IT_0152*IT_0153;
    const ccomplex_t IT_0155 = m_tau*IT_0130;
    const ccomplex_t IT_0156 = IT_0145 + -IT_0146 + -IT_0147 + IT_0150 + 
      -IT_0151 + IT_0154 + -IT_0155;
    const ccomplex_t IT_0157 = IT_0142 + -IT_0156;
    const ccomplex_t IT_0158 = -conj(IT_0156);
    const ccomplex_t IT_0159 = conj(IT_0142) + IT_0158;
    const ccomplex_t IT_0160 = IT_0003*IT_0088;
    const ccomplex_t IT_0161 = IT_0024*IT_0143*IT_0160;
    const ccomplex_t IT_0162 = IT_0046*IT_0127;
    const ccomplex_t IT_0163 = IT_0044*IT_0148*IT_0162;
    const ccomplex_t IT_0164 = -IT_0161 + -IT_0163;
    const ccomplex_t IT_0165 = 2*IT_0164;
    const ccomplex_t IT_0166 = IT_0062*IT_0114;
    const ccomplex_t IT_0167 = IT_0075*IT_0152*IT_0166;
    const ccomplex_t IT_0168 = IT_0029*IT_0099;
    const ccomplex_t IT_0169 = IT_0027*IT_0139*IT_0168;
    const ccomplex_t IT_0170 = m_tau*IT_0026;
    const ccomplex_t IT_0171 = m_tau*IT_0060;
    const ccomplex_t IT_0172 = m_tau*IT_0077;
    const ccomplex_t IT_0173 = m_tau*IT_0043;
    const ccomplex_t IT_0174 = IT_0167 + IT_0169 + -IT_0170 + -IT_0171 + 
      -IT_0172 + -IT_0173;
    const ccomplex_t IT_0175 = (-2)*IT_0174;
    const ccomplex_t IT_0176 = IT_0165 + IT_0175;
    const ccomplex_t IT_0177 = 2*conj(IT_0164);
    const ccomplex_t IT_0178 = (-2)*conj(IT_0174);
    const ccomplex_t IT_0179 = IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = m_mu*m_tau;
    const ccomplex_t IT_0181 = 2*IT_0174;
    const ccomplex_t IT_0182 = (-2)*IT_0164;
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = 2*conj(IT_0174);
    const ccomplex_t IT_0185 = (-2)*conj(IT_0164);
    const ccomplex_t IT_0186 = IT_0184 + IT_0185;
    const ccomplex_t IT_0187 = m_mu*s_24;
    const ccomplex_t IT_0188 = 2*s_34;
    const ccomplex_t IT_0189 = (-2)*s_34;
    const ccomplex_t IT_0190 = (-2)*conj(IT_0078)*(IT_0131*IT_0133 + (-2)
      *IT_0078*IT_0137) + (-2)*conj(IT_0131)*(IT_0078*IT_0133 + (-2)*IT_0131
      *IT_0137) + 2*IT_0138*(conj(IT_0131)*IT_0157 + IT_0131*IT_0159 + 0.5*conj
      (IT_0078)*IT_0176 + 0.5*IT_0078*IT_0179) + IT_0180*(conj(IT_0156)*IT_0176 
      + IT_0156*IT_0179 + conj(IT_0142)*IT_0183 + IT_0142*IT_0186) + (-2)*(conj
      (IT_0078)*IT_0157 + IT_0078*IT_0159 + (-0.5)*conj(IT_0131)*IT_0183 + (-0.5
      )*IT_0131*IT_0186)*IT_0187 + (IT_0142*conj(IT_0142) + IT_0156*conj(IT_0156
      ) + IT_0164*conj(IT_0164) + IT_0174*conj(IT_0174))*IT_0188 + (conj(IT_0142
      )*IT_0156 + IT_0142*conj(IT_0156) + conj(IT_0164)*IT_0174 + IT_0164*conj
      (IT_0174))*IT_0189;
    return create_ccomplex_return(IT_0190);
}

