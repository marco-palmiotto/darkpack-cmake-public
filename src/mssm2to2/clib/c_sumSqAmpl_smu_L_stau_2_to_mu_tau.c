#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_stau_2_to_mu_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_stau_2_to_mu_tau(
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
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_2 = param->m_stau_2;
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
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B4)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W4)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0006*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = N_d4*e_em*m_tau*IT_0018*U_stau_01;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0014 + 0.5*IT_0021);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0011*IT_0023;
    const ccomplex_t IT_0025 = pow(m_mu, 2);
    const ccomplex_t IT_0026 = pow(m_smu_L, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0025 + IT_0026 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0024*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0025 + IT_0026 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0031 = conj(N_B1)*e_em;
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W1)*e_em;
    const ccomplex_t IT_0035 = IT_0006*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0040 = IT_0001*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = N_d1*e_em*m_tau*IT_0018*U_stau_01;
    const ccomplex_t IT_0043 = IT_0017*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0041 + 0.5*IT_0044);
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = IT_0038*IT_0046;
    const ccomplex_t IT_0048 = IT_0030*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0025 + IT_0026 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0051 = conj(N_B2)*e_em;
    const ccomplex_t IT_0052 = IT_0001*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W2)*e_em;
    const ccomplex_t IT_0055 = IT_0006*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + IT_0056);
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0060 = IT_0001*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = N_d2*e_em*m_tau*IT_0018*U_stau_01;
    const ccomplex_t IT_0063 = IT_0017*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + 0.5*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0058*IT_0066;
    const ccomplex_t IT_0068 = IT_0050*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_13 + IT_0025 + IT_0026 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0071 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0072 = IT_0001*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = N_d3*e_em*m_tau*IT_0018*U_stau_01;
    const ccomplex_t IT_0075 = IT_0017*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + 0.5*IT_0076);
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = conj(N_B3)*e_em;
    const ccomplex_t IT_0080 = IT_0001*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W3)*e_em;
    const ccomplex_t IT_0083 = IT_0006*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + IT_0084);
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0078*IT_0086;
    const ccomplex_t IT_0088 = IT_0070*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = -IT_0029 + -IT_0049 + -IT_0069 + -IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*N_d3*e_em*m_mu
      *IT_0006*IT_0016*IT_0018;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0094 = IT_0001*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0097 = IT_0006*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = conj(N_d3)*e_em*m_tau*IT_0018*U_stau_11;
    const ccomplex_t IT_0100 = IT_0017*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0095 + IT_0098 + 
      -IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0092*IT_0103;
    const ccomplex_t IT_0105 = IT_0070*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0006*IT_0016*IT_0018;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0110 = IT_0001*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0113 = IT_0006*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = conj(N_d4)*e_em*m_tau*IT_0018*U_stau_11;
    const ccomplex_t IT_0116 = IT_0017*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0111 + IT_0114 + 
      -IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = IT_0108*IT_0119;
    const ccomplex_t IT_0121 = IT_0027*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0006*IT_0016*IT_0018;
    const ccomplex_t IT_0124 = (-0.5)*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0126 = IT_0001*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0129 = IT_0006*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = conj(N_d1)*e_em*m_tau*IT_0018*U_stau_11;
    const ccomplex_t IT_0132 = IT_0017*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0127 + IT_0130 + 
      -IT_0133);
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = IT_0124*IT_0135;
    const ccomplex_t IT_0137 = IT_0030*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*N_d2*e_em*m_mu
      *IT_0006*IT_0016*IT_0018;
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0142 = IT_0001*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0145 = IT_0006*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = conj(N_d2)*e_em*m_tau*IT_0018*U_stau_11;
    const ccomplex_t IT_0148 = IT_0017*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0143 + IT_0146 + 
      -IT_0149);
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = IT_0140*IT_0151;
    const ccomplex_t IT_0153 = IT_0050*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = -IT_0106 + -IT_0122 + -IT_0138 + -IT_0154;
    const ccomplex_t IT_0156 = pow(m_stau_2, 2);
    const ccomplex_t IT_0157 = m_mu*m_tau*IT_0156;
    const ccomplex_t IT_0158 = s_23*s_24;
    const ccomplex_t IT_0159 = s_34*IT_0156;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = IT_0158 + IT_0160;
    const ccomplex_t IT_0162 = m_mu*s_24;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0164 = IT_0086*IT_0103;
    const ccomplex_t IT_0165 = IT_0070*IT_0163*IT_0164;
    const ccomplex_t IT_0166 = m_tau*IT_0049;
    const ccomplex_t IT_0167 = m_tau*IT_0069;
    const ccomplex_t IT_0168 = m_tau*IT_0089;
    const ccomplex_t IT_0169 = m_tau*IT_0029;
    const ccomplex_t IT_0170 = IT_0165 + -IT_0166 + -IT_0167 + -IT_0168 + 
      -IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0172 = IT_0011*IT_0119;
    const ccomplex_t IT_0173 = IT_0027*IT_0171*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0175 = IT_0038*IT_0135;
    const ccomplex_t IT_0176 = IT_0030*IT_0174*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0178 = IT_0058*IT_0151;
    const ccomplex_t IT_0179 = IT_0050*IT_0177*IT_0178;
    const ccomplex_t IT_0180 = -IT_0173 + -IT_0176 + -IT_0179;
    const ccomplex_t IT_0181 = IT_0170 + -IT_0180;
    const ccomplex_t IT_0182 = -conj(IT_0180);
    const ccomplex_t IT_0183 = conj(IT_0170) + IT_0182;
    const ccomplex_t IT_0184 = IT_0046*IT_0124;
    const ccomplex_t IT_0185 = IT_0030*IT_0174*IT_0184;
    const ccomplex_t IT_0186 = IT_0023*IT_0108;
    const ccomplex_t IT_0187 = IT_0027*IT_0171*IT_0186;
    const ccomplex_t IT_0188 = m_tau*IT_0106;
    const ccomplex_t IT_0189 = m_tau*IT_0122;
    const ccomplex_t IT_0190 = m_tau*IT_0138;
    const ccomplex_t IT_0191 = m_tau*IT_0154;
    const ccomplex_t IT_0192 = IT_0185 + IT_0187 + -IT_0188 + -IT_0189 + 
      -IT_0190 + -IT_0191;
    const ccomplex_t IT_0193 = 2*IT_0192;
    const ccomplex_t IT_0194 = IT_0078*IT_0092;
    const ccomplex_t IT_0195 = IT_0070*IT_0163*IT_0194;
    const ccomplex_t IT_0196 = IT_0066*IT_0140;
    const ccomplex_t IT_0197 = IT_0050*IT_0177*IT_0196;
    const ccomplex_t IT_0198 = -IT_0195 + -IT_0197;
    const ccomplex_t IT_0199 = (-2)*IT_0198;
    const ccomplex_t IT_0200 = IT_0193 + IT_0199;
    const ccomplex_t IT_0201 = 2*conj(IT_0192);
    const ccomplex_t IT_0202 = (-2)*conj(IT_0198);
    const ccomplex_t IT_0203 = IT_0201 + IT_0202;
    const ccomplex_t IT_0204 = m_mu*m_tau;
    const ccomplex_t IT_0205 = 2*IT_0198;
    const ccomplex_t IT_0206 = (-2)*IT_0192;
    const ccomplex_t IT_0207 = IT_0205 + IT_0206;
    const ccomplex_t IT_0208 = 2*conj(IT_0198);
    const ccomplex_t IT_0209 = (-2)*conj(IT_0192);
    const ccomplex_t IT_0210 = IT_0208 + IT_0209;
    const ccomplex_t IT_0211 = s_23*m_tau;
    const ccomplex_t IT_0212 = 2*s_34;
    const ccomplex_t IT_0213 = (-2)*s_34;
    const ccomplex_t IT_0214 = (-2)*conj(IT_0090)*(IT_0155*IT_0157 + (-2)
      *IT_0090*IT_0161) + (-2)*conj(IT_0155)*(IT_0090*IT_0157 + (-2)*IT_0155
      *IT_0161) + 2*IT_0162*(conj(IT_0155)*IT_0181 + IT_0155*IT_0183 + 0.5*conj
      (IT_0090)*IT_0200 + 0.5*IT_0090*IT_0203) + IT_0204*(conj(IT_0180)*IT_0200 
      + IT_0180*IT_0203 + conj(IT_0170)*IT_0207 + IT_0170*IT_0210) + (-2)*(conj
      (IT_0090)*IT_0181 + IT_0090*IT_0183 + (-0.5)*conj(IT_0155)*IT_0207 + (-0.5
      )*IT_0155*IT_0210)*IT_0211 + (IT_0170*conj(IT_0170) + IT_0180*conj(IT_0180
      ) + IT_0192*conj(IT_0192) + IT_0198*conj(IT_0198))*IT_0212 + (conj(IT_0170
      )*IT_0180 + IT_0170*conj(IT_0180) + conj(IT_0192)*IT_0198 + IT_0192*conj
      (IT_0198))*IT_0213;
    return create_ccomplex_return(IT_0214);
}

