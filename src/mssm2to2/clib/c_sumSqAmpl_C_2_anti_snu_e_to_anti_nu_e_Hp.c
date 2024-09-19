#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_e_to_anti_nu_e_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_e_to_anti_nu_e_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_snu_mu = param->m_snu_mu;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_B2*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = N_W2*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = IT_0006*IT_0012;
    const ccomplex_t IT_0014 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 2*IT_0015;
    const ccomplex_t IT_0017 = IT_0001*IT_0012;
    const ccomplex_t IT_0018 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0022 = IT_0013*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = IT_0020 + IT_0023;
    const ccomplex_t IT_0025 = IT_0016 + IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0011*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0030 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0031 = cpow(s_23 + (-0.5)*IT_0030 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0028*IT_0029*IT_0031;
    const ccomplex_t IT_0033 = N_B3*e_em;
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = N_W3*e_em;
    const ccomplex_t IT_0037 = IT_0006*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + -IT_0038);
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0042 = IT_0013*IT_0041;
    const ccomplex_t IT_0043 = 2*IT_0042;
    const ccomplex_t IT_0044 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0045 = IT_0017*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0048 = IT_0013*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = IT_0046 + IT_0049;
    const ccomplex_t IT_0051 = IT_0043 + IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0040*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0056 = cpow(s_23 + (-0.5)*IT_0030 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0054*IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0059 = cpow(s_23 + (-0.5)*IT_0030 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0060 = N_B4*e_em;
    const ccomplex_t IT_0061 = IT_0001*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = N_W4*e_em;
    const ccomplex_t IT_0064 = IT_0006*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + -IT_0065);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0069 = IT_0013*IT_0068;
    const ccomplex_t IT_0070 = 2*IT_0069;
    const ccomplex_t IT_0071 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0072 = IT_0017*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0075 = IT_0013*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = IT_0073 + IT_0076;
    const ccomplex_t IT_0078 = IT_0070 + IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = IT_0067*IT_0080;
    const ccomplex_t IT_0082 = IT_0058*IT_0059*IT_0081;
    const ccomplex_t IT_0083 = pow(m_e, 2);
    const ccomplex_t IT_0084 = pow(m_C_2, 2);
    const ccomplex_t IT_0085 = cpow((-2)*s_12 + -IT_0030 + IT_0083 + -IT_0084 
      + -reg_prop, -1);
    const ccomplex_t IT_0086 = pow(m_W, -1);
    const ccomplex_t IT_0087 = cpow(IT_0012, -1);
    const ccomplex_t IT_0088 = sin(beta);
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0006*IT_0086*IT_0087*IT_0088;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0006*IT_0086*IT_0087;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = IT_0085*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = m_C_2*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_e;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0006;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0090*IT_0099;
    const ccomplex_t IT_0101 = IT_0085*IT_0097*IT_0100;
    const ccomplex_t IT_0102 = N_B1*e_em;
    const ccomplex_t IT_0103 = IT_0001*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = N_W1*e_em;
    const ccomplex_t IT_0106 = IT_0006*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + -IT_0107);
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0111 = IT_0013*IT_0110;
    const ccomplex_t IT_0112 = 2*IT_0111;
    const ccomplex_t IT_0113 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0114 = IT_0017*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0117 = IT_0013*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = IT_0115 + IT_0118;
    const ccomplex_t IT_0120 = IT_0112 + IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0109*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0125 = cpow(s_23 + (-0.5)*IT_0030 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0126 = IT_0123*IT_0124*IT_0125;
    const ccomplex_t IT_0127 = (-0.5)*IT_0032 + (-0.5)*IT_0057 + (-0.5)
      *IT_0082 + IT_0096 + -IT_0101 + (-0.5)*IT_0126;
    const ccomplex_t IT_0128 = IT_0006*IT_0088;
    const ccomplex_t IT_0129 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = IT_0001*IT_0088;
    const ccomplex_t IT_0133 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0134 = IT_0132*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = IT_0131 + IT_0135;
    const ccomplex_t IT_0137 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0138 = IT_0128*IT_0137;
    const ccomplex_t IT_0139 = (-2)*IT_0138;
    const ccomplex_t IT_0140 = IT_0136 + IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = IT_0109*IT_0142;
    const ccomplex_t IT_0144 = IT_0125*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0147 = IT_0128*IT_0146;
    const ccomplex_t IT_0148 = (-2)*IT_0147;
    const ccomplex_t IT_0149 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0150 = IT_0132*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0153 = IT_0128*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = IT_0151 + IT_0154;
    const ccomplex_t IT_0156 = IT_0148 + IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = 0.5*IT_0157;
    const ccomplex_t IT_0159 = IT_0040*IT_0158;
    const ccomplex_t IT_0160 = IT_0056*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0163 = IT_0128*IT_0162;
    const ccomplex_t IT_0164 = (-2)*IT_0163;
    const ccomplex_t IT_0165 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0166 = IT_0132*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0169 = IT_0128*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = IT_0167 + IT_0170;
    const ccomplex_t IT_0172 = IT_0164 + IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0173;
    const ccomplex_t IT_0175 = IT_0067*IT_0174;
    const ccomplex_t IT_0176 = IT_0059*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0179 = IT_0128*IT_0178;
    const ccomplex_t IT_0180 = (-2)*IT_0179;
    const ccomplex_t IT_0181 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0182 = IT_0132*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0185 = IT_0128*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = IT_0183 + IT_0186;
    const ccomplex_t IT_0188 = IT_0180 + IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = IT_0011*IT_0190;
    const ccomplex_t IT_0192 = IT_0031*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = -IT_0095 + 0.5*IT_0145 + 0.5*IT_0161 + 0.5
      *IT_0177 + 0.5*IT_0193;
    const ccomplex_t IT_0195 = s_23*m_C_2;
    const ccomplex_t IT_0196 = s_12*s_23;
    const ccomplex_t IT_0197 = s_13*IT_0030;
    return create_ccomplex_return(2*IT_0127*(s_13*conj(IT_0127) + -conj
      (IT_0194)*IT_0195) + (-2)*IT_0194*(conj(IT_0127)*IT_0195 + (-2)*conj
      (IT_0194)*(IT_0196 + (-0.5)*IT_0197)));
}

