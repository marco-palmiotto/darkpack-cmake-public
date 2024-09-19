#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_snu_mu_to_anti_nu_mu_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_snu_mu_to_anti_nu_mu_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0003 = 0.5*IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = pow(m_mu, 2);
    const ccomplex_t IT_0011 = pow(m_C_1, 2);
    const ccomplex_t IT_0012 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_mu;
    const ccomplex_t IT_0015 = IT_0009*IT_0013*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = N_B1*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = N_W1*e_em;
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0020 + -IT_0023);
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0029 = IT_0001*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0027 + 0.5*IT_0030);
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = IT_0025*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0035 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0034*IT_0035;
    const ccomplex_t IT_0037 = N_B2*e_em;
    const ccomplex_t IT_0038 = IT_0017*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = N_W2*e_em;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0046 = IT_0001*IT_0045;
    const ccomplex_t IT_0047 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0048 = IT_0001*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + 0.5*IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0044*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0054 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0052*IT_0053*IT_0054;
    const ccomplex_t IT_0056 = N_B3*e_em;
    const ccomplex_t IT_0057 = IT_0017*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = N_W3*e_em;
    const ccomplex_t IT_0060 = IT_0001*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + -IT_0061);
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0065 = IT_0001*IT_0064;
    const ccomplex_t IT_0066 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0067 = IT_0001*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + 0.5*IT_0068);
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = IT_0063*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0073 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0071*IT_0072*IT_0073;
    const ccomplex_t IT_0075 = N_B4*e_em;
    const ccomplex_t IT_0076 = IT_0017*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_W4*e_em;
    const ccomplex_t IT_0079 = IT_0001*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + -IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0084 = IT_0001*IT_0083;
    const ccomplex_t IT_0085 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0086 = IT_0001*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0084 + 0.5*IT_0087);
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = IT_0082*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0092 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0090*IT_0091*IT_0092;
    const ccomplex_t IT_0094 = IT_0015 + (-0.5)*IT_0036 + (-0.5)*IT_0055 + (
      -0.5)*IT_0074 + (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0096 = IT_0001*IT_0095;
    const ccomplex_t IT_0097 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0098 = IT_0001*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + (-0.5)*IT_0099);
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0025*IT_0101;
    const ccomplex_t IT_0103 = IT_0035*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0106 = IT_0001*IT_0105;
    const ccomplex_t IT_0107 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0108 = IT_0001*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + (-0.5)*IT_0109);
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = IT_0044*IT_0111;
    const ccomplex_t IT_0113 = IT_0054*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0116 = IT_0001*IT_0115;
    const ccomplex_t IT_0117 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0118 = IT_0001*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0116 + (-0.5)*IT_0119);
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = IT_0063*IT_0121;
    const ccomplex_t IT_0123 = IT_0073*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0126 = IT_0001*IT_0125;
    const ccomplex_t IT_0127 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0128 = IT_0001*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + (-0.5)*IT_0129);
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = IT_0082*IT_0131;
    const ccomplex_t IT_0133 = IT_0092*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0104 + IT_0114 + IT_0124 + IT_0134;
    const ccomplex_t IT_0136 = m_C_1*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0001;
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = IT_0003*IT_0139;
    const ccomplex_t IT_0141 = IT_0013*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = pow(m_W, -2);
    const ccomplex_t IT_0144 = pow(m_W, 2);
    const ccomplex_t IT_0145 = s_34*IT_0144;
    const ccomplex_t IT_0146 = IT_0143*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = s_34 + IT_0147;
    const ccomplex_t IT_0149 = m_C_1*IT_0148;
    const ccomplex_t IT_0150 = 2*IT_0149;
    const ccomplex_t IT_0151 = (-2)*conj(IT_0142)*IT_0150;
    const ccomplex_t IT_0152 = 2*IT_0142;
    const ccomplex_t IT_0153 = pow(s_34, 2);
    const ccomplex_t IT_0154 = m_C_1*IT_0143*IT_0153;
    const ccomplex_t IT_0155 = (-2)*IT_0154;
    const ccomplex_t IT_0156 = conj(IT_0152)*IT_0155;
    const ccomplex_t IT_0157 = -IT_0104 + -IT_0114 + -IT_0124 + -IT_0134;
    const ccomplex_t IT_0158 = s_14*s_34*IT_0143;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = s_13 + IT_0159;
    const ccomplex_t IT_0161 = m_C_1*IT_0160;
    const ccomplex_t IT_0162 = (-2)*IT_0161;
    const ccomplex_t IT_0163 = conj(IT_0157)*IT_0162;
    const ccomplex_t IT_0164 = conj(IT_0094) + conj(IT_0137);
    const ccomplex_t IT_0165 = s_13*IT_0144;
    const ccomplex_t IT_0166 = s_14*s_34;
    const ccomplex_t IT_0167 = (-2)*IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = IT_0143*IT_0168;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = s_13 + IT_0170;
    const ccomplex_t IT_0172 = 4*IT_0171;
    const ccomplex_t IT_0173 = IT_0164*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0104 + 0.5*IT_0114 + 0.5*IT_0124 + 0.5
      *IT_0134;
    const ccomplex_t IT_0175 = (-0.25)*IT_0146;
    const ccomplex_t IT_0176 = s_34 + IT_0175;
    const ccomplex_t IT_0177 = m_C_1*IT_0176;
    const ccomplex_t IT_0178 = (-4)*IT_0177;
    const ccomplex_t IT_0179 = 2*IT_0178;
    const ccomplex_t IT_0180 = conj(IT_0174)*IT_0179;
    const ccomplex_t IT_0181 = IT_0151 + IT_0156 + IT_0163 + IT_0173 + IT_0180;
    const ccomplex_t IT_0182 = s_13*s_34;
    const ccomplex_t IT_0183 = s_13*s_14;
    const ccomplex_t IT_0184 = s_34*IT_0011;
    const ccomplex_t IT_0185 = -IT_0184;
    const ccomplex_t IT_0186 = IT_0183 + IT_0185;
    const ccomplex_t IT_0187 = pow(m_W, 4);
    const ccomplex_t IT_0188 = s_13*IT_0187;
    const ccomplex_t IT_0189 = IT_0143*IT_0188;
    const ccomplex_t IT_0190 = (-0.25)*IT_0189;
    const ccomplex_t IT_0191 = IT_0166 + IT_0190;
    const ccomplex_t IT_0192 = IT_0165 + IT_0167 + (-0.5)*IT_0189;
    const ccomplex_t IT_0193 = (-2)*conj(IT_0142);
    const ccomplex_t IT_0194 = 2*IT_0186;
    const ccomplex_t IT_0195 = 2*IT_0182;
    const ccomplex_t IT_0196 = pow(s_14, 2);
    const ccomplex_t IT_0197 = IT_0143*IT_0196;
    const ccomplex_t IT_0198 = -IT_0011;
    const ccomplex_t IT_0199 = IT_0197 + IT_0198;
    const ccomplex_t IT_0200 = s_13*IT_0199;
    const ccomplex_t IT_0201 = s_13*IT_0160;
    const ccomplex_t IT_0202 = 2*IT_0201;
    const ccomplex_t IT_0203 = 0.5*IT_0164;
    const ccomplex_t IT_0204 = 0.5*conj(IT_0174);
    const ccomplex_t IT_0205 = s_13*IT_0143*IT_0153;
    const ccomplex_t IT_0206 = (IT_0094 + IT_0137)*IT_0181 + (-2)*IT_0142*
      (IT_0150*IT_0164 + -conj(IT_0152)*IT_0182 + -conj(IT_0157)*IT_0186 + (-4)
      *conj(IT_0174)*IT_0191 + -IT_0192*IT_0193) + (-4)*IT_0174*((-0.25)*IT_0164
      *IT_0179 + conj(IT_0174)*IT_0192 + IT_0191*IT_0193 + (-0.25)*conj(IT_0157)
      *IT_0194 + (-0.25)*conj(IT_0152)*IT_0195) + 2*IT_0157*(conj(IT_0142)
      *IT_0186 + conj(IT_0157)*IT_0200 + 0.5*conj(IT_0152)*IT_0202 + IT_0162
      *IT_0203 + IT_0194*IT_0204) + 2*IT_0152*(conj(IT_0142)*IT_0182 + 0.5*conj
      (IT_0157)*IT_0202 + IT_0155*IT_0203 + IT_0195*IT_0204 + conj(IT_0152)
      *IT_0205);
    return create_ccomplex_return(IT_0206);
}

