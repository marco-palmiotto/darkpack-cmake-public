#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_stau_1_to_anti_tau_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_stau_1_to_anti_tau_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
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
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = m_C_2*m_tau;
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_C_2, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = pow(m_tau, 2);
    const ccomplex_t IT_0009 = pow(m_stau_1, 2);
    const ccomplex_t IT_0010 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0014 + 0.5*IT_0017);
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = conj(N_B1)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = IT_0012*IT_0026;
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = conj(N_d1)*e_em*m_tau*IT_0028*conj(U_stau_00);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0024 + 0.5*IT_0031);
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = IT_0019*IT_0033;
    const ccomplex_t IT_0035 = IT_0010*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0038 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0039 = IT_0012*IT_0038;
    const ccomplex_t IT_0040 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0041 = IT_0012*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + 0.5*IT_0042);
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = conj(N_B2)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0046 = IT_0021*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = conj(N_d2)*e_em*m_tau*IT_0028*conj(U_stau_00);
    const ccomplex_t IT_0049 = IT_0027*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + 0.5*IT_0050);
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = IT_0044*IT_0052;
    const ccomplex_t IT_0054 = IT_0037*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0057 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0058 = IT_0012*IT_0057;
    const ccomplex_t IT_0059 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0060 = IT_0012*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + 0.5*IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = conj(N_B3)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0065 = IT_0021*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = conj(N_d3)*e_em*m_tau*IT_0028*conj(U_stau_00);
    const ccomplex_t IT_0068 = IT_0027*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + 0.5*IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0063*IT_0071;
    const ccomplex_t IT_0073 = IT_0056*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0076 = IT_0012*IT_0075;
    const ccomplex_t IT_0077 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0078 = IT_0012*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0076 + 0.5*IT_0079);
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = conj(N_B4)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0083 = IT_0021*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = conj(N_d4)*e_em*m_tau*IT_0028*conj(U_stau_00);
    const ccomplex_t IT_0086 = IT_0027*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0084 + 0.5*IT_0087);
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = IT_0081*IT_0089;
    const ccomplex_t IT_0091 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0012
      *conj(U_stau_00);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0012*IT_0026*IT_0028;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = cpow(s_13 + (-0.5)*IT_0004 + (-0.5)*IT_0008 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0095*IT_0100;
    const ccomplex_t IT_0102 = -IT_0036 + -IT_0055 + -IT_0074 + -IT_0093 + (
      -0.5)*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0012;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = IT_0098*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0095*IT_0106;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = 2*IT_0108;
    const ccomplex_t IT_0110 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0111 = IT_0012*IT_0110;
    const ccomplex_t IT_0112 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0113 = IT_0012*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + (-0.5)*IT_0114);
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = N_B1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0118 = IT_0021*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = N_W1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0121 = IT_0012*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = N_d1*e_em*m_tau*IT_0028*conj(U_stau_10);
    const ccomplex_t IT_0124 = IT_0027*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0119 + IT_0122 + 
      -IT_0125);
    const ccomplex_t IT_0127 = (-0.5)*IT_0126;
    const ccomplex_t IT_0128 = IT_0116*IT_0127;
    const ccomplex_t IT_0129 = IT_0010*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0132 = IT_0012*IT_0131;
    const ccomplex_t IT_0133 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0134 = IT_0012*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + (-0.5)*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = N_B2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0139 = IT_0021*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = N_W2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0142 = IT_0012*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = N_d2*e_em*m_tau*IT_0028*conj(U_stau_10);
    const ccomplex_t IT_0145 = IT_0027*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0140 + IT_0143 + 
      -IT_0146);
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0137*IT_0148;
    const ccomplex_t IT_0150 = IT_0037*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0153 = IT_0012*IT_0152;
    const ccomplex_t IT_0154 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0155 = IT_0012*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + (-0.5)*IT_0156);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = N_B3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0160 = IT_0021*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = N_W3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0163 = IT_0012*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = N_d3*e_em*m_tau*IT_0028*conj(U_stau_10);
    const ccomplex_t IT_0166 = IT_0027*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0161 + IT_0164 + 
      -IT_0167);
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0158*IT_0169;
    const ccomplex_t IT_0171 = IT_0056*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0174 = IT_0012*IT_0173;
    const ccomplex_t IT_0175 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0176 = IT_0012*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + (-0.5)*IT_0177);
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = N_B4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0181 = IT_0021*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = N_W4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0184 = IT_0012*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = N_d4*e_em*m_tau*IT_0028*conj(U_stau_10);
    const ccomplex_t IT_0187 = IT_0027*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0182 + IT_0185 + 
      -IT_0188);
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = IT_0179*IT_0190;
    const ccomplex_t IT_0192 = IT_0091*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = -IT_0130 + -IT_0151 + -IT_0172 + -IT_0193;
    const ccomplex_t IT_0195 = (-2)*IT_0194;
    const ccomplex_t IT_0196 = IT_0109 + IT_0195;
    const ccomplex_t IT_0197 = 2*conj(IT_0108);
    const ccomplex_t IT_0198 = (-2)*conj(IT_0194);
    const ccomplex_t IT_0199 = IT_0197 + IT_0198;
    const ccomplex_t IT_0200 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = s_13 + IT_0201;
    const ccomplex_t IT_0203 = m_C_2*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0205 = IT_0019*IT_0127;
    const ccomplex_t IT_0206 = IT_0010*IT_0204*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0208 = IT_0044*IT_0148;
    const ccomplex_t IT_0209 = IT_0037*IT_0207*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0211 = IT_0063*IT_0169;
    const ccomplex_t IT_0212 = IT_0056*IT_0210*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0214 = IT_0081*IT_0190;
    const ccomplex_t IT_0215 = IT_0091*IT_0213*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0206 + (-0.5)*IT_0209 + (-0.5)
      *IT_0212 + (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = IT_0130 + IT_0151 + IT_0172 + IT_0193;
    const ccomplex_t IT_0218 = m_C_2*IT_0217;
    const ccomplex_t IT_0219 = 0.5*IT_0218;
    const ccomplex_t IT_0220 = conj(IT_0216) + conj(IT_0219);
    const ccomplex_t IT_0221 = IT_0216 + IT_0219;
    const ccomplex_t IT_0222 = IT_0036 + IT_0055 + IT_0074 + IT_0093;
    const ccomplex_t IT_0223 = m_C_2*IT_0222;
    const ccomplex_t IT_0224 = (-0.5)*IT_0223;
    const ccomplex_t IT_0225 = IT_0033*IT_0116;
    const ccomplex_t IT_0226 = IT_0010*IT_0204*IT_0225;
    const ccomplex_t IT_0227 = IT_0052*IT_0137;
    const ccomplex_t IT_0228 = IT_0037*IT_0207*IT_0227;
    const ccomplex_t IT_0229 = IT_0071*IT_0158;
    const ccomplex_t IT_0230 = IT_0056*IT_0210*IT_0229;
    const ccomplex_t IT_0231 = IT_0089*IT_0179;
    const ccomplex_t IT_0232 = IT_0091*IT_0213*IT_0231;
    const ccomplex_t IT_0233 = 0.5*IT_0226 + 0.5*IT_0228 + 0.5*IT_0230 + 0.5
      *IT_0232;
    const ccomplex_t IT_0234 = conj(IT_0224) + conj(IT_0233);
    const ccomplex_t IT_0235 = 2*IT_0102;
    const ccomplex_t IT_0236 = IT_0224 + IT_0233;
    const ccomplex_t IT_0237 = 2*conj(IT_0102);
    const ccomplex_t IT_0238 = (-0.5)*IT_0107;
    const ccomplex_t IT_0239 = pow(s_34, 2);
    const ccomplex_t IT_0240 = IT_0001*IT_0239;
    const ccomplex_t IT_0241 = -IT_0008;
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = s_13*IT_0242;
    const ccomplex_t IT_0244 = 2*IT_0238;
    const ccomplex_t IT_0245 = 0.5*IT_0101;
    const ccomplex_t IT_0246 = IT_0000*IT_0242;
    const ccomplex_t IT_0247 = 2*conj(IT_0238);
    const ccomplex_t IT_0248 = 2*IT_0245;
    const ccomplex_t IT_0249 = m_tau*IT_0202;
    const ccomplex_t IT_0250 = 2*IT_0216;
    const ccomplex_t IT_0251 = 0.5*IT_0250;
    const ccomplex_t IT_0252 = IT_0219 + IT_0251;
    const ccomplex_t IT_0253 = 2*conj(IT_0216);
    const ccomplex_t IT_0254 = pow(m_W, 2);
    const ccomplex_t IT_0255 = s_13*IT_0254;
    const ccomplex_t IT_0256 = s_14*s_34;
    const ccomplex_t IT_0257 = (-2)*IT_0256;
    const ccomplex_t IT_0258 = IT_0255 + IT_0257;
    const ccomplex_t IT_0259 = IT_0001*IT_0258;
    const ccomplex_t IT_0260 = (-0.5)*IT_0259;
    const ccomplex_t IT_0261 = s_13 + IT_0260;
    const ccomplex_t IT_0262 = 4*IT_0224;
    const ccomplex_t IT_0263 = (-6)*conj(IT_0216);
    const ccomplex_t IT_0264 = s_13*IT_0006;
    const ccomplex_t IT_0265 = 2*IT_0194;
    const ccomplex_t IT_0266 = (-2)*IT_0108;
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = (-2)*conj(IT_0108);
    const ccomplex_t IT_0269 = s_13*IT_0202;
    const ccomplex_t IT_0270 = 2*conj(IT_0194);
    const ccomplex_t IT_0271 = IT_0268 + IT_0270;
    const ccomplex_t IT_0272 = (-2)*IT_0102;
    const ccomplex_t IT_0273 = (-2)*conj(IT_0102);
    const ccomplex_t IT_0274 = IT_0000*IT_0202;
    const ccomplex_t IT_0275 = (-2)*IT_0238;
    const ccomplex_t IT_0276 = (-2)*conj(IT_0238);
    const ccomplex_t IT_0277 = m_C_2*IT_0242;
    const ccomplex_t IT_0278 = (-2)*IT_0245;
    const ccomplex_t IT_0279 = (-2)*conj(IT_0245);
    const ccomplex_t IT_0280 = m_tau*IT_0006;
    const ccomplex_t IT_0281 = (-2)*IT_0216;
    const ccomplex_t IT_0282 = (-2)*conj(IT_0216);
    const ccomplex_t IT_0283 = 0.5*IT_0130 + 0.5*IT_0151 + 0.5*IT_0172 + 0.5
      *IT_0193;
    const ccomplex_t IT_0284 = s_14*IT_0254;
    const ccomplex_t IT_0285 = IT_0001*IT_0284;
    const ccomplex_t IT_0286 = 0.5*IT_0285;
    const ccomplex_t IT_0287 = s_14 + IT_0286;
    const ccomplex_t IT_0288 = m_tau*IT_0287;
    const ccomplex_t IT_0289 = 2*IT_0288;
    const ccomplex_t IT_0290 = s_34*IT_0254;
    const ccomplex_t IT_0291 = IT_0001*IT_0290;
    const ccomplex_t IT_0292 = (-0.25)*IT_0291;
    const ccomplex_t IT_0293 = s_34 + IT_0292;
    const ccomplex_t IT_0294 = m_C_2*IT_0293;
    const ccomplex_t IT_0295 = (-4)*IT_0294;
    const ccomplex_t IT_0296 = s_13*s_14;
    const ccomplex_t IT_0297 = s_34*IT_0004;
    const ccomplex_t IT_0298 = -IT_0297;
    const ccomplex_t IT_0299 = IT_0296 + IT_0298;
    const ccomplex_t IT_0300 = s_13*s_34;
    const ccomplex_t IT_0301 = s_14*IT_0008;
    const ccomplex_t IT_0302 = -IT_0301;
    const ccomplex_t IT_0303 = IT_0300 + IT_0302;
    const ccomplex_t IT_0304 = 0.5*IT_0303;
    const ccomplex_t IT_0305 = 0.5*IT_0036 + 0.5*IT_0055 + 0.5*IT_0074 + 0.5
      *IT_0093;
    const ccomplex_t IT_0306 = IT_0000*IT_0254;
    const ccomplex_t IT_0307 = (-3)*IT_0306;
    const ccomplex_t IT_0308 = pow(m_W, 4);
    const ccomplex_t IT_0309 = s_13*IT_0308;
    const ccomplex_t IT_0310 = IT_0001*IT_0309;
    const ccomplex_t IT_0311 = (-0.5)*IT_0310;
    const ccomplex_t IT_0312 = IT_0255 + IT_0257 + IT_0311;
    const ccomplex_t IT_0313 = 2*IT_0283*(IT_0234*IT_0289 + (conj(IT_0219) +
       0.5*IT_0253)*IT_0295 + 0.5*IT_0271*IT_0299 + IT_0247*IT_0304 + -conj
      (IT_0305)*IT_0307 + (-2)*conj(IT_0283)*IT_0312);
    const ccomplex_t IT_0314 = (-0.5)*IT_0303;
    const ccomplex_t IT_0315 = IT_0007*(conj(IT_0102)*IT_0196 + IT_0102
      *IT_0199) + IT_0203*(IT_0196*IT_0220 + IT_0199*IT_0221 + IT_0234*IT_0235 +
       IT_0236*IT_0237) + conj(IT_0238)*IT_0243*IT_0244 + IT_0245*IT_0246
      *IT_0247 + conj(IT_0245)*(IT_0244*IT_0246 + IT_0243*IT_0248) + 2*IT_0249*
      (0.5*IT_0234*IT_0244 + 0.5*IT_0236*IT_0247 + 0.5*conj(IT_0219)*IT_0248 +
       conj(IT_0245)*IT_0252 + 0.5*IT_0245*IT_0253) + 4*IT_0261*(IT_0220*IT_0221
       + IT_0234*(IT_0233 + 0.25*IT_0262)) + (-6)*IT_0000*(IT_0221*IT_0234 +
       IT_0236*(conj(IT_0219) + (-0.166666666666667)*IT_0263)) + IT_0264*(conj
      (IT_0108)*IT_0109 + conj(IT_0102)*IT_0235 + conj(IT_0194)*IT_0267 +
       IT_0194*IT_0268) + IT_0269*(conj(IT_0238)*IT_0267 + IT_0238*IT_0271 +
       conj(IT_0245)*IT_0272 + IT_0245*IT_0273) + IT_0274*(conj(IT_0245)*IT_0267
       + IT_0245*IT_0271 + conj(IT_0102)*IT_0275 + IT_0102*IT_0276) + IT_0277*
      (IT_0220*IT_0275 + IT_0221*IT_0276 + IT_0234*IT_0278 + IT_0236*IT_0279) +
       IT_0280*(IT_0234*IT_0267 + IT_0236*IT_0271 + conj(IT_0219)*IT_0272 +
       IT_0219*IT_0273 + conj(IT_0102)*IT_0281 + IT_0102*IT_0282) + 2*conj
      (IT_0283)*(IT_0236*IT_0289 + IT_0252*IT_0295 + 0.5*IT_0267*IT_0299 +
       IT_0244*IT_0304) + IT_0313 + (-2)*conj(IT_0305)*((IT_0219 + (-0.5)
      *IT_0281)*IT_0289 + IT_0236*IT_0295 + (-0.5)*IT_0235*IT_0299 + IT_0278
      *IT_0314) + (-2)*IT_0305*((conj(IT_0219) + (-0.5)*IT_0282)*IT_0289 +
       IT_0234*IT_0295 + (-0.5)*IT_0237*IT_0299 + conj(IT_0283)*IT_0307 + 2*conj
      (IT_0305)*IT_0312 + IT_0279*IT_0314);
    return create_ccomplex_return(IT_0315);
}

