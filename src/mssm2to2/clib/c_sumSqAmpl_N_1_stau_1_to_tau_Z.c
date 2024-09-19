#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_stau_1_to_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_stau_1_to_tau_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
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
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
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
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0006*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = N_d3*e_em*m_tau*IT_0020*U_stau_00;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0016 + 0.5*IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = pow(m_tau, 2);
    const ccomplex_t IT_0027 = pow(m_stau_1, 2);
    const ccomplex_t IT_0028 = cpow(s_23 + (-0.5)*IT_0026 + (-0.5)*IT_0027 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0013*IT_0030;
    const ccomplex_t IT_0032 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = IT_0007*IT_0032;
    const ccomplex_t IT_0035 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0036 = IT_0003*IT_0035;
    const ccomplex_t IT_0037 = IT_0007*IT_0035;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0033 + IT_0034 + 
      -IT_0036 + -IT_0037);
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = N_d1*e_em*m_tau*IT_0020*U_stau_00;
    const ccomplex_t IT_0044 = IT_0019*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 0.5*IT_0045);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = cpow(s_23 + (-0.5)*IT_0026 + (-0.5)*IT_0027 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0039*IT_0050;
    const ccomplex_t IT_0052 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0053 = IT_0003*IT_0052;
    const ccomplex_t IT_0054 = IT_0007*IT_0052;
    const ccomplex_t IT_0055 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0056 = IT_0003*IT_0055;
    const ccomplex_t IT_0057 = IT_0007*IT_0055;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0053 + IT_0054 + 
      -IT_0056 + -IT_0057);
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0061 = IT_0001*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = N_d4*e_em*m_tau*IT_0020*U_stau_00;
    const ccomplex_t IT_0064 = IT_0019*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + 0.5*IT_0065);
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = cpow(s_23 + (-0.5)*IT_0026 + (-0.5)*IT_0027 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0059*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0006)*U_stau_00*conj(U_stau_00) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_10));
    const ccomplex_t IT_0073 = (-2)*IT_0072;
    const ccomplex_t IT_0074 = pow(m_N_1, 2);
    const ccomplex_t IT_0075 = cpow(s_13 + (-0.5)*IT_0026 + (-0.5)*IT_0074 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0076 = IT_0047*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = IT_0073*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0006)*U_stau_00*conj(U_stau_01) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0080 = (-2)*IT_0079;
    const ccomplex_t IT_0081 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0082 = IT_0001*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = N_d1*e_em*m_tau*IT_0020*U_stau_01;
    const ccomplex_t IT_0085 = IT_0019*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + 0.5*IT_0086);
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = cpow(s_13 + (-0.5)*IT_0026 + (-0.5)*IT_0074 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0080*IT_0091;
    const ccomplex_t IT_0093 = 0.5*IT_0038;
    const ccomplex_t IT_0094 = IT_0050*IT_0093;
    const ccomplex_t IT_0095 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0096 = IT_0003*IT_0095;
    const ccomplex_t IT_0097 = IT_0007*IT_0095;
    const ccomplex_t IT_0098 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0099 = IT_0003*IT_0098;
    const ccomplex_t IT_0100 = IT_0007*IT_0098;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0096 + IT_0097 + 
      -IT_0099 + -IT_0100);
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0104 = IT_0001*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = N_d2*e_em*m_tau*IT_0020*U_stau_00;
    const ccomplex_t IT_0107 = IT_0019*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + 0.5*IT_0108);
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = cpow(s_23 + (-0.5)*IT_0026 + (-0.5)*IT_0027 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0102*IT_0113;
    const ccomplex_t IT_0115 = -IT_0031 + (-0.5)*IT_0051 + -IT_0071 + (-0.5)
      *IT_0078 + (-0.5)*IT_0092 + 0.5*IT_0094 + -IT_0114;
    const ccomplex_t IT_0116 = pow(m_Z, -2);
    const ccomplex_t IT_0117 = pow(s_14, 2);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = -IT_0074;
    const ccomplex_t IT_0120 = IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = s_13*IT_0120;
    const ccomplex_t IT_0122 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0123 = IT_0003*IT_0122;
    const ccomplex_t IT_0124 = IT_0007*IT_0122;
    const ccomplex_t IT_0125 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0126 = IT_0003*IT_0125;
    const ccomplex_t IT_0127 = IT_0007*IT_0125;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0123 + IT_0124 + 
      -IT_0126 + -IT_0127);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0131 = IT_0001*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0134 = IT_0006*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = conj(N_d3)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0137 = IT_0019*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0132 + IT_0135 + 
      -IT_0138);
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = IT_0028*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = IT_0129*IT_0142;
    const ccomplex_t IT_0144 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0145 = IT_0003*IT_0144;
    const ccomplex_t IT_0146 = IT_0007*IT_0144;
    const ccomplex_t IT_0147 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0148 = IT_0003*IT_0147;
    const ccomplex_t IT_0149 = IT_0007*IT_0147;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0145 + IT_0146 + 
      -IT_0148 + -IT_0149);
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0153 = IT_0001*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0156 = IT_0006*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = conj(N_d4)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0159 = IT_0019*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0154 + IT_0157 + 
      -IT_0160);
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0068*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0151*IT_0164;
    const ccomplex_t IT_0166 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0167 = IT_0001*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0170 = IT_0006*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = conj(N_d1)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0173 = IT_0019*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0168 + IT_0171 + 
      -IT_0174);
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = IT_0048*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = IT_0039*IT_0178;
    const ccomplex_t IT_0180 = IT_0093*IT_0178;
    const ccomplex_t IT_0181 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0182 = IT_0003*IT_0181;
    const ccomplex_t IT_0183 = IT_0007*IT_0181;
    const ccomplex_t IT_0184 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0185 = IT_0003*IT_0184;
    const ccomplex_t IT_0186 = IT_0007*IT_0184;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0182 + IT_0183 + 
      -IT_0185 + -IT_0186);
    const ccomplex_t IT_0188 = (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0190 = IT_0001*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0193 = IT_0006*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = conj(N_d2)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0196 = IT_0019*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0191 + IT_0194 + 
      -IT_0197);
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0111*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = IT_0188*IT_0201;
    const ccomplex_t IT_0203 = IT_0075*IT_0176;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = IT_0073*IT_0204;
    const ccomplex_t IT_0206 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0207 = IT_0001*IT_0206;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0210 = IT_0006*IT_0209;
    const ccomplex_t IT_0211 = 1.4142135623731*IT_0210;
    const ccomplex_t IT_0212 = conj(N_d1)*e_em*m_tau*IT_0020*U_stau_11;
    const ccomplex_t IT_0213 = IT_0019*IT_0212;
    const ccomplex_t IT_0214 = 1.4142135623731*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*(IT_0208 + IT_0211 + 
      -IT_0214);
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = IT_0089*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = IT_0080*IT_0218;
    const ccomplex_t IT_0220 = IT_0143 + IT_0165 + (-0.5)*IT_0179 + 0.5
      *IT_0180 + IT_0202 + 0.5*IT_0205 + 0.5*IT_0219;
    const ccomplex_t IT_0221 = m_N_1*m_tau;
    const ccomplex_t IT_0222 = IT_0120*IT_0221;
    const ccomplex_t IT_0223 = 0.5*IT_0031 + 0.25*IT_0051 + 0.5*IT_0071 + (
      -0.25)*IT_0094 + 0.5*IT_0114;
    const ccomplex_t IT_0224 = s_13*s_14;
    const ccomplex_t IT_0225 = s_34*IT_0074;
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = IT_0224 + IT_0226;
    const ccomplex_t IT_0228 = m_tau*IT_0120;
    const ccomplex_t IT_0229 = m_N_1*IT_0093;
    const ccomplex_t IT_0230 = IT_0178*IT_0229;
    const ccomplex_t IT_0231 = (-0.25)*IT_0230;
    const ccomplex_t IT_0232 = cpow(s_12 + (-0.5)*IT_0026 + 0.5*IT_0027 + 0.5
      *IT_0074 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0234 = e_em*IT_0003;
    const ccomplex_t IT_0235 = e_em*IT_0007;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*(IT_0234 + -IT_0235);
    const ccomplex_t IT_0237 = 0.5*IT_0236;
    const ccomplex_t IT_0238 = IT_0047*IT_0237;
    const ccomplex_t IT_0239 = IT_0232*IT_0233*IT_0238;
    const ccomplex_t IT_0240 = m_N_1*IT_0039;
    const ccomplex_t IT_0241 = IT_0178*IT_0240;
    const ccomplex_t IT_0242 = -IT_0229;
    const ccomplex_t IT_0243 = IT_0240 + IT_0242;
    const ccomplex_t IT_0244 = IT_0178*IT_0243;
    const ccomplex_t IT_0245 = m_N_2*IT_0102;
    const ccomplex_t IT_0246 = m_N_1*IT_0188;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = IT_0245 + IT_0247;
    const ccomplex_t IT_0249 = IT_0201*IT_0248;
    const ccomplex_t IT_0250 = m_N_3*IT_0013;
    const ccomplex_t IT_0251 = m_N_1*IT_0129;
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = IT_0250 + IT_0252;
    const ccomplex_t IT_0254 = IT_0142*IT_0253;
    const ccomplex_t IT_0255 = m_N_4*IT_0059;
    const ccomplex_t IT_0256 = m_N_1*IT_0151;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0255 + IT_0257;
    const ccomplex_t IT_0259 = IT_0164*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0261 = IT_0047*IT_0260;
    const ccomplex_t IT_0262 = IT_0232*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*IT_0262;
    const ccomplex_t IT_0264 = m_tau*IT_0263;
    const ccomplex_t IT_0265 = 0.5*IT_0239 + 0.25*IT_0241 + 0.25*IT_0244 + 0.5
      *IT_0249 + 0.5*IT_0254 + 0.5*IT_0259 + (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = conj(IT_0231) + conj(IT_0265);
    const ccomplex_t IT_0267 = pow(m_Z, 2);
    const ccomplex_t IT_0268 = s_13*IT_0267;
    const ccomplex_t IT_0269 = s_14*s_34;
    const ccomplex_t IT_0270 = (-2)*IT_0269;
    const ccomplex_t IT_0271 = IT_0268 + IT_0270;
    const ccomplex_t IT_0272 = IT_0116*IT_0271;
    const ccomplex_t IT_0273 = (-0.5)*IT_0272;
    const ccomplex_t IT_0274 = s_13 + IT_0273;
    const ccomplex_t IT_0275 = 4*IT_0231;
    const ccomplex_t IT_0276 = 0.5*IT_0078 + 0.5*IT_0092 + IT_0263;
    const ccomplex_t IT_0277 = s_14*s_34*IT_0116;
    const ccomplex_t IT_0278 = -IT_0277;
    const ccomplex_t IT_0279 = s_13 + IT_0278;
    const ccomplex_t IT_0280 = s_13*IT_0279;
    const ccomplex_t IT_0281 = IT_0176*IT_0237;
    const ccomplex_t IT_0282 = IT_0232*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = (-0.5)*IT_0205 + (-0.5)*IT_0219 + -IT_0283;
    const ccomplex_t IT_0285 = IT_0221*IT_0279;
    const ccomplex_t IT_0286 = 0.5*IT_0263;
    const ccomplex_t IT_0287 = m_N_1*IT_0279;
    const ccomplex_t IT_0288 = IT_0050*IT_0240;
    const ccomplex_t IT_0289 = 0.25*IT_0288;
    const ccomplex_t IT_0290 = IT_0176*IT_0260;
    const ccomplex_t IT_0291 = IT_0232*IT_0233*IT_0290;
    const ccomplex_t IT_0292 = m_N_4*IT_0151;
    const ccomplex_t IT_0293 = m_N_1*IT_0059;
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = IT_0292 + IT_0294;
    const ccomplex_t IT_0296 = IT_0070*IT_0295;
    const ccomplex_t IT_0297 = m_N_3*IT_0129;
    const ccomplex_t IT_0298 = m_N_1*IT_0013;
    const ccomplex_t IT_0299 = -IT_0298;
    const ccomplex_t IT_0300 = IT_0297 + IT_0299;
    const ccomplex_t IT_0301 = IT_0030*IT_0300;
    const ccomplex_t IT_0302 = m_tau*IT_0283;
    const ccomplex_t IT_0303 = IT_0050*IT_0229;
    const ccomplex_t IT_0304 = IT_0050*IT_0243;
    const ccomplex_t IT_0305 = m_N_2*IT_0188;
    const ccomplex_t IT_0306 = m_N_1*IT_0102;
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = IT_0305 + IT_0307;
    const ccomplex_t IT_0309 = IT_0113*IT_0308;
    const ccomplex_t IT_0310 = (-0.5)*IT_0291 + (-0.5)*IT_0296 + (-0.5)
      *IT_0301 + 0.5*IT_0302 + (-0.25)*IT_0303 + 0.25*IT_0304 + (-0.5)*IT_0309;
    const ccomplex_t IT_0311 = conj(IT_0289) + conj(IT_0310);
    const ccomplex_t IT_0312 = (-2)*IT_0115;
    const ccomplex_t IT_0313 = pow(s_34, 2);
    const ccomplex_t IT_0314 = IT_0116*IT_0313;
    const ccomplex_t IT_0315 = -IT_0026;
    const ccomplex_t IT_0316 = IT_0314 + IT_0315;
    const ccomplex_t IT_0317 = s_13*IT_0316;
    const ccomplex_t IT_0318 = IT_0221*IT_0316;
    const ccomplex_t IT_0319 = s_13*s_34;
    const ccomplex_t IT_0320 = s_14*IT_0026;
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = IT_0319 + IT_0321;
    const ccomplex_t IT_0323 = m_N_1*IT_0316;
    const ccomplex_t IT_0324 = (-2)*conj(IT_0115);
    const ccomplex_t IT_0325 = (-0.5)*IT_0283;
    const ccomplex_t IT_0326 = (-2)*IT_0220;
    const ccomplex_t IT_0327 = m_tau*IT_0279;
    const ccomplex_t IT_0328 = (-2)*IT_0276;
    const ccomplex_t IT_0329 = s_34*IT_0267;
    const ccomplex_t IT_0330 = IT_0116*IT_0329;
    const ccomplex_t IT_0331 = 0.5*IT_0330;
    const ccomplex_t IT_0332 = s_34 + IT_0331;
    const ccomplex_t IT_0333 = m_N_1*IT_0332;
    const ccomplex_t IT_0334 = 2*IT_0333;
    const ccomplex_t IT_0335 = (-0.5)*IT_0143 + (-0.5)*IT_0165 + 0.25*IT_0179 
      + (-0.25)*IT_0180 + (-0.5)*IT_0202;
    const ccomplex_t IT_0336 = (-0.25)*IT_0330;
    const ccomplex_t IT_0337 = s_34 + IT_0336;
    const ccomplex_t IT_0338 = m_N_1*IT_0337;
    const ccomplex_t IT_0339 = (-4)*IT_0338;
    const ccomplex_t IT_0340 = s_14*IT_0267;
    const ccomplex_t IT_0341 = IT_0116*IT_0340;
    const ccomplex_t IT_0342 = 0.5*IT_0341;
    const ccomplex_t IT_0343 = s_14 + IT_0342;
    const ccomplex_t IT_0344 = m_tau*IT_0343;
    const ccomplex_t IT_0345 = 2*IT_0344;
    const ccomplex_t IT_0346 = (-0.25)*IT_0341;
    const ccomplex_t IT_0347 = s_14 + IT_0346;
    const ccomplex_t IT_0348 = m_tau*IT_0347;
    const ccomplex_t IT_0349 = (-4)*IT_0348;
    const ccomplex_t IT_0350 = (-6)*conj(IT_0289);
    const ccomplex_t IT_0351 = (-0.166666666666667)*IT_0350;
    const ccomplex_t IT_0352 = conj(IT_0310) + IT_0351;
    const ccomplex_t IT_0353 = (-2)*conj(IT_0220);
    const ccomplex_t IT_0354 = (-2)*conj(IT_0276);
    const ccomplex_t IT_0355 = conj(IT_0115)*IT_0228 + conj(IT_0284)*IT_0323 +
       conj(IT_0325)*IT_0334 + conj(IT_0335)*IT_0339 + conj(IT_0223)*IT_0345 +
       conj(IT_0286)*IT_0349 + (-3)*IT_0221*IT_0352 + 0.5*IT_0287*IT_0353 + 0.5
      *IT_0327*IT_0354;
    const ccomplex_t IT_0356 = (-2)*IT_0284;
    const ccomplex_t IT_0357 = 0.5*IT_0285;
    const ccomplex_t IT_0358 = IT_0221*IT_0267;
    const ccomplex_t IT_0359 = 3*IT_0358;
    const ccomplex_t IT_0360 = (-3)*IT_0358;
    const ccomplex_t IT_0361 = pow(m_Z, 4);
    const ccomplex_t IT_0362 = s_13*IT_0361;
    const ccomplex_t IT_0363 = IT_0116*IT_0362;
    const ccomplex_t IT_0364 = (-0.5)*IT_0363;
    const ccomplex_t IT_0365 = IT_0268 + IT_0270 + IT_0364;
    const ccomplex_t IT_0366 = (-0.25)*IT_0363;
    const ccomplex_t IT_0367 = IT_0269 + IT_0366;
    const ccomplex_t IT_0368 = 0.5*IT_0227;
    const ccomplex_t IT_0369 = 0.5*IT_0322;
    const ccomplex_t IT_0370 = (-2)*conj(IT_0284);
    const ccomplex_t IT_0371 = (-0.333333333333333)*conj(IT_0276)*IT_0323;
    const ccomplex_t IT_0372 = (-0.333333333333333)*conj(IT_0220)*IT_0228;
    const ccomplex_t IT_0373 = IT_0221*(conj(IT_0231) + conj(IT_0265)) +
       IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = (-0.666666666666667)*IT_0274*IT_0311;
    const ccomplex_t IT_0375 = (-0.333333333333333)*conj(IT_0286);
    const ccomplex_t IT_0376 = (-0.333333333333333)*conj(IT_0223);
    const ccomplex_t IT_0377 = (-0.333333333333333)*conj(IT_0335);
    const ccomplex_t IT_0378 = IT_0345*IT_0377;
    const ccomplex_t IT_0379 = (-0.333333333333333)*conj(IT_0325);
    const ccomplex_t IT_0380 = IT_0349*IT_0379;
    const ccomplex_t IT_0381 = (-0.166666666666667)*IT_0287;
    const ccomplex_t IT_0382 = (-2)*IT_0289*(conj(IT_0284)*IT_0327 + 3*IT_0373
       + 3*IT_0374 + 3*IT_0334*IT_0375 + 3*IT_0339*IT_0376 + 3*IT_0378 + 3
      *IT_0380 + 3*IT_0324*IT_0381);
    const ccomplex_t IT_0383 = (-0.166666666666667)*IT_0327;
    const ccomplex_t IT_0384 = conj(IT_0115)*IT_0228 + 2*IT_0266*IT_0274 +
       conj(IT_0284)*IT_0323 + (-3)*IT_0221*IT_0352 + (-3)*IT_0349*IT_0375 + (-3
      )*IT_0345*IT_0376 + (-3)*IT_0339*IT_0377 + (-3)*IT_0334*IT_0379 + (-3)
      *IT_0353*IT_0381 + (-3)*IT_0354*IT_0383;
    const ccomplex_t IT_0385 = 2*IT_0115*(conj(IT_0115)*IT_0121 + conj(IT_0220
      )*IT_0222 + conj(IT_0223)*IT_0227 + IT_0228*IT_0266) + IT_0266*IT_0274
      *IT_0275 + (conj(IT_0276)*IT_0280 + conj(IT_0284)*IT_0285 + IT_0227*conj
      (IT_0286) + IT_0287*IT_0311)*IT_0312 + 2*IT_0276*(conj(IT_0276)*IT_0317 +
       conj(IT_0284)*IT_0318 + conj(IT_0286)*IT_0322 + IT_0311*IT_0323 + 0.5
      *IT_0280*IT_0324) + (IT_0280*conj(IT_0284) + IT_0266*IT_0287 + IT_0227
      *conj(IT_0325))*IT_0326 + (conj(IT_0220)*IT_0285 + conj(IT_0223)*IT_0322 +
       IT_0266*IT_0327)*IT_0328 + 2*IT_0231*IT_0355 + (conj(IT_0310)*IT_0327 +
       IT_0322*conj(IT_0335))*IT_0356 + 2*IT_0284*(conj(IT_0284)*IT_0317 + conj
      (IT_0276)*IT_0318 + IT_0266*IT_0323 + IT_0322*conj(IT_0325) + -conj
      (IT_0289)*IT_0327 + 0.5*IT_0280*IT_0353 + IT_0324*IT_0357) + 2*IT_0220*
      (IT_0121*conj(IT_0220) + conj(IT_0115)*IT_0222 + IT_0228*IT_0311 + IT_0227
      *conj(IT_0335) + IT_0354*IT_0357) + 2*IT_0286*(conj(IT_0276)*IT_0322 +
       IT_0311*IT_0334 + IT_0266*IT_0349 + conj(IT_0335)*IT_0359 + conj(IT_0325)
      *IT_0360 + (-2)*conj(IT_0286)*IT_0365 + (-4)*conj(IT_0223)*IT_0367 +
       IT_0324*IT_0368) + 2*IT_0325*(conj(IT_0284)*IT_0322 + IT_0266*IT_0334 +
       IT_0311*IT_0349 + conj(IT_0223)*IT_0359 + conj(IT_0286)*IT_0360 + (-2)
      *conj(IT_0325)*IT_0365 + (-4)*conj(IT_0335)*IT_0367 + IT_0353*IT_0368) + 2
      *IT_0223*(conj(IT_0115)*IT_0227 + IT_0311*IT_0339 + IT_0266*IT_0345 + conj
      (IT_0325)*IT_0359 + conj(IT_0335)*IT_0360 + (-2)*conj(IT_0223)*IT_0365 + (
      -4)*conj(IT_0286)*IT_0367 + IT_0354*IT_0369) + 2*IT_0335*(conj(IT_0220)
      *IT_0227 + IT_0266*IT_0339 + IT_0311*IT_0345 + conj(IT_0286)*IT_0359 +
       conj(IT_0223)*IT_0360 + (-2)*conj(IT_0335)*IT_0365 + (-4)*conj(IT_0325)
      *IT_0367 + IT_0369*IT_0370) + IT_0382 + (-6)*IT_0310*(IT_0373 + IT_0374 +
       IT_0334*IT_0375 + IT_0339*IT_0376 + IT_0378 + IT_0380 + IT_0324*IT_0381 +
       IT_0370*IT_0383) + 2*IT_0265*IT_0384;
    return create_ccomplex_return(IT_0385);
}

