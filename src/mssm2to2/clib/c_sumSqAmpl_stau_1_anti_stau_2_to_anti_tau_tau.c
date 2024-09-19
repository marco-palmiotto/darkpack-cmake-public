#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_anti_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_1 = param->m_stau_1;
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
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = mu_h*m_tau*IT_0000;
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = -IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0001*IT_0003;
    const ccomplex_t IT_0013 = IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -2);
    const ccomplex_t IT_0016 = pow(m_tau, 2);
    const ccomplex_t IT_0017 = (-2)*U_stau_00*(((-0.5)*m_W*IT_0013*(IT_0005 + 
      -IT_0004*IT_0015) + IT_0000*IT_0002*IT_0003*IT_0005*IT_0016)*conj
      (U_stau_01) + -1./2*IT_0002*IT_0005*IT_0008*IT_0009*conj(U_stau_11)) +
       U_stau_10*(IT_0002*IT_0005*IT_0008*IT_0009*conj(U_stau_01) + 2*(m_W
      *IT_0004*IT_0013*IT_0015 + -IT_0000*IT_0002*IT_0003*IT_0005*IT_0016)*conj
      (U_stau_11));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_tau*IT_0003*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = pow(m_stau_1, 2);
    const ccomplex_t IT_0022 = pow(m_stau_2, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + -IT_0021 + -IT_0022 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
      *IT_0005*IT_0009;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0009;
    const ccomplex_t IT_0029 = IT_0003*IT_0010;
    const ccomplex_t IT_0030 = IT_0028 + IT_0029;
    const ccomplex_t IT_0031 = 2*U_stau_00*((IT_0000*IT_0002*IT_0005*IT_0009
      *IT_0016 + (-0.5)*m_W*(IT_0005 + -IT_0004*IT_0015)*IT_0030)*conj(U_stau_01
      ) + 1./2*IT_0002*IT_0003*IT_0005*IT_0008*conj(U_stau_11)) + U_stau_10*
      (IT_0002*IT_0003*IT_0005*IT_0008*conj(U_stau_01) + 2*(IT_0000*IT_0002
      *IT_0005*IT_0009*IT_0016 + -m_W*IT_0004*IT_0015*IT_0030)*conj(U_stau_11));
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(e_em*IT_0031 +
       1.4142135623731*A_tau*IT_0009*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0021 + -IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = -IT_0025 + -IT_0037;
    const ccomplex_t IT_0039 = cpow(IT_0014, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0005
      *IT_0014 + -IT_0004*IT_0039)*U_stau_00*conj(U_stau_01) + IT_0004*IT_0039
      *U_stau_10*conj(U_stau_11));
    const ccomplex_t IT_0041 = (-2)*IT_0040;
    const ccomplex_t IT_0042 = IT_0004*IT_0039;
    const ccomplex_t IT_0043 = e_em*IT_0042;
    const ccomplex_t IT_0044 = IT_0005*IT_0014;
    const ccomplex_t IT_0045 = e_em*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0043 + -IT_0045);
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0021 + -IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0041*IT_0050;
    const ccomplex_t IT_0052 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0022 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0053 = N_B1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0054 = IT_0039*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_W1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0057 = IT_0005*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = IT_0002*IT_0005;
    const ccomplex_t IT_0060 = N_d1*e_em*m_tau*IT_0000*conj(U_stau_11);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0055 + IT_0058 + 
      -IT_0062);
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0066 = IT_0039*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0069 = IT_0005*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = conj(N_d1)*e_em*m_tau*IT_0000*U_stau_10;
    const ccomplex_t IT_0072 = IT_0059*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0067 + IT_0070 + 
      -IT_0073);
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = IT_0064*IT_0075;
    const ccomplex_t IT_0077 = IT_0052*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0080 = IT_0039*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0083 = IT_0005*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = conj(N_d2)*e_em*m_tau*IT_0000*U_stau_10;
    const ccomplex_t IT_0086 = IT_0059*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0081 + IT_0084 + 
      -IT_0087);
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = N_B2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0091 = IT_0039*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = N_W2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0094 = IT_0005*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = N_d2*e_em*m_tau*IT_0000*conj(U_stau_11);
    const ccomplex_t IT_0097 = IT_0059*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0092 + IT_0095 + 
      -IT_0098);
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0089*IT_0100;
    const ccomplex_t IT_0102 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0022 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0106 = IT_0039*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0109 = IT_0005*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = conj(N_d3)*e_em*m_tau*IT_0000*U_stau_10;
    const ccomplex_t IT_0112 = IT_0059*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0107 + IT_0110 + 
      -IT_0113);
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = N_B3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0117 = IT_0039*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_W3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0120 = IT_0005*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = N_d3*e_em*m_tau*IT_0000*conj(U_stau_11);
    const ccomplex_t IT_0123 = IT_0059*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0118 + IT_0121 + 
      -IT_0124);
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = IT_0115*IT_0126;
    const ccomplex_t IT_0128 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0022 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0129 = IT_0127*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0132 = IT_0039*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0135 = IT_0005*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = conj(N_d4)*e_em*m_tau*IT_0000*U_stau_10;
    const ccomplex_t IT_0138 = IT_0059*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0133 + IT_0136 + 
      -IT_0139);
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = N_B4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0143 = IT_0039*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = N_W4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0146 = IT_0005*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = N_d4*e_em*m_tau*IT_0000*conj(U_stau_11);
    const ccomplex_t IT_0149 = IT_0059*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0144 + IT_0147 + 
      -IT_0150);
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0141*IT_0152;
    const ccomplex_t IT_0154 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0022 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = -IT_0051 + 0.5*IT_0078 + 0.5*IT_0104 + 0.5
      *IT_0130 + 0.5*IT_0156;
    const ccomplex_t IT_0158 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0159 = IT_0039*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = N_d4*e_em*m_tau*IT_0000*U_stau_00;
    const ccomplex_t IT_0162 = IT_0059*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + 0.5*IT_0163);
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = conj(N_B4)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0167 = IT_0039*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = conj(N_d4)*e_em*m_tau*IT_0000*conj(U_stau_01);
    const ccomplex_t IT_0170 = IT_0059*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*(IT_0168 + 0.5*IT_0171);
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = IT_0165*IT_0173;
    const ccomplex_t IT_0175 = IT_0154*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0039;
    const ccomplex_t IT_0178 = IT_0048*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = IT_0041*IT_0179;
    const ccomplex_t IT_0181 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0182 = IT_0039*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = N_d2*e_em*m_tau*IT_0000*U_stau_00;
    const ccomplex_t IT_0185 = IT_0059*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0183 + 0.5*IT_0186);
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = conj(N_B2)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0190 = IT_0039*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = conj(N_d2)*e_em*m_tau*IT_0000*conj(U_stau_01);
    const ccomplex_t IT_0193 = IT_0059*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*(IT_0191 + 0.5*IT_0194);
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0188*IT_0196;
    const ccomplex_t IT_0198 = IT_0102*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0201 = IT_0039*IT_0200;
    const ccomplex_t IT_0202 = 1.4142135623731*IT_0201;
    const ccomplex_t IT_0203 = N_d3*e_em*m_tau*IT_0000*U_stau_00;
    const ccomplex_t IT_0204 = IT_0059*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*(IT_0202 + 0.5*IT_0205);
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = conj(N_B3)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0209 = IT_0039*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = conj(N_d3)*e_em*m_tau*IT_0000*conj(U_stau_01);
    const ccomplex_t IT_0212 = IT_0059*IT_0211;
    const ccomplex_t IT_0213 = 1.4142135623731*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*(IT_0210 + 0.5*IT_0213);
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = IT_0207*IT_0215;
    const ccomplex_t IT_0217 = IT_0128*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = conj(N_B1)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0220 = IT_0039*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = conj(N_d1)*e_em*m_tau*IT_0000*conj(U_stau_01);
    const ccomplex_t IT_0223 = IT_0059*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0221 + 0.5*IT_0224);
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0228 = IT_0039*IT_0227;
    const ccomplex_t IT_0229 = 1.4142135623731*IT_0228;
    const ccomplex_t IT_0230 = N_d1*e_em*m_tau*IT_0000*U_stau_00;
    const ccomplex_t IT_0231 = IT_0059*IT_0230;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*(IT_0229 + 0.5*IT_0232);
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0226*IT_0234;
    const ccomplex_t IT_0236 = IT_0052*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0176 + -IT_0180 + 0.5*IT_0199 + 0.5
      *IT_0218 + 0.5*IT_0237;
    const ccomplex_t IT_0239 = -s_14;
    const ccomplex_t IT_0240 = s_13 + IT_0239;
    const ccomplex_t IT_0241 = m_tau*IT_0240;
    const ccomplex_t IT_0242 = IT_0064*IT_0234;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0244 = IT_0052*IT_0242*IT_0243;
    const ccomplex_t IT_0245 = IT_0100*IT_0188;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0247 = IT_0102*IT_0245*IT_0246;
    const ccomplex_t IT_0248 = IT_0126*IT_0207;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0250 = IT_0128*IT_0248*IT_0249;
    const ccomplex_t IT_0251 = IT_0152*IT_0165;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0253 = IT_0154*IT_0251*IT_0252;
    const ccomplex_t IT_0254 = -IT_0040;
    const ccomplex_t IT_0255 = m_tau*IT_0254;
    const ccomplex_t IT_0256 = m_tau*IT_0041;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0255 + IT_0257;
    const ccomplex_t IT_0259 = IT_0050*IT_0258;
    const ccomplex_t IT_0260 = IT_0179*IT_0258;
    const ccomplex_t IT_0261 = m_tau*IT_0078;
    const ccomplex_t IT_0262 = m_tau*IT_0130;
    const ccomplex_t IT_0263 = m_tau*IT_0156;
    const ccomplex_t IT_0264 = m_tau*IT_0104;
    const ccomplex_t IT_0265 = 0.5*IT_0244 + 0.5*IT_0247 + 0.5*IT_0250 + 0.5
      *IT_0253 + -IT_0259 + IT_0260 + (-0.5)*IT_0261 + (-0.5)*IT_0262 + (-0.5)
      *IT_0263 + (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = IT_0089*IT_0196;
    const ccomplex_t IT_0267 = IT_0102*IT_0246*IT_0266;
    const ccomplex_t IT_0268 = IT_0115*IT_0215;
    const ccomplex_t IT_0269 = IT_0128*IT_0249*IT_0268;
    const ccomplex_t IT_0270 = IT_0141*IT_0173;
    const ccomplex_t IT_0271 = IT_0154*IT_0252*IT_0270;
    const ccomplex_t IT_0272 = IT_0075*IT_0226;
    const ccomplex_t IT_0273 = IT_0052*IT_0243*IT_0272;
    const ccomplex_t IT_0274 = m_tau*IT_0237;
    const ccomplex_t IT_0275 = m_tau*IT_0199;
    const ccomplex_t IT_0276 = m_tau*IT_0218;
    const ccomplex_t IT_0277 = m_tau*IT_0176;
    const ccomplex_t IT_0278 = IT_0259 + -IT_0260 + 0.5*IT_0267 + 0.5*IT_0269 
      + 0.5*IT_0271 + 0.5*IT_0273 + (-0.5)*IT_0274 + (-0.5)*IT_0275 + (-0.5)
      *IT_0276 + (-0.5)*IT_0277;
    const ccomplex_t IT_0279 = -IT_0016;
    const ccomplex_t IT_0280 = s_34 + IT_0279;
    const ccomplex_t IT_0281 = s_13*m_tau;
    const ccomplex_t IT_0282 = e_em*m_tau*IT_0000*IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0283 = (-0.5)*IT_0282;
    const ccomplex_t IT_0284 = e_em*mu_h*m_tau*IT_0000*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0285 = IT_0005*IT_0284;
    const ccomplex_t IT_0286 = e_em*mu_h*m_tau*IT_0000*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0287 = IT_0005*IT_0286;
    const ccomplex_t IT_0288 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0289 = IT_0010*IT_0288;
    const ccomplex_t IT_0290 = 1.4142135623731*IT_0289;
    const ccomplex_t IT_0291 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0292 = IT_0010*IT_0291;
    const ccomplex_t IT_0293 = 1.4142135623731*IT_0292;
    const ccomplex_t IT_0294 = (-0.5)*IT_0285 + 0.5*IT_0287 + (-0.5)*IT_0290 +
       0.5*IT_0293;
    const ccomplex_t IT_0295 = IT_0283*IT_0294;
    const ccomplex_t IT_0296 = cpow((-2)*s_12 + -IT_0021 + -IT_0022 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0297 = IT_0295*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = e_em*m_tau*IT_0000*IT_0005;
    const ccomplex_t IT_0300 = 0.5*IT_0299;
    const ccomplex_t IT_0301 = IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0302 = IT_0284*IT_0301;
    const ccomplex_t IT_0303 = IT_0286*IT_0301;
    const ccomplex_t IT_0304 = IT_0001*IT_0288;
    const ccomplex_t IT_0305 = 1.4142135623731*IT_0304;
    const ccomplex_t IT_0306 = IT_0001*IT_0291;
    const ccomplex_t IT_0307 = 1.4142135623731*IT_0306;
    const ccomplex_t IT_0308 = (-0.5)*IT_0302 + 0.5*IT_0303 + 0.5*IT_0305 + (
      -0.5)*IT_0307;
    const ccomplex_t IT_0309 = IT_0300*IT_0308;
    const ccomplex_t IT_0310 = pow(m_Z, 2);
    const ccomplex_t IT_0311 = cpow((-2)*s_12 + -IT_0021 + -IT_0022 + IT_0310 
      + -reg_prop, -1);
    const ccomplex_t IT_0312 = IT_0309*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*IT_0312;
    const ccomplex_t IT_0314 = -IT_0298 + -IT_0313;
    const ccomplex_t IT_0315 = s_13 + s_14;
    const ccomplex_t IT_0316 = m_tau*IT_0315;
    const ccomplex_t IT_0317 = s_13*s_14;
    const ccomplex_t IT_0318 = s_34*IT_0021;
    const ccomplex_t IT_0319 = (-0.5)*IT_0318;
    const ccomplex_t IT_0320 = IT_0317 + IT_0319;
    const ccomplex_t IT_0321 = s_14*m_tau;
    const ccomplex_t IT_0322 = (-2)*conj(IT_0265);
    const ccomplex_t IT_0323 = IT_0016*IT_0021;
    const ccomplex_t IT_0324 = -conj(IT_0157);
    const ccomplex_t IT_0325 = -conj(IT_0238);
    const ccomplex_t IT_0326 = s_34 + IT_0016;
    const ccomplex_t IT_0327 = (-2)*IT_0265;
    const ccomplex_t IT_0328 = conj(IT_0314)*IT_0326;
    const ccomplex_t IT_0329 = 2*s_34;
    const ccomplex_t IT_0330 = conj(IT_0038)*IT_0280;
    const ccomplex_t IT_0331 = 2*IT_0038*((conj(IT_0157) + conj(IT_0238))
      *IT_0241 + 2*(conj(IT_0038) + 1./2*conj(IT_0265) + 1./2*conj(IT_0278))
      *IT_0280) + 2*IT_0238*(conj(IT_0038)*IT_0241 + conj(IT_0278)*IT_0281 +
       conj(IT_0314)*IT_0316 + 2*conj(IT_0238)*IT_0320 + 0.5*IT_0321*IT_0322 +
       IT_0323*IT_0324) + 2*IT_0157*(conj(IT_0038)*IT_0241 + conj(IT_0265)
      *IT_0281 + -conj(IT_0314)*IT_0316 + 2*conj(IT_0157)*IT_0320 + -conj
      (IT_0278)*IT_0321 + IT_0323*IT_0325) + (-2)*IT_0314*(IT_0316*(conj(IT_0157
      ) + IT_0325) + -(conj(IT_0278) + 2*conj(IT_0314) + 0.5*IT_0322)*IT_0326) +
       IT_0327*(IT_0016*conj(IT_0278) + conj(IT_0238)*IT_0321 + IT_0328) + 2
      *IT_0265*(conj(IT_0157)*IT_0281 + 0.5*conj(IT_0265)*IT_0329 + IT_0330) + 2
      *IT_0278*(conj(IT_0238)*IT_0281 + 0.5*IT_0016*IT_0322 + IT_0321*IT_0324 +
       IT_0328 + 0.5*conj(IT_0278)*IT_0329 + IT_0330);
    return create_ccomplex_return(IT_0331);
}

