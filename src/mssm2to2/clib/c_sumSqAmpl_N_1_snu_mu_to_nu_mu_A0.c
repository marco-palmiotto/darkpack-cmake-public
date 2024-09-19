#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_snu_mu_to_nu_mu_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_snu_mu_to_nu_mu_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B3)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W3)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0013 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = cos(beta);
    const ccomplex_t IT_0017 = IT_0001*IT_0016;
    const ccomplex_t IT_0018 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = sin(beta);
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = IT_0006*IT_0016;
    const ccomplex_t IT_0025 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = IT_0006*IT_0020;
    const ccomplex_t IT_0028 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0031 = IT_0017*IT_0030;
    const ccomplex_t IT_0032 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0033 = IT_0021*IT_0032;
    const ccomplex_t IT_0034 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0035 = IT_0024*IT_0034;
    const ccomplex_t IT_0036 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0037 = IT_0027*IT_0036;
    const ccomplex_t IT_0038 = (-0.5)*IT_0019 + 0.5*IT_0023 + 0.5*IT_0026 + (
      -0.5)*IT_0029 + (-0.5)*IT_0031 + 0.5*IT_0033 + 0.5*IT_0035 + (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = m_N_1*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0041 = IT_0017*IT_0040;
    const ccomplex_t IT_0042 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0043 = IT_0021*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0045 = IT_0024*IT_0044;
    const ccomplex_t IT_0046 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0047 = IT_0027*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0049 = IT_0017*IT_0048;
    const ccomplex_t IT_0050 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0051 = IT_0021*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0053 = IT_0024*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0055 = IT_0027*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0041 + (-0.5)*IT_0043 + (-0.5)*IT_0045 +
       0.5*IT_0047 + 0.5*IT_0049 + (-0.5)*IT_0051 + (-0.5)*IT_0053 + 0.5*IT_0055;
    const ccomplex_t IT_0057 = m_N_3*IT_0056;
    const ccomplex_t IT_0058 = IT_0039 + IT_0057;
    const ccomplex_t IT_0059 = IT_0015*IT_0058;
    const ccomplex_t IT_0060 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0061 = IT_0017*IT_0060;
    const ccomplex_t IT_0062 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0063 = IT_0021*IT_0062;
    const ccomplex_t IT_0064 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0065 = IT_0024*IT_0064;
    const ccomplex_t IT_0066 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0067 = IT_0027*IT_0066;
    const ccomplex_t IT_0068 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0069 = IT_0017*IT_0068;
    const ccomplex_t IT_0070 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0071 = IT_0021*IT_0070;
    const ccomplex_t IT_0072 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0073 = IT_0024*IT_0072;
    const ccomplex_t IT_0074 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0075 = IT_0027*IT_0074;
    const ccomplex_t IT_0076 = (-0.5)*IT_0061 + 0.5*IT_0063 + 0.5*IT_0065 + (
      -0.5)*IT_0067 + (-0.5)*IT_0069 + 0.5*IT_0071 + 0.5*IT_0073 + (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = m_N_1*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0079 = IT_0017*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0081 = IT_0021*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0083 = IT_0024*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0085 = IT_0027*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0087 = IT_0017*IT_0086;
    const ccomplex_t IT_0088 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0089 = IT_0021*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0091 = IT_0024*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0093 = IT_0027*IT_0092;
    const ccomplex_t IT_0094 = 0.5*IT_0079 + (-0.5)*IT_0081 + (-0.5)*IT_0083 +
       0.5*IT_0085 + 0.5*IT_0087 + (-0.5)*IT_0089 + (-0.5)*IT_0091 + 0.5*IT_0093;
    const ccomplex_t IT_0095 = m_N_2*IT_0094;
    const ccomplex_t IT_0096 = IT_0077 + IT_0095;
    const ccomplex_t IT_0097 = conj(N_B2)*e_em;
    const ccomplex_t IT_0098 = IT_0001*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = conj(N_W2)*e_em;
    const ccomplex_t IT_0101 = IT_0006*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0099 + -IT_0102);
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0096*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B4)*e_em;
    const ccomplex_t IT_0110 = IT_0001*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W4)*e_em;
    const ccomplex_t IT_0113 = IT_0006*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + -IT_0114);
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0123 = IT_0021*IT_0122;
    const ccomplex_t IT_0124 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0125 = IT_0024*IT_0124;
    const ccomplex_t IT_0126 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0127 = IT_0027*IT_0126;
    const ccomplex_t IT_0128 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0129 = IT_0017*IT_0128;
    const ccomplex_t IT_0130 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0131 = IT_0021*IT_0130;
    const ccomplex_t IT_0132 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0133 = IT_0024*IT_0132;
    const ccomplex_t IT_0134 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0135 = IT_0027*IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0121 + 0.5*IT_0123 + 0.5*IT_0125 + (
      -0.5)*IT_0127 + (-0.5)*IT_0129 + 0.5*IT_0131 + 0.5*IT_0133 + (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = m_N_1*IT_0136;
    const ccomplex_t IT_0138 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0139 = IT_0017*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0141 = IT_0021*IT_0140;
    const ccomplex_t IT_0142 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0143 = IT_0024*IT_0142;
    const ccomplex_t IT_0144 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0145 = IT_0027*IT_0144;
    const ccomplex_t IT_0146 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0147 = IT_0017*IT_0146;
    const ccomplex_t IT_0148 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0149 = IT_0021*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0151 = IT_0024*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0153 = IT_0027*IT_0152;
    const ccomplex_t IT_0154 = 0.5*IT_0139 + (-0.5)*IT_0141 + (-0.5)*IT_0143 +
       0.5*IT_0145 + 0.5*IT_0147 + (-0.5)*IT_0149 + (-0.5)*IT_0151 + 0.5*IT_0153;
    const ccomplex_t IT_0155 = m_N_4*IT_0154;
    const ccomplex_t IT_0156 = IT_0137 + IT_0155;
    const ccomplex_t IT_0157 = IT_0119*IT_0156;
    const ccomplex_t IT_0158 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0159 = IT_0017*IT_0158;
    const ccomplex_t IT_0160 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0161 = IT_0021*IT_0160;
    const ccomplex_t IT_0162 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0163 = IT_0024*IT_0162;
    const ccomplex_t IT_0164 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0165 = IT_0027*IT_0164;
    const ccomplex_t IT_0166 = -IT_0159 + IT_0161 + IT_0163 + -IT_0165;
    const ccomplex_t IT_0167 = m_N_1*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0169 = IT_0017*IT_0168;
    const ccomplex_t IT_0170 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0171 = IT_0021*IT_0170;
    const ccomplex_t IT_0172 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0173 = IT_0024*IT_0172;
    const ccomplex_t IT_0174 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0175 = IT_0027*IT_0174;
    const ccomplex_t IT_0176 = IT_0169 + -IT_0171 + -IT_0173 + IT_0175;
    const ccomplex_t IT_0177 = m_N_1*IT_0176;
    const ccomplex_t IT_0178 = IT_0167 + IT_0177;
    const ccomplex_t IT_0179 = conj(N_B1)*e_em;
    const ccomplex_t IT_0180 = IT_0001*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = conj(N_W1)*e_em;
    const ccomplex_t IT_0183 = IT_0006*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0181 + -IT_0184);
    const ccomplex_t IT_0186 = (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0188 = IT_0186*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0178*IT_0189;
    const ccomplex_t IT_0191 = 0.5*IT_0059 + 0.5*IT_0108 + 0.5*IT_0157 + 0.5
      *IT_0190;
    const ccomplex_t IT_0192 = IT_0076*IT_0107;
    const ccomplex_t IT_0193 = IT_0015*IT_0038;
    const ccomplex_t IT_0194 = IT_0119*IT_0136;
    const ccomplex_t IT_0195 = IT_0166*IT_0189;
    const ccomplex_t IT_0196 = (-0.5)*IT_0192 + (-0.5)*IT_0193 + (-0.5)
      *IT_0194 + (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = s_34*m_N_1;
    const ccomplex_t IT_0198 = pow(m_A0, 2);
    const ccomplex_t IT_0199 = s_13*IT_0198;
    const ccomplex_t IT_0200 = s_14*s_34;
    return create_ccomplex_return(2*IT_0191*(s_13*conj(IT_0191) + conj(IT_0196
      )*IT_0197) + 2*IT_0196*(conj(IT_0191)*IT_0197 + -conj(IT_0196)*(IT_0199 + 
      (-2)*IT_0200)));
}

