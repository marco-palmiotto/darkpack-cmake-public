#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_snu_mu_to_nu_mu_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_snu_mu_to_nu_mu_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
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
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
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
    const ccomplex_t IT_0000 = pow(m_Z, -2);
    const ccomplex_t IT_0001 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0002 = -IT_0001;
    const ccomplex_t IT_0003 = s_13 + IT_0002;
    const ccomplex_t IT_0004 = m_N_3*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_B3)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = conj(N_W3)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0009 + -IT_0014);
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0018 = cpow(s_23 + (-0.5)*IT_0017 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0006*IT_0010;
    const ccomplex_t IT_0022 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = IT_0005*IT_0011;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0027 = IT_0021*IT_0026;
    const ccomplex_t IT_0028 = IT_0024*IT_0026;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0023 + IT_0025 + 
      -IT_0027 + -IT_0028);
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = m_N_3*IT_0030;
    const ccomplex_t IT_0032 = (-0.5)*IT_0029;
    const ccomplex_t IT_0033 = m_N_3*IT_0032;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = IT_0031 + IT_0034;
    const ccomplex_t IT_0036 = IT_0020*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B1)*e_em;
    const ccomplex_t IT_0038 = IT_0006*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W1)*e_em;
    const ccomplex_t IT_0041 = IT_0011*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = cpow(s_23 + (-0.5)*IT_0017 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0049 = IT_0021*IT_0048;
    const ccomplex_t IT_0050 = IT_0024*IT_0048;
    const ccomplex_t IT_0051 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0052 = IT_0021*IT_0051;
    const ccomplex_t IT_0053 = IT_0024*IT_0051;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0049 + IT_0050 + 
      -IT_0052 + -IT_0053);
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = m_N_1*IT_0055;
    const ccomplex_t IT_0057 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0058 = IT_0021*IT_0057;
    const ccomplex_t IT_0059 = IT_0024*IT_0057;
    const ccomplex_t IT_0060 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0061 = IT_0021*IT_0060;
    const ccomplex_t IT_0062 = IT_0024*IT_0060;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0058 + IT_0059 + 
      -IT_0061 + -IT_0062);
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = m_N_3*IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0056 + IT_0066;
    const ccomplex_t IT_0068 = IT_0047*IT_0067;
    const ccomplex_t IT_0069 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0070 = IT_0021*IT_0069;
    const ccomplex_t IT_0071 = IT_0024*IT_0069;
    const ccomplex_t IT_0072 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0073 = IT_0021*IT_0072;
    const ccomplex_t IT_0074 = IT_0024*IT_0072;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0070 + IT_0071 + 
      -IT_0073 + -IT_0074);
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = m_N_2*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B2)*e_em;
    const ccomplex_t IT_0079 = IT_0006*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = conj(N_W2)*e_em;
    const ccomplex_t IT_0082 = IT_0011*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = cpow(s_23 + (-0.5)*IT_0017 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0077*IT_0088;
    const ccomplex_t IT_0090 = (-0.25)*IT_0036 + (-0.5)*IT_0068 + (-0.5)
      *IT_0089;
    const ccomplex_t IT_0091 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0092 = IT_0021*IT_0091;
    const ccomplex_t IT_0093 = IT_0024*IT_0091;
    const ccomplex_t IT_0094 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0095 = IT_0021*IT_0094;
    const ccomplex_t IT_0096 = IT_0024*IT_0094;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0092 + IT_0093 + 
      -IT_0095 + -IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = m_N_3*IT_0098;
    const ccomplex_t IT_0100 = IT_0088*IT_0099;
    const ccomplex_t IT_0101 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0102 = IT_0021*IT_0101;
    const ccomplex_t IT_0103 = IT_0024*IT_0101;
    const ccomplex_t IT_0104 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0105 = IT_0021*IT_0104;
    const ccomplex_t IT_0106 = IT_0024*IT_0104;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0102 + IT_0103 + 
      -IT_0105 + -IT_0106);
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = m_N_4*IT_0108;
    const ccomplex_t IT_0110 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0111 = IT_0021*IT_0110;
    const ccomplex_t IT_0112 = IT_0024*IT_0110;
    const ccomplex_t IT_0113 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0114 = IT_0021*IT_0113;
    const ccomplex_t IT_0115 = IT_0024*IT_0113;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0111 + IT_0112 + 
      -IT_0114 + -IT_0115);
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = m_N_3*IT_0117;
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = IT_0109 + IT_0119;
    const ccomplex_t IT_0121 = conj(N_B4)*e_em;
    const ccomplex_t IT_0122 = IT_0006*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W4)*e_em;
    const ccomplex_t IT_0125 = IT_0011*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0123 + -IT_0126);
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = cpow(s_23 + (-0.5)*IT_0017 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = IT_0120*IT_0131;
    const ccomplex_t IT_0133 = (-0.25)*IT_0036 + 0.5*IT_0100 + 0.5*IT_0132;
    const ccomplex_t IT_0134 = conj(IT_0090) + conj(IT_0133);
    const ccomplex_t IT_0135 = IT_0020*IT_0032;
    const ccomplex_t IT_0136 = IT_0047*IT_0064;
    const ccomplex_t IT_0137 = 0.5*IT_0135 + IT_0136;
    const ccomplex_t IT_0138 = 2*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0011;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = pow(m_N_3, 2);
    const ccomplex_t IT_0142 = cpow(s_13 + (-0.5)*IT_0141 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0143 = IT_0016*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0140*IT_0144;
    const ccomplex_t IT_0146 = IT_0088*IT_0098;
    const ccomplex_t IT_0147 = IT_0117*IT_0131;
    const ccomplex_t IT_0148 = IT_0020*IT_0030;
    const ccomplex_t IT_0149 = 0.5*IT_0145 + -IT_0146 + IT_0147 + 0.5*IT_0148;
    const ccomplex_t IT_0150 = (-2)*IT_0149;
    const ccomplex_t IT_0151 = IT_0138 + IT_0150;
    const ccomplex_t IT_0152 = s_13*IT_0003;
    const ccomplex_t IT_0153 = cpow(s_12 + 0.5*IT_0017 + 0.5*IT_0141 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0154 = e_em*IT_0021;
    const ccomplex_t IT_0155 = e_em*IT_0024;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0154 + IT_0155);
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = IT_0016*IT_0157;
    const ccomplex_t IT_0159 = IT_0153*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (-0.5)*IT_0145 + -IT_0160;
    const ccomplex_t IT_0162 = 2*conj(IT_0137);
    const ccomplex_t IT_0163 = (-2)*conj(IT_0149);
    const ccomplex_t IT_0164 = IT_0162 + IT_0163;
    const ccomplex_t IT_0165 = pow(s_34, 2);
    const ccomplex_t IT_0166 = s_13*IT_0000*IT_0165;
    const ccomplex_t IT_0167 = m_N_3*IT_0000*IT_0165;
    const ccomplex_t IT_0168 = s_13*s_34;
    const ccomplex_t IT_0169 = (-0.5)*IT_0160;
    const ccomplex_t IT_0170 = (-0.25)*IT_0135 + (-0.5)*IT_0136;
    const ccomplex_t IT_0171 = conj(IT_0169) + conj(IT_0170);
    const ccomplex_t IT_0172 = 2*IT_0161;
    const ccomplex_t IT_0173 = pow(m_Z, 2);
    const ccomplex_t IT_0174 = s_34*IT_0173;
    const ccomplex_t IT_0175 = IT_0000*IT_0174;
    const ccomplex_t IT_0176 = 0.5*IT_0175;
    const ccomplex_t IT_0177 = s_34 + IT_0176;
    const ccomplex_t IT_0178 = m_N_3*IT_0177;
    const ccomplex_t IT_0179 = 2*IT_0178;
    const ccomplex_t IT_0180 = s_13*IT_0173;
    const ccomplex_t IT_0181 = s_14*s_34;
    const ccomplex_t IT_0182 = (-2)*IT_0181;
    const ccomplex_t IT_0183 = IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = IT_0000*IT_0183;
    const ccomplex_t IT_0185 = (-0.5)*IT_0184;
    const ccomplex_t IT_0186 = s_13 + IT_0185;
    const ccomplex_t IT_0187 = (-0.25)*IT_0175;
    const ccomplex_t IT_0188 = s_34 + IT_0187;
    const ccomplex_t IT_0189 = m_N_3*IT_0188;
    const ccomplex_t IT_0190 = (-4)*IT_0189;
    const ccomplex_t IT_0191 = 0.5*IT_0146 + (-0.5)*IT_0147 + (-0.25)*IT_0148;
    const ccomplex_t IT_0192 = -conj(IT_0191);
    const ccomplex_t IT_0193 = conj(IT_0170) + IT_0192;
    const ccomplex_t IT_0194 = 2*conj(IT_0161);
    const ccomplex_t IT_0195 = pow(m_Z, 4);
    const ccomplex_t IT_0196 = s_13*IT_0195;
    const ccomplex_t IT_0197 = IT_0000*IT_0196;
    const ccomplex_t IT_0198 = (-0.25)*IT_0197;
    const ccomplex_t IT_0199 = IT_0181 + IT_0198;
    const ccomplex_t IT_0200 = (-0.5)*IT_0197;
    const ccomplex_t IT_0201 = IT_0180 + IT_0182 + IT_0200;
    const ccomplex_t IT_0202 = 2*IT_0133;
    const ccomplex_t IT_0203 = pow(s_14, 2);
    const ccomplex_t IT_0204 = IT_0000*IT_0203;
    const ccomplex_t IT_0205 = -IT_0141;
    const ccomplex_t IT_0206 = IT_0204 + IT_0205;
    const ccomplex_t IT_0207 = s_13*IT_0206;
    const ccomplex_t IT_0208 = 2*IT_0149;
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = IT_0137 + IT_0209;
    const ccomplex_t IT_0211 = (-2)*conj(IT_0137);
    const ccomplex_t IT_0212 = s_13*s_14;
    const ccomplex_t IT_0213 = s_34*IT_0141;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = IT_0212 + IT_0214;
    const ccomplex_t IT_0216 = 2*conj(IT_0090);
    const ccomplex_t IT_0217 = 2*conj(IT_0133);
    const ccomplex_t IT_0218 = IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = (-2)*IT_0193;
    const ccomplex_t IT_0220 = (-0.5)*IT_0190;
    const ccomplex_t IT_0221 = IT_0004*IT_0134*IT_0151 + IT_0152*(IT_0151*conj
      (IT_0161) + IT_0161*IT_0164) + (conj(IT_0161)*IT_0166 + IT_0134*IT_0167 +
       IT_0168*IT_0171)*IT_0172 + 2*IT_0090*(0.5*IT_0004*IT_0164 + conj(IT_0169)
      *IT_0179 + 2*IT_0134*IT_0186 + -IT_0190*IT_0193 + 0.5*IT_0167*IT_0194) + (
      -2)*IT_0133*((-0.5)*IT_0004*IT_0164 + (-2)*IT_0134*IT_0186 + conj(IT_0170)
      *IT_0190 + (-0.5)*IT_0167*IT_0194) + 8*IT_0170*((-0.25)*IT_0134*IT_0190 +
       0.125*IT_0168*IT_0194 + conj(IT_0169)*IT_0199 + (-0.5)*IT_0193*IT_0201) +
       conj(IT_0169)*IT_0179*IT_0202 + (-2)*IT_0207*((-0.5)*conj(IT_0137)
      *IT_0138 + conj(IT_0149)*IT_0210 + (-0.5)*IT_0149*IT_0211) + (IT_0164*
      (IT_0169 + IT_0170) + IT_0151*IT_0171 + 2*IT_0192*IT_0210 + 2*IT_0191*
      (conj(IT_0149) + 0.5*IT_0211))*IT_0215 + (-4)*IT_0169*((-0.25)*IT_0168
      *IT_0194 + conj(IT_0169)*IT_0201 + (-0.25)*IT_0179*IT_0218 + IT_0199
      *IT_0219) + (-2)*conj(IT_0191)*(IT_0161*IT_0168 + IT_0202*IT_0220) + (-2)
      *IT_0191*(conj(IT_0161)*IT_0168 + 4*conj(IT_0169)*IT_0199 + IT_0201
      *IT_0219 + IT_0218*IT_0220);
    return create_ccomplex_return(IT_0221);
}

