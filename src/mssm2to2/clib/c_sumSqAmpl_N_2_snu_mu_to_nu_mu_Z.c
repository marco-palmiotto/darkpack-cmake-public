#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_snu_mu_to_nu_mu_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_snu_mu_to_nu_mu_Z(
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
    const ccomplex_t IT_0001 = pow(m_Z, 2);
    const ccomplex_t IT_0002 = s_13*IT_0001;
    const ccomplex_t IT_0003 = s_14*s_34;
    const ccomplex_t IT_0004 = (-2)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = s_13 + IT_0007;
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = conj(N_B2)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = conj(N_W2)*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0013 + -IT_0018);
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0022 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0010*IT_0014;
    const ccomplex_t IT_0026 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = IT_0009*IT_0015;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0031 = IT_0025*IT_0030;
    const ccomplex_t IT_0032 = IT_0028*IT_0030;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0027 + IT_0029 + 
      -IT_0031 + -IT_0032);
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = m_N_2*IT_0034;
    const ccomplex_t IT_0036 = IT_0024*IT_0035;
    const ccomplex_t IT_0037 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0038 = IT_0025*IT_0037;
    const ccomplex_t IT_0039 = IT_0028*IT_0037;
    const ccomplex_t IT_0040 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0041 = IT_0025*IT_0040;
    const ccomplex_t IT_0042 = IT_0028*IT_0040;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0038 + IT_0039 + 
      -IT_0041 + -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = m_N_1*IT_0044;
    const ccomplex_t IT_0046 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0047 = IT_0025*IT_0046;
    const ccomplex_t IT_0048 = IT_0028*IT_0046;
    const ccomplex_t IT_0049 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0050 = IT_0025*IT_0049;
    const ccomplex_t IT_0051 = IT_0028*IT_0049;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0047 + IT_0048 + 
      -IT_0050 + -IT_0051);
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = m_N_2*IT_0053;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0045 + IT_0055;
    const ccomplex_t IT_0057 = conj(N_B1)*e_em;
    const ccomplex_t IT_0058 = IT_0010*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W1)*e_em;
    const ccomplex_t IT_0061 = IT_0015*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + -IT_0062);
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0056*IT_0067;
    const ccomplex_t IT_0069 = (-0.25)*IT_0036 + (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0033;
    const ccomplex_t IT_0071 = m_N_2*IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0035 + IT_0072;
    const ccomplex_t IT_0074 = IT_0024*IT_0073;
    const ccomplex_t IT_0075 = IT_0024*IT_0071;
    const ccomplex_t IT_0076 = conj(N_B3)*e_em;
    const ccomplex_t IT_0077 = IT_0010*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = conj(N_W3)*e_em;
    const ccomplex_t IT_0080 = IT_0015*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0078 + -IT_0081);
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0088 = IT_0025*IT_0087;
    const ccomplex_t IT_0089 = IT_0028*IT_0087;
    const ccomplex_t IT_0090 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0091 = IT_0025*IT_0090;
    const ccomplex_t IT_0092 = IT_0028*IT_0090;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0088 + IT_0089 + 
      -IT_0091 + -IT_0092);
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = m_N_3*IT_0094;
    const ccomplex_t IT_0096 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0097 = IT_0025*IT_0096;
    const ccomplex_t IT_0098 = IT_0028*IT_0096;
    const ccomplex_t IT_0099 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0100 = IT_0025*IT_0099;
    const ccomplex_t IT_0101 = IT_0028*IT_0099;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0097 + IT_0098 + 
      -IT_0100 + -IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = m_N_2*IT_0103;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = IT_0095 + IT_0105;
    const ccomplex_t IT_0107 = IT_0086*IT_0106;
    const ccomplex_t IT_0108 = conj(N_B4)*e_em;
    const ccomplex_t IT_0109 = IT_0010*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = conj(N_W4)*e_em;
    const ccomplex_t IT_0112 = IT_0015*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + -IT_0113);
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0120 = IT_0025*IT_0119;
    const ccomplex_t IT_0121 = IT_0028*IT_0119;
    const ccomplex_t IT_0122 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0123 = IT_0025*IT_0122;
    const ccomplex_t IT_0124 = IT_0028*IT_0122;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0120 + IT_0121 + 
      -IT_0123 + -IT_0124);
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = m_N_4*IT_0126;
    const ccomplex_t IT_0128 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0129 = IT_0028*IT_0128;
    const ccomplex_t IT_0130 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0131 = IT_0025*IT_0130;
    const ccomplex_t IT_0132 = IT_0028*IT_0130;
    const ccomplex_t IT_0133 = IT_0025*IT_0128;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0129 + -IT_0131 + 
      -IT_0132 + IT_0133);
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = m_N_2*IT_0136;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0127 + IT_0138;
    const ccomplex_t IT_0140 = IT_0118*IT_0139;
    const ccomplex_t IT_0141 = (-0.25)*IT_0074 + 0.25*IT_0075 + 0.5*IT_0107 +
       0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0069 + IT_0141;
    const ccomplex_t IT_0143 = conj(IT_0069) + conj(IT_0141);
    const ccomplex_t IT_0144 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = s_13 + IT_0145;
    const ccomplex_t IT_0147 = s_13*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*e_em*IT_0010*IT_0015;
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = pow(m_N_2, 2);
    const ccomplex_t IT_0151 = cpow(s_13 + (-0.5)*IT_0150 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0152 = IT_0020*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = IT_0149*IT_0153;
    const ccomplex_t IT_0155 = cpow(s_12 + 0.5*IT_0021 + 0.5*IT_0150 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0156 = e_em*IT_0025;
    const ccomplex_t IT_0157 = e_em*IT_0028;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0156 + IT_0157);
    const ccomplex_t IT_0159 = 0.5*IT_0158;
    const ccomplex_t IT_0160 = IT_0020*IT_0159;
    const ccomplex_t IT_0161 = IT_0155*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = (-0.5)*IT_0154 + -IT_0162;
    const ccomplex_t IT_0164 = IT_0053*IT_0067;
    const ccomplex_t IT_0165 = 2*IT_0164;
    const ccomplex_t IT_0166 = IT_0086*IT_0103;
    const ccomplex_t IT_0167 = IT_0118*IT_0136;
    const ccomplex_t IT_0168 = IT_0024*IT_0070;
    const ccomplex_t IT_0169 = IT_0024*IT_0034;
    const ccomplex_t IT_0170 = 0.5*IT_0154 + IT_0166 + IT_0167 + (-0.5)
      *IT_0168 + 0.5*IT_0169;
    const ccomplex_t IT_0171 = (-2)*IT_0170;
    const ccomplex_t IT_0172 = IT_0165 + IT_0171;
    const ccomplex_t IT_0173 = 2*conj(IT_0164);
    const ccomplex_t IT_0174 = (-2)*conj(IT_0170);
    const ccomplex_t IT_0175 = IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = m_N_2*IT_0146;
    const ccomplex_t IT_0177 = s_34*IT_0001;
    const ccomplex_t IT_0178 = IT_0000*IT_0177;
    const ccomplex_t IT_0179 = (-0.25)*IT_0178;
    const ccomplex_t IT_0180 = s_34 + IT_0179;
    const ccomplex_t IT_0181 = m_N_2*IT_0180;
    const ccomplex_t IT_0182 = (-4)*IT_0181;
    const ccomplex_t IT_0183 = (-0.5)*IT_0164;
    const ccomplex_t IT_0184 = (-0.5)*IT_0166 + (-0.5)*IT_0167 + 0.25*IT_0168 
      + (-0.25)*IT_0169;
    const ccomplex_t IT_0185 = 2*IT_0069;
    const ccomplex_t IT_0186 = 2*IT_0141;
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = 2*conj(IT_0069);
    const ccomplex_t IT_0189 = 2*conj(IT_0141);
    const ccomplex_t IT_0190 = IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = pow(s_34, 2);
    const ccomplex_t IT_0192 = s_13*IT_0000*IT_0191;
    const ccomplex_t IT_0193 = 2*IT_0163;
    const ccomplex_t IT_0194 = s_13*s_34;
    const ccomplex_t IT_0195 = (-0.5)*IT_0162;
    const ccomplex_t IT_0196 = conj(IT_0183) + conj(IT_0195);
    const ccomplex_t IT_0197 = IT_0183 + IT_0195;
    const ccomplex_t IT_0198 = 2*conj(IT_0163);
    const ccomplex_t IT_0199 = m_N_2*IT_0000*IT_0191;
    const ccomplex_t IT_0200 = pow(s_14, 2);
    const ccomplex_t IT_0201 = IT_0000*IT_0200;
    const ccomplex_t IT_0202 = -IT_0150;
    const ccomplex_t IT_0203 = IT_0201 + IT_0202;
    const ccomplex_t IT_0204 = s_13*IT_0203;
    const ccomplex_t IT_0205 = 2*IT_0170;
    const ccomplex_t IT_0206 = (-0.5)*IT_0205;
    const ccomplex_t IT_0207 = IT_0164 + IT_0206;
    const ccomplex_t IT_0208 = (-2)*conj(IT_0164);
    const ccomplex_t IT_0209 = s_13*s_14;
    const ccomplex_t IT_0210 = s_34*IT_0150;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + IT_0211;
    const ccomplex_t IT_0213 = -conj(IT_0184);
    const ccomplex_t IT_0214 = pow(m_Z, 4);
    const ccomplex_t IT_0215 = s_13*IT_0214;
    const ccomplex_t IT_0216 = IT_0000*IT_0215;
    const ccomplex_t IT_0217 = (-0.5)*IT_0216;
    const ccomplex_t IT_0218 = IT_0002 + IT_0004 + IT_0217;
    const ccomplex_t IT_0219 = (-0.25)*IT_0216;
    const ccomplex_t IT_0220 = IT_0003 + IT_0219;
    const ccomplex_t IT_0221 = 0.5*IT_0178;
    const ccomplex_t IT_0222 = s_34 + IT_0221;
    const ccomplex_t IT_0223 = m_N_2*IT_0222;
    const ccomplex_t IT_0224 = 2*IT_0223;
    const ccomplex_t IT_0225 = 0.125*IT_0224;
    const ccomplex_t IT_0226 = 4*IT_0008*IT_0142*IT_0143 + IT_0147*(conj
      (IT_0163)*IT_0172 + IT_0163*IT_0175) + (IT_0143*IT_0172 + IT_0142*IT_0175)
      *IT_0176 + (-2)*IT_0182*(IT_0143*IT_0183 + IT_0142*conj(IT_0183) + (-0.5)
      *conj(IT_0184)*IT_0187 + (-0.5)*IT_0184*IT_0190) + conj(IT_0163)*IT_0192
      *IT_0193 + (-2)*IT_0194*(conj(IT_0163)*IT_0184 + IT_0163*conj(IT_0184) + (
      -0.5)*IT_0193*IT_0196 + (-0.5)*IT_0197*IT_0198) + (IT_0143*IT_0193 +
       IT_0142*IT_0198)*IT_0199 + (-2)*IT_0204*((-0.5)*conj(IT_0164)*IT_0165 +
       conj(IT_0170)*IT_0207 + (-0.5)*IT_0170*IT_0208) + IT_0212*(IT_0172
      *IT_0196 + IT_0175*IT_0197 + 2*IT_0184*(conj(IT_0170) + 0.5*IT_0208) + 2
      *IT_0207*IT_0213) + (-4)*(IT_0195*conj(IT_0195) + (IT_0183 + -IT_0184)*
      (conj(IT_0183) + IT_0213))*IT_0218 + 8*conj(IT_0195)*((IT_0183 + -IT_0184)
      *IT_0220 + IT_0187*IT_0225) + 8*IT_0195*((conj(IT_0183) + IT_0213)*IT_0220
       + IT_0190*IT_0225);
    return create_ccomplex_return(IT_0226);
}

