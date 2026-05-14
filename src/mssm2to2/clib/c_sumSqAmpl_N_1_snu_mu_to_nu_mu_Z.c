#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_snu_mu_to_nu_mu_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_snu_mu_to_nu_mu_Z(
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
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0011, -1);
    const ccomplex_t IT_0016 = IT_0009*IT_0015;
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0019 = IT_0012*IT_0018;
    const ccomplex_t IT_0020 = IT_0016*IT_0018;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0014 + IT_0017 + 
      -IT_0019 + -IT_0020);
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = m_N_1*IT_0022;
    const ccomplex_t IT_0024 = (-0.5)*IT_0021;
    const ccomplex_t IT_0025 = m_N_1*IT_0024;
    const ccomplex_t IT_0026 = conj(N_B1)*e_em;
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = conj(N_W1)*e_em;
    const ccomplex_t IT_0030 = IT_0015*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + -IT_0031);
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0035 = cpow(s_23 + (-0.5)*IT_0034 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (IT_0023 + -IT_0025)*IT_0037;
    const ccomplex_t IT_0039 = (-0.25)*IT_0038;
    const ccomplex_t IT_0040 = -IT_0025;
    const ccomplex_t IT_0041 = IT_0023 + IT_0040;
    const ccomplex_t IT_0042 = IT_0037*IT_0041;
    const ccomplex_t IT_0043 = conj(N_B3)*e_em;
    const ccomplex_t IT_0044 = IT_0010*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = conj(N_W3)*e_em;
    const ccomplex_t IT_0047 = IT_0015*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + -IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = cpow(s_23 + (-0.5)*IT_0034 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0055 = IT_0012*IT_0054;
    const ccomplex_t IT_0056 = IT_0016*IT_0054;
    const ccomplex_t IT_0057 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0058 = IT_0012*IT_0057;
    const ccomplex_t IT_0059 = IT_0016*IT_0057;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0055 + IT_0056 + 
      -IT_0058 + -IT_0059);
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = m_N_3*IT_0061;
    const ccomplex_t IT_0063 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0064 = IT_0012*IT_0063;
    const ccomplex_t IT_0065 = IT_0016*IT_0063;
    const ccomplex_t IT_0066 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0067 = IT_0012*IT_0066;
    const ccomplex_t IT_0068 = IT_0016*IT_0066;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0064 + IT_0065 + 
      -IT_0067 + -IT_0068);
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = m_N_1*IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0062 + IT_0072;
    const ccomplex_t IT_0074 = IT_0053*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B4)*e_em;
    const ccomplex_t IT_0076 = IT_0010*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W4)*e_em;
    const ccomplex_t IT_0079 = IT_0015*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + -IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = cpow(s_23 + (-0.5)*IT_0034 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0087 = IT_0012*IT_0086;
    const ccomplex_t IT_0088 = IT_0016*IT_0086;
    const ccomplex_t IT_0089 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0090 = IT_0012*IT_0089;
    const ccomplex_t IT_0091 = IT_0016*IT_0089;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0087 + IT_0088 + 
      -IT_0090 + -IT_0091);
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = m_N_4*IT_0093;
    const ccomplex_t IT_0095 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0096 = IT_0012*IT_0095;
    const ccomplex_t IT_0097 = IT_0016*IT_0095;
    const ccomplex_t IT_0098 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0099 = IT_0012*IT_0098;
    const ccomplex_t IT_0100 = IT_0016*IT_0098;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0096 + IT_0097 + 
      -IT_0099 + -IT_0100);
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = m_N_1*IT_0102;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = IT_0094 + IT_0104;
    const ccomplex_t IT_0106 = IT_0085*IT_0105;
    const ccomplex_t IT_0107 = conj(N_B2)*e_em;
    const ccomplex_t IT_0108 = IT_0010*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = conj(N_W2)*e_em;
    const ccomplex_t IT_0111 = IT_0015*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + -IT_0112);
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = cpow(s_23 + (-0.5)*IT_0034 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0119 = IT_0012*IT_0118;
    const ccomplex_t IT_0120 = IT_0016*IT_0118;
    const ccomplex_t IT_0121 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0122 = IT_0012*IT_0121;
    const ccomplex_t IT_0123 = IT_0016*IT_0121;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0119 + IT_0120 + 
      -IT_0122 + -IT_0123);
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = m_N_2*IT_0125;
    const ccomplex_t IT_0127 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0128 = IT_0012*IT_0127;
    const ccomplex_t IT_0129 = IT_0016*IT_0127;
    const ccomplex_t IT_0130 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0131 = IT_0012*IT_0130;
    const ccomplex_t IT_0132 = IT_0016*IT_0130;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0128 + IT_0129 + 
      -IT_0131 + -IT_0132);
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = m_N_1*IT_0134;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = IT_0126 + IT_0136;
    const ccomplex_t IT_0138 = IT_0117*IT_0137;
    const ccomplex_t IT_0139 = (-0.25)*IT_0042 + 0.5*IT_0074 + 0.5*IT_0106 +
       0.5*IT_0138;
    const ccomplex_t IT_0140 = IT_0039 + IT_0139;
    const ccomplex_t IT_0141 = conj(IT_0039) + conj(IT_0139);
    const ccomplex_t IT_0142 = pow(m_Z, 4);
    const ccomplex_t IT_0143 = s_13*IT_0142;
    const ccomplex_t IT_0144 = IT_0000*IT_0143;
    const ccomplex_t IT_0145 = (-0.5)*IT_0144;
    const ccomplex_t IT_0146 = IT_0002 + IT_0004 + IT_0145;
    const ccomplex_t IT_0147 = e_em*IT_0012;
    const ccomplex_t IT_0148 = e_em*IT_0016;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0147 + IT_0148);
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0033*IT_0150;
    const ccomplex_t IT_0152 = pow(m_N_1, 2);
    const ccomplex_t IT_0153 = cpow(s_12 + 0.5*IT_0034 + 0.5*IT_0152 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0154 = IT_0151*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = IT_0024*IT_0037;
    const ccomplex_t IT_0158 = (-0.25)*IT_0157;
    const ccomplex_t IT_0159 = IT_0117*IT_0134;
    const ccomplex_t IT_0160 = IT_0022*IT_0037;
    const ccomplex_t IT_0161 = IT_0053*IT_0070;
    const ccomplex_t IT_0162 = IT_0085*IT_0102;
    const ccomplex_t IT_0163 = (-0.5)*IT_0159 + (-0.25)*IT_0160 + (-0.5)
      *IT_0161 + (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = IT_0158 + -IT_0163;
    const ccomplex_t IT_0165 = -conj(IT_0163);
    const ccomplex_t IT_0166 = conj(IT_0158) + IT_0165;
    const ccomplex_t IT_0167 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = s_13 + IT_0168;
    const ccomplex_t IT_0170 = s_13*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*e_em*IT_0010*IT_0015;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = cpow(s_13 + (-0.5)*IT_0152 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0174 = IT_0033*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = IT_0172*IT_0175;
    const ccomplex_t IT_0177 = -IT_0155 + (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = 0.5*IT_0157;
    const ccomplex_t IT_0179 = 2*IT_0178;
    const ccomplex_t IT_0180 = IT_0159 + 0.5*IT_0160 + IT_0161 + IT_0162 + 0.5
      *IT_0176;
    const ccomplex_t IT_0181 = (-2)*IT_0180;
    const ccomplex_t IT_0182 = IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = 2*conj(IT_0178);
    const ccomplex_t IT_0184 = (-2)*conj(IT_0180);
    const ccomplex_t IT_0185 = IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = m_N_1*IT_0169;
    const ccomplex_t IT_0187 = s_34*IT_0001;
    const ccomplex_t IT_0188 = IT_0000*IT_0187;
    const ccomplex_t IT_0189 = (-0.25)*IT_0188;
    const ccomplex_t IT_0190 = s_34 + IT_0189;
    const ccomplex_t IT_0191 = m_N_1*IT_0190;
    const ccomplex_t IT_0192 = (-4)*IT_0191;
    const ccomplex_t IT_0193 = 2*IT_0039;
    const ccomplex_t IT_0194 = 2*IT_0139;
    const ccomplex_t IT_0195 = IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = 2*conj(IT_0039);
    const ccomplex_t IT_0197 = 2*conj(IT_0139);
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = pow(s_34, 2);
    const ccomplex_t IT_0200 = s_13*IT_0000*IT_0199;
    const ccomplex_t IT_0201 = 2*IT_0177;
    const ccomplex_t IT_0202 = m_N_1*IT_0000*IT_0199;
    const ccomplex_t IT_0203 = 2*conj(IT_0177);
    const ccomplex_t IT_0204 = s_13*s_34;
    const ccomplex_t IT_0205 = conj(IT_0156) + conj(IT_0158);
    const ccomplex_t IT_0206 = IT_0156 + IT_0158;
    const ccomplex_t IT_0207 = (-0.25)*IT_0144;
    const ccomplex_t IT_0208 = IT_0003 + IT_0207;
    const ccomplex_t IT_0209 = 0.5*IT_0188;
    const ccomplex_t IT_0210 = s_34 + IT_0209;
    const ccomplex_t IT_0211 = m_N_1*IT_0210;
    const ccomplex_t IT_0212 = 2*IT_0211;
    const ccomplex_t IT_0213 = 0.125*IT_0212;
    const ccomplex_t IT_0214 = pow(s_14, 2);
    const ccomplex_t IT_0215 = IT_0000*IT_0214;
    const ccomplex_t IT_0216 = -IT_0152;
    const ccomplex_t IT_0217 = IT_0215 + IT_0216;
    const ccomplex_t IT_0218 = s_13*IT_0217;
    const ccomplex_t IT_0219 = (-2)*conj(IT_0178);
    const ccomplex_t IT_0220 = 2*IT_0180;
    const ccomplex_t IT_0221 = (-0.5)*IT_0220;
    const ccomplex_t IT_0222 = s_13*s_14;
    const ccomplex_t IT_0223 = s_34*IT_0152;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0222 + IT_0224;
    const ccomplex_t IT_0226 = 4*IT_0008*IT_0140*IT_0141 + (-4)*IT_0146*
      (IT_0156*conj(IT_0156) + IT_0164*IT_0166) + IT_0170*(conj(IT_0177)*IT_0182
       + IT_0177*IT_0185) + (IT_0141*IT_0182 + IT_0140*IT_0185)*IT_0186 + (-2)
      *IT_0192*(IT_0141*IT_0158 + IT_0140*conj(IT_0158) + (-0.5)*conj(IT_0163)
      *IT_0195 + (-0.5)*IT_0163*IT_0198) + conj(IT_0177)*IT_0200*IT_0201 +
       IT_0202*(IT_0141*IT_0201 + IT_0140*IT_0203) + (-2)*IT_0204*(conj(IT_0163)
      *IT_0177 + IT_0163*conj(IT_0177) + (-0.5)*IT_0201*IT_0205 + (-0.5)*IT_0203
      *IT_0206) + 8*conj(IT_0156)*(IT_0164*IT_0208 + IT_0195*IT_0213) + 8
      *IT_0156*(IT_0166*IT_0208 + IT_0198*IT_0213) + IT_0218*(conj(IT_0178)
      *IT_0179 + IT_0180*IT_0219 + (-2)*conj(IT_0180)*(IT_0178 + IT_0221)) + 
      (IT_0182*IT_0205 + IT_0185*IT_0206 + 2*IT_0163*(conj(IT_0180) + 0.5
      *IT_0219) + 2*IT_0165*(IT_0178 + IT_0221))*IT_0225;
    return create_ccomplex_return(IT_0226);
}

