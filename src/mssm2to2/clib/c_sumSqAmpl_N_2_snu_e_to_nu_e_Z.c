#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_snu_e_to_nu_e_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_snu_e_to_nu_e_Z(
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
    const creal_t Gamma_nuel = param->Gamma_nuel;
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
    const ccomplex_t IT_0002 = conj(N_B1)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W1)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0013 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0001*IT_0005;
    const ccomplex_t IT_0017 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = IT_0000*IT_0006;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0022 = IT_0016*IT_0021;
    const ccomplex_t IT_0023 = IT_0019*IT_0021;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0018 + IT_0020 + 
      -IT_0022 + -IT_0023);
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = m_N_1*IT_0025;
    const ccomplex_t IT_0027 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0028 = IT_0016*IT_0027;
    const ccomplex_t IT_0029 = IT_0019*IT_0027;
    const ccomplex_t IT_0030 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0031 = IT_0016*IT_0030;
    const ccomplex_t IT_0032 = IT_0019*IT_0030;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0028 + IT_0029 + 
      -IT_0031 + -IT_0032);
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = m_N_2*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0026 + IT_0036;
    const ccomplex_t IT_0038 = IT_0015*IT_0037;
    const ccomplex_t IT_0039 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0040 = IT_0016*IT_0039;
    const ccomplex_t IT_0041 = IT_0019*IT_0039;
    const ccomplex_t IT_0042 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0043 = IT_0016*IT_0042;
    const ccomplex_t IT_0044 = IT_0019*IT_0042;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0040 + IT_0041 + 
      -IT_0043 + -IT_0044);
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = m_N_2*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B2)*e_em;
    const ccomplex_t IT_0049 = IT_0001*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W2)*e_em;
    const ccomplex_t IT_0052 = IT_0006*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + -IT_0053);
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0038 + (-0.25)*IT_0059;
    const ccomplex_t IT_0061 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0062 = IT_0016*IT_0061;
    const ccomplex_t IT_0063 = IT_0019*IT_0061;
    const ccomplex_t IT_0064 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0065 = IT_0016*IT_0064;
    const ccomplex_t IT_0066 = IT_0019*IT_0064;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0062 + IT_0063 + 
      -IT_0065 + -IT_0066);
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = m_N_3*IT_0068;
    const ccomplex_t IT_0070 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0071 = IT_0016*IT_0070;
    const ccomplex_t IT_0072 = IT_0019*IT_0070;
    const ccomplex_t IT_0073 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0074 = IT_0016*IT_0073;
    const ccomplex_t IT_0075 = IT_0019*IT_0073;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0071 + IT_0072 + 
      -IT_0074 + -IT_0075);
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = m_N_2*IT_0077;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = IT_0069 + IT_0079;
    const ccomplex_t IT_0081 = conj(N_B3)*e_em;
    const ccomplex_t IT_0082 = IT_0001*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W3)*e_em;
    const ccomplex_t IT_0085 = IT_0006*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + -IT_0086);
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0080*IT_0091;
    const ccomplex_t IT_0093 = (-0.5)*IT_0045;
    const ccomplex_t IT_0094 = m_N_2*IT_0093;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0047 + IT_0095;
    const ccomplex_t IT_0097 = IT_0058*IT_0096;
    const ccomplex_t IT_0098 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0099 = IT_0016*IT_0098;
    const ccomplex_t IT_0100 = IT_0019*IT_0098;
    const ccomplex_t IT_0101 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0102 = IT_0016*IT_0101;
    const ccomplex_t IT_0103 = IT_0019*IT_0101;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0099 + IT_0100 + 
      -IT_0102 + -IT_0103);
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = m_N_4*IT_0105;
    const ccomplex_t IT_0107 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0108 = IT_0016*IT_0107;
    const ccomplex_t IT_0109 = IT_0019*IT_0107;
    const ccomplex_t IT_0110 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0111 = IT_0016*IT_0110;
    const ccomplex_t IT_0112 = IT_0019*IT_0110;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0108 + IT_0109 + 
      -IT_0111 + -IT_0112);
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = m_N_2*IT_0114;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = IT_0106 + IT_0116;
    const ccomplex_t IT_0118 = conj(N_B4)*e_em;
    const ccomplex_t IT_0119 = IT_0001*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = conj(N_W4)*e_em;
    const ccomplex_t IT_0122 = IT_0006*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + -IT_0123);
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = cpow(s_23 + (-0.5)*IT_0012 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0117*IT_0128;
    const ccomplex_t IT_0130 = IT_0058*IT_0094;
    const ccomplex_t IT_0131 = 0.5*IT_0092 + (-0.25)*IT_0097 + 0.5*IT_0129 +
       0.25*IT_0130;
    const ccomplex_t IT_0132 = pow(m_Z, -2);
    const ccomplex_t IT_0133 = pow(m_Z, 2);
    const ccomplex_t IT_0134 = s_13*IT_0133;
    const ccomplex_t IT_0135 = s_14*s_34;
    const ccomplex_t IT_0136 = (-2)*IT_0135;
    const ccomplex_t IT_0137 = IT_0134 + IT_0136;
    const ccomplex_t IT_0138 = IT_0132*IT_0137;
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = s_13 + IT_0139;
    const ccomplex_t IT_0141 = conj(IT_0060) + conj(IT_0131);
    const ccomplex_t IT_0142 = pow(m_Z, 4);
    const ccomplex_t IT_0143 = s_13*IT_0142;
    const ccomplex_t IT_0144 = IT_0132*IT_0143;
    const ccomplex_t IT_0145 = (-0.5)*IT_0144;
    const ccomplex_t IT_0146 = IT_0134 + IT_0136 + IT_0145;
    const ccomplex_t IT_0147 = IT_0015*IT_0034;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0077*IT_0091;
    const ccomplex_t IT_0150 = IT_0114*IT_0128;
    const ccomplex_t IT_0151 = IT_0058*IT_0093;
    const ccomplex_t IT_0152 = IT_0046*IT_0058;
    const ccomplex_t IT_0153 = (-0.5)*IT_0149 + (-0.5)*IT_0150 + 0.25*IT_0151 
      + (-0.25)*IT_0152;
    const ccomplex_t IT_0154 = IT_0148 + -IT_0153;
    const ccomplex_t IT_0155 = -conj(IT_0153);
    const ccomplex_t IT_0156 = conj(IT_0148) + IT_0155;
    const ccomplex_t IT_0157 = e_em*IT_0016;
    const ccomplex_t IT_0158 = e_em*IT_0019;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0157 + IT_0158);
    const ccomplex_t IT_0160 = 0.5*IT_0159;
    const ccomplex_t IT_0161 = IT_0055*IT_0160;
    const ccomplex_t IT_0162 = pow(m_N_2, 2);
    const ccomplex_t IT_0163 = cpow(s_12 + 0.5*IT_0012 + 0.5*IT_0162 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0164 = IT_0161*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = (-0.25)*IT_0144;
    const ccomplex_t IT_0168 = IT_0135 + IT_0167;
    const ccomplex_t IT_0169 = s_34*IT_0133;
    const ccomplex_t IT_0170 = IT_0132*IT_0169;
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = s_34 + IT_0171;
    const ccomplex_t IT_0173 = m_N_2*IT_0172;
    const ccomplex_t IT_0174 = 2*IT_0173;
    const ccomplex_t IT_0175 = 2*conj(IT_0060);
    const ccomplex_t IT_0176 = 2*conj(IT_0131);
    const ccomplex_t IT_0177 = IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = s_14*s_34*IT_0132;
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = s_13 + IT_0179;
    const ccomplex_t IT_0181 = s_13*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0006;
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = cpow(s_13 + (-0.5)*IT_0162 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0185 = IT_0055*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = IT_0183*IT_0186;
    const ccomplex_t IT_0188 = -IT_0165 + (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = 2*IT_0147;
    const ccomplex_t IT_0190 = IT_0149 + IT_0150 + (-0.5)*IT_0151 + 0.5
      *IT_0152 + 0.5*IT_0187;
    const ccomplex_t IT_0191 = (-2)*IT_0190;
    const ccomplex_t IT_0192 = IT_0189 + IT_0191;
    const ccomplex_t IT_0193 = 2*conj(IT_0147);
    const ccomplex_t IT_0194 = (-2)*conj(IT_0190);
    const ccomplex_t IT_0195 = IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = m_N_2*IT_0180;
    const ccomplex_t IT_0197 = pow(s_34, 2);
    const ccomplex_t IT_0198 = s_13*IT_0132*IT_0197;
    const ccomplex_t IT_0199 = 2*IT_0188;
    const ccomplex_t IT_0200 = s_13*s_34;
    const ccomplex_t IT_0201 = 2*conj(IT_0188);
    const ccomplex_t IT_0202 = m_N_2*IT_0132*IT_0197;
    const ccomplex_t IT_0203 = (-0.25)*IT_0170;
    const ccomplex_t IT_0204 = s_34 + IT_0203;
    const ccomplex_t IT_0205 = m_N_2*IT_0204;
    const ccomplex_t IT_0206 = (-4)*IT_0205;
    const ccomplex_t IT_0207 = 2*IT_0131;
    const ccomplex_t IT_0208 = pow(s_14, 2);
    const ccomplex_t IT_0209 = IT_0132*IT_0208;
    const ccomplex_t IT_0210 = -IT_0162;
    const ccomplex_t IT_0211 = IT_0209 + IT_0210;
    const ccomplex_t IT_0212 = s_13*IT_0211;
    const ccomplex_t IT_0213 = 2*IT_0190;
    const ccomplex_t IT_0214 = (-0.5)*IT_0213;
    const ccomplex_t IT_0215 = IT_0147 + IT_0214;
    const ccomplex_t IT_0216 = (-2)*conj(IT_0147);
    const ccomplex_t IT_0217 = s_13*s_14;
    const ccomplex_t IT_0218 = s_34*IT_0162;
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = IT_0217 + IT_0219;
    const ccomplex_t IT_0221 = 4*(IT_0060 + IT_0131)*IT_0140*IT_0141 + (-4)
      *IT_0146*IT_0154*IT_0156 + 8*IT_0166*(IT_0156*IT_0168 + 0.125*IT_0174
      *IT_0177) + IT_0181*(conj(IT_0188)*IT_0192 + IT_0188*IT_0195) + (IT_0141
      *IT_0192 + (IT_0060 + IT_0131)*IT_0195)*IT_0196 + conj(IT_0188)*IT_0198
      *IT_0199 + (-2)*IT_0200*(conj(IT_0153)*IT_0188 + IT_0153*conj(IT_0188) + (
      -0.5)*(conj(IT_0148) + conj(IT_0166))*IT_0199 + (-0.5)*(IT_0148 + IT_0166)
      *IT_0201) + (IT_0141*IT_0199 + (IT_0060 + IT_0131)*IT_0201)*IT_0202 + (-2)
      *IT_0206*(IT_0141*IT_0148 + IT_0131*conj(IT_0148) + IT_0060*IT_0156 + (
      -0.5)*IT_0153*IT_0177 + (-0.5)*conj(IT_0153)*IT_0207) + (-4)*conj(IT_0166)
      *(IT_0146*IT_0166 + (-2)*IT_0154*IT_0168 + (-0.5)*IT_0174*(IT_0060 + 0.5
      *IT_0207)) + (-2)*IT_0212*((-0.5)*conj(IT_0147)*IT_0189 + conj(IT_0190)
      *IT_0215 + (-0.5)*IT_0190*IT_0216) + ((conj(IT_0148) + conj(IT_0166))
      *IT_0192 + (IT_0148 + IT_0166)*IT_0195 + 2*IT_0155*IT_0215 + 2*IT_0153*
      (conj(IT_0190) + 0.5*IT_0216))*IT_0220;
    return create_ccomplex_return(IT_0221);
}

