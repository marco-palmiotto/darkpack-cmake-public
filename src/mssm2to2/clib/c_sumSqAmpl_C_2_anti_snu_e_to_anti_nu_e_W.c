#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_snu_e_to_anti_nu_e_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_snu_e_to_anti_nu_e_W(
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
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
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
    const ccomplex_t IT_0000 = pow(m_W, -2);
    const ccomplex_t IT_0001 = pow(m_W, 2);
    const ccomplex_t IT_0002 = s_13*IT_0001;
    const ccomplex_t IT_0003 = s_14*s_34;
    const ccomplex_t IT_0004 = (-2)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = s_13 + IT_0007;
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = N_B3*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = N_W3*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0013 + -IT_0018);
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0022 = IT_0015*IT_0021;
    const ccomplex_t IT_0023 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0024 = IT_0015*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + (-0.5)*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0020*IT_0027;
    const ccomplex_t IT_0029 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0030 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = m_C_2*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = N_B1*e_em;
    const ccomplex_t IT_0036 = IT_0010*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = N_W1*e_em;
    const ccomplex_t IT_0039 = IT_0015*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0044 = IT_0015*IT_0043;
    const ccomplex_t IT_0045 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0046 = IT_0015*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + 0.5*IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = IT_0042*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0052 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0050*IT_0051*IT_0052;
    const ccomplex_t IT_0054 = N_B2*e_em;
    const ccomplex_t IT_0055 = IT_0010*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = N_W2*e_em;
    const ccomplex_t IT_0058 = IT_0015*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + -IT_0059);
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0063 = IT_0015*IT_0062;
    const ccomplex_t IT_0064 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0065 = IT_0015*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + 0.5*IT_0066);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0061*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0071 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0069*IT_0070*IT_0071;
    const ccomplex_t IT_0073 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0074 = IT_0015*IT_0073;
    const ccomplex_t IT_0075 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0076 = IT_0015*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + 0.5*IT_0077);
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = IT_0020*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0082 = IT_0030*IT_0080*IT_0081;
    const ccomplex_t IT_0083 = N_B4*e_em;
    const ccomplex_t IT_0084 = IT_0010*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = N_W4*e_em;
    const ccomplex_t IT_0087 = IT_0015*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0085 + -IT_0088);
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0092 = IT_0015*IT_0091;
    const ccomplex_t IT_0093 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0094 = IT_0015*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + 0.5*IT_0095);
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = IT_0090*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0100 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0101 = IT_0098*IT_0099*IT_0100;
    const ccomplex_t IT_0102 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0103 = IT_0015*IT_0102;
    const ccomplex_t IT_0104 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0105 = IT_0015*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0103 + (-0.5)*IT_0106);
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = IT_0042*IT_0108;
    const ccomplex_t IT_0110 = IT_0052*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = m_C_2*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0114 = IT_0015*IT_0113;
    const ccomplex_t IT_0115 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0116 = IT_0015*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0114 + (-0.5)*IT_0117);
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = IT_0090*IT_0119;
    const ccomplex_t IT_0121 = IT_0100*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = m_C_2*IT_0122;
    const ccomplex_t IT_0124 = pow(m_e, 2);
    const ccomplex_t IT_0125 = pow(m_C_2, 2);
    const ccomplex_t IT_0126 = cpow((-2)*s_12 + -IT_0029 + IT_0124 + -IT_0125 
      + -reg_prop, -1);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*m_e;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0015;
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = pow(m_W, -1);
    const ccomplex_t IT_0131 = cos(beta);
    const ccomplex_t IT_0132 = cpow(IT_0131, -1);
    const ccomplex_t IT_0133 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d2)
      *e_em*IT_0015*IT_0130*IT_0132;
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = IT_0129*IT_0134;
    const ccomplex_t IT_0136 = IT_0126*IT_0127*IT_0135;
    const ccomplex_t IT_0137 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0138 = IT_0015*IT_0137;
    const ccomplex_t IT_0139 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0140 = IT_0015*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + (-0.5)*IT_0141);
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = IT_0061*IT_0143;
    const ccomplex_t IT_0145 = IT_0071*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = m_C_2*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0053 + (-0.5)*IT_0072 + (-0.5)
      *IT_0082 + (-0.5)*IT_0101 + 0.5*IT_0112 + 0.5*IT_0123 + IT_0136 + 0.5
      *IT_0147;
    const ccomplex_t IT_0149 = IT_0034 + IT_0148;
    const ccomplex_t IT_0150 = conj(IT_0034) + conj(IT_0148);
    const ccomplex_t IT_0151 = pow(m_W, 4);
    const ccomplex_t IT_0152 = s_13*IT_0151;
    const ccomplex_t IT_0153 = IT_0000*IT_0152;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = IT_0002 + IT_0004 + IT_0154;
    const ccomplex_t IT_0156 = (-0.5)*IT_0111 + (-0.5)*IT_0122 + (-0.5)*IT_0146;
    const ccomplex_t IT_0157 = 0.5*IT_0032;
    const ccomplex_t IT_0158 = IT_0156 + -IT_0157;
    const ccomplex_t IT_0159 = -conj(IT_0157);
    const ccomplex_t IT_0160 = conj(IT_0156) + IT_0159;
    const ccomplex_t IT_0161 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = s_13 + IT_0162;
    const ccomplex_t IT_0164 = m_C_2*IT_0163;
    const ccomplex_t IT_0165 = IT_0111 + IT_0122 + IT_0146;
    const ccomplex_t IT_0166 = 2*IT_0165;
    const ccomplex_t IT_0167 = -IT_0032;
    const ccomplex_t IT_0168 = (-2)*IT_0167;
    const ccomplex_t IT_0169 = IT_0166 + IT_0168;
    const ccomplex_t IT_0170 = 2*conj(IT_0165);
    const ccomplex_t IT_0171 = (-2)*conj(IT_0167);
    const ccomplex_t IT_0172 = IT_0170 + IT_0171;
    const ccomplex_t IT_0173 = s_13*IT_0163;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0015;
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = IT_0129*IT_0175;
    const ccomplex_t IT_0177 = IT_0126*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = 2*IT_0178;
    const ccomplex_t IT_0180 = 2*IT_0167;
    const ccomplex_t IT_0181 = (-2)*IT_0165;
    const ccomplex_t IT_0182 = IT_0180 + IT_0181;
    const ccomplex_t IT_0183 = 2*conj(IT_0167);
    const ccomplex_t IT_0184 = (-2)*conj(IT_0165);
    const ccomplex_t IT_0185 = IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = s_13*s_14;
    const ccomplex_t IT_0187 = s_34*IT_0125;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = IT_0186 + IT_0188;
    const ccomplex_t IT_0190 = s_34*IT_0001;
    const ccomplex_t IT_0191 = IT_0000*IT_0190;
    const ccomplex_t IT_0192 = (-0.25)*IT_0191;
    const ccomplex_t IT_0193 = s_34 + IT_0192;
    const ccomplex_t IT_0194 = m_C_2*IT_0193;
    const ccomplex_t IT_0195 = (-4)*IT_0194;
    const ccomplex_t IT_0196 = (-2)*IT_0034;
    const ccomplex_t IT_0197 = (-2)*IT_0148;
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = (-2)*conj(IT_0034);
    const ccomplex_t IT_0200 = (-2)*conj(IT_0148);
    const ccomplex_t IT_0201 = IT_0199 + IT_0200;
    const ccomplex_t IT_0202 = pow(s_34, 2);
    const ccomplex_t IT_0203 = s_13*IT_0000*IT_0202;
    const ccomplex_t IT_0204 = 2*IT_0179;
    const ccomplex_t IT_0205 = (-0.25)*IT_0153;
    const ccomplex_t IT_0206 = IT_0003 + IT_0205;
    const ccomplex_t IT_0207 = 0.5*IT_0191;
    const ccomplex_t IT_0208 = s_34 + IT_0207;
    const ccomplex_t IT_0209 = m_C_2*IT_0208;
    const ccomplex_t IT_0210 = 2*IT_0209;
    const ccomplex_t IT_0211 = s_13*s_34;
    const ccomplex_t IT_0212 = 2*conj(IT_0179);
    const ccomplex_t IT_0213 = pow(s_14, 2);
    const ccomplex_t IT_0214 = IT_0000*IT_0213;
    const ccomplex_t IT_0215 = -IT_0125;
    const ccomplex_t IT_0216 = IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = s_13*IT_0216;
    const ccomplex_t IT_0218 = (-2)*IT_0179;
    const ccomplex_t IT_0219 = (-2)*conj(IT_0179);
    const ccomplex_t IT_0220 = m_C_2*IT_0000*IT_0202;
    const ccomplex_t IT_0221 = 4*IT_0008*IT_0149*IT_0150 + (-4)*IT_0155
      *IT_0158*IT_0160 + IT_0164*(IT_0150*IT_0169 + IT_0149*IT_0172) + IT_0173*
      (conj(IT_0179)*IT_0182 + IT_0179*IT_0185) + (conj(IT_0156)*IT_0169 +
       IT_0156*IT_0172 + (conj(IT_0157) + conj(IT_0178))*IT_0182 + (IT_0157 +
       IT_0178)*IT_0185)*IT_0189 + 2*IT_0195*(IT_0150*IT_0157 + IT_0149*conj
      (IT_0157) + 0.5*conj(IT_0156)*IT_0198 + 0.5*IT_0156*IT_0201) + conj
      (IT_0179)*IT_0203*IT_0204 + (-8)*conj(IT_0178)*(IT_0158*IT_0206 + (-0.125)
      *IT_0198*IT_0210 + (-0.125)*IT_0204*IT_0211) + (-4)*IT_0178*(IT_0155*conj
      (IT_0178) + 2*IT_0160*IT_0206 + (-0.25)*IT_0201*IT_0210 + (-0.25)*IT_0211
      *IT_0212) + (conj(IT_0165)*IT_0166 + conj(IT_0167)*IT_0182 + IT_0167
      *IT_0184)*IT_0217 + IT_0211*(conj(IT_0157)*IT_0204 + IT_0157*IT_0212 +
       conj(IT_0156)*IT_0218 + IT_0156*IT_0219) + (IT_0150*IT_0218 + IT_0149
      *IT_0219)*IT_0220;
    return create_ccomplex_return(IT_0221);
}

