#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_4_to_Hp_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_4_to_Hp_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0005 + IT_0008 + 
      -IT_0010 + -IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 2*IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0014;
    const ccomplex_t IT_0020 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0024 = IT_0015*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = IT_0022 + IT_0025;
    const ccomplex_t IT_0027 = IT_0018 + IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = pow(m_Hp, 2);
    const ccomplex_t IT_0031 = pow(m_C_1, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0030 + IT_0031 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0029*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0013*IT_0034;
    const ccomplex_t IT_0036 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = IT_0007*IT_0036;
    const ccomplex_t IT_0039 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0040 = IT_0003*IT_0039;
    const ccomplex_t IT_0041 = IT_0007*IT_0039;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0037 + IT_0038 + 
      -IT_0040 + -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0045 = IT_0015*IT_0044;
    const ccomplex_t IT_0046 = 2*IT_0045;
    const ccomplex_t IT_0047 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0048 = IT_0019*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0051 = IT_0015*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = IT_0049 + IT_0052;
    const ccomplex_t IT_0054 = IT_0046 + IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_13 + IT_0030 + IT_0031 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0043*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0042;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0064 = IT_0003*IT_0063;
    const ccomplex_t IT_0065 = IT_0007*IT_0063;
    const ccomplex_t IT_0066 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0067 = IT_0003*IT_0066;
    const ccomplex_t IT_0068 = IT_0007*IT_0066;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0064 + IT_0065 + 
      -IT_0067 + -IT_0068);
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0030 + IT_0031 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0072 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0073 = IT_0019*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0076 = IT_0015*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = IT_0074 + IT_0077;
    const ccomplex_t IT_0079 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0080 = IT_0015*IT_0079;
    const ccomplex_t IT_0081 = 2*IT_0080;
    const ccomplex_t IT_0082 = IT_0078 + IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0071*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0070*IT_0086;
    const ccomplex_t IT_0088 = IT_0035 + (-0.5)*IT_0060 + 0.5*IT_0062 + IT_0087;
    const ccomplex_t IT_0089 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0090 = IT_0003*IT_0089;
    const ccomplex_t IT_0091 = IT_0007*IT_0089;
    const ccomplex_t IT_0092 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0093 = IT_0003*IT_0092;
    const ccomplex_t IT_0094 = IT_0007*IT_0092;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0090 + IT_0091 + 
      -IT_0093 + -IT_0094);
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_13 + IT_0030 + IT_0031 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0098 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0099 = IT_0015*IT_0098;
    const ccomplex_t IT_0100 = 2*IT_0099;
    const ccomplex_t IT_0101 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0102 = IT_0019*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0105 = IT_0015*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = IT_0103 + IT_0106;
    const ccomplex_t IT_0108 = IT_0100 + IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0097*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0096*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0006);
    const ccomplex_t IT_0115 = pow(m_N_4, 2);
    const ccomplex_t IT_0116 = cpow((-2)*s_12 + -IT_0031 + -IT_0115 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0056*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0114*IT_0118;
    const ccomplex_t IT_0120 = (-2)*IT_0035 + IT_0060 + -IT_0062 + (-2)
      *IT_0087 + (-2)*IT_0113 + IT_0119;
    const ccomplex_t IT_0121 = s_12*s_24;
    const ccomplex_t IT_0122 = s_14*IT_0115;
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = IT_0121 + IT_0123;
    const ccomplex_t IT_0125 = 2*IT_0124;
    const ccomplex_t IT_0126 = conj(IT_0120)*IT_0125;
    const ccomplex_t IT_0127 = m_N_4*IT_0043;
    const ccomplex_t IT_0128 = sin(beta);
    const ccomplex_t IT_0129 = IT_0006*IT_0128;
    const ccomplex_t IT_0130 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0131 = IT_0129*IT_0130;
    const ccomplex_t IT_0132 = (-2)*IT_0131;
    const ccomplex_t IT_0133 = IT_0001*IT_0128;
    const ccomplex_t IT_0134 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0135 = IT_0133*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0138 = IT_0129*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = IT_0136 + IT_0139;
    const ccomplex_t IT_0141 = IT_0132 + IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = IT_0057*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = IT_0127*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = pow(m_Z, -2);
    const ccomplex_t IT_0149 = pow(m_Z, 2);
    const ccomplex_t IT_0150 = s_24*IT_0149;
    const ccomplex_t IT_0151 = IT_0148*IT_0150;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = s_24 + IT_0152;
    const ccomplex_t IT_0154 = m_C_1*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0154;
    const ccomplex_t IT_0156 = 2*IT_0155;
    const ccomplex_t IT_0157 = conj(IT_0147)*IT_0156;
    const ccomplex_t IT_0158 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0159 = IT_0003*IT_0158;
    const ccomplex_t IT_0160 = IT_0007*IT_0158;
    const ccomplex_t IT_0161 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0162 = IT_0007*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0159 + -IT_0160 + (-2)
      *IT_0162);
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = IT_0056*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0167 = cpow((-2)*s_23 + IT_0030 + IT_0115 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0168 = IT_0165*IT_0166*IT_0167;
    const ccomplex_t IT_0169 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0170 = IT_0003*IT_0169;
    const ccomplex_t IT_0171 = IT_0007*IT_0169;
    const ccomplex_t IT_0172 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0173 = IT_0007*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0170 + -IT_0171 + (-2)
      *IT_0173);
    const ccomplex_t IT_0175 = (-0.5)*IT_0174;
    const ccomplex_t IT_0176 = IT_0056*IT_0175;
    const ccomplex_t IT_0177 = IT_0167*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = m_C_1*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0181 = cpow((-2)*s_23 + IT_0030 + IT_0115 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0182 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0183 = IT_0003*IT_0182;
    const ccomplex_t IT_0184 = IT_0007*IT_0182;
    const ccomplex_t IT_0185 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0186 = IT_0007*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0183 + -IT_0184 + (-2)
      *IT_0186);
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0190 = IT_0015*IT_0189;
    const ccomplex_t IT_0191 = 2*IT_0190;
    const ccomplex_t IT_0192 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0193 = IT_0019*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0196 = IT_0015*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = IT_0194 + IT_0197;
    const ccomplex_t IT_0199 = IT_0191 + IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = 0.5*IT_0200;
    const ccomplex_t IT_0202 = IT_0188*IT_0201;
    const ccomplex_t IT_0203 = IT_0180*IT_0181*IT_0202;
    const ccomplex_t IT_0204 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0205 = IT_0003*IT_0204;
    const ccomplex_t IT_0206 = IT_0007*IT_0204;
    const ccomplex_t IT_0207 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0208 = IT_0007*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*(IT_0205 + -IT_0206 + (-2)
      *IT_0208);
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0201*IT_0210;
    const ccomplex_t IT_0212 = IT_0181*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = m_C_1*IT_0213;
    const ccomplex_t IT_0215 = m_N_4*IT_0061;
    const ccomplex_t IT_0216 = IT_0145*IT_0215;
    const ccomplex_t IT_0217 = -IT_0127;
    const ccomplex_t IT_0218 = IT_0215 + IT_0217;
    const ccomplex_t IT_0219 = IT_0145*IT_0218;
    const ccomplex_t IT_0220 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0221 = IT_0133*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0224 = IT_0129*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = IT_0222 + IT_0225;
    const ccomplex_t IT_0227 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0228 = IT_0129*IT_0227;
    const ccomplex_t IT_0229 = (-2)*IT_0228;
    const ccomplex_t IT_0230 = IT_0226 + IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = IT_0097*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = m_N_1*IT_0096;
    const ccomplex_t IT_0236 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0237 = IT_0003*IT_0236;
    const ccomplex_t IT_0238 = IT_0007*IT_0236;
    const ccomplex_t IT_0239 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0240 = IT_0003*IT_0239;
    const ccomplex_t IT_0241 = IT_0007*IT_0239;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*(IT_0237 + IT_0238 + 
      -IT_0240 + -IT_0241);
    const ccomplex_t IT_0243 = (-0.5)*IT_0242;
    const ccomplex_t IT_0244 = m_N_4*IT_0243;
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = IT_0235 + IT_0245;
    const ccomplex_t IT_0247 = IT_0234*IT_0246;
    const ccomplex_t IT_0248 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0249 = IT_0129*IT_0248;
    const ccomplex_t IT_0250 = (-2)*IT_0249;
    const ccomplex_t IT_0251 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0252 = IT_0133*IT_0251;
    const ccomplex_t IT_0253 = 1.4142135623731*IT_0252;
    const ccomplex_t IT_0254 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0255 = IT_0129*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = IT_0253 + IT_0256;
    const ccomplex_t IT_0258 = IT_0250 + IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = 0.5*IT_0259;
    const ccomplex_t IT_0261 = IT_0071*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = m_N_2*IT_0070;
    const ccomplex_t IT_0264 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0265 = IT_0003*IT_0264;
    const ccomplex_t IT_0266 = IT_0007*IT_0264;
    const ccomplex_t IT_0267 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0268 = IT_0003*IT_0267;
    const ccomplex_t IT_0269 = IT_0007*IT_0267;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*(IT_0265 + IT_0266 + 
      -IT_0268 + -IT_0269);
    const ccomplex_t IT_0271 = (-0.5)*IT_0270;
    const ccomplex_t IT_0272 = m_N_4*IT_0271;
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = IT_0263 + IT_0273;
    const ccomplex_t IT_0275 = IT_0262*IT_0274;
    const ccomplex_t IT_0276 = m_N_3*IT_0013;
    const ccomplex_t IT_0277 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0278 = IT_0003*IT_0277;
    const ccomplex_t IT_0279 = IT_0007*IT_0277;
    const ccomplex_t IT_0280 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0281 = IT_0003*IT_0280;
    const ccomplex_t IT_0282 = IT_0007*IT_0280;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*(IT_0278 + IT_0279 + 
      -IT_0281 + -IT_0282);
    const ccomplex_t IT_0284 = (-0.5)*IT_0283;
    const ccomplex_t IT_0285 = m_N_4*IT_0284;
    const ccomplex_t IT_0286 = -IT_0285;
    const ccomplex_t IT_0287 = IT_0276 + IT_0286;
    const ccomplex_t IT_0288 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0289 = IT_0129*IT_0288;
    const ccomplex_t IT_0290 = 1.4142135623731*IT_0289;
    const ccomplex_t IT_0291 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0292 = IT_0133*IT_0291;
    const ccomplex_t IT_0293 = 1.4142135623731*IT_0292;
    const ccomplex_t IT_0294 = IT_0290 + IT_0293;
    const ccomplex_t IT_0295 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0296 = IT_0129*IT_0295;
    const ccomplex_t IT_0297 = (-2)*IT_0296;
    const ccomplex_t IT_0298 = IT_0294 + IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*IT_0298;
    const ccomplex_t IT_0300 = 0.5*IT_0299;
    const ccomplex_t IT_0301 = IT_0032*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*IT_0301;
    const ccomplex_t IT_0303 = IT_0287*IT_0302;
    const ccomplex_t IT_0304 = -IT_0168 + -IT_0179 + -IT_0203 + -IT_0214 + (
      -0.5)*IT_0216 + (-0.5)*IT_0219 + -IT_0247 + -IT_0275 + -IT_0303;
    const ccomplex_t IT_0305 = IT_0156*conj(IT_0304);
    const ccomplex_t IT_0306 = IT_0059*IT_0215;
    const ccomplex_t IT_0307 = (-0.5)*IT_0306;
    const ccomplex_t IT_0308 = s_14*IT_0149;
    const ccomplex_t IT_0309 = IT_0148*IT_0308;
    const ccomplex_t IT_0310 = (-0.25)*IT_0309;
    const ccomplex_t IT_0311 = s_14 + IT_0310;
    const ccomplex_t IT_0312 = m_N_4*IT_0311;
    const ccomplex_t IT_0313 = (-4)*IT_0312;
    const ccomplex_t IT_0314 = (-2)*IT_0313;
    const ccomplex_t IT_0315 = conj(IT_0307)*IT_0314;
    const ccomplex_t IT_0316 = IT_0143*IT_0175;
    const ccomplex_t IT_0317 = IT_0166*IT_0167*IT_0316;
    const ccomplex_t IT_0318 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0319 = IT_0129*IT_0318;
    const ccomplex_t IT_0320 = (-2)*IT_0319;
    const ccomplex_t IT_0321 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0322 = IT_0133*IT_0321;
    const ccomplex_t IT_0323 = 1.4142135623731*IT_0322;
    const ccomplex_t IT_0324 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0325 = IT_0129*IT_0324;
    const ccomplex_t IT_0326 = 1.4142135623731*IT_0325;
    const ccomplex_t IT_0327 = IT_0323 + IT_0326;
    const ccomplex_t IT_0328 = IT_0320 + IT_0327;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*IT_0328;
    const ccomplex_t IT_0330 = 0.5*IT_0329;
    const ccomplex_t IT_0331 = IT_0210*IT_0330;
    const ccomplex_t IT_0332 = IT_0180*IT_0181*IT_0331;
    const ccomplex_t IT_0333 = IT_0143*IT_0164;
    const ccomplex_t IT_0334 = IT_0167*IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = m_C_1*IT_0335;
    const ccomplex_t IT_0337 = IT_0059*IT_0127;
    const ccomplex_t IT_0338 = IT_0059*IT_0218;
    const ccomplex_t IT_0339 = m_N_1*IT_0243;
    const ccomplex_t IT_0340 = m_N_4*IT_0096;
    const ccomplex_t IT_0341 = -IT_0340;
    const ccomplex_t IT_0342 = IT_0339 + IT_0341;
    const ccomplex_t IT_0343 = IT_0112*IT_0342;
    const ccomplex_t IT_0344 = m_N_2*IT_0271;
    const ccomplex_t IT_0345 = m_N_4*IT_0070;
    const ccomplex_t IT_0346 = -IT_0345;
    const ccomplex_t IT_0347 = IT_0344 + IT_0346;
    const ccomplex_t IT_0348 = IT_0086*IT_0347;
    const ccomplex_t IT_0349 = m_N_3*IT_0284;
    const ccomplex_t IT_0350 = m_N_4*IT_0013;
    const ccomplex_t IT_0351 = -IT_0350;
    const ccomplex_t IT_0352 = IT_0349 + IT_0351;
    const ccomplex_t IT_0353 = IT_0034*IT_0352;
    const ccomplex_t IT_0354 = IT_0188*IT_0330;
    const ccomplex_t IT_0355 = IT_0181*IT_0354;
    const ccomplex_t IT_0356 = (0 + _Complex_I*1)*IT_0355;
    const ccomplex_t IT_0357 = m_C_1*IT_0356;
    const ccomplex_t IT_0358 = -IT_0317 + -IT_0332 + -IT_0336 + 0.5*IT_0337 + 
      (-0.5)*IT_0338 + IT_0343 + IT_0348 + IT_0353 + -IT_0357;
    const ccomplex_t IT_0359 = IT_0314*conj(IT_0358);
    const ccomplex_t IT_0360 = -IT_0335 + -IT_0356;
    const ccomplex_t IT_0361 = m_C_1*m_N_4;
    const ccomplex_t IT_0362 = IT_0149*IT_0361;
    const ccomplex_t IT_0363 = 3*IT_0362;
    const ccomplex_t IT_0364 = (-2)*IT_0363;
    const ccomplex_t IT_0365 = conj(IT_0360)*IT_0364;
    const ccomplex_t IT_0366 = IT_0234*IT_0243;
    const ccomplex_t IT_0367 = IT_0262*IT_0271;
    const ccomplex_t IT_0368 = IT_0061*IT_0145;
    const ccomplex_t IT_0369 = IT_0043*IT_0145;
    const ccomplex_t IT_0370 = IT_0284*IT_0302;
    const ccomplex_t IT_0371 = -IT_0366 + -IT_0367 + 0.5*IT_0368 + (-0.5)
      *IT_0369 + -IT_0370;
    const ccomplex_t IT_0372 = (-3)*IT_0362;
    const ccomplex_t IT_0373 = (-2)*IT_0372;
    const ccomplex_t IT_0374 = conj(IT_0371)*IT_0373;
    const ccomplex_t IT_0375 = IT_0119 + 2*IT_0178 + 2*IT_0213;
    const ccomplex_t IT_0376 = s_12*s_14;
    const ccomplex_t IT_0377 = s_24*IT_0031;
    const ccomplex_t IT_0378 = -IT_0377;
    const ccomplex_t IT_0379 = IT_0376 + IT_0378;
    const ccomplex_t IT_0380 = (-2)*IT_0379;
    const ccomplex_t IT_0381 = conj(IT_0375)*IT_0380;
    const ccomplex_t IT_0382 = conj(IT_0088) + conj(IT_0113);
    const ccomplex_t IT_0383 = s_12*IT_0149;
    const ccomplex_t IT_0384 = pow(m_Z, 4);
    const ccomplex_t IT_0385 = s_12*IT_0384;
    const ccomplex_t IT_0386 = IT_0148*IT_0385;
    const ccomplex_t IT_0387 = (-0.5)*IT_0386;
    const ccomplex_t IT_0388 = s_14*s_24;
    const ccomplex_t IT_0389 = (-2)*IT_0388;
    const ccomplex_t IT_0390 = IT_0383 + IT_0387 + IT_0389;
    const ccomplex_t IT_0391 = (-4)*IT_0390;
    const ccomplex_t IT_0392 = IT_0382*IT_0391;
    const ccomplex_t IT_0393 = -IT_0178 + -IT_0213;
    const ccomplex_t IT_0394 = (-0.25)*IT_0386;
    const ccomplex_t IT_0395 = IT_0388 + IT_0394;
    const ccomplex_t IT_0396 = (-8)*IT_0395;
    const ccomplex_t IT_0397 = conj(IT_0393)*IT_0396;
    const ccomplex_t IT_0398 = IT_0126 + IT_0157 + IT_0305 + IT_0315 + IT_0359
       + IT_0365 + IT_0374 + IT_0381 + IT_0392 + IT_0397;
    const ccomplex_t IT_0399 = conj(IT_0147) + conj(IT_0304);
    const ccomplex_t IT_0400 = IT_0383 + IT_0389;
    const ccomplex_t IT_0401 = IT_0148*IT_0400;
    const ccomplex_t IT_0402 = (-0.5)*IT_0401;
    const ccomplex_t IT_0403 = s_12 + IT_0402;
    const ccomplex_t IT_0404 = 4*IT_0403;
    const ccomplex_t IT_0405 = IT_0147*IT_0404;
    const ccomplex_t IT_0406 = conj(IT_0307) + conj(IT_0358);
    const ccomplex_t IT_0407 = IT_0116*IT_0143;
    const ccomplex_t IT_0408 = (0 + _Complex_I*1)*IT_0407;
    const ccomplex_t IT_0409 = IT_0114*IT_0408;
    const ccomplex_t IT_0410 = 2*IT_0366 + 2*IT_0367 + -IT_0368 + IT_0369 + 2
      *IT_0370 + -IT_0409;
    const ccomplex_t IT_0411 = pow(s_24, 2);
    const ccomplex_t IT_0412 = IT_0148*IT_0411;
    const ccomplex_t IT_0413 = -IT_0115;
    const ccomplex_t IT_0414 = IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = m_C_1*IT_0414;
    const ccomplex_t IT_0416 = 2*IT_0415;
    const ccomplex_t IT_0417 = IT_0410*IT_0416;
    const ccomplex_t IT_0418 = 2*IT_0335 + 2*IT_0356 + -IT_0409;
    const ccomplex_t IT_0419 = s_14*s_24*IT_0148;
    const ccomplex_t IT_0420 = -IT_0419;
    const ccomplex_t IT_0421 = s_12 + IT_0420;
    const ccomplex_t IT_0422 = m_C_1*IT_0421;
    const ccomplex_t IT_0423 = (-2)*IT_0422;
    const ccomplex_t IT_0424 = IT_0418*IT_0423;
    const ccomplex_t IT_0425 = pow(s_14, 2);
    const ccomplex_t IT_0426 = IT_0148*IT_0425;
    const ccomplex_t IT_0427 = -IT_0031;
    const ccomplex_t IT_0428 = IT_0426 + IT_0427;
    const ccomplex_t IT_0429 = s_12*IT_0428;
    const ccomplex_t IT_0430 = 2*IT_0379;
    const ccomplex_t IT_0431 = IT_0361*IT_0421;
    const ccomplex_t IT_0432 = 2*IT_0431;
    const ccomplex_t IT_0433 = s_12*IT_0421;
    const ccomplex_t IT_0434 = (-2)*IT_0433;
    const ccomplex_t IT_0435 = IT_0361*IT_0428;
    const ccomplex_t IT_0436 = (-2)*IT_0435;
    const ccomplex_t IT_0437 = m_N_4*IT_0428;
    const ccomplex_t IT_0438 = (-2)*IT_0437;
    const ccomplex_t IT_0439 = 0.5*IT_0380*IT_0382 + 0.5*IT_0399*IT_0423 +
       conj(IT_0375)*IT_0429 + 0.5*conj(IT_0393)*IT_0430 + 0.5*conj(IT_0410)
      *IT_0432 + 0.5*conj(IT_0120)*IT_0434 + 0.5*conj(IT_0418)*IT_0436 + 0.5
      *IT_0406*IT_0438;
    const ccomplex_t IT_0440 = 2*IT_0418*(0.5*conj(IT_0371)*IT_0380 + conj
      (IT_0418)*IT_0429 + 0.5*conj(IT_0360)*IT_0430 + 0.5*conj(IT_0120)*IT_0432 
      + 0.5*conj(IT_0410)*IT_0434 + 0.5*conj(IT_0375)*IT_0436 + 0.5*IT_0399
      *IT_0438);
    const ccomplex_t IT_0441 = s_12*IT_0414;
    const ccomplex_t IT_0442 = 2*IT_0441;
    const ccomplex_t IT_0443 = m_N_4*IT_0421;
    const ccomplex_t IT_0444 = 2*IT_0443;
    const ccomplex_t IT_0445 = (-2)*IT_0124;
    const ccomplex_t IT_0446 = IT_0361*IT_0414;
    const ccomplex_t IT_0447 = (-2)*IT_0446;
    const ccomplex_t IT_0448 = IT_0156*IT_0304;
    const ccomplex_t IT_0449 = IT_0314*IT_0358;
    const ccomplex_t IT_0450 = (-0.25)*IT_0151;
    const ccomplex_t IT_0451 = s_24 + IT_0450;
    const ccomplex_t IT_0452 = m_C_1*IT_0451;
    const ccomplex_t IT_0453 = (-4)*IT_0452;
    const ccomplex_t IT_0454 = 2*IT_0453;
    const ccomplex_t IT_0455 = 0.5*IT_0309;
    const ccomplex_t IT_0456 = s_14 + IT_0455;
    const ccomplex_t IT_0457 = m_N_4*IT_0456;
    const ccomplex_t IT_0458 = 2*IT_0457;
    const ccomplex_t IT_0459 = (-2)*IT_0458;
    const ccomplex_t IT_0460 = conj(IT_0410)*IT_0416;
    const ccomplex_t IT_0461 = conj(IT_0120)*IT_0444;
    const ccomplex_t IT_0462 = 6*IT_0361;
    const ccomplex_t IT_0463 = conj(IT_0147)*IT_0462;
    const ccomplex_t IT_0464 = conj(IT_0418)*IT_0423;
    const ccomplex_t IT_0465 = conj(IT_0375)*IT_0438;
    const ccomplex_t IT_0466 = IT_0404*IT_0406;
    const ccomplex_t IT_0467 = conj(IT_0304)*IT_0462;
    const ccomplex_t IT_0468 = conj(IT_0393)*IT_0459;
    const ccomplex_t IT_0469 = IT_0156*conj(IT_0371);
    const ccomplex_t IT_0470 = conj(IT_0120)*IT_0416;
    const ccomplex_t IT_0471 = conj(IT_0410)*IT_0444;
    const ccomplex_t IT_0472 = IT_0470 + IT_0471;
    const ccomplex_t IT_0473 = conj(IT_0375)*IT_0423;
    const ccomplex_t IT_0474 = conj(IT_0418)*IT_0438;
    const ccomplex_t IT_0475 = IT_0406*IT_0462;
    const ccomplex_t IT_0476 = IT_0314*conj(IT_0371);
    const ccomplex_t IT_0477 = conj(IT_0360)*IT_0459;
    const ccomplex_t IT_0478 = (IT_0088 + IT_0113)*IT_0398 + IT_0399*IT_0405 +
       IT_0371*(IT_0373*IT_0382 + conj(IT_0371)*IT_0391 + IT_0364*conj(IT_0393) 
      + conj(IT_0360)*IT_0396 + IT_0314*IT_0399 + IT_0156*IT_0406 + IT_0125*conj
      (IT_0410) + IT_0380*conj(IT_0418)) + IT_0406*(IT_0417 + IT_0424) + 2
      *IT_0375*IT_0439 + IT_0440 + IT_0410*(IT_0125*conj(IT_0371) + conj(IT_0375
      )*IT_0432 + conj(IT_0418)*IT_0434 + conj(IT_0410)*IT_0442 + IT_0399
      *IT_0444 + conj(IT_0360)*IT_0445 + conj(IT_0120)*IT_0447) + IT_0120*
      (IT_0125*IT_0382 + IT_0399*IT_0416 + conj(IT_0418)*IT_0432 + conj(IT_0375)
      *IT_0434 + conj(IT_0120)*IT_0442 + IT_0406*IT_0444 + conj(IT_0393)*IT_0445
       + conj(IT_0410)*IT_0447) + IT_0382*(IT_0448 + IT_0449) + IT_0360*(IT_0364
      *IT_0382 + conj(IT_0360)*IT_0391 + IT_0373*conj(IT_0393) + conj(IT_0371)
      *IT_0396 + conj(IT_0418)*IT_0430 + conj(IT_0410)*IT_0445 + IT_0406*IT_0454
       + IT_0399*IT_0459) + IT_0393*(IT_0364*conj(IT_0371) + conj(IT_0360)
      *IT_0373 + IT_0391*conj(IT_0393) + IT_0382*IT_0396 + conj(IT_0375)*IT_0430
       + conj(IT_0120)*IT_0445 + IT_0399*IT_0454 + IT_0406*IT_0459) + IT_0358*
      (conj(IT_0360)*IT_0454 + IT_0460 + IT_0461 + IT_0463 + IT_0464 + IT_0465 +
       IT_0466 + IT_0467 + IT_0468 + IT_0469) + IT_0307*(IT_0314*IT_0382 + conj
      (IT_0360)*IT_0454 + IT_0460 + IT_0461 + IT_0463 + IT_0464 + IT_0465 +
       IT_0466 + IT_0467 + IT_0468 + IT_0469) + IT_0147*(IT_0156*IT_0382 + conj
      (IT_0393)*IT_0454 + IT_0472 + IT_0473 + IT_0474 + IT_0475 + IT_0476 +
       IT_0477) + IT_0304*(IT_0399*IT_0404 + conj(IT_0393)*IT_0454 + IT_0472 +
       IT_0473 + IT_0474 + IT_0475 + IT_0476 + IT_0477);
    return create_ccomplex_return(IT_0478);
}

