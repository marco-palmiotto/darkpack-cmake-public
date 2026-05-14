#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_3_to_Z_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_3_to_Z_W(
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
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
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
    const ccomplex_t IT_0000 = pow(m_Z, -2);
    const ccomplex_t IT_0001 = pow(m_Z, 2);
    const ccomplex_t IT_0002 = s_23*IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -2);
    const ccomplex_t IT_0005 = s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0006 = s_24*IT_0001;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, 2);
    const ccomplex_t IT_0009 = s_23*IT_0008;
    const ccomplex_t IT_0010 = s_24*s_34;
    const ccomplex_t IT_0011 = (-2)*IT_0010;
    const ccomplex_t IT_0012 = IT_0009 + IT_0011;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0003 + (-0.5)*IT_0007 + (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = s_23 + IT_0014;
    const ccomplex_t IT_0016 = m_C_1*IT_0015;
    const ccomplex_t IT_0017 = 2*IT_0016;
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = sin(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*IT_0018*IT_0020;
    const ccomplex_t IT_0022 = m_N_3*IT_0021;
    const ccomplex_t IT_0023 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + 0.5*IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = pow(m_C_1, 2);
    const ccomplex_t IT_0031 = pow(m_N_3, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0030 + -IT_0031 + -reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0029*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0022*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = m_C_1*m_N_3;
    const ccomplex_t IT_0038 = s_34*IT_0001;
    const ccomplex_t IT_0039 = IT_0005*IT_0038;
    const ccomplex_t IT_0040 = IT_0001*IT_0008;
    const ccomplex_t IT_0041 = pow(s_34, 2);
    const ccomplex_t IT_0042 = (-2)*IT_0041;
    const ccomplex_t IT_0043 = IT_0040 + IT_0042;
    const ccomplex_t IT_0044 = IT_0004*IT_0043;
    const ccomplex_t IT_0045 = (-0.333333333333333)*IT_0039 + (
      -0.333333333333333)*IT_0044;
    const ccomplex_t IT_0046 = IT_0001 + IT_0045;
    const ccomplex_t IT_0047 = IT_0037*IT_0046;
    const ccomplex_t IT_0048 = 3*IT_0047;
    const ccomplex_t IT_0049 = 2*IT_0021;
    const ccomplex_t IT_0050 = IT_0034*IT_0049;
    const ccomplex_t IT_0051 = s_13*IT_0001;
    const ccomplex_t IT_0052 = IT_0000*IT_0051;
    const ccomplex_t IT_0053 = s_14*IT_0001;
    const ccomplex_t IT_0054 = IT_0005*IT_0053;
    const ccomplex_t IT_0055 = s_13*IT_0008;
    const ccomplex_t IT_0056 = s_14*s_34;
    const ccomplex_t IT_0057 = (-2)*IT_0056;
    const ccomplex_t IT_0058 = IT_0055 + IT_0057;
    const ccomplex_t IT_0059 = IT_0004*IT_0058;
    const ccomplex_t IT_0060 = 0.5*IT_0052 + (-0.5)*IT_0054 + (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = s_13 + IT_0060;
    const ccomplex_t IT_0062 = m_N_3*IT_0061;
    const ccomplex_t IT_0063 = 2*IT_0062;
    const ccomplex_t IT_0064 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0065 = IT_0020*IT_0064;
    const ccomplex_t IT_0066 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0067 = IT_0020*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + (-0.5)*IT_0068);
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = IT_0032*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0022*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_12 + IT_0008 + -IT_0030 + -IT_0031 
      + -reg_prop, -1);
    const ccomplex_t IT_0075 = sin(beta);
    const ccomplex_t IT_0076 = IT_0020*IT_0075;
    const ccomplex_t IT_0077 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = 2*IT_0078;
    const ccomplex_t IT_0080 = cpow(IT_0018, -1);
    const ccomplex_t IT_0081 = IT_0075*IT_0080;
    const ccomplex_t IT_0082 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0086 = IT_0076*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = IT_0079 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*m_W*e_em*IT_0019*IT_0080;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = IT_0091*IT_0093;
    const ccomplex_t IT_0095 = IT_0074*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = IT_0073 + -IT_0096;
    const ccomplex_t IT_0098 = IT_0004*IT_0041;
    const ccomplex_t IT_0099 = IT_0001*IT_0004*IT_0041;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = pow(m_Z, 4);
    const ccomplex_t IT_0102 = IT_0100 + IT_0101;
    const ccomplex_t IT_0103 = IT_0000*IT_0102;
    const ccomplex_t IT_0104 = -IT_0001 + 0.25*IT_0103;
    const ccomplex_t IT_0105 = IT_0098 + IT_0104;
    const ccomplex_t IT_0106 = IT_0037*IT_0105;
    const ccomplex_t IT_0107 = 4*IT_0106;
    const ccomplex_t IT_0108 = -IT_0021;
    const ccomplex_t IT_0109 = IT_0034*IT_0108;
    const ccomplex_t IT_0110 = m_C_1*IT_0004;
    const ccomplex_t IT_0111 = s_23*IT_0001*IT_0008;
    const ccomplex_t IT_0112 = IT_0000*IT_0111;
    const ccomplex_t IT_0113 = (-2)*IT_0010 + (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0009 + IT_0113;
    const ccomplex_t IT_0115 = IT_0110*IT_0114;
    const ccomplex_t IT_0116 = (-2)*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_23 + IT_0001 + IT_0031 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0118 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0119 = IT_0020*IT_0118;
    const ccomplex_t IT_0120 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0121 = IT_0020*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0119 + 0.5*IT_0122);
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0117*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = IT_0019*IT_0080;
    const ccomplex_t IT_0128 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0129 = IT_0127*IT_0128;
    const ccomplex_t IT_0130 = IT_0018*IT_0020;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0133 = IT_0127*IT_0132;
    const ccomplex_t IT_0134 = IT_0130*IT_0132;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0129 + IT_0131 + 
      -IT_0133 + -IT_0134);
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = m_N_1*IT_0136;
    const ccomplex_t IT_0138 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0139 = IT_0127*IT_0138;
    const ccomplex_t IT_0140 = IT_0130*IT_0138;
    const ccomplex_t IT_0141 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0142 = IT_0127*IT_0141;
    const ccomplex_t IT_0143 = IT_0130*IT_0141;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0139 + IT_0140 + 
      -IT_0142 + -IT_0143);
    const ccomplex_t IT_0145 = (-0.5)*IT_0144;
    const ccomplex_t IT_0146 = m_N_3*IT_0145;
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = IT_0137 + IT_0147;
    const ccomplex_t IT_0149 = IT_0126*IT_0148;
    const ccomplex_t IT_0150 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0151 = IT_0020*IT_0150;
    const ccomplex_t IT_0152 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0153 = IT_0020*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0151 + 0.5*IT_0154);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = cpow((-2)*s_23 + IT_0001 + IT_0031 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0158 = IT_0156*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0161 = IT_0127*IT_0160;
    const ccomplex_t IT_0162 = IT_0130*IT_0160;
    const ccomplex_t IT_0163 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0164 = IT_0127*IT_0163;
    const ccomplex_t IT_0165 = IT_0130*IT_0163;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*(IT_0161 + IT_0162 + 
      -IT_0164 + -IT_0165);
    const ccomplex_t IT_0167 = (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = m_N_2*IT_0167;
    const ccomplex_t IT_0169 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0170 = IT_0127*IT_0169;
    const ccomplex_t IT_0171 = IT_0130*IT_0169;
    const ccomplex_t IT_0172 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0173 = IT_0127*IT_0172;
    const ccomplex_t IT_0174 = IT_0130*IT_0172;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0170 + IT_0171 + 
      -IT_0173 + -IT_0174);
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = m_N_3*IT_0176;
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = IT_0168 + IT_0178;
    const ccomplex_t IT_0180 = IT_0159*IT_0179;
    const ccomplex_t IT_0181 = -IT_0149 + -IT_0180;
    const ccomplex_t IT_0182 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0183 = IT_0020*IT_0182;
    const ccomplex_t IT_0184 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0185 = IT_0020*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0183 + (-0.5)*IT_0186);
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0190 = IT_0127*IT_0189;
    const ccomplex_t IT_0191 = IT_0130*IT_0189;
    const ccomplex_t IT_0192 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0193 = IT_0130*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*(IT_0190 + -IT_0191 + (-2)
      *IT_0193);
    const ccomplex_t IT_0195 = (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = IT_0188*IT_0195;
    const ccomplex_t IT_0197 = cpow((-2)*s_13 + IT_0001 + IT_0030 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0198 = IT_0196*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = 2*IT_0199;
    const ccomplex_t IT_0201 = cpow((-2)*s_13 + IT_0001 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0202 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0203 = IT_0127*IT_0202;
    const ccomplex_t IT_0204 = IT_0130*IT_0202;
    const ccomplex_t IT_0205 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0206 = IT_0130*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*(IT_0203 + -IT_0204 + (-2)
      *IT_0206);
    const ccomplex_t IT_0208 = (-0.5)*IT_0207;
    const ccomplex_t IT_0209 = IT_0070*IT_0208;
    const ccomplex_t IT_0210 = IT_0201*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = 2*IT_0211;
    const ccomplex_t IT_0213 = IT_0200 + IT_0212;
    const ccomplex_t IT_0214 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0215 = IT_0127*IT_0214;
    const ccomplex_t IT_0216 = IT_0130*IT_0214;
    const ccomplex_t IT_0217 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0218 = IT_0127*IT_0217;
    const ccomplex_t IT_0219 = IT_0130*IT_0217;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*(IT_0215 + IT_0216 + 
      -IT_0218 + -IT_0219);
    const ccomplex_t IT_0221 = (-0.5)*IT_0220;
    const ccomplex_t IT_0222 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0223 = IT_0020*IT_0222;
    const ccomplex_t IT_0224 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0225 = IT_0020*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + (-0.5)*IT_0226);
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = cpow((-2)*s_23 + IT_0001 + IT_0031 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0230 = IT_0228*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = IT_0221*IT_0231;
    const ccomplex_t IT_0233 = (-2)*IT_0232;
    const ccomplex_t IT_0234 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0235 = IT_0020*IT_0234;
    const ccomplex_t IT_0236 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0237 = IT_0020*IT_0236;
    const ccomplex_t IT_0238 = 1.4142135623731*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*(IT_0235 + (-0.5)*IT_0238);
    const ccomplex_t IT_0240 = -IT_0239;
    const ccomplex_t IT_0241 = IT_0117*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = IT_0136*IT_0242;
    const ccomplex_t IT_0244 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0245 = IT_0127*IT_0244;
    const ccomplex_t IT_0246 = IT_0130*IT_0244;
    const ccomplex_t IT_0247 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0248 = IT_0127*IT_0247;
    const ccomplex_t IT_0249 = IT_0130*IT_0247;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0245 + IT_0246 + 
      -IT_0248 + -IT_0249);
    const ccomplex_t IT_0251 = 0.5*IT_0250;
    const ccomplex_t IT_0252 = cpow((-2)*s_23 + IT_0001 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0253 = IT_0070*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0251*IT_0254;
    const ccomplex_t IT_0256 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0257 = IT_0020*IT_0256;
    const ccomplex_t IT_0258 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0259 = IT_0020*IT_0258;
    const ccomplex_t IT_0260 = 1.4142135623731*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*(IT_0257 + (-0.5)*IT_0260);
    const ccomplex_t IT_0262 = -IT_0261;
    const ccomplex_t IT_0263 = IT_0157*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*IT_0263;
    const ccomplex_t IT_0265 = IT_0167*IT_0264;
    const ccomplex_t IT_0266 = (-2)*IT_0021;
    const ccomplex_t IT_0267 = IT_0072*IT_0266;
    const ccomplex_t IT_0268 = 2*IT_0243 + (-2)*IT_0255 + 2*IT_0265 + -IT_0267;
    const ccomplex_t IT_0269 = pow(s_23, 2);
    const ccomplex_t IT_0270 = IT_0000*IT_0269;
    const ccomplex_t IT_0271 = pow(s_24, 2);
    const ccomplex_t IT_0272 = IT_0004*IT_0271;
    const ccomplex_t IT_0273 = IT_0270 + IT_0272;
    const ccomplex_t IT_0274 = s_23*s_24*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0275 = -IT_0031 + -IT_0274;
    const ccomplex_t IT_0276 = IT_0273 + IT_0275;
    const ccomplex_t IT_0277 = m_C_1*IT_0276;
    const ccomplex_t IT_0278 = 2*IT_0277;
    const ccomplex_t IT_0279 = IT_0036*IT_0278;
    const ccomplex_t IT_0280 = m_C_1*IT_0021;
    const ccomplex_t IT_0281 = IT_0034*IT_0280;
    const ccomplex_t IT_0282 = s_14*s_23*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0283 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0284 = s_14*s_24*IT_0004;
    const ccomplex_t IT_0285 = -IT_0283 + -IT_0284;
    const ccomplex_t IT_0286 = s_12 + IT_0282 + IT_0285;
    const ccomplex_t IT_0287 = m_N_3*IT_0286;
    const ccomplex_t IT_0288 = 2*IT_0287;
    const ccomplex_t IT_0289 = IT_0281*IT_0288;
    const ccomplex_t IT_0290 = s_34*m_C_1*m_N_3*IT_0004;
    const ccomplex_t IT_0291 = s_23*s_34*IT_0000;
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = s_24 + IT_0292;
    const ccomplex_t IT_0294 = IT_0290*IT_0293;
    const ccomplex_t IT_0295 = 2*IT_0294;
    const ccomplex_t IT_0296 = IT_0109*IT_0295;
    const ccomplex_t IT_0297 = cos(beta);
    const ccomplex_t IT_0298 = IT_0020*IT_0297;
    const ccomplex_t IT_0299 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0300 = IT_0298*IT_0299;
    const ccomplex_t IT_0301 = (-2)*IT_0300;
    const ccomplex_t IT_0302 = IT_0080*IT_0297;
    const ccomplex_t IT_0303 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0304 = IT_0302*IT_0303;
    const ccomplex_t IT_0305 = 1.4142135623731*IT_0304;
    const ccomplex_t IT_0306 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0307 = IT_0298*IT_0306;
    const ccomplex_t IT_0308 = 1.4142135623731*IT_0307;
    const ccomplex_t IT_0309 = IT_0305 + IT_0308;
    const ccomplex_t IT_0310 = IT_0301 + IT_0309;
    const ccomplex_t IT_0311 = (0 + _Complex_I*1)*IT_0310;
    const ccomplex_t IT_0312 = (-0.5)*IT_0311;
    const ccomplex_t IT_0313 = IT_0093*IT_0312;
    const ccomplex_t IT_0314 = IT_0074*IT_0313;
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*IT_0314;
    const ccomplex_t IT_0316 = IT_0072*IT_0280;
    const ccomplex_t IT_0317 = IT_0315 + -IT_0316;
    const ccomplex_t IT_0318 = (-2)*IT_0277;
    const ccomplex_t IT_0319 = IT_0317*IT_0318;
    const ccomplex_t IT_0320 = s_12*IT_0008;
    const ccomplex_t IT_0321 = s_14*s_24;
    const ccomplex_t IT_0322 = (-2)*IT_0321;
    const ccomplex_t IT_0323 = IT_0320 + IT_0322;
    const ccomplex_t IT_0324 = IT_0004*IT_0323;
    const ccomplex_t IT_0325 = -(IT_0031 + -IT_0270)*(s_12 + (-0.5)*IT_0324);
    const ccomplex_t IT_0326 = (-2)*IT_0325;
    const ccomplex_t IT_0327 = (-2)*IT_0326;
    const ccomplex_t IT_0328 = IT_0233*IT_0327;
    const ccomplex_t IT_0329 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0330 = IT_0127*IT_0329;
    const ccomplex_t IT_0331 = IT_0130*IT_0329;
    const ccomplex_t IT_0332 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0333 = IT_0130*IT_0332;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*(IT_0330 + -IT_0331 + (-2)
      *IT_0333);
    const ccomplex_t IT_0335 = 0.5*IT_0334;
    const ccomplex_t IT_0336 = IT_0029*IT_0335;
    const ccomplex_t IT_0337 = IT_0201*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*IT_0337;
    const ccomplex_t IT_0339 = m_C_1*IT_0338;
    const ccomplex_t IT_0340 = s_13*s_23*IT_0000*IT_0008;
    const ccomplex_t IT_0341 = -IT_0340;
    const ccomplex_t IT_0342 = IT_0320 + IT_0341;
    const ccomplex_t IT_0343 = IT_0004*IT_0342;
    const ccomplex_t IT_0344 = -IT_0283 + (-0.25)*IT_0343;
    const ccomplex_t IT_0345 = s_12 + IT_0344;
    const ccomplex_t IT_0346 = m_N_3*IT_0345;
    const ccomplex_t IT_0347 = (-4)*IT_0346;
    const ccomplex_t IT_0348 = (-2)*IT_0347;
    const ccomplex_t IT_0349 = IT_0339*IT_0348;
    const ccomplex_t IT_0350 = s_12*s_23;
    const ccomplex_t IT_0351 = s_13*IT_0031;
    const ccomplex_t IT_0352 = -IT_0351;
    const ccomplex_t IT_0353 = IT_0350 + IT_0352;
    const ccomplex_t IT_0354 = s_12*s_23*IT_0008;
    const ccomplex_t IT_0355 = s_13*IT_0008*IT_0031;
    const ccomplex_t IT_0356 = s_13*IT_0271;
    const ccomplex_t IT_0357 = s_12*s_24*s_34;
    const ccomplex_t IT_0358 = -IT_0355 + 2*IT_0356 + (-2)*IT_0357;
    const ccomplex_t IT_0359 = IT_0354 + IT_0358;
    const ccomplex_t IT_0360 = IT_0004*IT_0359;
    const ccomplex_t IT_0361 = (-0.5)*IT_0360;
    const ccomplex_t IT_0362 = IT_0353 + IT_0361;
    const ccomplex_t IT_0363 = (-4)*IT_0199;
    const ccomplex_t IT_0364 = (-4)*IT_0211;
    const ccomplex_t IT_0365 = IT_0363 + IT_0364;
    const ccomplex_t IT_0366 = IT_0362*IT_0365;
    const ccomplex_t IT_0367 = IT_0279 + IT_0289 + IT_0296 + IT_0319 + IT_0328
       + IT_0349 + IT_0366;
    const ccomplex_t IT_0368 = IT_0072*IT_0108;
    const ccomplex_t IT_0369 = s_12*IT_0001;
    const ccomplex_t IT_0370 = s_14*s_23*s_34*IT_0004;
    const ccomplex_t IT_0371 = s_12*IT_0004*IT_0041;
    const ccomplex_t IT_0372 = -IT_0371;
    const ccomplex_t IT_0373 = s_12*IT_0101;
    const ccomplex_t IT_0374 = s_14*s_23*s_34*IT_0001*IT_0004;
    const ccomplex_t IT_0375 = s_13*s_24*s_34*IT_0001*IT_0004;
    const ccomplex_t IT_0376 = IT_0373 + IT_0374 + IT_0375;
    const ccomplex_t IT_0377 = s_13*s_23*IT_0001;
    const ccomplex_t IT_0378 = s_12*IT_0001*IT_0004*IT_0041;
    const ccomplex_t IT_0379 = (-2)*IT_0377 + -IT_0378;
    const ccomplex_t IT_0380 = IT_0376 + IT_0379;
    const ccomplex_t IT_0381 = IT_0000*IT_0380;
    const ccomplex_t IT_0382 = (-0.5)*IT_0381;
    const ccomplex_t IT_0383 = s_13*s_24*s_34*IT_0004;
    const ccomplex_t IT_0384 = -IT_0383;
    const ccomplex_t IT_0385 = IT_0369 + IT_0370 + IT_0372 + IT_0382 + IT_0384;
    const ccomplex_t IT_0386 = (-0.5)*IT_0250;
    const ccomplex_t IT_0387 = IT_0029*IT_0252;
    const ccomplex_t IT_0388 = (0 + _Complex_I*1)*IT_0387;
    const ccomplex_t IT_0389 = IT_0386*IT_0388;
    const ccomplex_t IT_0390 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0391 = IT_0127*IT_0390;
    const ccomplex_t IT_0392 = IT_0130*IT_0390;
    const ccomplex_t IT_0393 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0394 = IT_0127*IT_0393;
    const ccomplex_t IT_0395 = IT_0130*IT_0393;
    const ccomplex_t IT_0396 = (0 + _Complex_I*1)*(IT_0391 + IT_0392 + 
      -IT_0394 + -IT_0395);
    const ccomplex_t IT_0397 = (-0.5)*IT_0396;
    const ccomplex_t IT_0398 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0399 = IT_0020*IT_0398;
    const ccomplex_t IT_0400 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0401 = IT_0020*IT_0400;
    const ccomplex_t IT_0402 = 1.4142135623731*IT_0401;
    const ccomplex_t IT_0403 = (0 + _Complex_I*1)*(IT_0399 + 0.5*IT_0402);
    const ccomplex_t IT_0404 = -IT_0403;
    const ccomplex_t IT_0405 = IT_0229*IT_0404;
    const ccomplex_t IT_0406 = (0 + _Complex_I*1)*IT_0405;
    const ccomplex_t IT_0407 = IT_0397*IT_0406;
    const ccomplex_t IT_0408 = IT_0034*IT_0266;
    const ccomplex_t IT_0409 = (-2)*IT_0389 + (-2)*IT_0407 + -IT_0408;
    const ccomplex_t IT_0410 = IT_0295*IT_0409;
    const ccomplex_t IT_0411 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0412 = IT_0020*IT_0411;
    const ccomplex_t IT_0413 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0414 = IT_0020*IT_0413;
    const ccomplex_t IT_0415 = 1.4142135623731*IT_0414;
    const ccomplex_t IT_0416 = (0 + _Complex_I*1)*(IT_0412 + 0.5*IT_0415);
    const ccomplex_t IT_0417 = -IT_0416;
    const ccomplex_t IT_0418 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0419 = IT_0127*IT_0418;
    const ccomplex_t IT_0420 = IT_0130*IT_0418;
    const ccomplex_t IT_0421 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0422 = IT_0130*IT_0421;
    const ccomplex_t IT_0423 = (0 + _Complex_I*1)*(IT_0419 + -IT_0420 + (-2)
      *IT_0422);
    const ccomplex_t IT_0424 = 0.5*IT_0423;
    const ccomplex_t IT_0425 = IT_0417*IT_0424;
    const ccomplex_t IT_0426 = IT_0197*IT_0425;
    const ccomplex_t IT_0427 = (0 + _Complex_I*1)*IT_0426;
    const ccomplex_t IT_0428 = 2*IT_0338 + -IT_0408 + 2*IT_0427;
    const ccomplex_t IT_0429 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0430 = -IT_0429;
    const ccomplex_t IT_0431 = s_14 + IT_0430;
    const ccomplex_t IT_0432 = IT_0290*IT_0431;
    const ccomplex_t IT_0433 = 2*IT_0432;
    const ccomplex_t IT_0434 = IT_0428*IT_0433;
    const ccomplex_t IT_0435 = IT_0410 + IT_0434;
    const ccomplex_t IT_0436 = s_14*s_34*IT_0004;
    const ccomplex_t IT_0437 = IT_0000*IT_0004;
    const ccomplex_t IT_0438 = s_13*IT_0001*IT_0008;
    const ccomplex_t IT_0439 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0440 = (-2)*IT_0439;
    const ccomplex_t IT_0441 = IT_0438 + IT_0440;
    const ccomplex_t IT_0442 = IT_0437*IT_0441;
    const ccomplex_t IT_0443 = (-0.666666666666667)*IT_0436 + (
      -0.166666666666667)*IT_0442;
    const ccomplex_t IT_0444 = s_13 + IT_0443;
    const ccomplex_t IT_0445 = m_N_3*IT_0444;
    const ccomplex_t IT_0446 = (-6)*IT_0445;
    const ccomplex_t IT_0447 = m_N_4*IT_0397;
    const ccomplex_t IT_0448 = m_N_3*IT_0221;
    const ccomplex_t IT_0449 = IT_0231*(IT_0447 + -IT_0448);
    const ccomplex_t IT_0450 = -IT_0449;
    const ccomplex_t IT_0451 = m_N_3*IT_0251;
    const ccomplex_t IT_0452 = m_N_3*IT_0386;
    const ccomplex_t IT_0453 = -IT_0452;
    const ccomplex_t IT_0454 = IT_0451 + IT_0453;
    const ccomplex_t IT_0455 = IT_0388*IT_0454;
    const ccomplex_t IT_0456 = m_N_4*IT_0221;
    const ccomplex_t IT_0457 = m_N_3*IT_0397;
    const ccomplex_t IT_0458 = -IT_0457;
    const ccomplex_t IT_0459 = IT_0456 + IT_0458;
    const ccomplex_t IT_0460 = IT_0406*IT_0459;
    const ccomplex_t IT_0461 = -IT_0455 + -IT_0460;
    const ccomplex_t IT_0462 = IT_0437*IT_0438;
    const ccomplex_t IT_0463 = 0.666666666666667*IT_0052 + (-0.166666666666667
      )*IT_0462;
    const ccomplex_t IT_0464 = s_13 + IT_0463;
    const ccomplex_t IT_0465 = m_N_3*IT_0464;
    const ccomplex_t IT_0466 = (-6)*IT_0465;
    const ccomplex_t IT_0467 = IT_0029*IT_0208;
    const ccomplex_t IT_0468 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0469 = IT_0201*IT_0467*IT_0468;
    const ccomplex_t IT_0470 = IT_0195*IT_0417;
    const ccomplex_t IT_0471 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0472 = IT_0197*IT_0470*IT_0471;
    const ccomplex_t IT_0473 = -IT_0469 + -IT_0472;
    const ccomplex_t IT_0474 = s_24*s_34*IT_0001;
    const ccomplex_t IT_0475 = (-2)*IT_0474;
    const ccomplex_t IT_0476 = IT_0111 + IT_0475;
    const ccomplex_t IT_0477 = IT_0437*IT_0476;
    const ccomplex_t IT_0478 = (-0.666666666666667)*IT_0013 +
       0.166666666666667*IT_0477;
    const ccomplex_t IT_0479 = s_23 + IT_0478;
    const ccomplex_t IT_0480 = m_C_1*IT_0479;
    const ccomplex_t IT_0481 = 6*IT_0480;
    const ccomplex_t IT_0482 = IT_0070*IT_0335;
    const ccomplex_t IT_0483 = IT_0201*IT_0468*IT_0482;
    const ccomplex_t IT_0484 = IT_0188*IT_0424;
    const ccomplex_t IT_0485 = IT_0197*IT_0471*IT_0484;
    const ccomplex_t IT_0486 = m_C_1*IT_0199;
    const ccomplex_t IT_0487 = m_C_1*IT_0211;
    const ccomplex_t IT_0488 = IT_0483 + IT_0485 + IT_0486 + IT_0487;
    const ccomplex_t IT_0489 = (-2)*IT_0199;
    const ccomplex_t IT_0490 = (-2)*IT_0211;
    const ccomplex_t IT_0491 = IT_0489 + IT_0490;
    const ccomplex_t IT_0492 = s_13*s_23;
    const ccomplex_t IT_0493 = s_12*IT_0001*IT_0008;
    const ccomplex_t IT_0494 = s_13*s_23*IT_0008;
    const ccomplex_t IT_0495 = s_13*s_24*s_34;
    const ccomplex_t IT_0496 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0497 = (-2)*IT_0494 + 4*IT_0495 + (-2)*IT_0496;
    const ccomplex_t IT_0498 = IT_0493 + IT_0497;
    const ccomplex_t IT_0499 = IT_0004*IT_0498;
    const ccomplex_t IT_0500 = -IT_0499;
    const ccomplex_t IT_0501 = s_12*IT_0008*IT_0101;
    const ccomplex_t IT_0502 = s_14*s_24*IT_0101;
    const ccomplex_t IT_0503 = (-2)*IT_0502;
    const ccomplex_t IT_0504 = IT_0501 + IT_0503;
    const ccomplex_t IT_0505 = IT_0437*IT_0504;
    const ccomplex_t IT_0506 = 0.5*IT_0505;
    const ccomplex_t IT_0507 = IT_0369 + (-4)*IT_0492 + IT_0500 + IT_0506;
    const ccomplex_t IT_0508 = (-4)*conj(IT_0211);
    const ccomplex_t IT_0509 = -IT_0338 + -IT_0427;
    const ccomplex_t IT_0510 = IT_0008*IT_0101;
    const ccomplex_t IT_0511 = IT_0437*IT_0510;
    const ccomplex_t IT_0512 = 0.125*IT_0511;
    const ccomplex_t IT_0513 = IT_0001 + IT_0512;
    const ccomplex_t IT_0514 = IT_0037*IT_0513;
    const ccomplex_t IT_0515 = 8*IT_0514;
    const ccomplex_t IT_0516 = 2*IT_0063;
    const ccomplex_t IT_0517 = IT_0317*IT_0516;
    const ccomplex_t IT_0518 = 2*IT_0017;
    const ccomplex_t IT_0519 = IT_0281*IT_0518;
    const ccomplex_t IT_0520 = IT_0021*IT_0072;
    const ccomplex_t IT_0521 = -IT_0520;
    const ccomplex_t IT_0522 = s_14*s_34*IT_0001*IT_0004;
    const ccomplex_t IT_0523 = -IT_0522;
    const ccomplex_t IT_0524 = IT_0051 + IT_0523;
    const ccomplex_t IT_0525 = IT_0000*IT_0524;
    const ccomplex_t IT_0526 = -IT_0436 + (-0.25)*IT_0525;
    const ccomplex_t IT_0527 = s_13 + IT_0526;
    const ccomplex_t IT_0528 = IT_0037*IT_0527;
    const ccomplex_t IT_0529 = (-4)*IT_0528;
    const ccomplex_t IT_0530 = 2*IT_0529;
    const ccomplex_t IT_0531 = IT_0521*IT_0530;
    const ccomplex_t IT_0532 = s_24*s_34*IT_0004;
    const ccomplex_t IT_0533 = s_24*s_34*IT_0001*IT_0004;
    const ccomplex_t IT_0534 = -IT_0533;
    const ccomplex_t IT_0535 = IT_0002 + IT_0534;
    const ccomplex_t IT_0536 = IT_0000*IT_0535;
    const ccomplex_t IT_0537 = -IT_0532 + (-0.25)*IT_0536;
    const ccomplex_t IT_0538 = s_23 + IT_0537;
    const ccomplex_t IT_0539 = IT_0037*IT_0538;
    const ccomplex_t IT_0540 = (-4)*IT_0539;
    const ccomplex_t IT_0541 = 2*IT_0540;
    const ccomplex_t IT_0542 = IT_0521*IT_0541;
    const ccomplex_t IT_0543 = 2*IT_0107;
    const ccomplex_t IT_0544 = IT_0368*IT_0543;
    const ccomplex_t IT_0545 = IT_0049*IT_0072;
    const ccomplex_t IT_0546 = 2*IT_0048;
    const ccomplex_t IT_0547 = IT_0545*IT_0546;
    const ccomplex_t IT_0548 = (-2)*IT_0377;
    const ccomplex_t IT_0549 = IT_0373 + IT_0548;
    const ccomplex_t IT_0550 = IT_0000*IT_0549;
    const ccomplex_t IT_0551 = 0.5*IT_0550;
    const ccomplex_t IT_0552 = s_14*s_23*s_34;
    const ccomplex_t IT_0553 = s_12*IT_0041;
    const ccomplex_t IT_0554 = (-2)*IT_0494 + 2*IT_0495 + 2*IT_0552 + (-2)
      *IT_0553;
    const ccomplex_t IT_0555 = IT_0493 + IT_0554;
    const ccomplex_t IT_0556 = IT_0004*IT_0555;
    const ccomplex_t IT_0557 = (-0.5)*IT_0556;
    const ccomplex_t IT_0558 = s_14*s_23*IT_0001;
    const ccomplex_t IT_0559 = s_13*s_24*IT_0001;
    const ccomplex_t IT_0560 = IT_0558 + IT_0559;
    const ccomplex_t IT_0561 = s_12*s_34*IT_0001;
    const ccomplex_t IT_0562 = -IT_0561;
    const ccomplex_t IT_0563 = IT_0560 + IT_0562;
    const ccomplex_t IT_0564 = IT_0005*IT_0563;
    const ccomplex_t IT_0565 = 0.5*IT_0564;
    const ccomplex_t IT_0566 = (-2)*IT_0492;
    const ccomplex_t IT_0567 = IT_0369 + IT_0551 + IT_0557 + IT_0565 + IT_0566;
    const ccomplex_t IT_0568 = 4*IT_0567;
    const ccomplex_t IT_0569 = IT_0050*IT_0568;
    const ccomplex_t IT_0570 = (-2)*IT_0063;
    const ccomplex_t IT_0571 = IT_0036*IT_0570;
    const ccomplex_t IT_0572 = (-2)*IT_0017;
    const ccomplex_t IT_0573 = IT_0097*IT_0572;
    const ccomplex_t IT_0574 = IT_0021*IT_0034;
    const ccomplex_t IT_0575 = -IT_0574;
    const ccomplex_t IT_0576 = IT_0001*IT_0004*IT_0271;
    const ccomplex_t IT_0577 = IT_0001*IT_0031;
    const ccomplex_t IT_0578 = s_23*s_24*s_34*IT_0004;
    const ccomplex_t IT_0579 = 0.5*IT_0576 + (-0.5)*IT_0577 + -IT_0578;
    const ccomplex_t IT_0580 = IT_0269 + IT_0579;
    const ccomplex_t IT_0581 = IT_0000*IT_0580;
    const ccomplex_t IT_0582 = -IT_0272;
    const ccomplex_t IT_0583 = IT_0031 + IT_0581 + IT_0582;
    const ccomplex_t IT_0584 = m_C_1*IT_0583;
    const ccomplex_t IT_0585 = IT_0254*(IT_0451 + -IT_0452);
    const ccomplex_t IT_0586 = -IT_0585;
    const ccomplex_t IT_0587 = m_N_1*IT_0145;
    const ccomplex_t IT_0588 = m_N_3*IT_0136;
    const ccomplex_t IT_0589 = -IT_0588;
    const ccomplex_t IT_0590 = IT_0587 + IT_0589;
    const ccomplex_t IT_0591 = IT_0242*IT_0590;
    const ccomplex_t IT_0592 = m_N_2*IT_0176;
    const ccomplex_t IT_0593 = m_N_3*IT_0167;
    const ccomplex_t IT_0594 = -IT_0593;
    const ccomplex_t IT_0595 = IT_0592 + IT_0594;
    const ccomplex_t IT_0596 = IT_0264*IT_0595;
    const ccomplex_t IT_0597 = -IT_0591 + -IT_0596;
    const ccomplex_t IT_0598 = IT_0389 + IT_0407;
    const ccomplex_t IT_0599 = IT_0126*IT_0145;
    const ccomplex_t IT_0600 = IT_0159*IT_0176;
    const ccomplex_t IT_0601 = -IT_0599 + -IT_0600;
    const ccomplex_t IT_0602 = s_12*s_24*s_34*IT_0004;
    const ccomplex_t IT_0603 = (-0.5)*IT_0602;
    const ccomplex_t IT_0604 = s_14*s_23*s_24*IT_0001*IT_0004;
    const ccomplex_t IT_0605 = s_13*IT_0001*IT_0004*IT_0271;
    const ccomplex_t IT_0606 = IT_0604 + IT_0605;
    const ccomplex_t IT_0607 = s_12*s_24*s_34*IT_0001*IT_0004;
    const ccomplex_t IT_0608 = s_13*IT_0001*IT_0031;
    const ccomplex_t IT_0609 = -IT_0607 + -IT_0608;
    const ccomplex_t IT_0610 = IT_0606 + IT_0609;
    const ccomplex_t IT_0611 = IT_0000*IT_0610;
    const ccomplex_t IT_0612 = (-0.25)*IT_0611;
    const ccomplex_t IT_0613 = (-0.5)*IT_0351;
    const ccomplex_t IT_0614 = s_13*IT_0004*IT_0271;
    const ccomplex_t IT_0615 = 0.5*IT_0614;
    const ccomplex_t IT_0616 = s_14*s_23*s_24*IT_0004;
    const ccomplex_t IT_0617 = (-0.5)*IT_0616;
    const ccomplex_t IT_0618 = IT_0350 + IT_0603 + IT_0612 + IT_0613 + IT_0615
       + IT_0617;
    const ccomplex_t IT_0619 = s_12*s_13;
    const ccomplex_t IT_0620 = s_12*s_14*s_34*IT_0004;
    const ccomplex_t IT_0621 = (-0.5)*IT_0620;
    const ccomplex_t IT_0622 = pow(s_14, 2);
    const ccomplex_t IT_0623 = s_23*IT_0001*IT_0004*IT_0622;
    const ccomplex_t IT_0624 = s_13*s_14*s_24*IT_0001*IT_0004;
    const ccomplex_t IT_0625 = IT_0623 + IT_0624;
    const ccomplex_t IT_0626 = s_12*s_14*s_34*IT_0001*IT_0004;
    const ccomplex_t IT_0627 = s_23*IT_0001*IT_0030;
    const ccomplex_t IT_0628 = -IT_0626 + -IT_0627;
    const ccomplex_t IT_0629 = IT_0625 + IT_0628;
    const ccomplex_t IT_0630 = IT_0000*IT_0629;
    const ccomplex_t IT_0631 = (-0.25)*IT_0630;
    const ccomplex_t IT_0632 = s_23*IT_0030;
    const ccomplex_t IT_0633 = (-0.5)*IT_0632;
    const ccomplex_t IT_0634 = s_13*s_14*s_24*IT_0004;
    const ccomplex_t IT_0635 = (-0.5)*IT_0634;
    const ccomplex_t IT_0636 = s_23*IT_0004*IT_0622;
    const ccomplex_t IT_0637 = 0.5*IT_0636;
    const ccomplex_t IT_0638 = IT_0619 + IT_0621 + IT_0631 + IT_0633 + IT_0635
       + IT_0637;
    const ccomplex_t IT_0639 = -IT_0614;
    const ccomplex_t IT_0640 = -IT_0602;
    const ccomplex_t IT_0641 = (-0.5)*IT_0611;
    const ccomplex_t IT_0642 = IT_0351 + IT_0616 + IT_0639 + IT_0640 + IT_0641;
    const ccomplex_t IT_0643 = 0.5*IT_0642;
    const ccomplex_t IT_0644 = IT_0638 + IT_0643;
    const ccomplex_t IT_0645 = IT_0529 + IT_0540;
    const ccomplex_t IT_0646 = s_13*s_23*IT_0000*IT_0004*IT_0622;
    const ccomplex_t IT_0647 = pow(s_13, 2);
    const ccomplex_t IT_0648 = s_14*s_24*IT_0000*IT_0004*IT_0647;
    const ccomplex_t IT_0649 = s_12*IT_0030;
    const ccomplex_t IT_0650 = s_14*s_24*IT_0004*IT_0030;
    const ccomplex_t IT_0651 = -IT_0650;
    const ccomplex_t IT_0652 = s_13*s_23*IT_0000*IT_0030;
    const ccomplex_t IT_0653 = -IT_0652;
    const ccomplex_t IT_0654 = s_12*s_13*s_14*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0655 = -IT_0654;
    const ccomplex_t IT_0656 = IT_0646 + IT_0648 + IT_0649 + IT_0651 + IT_0653
       + IT_0655;
    const ccomplex_t IT_0657 = 2*IT_0656;
    const ccomplex_t IT_0658 = s_13*s_14*s_23*s_24*IT_0000*IT_0004;
    const ccomplex_t IT_0659 = IT_0030*IT_0031;
    const ccomplex_t IT_0660 = IT_0000*IT_0004*IT_0271*IT_0647;
    const ccomplex_t IT_0661 = IT_0000*IT_0031*IT_0647;
    const ccomplex_t IT_0662 = -IT_0661;
    const ccomplex_t IT_0663 = IT_0004*IT_0030*IT_0271;
    const ccomplex_t IT_0664 = -IT_0663;
    const ccomplex_t IT_0665 = s_12*s_13*s_24*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0666 = -IT_0665;
    const ccomplex_t IT_0667 = IT_0658 + IT_0659 + IT_0660 + IT_0662 + IT_0664
       + IT_0666;
    const ccomplex_t IT_0668 = 2*IT_0667;
    const ccomplex_t IT_0669 = IT_0657 + IT_0668;
    const ccomplex_t IT_0670 = s_13*s_14*s_23*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0671 = s_23*IT_0000*IT_0647;
    const ccomplex_t IT_0672 = (-2)*IT_0671;
    const ccomplex_t IT_0673 = IT_0619 + IT_0672;
    const ccomplex_t IT_0674 = -IT_0634;
    const ccomplex_t IT_0675 = -IT_0636;
    const ccomplex_t IT_0676 = IT_0673 + IT_0674 + IT_0675;
    const ccomplex_t IT_0677 = 2*IT_0632 + 4*IT_0670 + 2*IT_0676;
    const ccomplex_t IT_0678 = s_13*s_23*s_24*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0679 = s_13*IT_0000*IT_0269;
    const ccomplex_t IT_0680 = (-2)*IT_0679;
    const ccomplex_t IT_0681 = IT_0350 + IT_0680;
    const ccomplex_t IT_0682 = -IT_0616;
    const ccomplex_t IT_0683 = IT_0639 + IT_0681 + IT_0682;
    const ccomplex_t IT_0684 = 2*IT_0351 + 4*IT_0678 + 2*IT_0683;
    const ccomplex_t IT_0685 = IT_0677 + IT_0684;
    const ccomplex_t IT_0686 = IT_0000*IT_0647;
    const ccomplex_t IT_0687 = IT_0004*IT_0622;
    const ccomplex_t IT_0688 = IT_0686 + IT_0687;
    const ccomplex_t IT_0689 = s_13*s_14*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0690 = -IT_0030 + -IT_0689;
    const ccomplex_t IT_0691 = IT_0688 + IT_0690;
    const ccomplex_t IT_0692 = m_N_3*IT_0691;
    const ccomplex_t IT_0693 = 2*IT_0692;
    const ccomplex_t IT_0694 = s_13*s_24*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0695 = s_12 + IT_0285 + IT_0694;
    const ccomplex_t IT_0696 = m_N_3*IT_0695;
    const ccomplex_t IT_0697 = (-2)*IT_0696;
    const ccomplex_t IT_0698 = IT_0693 + IT_0697;
    const ccomplex_t IT_0699 = IT_0369 + IT_0566;
    const ccomplex_t IT_0700 = IT_0000*IT_0699;
    const ccomplex_t IT_0701 = (s_13 + -IT_0436)*(s_12 + (-0.5)*IT_0700);
    const ccomplex_t IT_0702 = (-2)*IT_0701;
    const ccomplex_t IT_0703 = 2*IT_0702;
    const ccomplex_t IT_0704 = (s_23 + -IT_0532)*(s_12 + (-0.5)*IT_0700);
    const ccomplex_t IT_0705 = (-2)*IT_0704;
    const ccomplex_t IT_0706 = 2*IT_0705;
    const ccomplex_t IT_0707 = IT_0703 + IT_0706;
    const ccomplex_t IT_0708 = -IT_0577;
    const ccomplex_t IT_0709 = IT_0576 + IT_0708;
    const ccomplex_t IT_0710 = IT_0000*IT_0709;
    const ccomplex_t IT_0711 = -IT_0031;
    const ccomplex_t IT_0712 = IT_0272 + IT_0711;
    const ccomplex_t IT_0713 = (-4)*IT_0712;
    const ccomplex_t IT_0714 = IT_0710 + IT_0713;
    const ccomplex_t IT_0715 = m_C_1*IT_0714;
    const ccomplex_t IT_0716 = 2*IT_0715;
    const ccomplex_t IT_0717 = s_14*s_24*IT_0001*IT_0004;
    const ccomplex_t IT_0718 = -IT_0717;
    const ccomplex_t IT_0719 = IT_0369 + IT_0718;
    const ccomplex_t IT_0720 = IT_0000*IT_0719;
    const ccomplex_t IT_0721 = -IT_0284 + (-0.25)*IT_0720;
    const ccomplex_t IT_0722 = s_12 + IT_0721;
    const ccomplex_t IT_0723 = m_C_1*IT_0722;
    const ccomplex_t IT_0724 = (-4)*IT_0723;
    const ccomplex_t IT_0725 = (-2)*IT_0724;
    const ccomplex_t IT_0726 = IT_0716 + IT_0725;
    const ccomplex_t IT_0727 = IT_0001*IT_0004*IT_0622;
    const ccomplex_t IT_0728 = IT_0001*IT_0030;
    const ccomplex_t IT_0729 = -IT_0728;
    const ccomplex_t IT_0730 = IT_0727 + IT_0729;
    const ccomplex_t IT_0731 = IT_0000*IT_0730;
    const ccomplex_t IT_0732 = -IT_0030;
    const ccomplex_t IT_0733 = IT_0687 + IT_0732;
    const ccomplex_t IT_0734 = (-4)*IT_0733;
    const ccomplex_t IT_0735 = IT_0731 + IT_0734;
    const ccomplex_t IT_0736 = m_N_3*IT_0735;
    const ccomplex_t IT_0737 = 2*IT_0736;
    const ccomplex_t IT_0738 = m_N_3*IT_0722;
    const ccomplex_t IT_0739 = (-4)*IT_0738;
    const ccomplex_t IT_0740 = (-2)*IT_0739;
    const ccomplex_t IT_0741 = IT_0737 + IT_0740;
    const ccomplex_t IT_0742 = 2*IT_0739;
    const ccomplex_t IT_0743 = (-2)*IT_0736;
    const ccomplex_t IT_0744 = IT_0742 + IT_0743;
    const ccomplex_t IT_0745 = (-2)*IT_0199 + (-2)*IT_0211 + -IT_0267;
    const ccomplex_t IT_0746 = IT_0037*IT_0691;
    const ccomplex_t IT_0747 = 2*IT_0746;
    const ccomplex_t IT_0748 = IT_0037*IT_0695;
    const ccomplex_t IT_0749 = (-2)*IT_0748;
    const ccomplex_t IT_0750 = IT_0747 + IT_0749;
    const ccomplex_t IT_0751 = 4*IT_0584;
    const ccomplex_t IT_0752 = 2*IT_0370 + (-2)*IT_0492 + -IT_0717;
    const ccomplex_t IT_0753 = IT_0369 + IT_0752;
    const ccomplex_t IT_0754 = IT_0000*IT_0753;
    const ccomplex_t IT_0755 = -IT_0284 + (-0.5)*IT_0754;
    const ccomplex_t IT_0756 = s_12 + IT_0755;
    const ccomplex_t IT_0757 = m_C_1*IT_0756;
    const ccomplex_t IT_0758 = (-2)*IT_0757;
    const ccomplex_t IT_0759 = (-2)*IT_0758;
    const ccomplex_t IT_0760 = IT_0751 + IT_0759;
    const ccomplex_t IT_0761 = s_13*s_14*s_34*IT_0004;
    const ccomplex_t IT_0762 = 0.5*IT_0727 + (-0.5)*IT_0728 + -IT_0761;
    const ccomplex_t IT_0763 = IT_0647 + IT_0762;
    const ccomplex_t IT_0764 = IT_0000*IT_0763;
    const ccomplex_t IT_0765 = -IT_0687;
    const ccomplex_t IT_0766 = IT_0030 + IT_0764 + IT_0765;
    const ccomplex_t IT_0767 = m_N_3*IT_0766;
    const ccomplex_t IT_0768 = 4*IT_0767;
    const ccomplex_t IT_0769 = 2*IT_0383 + (-2)*IT_0492 + -IT_0717;
    const ccomplex_t IT_0770 = IT_0369 + IT_0769;
    const ccomplex_t IT_0771 = IT_0000*IT_0770;
    const ccomplex_t IT_0772 = -IT_0284 + (-0.5)*IT_0771;
    const ccomplex_t IT_0773 = s_12 + IT_0772;
    const ccomplex_t IT_0774 = m_N_3*IT_0773;
    const ccomplex_t IT_0775 = (-2)*IT_0774;
    const ccomplex_t IT_0776 = (-2)*IT_0775;
    const ccomplex_t IT_0777 = IT_0768 + IT_0776;
    const ccomplex_t IT_0778 = -IT_0436;
    const ccomplex_t IT_0779 = s_13 + IT_0778;
    const ccomplex_t IT_0780 = IT_0037*IT_0779;
    const ccomplex_t IT_0781 = (-6)*IT_0780;
    const ccomplex_t IT_0782 = -IT_0532;
    const ccomplex_t IT_0783 = s_23 + IT_0782;
    const ccomplex_t IT_0784 = IT_0037*IT_0783;
    const ccomplex_t IT_0785 = (-6)*IT_0784;
    const ccomplex_t IT_0786 = IT_0781 + IT_0785;
    const ccomplex_t IT_0787 = m_C_1*IT_0286;
    const ccomplex_t IT_0788 = 2*IT_0787;
    const ccomplex_t IT_0789 = IT_0278*IT_0281;
    const ccomplex_t IT_0790 = (-2)*IT_0787;
    const ccomplex_t IT_0791 = IT_0097*IT_0318;
    const ccomplex_t IT_0792 = 2*IT_0696;
    const ccomplex_t IT_0793 = (-2)*IT_0692;
    const ccomplex_t IT_0794 = IT_0584*(IT_0586 + IT_0597) + 2*(IT_0598 +
       IT_0601)*IT_0618 + 2*IT_0509*IT_0644 + (-0.25)*IT_0213*IT_0645 + (-0.25)
      *IT_0428*IT_0669 + (-0.25)*IT_0050*IT_0685 + (-0.25)*IT_0317*IT_0698 + (
      -0.25)*IT_0109*IT_0707 + (-0.25)*IT_0473*IT_0726 + (-0.25)*IT_0461*IT_0741
       + (-0.25)*IT_0181*IT_0744 + (-0.25)*IT_0745*IT_0750 + (-0.25)*IT_0450
      *IT_0760 + (-0.25)*IT_0488*IT_0777 + (-0.25)*IT_0368*IT_0786 + (-0.25)
      *IT_0097*IT_0788 + (-0.25)*IT_0789 + (-0.25)*IT_0281*IT_0790 + (-0.25)
      *IT_0791 + (-0.25)*IT_0036*(IT_0792 + IT_0793);
    const ccomplex_t IT_0795 = 2*IT_0599 + 2*IT_0600;
    const ccomplex_t IT_0796 = conj(IT_0409) + conj(IT_0795);
    const ccomplex_t IT_0797 = (-2)*IT_0287;
    const ccomplex_t IT_0798 = IT_0317*IT_0797;
    const ccomplex_t IT_0799 = conj(IT_0598) + conj(IT_0601);
    const ccomplex_t IT_0800 = (-0.166666666666667)*IT_0477 + (
      -0.666666666666667)*IT_0532;
    const ccomplex_t IT_0801 = s_23 + IT_0800;
    const ccomplex_t IT_0802 = m_C_1*IT_0801;
    const ccomplex_t IT_0803 = (-6)*IT_0802;
    const ccomplex_t IT_0804 = 2*IT_0803;
    const ccomplex_t IT_0805 = IT_0339*IT_0804;
    const ccomplex_t IT_0806 = (-0.666666666666667)*IT_0059 +
       0.166666666666667*IT_0442;
    const ccomplex_t IT_0807 = s_13 + IT_0806;
    const ccomplex_t IT_0808 = m_N_3*IT_0807;
    const ccomplex_t IT_0809 = 6*IT_0808;
    const ccomplex_t IT_0810 = 2*IT_0809;
    const ccomplex_t IT_0811 = IT_0181*IT_0810;
    const ccomplex_t IT_0812 = IT_0111*IT_0437;
    const ccomplex_t IT_0813 = (-10)*s_23;
    const ccomplex_t IT_0814 = IT_0812 + IT_0813;
    const ccomplex_t IT_0815 = m_C_1*IT_0814;
    const ccomplex_t IT_0816 = 2*IT_0815;
    const ccomplex_t IT_0817 = IT_0597*IT_0816;
    const ccomplex_t IT_0818 = s_14*s_23*s_24;
    const ccomplex_t IT_0819 = s_14*s_34*IT_0031;
    const ccomplex_t IT_0820 = -IT_0355 + (-2)*IT_0818 + 2*IT_0819;
    const ccomplex_t IT_0821 = IT_0354 + IT_0820;
    const ccomplex_t IT_0822 = IT_0004*IT_0821;
    const ccomplex_t IT_0823 = (-0.5)*IT_0822;
    const ccomplex_t IT_0824 = IT_0353 + IT_0823;
    const ccomplex_t IT_0825 = 4*IT_0824;
    const ccomplex_t IT_0826 = IT_0795*IT_0825;
    const ccomplex_t IT_0827 = m_C_1*m_N_3*IT_0004;
    const ccomplex_t IT_0828 = IT_0000*IT_0510;
    const ccomplex_t IT_0829 = (-0.25)*IT_0828;
    const ccomplex_t IT_0830 = IT_0041 + IT_0829;
    const ccomplex_t IT_0831 = IT_0827*IT_0830;
    const ccomplex_t IT_0832 = (-4)*IT_0831;
    const ccomplex_t IT_0833 = IT_0489*IT_0832;
    const ccomplex_t IT_0834 = IT_0490*IT_0832;
    const ccomplex_t IT_0835 = (-2)*IT_0809;
    const ccomplex_t IT_0836 = IT_0461*IT_0835;
    const ccomplex_t IT_0837 = m_N_3*IT_0004;
    const ccomplex_t IT_0838 = IT_0000*IT_0438;
    const ccomplex_t IT_0839 = (-2)*IT_0056 + (-0.5)*IT_0838;
    const ccomplex_t IT_0840 = IT_0055 + IT_0839;
    const ccomplex_t IT_0841 = IT_0837*IT_0840;
    const ccomplex_t IT_0842 = (-2)*IT_0841;
    const ccomplex_t IT_0843 = (-2)*IT_0842;
    const ccomplex_t IT_0844 = IT_0488*IT_0843;
    const ccomplex_t IT_0845 = IT_0097*IT_0516;
    const ccomplex_t IT_0846 = IT_0036*IT_0518;
    const ccomplex_t IT_0847 = IT_0530*IT_0575;
    const ccomplex_t IT_0848 = IT_0541*IT_0575;
    const ccomplex_t IT_0849 = IT_0109*IT_0543;
    const ccomplex_t IT_0850 = IT_0050*IT_0546;
    const ccomplex_t IT_0851 = IT_0586*IT_0810;
    const ccomplex_t IT_0852 = IT_0597*IT_0810;
    const ccomplex_t IT_0853 = IT_0181*IT_0816;
    const ccomplex_t IT_0854 = m_C_1*IT_0427;
    const ccomplex_t IT_0855 = 2*IT_0842;
    const ccomplex_t IT_0856 = IT_0854*IT_0855;
    const ccomplex_t IT_0857 = IT_0339*IT_0855;
    const ccomplex_t IT_0858 = -IT_0632;
    const ccomplex_t IT_0859 = 0.5*IT_0630;
    const ccomplex_t IT_0860 = IT_0620 + IT_0636 + IT_0674 + IT_0858 + IT_0859;
    const ccomplex_t IT_0861 = 4*IT_0860;
    const ccomplex_t IT_0862 = IT_0521*IT_0861;
    const ccomplex_t IT_0863 = IT_0545*IT_0568;
    const ccomplex_t IT_0864 = IT_0233*IT_0825;
    const ccomplex_t IT_0865 = IT_0268*IT_0825;
    const ccomplex_t IT_0866 = 0.25*IT_0505;
    const ccomplex_t IT_0867 = IT_0369 + IT_0866;
    const ccomplex_t IT_0868 = IT_0199*IT_0867;
    const ccomplex_t IT_0869 = (-8)*IT_0868;
    const ccomplex_t IT_0870 = (-8)*IT_0211;
    const ccomplex_t IT_0871 = IT_0867*IT_0870;
    const ccomplex_t IT_0872 = IT_0281*IT_0570;
    const ccomplex_t IT_0873 = IT_0317*IT_0572;
    const ccomplex_t IT_0874 = (-2)*IT_0832;
    const ccomplex_t IT_0875 = IT_0509*IT_0874;
    const ccomplex_t IT_0876 = (-2)*IT_0803;
    const ccomplex_t IT_0877 = IT_0488*IT_0876;
    const ccomplex_t IT_0878 = IT_0450*IT_0835;
    const ccomplex_t IT_0879 = IT_0473*IT_0843;
    const ccomplex_t IT_0880 = (-2)*IT_0515;
    const ccomplex_t IT_0881 = IT_0598*IT_0880;
    const ccomplex_t IT_0882 = IT_0845 + IT_0846 + IT_0847 + IT_0848 + IT_0849
       + IT_0850 + IT_0851 + IT_0852 + IT_0853 + IT_0856 + IT_0857 + IT_0862 +
       IT_0863 + IT_0864 + IT_0865 + IT_0869 + IT_0871 + IT_0872 + IT_0873 +
       IT_0875 + IT_0877 + IT_0878 + IT_0879 + IT_0881;
    const ccomplex_t IT_0883 = -IT_0265;
    const ccomplex_t IT_0884 = -IT_0243;
    const ccomplex_t IT_0885 = conj(IT_0255) + conj(IT_0883) + conj(IT_0884);
    const ccomplex_t IT_0886 = (-2)*IT_0815;
    const ccomplex_t IT_0887 = IT_0461*IT_0886;
    const ccomplex_t IT_0888 = IT_0845 + IT_0846 + IT_0847 + IT_0848 + IT_0849
       + IT_0850 + IT_0851 + IT_0852 + IT_0853 + IT_0856 + IT_0857 + IT_0862 +
       IT_0863 + IT_0864 + IT_0865 + IT_0869 + IT_0871 + IT_0872 + IT_0873 +
       IT_0875 + IT_0877 + IT_0878 + IT_0879 + IT_0881 + IT_0887;
    const ccomplex_t IT_0889 = IT_0000*IT_0004*IT_0041;
    const ccomplex_t IT_0890 = 0.5*IT_0889;
    const ccomplex_t IT_0891 = 1 + IT_0890;
    const ccomplex_t IT_0892 = IT_0037*IT_0891;
    const ccomplex_t IT_0893 = (-4)*IT_0892;
    const ccomplex_t IT_0894 = (-2)*conj(IT_0199);
    const ccomplex_t IT_0895 = (-2)*conj(IT_0211);
    const ccomplex_t IT_0896 = IT_0894 + IT_0895;
    const ccomplex_t IT_0897 = IT_0002 + IT_0533;
    const ccomplex_t IT_0898 = s_23*IT_0004*IT_0041;
    const ccomplex_t IT_0899 = (-2)*IT_0898;
    const ccomplex_t IT_0900 = IT_0897 + IT_0899;
    const ccomplex_t IT_0901 = IT_0000*IT_0900;
    const ccomplex_t IT_0902 = -IT_0532 + 0.5*IT_0901;
    const ccomplex_t IT_0903 = s_23 + IT_0902;
    const ccomplex_t IT_0904 = m_C_1*IT_0903;
    const ccomplex_t IT_0905 = 2*IT_0904;
    const ccomplex_t IT_0906 = 2*IT_0905;
    const ccomplex_t IT_0907 = m_C_1*IT_0345;
    const ccomplex_t IT_0908 = (-4)*IT_0907;
    const ccomplex_t IT_0909 = 2*IT_0908;
    const ccomplex_t IT_0910 = IT_0004*IT_0009;
    const ccomplex_t IT_0911 = s_23*s_34;
    const ccomplex_t IT_0912 = (-2)*IT_0911;
    const ccomplex_t IT_0913 = IT_0006 + IT_0912;
    const ccomplex_t IT_0914 = IT_0005*IT_0913;
    const ccomplex_t IT_0915 = (-0.25)*IT_0003 + (-0.25)*IT_0910 + (-0.25)
      *IT_0914;
    const ccomplex_t IT_0916 = s_23 + IT_0915;
    const ccomplex_t IT_0917 = m_C_1*IT_0916;
    const ccomplex_t IT_0918 = (-4)*IT_0917;
    const ccomplex_t IT_0919 = 2*IT_0918;
    const ccomplex_t IT_0920 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_0921 = (-2)*IT_0321 + -IT_0340 + 2*IT_0920;
    const ccomplex_t IT_0922 = IT_0320 + IT_0921;
    const ccomplex_t IT_0923 = IT_0004*IT_0922;
    const ccomplex_t IT_0924 = -IT_0283 + (-0.5)*IT_0923;
    const ccomplex_t IT_0925 = s_12 + IT_0924;
    const ccomplex_t IT_0926 = m_N_3*IT_0925;
    const ccomplex_t IT_0927 = (-2)*IT_0926;
    const ccomplex_t IT_0928 = 2*IT_0927;
    const ccomplex_t IT_0929 = IT_0000*IT_0369;
    const ccomplex_t IT_0930 = (-0.25)*IT_0929;
    const ccomplex_t IT_0931 = IT_0004*IT_0320;
    const ccomplex_t IT_0932 = (-0.25)*IT_0931;
    const ccomplex_t IT_0933 = s_12 + IT_0930 + IT_0932;
    const ccomplex_t IT_0934 = s_14*s_23;
    const ccomplex_t IT_0935 = s_12*s_34;
    const ccomplex_t IT_0936 = IT_0934 + IT_0935;
    const ccomplex_t IT_0937 = s_13*s_24;
    const ccomplex_t IT_0938 = -IT_0937;
    const ccomplex_t IT_0939 = IT_0936 + IT_0938;
    const ccomplex_t IT_0940 = IT_0005*IT_0939;
    const ccomplex_t IT_0941 = 0.25*IT_0940;
    const ccomplex_t IT_0942 = IT_0933 + IT_0941;
    const ccomplex_t IT_0943 = 8*IT_0942;
    const ccomplex_t IT_0944 = (-0.5)*IT_0324;
    const ccomplex_t IT_0945 = (-0.5)*IT_0700;
    const ccomplex_t IT_0946 = s_12 + IT_0944 + IT_0945;
    const ccomplex_t IT_0947 = (-2)*IT_0494 + (-2)*IT_0496 + 4*IT_0552;
    const ccomplex_t IT_0948 = IT_0493 + IT_0947;
    const ccomplex_t IT_0949 = IT_0437*IT_0948;
    const ccomplex_t IT_0950 = 0.25*IT_0949;
    const ccomplex_t IT_0951 = IT_0946 + IT_0950;
    const ccomplex_t IT_0952 = 8*IT_0951;
    const ccomplex_t IT_0953 = conj(IT_0339) + conj(IT_0854);
    const ccomplex_t IT_0954 = IT_0043*IT_0437;
    const ccomplex_t IT_0955 = 0.25*IT_0954;
    const ccomplex_t IT_0956 = 1 + IT_0955;
    const ccomplex_t IT_0957 = IT_0037*IT_0956;
    const ccomplex_t IT_0958 = 8*IT_0957;
    const ccomplex_t IT_0959 = conj(IT_0586) + conj(IT_0597);
    const ccomplex_t IT_0960 = IT_0040*IT_0437;
    const ccomplex_t IT_0961 = 0.125*IT_0960;
    const ccomplex_t IT_0962 = 1 + IT_0961;
    const ccomplex_t IT_0963 = IT_0037*IT_0962;
    const ccomplex_t IT_0964 = 16*IT_0963;
    const ccomplex_t IT_0965 = (-16)*IT_0963;
    const ccomplex_t IT_0966 = (-8)*IT_0942;
    const ccomplex_t IT_0967 = (-2)*IT_0494 + 2*IT_0495 + (-2)*IT_0496 + 2
      *IT_0552 + (-2)*IT_0553;
    const ccomplex_t IT_0968 = IT_0493 + IT_0967;
    const ccomplex_t IT_0969 = IT_0437*IT_0968;
    const ccomplex_t IT_0970 = (-0.25)*IT_0969;
    const ccomplex_t IT_0971 = -IT_0284;
    const ccomplex_t IT_0972 = -IT_0283;
    const ccomplex_t IT_0973 = s_12 + IT_0970 + IT_0971 + IT_0972;
    const ccomplex_t IT_0974 = (-8)*IT_0973;
    const ccomplex_t IT_0975 = (-8)*IT_0951;
    const ccomplex_t IT_0976 = (-8)*IT_0957;
    const ccomplex_t IT_0977 = (-2)*IT_0446;
    const ccomplex_t IT_0978 = conj(IT_0233) + conj(IT_0268);
    const ccomplex_t IT_0979 = IT_0000*IT_0008*IT_0269;
    const ccomplex_t IT_0980 = IT_0008*IT_0031;
    const ccomplex_t IT_0981 = -IT_0980;
    const ccomplex_t IT_0982 = IT_0979 + IT_0981;
    const ccomplex_t IT_0983 = IT_0004*IT_0982;
    const ccomplex_t IT_0984 = IT_0270 + IT_0711;
    const ccomplex_t IT_0985 = (-4)*IT_0984;
    const ccomplex_t IT_0986 = IT_0983 + IT_0985;
    const ccomplex_t IT_0987 = m_C_1*IT_0986;
    const ccomplex_t IT_0988 = (-2)*IT_0987;
    const ccomplex_t IT_0989 = m_N_3*IT_0527;
    const ccomplex_t IT_0990 = (-4)*IT_0989;
    const ccomplex_t IT_0991 = (-2)*IT_0990;
    const ccomplex_t IT_0992 = IT_0000*IT_0008*IT_0647;
    const ccomplex_t IT_0993 = s_13*s_14*s_34*IT_0000;
    const ccomplex_t IT_0994 = IT_0008*IT_0030;
    const ccomplex_t IT_0995 = 2*IT_0622 + (-2)*IT_0993 + -IT_0994;
    const ccomplex_t IT_0996 = IT_0992 + IT_0995;
    const ccomplex_t IT_0997 = IT_0004*IT_0996;
    const ccomplex_t IT_0998 = -IT_0030 + (-0.5)*IT_0997;
    const ccomplex_t IT_0999 = IT_0686 + IT_0998;
    const ccomplex_t IT_1000 = m_N_3*IT_0999;
    const ccomplex_t IT_1001 = (-2)*IT_1000;
    const ccomplex_t IT_1002 = (-2)*IT_1001;
    const ccomplex_t IT_1003 = IT_0004*IT_0055;
    const ccomplex_t IT_1004 = (-0.25)*IT_0052 + 0.25*IT_0054 + (-0.25)*IT_1003;
    const ccomplex_t IT_1005 = s_13 + IT_1004;
    const ccomplex_t IT_1006 = m_N_3*IT_1005;
    const ccomplex_t IT_1007 = (-4)*IT_1006;
    const ccomplex_t IT_1008 = (-2)*IT_1007;
    const ccomplex_t IT_1009 = 0.666666666666667*IT_0003 + (-0.166666666666667
      )*IT_0812;
    const ccomplex_t IT_1010 = s_23 + IT_1009;
    const ccomplex_t IT_1011 = m_C_1*IT_1010;
    const ccomplex_t IT_1012 = (-6)*IT_1011;
    const ccomplex_t IT_1013 = (-2)*IT_1012;
    const ccomplex_t IT_1014 = -IT_0893;
    const ccomplex_t IT_1015 = conj(IT_0281)*IT_0893 + IT_0116*IT_0896 + conj
      (IT_0368)*IT_0906 + conj(IT_0745)*IT_0909 + conj(IT_0545)*IT_0919 +
       IT_0796*IT_0928 + conj(IT_0036)*IT_0943 + conj(IT_0461)*IT_0952 + IT_0953
      *IT_0958 + IT_0959*IT_0964 + conj(IT_0450)*IT_0965 + conj(IT_0317)*IT_0966
       + conj(IT_0488)*IT_0974 + conj(IT_0181)*IT_0975 + conj(IT_0473)*IT_0976 +
       conj(IT_0509)*IT_0977 + IT_0978*IT_0988 + conj(IT_0109)*IT_0991 + conj
      (IT_0428)*IT_1002 + conj(IT_0050)*IT_1008 + conj(IT_0232)*IT_1013 + conj
      (IT_0097)*IT_1014;
    const ccomplex_t IT_1016 = IT_0317*IT_1014;
    const ccomplex_t IT_1017 = IT_0036*IT_0893;
    const ccomplex_t IT_1018 = (-10)*s_13;
    const ccomplex_t IT_1019 = IT_0462 + IT_1018;
    const ccomplex_t IT_1020 = m_N_3*IT_1019;
    const ccomplex_t IT_1021 = IT_0200*IT_1020;
    const ccomplex_t IT_1022 = IT_0212*IT_1020;
    const ccomplex_t IT_1023 = 8*IT_0973;
    const ccomplex_t IT_1024 = (-2)*IT_0116;
    const ccomplex_t IT_1025 = conj(IT_0317)*IT_1014;
    const ccomplex_t IT_1026 = conj(IT_0036)*IT_0893;
    const ccomplex_t IT_1027 = conj(IT_0601)*IT_0886 + IT_0446*IT_0896 + conj
      (IT_0109)*IT_0906 + conj(IT_0428)*IT_0909 + conj(IT_0050)*IT_0919 + conj
      (IT_0281)*IT_0943 + conj(IT_0450)*IT_0952 + conj(IT_0181)*IT_0964 + conj
      (IT_0461)*IT_0965 + conj(IT_0097)*IT_0966 + conj(IT_0473)*IT_0974 +
       IT_0959*IT_0975 + conj(IT_0488)*IT_0976 + IT_0928*IT_0978 + IT_0796
      *IT_0988 + conj(IT_0368)*IT_0991 + conj(IT_0745)*IT_1002 + conj(IT_0545)
      *IT_1008 + conj(IT_0598)*IT_1013 + IT_0953*IT_1023 + conj(IT_0509)*IT_1024
       + IT_1025 + IT_1026;
    const ccomplex_t IT_1028 = IT_0281*IT_1014;
    const ccomplex_t IT_1029 = IT_0097*IT_0893;
    const ccomplex_t IT_1030 = 2*conj(IT_0199);
    const ccomplex_t IT_1031 = 2*conj(IT_0211);
    const ccomplex_t IT_1032 = IT_1030 + IT_1031;
    const ccomplex_t IT_1033 = 2*IT_0446;
    const ccomplex_t IT_1034 = 2*IT_0987;
    const ccomplex_t IT_1035 = 2*IT_0990;
    const ccomplex_t IT_1036 = 2*IT_1001;
    const ccomplex_t IT_1037 = 2*IT_1007;
    const ccomplex_t IT_1038 = (-2)*IT_0905;
    const ccomplex_t IT_1039 = (-2)*IT_0908;
    const ccomplex_t IT_1040 = (-2)*IT_0918;
    const ccomplex_t IT_1041 = (-2)*IT_0927;
    const ccomplex_t IT_1042 = conj(IT_0281)*IT_1014;
    const ccomplex_t IT_1043 = conj(IT_0097)*IT_0893;
    const ccomplex_t IT_1044 = conj(IT_0317)*IT_0943 + conj(IT_0181)*IT_0952 +
       conj(IT_0473)*IT_0958 + conj(IT_0450)*IT_0964 + IT_0959*IT_0965 + conj
      (IT_0036)*IT_0966 + conj(IT_0461)*IT_0975 + IT_0953*IT_0976 + conj(IT_0488
      )*IT_1023 + IT_0116*IT_1032 + conj(IT_0509)*IT_1033 + IT_0978*IT_1034 +
       conj(IT_0109)*IT_1035 + conj(IT_0428)*IT_1036 + conj(IT_0050)*IT_1037 +
       conj(IT_0368)*IT_1038 + conj(IT_0745)*IT_1039 + conj(IT_0545)*IT_1040 +
       IT_0796*IT_1041 + IT_1042 + IT_1043;
    const ccomplex_t IT_1045 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_1046 = (-2)*IT_0321 + -IT_0340 + 2*IT_1045;
    const ccomplex_t IT_1047 = IT_0320 + IT_1046;
    const ccomplex_t IT_1048 = IT_0004*IT_1047;
    const ccomplex_t IT_1049 = -IT_0283 + (-0.5)*IT_1048;
    const ccomplex_t IT_1050 = s_12 + IT_1049;
    const ccomplex_t IT_1051 = m_C_1*IT_1050;
    const ccomplex_t IT_1052 = (-2)*IT_1051;
    const ccomplex_t IT_1053 = 2*IT_1052;
    const ccomplex_t IT_1054 = 2*IT_0347;
    const ccomplex_t IT_1055 = IT_0051 + IT_0522;
    const ccomplex_t IT_1056 = s_13*IT_0004*IT_0041;
    const ccomplex_t IT_1057 = (-2)*IT_1056;
    const ccomplex_t IT_1058 = IT_1055 + IT_1057;
    const ccomplex_t IT_1059 = IT_0000*IT_1058;
    const ccomplex_t IT_1060 = -IT_0436 + 0.5*IT_1059;
    const ccomplex_t IT_1061 = s_13 + IT_1060;
    const ccomplex_t IT_1062 = m_N_3*IT_1061;
    const ccomplex_t IT_1063 = 2*IT_1062;
    const ccomplex_t IT_1064 = 2*IT_1063;
    const ccomplex_t IT_1065 = s_13*s_34;
    const ccomplex_t IT_1066 = (-2)*IT_1065;
    const ccomplex_t IT_1067 = IT_0053 + IT_1066;
    const ccomplex_t IT_1068 = IT_0005*IT_1067;
    const ccomplex_t IT_1069 = (-0.25)*IT_0052 + (-0.25)*IT_1003 + (-0.25)
      *IT_1068;
    const ccomplex_t IT_1070 = s_13 + IT_1069;
    const ccomplex_t IT_1071 = m_N_3*IT_1070;
    const ccomplex_t IT_1072 = (-4)*IT_1071;
    const ccomplex_t IT_1073 = 2*IT_1072;
    const ccomplex_t IT_1074 = -IT_0935 + -IT_0937;
    const ccomplex_t IT_1075 = IT_0934 + IT_1074;
    const ccomplex_t IT_1076 = IT_0005*IT_1075;
    const ccomplex_t IT_1077 = (-0.25)*IT_1076;
    const ccomplex_t IT_1078 = IT_0933 + IT_1077;
    const ccomplex_t IT_1079 = 8*IT_1078;
    const ccomplex_t IT_1080 = IT_0437*IT_0498;
    const ccomplex_t IT_1081 = 0.25*IT_1080;
    const ccomplex_t IT_1082 = IT_0946 + IT_1081;
    const ccomplex_t IT_1083 = 8*IT_1082;
    const ccomplex_t IT_1084 = (-8)*IT_1078;
    const ccomplex_t IT_1085 = (-2)*IT_1020;
    const ccomplex_t IT_1086 = m_C_1*IT_0538;
    const ccomplex_t IT_1087 = (-4)*IT_1086;
    const ccomplex_t IT_1088 = (-2)*IT_1087;
    const ccomplex_t IT_1089 = s_23*s_24*s_34*IT_0000;
    const ccomplex_t IT_1090 = 2*IT_0271 + -IT_0980 + (-2)*IT_1089;
    const ccomplex_t IT_1091 = IT_0979 + IT_1090;
    const ccomplex_t IT_1092 = IT_0004*IT_1091;
    const ccomplex_t IT_1093 = -IT_0031 + (-0.5)*IT_1092;
    const ccomplex_t IT_1094 = IT_0270 + IT_1093;
    const ccomplex_t IT_1095 = m_C_1*IT_1094;
    const ccomplex_t IT_1096 = (-2)*IT_1095;
    const ccomplex_t IT_1097 = (-2)*IT_1096;
    const ccomplex_t IT_1098 = (-0.25)*IT_0003 + 0.25*IT_0007 + (-0.25)*IT_0910;
    const ccomplex_t IT_1099 = s_23 + IT_1098;
    const ccomplex_t IT_1100 = m_C_1*IT_1099;
    const ccomplex_t IT_1101 = (-4)*IT_1100;
    const ccomplex_t IT_1102 = (-2)*IT_1101;
    const ccomplex_t IT_1103 = -IT_0994;
    const ccomplex_t IT_1104 = IT_0992 + IT_1103;
    const ccomplex_t IT_1105 = IT_0004*IT_1104;
    const ccomplex_t IT_1106 = IT_0686 + IT_0732;
    const ccomplex_t IT_1107 = (-4)*IT_1106;
    const ccomplex_t IT_1108 = IT_1105 + IT_1107;
    const ccomplex_t IT_1109 = m_N_3*IT_1108;
    const ccomplex_t IT_1110 = (-2)*IT_1109;
    const ccomplex_t IT_1111 = IT_0481*IT_0896 + IT_0958*IT_0959 + IT_0953
      *IT_0964 + conj(IT_0473)*IT_0965 + conj(IT_0461)*IT_0974 + conj(IT_0450)
      *IT_0976 + conj(IT_0181)*IT_1023 + IT_1042 + IT_1043 + conj(IT_0745)
      *IT_1053 + IT_0796*IT_1054 + conj(IT_0109)*IT_1064 + conj(IT_0050)*IT_1073
       + conj(IT_0317)*IT_1079 + conj(IT_0488)*IT_1083 + conj(IT_0036)*IT_1084 +
       conj(IT_0509)*IT_1085 + conj(IT_0368)*IT_1088 + IT_0978*IT_1097 + conj
      (IT_0545)*IT_1102 + conj(IT_0428)*IT_1110;
    const ccomplex_t IT_1112 = IT_0036*IT_1014;
    const ccomplex_t IT_1113 = IT_0317*IT_0893;
    const ccomplex_t IT_1114 = IT_0200*IT_0446;
    const ccomplex_t IT_1115 = IT_0212*IT_0446;
    const ccomplex_t IT_1116 = IT_0795*IT_1034;
    const ccomplex_t IT_1117 = IT_0097*IT_0943;
    const ccomplex_t IT_1118 = IT_0586*IT_0952;
    const ccomplex_t IT_1119 = IT_0854*IT_0974;
    const ccomplex_t IT_1120 = IT_0339*IT_0974;
    const ccomplex_t IT_1121 = IT_0268*IT_1041;
    const ccomplex_t IT_1122 = (-8)*IT_1082;
    const ccomplex_t IT_1123 = (-2)*IT_0481;
    const ccomplex_t IT_1124 = conj(IT_0036)*IT_1014;
    const ccomplex_t IT_1125 = conj(IT_0317)*IT_0893;
    const ccomplex_t IT_1126 = IT_0799*IT_0876 + IT_0466*IT_0896 + conj
      (IT_0181)*IT_0958 + conj(IT_0488)*IT_0965 + conj(IT_0450)*IT_0974 + conj
      (IT_0461)*IT_0976 + IT_0959*IT_1023 + conj(IT_0428)*IT_1053 + IT_0978
      *IT_1054 + conj(IT_0368)*IT_1064 + conj(IT_0545)*IT_1073 + conj(IT_0097)
      *IT_1079 + conj(IT_0473)*IT_1083 + conj(IT_0281)*IT_1084 + conj(IT_0109)
      *IT_1088 + IT_0796*IT_1097 + conj(IT_0050)*IT_1102 + conj(IT_0745)*IT_1110
       + IT_0953*IT_1122 + conj(IT_0509)*IT_1123 + IT_1124 + IT_1125;
    const ccomplex_t IT_1127 = 4*IT_0211;
    const ccomplex_t IT_1128 = 0.25*IT_1127;
    const ccomplex_t IT_1129 = IT_0281*IT_0793;
    const ccomplex_t IT_1130 = IT_0619 + IT_0858;
    const ccomplex_t IT_1131 = s_12*s_13*IT_0008;
    const ccomplex_t IT_1132 = s_23*IT_0008*IT_0030;
    const ccomplex_t IT_1133 = s_13*s_14*s_24;
    const ccomplex_t IT_1134 = s_24*s_34*IT_0030;
    const ccomplex_t IT_1135 = -IT_1132 + (-2)*IT_1133 + 2*IT_1134;
    const ccomplex_t IT_1136 = IT_1131 + IT_1135;
    const ccomplex_t IT_1137 = IT_0004*IT_1136;
    const ccomplex_t IT_1138 = (-0.5)*IT_1137;
    const ccomplex_t IT_1139 = IT_1130 + IT_1138;
    const ccomplex_t IT_1140 = s_14*s_34*IT_0000*IT_0004*IT_0269;
    const ccomplex_t IT_1141 = s_12*s_23*IT_0000*IT_0004*IT_0041;
    const ccomplex_t IT_1142 = -IT_1141;
    const ccomplex_t IT_1143 = IT_0351 + IT_0602 + IT_0678 + IT_0683 + IT_1140
       + IT_1142;
    const ccomplex_t IT_1144 = s_12 + IT_0972;
    const ccomplex_t IT_1145 = IT_0037*IT_1144;
    const ccomplex_t IT_1146 = (-6)*IT_1145;
    const ccomplex_t IT_1147 = conj(IT_0461)*IT_0988;
    const ccomplex_t IT_1148 = 0.5*IT_1147;
    const ccomplex_t IT_1149 = s_14*s_34*IT_0004*IT_0031;
    const ccomplex_t IT_1150 = -IT_1149;
    const ccomplex_t IT_1151 = IT_0351 + IT_0678 + IT_0681 + IT_1140 + IT_1142
       + IT_1150;
    const ccomplex_t IT_1152 = 2*IT_1151;
    const ccomplex_t IT_1153 = conj(IT_0368)*IT_1152;
    const ccomplex_t IT_1154 = 0.5*IT_1153;
    const ccomplex_t IT_1155 = (s_12 + -IT_0283)*(s_12 + (-0.5)*IT_0324);
    const ccomplex_t IT_1156 = (-2)*IT_1155;
    const ccomplex_t IT_1157 = 2*IT_1156;
    const ccomplex_t IT_1158 = conj(IT_0745)*IT_1157;
    const ccomplex_t IT_1159 = 0.5*IT_1158;
    const ccomplex_t IT_1160 = IT_0037*IT_0984;
    const ccomplex_t IT_1161 = 6*IT_1160;
    const ccomplex_t IT_1162 = IT_0796*IT_1161;
    const ccomplex_t IT_1163 = 0.5*IT_1162;
    const ccomplex_t IT_1164 = conj(IT_0473)*IT_1054;
    const ccomplex_t IT_1165 = 0.5*IT_1164;
    const ccomplex_t IT_1166 = conj(IT_0488)*IT_1097;
    const ccomplex_t IT_1167 = 0.5*IT_1166;
    const ccomplex_t IT_1168 = 0.5*conj(IT_0181);
    const ccomplex_t IT_1169 = IT_1034*IT_1168;
    const ccomplex_t IT_1170 = conj(IT_0036)*IT_0278;
    const ccomplex_t IT_1171 = 0.5*IT_1170;
    const ccomplex_t IT_1172 = conj(IT_0281)*IT_0288;
    const ccomplex_t IT_1173 = 0.5*IT_1172;
    const ccomplex_t IT_1174 = IT_0037*IT_0276;
    const ccomplex_t IT_1175 = 2*IT_1174;
    const ccomplex_t IT_1176 = conj(IT_0575)*IT_1175;
    const ccomplex_t IT_1177 = 0.5*IT_1176;
    const ccomplex_t IT_1178 = conj(IT_0050)*IT_0295;
    const ccomplex_t IT_1179 = 0.5*IT_1178;
    const ccomplex_t IT_1180 = conj(IT_0109)*IT_0295;
    const ccomplex_t IT_1181 = 0.5*IT_1180;
    const ccomplex_t IT_1182 = conj(IT_0450)*IT_0928;
    const ccomplex_t IT_1183 = 0.5*IT_1182;
    const ccomplex_t IT_1184 = (-4)*conj(IT_0199);
    const ccomplex_t IT_1185 = IT_0362*IT_1184;
    const ccomplex_t IT_1186 = 0.5*IT_1185;
    const ccomplex_t IT_1187 = IT_0362*IT_0508;
    const ccomplex_t IT_1188 = 0.5*IT_1187;
    const ccomplex_t IT_1189 = conj(IT_0317)*IT_0318;
    const ccomplex_t IT_1190 = 0.5*IT_1189;
    const ccomplex_t IT_1191 = conj(IT_0097)*IT_0797;
    const ccomplex_t IT_1192 = 0.5*IT_1191;
    const ccomplex_t IT_1193 = IT_0037*IT_0286;
    const ccomplex_t IT_1194 = (-2)*IT_1193;
    const ccomplex_t IT_1195 = conj(IT_0575)*IT_1194;
    const ccomplex_t IT_1196 = 0.5*IT_1195;
    const ccomplex_t IT_1197 = IT_0348*conj(IT_0854);
    const ccomplex_t IT_1198 = 0.5*IT_1197;
    const ccomplex_t IT_1199 = conj(IT_0339)*IT_0348;
    const ccomplex_t IT_1200 = 0.5*IT_1199;
    const ccomplex_t IT_1201 = 0.5*IT_0327*IT_0978 + conj(IT_0545)*IT_1143 +
       0.5*conj(IT_0428)*IT_1146 + IT_1148 + IT_1154 + IT_1159 + IT_1163 +
       IT_1165 + IT_1167 + IT_1169 + IT_1171 + IT_1173 + IT_1177 + IT_1179 +
       IT_1181 + IT_1183 + IT_1186 + IT_1188 + IT_1190 + IT_1192 + IT_1196 +
       IT_1198 + IT_1200;
    const ccomplex_t IT_1202 = 0.5*IT_0959;
    const ccomplex_t IT_1203 = conj(IT_0545)*IT_1143 + 0.5*conj(IT_0428)
      *IT_1146 + IT_1148 + IT_1154 + IT_1159 + IT_1163 + IT_1165 + IT_1167 +
       IT_1169 + IT_1171 + IT_1173 + IT_1177 + IT_1179 + IT_1181 + IT_1183 +
       IT_1186 + IT_1188 + IT_1190 + IT_1192 + IT_1196 + IT_1198 + IT_1200 +
       IT_1041*IT_1202;
    const ccomplex_t IT_1204 = s_12*IT_0891;
    const ccomplex_t IT_1205 = 4*IT_0892;
    const ccomplex_t IT_1206 = (-0.25)*conj(IT_0036);
    const ccomplex_t IT_1207 = s_23*m_C_1;
    const ccomplex_t IT_1208 = IT_0891*IT_1207;
    const ccomplex_t IT_1209 = 4*IT_1208;
    const ccomplex_t IT_1210 = (-0.25)*conj(IT_0050);
    const ccomplex_t IT_1211 = (-0.25)*conj(IT_0450);
    const ccomplex_t IT_1212 = 4*IT_1204;
    const ccomplex_t IT_1213 = (-0.25)*conj(IT_0097);
    const ccomplex_t IT_1214 = s_34*m_C_1*IT_0004;
    const ccomplex_t IT_1215 = IT_0293*IT_1214;
    const ccomplex_t IT_1216 = (-2)*IT_1215;
    const ccomplex_t IT_1217 = (-0.25)*conj(IT_0109);
    const ccomplex_t IT_1218 = m_C_1*IT_0695;
    const ccomplex_t IT_1219 = 2*IT_1218;
    const ccomplex_t IT_1220 = (-0.25)*conj(IT_0428);
    const ccomplex_t IT_1221 = (-0.25)*IT_0953;
    const ccomplex_t IT_1222 = (-0.25)*conj(IT_0461);
    const ccomplex_t IT_1223 = s_34*m_N_3*IT_0004;
    const ccomplex_t IT_1224 = IT_0431*IT_1223;
    const ccomplex_t IT_1225 = 2*IT_1224;
    const ccomplex_t IT_1226 = (-0.25)*conj(IT_0368);
    const ccomplex_t IT_1227 = (-0.25)*conj(IT_0745);
    const ccomplex_t IT_1228 = (-0.25)*IT_0978;
    const ccomplex_t IT_1229 = s_13*m_N_3;
    const ccomplex_t IT_1230 = IT_0891*IT_1229;
    const ccomplex_t IT_1231 = (-4)*IT_1230;
    const ccomplex_t IT_1232 = (-0.25)*conj(IT_0545);
    const ccomplex_t IT_1233 = (-0.25)*conj(IT_0473);
    const ccomplex_t IT_1234 = conj(IT_0317)*IT_0892 + (-0.25)*IT_0063*IT_1032
       + conj(IT_0281)*IT_1204 + IT_1205*IT_1206 + IT_1209*IT_1210 + IT_0966
      *IT_1211 + IT_1212*IT_1213 + IT_1216*IT_1217 + IT_1219*IT_1220 + IT_1084
      *IT_1221 + IT_1014*IT_1222 + IT_1225*IT_1226 + IT_0693*IT_1227 + IT_0797
      *IT_1228 + IT_1231*IT_1232 + IT_1079*IT_1233;
    const ccomplex_t IT_1235 = (-4)*IT_1208;
    const ccomplex_t IT_1236 = (-0.25)*conj(IT_0281);
    const ccomplex_t IT_1237 = (-0.25)*conj(IT_0317);
    const ccomplex_t IT_1238 = 2*IT_1215;
    const ccomplex_t IT_1239 = (-2)*IT_1218;
    const ccomplex_t IT_1240 = (-2)*IT_1224;
    const ccomplex_t IT_1241 = 4*IT_1230;
    const ccomplex_t IT_1242 = (-0.25)*IT_0959;
    const ccomplex_t IT_1243 = conj(IT_0036)*IT_0892 + (-0.25)*conj(IT_0461)
      *IT_0893 + (-0.25)*IT_0063*IT_0896 + conj(IT_0097)*IT_1204 + IT_0943
      *IT_1211 + IT_1079*IT_1221 + IT_1084*IT_1233 + IT_1210*IT_1235 + IT_1212
      *IT_1236 + IT_1205*IT_1237 + IT_1217*IT_1238 + IT_1220*IT_1239 + IT_1226
      *IT_1240 + IT_1232*IT_1241 + IT_0966*IT_1242;
    const ccomplex_t IT_1244 = -IT_0370;
    const ccomplex_t IT_1245 = IT_0369 + IT_0372 + IT_0382 + IT_0383 + IT_1244;
    const ccomplex_t IT_1246 = s_34*IT_0004;
    const ccomplex_t IT_1247 = IT_0934 + IT_0937;
    const ccomplex_t IT_1248 = s_13*s_23*s_34*IT_0000;
    const ccomplex_t IT_1249 = (-2)*IT_1248;
    const ccomplex_t IT_1250 = IT_1247 + IT_1249;
    const ccomplex_t IT_1251 = IT_1246*IT_1250;
    const ccomplex_t IT_1252 = (-2)*IT_1251;
    const ccomplex_t IT_1253 = (IT_0001 + -IT_0098)*(s_12 + (-0.5)*IT_0700);
    const ccomplex_t IT_1254 = (-2)*IT_1253;
    const ccomplex_t IT_1255 = 2*IT_1254;
    const ccomplex_t IT_1256 = s_24*s_34*IT_0000*IT_0004*IT_0647;
    const ccomplex_t IT_1257 = s_24*s_34*IT_0004*IT_0030;
    const ccomplex_t IT_1258 = -IT_1257;
    const ccomplex_t IT_1259 = s_12*s_13*IT_0000*IT_0004*IT_0041;
    const ccomplex_t IT_1260 = -IT_1259;
    const ccomplex_t IT_1261 = IT_0632 + IT_0670 + IT_0673 + IT_1256 + IT_1258
       + IT_1260;
    const ccomplex_t IT_1262 = 2*IT_1261;
    const ccomplex_t IT_1263 = 2*IT_1087;
    const ccomplex_t IT_1264 = -IT_0098;
    const ccomplex_t IT_1265 = IT_0001 + IT_1264;
    const ccomplex_t IT_1266 = IT_0037*IT_1265;
    const ccomplex_t IT_1267 = (-6)*IT_1266;
    const ccomplex_t IT_1268 = IT_0433*IT_1227;
    const ccomplex_t IT_1269 = (-0.25)*IT_0796;
    const ccomplex_t IT_1270 = (-0.25)*conj(IT_0488);
    const ccomplex_t IT_1271 = (-0.25)*conj(IT_0181);
    const ccomplex_t IT_1272 = IT_0385*conj(IT_0509) + (-0.25)*IT_0107*IT_1032
       + IT_0906*IT_1211 + IT_1213*IT_1216 + IT_0991*IT_1222 + IT_1088*IT_1233 +
       IT_1225*IT_1237 + IT_1236*IT_1238 + IT_1206*IT_1240 + IT_1038*IT_1242 +
       IT_0799*IT_1245 + IT_1210*IT_1252 + IT_1217*IT_1255 + IT_1220*IT_1262 +
       IT_1221*IT_1263 + IT_1226*IT_1267 + IT_1268 + IT_1152*IT_1269 + IT_1064
      *IT_1270 + IT_1035*IT_1271;
    const ccomplex_t IT_1273 = (IT_0369 + (-2)*IT_0492)*(1 + 0.5*IT_0889);
    const ccomplex_t IT_1274 = IT_0620 + IT_0632 + IT_0670 + IT_0676 + IT_1256
       + IT_1260;
    const ccomplex_t IT_1275 = 2*IT_1101;
    const ccomplex_t IT_1276 = m_C_1*m_N_3*IT_0001;
    const ccomplex_t IT_1277 = IT_0891*IT_1276;
    const ccomplex_t IT_1278 = (-4)*IT_1277;
    const ccomplex_t IT_1279 = (-0.25)*IT_0048*IT_1032 + (-0.5)*IT_0796
      *IT_1143 + IT_0919*IT_1211 + IT_1209*IT_1213 + IT_1008*IT_1222 + IT_0295
      *IT_1228 + IT_1102*IT_1233 + IT_1235*IT_1236 + IT_1231*IT_1237 + IT_1206
      *IT_1241 + IT_1040*IT_1242 + IT_1217*IT_1252 + IT_1268 + IT_1073*IT_1270 +
       IT_1037*IT_1271 + conj(IT_0050)*IT_1273 + (-0.5)*conj(IT_0428)*IT_1274 +
       IT_1221*IT_1275 + IT_1232*IT_1278;
    const ccomplex_t IT_1280 = (-0.25)*IT_0017;
    const ccomplex_t IT_1281 = conj(IT_0097)*IT_0892 + conj(IT_0036)*IT_1204 +
       IT_0693*IT_1220 + IT_0966*IT_1222 + IT_1217*IT_1225 + IT_1216*IT_1226 +
       IT_1219*IT_1227 + IT_1210*IT_1231 + IT_1209*IT_1232 + IT_1205*IT_1236 +
       IT_1212*IT_1237 + IT_1079*IT_1270 + IT_0943*IT_1271 + IT_0896*IT_1280;
    const ccomplex_t IT_1282 = conj(IT_0281)*IT_0892 + conj(IT_0317)*IT_1204 +
       IT_1206*IT_1212 + IT_1205*IT_1213 + IT_0793*IT_1220 + IT_0943*IT_1222 +
       IT_1232*IT_1235 + IT_1226*IT_1238 + IT_1227*IT_1239 + IT_1217*IT_1240 +
       IT_1210*IT_1241 + IT_0288*IT_1269 + IT_1084*IT_1270 + IT_0966*IT_1271 +
       IT_1032*IT_1280;
    const ccomplex_t IT_1283 = s_23*IT_0622;
    const ccomplex_t IT_1284 = s_12*s_14*s_34;
    const ccomplex_t IT_1285 = -IT_1132 + 2*IT_1283 + (-2)*IT_1284;
    const ccomplex_t IT_1286 = IT_1131 + IT_1285;
    const ccomplex_t IT_1287 = IT_0004*IT_1286;
    const ccomplex_t IT_1288 = (-0.5)*IT_1287;
    const ccomplex_t IT_1289 = IT_1130 + IT_1288;
    const ccomplex_t IT_1290 = conj(IT_0232) + conj(IT_0255) + conj(IT_0883) +
       conj(IT_0884);
    const ccomplex_t IT_1291 = 4*conj(IT_0199);
    const ccomplex_t IT_1292 = 4*conj(IT_0211);
    const ccomplex_t IT_1293 = IT_1291 + IT_1292;
    const ccomplex_t IT_1294 = conj(IT_0097)*IT_0693;
    const ccomplex_t IT_1295 = conj(IT_0281)*IT_0793;
    const ccomplex_t IT_1296 = IT_0037*IT_1106;
    const ccomplex_t IT_1297 = 6*IT_1296;
    const ccomplex_t IT_1298 = 2*IT_1109;
    const ccomplex_t IT_1299 = -(s_12 + (-0.5)*IT_0324)*(IT_0030 + -IT_0686);
    const ccomplex_t IT_1300 = (-2)*IT_1299;
    const ccomplex_t IT_1301 = (-2)*IT_1300;
    const ccomplex_t IT_1302 = 0.5*(conj(IT_0050) + conj(IT_0109))*IT_0433 + (
      -2)*IT_1002*IT_1211 + (-2)*IT_0909*IT_1222 + (-2)*IT_1157*IT_1228 + (-2)
      *IT_1110*IT_1233 + (-2)*IT_1219*IT_1237 + (-2)*IT_1206*IT_1239 + (-2)
      *IT_1036*IT_1242 + (-2)*IT_1226*IT_1262 + (-2)*IT_1146*IT_1269 + (-2)
      *IT_1053*IT_1270 + (-2)*IT_1039*IT_1271 + conj(IT_0545)*IT_1274 + (-2)
      *IT_1289*IT_1290 + 0.5*IT_1139*IT_1293 + 0.5*IT_1294 + 0.5*IT_1295 + (-2)
      *IT_1220*IT_1297 + (-2)*IT_1221*IT_1298 + (-2)*IT_1227*IT_1301;
    const ccomplex_t IT_1303 = conj(IT_0317)*IT_0693;
    const ccomplex_t IT_1304 = conj(IT_0036)*IT_0793;
    const ccomplex_t IT_1305 = (-2)*IT_1052;
    const ccomplex_t IT_1306 = 4*IT_1139;
    const ccomplex_t IT_1307 = (-0.25)*conj(IT_0509);
    const ccomplex_t IT_1308 = (-2)*IT_0909*IT_1211 + (-2)*IT_1213*IT_1219 + (
      -2)*IT_1002*IT_1222 + (-2)*IT_1146*IT_1228 + (-2)*IT_1053*IT_1233 + (-2)
      *IT_1236*IT_1239 + (-2)*IT_1039*IT_1242 + (-2)*IT_1217*IT_1262 + (-2)
      *IT_1157*IT_1269 + (-2)*IT_1110*IT_1270 + (-2)*IT_1036*IT_1271 + conj
      (IT_0050)*IT_1274 + (-2)*IT_0799*IT_1289 + (-2)*IT_1227*IT_1297 + (-2)
      *IT_1220*IT_1301 + 0.5*IT_1303 + 0.5*IT_1304 + (-2)*IT_1221*IT_1305 + (-2)
      *IT_1306*IT_1307;
    const ccomplex_t IT_1309 = IT_0508 + IT_1184;
    const ccomplex_t IT_1310 = (-2)*IT_1063;
    const ccomplex_t IT_1311 = IT_0295*conj(IT_0795);
    const ccomplex_t IT_1312 = (-0.25)*IT_1311;
    const ccomplex_t IT_1313 = IT_0295*conj(IT_0409);
    const ccomplex_t IT_1314 = (-0.25)*IT_1313;
    const ccomplex_t IT_1315 = conj(IT_0428)*IT_0433;
    const ccomplex_t IT_1316 = (-0.25)*IT_1315;
    const ccomplex_t IT_1317 = IT_0991*IT_1211 + IT_0906*IT_1222 + IT_1213
      *IT_1225 + IT_1152*IT_1228 + IT_1064*IT_1233 + IT_1216*IT_1237 + IT_1206
      *IT_1238 + IT_1236*IT_1240 + IT_1035*IT_1242 + IT_1232*IT_1252 + IT_1226
      *IT_1255 + IT_1227*IT_1262 + IT_1217*IT_1267 + IT_1088*IT_1270 + IT_1038
      *IT_1271 + IT_1245*IT_1290 + (-0.25)*IT_0385*IT_1309 + IT_1221*IT_1310 +
       IT_1312 + IT_1314 + IT_1316;
    const ccomplex_t IT_1318 = (-2)*IT_1072;
    const ccomplex_t IT_1319 = (-0.5)*IT_0978*IT_1143 + IT_1008*IT_1211 +
       IT_0919*IT_1222 + IT_1213*IT_1231 + IT_1073*IT_1233 + IT_1206*IT_1235 +
       IT_1209*IT_1237 + IT_1236*IT_1241 + IT_1037*IT_1242 + IT_1226*IT_1252 +
       IT_1102*IT_1270 + IT_1040*IT_1271 + conj(IT_0545)*IT_1273 + (-0.5)*conj
      (IT_0745)*IT_1274 + IT_1210*IT_1278 + (-0.25)*IT_0567*IT_1293 + IT_1312 +
       IT_1314 + IT_1316 + IT_1221*IT_1318;
    const ccomplex_t IT_1320 = 2*IT_0481;
    const ccomplex_t IT_1321 = 2*IT_0116;
    const ccomplex_t IT_1322 = (-8)*IT_0867;
    const ccomplex_t IT_1323 = (-0.25)*IT_0799;
    const ccomplex_t IT_1324 = conj(IT_0317)*IT_0516;
    const ccomplex_t IT_1325 = (-0.25)*IT_1324;
    const ccomplex_t IT_1326 = conj(IT_0281)*IT_0518;
    const ccomplex_t IT_1327 = (-0.25)*IT_1326;
    const ccomplex_t IT_1328 = conj(IT_0521)*IT_0530;
    const ccomplex_t IT_1329 = (-0.25)*IT_1328;
    const ccomplex_t IT_1330 = conj(IT_0521)*IT_0541;
    const ccomplex_t IT_1331 = (-0.25)*IT_1330;
    const ccomplex_t IT_1332 = conj(IT_0368)*IT_0543;
    const ccomplex_t IT_1333 = (-0.25)*IT_1332;
    const ccomplex_t IT_1334 = conj(IT_0545)*IT_0546;
    const ccomplex_t IT_1335 = (-0.25)*IT_1334;
    const ccomplex_t IT_1336 = conj(IT_0050)*IT_0568;
    const ccomplex_t IT_1337 = (-0.25)*IT_1336;
    const ccomplex_t IT_1338 = conj(IT_0036)*IT_0570;
    const ccomplex_t IT_1339 = (-0.25)*IT_1338;
    const ccomplex_t IT_1340 = conj(IT_0097)*IT_0572;
    const ccomplex_t IT_1341 = (-0.25)*IT_1340;
    const ccomplex_t IT_1342 = conj(IT_0109)*IT_0385 + IT_0507*conj(IT_0509) +
       IT_0362*IT_0796 + (-0.25)*IT_0515*IT_0896 + IT_1024*IT_1211 + IT_0977
      *IT_1222 + IT_1123*IT_1233 + IT_1085*IT_1270 + IT_1033*IT_1271 + IT_1220
      *IT_1306 + IT_1221*IT_1320 + IT_1242*IT_1321 + IT_1322*IT_1323 + IT_1325 +
       IT_1327 + IT_1329 + IT_1331 + IT_1333 + IT_1335 + IT_1337 + IT_1339 +
       IT_1341;
    const ccomplex_t IT_1343 = IT_0278*conj(IT_0281);
    const ccomplex_t IT_1344 = conj(IT_0097)*IT_0318;
    const ccomplex_t IT_1345 = 2*IT_0758;
    const ccomplex_t IT_1346 = (-0.25)*IT_1345;
    const ccomplex_t IT_1347 = IT_0584 + IT_1346;
    const ccomplex_t IT_1348 = (-0.125)*conj(IT_0036);
    const ccomplex_t IT_1349 = (s_12 + -IT_0284)*(s_12 + (-0.5)*IT_0700);
    const ccomplex_t IT_1350 = (-2)*IT_1349;
    const ccomplex_t IT_1351 = 4*IT_1350;
    const ccomplex_t IT_1352 = -(IT_0030 + -IT_0687)*(s_12 + (-0.5)*IT_0700);
    const ccomplex_t IT_1353 = (-2)*IT_1352;
    const ccomplex_t IT_1354 = (-2)*IT_1353;
    const ccomplex_t IT_1355 = -(IT_0031 + -IT_0272)*(s_12 + (-0.5)*IT_0700);
    const ccomplex_t IT_1356 = (-2)*IT_1355;
    const ccomplex_t IT_1357 = (-2)*IT_1356;
    const ccomplex_t IT_1358 = IT_1351 + IT_1354 + IT_1357;
    const ccomplex_t IT_1359 = (-0.125)*conj(IT_0575);
    const ccomplex_t IT_1360 = (-0.125)*conj(IT_0450);
    const ccomplex_t IT_1361 = (-0.125)*conj(IT_0281);
    const ccomplex_t IT_1362 = (-0.125)*conj(IT_0097);
    const ccomplex_t IT_1363 = (-0.125)*conj(IT_0317);
    const ccomplex_t IT_1364 = (-0.125)*conj(IT_0109);
    const ccomplex_t IT_1365 = (-0.125)*conj(IT_0428);
    const ccomplex_t IT_1366 = 2*IT_0724;
    const ccomplex_t IT_1367 = (-2)*IT_0715;
    const ccomplex_t IT_1368 = (-0.125)*IT_0953;
    const ccomplex_t IT_1369 = (-0.125)*conj(IT_0461);
    const ccomplex_t IT_1370 = IT_0037*IT_0733;
    const ccomplex_t IT_1371 = 6*IT_1370;
    const ccomplex_t IT_1372 = IT_0037*IT_0712;
    const ccomplex_t IT_1373 = 6*IT_1372;
    const ccomplex_t IT_1374 = s_12 + IT_0971;
    const ccomplex_t IT_1375 = IT_0037*IT_1374;
    const ccomplex_t IT_1376 = (-12)*IT_1375;
    const ccomplex_t IT_1377 = IT_1371 + IT_1373 + IT_1376;
    const ccomplex_t IT_1378 = (-0.125)*conj(IT_0521);
    const ccomplex_t IT_1379 = (-0.125)*conj(IT_0368);
    const ccomplex_t IT_1380 = (-0.125)*conj(IT_0745);
    const ccomplex_t IT_1381 = IT_1175 + IT_1194;
    const ccomplex_t IT_1382 = (-0.125)*IT_0978;
    const ccomplex_t IT_1383 = IT_0000*IT_0004*IT_0269*IT_0622;
    const ccomplex_t IT_1384 = IT_0004*IT_0031*IT_0622;
    const ccomplex_t IT_1385 = -IT_1384;
    const ccomplex_t IT_1386 = IT_0000*IT_0030*IT_0269;
    const ccomplex_t IT_1387 = -IT_1386;
    const ccomplex_t IT_1388 = s_12*s_14*s_23*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_1389 = -IT_1388;
    const ccomplex_t IT_1390 = IT_0658 + IT_0659 + IT_1383 + IT_1385 + IT_1387
       + IT_1389;
    const ccomplex_t IT_1391 = 2*IT_1390;
    const ccomplex_t IT_1392 = s_14*s_24*IT_0000*IT_0004*IT_0269;
    const ccomplex_t IT_1393 = s_13*s_23*IT_0000*IT_0004*IT_0271;
    const ccomplex_t IT_1394 = s_12*IT_0031;
    const ccomplex_t IT_1395 = s_14*s_24*IT_0004*IT_0031;
    const ccomplex_t IT_1396 = -IT_1395;
    const ccomplex_t IT_1397 = s_13*s_23*IT_0000*IT_0031;
    const ccomplex_t IT_1398 = -IT_1397;
    const ccomplex_t IT_1399 = s_12*s_23*s_24*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_1400 = -IT_1399;
    const ccomplex_t IT_1401 = IT_1392 + IT_1393 + IT_1394 + IT_1396 + IT_1398
       + IT_1400;
    const ccomplex_t IT_1402 = 2*IT_1401;
    const ccomplex_t IT_1403 = IT_1391 + IT_1402;
    const ccomplex_t IT_1404 = (-0.125)*IT_0796;
    const ccomplex_t IT_1405 = (-0.125)*conj(IT_0473);
    const ccomplex_t IT_1406 = (-0.125)*conj(IT_0488);
    const ccomplex_t IT_1407 = conj(IT_0509)*IT_0644 + (-0.125)*conj(IT_0050)
      *IT_0685 + (-0.125)*conj(IT_0181)*IT_0744 + IT_0799*(IT_0618 + (-0.125)
      *IT_0861) + (-0.125)*IT_0645*IT_1032 + (-0.125)*IT_1303 + (-0.125)*IT_1304
       + (-0.125)*IT_1343 + (-0.125)*IT_1344 + IT_1202*IT_1347 + IT_0792*IT_1348
       + IT_1358*IT_1359 + IT_0760*IT_1360 + IT_0790*IT_1361 + IT_0788*IT_1362 +
       IT_0697*IT_1363 + IT_0707*IT_1364 + IT_0669*IT_1365 + (IT_1366 + IT_1367)
      *IT_1368 + IT_0741*IT_1369 + IT_1377*IT_1378 + IT_0786*IT_1379 + IT_0750
      *IT_1380 + IT_1381*IT_1382 + IT_1403*IT_1404 + IT_0726*IT_1405 + IT_0777
      *IT_1406;
    const ccomplex_t IT_1408 = 2*IT_0775;
    const ccomplex_t IT_1409 = (-4)*IT_0767;
    const ccomplex_t IT_1410 = IT_1408 + IT_1409;
    const ccomplex_t IT_1411 = (-0.125)*IT_0642;
    const ccomplex_t IT_1412 = (-0.125)*IT_0279 + (-0.125)*IT_0319 + (-0.125)
      *IT_0545*IT_0685 + (-0.125)*IT_0097*IT_0698 + (-0.125)*IT_0368*IT_0707 + (
      -0.125)*IT_0488*IT_0726 + (-0.125)*IT_0450*IT_0741 + (-0.125)*IT_0669
      *IT_0745 + (-0.125)*IT_0428*IT_0750 + (-0.125)*IT_0461*IT_0760 + (-0.125)
      *IT_0473*IT_0777 + (-0.125)*IT_0109*IT_0786 + (-0.125)*IT_0317*IT_0788 + (
      -0.125)*IT_0036*IT_0790 + (-0.125)*IT_0281*IT_0792 + IT_0638*(IT_0199 + (
      -0.125)*IT_0870) + IT_0618*(IT_0232 + IT_0255 + IT_0883 + IT_0884) + (
      -0.125)*IT_1129 + 0.5*IT_0181*IT_1347 + (-0.125)*(IT_0409 + IT_0795)
      *IT_1381 + (-0.125)*(IT_0233 + IT_0268)*IT_1403 + (-0.125)*(IT_0339 +
       IT_0854)*IT_1410 + IT_0365*IT_1411;
    const ccomplex_t IT_1413 = (-8)*conj(IT_0199);
    const ccomplex_t IT_1414 = (-8)*conj(IT_0211);
    const ccomplex_t IT_1415 = (-0.125)*conj(IT_0545)*IT_0685 + (-0.125)
      *IT_0744*IT_0959 + (-0.125)*IT_1170 + (-0.125)*IT_1189 + IT_0618*IT_1290 +
       (-0.125)*IT_1294 + (-0.125)*IT_1295 + IT_1168*IT_1347 + IT_0790*IT_1348 +
       IT_0741*IT_1360 + IT_0792*IT_1361 + IT_0697*IT_1362 + IT_0788*IT_1363 +
       IT_0786*IT_1364 + IT_0750*IT_1365 + IT_0760*IT_1369 + IT_1359*IT_1377 +
       IT_1358*IT_1378 + IT_0707*IT_1379 + IT_0669*IT_1380 + IT_1382*IT_1403 +
       IT_1381*IT_1404 + IT_0777*IT_1405 + IT_0726*IT_1406 + IT_1368*IT_1410 +
       IT_1309*IT_1411 + (-0.125)*IT_0638*(IT_1413 + IT_1414);
    const ccomplex_t IT_1416 = conj(IT_0050)*IT_1143;
    const ccomplex_t IT_1417 = (-0.5)*IT_1416;
    const ccomplex_t IT_1418 = IT_0988*IT_1211;
    const ccomplex_t IT_1419 = 2*IT_1096;
    const ccomplex_t IT_1420 = IT_1146*IT_1227;
    const ccomplex_t IT_1421 = IT_1161*IT_1228;
    const ccomplex_t IT_1422 = IT_0327*IT_1269;
    const ccomplex_t IT_1423 = IT_1097*IT_1233;
    const ccomplex_t IT_1424 = IT_1041*IT_1271;
    const ccomplex_t IT_1425 = conj(IT_0575)*IT_1391;
    const ccomplex_t IT_1426 = (-0.25)*IT_1425;
    const ccomplex_t IT_1427 = conj(IT_0575)*IT_1402;
    const ccomplex_t IT_1428 = (-0.25)*IT_1427;
    const ccomplex_t IT_1429 = conj(IT_0109)*IT_1152;
    const ccomplex_t IT_1430 = (-0.25)*IT_1429;
    const ccomplex_t IT_1431 = (-0.25)*IT_1343;
    const ccomplex_t IT_1432 = conj(IT_0036)*IT_0288;
    const ccomplex_t IT_1433 = (-0.25)*IT_1432;
    const ccomplex_t IT_1434 = conj(IT_0428)*IT_1157;
    const ccomplex_t IT_1435 = (-0.25)*IT_1434;
    const ccomplex_t IT_1436 = (-0.25)*IT_1344;
    const ccomplex_t IT_1437 = conj(IT_0317)*IT_0797;
    const ccomplex_t IT_1438 = (-0.25)*IT_1437;
    const ccomplex_t IT_1439 = IT_0362*conj(IT_0509);
    const ccomplex_t IT_1440 = IT_0928*IT_1222 + IT_1034*IT_1242 + IT_1054
      *IT_1270 + IT_0825*IT_1323 + IT_1417 + IT_1418 + IT_1221*IT_1419 + IT_1420
       + IT_1421 + IT_1422 + IT_1423 + IT_1424 + IT_1426 + IT_1428 + IT_1430 +
       IT_1431 + IT_1433 + IT_1435 + IT_1436 + IT_1438 + IT_1439;
    const ccomplex_t IT_1441 = IT_0295*(conj(IT_0368) + conj(IT_0545)) + (-4)
      *IT_0928*IT_1222 + (-4)*IT_1054*IT_1270 + (-4)*IT_1417 + (-4)*IT_1418 + (
      -4)*IT_1221*IT_1419 + (-4)*IT_1420 + (-4)*IT_1421 + (-4)*IT_1422 + (-4)
      *IT_1423 + (-4)*IT_1424 + (-4)*IT_1426 + (-4)*IT_1428 + (-4)*IT_1430 + (-4
      )*IT_1431 + (-4)*IT_1433 + (-4)*IT_1435 + (-4)*IT_1436 + (-4)*IT_1438 + (
      -4)*IT_1439;
    const ccomplex_t IT_1442 = IT_0446*IT_1030;
    const ccomplex_t IT_1443 = IT_0446*IT_1031;
    const ccomplex_t IT_1444 = conj(IT_0575)*IT_1345;
    const ccomplex_t IT_1445 = conj(IT_0409)*IT_1034;
    const ccomplex_t IT_1446 = conj(IT_0795)*IT_1034;
    const ccomplex_t IT_1447 = conj(IT_0521)*IT_0742;
    const ccomplex_t IT_1448 = conj(IT_0368)*IT_1035;
    const ccomplex_t IT_1449 = conj(IT_0745)*IT_1036;
    const ccomplex_t IT_1450 = conj(IT_0545)*IT_1037;
    const ccomplex_t IT_1451 = conj(IT_0097)*IT_0943;
    const ccomplex_t IT_1452 = conj(IT_0461)*IT_0964;
    const ccomplex_t IT_1453 = conj(IT_0281)*IT_0966;
    const ccomplex_t IT_1454 = conj(IT_0854)*IT_0974;
    const ccomplex_t IT_1455 = conj(IT_0339)*IT_0974;
    const ccomplex_t IT_1456 = conj(IT_0450)*IT_0975;
    const ccomplex_t IT_1457 = conj(IT_0109)*IT_1038;
    const ccomplex_t IT_1458 = conj(IT_0428)*IT_1039;
    const ccomplex_t IT_1459 = conj(IT_0050)*IT_1040;
    const ccomplex_t IT_1460 = conj(IT_0233)*IT_1041;
    const ccomplex_t IT_1461 = conj(IT_0268)*IT_1041;
    const ccomplex_t IT_1462 = conj(IT_0521)*IT_0743;
    const ccomplex_t IT_1463 = conj(IT_0488)*IT_0958;
    const ccomplex_t IT_1464 = conj(IT_0181)*IT_0965;
    const ccomplex_t IT_1465 = conj(IT_0509)*IT_1321;
    const ccomplex_t IT_1466 = IT_0799*IT_0816 + conj(IT_0473)*IT_1023 +
       IT_1124 + IT_1125 + IT_1442 + IT_1443 + IT_1444 + IT_1445 + IT_1446 +
       IT_1447 + IT_1448 + IT_1449 + IT_1450 + IT_1451 + IT_1452 + IT_1453 +
       IT_1454 + IT_1455 + IT_1456 + IT_1457 + IT_1458 + IT_1459 + IT_1460 +
       IT_1461 + IT_1462 + IT_1463 + IT_1464 + IT_1465;
    const ccomplex_t IT_1467 = IT_0952*IT_0959 + conj(IT_0473)*IT_1023 +
       IT_1124 + IT_1125 + IT_1442 + IT_1443 + IT_1444 + IT_1445 + IT_1446 +
       IT_1447 + IT_1448 + IT_1449 + IT_1450 + IT_1451 + IT_1452 + IT_1453 +
       IT_1454 + IT_1455 + IT_1456 + IT_1457 + IT_1458 + IT_1459 + IT_1460 +
       IT_1461 + IT_1462 + IT_1463 + IT_1464 + IT_1465;
    const ccomplex_t IT_1468 = IT_1020*IT_1030;
    const ccomplex_t IT_1469 = IT_1020*IT_1031;
    const ccomplex_t IT_1470 = conj(IT_0575)*IT_1366;
    const ccomplex_t IT_1471 = conj(IT_0109)*IT_1263;
    const ccomplex_t IT_1472 = conj(IT_0409)*IT_1419;
    const ccomplex_t IT_1473 = conj(IT_0050)*IT_1275;
    const ccomplex_t IT_1474 = conj(IT_0575)*IT_1367;
    const ccomplex_t IT_1475 = conj(IT_0428)*IT_1305;
    const ccomplex_t IT_1476 = conj(IT_0281)*IT_1079;
    const ccomplex_t IT_1477 = conj(IT_0461)*IT_0958;
    const ccomplex_t IT_1478 = conj(IT_0368)*IT_1310;
    const ccomplex_t IT_1479 = conj(IT_0745)*IT_1298;
    const ccomplex_t IT_1480 = conj(IT_0795)*IT_1419;
    const ccomplex_t IT_1481 = conj(IT_0545)*IT_1318;
    const ccomplex_t IT_1482 = conj(IT_0181)*IT_0976;
    const ccomplex_t IT_1483 = conj(IT_0509)*IT_1320;
    const ccomplex_t IT_1484 = conj(IT_0488)*IT_0964 + conj(IT_0450)*IT_1023 +
       IT_1025 + IT_1026 + IT_0953*IT_1083 + conj(IT_0097)*IT_1084 + conj
      (IT_0473)*IT_1122 + IT_1468 + IT_1469 + IT_1470 + IT_1471 + IT_1472 +
       IT_1473 + IT_1474 + IT_1475 + IT_1476 + IT_1477 + IT_1478 + IT_1479 +
       IT_1480 + IT_1481 + IT_1482 + IT_1483;
    const ccomplex_t IT_1485 = IT_0799*IT_0804 + conj(IT_0488)*IT_0964 +
       IT_0348*IT_0978 + conj(IT_0450)*IT_1023 + IT_1025 + IT_1026 + IT_0953
      *IT_1083 + conj(IT_0097)*IT_1084 + conj(IT_0473)*IT_1122 + IT_1468 +
       IT_1469 + IT_1470 + IT_1471 + IT_1472 + IT_1473 + IT_1474 + IT_1475 +
       IT_1476 + IT_1477 + IT_1478 + IT_1479 + IT_1480 + IT_1481 + IT_1482 +
       IT_1483;
    const ccomplex_t IT_1486 = (-0.25)*conj(IT_0601);
    const ccomplex_t IT_1487 = conj(IT_0097)*IT_0516;
    const ccomplex_t IT_1488 = (-0.25)*IT_1487;
    const ccomplex_t IT_1489 = conj(IT_0036)*IT_0518;
    const ccomplex_t IT_1490 = (-0.25)*IT_1489;
    const ccomplex_t IT_1491 = IT_0530*conj(IT_0575);
    const ccomplex_t IT_1492 = (-0.25)*IT_1491;
    const ccomplex_t IT_1493 = IT_0541*conj(IT_0575);
    const ccomplex_t IT_1494 = (-0.25)*IT_1493;
    const ccomplex_t IT_1495 = conj(IT_0109)*IT_0543;
    const ccomplex_t IT_1496 = (-0.25)*IT_1495;
    const ccomplex_t IT_1497 = conj(IT_0050)*IT_0546;
    const ccomplex_t IT_1498 = (-0.25)*IT_1497;
    const ccomplex_t IT_1499 = conj(IT_0586)*IT_0810;
    const ccomplex_t IT_1500 = (-0.25)*IT_1499;
    const ccomplex_t IT_1501 = conj(IT_0597)*IT_0810;
    const ccomplex_t IT_1502 = (-0.25)*IT_1501;
    const ccomplex_t IT_1503 = conj(IT_0181)*IT_0816;
    const ccomplex_t IT_1504 = (-0.25)*IT_1503;
    const ccomplex_t IT_1505 = conj(IT_0854)*IT_0855;
    const ccomplex_t IT_1506 = (-0.25)*IT_1505;
    const ccomplex_t IT_1507 = conj(IT_0339)*IT_0855;
    const ccomplex_t IT_1508 = (-0.25)*IT_1507;
    const ccomplex_t IT_1509 = conj(IT_0521)*IT_0861;
    const ccomplex_t IT_1510 = (-0.25)*IT_1509;
    const ccomplex_t IT_1511 = conj(IT_0545)*IT_0568;
    const ccomplex_t IT_1512 = (-0.25)*IT_1511;
    const ccomplex_t IT_1513 = conj(IT_0233)*IT_0825;
    const ccomplex_t IT_1514 = (-0.25)*IT_1513;
    const ccomplex_t IT_1515 = conj(IT_0268)*IT_0825;
    const ccomplex_t IT_1516 = (-0.25)*IT_1515;
    const ccomplex_t IT_1517 = IT_0867*IT_1413;
    const ccomplex_t IT_1518 = (-0.25)*IT_1517;
    const ccomplex_t IT_1519 = IT_0867*IT_1414;
    const ccomplex_t IT_1520 = (-0.25)*IT_1519;
    const ccomplex_t IT_1521 = conj(IT_0281)*IT_0570;
    const ccomplex_t IT_1522 = (-0.25)*IT_1521;
    const ccomplex_t IT_1523 = conj(IT_0317)*IT_0572;
    const ccomplex_t IT_1524 = (-0.25)*IT_1523;
    const ccomplex_t IT_1525 = conj(IT_0509)*IT_0874;
    const ccomplex_t IT_1526 = (-0.25)*IT_1525;
    const ccomplex_t IT_1527 = conj(IT_0488)*IT_0876;
    const ccomplex_t IT_1528 = (-0.25)*IT_1527;
    const ccomplex_t IT_1529 = conj(IT_0450)*IT_0835;
    const ccomplex_t IT_1530 = (-0.25)*IT_1529;
    const ccomplex_t IT_1531 = conj(IT_0473)*IT_0843;
    const ccomplex_t IT_1532 = (-0.25)*IT_1531;
    const ccomplex_t IT_1533 = conj(IT_0598)*IT_0880;
    const ccomplex_t IT_1534 = (-0.25)*IT_1533;
    const ccomplex_t IT_1535 = conj(IT_0745)*IT_1289;
    const ccomplex_t IT_1536 = IT_0004*IT_0948;
    const ccomplex_t IT_1537 = -IT_1536;
    const ccomplex_t IT_1538 = (-4)*IT_0492;
    const ccomplex_t IT_1539 = IT_0369 + IT_0506 + IT_1537 + IT_1538;
    const ccomplex_t IT_1540 = IT_1290*IT_1539;
    const ccomplex_t IT_1541 = IT_1013*IT_1222 + conj(IT_0368)*IT_1245 +
       IT_0880*IT_1486 + IT_1488 + IT_1490 + IT_1492 + IT_1494 + IT_1496 +
       IT_1498 + IT_1500 + IT_1502 + IT_1504 + IT_1506 + IT_1508 + IT_1510 +
       IT_1512 + IT_1514 + IT_1516 + IT_1518 + IT_1520 + IT_1522 + IT_1524 +
       IT_1526 + IT_1528 + IT_1530 + IT_1532 + IT_1534 + IT_1535 + IT_1540;
    const ccomplex_t IT_1542 = conj(IT_0461)*IT_0886;
    const ccomplex_t IT_1543 = IT_0004*IT_0040;
    const ccomplex_t IT_1544 = 12*IT_0001 + (-4)*IT_1543;
    const ccomplex_t IT_1545 = IT_0511 + IT_1544;
    const ccomplex_t IT_1546 = IT_0037*IT_1545;
    const ccomplex_t IT_1547 = (-2)*IT_1546;
    const ccomplex_t IT_1548 = conj(IT_0368)*IT_1245 + IT_1488 + IT_1490 +
       IT_1492 + IT_1494 + IT_1496 + IT_1498 + IT_1500 + IT_1502 + IT_1504 +
       IT_1506 + IT_1508 + IT_1510 + IT_1512 + IT_1514 + IT_1516 + IT_1518 +
       IT_1520 + IT_1522 + IT_1524 + IT_1526 + IT_1528 + IT_1530 + IT_1532 +
       IT_1534 + IT_1535 + IT_1540 + (-0.25)*IT_1542 + IT_1486*IT_1547;
    const ccomplex_t IT_1549 = IT_1307*IT_1322;
    const ccomplex_t IT_1550 = IT_0804*conj(IT_0854);
    const ccomplex_t IT_1551 = (-0.25)*IT_1550;
    const ccomplex_t IT_1552 = conj(IT_0339)*IT_0804;
    const ccomplex_t IT_1553 = (-0.25)*IT_1552;
    const ccomplex_t IT_1554 = conj(IT_0181)*IT_0810;
    const ccomplex_t IT_1555 = (-0.25)*IT_1554;
    const ccomplex_t IT_1556 = conj(IT_0586)*IT_0816;
    const ccomplex_t IT_1557 = (-0.25)*IT_1556;
    const ccomplex_t IT_1558 = conj(IT_0597)*IT_0816;
    const ccomplex_t IT_1559 = (-0.25)*IT_1558;
    const ccomplex_t IT_1560 = conj(IT_0575)*IT_0861;
    const ccomplex_t IT_1561 = (-0.25)*IT_1560;
    const ccomplex_t IT_1562 = conj(IT_0409)*IT_0825;
    const ccomplex_t IT_1563 = (-0.25)*IT_1562;
    const ccomplex_t IT_1564 = conj(IT_0795)*IT_0825;
    const ccomplex_t IT_1565 = (-0.25)*IT_1564;
    const ccomplex_t IT_1566 = IT_0832*IT_0894;
    const ccomplex_t IT_1567 = (-0.25)*IT_1566;
    const ccomplex_t IT_1568 = IT_0832*IT_0895;
    const ccomplex_t IT_1569 = (-0.25)*IT_1568;
    const ccomplex_t IT_1570 = conj(IT_0473)*IT_0876;
    const ccomplex_t IT_1571 = (-0.25)*IT_1570;
    const ccomplex_t IT_1572 = conj(IT_0461)*IT_0835;
    const ccomplex_t IT_1573 = (-0.25)*IT_1572;
    const ccomplex_t IT_1574 = conj(IT_0488)*IT_0843;
    const ccomplex_t IT_1575 = (-0.25)*IT_1574;
    const ccomplex_t IT_1576 = conj(IT_0109)*IT_1245;
    const ccomplex_t IT_1577 = conj(IT_0428)*IT_1289;
    const ccomplex_t IT_1578 = IT_0799*IT_1539;
    const ccomplex_t IT_1579 = IT_1013*IT_1211 + IT_1325 + IT_1327 + IT_1329 +
       IT_1331 + IT_1333 + IT_1335 + IT_1337 + IT_1339 + IT_1341 + IT_1549 +
       IT_1551 + IT_1553 + IT_1555 + IT_1557 + IT_1559 + IT_1561 + IT_1563 +
       IT_1565 + IT_1567 + IT_1569 + IT_1571 + IT_1573 + IT_1575 + IT_1576 +
       IT_1577 + IT_1578;
    const ccomplex_t IT_1580 = conj(IT_0232)*IT_0880 + (-4)*IT_0886*IT_1211 + 
      (-4)*IT_1325 + (-4)*IT_1327 + (-4)*IT_1329 + (-4)*IT_1331 + (-4)*IT_1333 +
       (-4)*IT_1335 + (-4)*IT_1337 + (-4)*IT_1339 + (-4)*IT_1341 + IT_0885
      *IT_1547 + (-4)*IT_1549 + (-4)*IT_1551 + (-4)*IT_1553 + (-4)*IT_1555 + (-4
      )*IT_1557 + (-4)*IT_1559 + (-4)*IT_1561 + (-4)*IT_1563 + (-4)*IT_1565 + (
      -4)*IT_1567 + (-4)*IT_1569 + (-4)*IT_1571 + (-4)*IT_1573 + (-4)*IT_1575 + 
      (-4)*IT_1576 + (-4)*IT_1577 + (-4)*IT_1578;
    const ccomplex_t IT_1581 = (IT_0017*conj(IT_0036) + IT_0048*conj(IT_0050) 
      + IT_0063*conj(IT_0097) + IT_0107*conj(IT_0109) + IT_0116*conj(IT_0181))
      *IT_0213 + (conj(IT_0233) + conj(IT_0268))*IT_0367 + conj(IT_0368)*
      (IT_0365*IT_0385 + IT_0435) + (IT_0063*conj(IT_0281) + IT_0017*conj
      (IT_0317) + IT_0446*conj(IT_0450) + IT_0116*conj(IT_0461) + IT_0466*conj
      (IT_0473) + IT_0481*conj(IT_0488))*IT_0491 + IT_0507*(conj(IT_0211)
      *IT_0364 + conj(IT_0199)*IT_0365 + IT_0199*IT_0508) + conj(IT_0509)*
      (IT_0491*IT_0515 + IT_0517 + IT_0519 + IT_0531 + IT_0542 + IT_0544 +
       IT_0547 + IT_0569 + IT_0571 + IT_0573) + (-4)*conj(IT_0575)*IT_0794 +
       IT_0796*(IT_0789 + IT_0791 + IT_0798) + IT_0799*(IT_0517 + IT_0519 +
       IT_0531 + IT_0542 + IT_0544 + IT_0547 + IT_0569 + IT_0571 + IT_0573 +
       IT_0805 + IT_0811 + IT_0817 + IT_0826 + IT_0833 + IT_0834 + IT_0836 +
       IT_0844) + conj(IT_0232)*IT_0882 + IT_0885*IT_0888 + IT_0461*IT_1015 +
       conj(IT_0450)*(IT_1016 + IT_1017) + IT_0953*(IT_1016 + IT_1017 + IT_1021 
      + IT_1022) + IT_0450*IT_1027 + (conj(IT_0181) + conj(IT_0488))*(IT_1028 +
       IT_1029) + IT_0181*IT_1044 + IT_0488*IT_1111 + conj(IT_0473)*(IT_1112 +
       IT_1113) + IT_0959*(IT_1112 + IT_1113 + IT_1114 + IT_1115 + IT_1116 +
       IT_1117 + IT_1118 + IT_1119 + IT_1120 + IT_1121) + IT_0473*IT_1126 + conj
      (IT_0545)*(IT_0435 + 4*IT_0567*(IT_0199 + IT_1128)) + conj(IT_0745)*
      (IT_1129 + 4*(IT_0199 + IT_1128)*IT_1139) + 2*IT_0268*IT_1201 + 2*IT_0233
      *IT_1203 + (-4)*IT_0097*IT_1234 + (-4)*IT_0281*IT_1243 + (-4)*IT_0109
      *IT_1272 + (-4)*IT_0050*IT_1279 + (-4)*IT_0317*IT_1281 + (-4)*IT_0036
      *IT_1282 + 2*IT_0745*IT_1302 + 2*IT_0428*IT_1308 + (-4)*IT_0368*IT_1317 + 
      (-4)*IT_0545*IT_1319 + (-4)*IT_0509*IT_1342 + (-8)*IT_0575*IT_1407 + (-8)
      *conj(IT_0521)*IT_1412 + (-8)*IT_0521*IT_1415 + (-4)*IT_0409*IT_1440 +
       IT_0795*IT_1441 + IT_0586*IT_1466 + IT_0597*IT_1467 + IT_0339*IT_1484 +
       IT_0854*IT_1485 + (-4)*IT_0232*IT_1541 + (-4)*(IT_0255 + IT_0883 +
       IT_0884)*IT_1548 + (-4)*IT_0598*IT_1579 + IT_0601*IT_1580;
    return create_ccomplex_return(IT_1581);
}

