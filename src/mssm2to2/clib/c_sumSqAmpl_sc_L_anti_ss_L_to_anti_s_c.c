#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_anti_s_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_anti_s_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
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
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B2)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W2)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + 3*IT_0009);
    const ccomplex_t IT_0011 = 0.166666666666667*IT_0010;
    const ccomplex_t IT_0012 = N_B2*e_em;
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = N_W2*e_em;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0014 + (-3)*IT_0017);
    const ccomplex_t IT_0019 = 0.166666666666667*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0019;
    const ccomplex_t IT_0021 = pow(m_s, 2);
    const ccomplex_t IT_0022 = pow(m_ss_L, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = conj(N_B3)*e_em;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = conj(N_W3)*e_em;
    const ccomplex_t IT_0030 = IT_0006*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + 3*IT_0031);
    const ccomplex_t IT_0033 = 0.166666666666667*IT_0032;
    const ccomplex_t IT_0034 = N_B3*e_em;
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_W3*e_em;
    const ccomplex_t IT_0038 = IT_0006*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + (-3)*IT_0039);
    const ccomplex_t IT_0041 = 0.166666666666667*IT_0040;
    const ccomplex_t IT_0042 = IT_0033*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = conj(N_B4)*e_em;
    const ccomplex_t IT_0047 = IT_0001*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W4)*e_em;
    const ccomplex_t IT_0050 = IT_0006*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + 3*IT_0051);
    const ccomplex_t IT_0053 = 0.166666666666667*IT_0052;
    const ccomplex_t IT_0054 = N_B4*e_em;
    const ccomplex_t IT_0055 = IT_0001*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = N_W4*e_em;
    const ccomplex_t IT_0058 = IT_0006*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + (-3)*IT_0059);
    const ccomplex_t IT_0061 = 0.166666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = pow(m_sG, 2);
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -IT_0066 +
       reg_prop, -1);
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0069 = cpow(IT_0068, 2);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0067*IT_0070;
    const ccomplex_t IT_0072 = conj(N_B1)*e_em;
    const ccomplex_t IT_0073 = IT_0001*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W1)*e_em;
    const ccomplex_t IT_0076 = IT_0006*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + 3*IT_0077);
    const ccomplex_t IT_0079 = 0.166666666666667*IT_0078;
    const ccomplex_t IT_0080 = N_B1*e_em;
    const ccomplex_t IT_0081 = IT_0001*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = N_W1*e_em;
    const ccomplex_t IT_0084 = IT_0006*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + (-3)*IT_0085);
    const ccomplex_t IT_0087 = 0.166666666666667*IT_0086;
    const ccomplex_t IT_0088 = IT_0079*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0021 + IT_0022 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = -IT_0025 + -IT_0045 + -IT_0065 +
       0.166666666666667*IT_0071 + -IT_0091;
    const ccomplex_t IT_0093 = pow(m_W, -1);
    const ccomplex_t IT_0094 = cos(beta);
    const ccomplex_t IT_0095 = cpow(IT_0094, -1);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0006*IT_0093*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = sin(beta);
    const ccomplex_t IT_0099 = cpow(IT_0098, -1);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0006*IT_0093*IT_0099;
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0097*IT_0101;
    const ccomplex_t IT_0103 = IT_0089*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0006*IT_0093*IT_0095;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0006*IT_0093*IT_0099;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0106*IT_0108;
    const ccomplex_t IT_0110 = IT_0023*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0006*IT_0093*IT_0095;
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0006*IT_0093*IT_0099;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = IT_0113*IT_0115;
    const ccomplex_t IT_0117 = IT_0043*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0006*IT_0093*IT_0095;
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0006*IT_0093*IT_0099;
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = IT_0063*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = -IT_0104 + -IT_0111 + -IT_0118 + -IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0006;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0006;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = pow(m_sc_L, 2);
    const ccomplex_t IT_0132 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0022 + -IT_0131 + -reg_prop, -1);
    const ccomplex_t IT_0133 = IT_0130*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0128*IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0071 + -IT_0135;
    const ccomplex_t IT_0137 = s_34*IT_0131;
    const ccomplex_t IT_0138 = (-18)*IT_0137;
    const ccomplex_t IT_0139 = s_13*s_14;
    const ccomplex_t IT_0140 = 36*IT_0139;
    const ccomplex_t IT_0141 = IT_0138 + IT_0140;
    const ccomplex_t IT_0142 = m_c*m_s*IT_0131;
    const ccomplex_t IT_0143 = (-6)*IT_0136;
    const ccomplex_t IT_0144 = (-18)*IT_0092 + IT_0143;
    const ccomplex_t IT_0145 = (-6)*conj(IT_0136);
    const ccomplex_t IT_0146 = (-18)*conj(IT_0092);
    const ccomplex_t IT_0147 = IT_0145 + IT_0146;
    const ccomplex_t IT_0148 = pow(m_c, 2);
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0006*(IT_0021*IT_0093*IT_0095*IT_0098 + (-2)*IT_0094*(m_W*IT_0098
       + (-0.5)*IT_0093*IT_0099*IT_0148));
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0006*IT_0093*IT_0094*IT_0099;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = IT_0150*IT_0152;
    const ccomplex_t IT_0154 = cpow((-2)*s_12 + -IT_0022 + -IT_0131 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0155 = IT_0153*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = m_c*IT_0071;
    const ccomplex_t IT_0158 = (-0.5)*IT_0127;
    const ccomplex_t IT_0159 = m_c*IT_0158;
    const ccomplex_t IT_0160 = m_c*IT_0128;
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = IT_0159 + IT_0161;
    const ccomplex_t IT_0163 = IT_0134*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0006*IT_0093;
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = cpow(IT_0094, 2);
    const ccomplex_t IT_0167 = cpow(IT_0098, 2);
    const ccomplex_t IT_0168 = (0 + _Complex_I*(-1.4142135623731))*conj(V_cs)
      *e_em*IT_0006*(IT_0093*(IT_0021 + -IT_0148) + -m_W*(IT_0166 + -IT_0167));
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = IT_0165*IT_0169;
    const ccomplex_t IT_0171 = pow(m_W, 2);
    const ccomplex_t IT_0172 = cpow((-2)*s_12 + -IT_0022 + -IT_0131 + IT_0171 
      + -reg_prop, -1);
    const ccomplex_t IT_0173 = IT_0170*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = -IT_0156 + 0.5*IT_0157 + -IT_0163 + -IT_0174;
    const ccomplex_t IT_0176 = IT_0087*IT_0101;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0178 = IT_0089*IT_0176*IT_0177;
    const ccomplex_t IT_0179 = IT_0019*IT_0108;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0181 = IT_0023*IT_0179*IT_0180;
    const ccomplex_t IT_0182 = IT_0041*IT_0115;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0184 = IT_0043*IT_0182*IT_0183;
    const ccomplex_t IT_0185 = IT_0061*IT_0122;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0187 = IT_0063*IT_0185*IT_0186;
    const ccomplex_t IT_0188 = m_c*IT_0025;
    const ccomplex_t IT_0189 = m_c*IT_0045;
    const ccomplex_t IT_0190 = m_c*IT_0065;
    const ccomplex_t IT_0191 = m_c*IT_0091;
    const ccomplex_t IT_0192 = (-0.166666666666667)*IT_0157 + -IT_0178 + 
      -IT_0181 + -IT_0184 + -IT_0187 + IT_0188 + IT_0189 + IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = IT_0079*IT_0097;
    const ccomplex_t IT_0194 = IT_0089*IT_0177*IT_0193;
    const ccomplex_t IT_0195 = IT_0011*IT_0106;
    const ccomplex_t IT_0196 = IT_0023*IT_0180*IT_0195;
    const ccomplex_t IT_0197 = IT_0033*IT_0113;
    const ccomplex_t IT_0198 = IT_0043*IT_0183*IT_0197;
    const ccomplex_t IT_0199 = IT_0053*IT_0120;
    const ccomplex_t IT_0200 = IT_0063*IT_0186*IT_0199;
    const ccomplex_t IT_0201 = m_c*IT_0104;
    const ccomplex_t IT_0202 = m_c*IT_0111;
    const ccomplex_t IT_0203 = m_c*IT_0118;
    const ccomplex_t IT_0204 = m_c*IT_0125;
    const ccomplex_t IT_0205 = -IT_0194 + -IT_0196 + -IT_0198 + -IT_0200 +
       IT_0201 + IT_0202 + IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0006*IT_0093;
    const ccomplex_t IT_0207 = (-0.5)*IT_0206;
    const ccomplex_t IT_0208 = IT_0169*IT_0207;
    const ccomplex_t IT_0209 = IT_0172*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*IT_0209;
    const ccomplex_t IT_0211 = m_s*IT_0158;
    const ccomplex_t IT_0212 = m_s*IT_0128;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = IT_0211 + IT_0213;
    const ccomplex_t IT_0215 = IT_0134*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0006*IT_0093*IT_0095*IT_0098;
    const ccomplex_t IT_0217 = 0.5*IT_0216;
    const ccomplex_t IT_0218 = IT_0150*IT_0217;
    const ccomplex_t IT_0219 = IT_0154*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = -IT_0210 + IT_0215 + -IT_0220;
    const ccomplex_t IT_0222 = 6*s_34;
    const ccomplex_t IT_0223 = 18*s_34;
    const ccomplex_t IT_0224 = m_c*m_s;
    const ccomplex_t IT_0225 = (-18)*conj(IT_0192);
    const ccomplex_t IT_0226 = (-6)*conj(IT_0175);
    const ccomplex_t IT_0227 = IT_0225 + IT_0226;
    const ccomplex_t IT_0228 = (-18)*IT_0192;
    const ccomplex_t IT_0229 = 12*IT_0139;
    const ccomplex_t IT_0230 = (-6)*IT_0137;
    const ccomplex_t IT_0231 = m_c*s_13;
    const ccomplex_t IT_0232 = 3*IT_0136;
    const ccomplex_t IT_0233 = 3*conj(IT_0136);
    const ccomplex_t IT_0234 = (conj(IT_0092) + 0.333333333333333*conj(IT_0136
      ))*IT_0192 + (IT_0092 + 0.333333333333333*IT_0136)*conj(IT_0192) + conj
      (IT_0126)*(IT_0205 + 0.333333333333333*IT_0221) + IT_0126*(conj(IT_0205) +
       0.333333333333333*conj(IT_0221)) + 0.333333333333333*conj(IT_0175)*
      (IT_0092 + IT_0232) + 0.333333333333333*IT_0175*(conj(IT_0092) + IT_0233);
    const ccomplex_t IT_0235 = m_s*s_14;
    const ccomplex_t IT_0236 = (IT_0092*conj(IT_0092) + IT_0126*conj(IT_0126) 
      + IT_0136*conj(IT_0136))*IT_0141 + IT_0142*(conj(IT_0126)*IT_0144 +
       IT_0126*IT_0147) + (conj(IT_0175)*IT_0192 + IT_0175*conj(IT_0192) + conj
      (IT_0205)*IT_0221 + IT_0205*conj(IT_0221))*IT_0222 + (IT_0175*conj(IT_0175
      ) + IT_0192*conj(IT_0192) + IT_0205*conj(IT_0205) + IT_0221*conj(IT_0221))
      *IT_0223 + (-18)*IT_0224*((conj(IT_0175) + 0.333333333333333*conj(IT_0192)
      )*IT_0221 + 0.333333333333333*IT_0192*conj(IT_0221) + 0.333333333333333
      *IT_0175*(conj(IT_0205) + 3*conj(IT_0221)) + (-0.0555555555555556)*IT_0205
      *IT_0227 + (-0.0555555555555556)*conj(IT_0205)*IT_0228) + (conj(IT_0092)
      *IT_0136 + IT_0092*conj(IT_0136))*(IT_0229 + IT_0230) + 18*IT_0231*IT_0234
       + (IT_0147*IT_0205 + IT_0144*conj(IT_0205) + IT_0126*IT_0227 + (-6)*conj
      (IT_0126)*(IT_0175 + (-0.166666666666667)*IT_0228) + (-6)*conj(IT_0221)*
      (IT_0092 + IT_0232) + (-6)*IT_0221*(conj(IT_0092) + IT_0233))*IT_0235;
    return create_ccomplex_return(IT_0236);
}

