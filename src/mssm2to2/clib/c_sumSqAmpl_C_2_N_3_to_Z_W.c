#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_Z_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_Z_W(
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0005 + 0.5*IT_0008);
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = pow(m_C_2, 2);
    const ccomplex_t IT_0012 = pow(m_N_3, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -2);
    const ccomplex_t IT_0019 = pow(m_Z, -2);
    const ccomplex_t IT_0020 = s_13*s_24*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0021 = s_13*s_23*IT_0019;
    const ccomplex_t IT_0022 = s_14*s_24*IT_0018;
    const ccomplex_t IT_0023 = -IT_0021 + -IT_0022;
    const ccomplex_t IT_0024 = s_12 + IT_0020 + IT_0023;
    const ccomplex_t IT_0025 = m_N_3*IT_0024;
    const ccomplex_t IT_0026 = m_N_3*IT_0003;
    const ccomplex_t IT_0027 = IT_0015*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = 2*conj(IT_0028);
    const ccomplex_t IT_0030 = pow(m_W, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + IT_0030 
      + -reg_prop, -1);
    const ccomplex_t IT_0032 = cpow(IT_0000, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_W*e_em*IT_0001*IT_0032;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = cos(beta);
    const ccomplex_t IT_0036 = IT_0002*IT_0035;
    const ccomplex_t IT_0037 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (-2)*IT_0038;
    const ccomplex_t IT_0040 = IT_0032*IT_0035;
    const ccomplex_t IT_0041 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0045 = IT_0036*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = IT_0043 + IT_0046;
    const ccomplex_t IT_0048 = IT_0039 + IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0034*IT_0050;
    const ccomplex_t IT_0052 = IT_0031*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = m_C_2*IT_0003;
    const ccomplex_t IT_0055 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0056 = IT_0002*IT_0055;
    const ccomplex_t IT_0057 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0058 = IT_0002*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + (-0.5)*IT_0059);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0013*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0054*IT_0063;
    const ccomplex_t IT_0065 = IT_0053 + -IT_0064;
    const ccomplex_t IT_0066 = (-2)*conj(IT_0065);
    const ccomplex_t IT_0067 = IT_0029 + IT_0066;
    const ccomplex_t IT_0068 = pow(s_23, 2);
    const ccomplex_t IT_0069 = IT_0019*IT_0068;
    const ccomplex_t IT_0070 = pow(s_24, 2);
    const ccomplex_t IT_0071 = IT_0018*IT_0070;
    const ccomplex_t IT_0072 = IT_0069 + IT_0071;
    const ccomplex_t IT_0073 = s_23*s_24*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0074 = -IT_0012 + -IT_0073;
    const ccomplex_t IT_0075 = IT_0072 + IT_0074;
    const ccomplex_t IT_0076 = m_C_2*IT_0075;
    const ccomplex_t IT_0077 = IT_0001*IT_0032;
    const ccomplex_t IT_0078 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = IT_0000*IT_0002;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0083 = IT_0077*IT_0082;
    const ccomplex_t IT_0084 = IT_0080*IT_0082;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0079 + IT_0081 + 
      -IT_0083 + -IT_0084);
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0088 = IT_0002*IT_0087;
    const ccomplex_t IT_0089 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0090 = IT_0002*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + (-0.5)*IT_0091);
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = pow(m_Z, 2);
    const ccomplex_t IT_0095 = cpow((-2)*s_23 + IT_0012 + IT_0094 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0086*IT_0097;
    const ccomplex_t IT_0099 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0100 = IT_0077*IT_0099;
    const ccomplex_t IT_0101 = IT_0080*IT_0099;
    const ccomplex_t IT_0102 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0103 = IT_0077*IT_0102;
    const ccomplex_t IT_0104 = IT_0080*IT_0102;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0100 + IT_0101 + 
      -IT_0103 + -IT_0104);
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0094 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0061*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0106*IT_0109;
    const ccomplex_t IT_0111 = (-2)*IT_0098 + -IT_0110;
    const ccomplex_t IT_0112 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0113 = IT_0077*IT_0112;
    const ccomplex_t IT_0114 = IT_0080*IT_0112;
    const ccomplex_t IT_0115 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0116 = IT_0077*IT_0115;
    const ccomplex_t IT_0117 = IT_0080*IT_0115;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0113 + IT_0114 + 
      -IT_0116 + -IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0121 = IT_0002*IT_0120;
    const ccomplex_t IT_0122 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0123 = IT_0002*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + (-0.5)*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_23 + IT_0012 + IT_0094 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0126*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0119*IT_0129;
    const ccomplex_t IT_0131 = (-2)*IT_0003;
    const ccomplex_t IT_0132 = IT_0063*IT_0131;
    const ccomplex_t IT_0133 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0134 = IT_0077*IT_0133;
    const ccomplex_t IT_0135 = IT_0080*IT_0133;
    const ccomplex_t IT_0136 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0137 = IT_0077*IT_0136;
    const ccomplex_t IT_0138 = IT_0080*IT_0136;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0134 + IT_0135 + 
      -IT_0137 + -IT_0138);
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0142 = IT_0002*IT_0141;
    const ccomplex_t IT_0143 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0144 = IT_0002*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0142 + (-0.5)*IT_0145);
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = cpow((-2)*s_23 + IT_0012 + IT_0094 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0149 = IT_0147*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = IT_0140*IT_0150;
    const ccomplex_t IT_0152 = (-0.5)*IT_0105;
    const ccomplex_t IT_0153 = IT_0109*IT_0152;
    const ccomplex_t IT_0154 = 2*IT_0130 + -IT_0132 + 2*IT_0151 + IT_0153;
    const ccomplex_t IT_0155 = conj(IT_0111) + conj(IT_0154);
    const ccomplex_t IT_0156 = 2*IT_0028;
    const ccomplex_t IT_0157 = (-2)*IT_0065;
    const ccomplex_t IT_0158 = IT_0156 + IT_0157;
    const ccomplex_t IT_0159 = IT_0003*IT_0063;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0010*IT_0107;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = IT_0152*IT_0162;
    const ccomplex_t IT_0164 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0165 = IT_0077*IT_0164;
    const ccomplex_t IT_0166 = IT_0080*IT_0164;
    const ccomplex_t IT_0167 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0168 = IT_0077*IT_0167;
    const ccomplex_t IT_0169 = IT_0080*IT_0167;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0165 + IT_0166 + 
      -IT_0168 + -IT_0169);
    const ccomplex_t IT_0171 = (-0.5)*IT_0170;
    const ccomplex_t IT_0172 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0173 = IT_0002*IT_0172;
    const ccomplex_t IT_0174 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0175 = IT_0002*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0173 + 0.5*IT_0176);
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = IT_0095*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = IT_0171*IT_0180;
    const ccomplex_t IT_0182 = IT_0015*IT_0131;
    const ccomplex_t IT_0183 = -IT_0163 + (-2)*IT_0181 + -IT_0182;
    const ccomplex_t IT_0184 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0185 = IT_0077*IT_0184;
    const ccomplex_t IT_0186 = IT_0080*IT_0184;
    const ccomplex_t IT_0187 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0188 = IT_0077*IT_0187;
    const ccomplex_t IT_0189 = IT_0080*IT_0187;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*(IT_0185 + IT_0186 + 
      -IT_0188 + -IT_0189);
    const ccomplex_t IT_0191 = (-0.5)*IT_0190;
    const ccomplex_t IT_0192 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0193 = IT_0002*IT_0192;
    const ccomplex_t IT_0194 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0195 = IT_0002*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + 0.5*IT_0196);
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0148*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = IT_0191*IT_0200;
    const ccomplex_t IT_0202 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0203 = IT_0077*IT_0202;
    const ccomplex_t IT_0204 = IT_0080*IT_0202;
    const ccomplex_t IT_0205 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0206 = IT_0077*IT_0205;
    const ccomplex_t IT_0207 = IT_0080*IT_0205;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0203 + IT_0204 + 
      -IT_0206 + -IT_0207);
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0211 = IT_0002*IT_0210;
    const ccomplex_t IT_0212 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0213 = IT_0002*IT_0212;
    const ccomplex_t IT_0214 = 1.4142135623731*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*(IT_0211 + 0.5*IT_0214);
    const ccomplex_t IT_0216 = -IT_0215;
    const ccomplex_t IT_0217 = IT_0127*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = IT_0209*IT_0218;
    const ccomplex_t IT_0220 = IT_0106*IT_0162;
    const ccomplex_t IT_0221 = 2*IT_0201 + 2*IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = IT_0015*IT_0054;
    const ccomplex_t IT_0223 = 2*IT_0222;
    const ccomplex_t IT_0224 = IT_0026*IT_0063;
    const ccomplex_t IT_0225 = sin(beta);
    const ccomplex_t IT_0226 = IT_0002*IT_0225;
    const ccomplex_t IT_0227 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0228 = IT_0226*IT_0227;
    const ccomplex_t IT_0229 = 2*IT_0228;
    const ccomplex_t IT_0230 = IT_0032*IT_0225;
    const ccomplex_t IT_0231 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0232 = IT_0230*IT_0231;
    const ccomplex_t IT_0233 = 1.4142135623731*IT_0232;
    const ccomplex_t IT_0234 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0235 = IT_0226*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = IT_0233 + IT_0236;
    const ccomplex_t IT_0238 = IT_0229 + IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = 0.5*IT_0239;
    const ccomplex_t IT_0241 = IT_0034*IT_0240;
    const ccomplex_t IT_0242 = IT_0031*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = IT_0224 + -IT_0243;
    const ccomplex_t IT_0245 = (-2)*IT_0244;
    const ccomplex_t IT_0246 = IT_0223 + IT_0245;
    const ccomplex_t IT_0247 = 2*conj(IT_0222);
    const ccomplex_t IT_0248 = (-2)*conj(IT_0244);
    const ccomplex_t IT_0249 = IT_0247 + IT_0248;
    const ccomplex_t IT_0250 = s_23*IT_0094;
    const ccomplex_t IT_0251 = IT_0019*IT_0250;
    const ccomplex_t IT_0252 = s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0253 = s_24*IT_0094;
    const ccomplex_t IT_0254 = IT_0252*IT_0253;
    const ccomplex_t IT_0255 = s_23*IT_0030;
    const ccomplex_t IT_0256 = s_24*s_34;
    const ccomplex_t IT_0257 = (-2)*IT_0256;
    const ccomplex_t IT_0258 = IT_0255 + IT_0257;
    const ccomplex_t IT_0259 = IT_0018*IT_0258;
    const ccomplex_t IT_0260 = 0.5*IT_0251 + (-0.5)*IT_0254 + (-0.5)*IT_0259;
    const ccomplex_t IT_0261 = s_23 + IT_0260;
    const ccomplex_t IT_0262 = m_C_2*IT_0261;
    const ccomplex_t IT_0263 = 2*IT_0262;
    const ccomplex_t IT_0264 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0265 = IT_0002*IT_0264;
    const ccomplex_t IT_0266 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0267 = IT_0002*IT_0266;
    const ccomplex_t IT_0268 = 1.4142135623731*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*(IT_0265 + (-0.5)*IT_0268);
    const ccomplex_t IT_0270 = -IT_0269;
    const ccomplex_t IT_0271 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0272 = IT_0077*IT_0271;
    const ccomplex_t IT_0273 = IT_0080*IT_0271;
    const ccomplex_t IT_0274 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0275 = IT_0080*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*(IT_0272 + -IT_0273 + (-2)
      *IT_0275);
    const ccomplex_t IT_0277 = (-0.5)*IT_0276;
    const ccomplex_t IT_0278 = IT_0270*IT_0277;
    const ccomplex_t IT_0279 = cpow((-2)*s_13 + IT_0011 + IT_0094 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0280 = IT_0278*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = cpow((-2)*s_13 + IT_0094 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0283 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0284 = IT_0077*IT_0283;
    const ccomplex_t IT_0285 = IT_0080*IT_0283;
    const ccomplex_t IT_0286 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0287 = IT_0080*IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*(IT_0284 + -IT_0285 + (-2)
      *IT_0287);
    const ccomplex_t IT_0289 = (-0.5)*IT_0288;
    const ccomplex_t IT_0290 = IT_0061*IT_0289;
    const ccomplex_t IT_0291 = IT_0282*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = IT_0098 + 0.5*IT_0110;
    const ccomplex_t IT_0294 = -IT_0130 + (-0.5)*IT_0153;
    const ccomplex_t IT_0295 = -IT_0151;
    const ccomplex_t IT_0296 = conj(IT_0293) + conj(IT_0294) + conj(IT_0295);
    const ccomplex_t IT_0297 = conj(IT_0281) + conj(IT_0292) + IT_0296;
    const ccomplex_t IT_0298 = IT_0293 + IT_0294 + IT_0295;
    const ccomplex_t IT_0299 = IT_0281 + IT_0292 + IT_0298;
    const ccomplex_t IT_0300 = 0.5*IT_0163 + IT_0181;
    const ccomplex_t IT_0301 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0302 = IT_0002*IT_0301;
    const ccomplex_t IT_0303 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0304 = IT_0002*IT_0303;
    const ccomplex_t IT_0305 = 1.4142135623731*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*(IT_0302 + 0.5*IT_0305);
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0309 = IT_0077*IT_0308;
    const ccomplex_t IT_0310 = IT_0080*IT_0308;
    const ccomplex_t IT_0311 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0312 = IT_0080*IT_0311;
    const ccomplex_t IT_0313 = (0 + _Complex_I*1)*(IT_0309 + -IT_0310 + (-2)
      *IT_0312);
    const ccomplex_t IT_0314 = 0.5*IT_0313;
    const ccomplex_t IT_0315 = IT_0307*IT_0314;
    const ccomplex_t IT_0316 = IT_0279*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*IT_0316;
    const ccomplex_t IT_0318 = -IT_0317;
    const ccomplex_t IT_0319 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0320 = IT_0077*IT_0319;
    const ccomplex_t IT_0321 = IT_0080*IT_0319;
    const ccomplex_t IT_0322 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0323 = IT_0080*IT_0322;
    const ccomplex_t IT_0324 = (0 + _Complex_I*1)*(IT_0320 + -IT_0321 + (-2)
      *IT_0323);
    const ccomplex_t IT_0325 = 0.5*IT_0324;
    const ccomplex_t IT_0326 = IT_0010*IT_0325;
    const ccomplex_t IT_0327 = IT_0282*IT_0326;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*IT_0327;
    const ccomplex_t IT_0329 = -IT_0328;
    const ccomplex_t IT_0330 = -IT_0201 + -IT_0219 + (-0.5)*IT_0220;
    const ccomplex_t IT_0331 = s_13*IT_0094;
    const ccomplex_t IT_0332 = IT_0019*IT_0331;
    const ccomplex_t IT_0333 = s_14*IT_0094;
    const ccomplex_t IT_0334 = IT_0252*IT_0333;
    const ccomplex_t IT_0335 = s_13*IT_0030;
    const ccomplex_t IT_0336 = s_14*s_34;
    const ccomplex_t IT_0337 = (-2)*IT_0336;
    const ccomplex_t IT_0338 = IT_0335 + IT_0337;
    const ccomplex_t IT_0339 = IT_0018*IT_0338;
    const ccomplex_t IT_0340 = 0.5*IT_0332 + (-0.5)*IT_0334 + (-0.5)*IT_0339;
    const ccomplex_t IT_0341 = s_13 + IT_0340;
    const ccomplex_t IT_0342 = m_N_3*IT_0341;
    const ccomplex_t IT_0343 = 2*IT_0342;
    const ccomplex_t IT_0344 = 2*IT_0244;
    const ccomplex_t IT_0345 = (-2)*IT_0222;
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = 2*conj(IT_0244);
    const ccomplex_t IT_0348 = (-2)*conj(IT_0222);
    const ccomplex_t IT_0349 = IT_0347 + IT_0348;
    const ccomplex_t IT_0350 = 2*IT_0065;
    const ccomplex_t IT_0351 = (-2)*IT_0028;
    const ccomplex_t IT_0352 = IT_0350 + IT_0351;
    const ccomplex_t IT_0353 = 2*conj(IT_0065);
    const ccomplex_t IT_0354 = (-2)*conj(IT_0028);
    const ccomplex_t IT_0355 = IT_0353 + IT_0354;
    const ccomplex_t IT_0356 = IT_0018*IT_0019;
    const ccomplex_t IT_0357 = s_13*IT_0030*IT_0094;
    const ccomplex_t IT_0358 = s_14*s_34*IT_0094;
    const ccomplex_t IT_0359 = (-2)*IT_0358;
    const ccomplex_t IT_0360 = IT_0357 + IT_0359;
    const ccomplex_t IT_0361 = IT_0356*IT_0360;
    const ccomplex_t IT_0362 = (-0.666666666666667)*IT_0339 +
       0.166666666666667*IT_0361;
    const ccomplex_t IT_0363 = s_13 + IT_0362;
    const ccomplex_t IT_0364 = m_N_3*IT_0363;
    const ccomplex_t IT_0365 = 6*IT_0364;
    const ccomplex_t IT_0366 = m_N_3*IT_0152;
    const ccomplex_t IT_0367 = m_N_3*IT_0106;
    const ccomplex_t IT_0368 = IT_0109*(IT_0366 + -IT_0367);
    const ccomplex_t IT_0369 = 0.5*IT_0368;
    const ccomplex_t IT_0370 = 2*IT_0369;
    const ccomplex_t IT_0371 = m_N_2*IT_0209;
    const ccomplex_t IT_0372 = m_N_3*IT_0119;
    const ccomplex_t IT_0373 = IT_0129*(IT_0371 + -IT_0372);
    const ccomplex_t IT_0374 = -IT_0373;
    const ccomplex_t IT_0375 = 2*IT_0374;
    const ccomplex_t IT_0376 = m_N_1*IT_0191;
    const ccomplex_t IT_0377 = m_N_3*IT_0140;
    const ccomplex_t IT_0378 = IT_0150*(IT_0376 + -IT_0377);
    const ccomplex_t IT_0379 = -IT_0378;
    const ccomplex_t IT_0380 = 2*IT_0379;
    const ccomplex_t IT_0381 = IT_0370 + IT_0375 + IT_0380;
    const ccomplex_t IT_0382 = -IT_0367;
    const ccomplex_t IT_0383 = IT_0366 + IT_0382;
    const ccomplex_t IT_0384 = IT_0109*IT_0383;
    const ccomplex_t IT_0385 = m_N_4*IT_0171;
    const ccomplex_t IT_0386 = m_N_3*IT_0086;
    const ccomplex_t IT_0387 = -IT_0386;
    const ccomplex_t IT_0388 = IT_0385 + IT_0387;
    const ccomplex_t IT_0389 = IT_0097*IT_0388;
    const ccomplex_t IT_0390 = (-0.5)*IT_0384 + -IT_0389;
    const ccomplex_t IT_0391 = (-2)*IT_0390;
    const ccomplex_t IT_0392 = IT_0381 + IT_0391;
    const ccomplex_t IT_0393 = 2*conj(IT_0369);
    const ccomplex_t IT_0394 = 2*conj(IT_0374);
    const ccomplex_t IT_0395 = 2*conj(IT_0379);
    const ccomplex_t IT_0396 = IT_0393 + IT_0394 + IT_0395;
    const ccomplex_t IT_0397 = (-2)*conj(IT_0390);
    const ccomplex_t IT_0398 = IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = conj(IT_0300) + conj(IT_0330);
    const ccomplex_t IT_0400 = IT_0162*IT_0383;
    const ccomplex_t IT_0401 = m_N_2*IT_0119;
    const ccomplex_t IT_0402 = m_N_3*IT_0209;
    const ccomplex_t IT_0403 = -IT_0402;
    const ccomplex_t IT_0404 = IT_0401 + IT_0403;
    const ccomplex_t IT_0405 = IT_0218*IT_0404;
    const ccomplex_t IT_0406 = (-0.5)*IT_0400 + -IT_0405;
    const ccomplex_t IT_0407 = 2*IT_0406;
    const ccomplex_t IT_0408 = m_N_1*IT_0140;
    const ccomplex_t IT_0409 = m_N_3*IT_0191;
    const ccomplex_t IT_0410 = IT_0200*(IT_0408 + -IT_0409);
    const ccomplex_t IT_0411 = -IT_0410;
    const ccomplex_t IT_0412 = 2*IT_0411;
    const ccomplex_t IT_0413 = IT_0407 + IT_0412;
    const ccomplex_t IT_0414 = m_N_4*IT_0086;
    const ccomplex_t IT_0415 = m_N_3*IT_0171;
    const ccomplex_t IT_0416 = -IT_0415;
    const ccomplex_t IT_0417 = IT_0414 + IT_0416;
    const ccomplex_t IT_0418 = IT_0180*IT_0417;
    const ccomplex_t IT_0419 = 0.5*IT_0400 + -IT_0418;
    const ccomplex_t IT_0420 = (-2)*IT_0419;
    const ccomplex_t IT_0421 = IT_0413 + IT_0420;
    const ccomplex_t IT_0422 = IT_0300 + IT_0330;
    const ccomplex_t IT_0423 = 2*conj(IT_0406);
    const ccomplex_t IT_0424 = 2*conj(IT_0411);
    const ccomplex_t IT_0425 = IT_0423 + IT_0424;
    const ccomplex_t IT_0426 = (-2)*conj(IT_0419);
    const ccomplex_t IT_0427 = IT_0425 + IT_0426;
    const ccomplex_t IT_0428 = m_C_2*m_N_3;
    const ccomplex_t IT_0429 = pow(s_34, 2);
    const ccomplex_t IT_0430 = IT_0018*IT_0019*IT_0429;
    const ccomplex_t IT_0431 = 0.5*IT_0430;
    const ccomplex_t IT_0432 = 1 + IT_0431;
    const ccomplex_t IT_0433 = IT_0428*IT_0432;
    const ccomplex_t IT_0434 = (-4)*conj(IT_0244);
    const ccomplex_t IT_0435 = 0.25*IT_0434;
    const ccomplex_t IT_0436 = conj(IT_0222) + IT_0435;
    const ccomplex_t IT_0437 = (-4)*IT_0222;
    const ccomplex_t IT_0438 = 0.25*IT_0437;
    const ccomplex_t IT_0439 = IT_0244 + IT_0438;
    const ccomplex_t IT_0440 = (-4)*conj(IT_0222);
    const ccomplex_t IT_0441 = 0.25*IT_0440;
    const ccomplex_t IT_0442 = conj(IT_0244) + IT_0441;
    const ccomplex_t IT_0443 = 4*IT_0222;
    const ccomplex_t IT_0444 = (-4)*IT_0244;
    const ccomplex_t IT_0445 = IT_0443 + IT_0444;
    const ccomplex_t IT_0446 = -IT_0003;
    const ccomplex_t IT_0447 = IT_0015*IT_0446;
    const ccomplex_t IT_0448 = s_34*m_C_2*IT_0018;
    const ccomplex_t IT_0449 = s_23*s_34*IT_0019;
    const ccomplex_t IT_0450 = -IT_0449;
    const ccomplex_t IT_0451 = s_24 + IT_0450;
    const ccomplex_t IT_0452 = IT_0448*IT_0451;
    const ccomplex_t IT_0453 = s_34*m_N_3*IT_0018;
    const ccomplex_t IT_0454 = s_13*s_34*IT_0019;
    const ccomplex_t IT_0455 = -IT_0454;
    const ccomplex_t IT_0456 = s_14 + IT_0455;
    const ccomplex_t IT_0457 = IT_0453*IT_0456;
    const ccomplex_t IT_0458 = s_24*s_34*IT_0018;
    const ccomplex_t IT_0459 = s_12*IT_0094;
    const ccomplex_t IT_0460 = s_13*s_23;
    const ccomplex_t IT_0461 = (-2)*IT_0460;
    const ccomplex_t IT_0462 = IT_0459 + IT_0461;
    const ccomplex_t IT_0463 = IT_0019*IT_0462;
    const ccomplex_t IT_0464 = (s_23 + -IT_0458)*(s_12 + (-0.5)*IT_0463);
    const ccomplex_t IT_0465 = (-2)*IT_0464;
    const ccomplex_t IT_0466 = 2*conj(IT_0017);
    const ccomplex_t IT_0467 = (s_12 + -IT_0022)*(s_12 + (-0.5)*IT_0463);
    const ccomplex_t IT_0468 = (-2)*IT_0467;
    const ccomplex_t IT_0469 = s_14*s_34*IT_0018;
    const ccomplex_t IT_0470 = (s_12 + (-0.5)*IT_0463)*(s_13 + -IT_0469);
    const ccomplex_t IT_0471 = (-2)*IT_0470;
    const ccomplex_t IT_0472 = m_C_2*IT_0317;
    const ccomplex_t IT_0473 = s_14*s_24*IT_0018*IT_0094;
    const ccomplex_t IT_0474 = -IT_0473;
    const ccomplex_t IT_0475 = IT_0459 + IT_0474;
    const ccomplex_t IT_0476 = IT_0019*IT_0475;
    const ccomplex_t IT_0477 = -IT_0022 + (-0.25)*IT_0476;
    const ccomplex_t IT_0478 = s_12 + IT_0477;
    const ccomplex_t IT_0479 = m_C_2*IT_0478;
    const ccomplex_t IT_0480 = (-4)*IT_0479;
    const ccomplex_t IT_0481 = pow(s_14, 2);
    const ccomplex_t IT_0482 = IT_0018*IT_0094*IT_0481;
    const ccomplex_t IT_0483 = IT_0011*IT_0094;
    const ccomplex_t IT_0484 = -IT_0483;
    const ccomplex_t IT_0485 = IT_0482 + IT_0484;
    const ccomplex_t IT_0486 = IT_0019*IT_0485;
    const ccomplex_t IT_0487 = IT_0018*IT_0481;
    const ccomplex_t IT_0488 = -IT_0011;
    const ccomplex_t IT_0489 = IT_0487 + IT_0488;
    const ccomplex_t IT_0490 = (-4)*IT_0489;
    const ccomplex_t IT_0491 = IT_0486 + IT_0490;
    const ccomplex_t IT_0492 = m_N_3*IT_0491;
    const ccomplex_t IT_0493 = -IT_0132 + (-2)*IT_0281 + (-2)*IT_0292;
    const ccomplex_t IT_0494 = pow(s_13, 2);
    const ccomplex_t IT_0495 = IT_0019*IT_0494;
    const ccomplex_t IT_0496 = IT_0487 + IT_0495;
    const ccomplex_t IT_0497 = s_13*s_14*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0498 = -IT_0011 + -IT_0497;
    const ccomplex_t IT_0499 = IT_0496 + IT_0498;
    const ccomplex_t IT_0500 = IT_0428*IT_0499;
    const ccomplex_t IT_0501 = -IT_0182 + 2*IT_0317;
    const ccomplex_t IT_0502 = s_13*s_23*IT_0018*IT_0019*IT_0481;
    const ccomplex_t IT_0503 = s_14*s_24*IT_0018*IT_0019*IT_0494;
    const ccomplex_t IT_0504 = s_12*IT_0011;
    const ccomplex_t IT_0505 = s_14*s_24*IT_0011*IT_0018;
    const ccomplex_t IT_0506 = -IT_0505;
    const ccomplex_t IT_0507 = s_13*s_23*IT_0011*IT_0019;
    const ccomplex_t IT_0508 = -IT_0507;
    const ccomplex_t IT_0509 = s_12*s_13*s_14*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0510 = -IT_0509;
    const ccomplex_t IT_0511 = IT_0502 + IT_0503 + IT_0504 + IT_0506 + IT_0508
       + IT_0510;
    const ccomplex_t IT_0512 = IT_0018*IT_0019*IT_0068*IT_0481;
    const ccomplex_t IT_0513 = s_13*s_14*s_23*s_24*IT_0018*IT_0019;
    const ccomplex_t IT_0514 = IT_0011*IT_0012;
    const ccomplex_t IT_0515 = IT_0012*IT_0018*IT_0481;
    const ccomplex_t IT_0516 = -IT_0515;
    const ccomplex_t IT_0517 = IT_0011*IT_0019*IT_0068;
    const ccomplex_t IT_0518 = -IT_0517;
    const ccomplex_t IT_0519 = s_12*s_14*s_23*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0520 = -IT_0519;
    const ccomplex_t IT_0521 = IT_0512 + IT_0513 + IT_0514 + IT_0516 + IT_0518
       + IT_0520;
    const ccomplex_t IT_0522 = 2*IT_0003;
    const ccomplex_t IT_0523 = IT_0015*IT_0522;
    const ccomplex_t IT_0524 = s_12*s_13;
    const ccomplex_t IT_0525 = s_23*IT_0011;
    const ccomplex_t IT_0526 = s_13*s_14*s_23*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0527 = 2*IT_0526;
    const ccomplex_t IT_0528 = s_13*s_14*s_24*IT_0018;
    const ccomplex_t IT_0529 = -IT_0528;
    const ccomplex_t IT_0530 = s_23*IT_0018*IT_0481;
    const ccomplex_t IT_0531 = -IT_0530;
    const ccomplex_t IT_0532 = s_23*IT_0019*IT_0494;
    const ccomplex_t IT_0533 = (-2)*IT_0532;
    const ccomplex_t IT_0534 = IT_0524 + IT_0525 + IT_0527 + IT_0529 + IT_0531
       + IT_0533;
    const ccomplex_t IT_0535 = s_14*s_23*s_34*IT_0018;
    const ccomplex_t IT_0536 = (-2)*IT_0460 + -IT_0473 + 2*IT_0535;
    const ccomplex_t IT_0537 = IT_0459 + IT_0536;
    const ccomplex_t IT_0538 = IT_0019*IT_0537;
    const ccomplex_t IT_0539 = -IT_0022 + (-0.5)*IT_0538;
    const ccomplex_t IT_0540 = s_12 + IT_0539;
    const ccomplex_t IT_0541 = m_C_2*IT_0540;
    const ccomplex_t IT_0542 = (-2)*IT_0541;
    const ccomplex_t IT_0543 = IT_0369 + IT_0374 + IT_0379;
    const ccomplex_t IT_0544 = s_14*s_34*IT_0018*IT_0094;
    const ccomplex_t IT_0545 = -IT_0544;
    const ccomplex_t IT_0546 = IT_0331 + IT_0545;
    const ccomplex_t IT_0547 = IT_0019*IT_0546;
    const ccomplex_t IT_0548 = -IT_0469 + (-0.25)*IT_0547;
    const ccomplex_t IT_0549 = s_13 + IT_0548;
    const ccomplex_t IT_0550 = IT_0428*IT_0549;
    const ccomplex_t IT_0551 = (-4)*IT_0550;
    const ccomplex_t IT_0552 = IT_0466*(IT_0017*IT_0468 + IT_0447*IT_0471 +
       IT_0472*IT_0480 + IT_0419*IT_0492 + IT_0493*IT_0500 + IT_0501*IT_0511 +
       IT_0183*IT_0521 + IT_0523*IT_0534 + IT_0542*IT_0543 + IT_0298*IT_0551);
    const ccomplex_t IT_0553 = IT_0018*IT_0019*IT_0070*IT_0494;
    const ccomplex_t IT_0554 = IT_0012*IT_0019*IT_0494;
    const ccomplex_t IT_0555 = -IT_0554;
    const ccomplex_t IT_0556 = IT_0011*IT_0018*IT_0070;
    const ccomplex_t IT_0557 = -IT_0556;
    const ccomplex_t IT_0558 = s_12*s_13*s_24*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0559 = -IT_0558;
    const ccomplex_t IT_0560 = IT_0513 + IT_0514 + IT_0553 + IT_0555 + IT_0557
       + IT_0559;
    const ccomplex_t IT_0561 = 2*IT_0160;
    const ccomplex_t IT_0562 = conj(IT_0493)*IT_0561;
    const ccomplex_t IT_0563 = 2*conj(IT_0160);
    const ccomplex_t IT_0564 = IT_0493*IT_0563;
    const ccomplex_t IT_0565 = s_12*s_23;
    const ccomplex_t IT_0566 = s_13*IT_0012;
    const ccomplex_t IT_0567 = s_13*s_23*s_24*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0568 = 2*IT_0567;
    const ccomplex_t IT_0569 = s_13*IT_0018*IT_0070;
    const ccomplex_t IT_0570 = -IT_0569;
    const ccomplex_t IT_0571 = s_14*s_23*s_24*IT_0018;
    const ccomplex_t IT_0572 = -IT_0571;
    const ccomplex_t IT_0573 = s_13*IT_0019*IT_0068;
    const ccomplex_t IT_0574 = (-2)*IT_0573;
    const ccomplex_t IT_0575 = IT_0565 + IT_0566 + IT_0568 + IT_0570 + IT_0572
       + IT_0574;
    const ccomplex_t IT_0576 = IT_0063*IT_0522;
    const ccomplex_t IT_0577 = IT_0561*conj(IT_0576);
    const ccomplex_t IT_0578 = IT_0563*IT_0576;
    const ccomplex_t IT_0579 = conj(IT_0300)*IT_0561;
    const ccomplex_t IT_0580 = conj(IT_0318)*IT_0561;
    const ccomplex_t IT_0581 = conj(IT_0329)*IT_0561;
    const ccomplex_t IT_0582 = conj(IT_0330)*IT_0561;
    const ccomplex_t IT_0583 = IT_0300*IT_0563;
    const ccomplex_t IT_0584 = IT_0318*IT_0563;
    const ccomplex_t IT_0585 = IT_0329*IT_0563;
    const ccomplex_t IT_0586 = IT_0330*IT_0563;
    const ccomplex_t IT_0587 = s_24*s_34*IT_0018*IT_0094;
    const ccomplex_t IT_0588 = -IT_0587;
    const ccomplex_t IT_0589 = IT_0250 + IT_0588;
    const ccomplex_t IT_0590 = IT_0019*IT_0589;
    const ccomplex_t IT_0591 = -IT_0458 + (-0.25)*IT_0590;
    const ccomplex_t IT_0592 = s_23 + IT_0591;
    const ccomplex_t IT_0593 = IT_0428*IT_0592;
    const ccomplex_t IT_0594 = (-4)*IT_0593;
    const ccomplex_t IT_0595 = s_14*s_24*IT_0018*IT_0019*IT_0068;
    const ccomplex_t IT_0596 = s_13*s_23*IT_0018*IT_0019*IT_0070;
    const ccomplex_t IT_0597 = s_12*IT_0012;
    const ccomplex_t IT_0598 = s_14*s_24*IT_0012*IT_0018;
    const ccomplex_t IT_0599 = -IT_0598;
    const ccomplex_t IT_0600 = s_13*s_23*IT_0012*IT_0019;
    const ccomplex_t IT_0601 = -IT_0600;
    const ccomplex_t IT_0602 = s_12*s_23*s_24*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0603 = -IT_0602;
    const ccomplex_t IT_0604 = IT_0595 + IT_0596 + IT_0597 + IT_0599 + IT_0601
       + IT_0603;
    const ccomplex_t IT_0605 = conj(IT_0154)*IT_0561;
    const ccomplex_t IT_0606 = IT_0154*IT_0563;
    const ccomplex_t IT_0607 = 2*IT_0111;
    const ccomplex_t IT_0608 = conj(IT_0160)*IT_0607;
    const ccomplex_t IT_0609 = 2*conj(IT_0111);
    const ccomplex_t IT_0610 = IT_0160*IT_0609;
    const ccomplex_t IT_0611 = 2*IT_0017;
    const ccomplex_t IT_0612 = IT_0018*IT_0429;
    const ccomplex_t IT_0613 = (s_12 + (-0.5)*IT_0463)*(IT_0094 + -IT_0612);
    const ccomplex_t IT_0614 = (-2)*IT_0613;
    const ccomplex_t IT_0615 = 2*IT_0447;
    const ccomplex_t IT_0616 = s_14*s_23*s_34*IT_0018*IT_0019;
    const ccomplex_t IT_0617 = s_12 + IT_0023 + IT_0616;
    const ccomplex_t IT_0618 = m_N_3*IT_0617;
    const ccomplex_t IT_0619 = s_14*s_34*IT_0018*IT_0019*IT_0068;
    const ccomplex_t IT_0620 = IT_0565 + IT_0566 + IT_0567;
    const ccomplex_t IT_0621 = s_12*s_23*IT_0018*IT_0019*IT_0429;
    const ccomplex_t IT_0622 = -IT_0621;
    const ccomplex_t IT_0623 = s_14*s_34*IT_0012*IT_0018;
    const ccomplex_t IT_0624 = -IT_0623;
    const ccomplex_t IT_0625 = IT_0574 + IT_0619 + IT_0620 + IT_0622 + IT_0624;
    const ccomplex_t IT_0626 = 2*conj(IT_0447);
    const ccomplex_t IT_0627 = s_12*IT_0030;
    const ccomplex_t IT_0628 = s_14*s_24;
    const ccomplex_t IT_0629 = (-2)*IT_0628;
    const ccomplex_t IT_0630 = IT_0627 + IT_0629;
    const ccomplex_t IT_0631 = IT_0018*IT_0630;
    const ccomplex_t IT_0632 = (s_12 + -IT_0021)*(s_12 + (-0.5)*IT_0631);
    const ccomplex_t IT_0633 = (-2)*IT_0632;
    const ccomplex_t IT_0634 = 2*conj(IT_0501);
    const ccomplex_t IT_0635 = s_13*s_23*IT_0019*IT_0030;
    const ccomplex_t IT_0636 = s_14*s_23*s_34*IT_0019;
    const ccomplex_t IT_0637 = (-2)*IT_0628 + -IT_0635 + 2*IT_0636;
    const ccomplex_t IT_0638 = IT_0627 + IT_0637;
    const ccomplex_t IT_0639 = IT_0018*IT_0638;
    const ccomplex_t IT_0640 = -IT_0021 + (-0.5)*IT_0639;
    const ccomplex_t IT_0641 = s_12 + IT_0640;
    const ccomplex_t IT_0642 = m_N_3*IT_0641;
    const ccomplex_t IT_0643 = (-2)*IT_0642;
    const ccomplex_t IT_0644 = conj(IT_0406) + conj(IT_0411);
    const ccomplex_t IT_0645 = (-2)*IT_0183;
    const ccomplex_t IT_0646 = (-2)*IT_0221;
    const ccomplex_t IT_0647 = IT_0645 + IT_0646;
    const ccomplex_t IT_0648 = IT_0406 + IT_0411;
    const ccomplex_t IT_0649 = (-2)*conj(IT_0183);
    const ccomplex_t IT_0650 = (-2)*conj(IT_0221);
    const ccomplex_t IT_0651 = IT_0649 + IT_0650;
    const ccomplex_t IT_0652 = conj(IT_0369) + conj(IT_0374) + conj(IT_0379);
    const ccomplex_t IT_0653 = (-2)*IT_0111;
    const ccomplex_t IT_0654 = (-2)*IT_0154;
    const ccomplex_t IT_0655 = IT_0653 + IT_0654;
    const ccomplex_t IT_0656 = (-2)*conj(IT_0111);
    const ccomplex_t IT_0657 = (-2)*conj(IT_0154);
    const ccomplex_t IT_0658 = IT_0656 + IT_0657;
    const ccomplex_t IT_0659 = 2*IT_0183;
    const ccomplex_t IT_0660 = 2*conj(IT_0183);
    const ccomplex_t IT_0661 = 2*IT_0390;
    const ccomplex_t IT_0662 = IT_0111 + IT_0154;
    const ccomplex_t IT_0663 = 2*conj(IT_0390);
    const ccomplex_t IT_0664 = 2*IT_0419;
    const ccomplex_t IT_0665 = 2*conj(IT_0419);
    const ccomplex_t IT_0666 = IT_0643*(IT_0644*IT_0647 + IT_0648*IT_0651 +
       IT_0652*IT_0655 + IT_0543*IT_0658 + conj(IT_0419)*IT_0659 + IT_0419
      *IT_0660 + IT_0155*IT_0661 + IT_0662*IT_0663 + conj(IT_0221)*IT_0664 +
       IT_0221*IT_0665);
    const ccomplex_t IT_0667 = -IT_0566;
    const ccomplex_t IT_0668 = IT_0565 + IT_0667;
    const ccomplex_t IT_0669 = s_12*s_23*IT_0030;
    const ccomplex_t IT_0670 = s_14*s_23*s_24;
    const ccomplex_t IT_0671 = s_13*IT_0012*IT_0030;
    const ccomplex_t IT_0672 = s_14*s_34*IT_0012;
    const ccomplex_t IT_0673 = (-2)*IT_0670 + -IT_0671 + 2*IT_0672;
    const ccomplex_t IT_0674 = IT_0669 + IT_0673;
    const ccomplex_t IT_0675 = IT_0018*IT_0674;
    const ccomplex_t IT_0676 = (-0.5)*IT_0675;
    const ccomplex_t IT_0677 = IT_0668 + IT_0676;
    const ccomplex_t IT_0678 = 4*IT_0111;
    const ccomplex_t IT_0679 = 4*IT_0154;
    const ccomplex_t IT_0680 = 4*conj(IT_0111);
    const ccomplex_t IT_0681 = 4*conj(IT_0154);
    const ccomplex_t IT_0682 = 4*IT_0183;
    const ccomplex_t IT_0683 = 4*IT_0221;
    const ccomplex_t IT_0684 = 4*conj(IT_0183);
    const ccomplex_t IT_0685 = 4*conj(IT_0221);
    const ccomplex_t IT_0686 = IT_0030*IT_0094;
    const ccomplex_t IT_0687 = IT_0356*IT_0686;
    const ccomplex_t IT_0688 = 0.125*IT_0687;
    const ccomplex_t IT_0689 = 1 + IT_0688;
    const ccomplex_t IT_0690 = IT_0428*IT_0689;
    const ccomplex_t IT_0691 = m_C_2*IT_0328;
    const ccomplex_t IT_0692 = -IT_0691;
    const ccomplex_t IT_0693 = IT_0010*IT_0289;
    const ccomplex_t IT_0694 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0695 = IT_0282*IT_0693*IT_0694;
    const ccomplex_t IT_0696 = -IT_0695;
    const ccomplex_t IT_0697 = IT_0277*IT_0307;
    const ccomplex_t IT_0698 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0699 = IT_0279*IT_0697*IT_0698;
    const ccomplex_t IT_0700 = -IT_0699;
    const ccomplex_t IT_0701 = IT_0061*IT_0325;
    const ccomplex_t IT_0702 = IT_0282*IT_0694*IT_0701;
    const ccomplex_t IT_0703 = IT_0270*IT_0314;
    const ccomplex_t IT_0704 = IT_0279*IT_0698*IT_0703;
    const ccomplex_t IT_0705 = m_C_2*IT_0281;
    const ccomplex_t IT_0706 = m_C_2*IT_0292;
    const ccomplex_t IT_0707 = IT_0702 + IT_0704 + IT_0705 + IT_0706;
    const ccomplex_t IT_0708 = 16*IT_0419;
    const ccomplex_t IT_0709 = 16*conj(IT_0419);
    const ccomplex_t IT_0710 = (-16)*IT_0369;
    const ccomplex_t IT_0711 = (-16)*IT_0374;
    const ccomplex_t IT_0712 = (-16)*IT_0379;
    const ccomplex_t IT_0713 = 16*conj(IT_0390);
    const ccomplex_t IT_0714 = (-16)*conj(IT_0369);
    const ccomplex_t IT_0715 = (-16)*conj(IT_0374);
    const ccomplex_t IT_0716 = (-16)*conj(IT_0379);
    const ccomplex_t IT_0717 = (-16)*conj(IT_0390)*IT_0419 + (-16)*IT_0390
      *conj(IT_0419) + (16*conj(IT_0472) + (-16)*conj(IT_0692) + (-16)*conj
      (IT_0696) + (-16)*conj(IT_0700))*IT_0707 + (16*IT_0472 + (-16)*IT_0692 + (
      -16)*IT_0696 + (-16)*IT_0700)*conj(IT_0707) + IT_0652*IT_0708 + IT_0543
      *IT_0709 + IT_0644*(16*IT_0390 + IT_0710 + IT_0711 + IT_0712) + IT_0648*
      (IT_0713 + IT_0714 + IT_0715 + IT_0716);
    const ccomplex_t IT_0718 = -IT_0469;
    const ccomplex_t IT_0719 = s_13 + IT_0718;
    const ccomplex_t IT_0720 = IT_0428*IT_0719;
    const ccomplex_t IT_0721 = -IT_0458;
    const ccomplex_t IT_0722 = s_23 + IT_0721;
    const ccomplex_t IT_0723 = IT_0428*IT_0722;
    const ccomplex_t IT_0724 = (-6)*IT_0447;
    const ccomplex_t IT_0725 = conj(IT_0160)*IT_0724;
    const ccomplex_t IT_0726 = (-6)*conj(IT_0447);
    const ccomplex_t IT_0727 = IT_0160*IT_0726;
    const ccomplex_t IT_0728 = s_12*IT_0432;
    const ccomplex_t IT_0729 = -conj(IT_0065);
    const ccomplex_t IT_0730 = conj(IT_0028) + IT_0729;
    const ccomplex_t IT_0731 = (-4)*conj(IT_0028);
    const ccomplex_t IT_0732 = 0.25*IT_0731;
    const ccomplex_t IT_0733 = conj(IT_0065) + IT_0732;
    const ccomplex_t IT_0734 = IT_0018*IT_0070*IT_0094;
    const ccomplex_t IT_0735 = s_23*s_24*s_34*IT_0018;
    const ccomplex_t IT_0736 = IT_0012*IT_0094;
    const ccomplex_t IT_0737 = 0.5*IT_0734 + -IT_0735 + (-0.5)*IT_0736;
    const ccomplex_t IT_0738 = IT_0068 + IT_0737;
    const ccomplex_t IT_0739 = IT_0019*IT_0738;
    const ccomplex_t IT_0740 = -IT_0071;
    const ccomplex_t IT_0741 = IT_0012 + IT_0739 + IT_0740;
    const ccomplex_t IT_0742 = m_C_2*IT_0741;
    const ccomplex_t IT_0743 = (-4)*IT_0017;
    const ccomplex_t IT_0744 = (-4)*conj(IT_0017);
    const ccomplex_t IT_0745 = (-4)*IT_0160;
    const ccomplex_t IT_0746 = (-4)*conj(IT_0160);
    const ccomplex_t IT_0747 = -(IT_0012 + -IT_0071)*(s_12 + (-0.5)*IT_0463);
    const ccomplex_t IT_0748 = (-2)*IT_0747;
    const ccomplex_t IT_0749 = (-2)*IT_0017;
    const ccomplex_t IT_0750 = m_N_3*IT_0478;
    const ccomplex_t IT_0751 = (-4)*IT_0750;
    const ccomplex_t IT_0752 = (-2)*conj(IT_0017);
    const ccomplex_t IT_0753 = m_C_2*IT_0617;
    const ccomplex_t IT_0754 = m_N_3*IT_0499;
    const ccomplex_t IT_0755 = -(s_12 + (-0.5)*IT_0463)*(IT_0011 + -IT_0487);
    const ccomplex_t IT_0756 = (-2)*IT_0755;
    const ccomplex_t IT_0757 = (-2)*IT_0160;
    const ccomplex_t IT_0758 = (-2)*conj(IT_0160);
    const ccomplex_t IT_0759 = 2*IT_0501;
    const ccomplex_t IT_0760 = -(IT_0012 + -IT_0069)*(s_12 + (-0.5)*IT_0631);
    const ccomplex_t IT_0761 = (-2)*IT_0760;
    const ccomplex_t IT_0762 = m_N_3*IT_0549;
    const ccomplex_t IT_0763 = (-4)*IT_0762;
    const ccomplex_t IT_0764 = IT_0063*IT_0446;
    const ccomplex_t IT_0765 = 2*IT_0764;
    const ccomplex_t IT_0766 = 2*conj(IT_0764);
    const ccomplex_t IT_0767 = (-2)*IT_0447;
    const ccomplex_t IT_0768 = (-2)*conj(IT_0447);
    const ccomplex_t IT_0769 = IT_0250 + IT_0587;
    const ccomplex_t IT_0770 = s_23*IT_0018*IT_0429;
    const ccomplex_t IT_0771 = (-2)*IT_0770;
    const ccomplex_t IT_0772 = IT_0769 + IT_0771;
    const ccomplex_t IT_0773 = IT_0019*IT_0772;
    const ccomplex_t IT_0774 = -IT_0458 + 0.5*IT_0773;
    const ccomplex_t IT_0775 = s_23 + IT_0774;
    const ccomplex_t IT_0776 = m_C_2*IT_0775;
    const ccomplex_t IT_0777 = 2*IT_0776;
    const ccomplex_t IT_0778 = (-2)*IT_0764;
    const ccomplex_t IT_0779 = (-2)*conj(IT_0764);
    const ccomplex_t IT_0780 = m_C_2*IT_0024;
    const ccomplex_t IT_0781 = s_24*s_34*IT_0018*IT_0019*IT_0494;
    const ccomplex_t IT_0782 = IT_0524 + IT_0525 + IT_0526;
    const ccomplex_t IT_0783 = s_12*s_13*IT_0018*IT_0019*IT_0429;
    const ccomplex_t IT_0784 = -IT_0783;
    const ccomplex_t IT_0785 = s_24*s_34*IT_0011*IT_0018;
    const ccomplex_t IT_0786 = -IT_0785;
    const ccomplex_t IT_0787 = IT_0533 + IT_0781 + IT_0782 + IT_0784 + IT_0786;
    const ccomplex_t IT_0788 = -(IT_0011 + -IT_0495)*(s_12 + (-0.5)*IT_0631);
    const ccomplex_t IT_0789 = (-2)*IT_0788;
    const ccomplex_t IT_0790 = (-2)*IT_0501;
    const ccomplex_t IT_0791 = IT_0019*IT_0459;
    const ccomplex_t IT_0792 = (-0.25)*IT_0791;
    const ccomplex_t IT_0793 = IT_0018*IT_0627;
    const ccomplex_t IT_0794 = (-0.25)*IT_0793;
    const ccomplex_t IT_0795 = s_12 + IT_0792 + IT_0794;
    const ccomplex_t IT_0796 = s_14*s_23;
    const ccomplex_t IT_0797 = s_13*s_24;
    const ccomplex_t IT_0798 = s_12*s_34;
    const ccomplex_t IT_0799 = -IT_0797 + -IT_0798;
    const ccomplex_t IT_0800 = IT_0796 + IT_0799;
    const ccomplex_t IT_0801 = IT_0252*IT_0800;
    const ccomplex_t IT_0802 = (-0.25)*IT_0801;
    const ccomplex_t IT_0803 = IT_0795 + IT_0802;
    const ccomplex_t IT_0804 = 8*IT_0222;
    const ccomplex_t IT_0805 = (-8)*IT_0244;
    const ccomplex_t IT_0806 = IT_0804 + IT_0805;
    const ccomplex_t IT_0807 = IT_0019*IT_0030*IT_0068;
    const ccomplex_t IT_0808 = IT_0012*IT_0030;
    const ccomplex_t IT_0809 = s_23*s_24*s_34*IT_0019;
    const ccomplex_t IT_0810 = 2*IT_0070 + -IT_0808 + (-2)*IT_0809;
    const ccomplex_t IT_0811 = IT_0807 + IT_0810;
    const ccomplex_t IT_0812 = IT_0018*IT_0811;
    const ccomplex_t IT_0813 = -IT_0012 + (-0.5)*IT_0812;
    const ccomplex_t IT_0814 = IT_0069 + IT_0813;
    const ccomplex_t IT_0815 = m_C_2*IT_0814;
    const ccomplex_t IT_0816 = (-2)*IT_0815;
    const ccomplex_t IT_0817 = m_C_2*IT_0592;
    const ccomplex_t IT_0818 = (-4)*IT_0817;
    const ccomplex_t IT_0819 = IT_0018*IT_0255;
    const ccomplex_t IT_0820 = (-0.25)*IT_0251 + 0.25*IT_0254 + (-0.25)*IT_0819;
    const ccomplex_t IT_0821 = s_23 + IT_0820;
    const ccomplex_t IT_0822 = m_C_2*IT_0821;
    const ccomplex_t IT_0823 = (-4)*IT_0822;
    const ccomplex_t IT_0824 = 2*IT_0523;
    const ccomplex_t IT_0825 = (-0.5)*IT_0631;
    const ccomplex_t IT_0826 = (-0.5)*IT_0463;
    const ccomplex_t IT_0827 = s_12*IT_0030*IT_0094;
    const ccomplex_t IT_0828 = s_13*s_23*IT_0030;
    const ccomplex_t IT_0829 = s_13*s_24*s_34;
    const ccomplex_t IT_0830 = s_14*s_24*IT_0094;
    const ccomplex_t IT_0831 = (-2)*IT_0828 + 4*IT_0829 + (-2)*IT_0830;
    const ccomplex_t IT_0832 = IT_0827 + IT_0831;
    const ccomplex_t IT_0833 = IT_0356*IT_0832;
    const ccomplex_t IT_0834 = 0.25*IT_0833;
    const ccomplex_t IT_0835 = s_12 + IT_0825 + IT_0826 + IT_0834;
    const ccomplex_t IT_0836 = 8*IT_0472;
    const ccomplex_t IT_0837 = -IT_0736;
    const ccomplex_t IT_0838 = IT_0734 + IT_0837;
    const ccomplex_t IT_0839 = IT_0019*IT_0838;
    const ccomplex_t IT_0840 = -IT_0012;
    const ccomplex_t IT_0841 = IT_0071 + IT_0840;
    const ccomplex_t IT_0842 = (-4)*IT_0841;
    const ccomplex_t IT_0843 = IT_0839 + IT_0842;
    const ccomplex_t IT_0844 = m_C_2*IT_0843;
    const ccomplex_t IT_0845 = s_13*s_24*s_34*IT_0019;
    const ccomplex_t IT_0846 = (-2)*IT_0628 + -IT_0635 + 2*IT_0845;
    const ccomplex_t IT_0847 = IT_0627 + IT_0846;
    const ccomplex_t IT_0848 = IT_0018*IT_0847;
    const ccomplex_t IT_0849 = -IT_0021 + (-0.5)*IT_0848;
    const ccomplex_t IT_0850 = s_12 + IT_0849;
    const ccomplex_t IT_0851 = m_C_2*IT_0850;
    const ccomplex_t IT_0852 = (-2)*IT_0851;
    const ccomplex_t IT_0853 = 8*conj(IT_0222);
    const ccomplex_t IT_0854 = (-8)*conj(IT_0244);
    const ccomplex_t IT_0855 = IT_0853 + IT_0854;
    const ccomplex_t IT_0856 = 2*conj(IT_0523);
    const ccomplex_t IT_0857 = (-2)*conj(IT_0501);
    const ccomplex_t IT_0858 = -IT_0635;
    const ccomplex_t IT_0859 = IT_0627 + IT_0858;
    const ccomplex_t IT_0860 = IT_0018*IT_0859;
    const ccomplex_t IT_0861 = -IT_0021 + (-0.25)*IT_0860;
    const ccomplex_t IT_0862 = s_12 + IT_0861;
    const ccomplex_t IT_0863 = m_C_2*IT_0862;
    const ccomplex_t IT_0864 = (-4)*IT_0863;
    const ccomplex_t IT_0865 = 2*IT_0328;
    const ccomplex_t IT_0866 = (-2)*IT_0865;
    const ccomplex_t IT_0867 = IT_0790 + IT_0866;
    const ccomplex_t IT_0868 = (-2)*conj(IT_0865);
    const ccomplex_t IT_0869 = IT_0857 + IT_0868;
    const ccomplex_t IT_0870 = (-2)*IT_0493;
    const ccomplex_t IT_0871 = (-2)*conj(IT_0493);
    const ccomplex_t IT_0872 = IT_0864*(IT_0390*IT_0634 + conj(IT_0493)
      *IT_0664 + IT_0493*IT_0665 + conj(IT_0390)*IT_0759 + IT_0663*IT_0865 +
       IT_0661*conj(IT_0865) + IT_0652*IT_0867 + IT_0543*IT_0869 + IT_0644
      *IT_0870 + IT_0648*IT_0871);
    const ccomplex_t IT_0873 = IT_0018*IT_0335;
    const ccomplex_t IT_0874 = (-0.25)*IT_0332 + 0.25*IT_0334 + (-0.25)*IT_0873;
    const ccomplex_t IT_0875 = s_13 + IT_0874;
    const ccomplex_t IT_0876 = m_N_3*IT_0875;
    const ccomplex_t IT_0877 = (-4)*IT_0876;
    const ccomplex_t IT_0878 = 2*IT_0576;
    const ccomplex_t IT_0879 = 2*conj(IT_0576);
    const ccomplex_t IT_0880 = (-2)*IT_0523;
    const ccomplex_t IT_0881 = (-2)*conj(IT_0523);
    const ccomplex_t IT_0882 = m_N_3*IT_0862;
    const ccomplex_t IT_0883 = (-4)*IT_0882;
    const ccomplex_t IT_0884 = 2*IT_0221;
    const ccomplex_t IT_0885 = IT_0659 + IT_0884;
    const ccomplex_t IT_0886 = 8*IT_0065;
    const ccomplex_t IT_0887 = (-8)*IT_0028;
    const ccomplex_t IT_0888 = IT_0886 + IT_0887;
    const ccomplex_t IT_0889 = IT_0019*IT_0030*IT_0494;
    const ccomplex_t IT_0890 = IT_0011*IT_0030;
    const ccomplex_t IT_0891 = -IT_0890;
    const ccomplex_t IT_0892 = IT_0889 + IT_0891;
    const ccomplex_t IT_0893 = IT_0018*IT_0892;
    const ccomplex_t IT_0894 = IT_0488 + IT_0495;
    const ccomplex_t IT_0895 = (-4)*IT_0894;
    const ccomplex_t IT_0896 = IT_0893 + IT_0895;
    const ccomplex_t IT_0897 = m_N_3*IT_0896;
    const ccomplex_t IT_0898 = IT_0331 + IT_0544;
    const ccomplex_t IT_0899 = s_13*IT_0018*IT_0429;
    const ccomplex_t IT_0900 = (-2)*IT_0899;
    const ccomplex_t IT_0901 = IT_0898 + IT_0900;
    const ccomplex_t IT_0902 = IT_0019*IT_0901;
    const ccomplex_t IT_0903 = -IT_0469 + 0.5*IT_0902;
    const ccomplex_t IT_0904 = s_13 + IT_0903;
    const ccomplex_t IT_0905 = m_N_3*IT_0904;
    const ccomplex_t IT_0906 = 2*IT_0905;
    const ccomplex_t IT_0907 = 2*IT_0493;
    const ccomplex_t IT_0908 = s_13*s_34;
    const ccomplex_t IT_0909 = (-2)*IT_0908;
    const ccomplex_t IT_0910 = IT_0333 + IT_0909;
    const ccomplex_t IT_0911 = IT_0252*IT_0910;
    const ccomplex_t IT_0912 = (-0.25)*IT_0332 + (-0.25)*IT_0873 + (-0.25)
      *IT_0911;
    const ccomplex_t IT_0913 = s_13 + IT_0912;
    const ccomplex_t IT_0914 = m_N_3*IT_0913;
    const ccomplex_t IT_0915 = (-4)*IT_0914;
    const ccomplex_t IT_0916 = s_13*s_14*s_34*IT_0018;
    const ccomplex_t IT_0917 = 0.5*IT_0482 + (-0.5)*IT_0483 + -IT_0916;
    const ccomplex_t IT_0918 = IT_0494 + IT_0917;
    const ccomplex_t IT_0919 = IT_0019*IT_0918;
    const ccomplex_t IT_0920 = -IT_0487;
    const ccomplex_t IT_0921 = IT_0011 + IT_0919 + IT_0920;
    const ccomplex_t IT_0922 = m_N_3*IT_0921;
    const ccomplex_t IT_0923 = 4*IT_0017;
    const ccomplex_t IT_0924 = 8*IT_0707;
    const ccomplex_t IT_0925 = s_13*s_24*s_34*IT_0018;
    const ccomplex_t IT_0926 = (-2)*IT_0460 + -IT_0473 + 2*IT_0925;
    const ccomplex_t IT_0927 = IT_0459 + IT_0926;
    const ccomplex_t IT_0928 = IT_0019*IT_0927;
    const ccomplex_t IT_0929 = -IT_0022 + (-0.5)*IT_0928;
    const ccomplex_t IT_0930 = s_12 + IT_0929;
    const ccomplex_t IT_0931 = m_N_3*IT_0930;
    const ccomplex_t IT_0932 = (-2)*IT_0931;
    const ccomplex_t IT_0933 = (-2)*IT_0576;
    const ccomplex_t IT_0934 = IT_0480*IT_0757 + IT_0655*IT_0816 + IT_0778
      *IT_0818 + IT_0561*IT_0844 + IT_0883*IT_0885 + IT_0803*IT_0888 + IT_0867
      *IT_0897 + IT_0615*IT_0906 + IT_0852*IT_0907 + IT_0824*IT_0915 + IT_0922
      *IT_0923 + IT_0835*IT_0924 + IT_0749*IT_0932 + IT_0823*IT_0933;
    const ccomplex_t IT_0935 = 2*conj(IT_0221);
    const ccomplex_t IT_0936 = IT_0660 + IT_0935;
    const ccomplex_t IT_0937 = 8*conj(IT_0065);
    const ccomplex_t IT_0938 = (-8)*conj(IT_0028);
    const ccomplex_t IT_0939 = IT_0937 + IT_0938;
    const ccomplex_t IT_0940 = 2*conj(IT_0493);
    const ccomplex_t IT_0941 = 4*conj(IT_0017);
    const ccomplex_t IT_0942 = (-2)*conj(IT_0576);
    const ccomplex_t IT_0943 = IT_0480*IT_0758 + IT_0658*IT_0816 + IT_0779
      *IT_0818 + IT_0563*IT_0844 + IT_0869*IT_0897 + IT_0626*IT_0906 + IT_0856
      *IT_0915 + IT_0752*IT_0932 + IT_0883*IT_0936 + IT_0803*IT_0939 + IT_0852
      *IT_0940 + IT_0922*IT_0941 + IT_0823*IT_0942;
    const ccomplex_t IT_0944 = s_23*s_34;
    const ccomplex_t IT_0945 = (-2)*IT_0944;
    const ccomplex_t IT_0946 = IT_0253 + IT_0945;
    const ccomplex_t IT_0947 = IT_0252*IT_0946;
    const ccomplex_t IT_0948 = (-0.25)*IT_0251 + (-0.25)*IT_0819 + (-0.25)
      *IT_0947;
    const ccomplex_t IT_0949 = s_23 + IT_0948;
    const ccomplex_t IT_0950 = m_C_2*IT_0949;
    const ccomplex_t IT_0951 = (-4)*IT_0950;
    const ccomplex_t IT_0952 = IT_0018*IT_0094*IT_0429;
    const ccomplex_t IT_0953 = -IT_0952;
    const ccomplex_t IT_0954 = pow(m_Z, 4);
    const ccomplex_t IT_0955 = IT_0953 + IT_0954;
    const ccomplex_t IT_0956 = IT_0019*IT_0955;
    const ccomplex_t IT_0957 = -IT_0094 + 0.25*IT_0956;
    const ccomplex_t IT_0958 = IT_0612 + IT_0957;
    const ccomplex_t IT_0959 = IT_0428*IT_0958;
    const ccomplex_t IT_0960 = 4*IT_0959;
    const ccomplex_t IT_0961 = s_12*s_14*s_34*IT_0018;
    const ccomplex_t IT_0962 = -IT_0525;
    const ccomplex_t IT_0963 = s_23*IT_0018*IT_0094*IT_0481;
    const ccomplex_t IT_0964 = s_13*s_14*s_24*IT_0018*IT_0094;
    const ccomplex_t IT_0965 = IT_0963 + IT_0964;
    const ccomplex_t IT_0966 = s_12*s_14*s_34*IT_0018*IT_0094;
    const ccomplex_t IT_0967 = s_23*IT_0011*IT_0094;
    const ccomplex_t IT_0968 = -IT_0966 + -IT_0967;
    const ccomplex_t IT_0969 = IT_0965 + IT_0968;
    const ccomplex_t IT_0970 = IT_0019*IT_0969;
    const ccomplex_t IT_0971 = 0.5*IT_0970;
    const ccomplex_t IT_0972 = IT_0529 + IT_0530 + IT_0961 + IT_0962 + IT_0971;
    const ccomplex_t IT_0973 = s_12*IT_0954;
    const ccomplex_t IT_0974 = s_13*s_23*IT_0094;
    const ccomplex_t IT_0975 = (-2)*IT_0974;
    const ccomplex_t IT_0976 = IT_0973 + IT_0975;
    const ccomplex_t IT_0977 = IT_0019*IT_0976;
    const ccomplex_t IT_0978 = 0.5*IT_0977;
    const ccomplex_t IT_0979 = s_14*s_23*s_34;
    const ccomplex_t IT_0980 = s_12*IT_0429;
    const ccomplex_t IT_0981 = (-2)*IT_0828 + 2*IT_0829 + 2*IT_0979 + (-2)
      *IT_0980;
    const ccomplex_t IT_0982 = IT_0827 + IT_0981;
    const ccomplex_t IT_0983 = IT_0018*IT_0982;
    const ccomplex_t IT_0984 = (-0.5)*IT_0983;
    const ccomplex_t IT_0985 = s_14*s_23*IT_0094;
    const ccomplex_t IT_0986 = s_13*s_24*IT_0094;
    const ccomplex_t IT_0987 = IT_0985 + IT_0986;
    const ccomplex_t IT_0988 = s_12*s_34*IT_0094;
    const ccomplex_t IT_0989 = -IT_0988;
    const ccomplex_t IT_0990 = IT_0987 + IT_0989;
    const ccomplex_t IT_0991 = IT_0252*IT_0990;
    const ccomplex_t IT_0992 = 0.5*IT_0991;
    const ccomplex_t IT_0993 = IT_0459 + IT_0461 + IT_0978 + IT_0984 + IT_0992;
    const ccomplex_t IT_0994 = 4*IT_0523;
    const ccomplex_t IT_0995 = s_12*s_24*s_34*IT_0018;
    const ccomplex_t IT_0996 = (-0.5)*IT_0995;
    const ccomplex_t IT_0997 = s_14*s_23*s_24*IT_0018*IT_0094;
    const ccomplex_t IT_0998 = s_13*IT_0018*IT_0070*IT_0094;
    const ccomplex_t IT_0999 = IT_0997 + IT_0998;
    const ccomplex_t IT_1000 = s_12*s_24*s_34*IT_0018*IT_0094;
    const ccomplex_t IT_1001 = s_13*IT_0012*IT_0094;
    const ccomplex_t IT_1002 = -IT_1000 + -IT_1001;
    const ccomplex_t IT_1003 = IT_0999 + IT_1002;
    const ccomplex_t IT_1004 = IT_0019*IT_1003;
    const ccomplex_t IT_1005 = (-0.25)*IT_1004;
    const ccomplex_t IT_1006 = (-0.5)*IT_0566;
    const ccomplex_t IT_1007 = 0.5*IT_0569;
    const ccomplex_t IT_1008 = (-0.5)*IT_0571;
    const ccomplex_t IT_1009 = IT_0565 + IT_0996 + IT_1005 + IT_1006 + IT_1007
       + IT_1008;
    const ccomplex_t IT_1010 = (-8)*IT_0017;
    const ccomplex_t IT_1011 = IT_0524 + IT_0962;
    const ccomplex_t IT_1012 = s_12*s_13*IT_0030;
    const ccomplex_t IT_1013 = s_23*IT_0481;
    const ccomplex_t IT_1014 = s_12*s_14*s_34;
    const ccomplex_t IT_1015 = s_23*IT_0011*IT_0030;
    const ccomplex_t IT_1016 = 2*IT_1013 + (-2)*IT_1014 + -IT_1015;
    const ccomplex_t IT_1017 = IT_1012 + IT_1016;
    const ccomplex_t IT_1018 = IT_0018*IT_1017;
    const ccomplex_t IT_1019 = (-0.5)*IT_1018;
    const ccomplex_t IT_1020 = IT_1011 + IT_1019;
    const ccomplex_t IT_1021 = (-4)*IT_0865;
    const ccomplex_t IT_1022 = (-4)*IT_0501;
    const ccomplex_t IT_1023 = s_12*IT_0030*IT_0954;
    const ccomplex_t IT_1024 = s_14*s_24*IT_0954;
    const ccomplex_t IT_1025 = (-2)*IT_1024;
    const ccomplex_t IT_1026 = IT_1023 + IT_1025;
    const ccomplex_t IT_1027 = IT_0356*IT_1026;
    const ccomplex_t IT_1028 = 0.5*IT_1027;
    const ccomplex_t IT_1029 = (-4)*IT_0460;
    const ccomplex_t IT_1030 = IT_0459 + IT_1028 + IT_1029;
    const ccomplex_t IT_1031 = (-2)*IT_0828 + (-2)*IT_0830 + 4*IT_0979;
    const ccomplex_t IT_1032 = IT_0827 + IT_1031;
    const ccomplex_t IT_1033 = IT_0018*IT_1032;
    const ccomplex_t IT_1034 = -IT_1033;
    const ccomplex_t IT_1035 = IT_1030 + IT_1034;
    const ccomplex_t IT_1036 = (-4)*IT_0300;
    const ccomplex_t IT_1037 = s_12*IT_0018*IT_0429;
    const ccomplex_t IT_1038 = -IT_1037;
    const ccomplex_t IT_1039 = s_14*s_23*s_34*IT_0018*IT_0094;
    const ccomplex_t IT_1040 = s_13*s_24*s_34*IT_0018*IT_0094;
    const ccomplex_t IT_1041 = IT_0973 + IT_1039 + IT_1040;
    const ccomplex_t IT_1042 = s_12*IT_0018*IT_0094*IT_0429;
    const ccomplex_t IT_1043 = (-2)*IT_0974 + -IT_1042;
    const ccomplex_t IT_1044 = IT_1041 + IT_1043;
    const ccomplex_t IT_1045 = IT_0019*IT_1044;
    const ccomplex_t IT_1046 = (-0.5)*IT_1045;
    const ccomplex_t IT_1047 = -IT_0535;
    const ccomplex_t IT_1048 = IT_0459 + IT_0925 + IT_1038 + IT_1046 + IT_1047;
    const ccomplex_t IT_1049 = (-4)*IT_0447;
    const ccomplex_t IT_1050 = IT_1048*IT_1049;
    const ccomplex_t IT_1051 = s_34*IT_0094;
    const ccomplex_t IT_1052 = IT_0252*IT_1051;
    const ccomplex_t IT_1053 = (-2)*IT_0429;
    const ccomplex_t IT_1054 = IT_0686 + IT_1053;
    const ccomplex_t IT_1055 = IT_0018*IT_1054;
    const ccomplex_t IT_1056 = (-0.333333333333333)*IT_1052 + (
      -0.333333333333333)*IT_1055;
    const ccomplex_t IT_1057 = IT_0094 + IT_1056;
    const ccomplex_t IT_1058 = IT_0428*IT_1057;
    const ccomplex_t IT_1059 = 3*IT_1058;
    const ccomplex_t IT_1060 = IT_0878*IT_1059;
    const ccomplex_t IT_1061 = IT_0765*IT_0960 + IT_0923*IT_0972 + IT_0993
      *IT_0994 + IT_1009*IT_1010 + IT_1020*(IT_1021 + IT_1022) + IT_1035*IT_1036
       + IT_1050 + IT_1060;
    const ccomplex_t IT_1062 = s_23*IT_0030*IT_0094;
    const ccomplex_t IT_1063 = IT_0356*IT_1062;
    const ccomplex_t IT_1064 = 0.666666666666667*IT_0251 + (-0.166666666666667
      )*IT_1063;
    const ccomplex_t IT_1065 = s_23 + IT_1064;
    const ccomplex_t IT_1066 = m_C_2*IT_1065;
    const ccomplex_t IT_1067 = (-6)*IT_1066;
    const ccomplex_t IT_1068 = s_13*s_14*s_34*IT_0019;
    const ccomplex_t IT_1069 = 2*IT_0481 + -IT_0890 + (-2)*IT_1068;
    const ccomplex_t IT_1070 = IT_0889 + IT_1069;
    const ccomplex_t IT_1071 = IT_0018*IT_1070;
    const ccomplex_t IT_1072 = -IT_0011 + (-0.5)*IT_1071;
    const ccomplex_t IT_1073 = IT_0495 + IT_1072;
    const ccomplex_t IT_1074 = m_N_3*IT_1073;
    const ccomplex_t IT_1075 = (-2)*IT_1074;
    const ccomplex_t IT_1076 = 2*IT_0865;
    const ccomplex_t IT_1077 = IT_0759 + IT_1076;
    const ccomplex_t IT_1078 = 2*conj(IT_0865);
    const ccomplex_t IT_1079 = IT_0634 + IT_1078;
    const ccomplex_t IT_1080 = IT_1075*(IT_0397*IT_0493 + IT_0391*conj(IT_0493
      ) + conj(IT_0419)*IT_0790 + IT_0419*IT_0857 + IT_0426*IT_0865 + IT_0420
      *conj(IT_0865) + IT_0652*IT_0907 + IT_0543*IT_0940 + IT_0644*IT_1077 +
       IT_0648*IT_1079);
    const ccomplex_t IT_1081 = -IT_0808;
    const ccomplex_t IT_1082 = IT_0807 + IT_1081;
    const ccomplex_t IT_1083 = IT_0018*IT_1082;
    const ccomplex_t IT_1084 = IT_0069 + IT_0840;
    const ccomplex_t IT_1085 = (-4)*IT_1084;
    const ccomplex_t IT_1086 = IT_1083 + IT_1085;
    const ccomplex_t IT_1087 = m_C_2*IT_1086;
    const ccomplex_t IT_1088 = 2*IT_0154;
    const ccomplex_t IT_1089 = IT_0607 + IT_1088;
    const ccomplex_t IT_1090 = 2*conj(IT_0154);
    const ccomplex_t IT_1091 = IT_0609 + IT_1090;
    const ccomplex_t IT_1092 = conj(IT_0221)*IT_0391 + IT_0221*IT_0397 + conj
      (IT_0154)*IT_0420 + IT_0154*IT_0426 + conj(IT_0390)*IT_0645 + IT_0390
      *IT_0649 + conj(IT_0419)*IT_0653 + IT_0419*IT_0656 + IT_0652*IT_0885 +
       IT_0543*IT_0936 + IT_0644*IT_1089 + IT_0648*IT_1091;
    const ccomplex_t IT_1093 = (-10)*s_23;
    const ccomplex_t IT_1094 = IT_1063 + IT_1093;
    const ccomplex_t IT_1095 = m_C_2*IT_1094;
    const ccomplex_t IT_1096 = conj(IT_0294) + conj(IT_0295);
    const ccomplex_t IT_1097 = IT_0294 + IT_0295;
    const ccomplex_t IT_1098 = s_23*m_C_2;
    const ccomplex_t IT_1099 = IT_0432*IT_1098;
    const ccomplex_t IT_1100 = -IT_0065;
    const ccomplex_t IT_1101 = IT_0028 + IT_1100;
    const ccomplex_t IT_1102 = s_34*m_C_2*m_N_3*IT_0018;
    const ccomplex_t IT_1103 = IT_0456*IT_1102;
    const ccomplex_t IT_1104 = IT_0451*IT_1102;
    const ccomplex_t IT_1105 = IT_0570 + IT_0572 + IT_0574 + IT_0619 + IT_0620
       + IT_0622 + IT_0995;
    const ccomplex_t IT_1106 = s_13*m_N_3;
    const ccomplex_t IT_1107 = IT_0432*IT_1106;
    const ccomplex_t IT_1108 = IT_0529 + IT_0531 + IT_0533 + IT_0781 + IT_0782
       + IT_0784 + IT_0961;
    const ccomplex_t IT_1109 = m_C_2*m_N_3*IT_0094;
    const ccomplex_t IT_1110 = IT_0432*IT_1109;
    const ccomplex_t IT_1111 = (-4)*IT_0523;
    const ccomplex_t IT_1112 = s_34*IT_0018;
    const ccomplex_t IT_1113 = IT_0796 + IT_0797;
    const ccomplex_t IT_1114 = s_13*s_23*s_34*IT_0019;
    const ccomplex_t IT_1115 = (-2)*IT_1114;
    const ccomplex_t IT_1116 = IT_1113 + IT_1115;
    const ccomplex_t IT_1117 = IT_1112*IT_1116;
    const ccomplex_t IT_1118 = (-2)*IT_0472;
    const ccomplex_t IT_1119 = (-0.25)*IT_0391*IT_0877 + (-0.25)*IT_0664
      *IT_0951 + IT_1099*IT_1101 + (-0.25)*IT_1077*IT_1103 + (-0.25)*IT_0885
      *IT_1104 + (-0.25)*IT_1089*IT_1105 + (-0.25)*IT_0445*IT_1107 + (-0.25)
      *IT_0907*IT_1108 + (-0.25)*IT_1110*IT_1111 + (-0.25)*IT_0778*IT_1117 + (
      -0.25)*IT_0915*IT_1118;
    const ccomplex_t IT_1120 = IT_0465 + IT_0471;
    const ccomplex_t IT_1121 = (-6)*conj(IT_0017);
    const ccomplex_t IT_1122 = -IT_0612;
    const ccomplex_t IT_1123 = IT_0094 + IT_1122;
    const ccomplex_t IT_1124 = IT_0428*IT_1123;
    const ccomplex_t IT_1125 = (-2)*conj(IT_0472);
    const ccomplex_t IT_1126 = IT_0067*IT_0452 + IT_0349*IT_0457 + IT_0609
      *IT_0625 + IT_0397*IT_0763 + IT_0665*IT_0777 + IT_0634*IT_1103 + IT_0660
      *IT_1104 + IT_0563*IT_1120 + IT_0720*IT_1121 + IT_0723*IT_1121 + IT_0726
      *IT_1124 + IT_0906*IT_1125;
    const ccomplex_t IT_1127 = m_N_3*IT_0018;
    const ccomplex_t IT_1128 = IT_0019*IT_0357;
    const ccomplex_t IT_1129 = (-2)*IT_0336 + (-0.5)*IT_1128;
    const ccomplex_t IT_1130 = IT_0335 + IT_1129;
    const ccomplex_t IT_1131 = IT_1127*IT_1130;
    const ccomplex_t IT_1132 = (-2)*IT_1131;
    const ccomplex_t IT_1133 = 2*IT_0472;
    const ccomplex_t IT_1134 = (-2)*IT_0692;
    const ccomplex_t IT_1135 = (-2)*IT_0696;
    const ccomplex_t IT_1136 = (-2)*IT_0700;
    const ccomplex_t IT_1137 = IT_1133 + IT_1134 + IT_1135 + IT_1136;
    const ccomplex_t IT_1138 = 2*conj(IT_0472);
    const ccomplex_t IT_1139 = (-2)*conj(IT_0692);
    const ccomplex_t IT_1140 = (-2)*conj(IT_0696);
    const ccomplex_t IT_1141 = (-2)*conj(IT_0700);
    const ccomplex_t IT_1142 = IT_1138 + IT_1139 + IT_1140 + IT_1141;
    const ccomplex_t IT_1143 = (-2)*IT_0707;
    const ccomplex_t IT_1144 = (-2)*conj(IT_0707);
    const ccomplex_t IT_1145 = s_24*s_34*IT_0094;
    const ccomplex_t IT_1146 = (-2)*IT_1145;
    const ccomplex_t IT_1147 = IT_1062 + IT_1146;
    const ccomplex_t IT_1148 = IT_0356*IT_1147;
    const ccomplex_t IT_1149 = (-0.666666666666667)*IT_0458 + (
      -0.166666666666667)*IT_1148;
    const ccomplex_t IT_1150 = s_23 + IT_1149;
    const ccomplex_t IT_1151 = m_C_2*IT_1150;
    const ccomplex_t IT_1152 = (-6)*IT_1151;
    const ccomplex_t IT_1153 = 4*IT_0160;
    const ccomplex_t IT_1154 = 4*IT_0576;
    const ccomplex_t IT_1155 = (-8)*IT_0160;
    const ccomplex_t IT_1156 = (-4)*IT_0764;
    const ccomplex_t IT_1157 = (-4)*IT_0493;
    const ccomplex_t IT_1158 = (-4)*IT_0293;
    const ccomplex_t IT_1159 = IT_0030*IT_0954;
    const ccomplex_t IT_1160 = IT_0356*IT_1159;
    const ccomplex_t IT_1161 = 0.125*IT_1160;
    const ccomplex_t IT_1162 = IT_0094 + IT_1161;
    const ccomplex_t IT_1163 = IT_0428*IT_1162;
    const ccomplex_t IT_1164 = 8*IT_1163;
    const ccomplex_t IT_1165 = (-2)*IT_0300;
    const ccomplex_t IT_1166 = m_C_2*m_N_3*IT_0018;
    const ccomplex_t IT_1167 = IT_0019*IT_1159;
    const ccomplex_t IT_1168 = (-0.25)*IT_1167;
    const ccomplex_t IT_1169 = IT_0429 + IT_1168;
    const ccomplex_t IT_1170 = IT_1166*IT_1169;
    const ccomplex_t IT_1171 = (-4)*IT_1170;
    const ccomplex_t IT_1172 = (-2)*IT_0318;
    const ccomplex_t IT_1173 = IT_0551*IT_0611 + IT_0594*IT_0611 + IT_0615
      *IT_0960 + IT_0824*IT_1059 + IT_0420*IT_1067 + IT_0972*IT_1153 + IT_0993
      *IT_1154 + IT_1009*IT_1155 + IT_1048*IT_1156 + IT_1020*IT_1157 + IT_1035
      *IT_1158 + IT_1164*IT_1165 + IT_1171*IT_1172;
    const ccomplex_t IT_1174 = (-2)*IT_0293;
    const ccomplex_t IT_1175 = IT_0018*IT_0686;
    const ccomplex_t IT_1176 = 12*IT_0094 + (-4)*IT_1175;
    const ccomplex_t IT_1177 = IT_1160 + IT_1176;
    const ccomplex_t IT_1178 = IT_0428*IT_1177;
    const ccomplex_t IT_1179 = (-2)*IT_0294;
    const ccomplex_t IT_1180 = (-2)*IT_0295;
    const ccomplex_t IT_1181 = (-4)*IT_0433;
    const ccomplex_t IT_1182 = -IT_0028;
    const ccomplex_t IT_1183 = -conj(IT_0028);
    const ccomplex_t IT_1184 = -IT_0222;
    const ccomplex_t IT_1185 = -conj(IT_0222);
    const ccomplex_t IT_1186 = conj(IT_0222)*IT_0419 + IT_0222*conj(IT_0419) +
       IT_0028*(conj(IT_0390) + conj(IT_0472)) + conj(IT_0244)*(IT_0406 +
       IT_0411 + -IT_0419 + IT_0707) + IT_0244*(conj(IT_0406) + conj(IT_0411) + 
      -conj(IT_0419) + conj(IT_0707)) + (IT_0390 + IT_0472)*(conj(IT_0028) +
       IT_0729) + -(conj(IT_0369) + conj(IT_0374) + conj(IT_0379) + -conj
      (IT_0390) + -conj(IT_0472) + conj(IT_0692) + conj(IT_0696) + conj(IT_0700)
      )*IT_1100 + (IT_0652 + conj(IT_0692) + conj(IT_0696) + conj(IT_0700))
      *IT_1182 + (IT_0543 + IT_0692 + IT_0696 + IT_0700)*(conj(IT_0065) +
       IT_1183) + (IT_0644 + conj(IT_0707))*IT_1184 + (IT_0648 + IT_0707)*IT_1185;
    const ccomplex_t IT_1187 = (-0.166666666666667)*IT_0361 + (
      -0.666666666666667)*IT_0469;
    const ccomplex_t IT_1188 = s_13 + IT_1187;
    const ccomplex_t IT_1189 = m_N_3*IT_1188;
    const ccomplex_t IT_1190 = (-6)*IT_1189;
    const ccomplex_t IT_1191 = m_C_2*IT_0018;
    const ccomplex_t IT_1192 = IT_0019*IT_1062;
    const ccomplex_t IT_1193 = (-2)*IT_0256 + (-0.5)*IT_1192;
    const ccomplex_t IT_1194 = IT_0255 + IT_1193;
    const ccomplex_t IT_1195 = IT_1191*IT_1194;
    const ccomplex_t IT_1196 = (-2)*IT_1195;
    const ccomplex_t IT_1197 = IT_0356*IT_0357;
    const ccomplex_t IT_1198 = (-10)*s_13;
    const ccomplex_t IT_1199 = IT_1197 + IT_1198;
    const ccomplex_t IT_1200 = m_N_3*IT_1199;
    const ccomplex_t IT_1201 = s_13*s_14*s_24;
    const ccomplex_t IT_1202 = s_24*s_34*IT_0011;
    const ccomplex_t IT_1203 = -IT_1015 + (-2)*IT_1201 + 2*IT_1202;
    const ccomplex_t IT_1204 = IT_1012 + IT_1203;
    const ccomplex_t IT_1205 = IT_0018*IT_1204;
    const ccomplex_t IT_1206 = (-0.5)*IT_1205;
    const ccomplex_t IT_1207 = IT_1011 + IT_1206;
    const ccomplex_t IT_1208 = 4*conj(IT_0493);
    const ccomplex_t IT_1209 = (-0.5)*IT_0961;
    const ccomplex_t IT_1210 = (-0.25)*IT_0970;
    const ccomplex_t IT_1211 = (-0.5)*IT_0525;
    const ccomplex_t IT_1212 = (-0.5)*IT_0528;
    const ccomplex_t IT_1213 = 0.5*IT_0530;
    const ccomplex_t IT_1214 = IT_0524 + IT_1209 + IT_1210 + IT_1211 + IT_1212
       + IT_1213;
    const ccomplex_t IT_1215 = (-8)*conj(IT_0160);
    const ccomplex_t IT_1216 = 0.25*IT_1027;
    const ccomplex_t IT_1217 = IT_0459 + IT_1216;
    const ccomplex_t IT_1218 = (-8)*conj(IT_0293);
    const ccomplex_t IT_1219 = (-8)*conj(IT_0294);
    const ccomplex_t IT_1220 = (-8)*conj(IT_0295);
    const ccomplex_t IT_1221 = -IT_0925;
    const ccomplex_t IT_1222 = IT_0459 + IT_0535 + IT_1038 + IT_1046 + IT_1221;
    const ccomplex_t IT_1223 = (-4)*conj(IT_0764);
    const ccomplex_t IT_1224 = -IT_0995;
    const ccomplex_t IT_1225 = (-0.5)*IT_1004;
    const ccomplex_t IT_1226 = IT_0566 + IT_0570 + IT_0571 + IT_1224 + IT_1225;
    const ccomplex_t IT_1227 = s_13*IT_0070;
    const ccomplex_t IT_1228 = s_12*s_24*s_34;
    const ccomplex_t IT_1229 = -IT_0671 + 2*IT_1227 + (-2)*IT_1228;
    const ccomplex_t IT_1230 = IT_0669 + IT_1229;
    const ccomplex_t IT_1231 = IT_0018*IT_1230;
    const ccomplex_t IT_1232 = (-0.5)*IT_1231;
    const ccomplex_t IT_1233 = IT_0668 + IT_1232;
    const ccomplex_t IT_1234 = (-4)*conj(IT_0111);
    const ccomplex_t IT_1235 = (-4)*conj(IT_0154);
    const ccomplex_t IT_1236 = 0.666666666666667*IT_0332 + (-0.166666666666667
      )*IT_1197;
    const ccomplex_t IT_1237 = s_13 + IT_1236;
    const ccomplex_t IT_1238 = m_N_3*IT_1237;
    const ccomplex_t IT_1239 = (-6)*IT_1238;
    const ccomplex_t IT_1240 = (-0.666666666666667)*IT_0259 +
       0.166666666666667*IT_1148;
    const ccomplex_t IT_1241 = s_23 + IT_1240;
    const ccomplex_t IT_1242 = m_C_2*IT_1241;
    const ccomplex_t IT_1243 = 6*IT_1242;
    const ccomplex_t IT_1244 = (-2)*conj(IT_0318);
    const ccomplex_t IT_1245 = (-2)*conj(IT_0329);
    const ccomplex_t IT_1246 = (-2)*conj(IT_0300);
    const ccomplex_t IT_1247 = (-2)*conj(IT_0330);
    const ccomplex_t IT_1248 = 4*conj(IT_0576);
    const ccomplex_t IT_1249 = IT_0993*IT_1248;
    const ccomplex_t IT_1250 = IT_0856*IT_1059;
    const ccomplex_t IT_1251 = IT_0626*IT_0960;
    const ccomplex_t IT_1252 = IT_0466*IT_0551 + IT_0466*IT_0594 + IT_0398
      *IT_1190 + IT_0427*IT_1196 + IT_1138*IT_1200 + IT_1207*IT_1208 + IT_1214
      *IT_1215 + IT_1217*(IT_1218 + IT_1219 + IT_1220) + IT_1222*IT_1223 +
       IT_0746*IT_1226 + IT_1233*(IT_1234 + IT_1235) + (IT_1139 + IT_1140 +
       IT_1141)*IT_1239 + IT_1144*IT_1243 + IT_1164*(IT_1244 + IT_1245) +
       IT_1171*(IT_1246 + IT_1247) + IT_1249 + IT_1250 + IT_1251;
    const ccomplex_t IT_1253 = (-2)*IT_0828 + 2*IT_0829 + (-2)*IT_0830 + 2
      *IT_0979 + (-2)*IT_0980;
    const ccomplex_t IT_1254 = IT_0827 + IT_1253;
    const ccomplex_t IT_1255 = IT_0356*IT_1254;
    const ccomplex_t IT_1256 = (-0.25)*IT_1255;
    const ccomplex_t IT_1257 = -IT_0022;
    const ccomplex_t IT_1258 = -IT_0021;
    const ccomplex_t IT_1259 = 8*IT_0369;
    const ccomplex_t IT_1260 = (-8)*IT_0390;
    const ccomplex_t IT_1261 = IT_1259 + IT_1260;
    const ccomplex_t IT_1262 = 8*conj(IT_0369);
    const ccomplex_t IT_1263 = (-8)*conj(IT_0390);
    const ccomplex_t IT_1264 = IT_1262 + IT_1263;
    const ccomplex_t IT_1265 = 8*IT_0390;
    const ccomplex_t IT_1266 = 8*IT_0374;
    const ccomplex_t IT_1267 = IT_1259 + IT_1266;
    const ccomplex_t IT_1268 = 8*conj(IT_0374);
    const ccomplex_t IT_1269 = 8*IT_0406;
    const ccomplex_t IT_1270 = (-8)*IT_0419;
    const ccomplex_t IT_1271 = 8*conj(IT_0406);
    const ccomplex_t IT_1272 = (-8)*conj(IT_0419);
    const ccomplex_t IT_1273 = 8*conj(IT_0692);
    const ccomplex_t IT_1274 = (-8)*conj(IT_0472);
    const ccomplex_t IT_1275 = IT_1273 + IT_1274;
    const ccomplex_t IT_1276 = 0.125*IT_0369;
    const ccomplex_t IT_1277 = 8*IT_0692;
    const ccomplex_t IT_1278 = (-8)*IT_0472;
    const ccomplex_t IT_1279 = IT_1277 + IT_1278;
    const ccomplex_t IT_1280 = 0.125*conj(IT_0369);
    const ccomplex_t IT_1281 = 8*conj(IT_0696);
    const ccomplex_t IT_1282 = IT_1273 + IT_1281;
    const ccomplex_t IT_1283 = IT_1274 + IT_1282;
    const ccomplex_t IT_1284 = 0.125*IT_0374;
    const ccomplex_t IT_1285 = 8*IT_0696;
    const ccomplex_t IT_1286 = IT_1277 + IT_1285;
    const ccomplex_t IT_1287 = IT_1278 + IT_1286;
    const ccomplex_t IT_1288 = 0.125*conj(IT_0374);
    const ccomplex_t IT_1289 = 8*conj(IT_0700);
    const ccomplex_t IT_1290 = 0.125*IT_0379;
    const ccomplex_t IT_1291 = 8*IT_0700;
    const ccomplex_t IT_1292 = 0.125*conj(IT_0379);
    const ccomplex_t IT_1293 = conj(IT_0390)*IT_0472 + 0.125*conj(IT_0411)
      *IT_0924 + 0.125*conj(IT_0692)*IT_1260 + 0.125*conj(IT_0696)*IT_1261 +
       0.125*IT_0692*IT_1263 + 0.125*IT_0696*IT_1264 + 0.125*conj(IT_0472)
      *IT_1265 + 0.125*conj(IT_0700)*(IT_1260 + IT_1267) + 0.125*IT_0700*
      (IT_1262 + IT_1263 + IT_1268) + conj(IT_0707)*(IT_0411 + 0.125*IT_1269 +
       0.125*IT_1270) + 0.125*IT_0707*(IT_1271 + IT_1272) + IT_1275*IT_1276 +
       IT_1279*IT_1280 + IT_1283*IT_1284 + IT_1287*IT_1288 + (IT_1274 + IT_1282 
      + IT_1289)*IT_1290 + (IT_1278 + IT_1286 + IT_1291)*IT_1292;
    const ccomplex_t IT_1294 = IT_0796 + IT_0798;
    const ccomplex_t IT_1295 = -IT_0797;
    const ccomplex_t IT_1296 = IT_1294 + IT_1295;
    const ccomplex_t IT_1297 = IT_0252*IT_1296;
    const ccomplex_t IT_1298 = 0.25*IT_1297;
    const ccomplex_t IT_1299 = IT_0795 + IT_1298;
    const ccomplex_t IT_1300 = 8*IT_0244;
    const ccomplex_t IT_1301 = (-8)*IT_0222;
    const ccomplex_t IT_1302 = IT_1300 + IT_1301;
    const ccomplex_t IT_1303 = 8*conj(IT_0244);
    const ccomplex_t IT_1304 = (-8)*conj(IT_0222);
    const ccomplex_t IT_1305 = IT_1303 + IT_1304;
    const ccomplex_t IT_1306 = 0.125*conj(IT_0390);
    const ccomplex_t IT_1307 = 0.125*IT_0648;
    const ccomplex_t IT_1308 = 0.125*IT_0644;
    const ccomplex_t IT_1309 = 8*IT_1299*(IT_0419*IT_0730 + 0.125*IT_0390
      *IT_0855 + conj(IT_0419)*IT_1101 + 0.125*IT_0652*IT_1302 + 0.125*IT_0543
      *IT_1305 + IT_0806*IT_1306 + IT_0939*IT_1307 + IT_0888*IT_1308);
    const ccomplex_t IT_1310 = IT_0356*IT_1032;
    const ccomplex_t IT_1311 = 0.25*IT_1310;
    const ccomplex_t IT_1312 = 8*IT_0379;
    const ccomplex_t IT_1313 = IT_1260 + IT_1267 + IT_1312;
    const ccomplex_t IT_1314 = IT_0419*conj(IT_0419) + 0.125*IT_0406*IT_1272 +
       IT_0411*(conj(IT_0411) + 0.125*IT_1271 + 0.125*IT_1272) + IT_1263*IT_1276
       + IT_1261*IT_1280 + IT_1264*IT_1284 + (IT_1260 + IT_1267)*IT_1288 + 
      (IT_1262 + IT_1263 + IT_1268)*IT_1290 + IT_1265*IT_1306 + (IT_1269 +
       IT_1270)*IT_1308 + IT_1292*IT_1313;
    const ccomplex_t IT_1315 = IT_0356*IT_1054;
    const ccomplex_t IT_1316 = 0.25*IT_1315;
    const ccomplex_t IT_1317 = 1 + IT_1316;
    const ccomplex_t IT_1318 = IT_0428*IT_1317;
    const ccomplex_t IT_1319 = IT_0419*conj(IT_0472) + 0.125*conj(IT_0419)
      *IT_0836 + IT_0707*(conj(IT_0379) + 0.125*IT_1262 + 0.125*IT_1263 + 0.125
      *IT_1268) + 0.125*(conj(IT_0692) + conj(IT_0696) + conj(IT_0700))*IT_1270 
      + 0.125*(IT_0692 + IT_0696 + IT_0700)*IT_1272 + (IT_1274 + IT_1282 +
       IT_1289)*IT_1307 + (IT_1278 + IT_1286 + IT_1291)*IT_1308 + 0.125*conj
      (IT_0707)*IT_1313;
    const ccomplex_t IT_1320 = IT_0428*IT_0489;
    const ccomplex_t IT_1321 = s_12 + IT_1257;
    const ccomplex_t IT_1322 = IT_0428*IT_1321;
    const ccomplex_t IT_1323 = IT_0428*IT_0841;
    const ccomplex_t IT_1324 = 0.166666666666667*IT_0754;
    const ccomplex_t IT_1325 = 0.166666666666667*IT_0025;
    const ccomplex_t IT_1326 = 0.166666666666667*IT_0542;
    const ccomplex_t IT_1327 = IT_0428*IT_0617;
    const ccomplex_t IT_1328 = 0.166666666666667*IT_1327;
    const ccomplex_t IT_1329 = 0.166666666666667*IT_0492;
    const ccomplex_t IT_1330 = (-6)*IT_0017;
    const ccomplex_t IT_1331 = 0.166666666666667*IT_1322;
    const ccomplex_t IT_1332 = 0.166666666666667*IT_0500;
    const ccomplex_t IT_1333 = IT_0075*IT_0428;
    const ccomplex_t IT_1334 = 0.166666666666667*IT_1333;
    const ccomplex_t IT_1335 = 0.166666666666667*IT_0751;
    const ccomplex_t IT_1336 = IT_0024*IT_0428;
    const ccomplex_t IT_1337 = 0.166666666666667*IT_1336;
    const ccomplex_t IT_1338 = 0.166666666666667*IT_0932;
    const ccomplex_t IT_1339 = 0.166666666666667*IT_0076*IT_0158 +
       0.666666666666667*IT_0419*IT_0742 + 0.166666666666667*IT_0352*IT_0753 +
       0.166666666666667*(IT_0748 + IT_0756)*IT_0757 + (-0.666666666666667)
      *IT_0472*IT_0922 + 0.166666666666667*IT_0468*IT_1153 + IT_0017*IT_1320 + 
      -IT_0017*(IT_1322 + -IT_1323) + IT_0346*IT_1324 + IT_0246*IT_1325 +
       IT_0420*IT_1326 + IT_0645*IT_1328 + IT_0661*IT_1329 + IT_1330*IT_1331 +
       IT_0759*IT_1332 + IT_0659*IT_1334 + IT_0391*IT_1335 + IT_0790*IT_1337 +
       IT_1133*IT_1338;
    const ccomplex_t IT_1340 = 0.666666666666667*conj(IT_0419)*IT_0742 +
       0.166666666666667*IT_0355*IT_0753 + (-0.666666666666667)*conj(IT_0472)
      *IT_0922 + conj(IT_0017)*IT_1320 + conj(IT_0017)*IT_1323 + IT_0349*IT_1324
       + IT_0249*IT_1325 + IT_0426*IT_1326 + IT_0649*IT_1328 + IT_0663*IT_1329 +
       2*IT_1121*IT_1331 + IT_0634*IT_1332 + IT_0660*IT_1334 + IT_0397*IT_1335 +
       IT_0857*IT_1337 + IT_1138*IT_1338;
    const ccomplex_t IT_1341 = (1 + 0.5*IT_0430)*(IT_0459 + (-2)*IT_0460);
    const ccomplex_t IT_1342 = 0.25*IT_1117;
    const ccomplex_t IT_1343 = 0.25*IT_1108;
    const ccomplex_t IT_1344 = 0.25*IT_1105;
    const ccomplex_t IT_1345 = (-0.25)*IT_0397*IT_0877 + (-0.25)*IT_0665
      *IT_0951 + -IT_0733*IT_1099 + (-0.25)*IT_1079*IT_1103 + (-0.25)*IT_0936
      *IT_1104 + -IT_0436*IT_1107 + conj(IT_0523)*IT_1110 + (-0.25)*IT_0915
      *IT_1125 + conj(IT_0576)*IT_1341 + -IT_0779*IT_1342 + -IT_0940*IT_1343 + 
      -IT_1091*IT_1344;
    const ccomplex_t IT_1346 = 0.25*IT_0575;
    const ccomplex_t IT_1347 = (-2)*conj(IT_0294);
    const ccomplex_t IT_1348 = (-2)*conj(IT_0295);
    const ccomplex_t IT_1349 = (-0.25)*IT_0960;
    const ccomplex_t IT_1350 = (-0.25)*IT_1059;
    const ccomplex_t IT_1351 = IT_0879*IT_1350;
    const ccomplex_t IT_1352 = (-2)*conj(IT_0293);
    const ccomplex_t IT_1353 = (-0.25)*IT_1164;
    const ccomplex_t IT_1354 = 4*conj(IT_0523);
    const ccomplex_t IT_1355 = (-0.25)*IT_0993;
    const ccomplex_t IT_1356 = (-4)*conj(IT_0447);
    const ccomplex_t IT_1357 = (-0.25)*IT_1048;
    const ccomplex_t IT_1358 = (-4)*conj(IT_0865);
    const ccomplex_t IT_1359 = (-4)*conj(IT_0501);
    const ccomplex_t IT_1360 = IT_1358 + IT_1359;
    const ccomplex_t IT_1361 = (-0.25)*IT_1020;
    const ccomplex_t IT_1362 = (-0.25)*IT_0972;
    const ccomplex_t IT_1363 = (-8)*conj(IT_0017);
    const ccomplex_t IT_1364 = (-0.25)*IT_1009;
    const ccomplex_t IT_1365 = IT_0399*IT_1035 + (-0.25)*IT_1178*(IT_1347 +
       IT_1348) + IT_0766*IT_1349 + IT_1351 + IT_1352*IT_1353 + IT_1354*IT_1355 
      + IT_1356*IT_1357 + IT_1360*IT_1361 + IT_0941*IT_1362 + IT_1363*IT_1364;
    const ccomplex_t IT_1366 = (-4)*IT_0294;
    const ccomplex_t IT_1367 = IT_0551*IT_0611;
    const ccomplex_t IT_1368 = (-0.25)*IT_1367;
    const ccomplex_t IT_1369 = IT_0594*IT_0611;
    const ccomplex_t IT_1370 = (-0.25)*IT_1369;
    const ccomplex_t IT_1371 = IT_0615*IT_1349;
    const ccomplex_t IT_1372 = (-0.25)*IT_1171;
    const ccomplex_t IT_1373 = IT_1154*IT_1355;
    const ccomplex_t IT_1374 = IT_0824*IT_1350 + IT_1165*IT_1353 + IT_1156
      *IT_1357 + IT_1157*IT_1361 + IT_1153*IT_1362 + IT_1155*IT_1364 + IT_1035*
      (IT_0295 + (-0.25)*IT_1158 + (-0.25)*IT_1366) + IT_1368 + IT_1370 +
       IT_1371 + IT_1172*IT_1372 + IT_1373;
    const ccomplex_t IT_1375 = IT_0018*IT_0832;
    const ccomplex_t IT_1376 = -IT_1375;
    const ccomplex_t IT_1377 = IT_1030 + IT_1376;
    const ccomplex_t IT_1378 = (-8)*conj(IT_0300);
    const ccomplex_t IT_1379 = (-0.125)*IT_1378;
    const ccomplex_t IT_1380 = conj(IT_0330) + IT_1379;
    const ccomplex_t IT_1381 = (-0.25)*IT_1196;
    const ccomplex_t IT_1382 = (-0.25)*IT_1190;
    const ccomplex_t IT_1383 = (-0.25)*IT_1243;
    const ccomplex_t IT_1384 = (-0.25)*IT_1239;
    const ccomplex_t IT_1385 = (-0.25)*IT_1222;
    const ccomplex_t IT_1386 = (-0.25)*IT_1214;
    const ccomplex_t IT_1387 = 4*conj(IT_0865);
    const ccomplex_t IT_1388 = 4*conj(IT_0501);
    const ccomplex_t IT_1389 = IT_1387 + IT_1388;
    const ccomplex_t IT_1390 = (-0.25)*IT_1207;
    const ccomplex_t IT_1391 = (-4)*conj(IT_0183);
    const ccomplex_t IT_1392 = (-4)*conj(IT_0221);
    const ccomplex_t IT_1393 = IT_1391 + IT_1392;
    const ccomplex_t IT_1394 = (-0.25)*IT_1233;
    const ccomplex_t IT_1395 = (-0.25)*IT_1226;
    const ccomplex_t IT_1396 = IT_0766*IT_1349 + IT_1351 + IT_1354*IT_1355 + 
      (IT_1347 + IT_1348 + IT_1352)*IT_1372 + (conj(IT_0318) + conj(IT_0329))
      *IT_1377 + 2*IT_1217*IT_1380 + IT_0398*IT_1381 + IT_0427*IT_1382 + IT_1142
      *IT_1383 + IT_1144*IT_1384 + IT_1356*IT_1385 + IT_1363*IT_1386 + IT_1389
      *IT_1390 + IT_1393*IT_1394 + IT_0744*IT_1395;
    const ccomplex_t IT_1397 = (-8)*IT_0293;
    const ccomplex_t IT_1398 = (-8)*IT_0294;
    const ccomplex_t IT_1399 = (-8)*IT_0295;
    const ccomplex_t IT_1400 = (-4)*IT_0281;
    const ccomplex_t IT_1401 = (-2)*IT_0329;
    const ccomplex_t IT_1402 = (-2)*IT_0330;
    const ccomplex_t IT_1403 = 4*IT_0493;
    const ccomplex_t IT_1404 = (-4)*IT_0111;
    const ccomplex_t IT_1405 = (-4)*IT_0154;
    const ccomplex_t IT_1406 = IT_1133*IT_1200 + (-4)*IT_0824*IT_1350 + (-4)
      *IT_1368 + (-4)*IT_1370 + (-4)*IT_1371 + (-4)*IT_1373 + (-4)*IT_0421
      *IT_1381 + (-4)*IT_0392*IT_1382 + (-4)*IT_1143*IT_1383 + (-4)*(IT_1134 +
       IT_1135 + IT_1136)*IT_1384 + (-4)*IT_1156*IT_1385 + (-4)*IT_1155*IT_1386 
      + (-4)*IT_0745*IT_1395 + IT_1217*(IT_1397 + IT_1398 + IT_1399) + (-4)
      *IT_1377*(IT_0292 + (-0.25)*IT_1400) + (-4)*IT_1353*(IT_1172 + IT_1401) + 
      (-4)*IT_1372*(IT_1165 + IT_1402) + (-4)*IT_1390*IT_1403 + (-4)*IT_1394*
      (IT_1404 + IT_1405);
    const ccomplex_t IT_1407 = s_12 + IT_1258;
    const ccomplex_t IT_1408 = IT_0428*IT_1407;
    const ccomplex_t IT_1409 = IT_0428*IT_0894;
    const ccomplex_t IT_1410 = (-0.166666666666667)*IT_0780;
    const ccomplex_t IT_1411 = (-0.166666666666667)*IT_0754;
    const ccomplex_t IT_1412 = (-0.166666666666667)*IT_0500;
    const ccomplex_t IT_1413 = (-0.166666666666667)*IT_1336;
    const ccomplex_t IT_1414 = IT_0615 + IT_0824;
    const ccomplex_t IT_1415 = (-0.166666666666667)*IT_1103;
    const ccomplex_t IT_1416 = (-0.166666666666667)*IT_0789;
    const ccomplex_t IT_1417 = (-0.166666666666667)*IT_0633;
    const ccomplex_t IT_1418 = (-0.166666666666667)*IT_0897;
    const ccomplex_t IT_1419 = (-0.166666666666667)*IT_0787;
    const ccomplex_t IT_1420 = IT_0183*IT_1408 + -IT_0501*IT_1409 + IT_0352
      *IT_1410 + IT_0346*IT_1411 + IT_0611*IT_1412 + IT_0749*IT_1413 + IT_1414
      *IT_1415 + IT_0870*IT_1416 + IT_0607*IT_1417 + IT_1133*IT_1418 + IT_0765
      *IT_1419;
    const ccomplex_t IT_1421 = IT_0626 + IT_0856;
    const ccomplex_t IT_1422 = (-6)*IT_0493*(conj(IT_0183)*IT_1408 + -conj
      (IT_0501)*IT_1409 + IT_0355*IT_1410 + IT_0349*IT_1411 + IT_0752*IT_1413 +
       IT_0609*IT_1417 + IT_1138*IT_1418 + IT_0766*IT_1419 + IT_1415*IT_1421);
    const ccomplex_t IT_1423 = (-0.166666666666667)*IT_0852;
    const ccomplex_t IT_1424 = (-0.166666666666667)*IT_0511;
    const ccomplex_t IT_1425 = (-0.166666666666667)*IT_0560;
    const ccomplex_t IT_1426 = (-0.166666666666667)*IT_1108;
    const ccomplex_t IT_1427 = -IT_0662*IT_1408 + IT_0493*IT_1409 + -IT_0346
      *IT_1410 + -IT_0352*IT_1411 + -IT_0561*IT_1412 + -IT_0757*IT_1413 + 
      -IT_0765*IT_1415 + -IT_0867*IT_1416 + -IT_0659*IT_1417 + -IT_0615*IT_1419 
      + -IT_1118*IT_1423 + -IT_0611*IT_1424 + -IT_0611*IT_1425 + -IT_0824*IT_1426;
    const ccomplex_t IT_1428 = -IT_0155*IT_1408 + conj(IT_0493)*IT_1409 + 
      -IT_0349*IT_1410 + -IT_0355*IT_1411 + -IT_0563*IT_1412 + -IT_0758*IT_1413 
      + -IT_0766*IT_1415 + -IT_0857*IT_1416 + -IT_0660*IT_1417 + -IT_0626
      *IT_1419 + -IT_1125*IT_1423 + -IT_0466*IT_1424 + -IT_0466*IT_1425 + 
      -IT_0856*IT_1426;
    const ccomplex_t IT_1429 = (-0.166666666666667)*IT_1104;
    const ccomplex_t IT_1430 = (-0.166666666666667)*IT_0625;
    const ccomplex_t IT_1431 = (-0.166666666666667)*IT_0158*IT_0452 + (
      -0.166666666666667)*IT_0346*IT_0457 + IT_0017*IT_0723 + (
      -0.166666666666667)*IT_0391*IT_0763 + (-0.166666666666667)*IT_0614*IT_0765
       + (-0.166666666666667)*IT_0664*IT_0777 + (-0.166666666666667)*IT_0906
      *IT_1118 + (-0.166666666666667)*IT_0561*IT_1120 + (-0.166666666666667)
      *IT_0724*IT_1124 + (-0.166666666666667)*IT_0720*IT_1330 + IT_0759*IT_1415 
      + IT_0659*IT_1429 + IT_0607*IT_1430;
    const ccomplex_t IT_1432 = IT_0428*IT_1084;
    const ccomplex_t IT_1433 = (-0.166666666666667)*IT_0618;
    const ccomplex_t IT_1434 = (-0.166666666666667)*IT_1327;
    const ccomplex_t IT_1435 = (-0.166666666666667)*IT_1333;
    const ccomplex_t IT_1436 = (-0.166666666666667)*IT_0761;
    const ccomplex_t IT_1437 = (-0.166666666666667)*IT_0816;
    const ccomplex_t IT_1438 = (-0.166666666666667)*IT_0521;
    const ccomplex_t IT_1439 = (-0.166666666666667)*IT_0604;
    const ccomplex_t IT_1440 = (-0.166666666666667)*IT_1105;
    const ccomplex_t IT_1441 = IT_0493*IT_1408 + IT_1077*IT_1417 + IT_0765
      *IT_1429 + IT_0615*IT_1430 + -IT_0662*IT_1432 + IT_0158*IT_1433 + IT_0757
      *IT_1434 + IT_0561*IT_1435 + IT_0647*IT_1436 + IT_1133*IT_1437 + IT_0611
      *IT_1438 + IT_0611*IT_1439 + IT_0824*IT_1440;
    const ccomplex_t IT_1442 = conj(IT_0493)*IT_1408 + IT_1079*IT_1417 +
       IT_0766*IT_1429 + IT_0626*IT_1430 + -IT_0155*IT_1432 + IT_0067*IT_1433 +
       IT_0758*IT_1434 + IT_0563*IT_1435 + IT_0649*IT_1436 + IT_1138*IT_1437 +
       IT_0466*IT_1438 + IT_0466*IT_1439 + IT_0856*IT_1440;
    const ccomplex_t IT_1443 = IT_0355*IT_0754;
    const ccomplex_t IT_1444 = IT_0611*IT_1333;
    const ccomplex_t IT_1445 = IT_0749*IT_1327;
    const ccomplex_t IT_1446 = IT_1104*IT_1414;
    const ccomplex_t IT_1447 = (-6)*IT_0501;
    const ccomplex_t IT_1448 = IT_1408*IT_1447;
    const ccomplex_t IT_1449 = 6*IT_0183;
    const ccomplex_t IT_1450 = IT_1432*IT_1449;
    const ccomplex_t IT_1451 = IT_0883*IT_1118;
    const ccomplex_t IT_1452 = IT_0466*IT_1333;
    const ccomplex_t IT_1453 = IT_0249*IT_0618;
    const ccomplex_t IT_1454 = IT_0752*IT_1327;
    const ccomplex_t IT_1455 = IT_1104*IT_1421;
    const ccomplex_t IT_1456 = (-6)*conj(IT_0501);
    const ccomplex_t IT_1457 = IT_1408*IT_1456;
    const ccomplex_t IT_1458 = 6*conj(IT_0183);
    const ccomplex_t IT_1459 = IT_1432*IT_1458;
    const ccomplex_t IT_1460 = IT_0883*IT_1125;
    const ccomplex_t IT_1461 = (-4)*conj(IT_0493);
    const ccomplex_t IT_1462 = IT_1020*IT_1461;
    const ccomplex_t IT_1463 = IT_1048*IT_1223;
    const ccomplex_t IT_1464 = IT_1009*IT_1215;
    const ccomplex_t IT_1465 = 4*conj(IT_0160);
    const ccomplex_t IT_1466 = IT_0972*IT_1465;
    const ccomplex_t IT_1467 = IT_1249 + IT_1250 + IT_1251 + IT_1462 + IT_1463
       + IT_1464 + IT_1466;
    const ccomplex_t IT_1468 = (-4)*conj(IT_0293);
    const ccomplex_t IT_1469 = IT_1164*IT_1246;
    const ccomplex_t IT_1470 = IT_1171*IT_1244;
    const ccomplex_t IT_1471 = IT_0480*IT_0749;
    const ccomplex_t IT_1472 = IT_0922*IT_1153;
    const ccomplex_t IT_1473 = IT_0611*IT_0844;
    const ccomplex_t IT_1474 = IT_0757*IT_0932;
    const ccomplex_t IT_1475 = IT_0767*IT_0818;
    const ccomplex_t IT_1476 = IT_0765*IT_0906;
    const ccomplex_t IT_1477 = IT_0852*IT_1077;
    const ccomplex_t IT_1478 = IT_0870*IT_0897;
    const ccomplex_t IT_1479 = IT_0647*IT_0816;
    const ccomplex_t IT_1480 = IT_0883*IT_1089;
    const ccomplex_t IT_1481 = IT_0823*IT_0880;
    const ccomplex_t IT_1482 = IT_0878*IT_0915;
    const ccomplex_t IT_1483 = IT_0803*IT_1302;
    const ccomplex_t IT_1484 = IT_0480*IT_0752;
    const ccomplex_t IT_1485 = IT_0922*IT_1465;
    const ccomplex_t IT_1486 = IT_0466*IT_0844;
    const ccomplex_t IT_1487 = IT_0758*IT_0932;
    const ccomplex_t IT_1488 = IT_0768*IT_0818;
    const ccomplex_t IT_1489 = IT_0766*IT_0906;
    const ccomplex_t IT_1490 = IT_0852*IT_1079;
    const ccomplex_t IT_1491 = IT_0871*IT_0897;
    const ccomplex_t IT_1492 = IT_0651*IT_0816;
    const ccomplex_t IT_1493 = IT_0883*IT_1091;
    const ccomplex_t IT_1494 = IT_0823*IT_0881;
    const ccomplex_t IT_1495 = IT_0879*IT_0915;
    const ccomplex_t IT_1496 = IT_0803*IT_1305;
    const ccomplex_t IT_1497 = IT_0766*IT_0960;
    const ccomplex_t IT_1498 = IT_0879*IT_1059;
    const ccomplex_t IT_1499 = IT_0993*IT_1354;
    const ccomplex_t IT_1500 = (-0.125)*IT_0427*IT_1190 + (-0.125)*IT_0398
      *IT_1196 + (-0.125)*IT_1144*IT_1200 + (-0.125)*IT_0744*IT_1226 + (-0.125)
      *IT_1142*IT_1243 + (-0.125)*IT_1222*IT_1356 + (-0.125)*IT_1214*IT_1363 +
       IT_1217*IT_1380 + (-0.125)*IT_1207*IT_1389 + (-0.125)*IT_1233*IT_1393 + (
      -0.125)*IT_1497 + (-0.125)*IT_1498 + (-0.125)*IT_1499;
    const ccomplex_t IT_1501 = (-0.125)*IT_1060;
    const ccomplex_t IT_1502 = IT_0392*IT_1196;
    const ccomplex_t IT_1503 = (-0.125)*IT_1502;
    const ccomplex_t IT_1504 = IT_0421*IT_1190;
    const ccomplex_t IT_1505 = (-0.125)*IT_1504;
    const ccomplex_t IT_1506 = IT_1137*IT_1243;
    const ccomplex_t IT_1507 = (-0.125)*IT_1506;
    const ccomplex_t IT_1508 = (-4)*IT_0318;
    const ccomplex_t IT_1509 = IT_1377*IT_1508;
    const ccomplex_t IT_1510 = (-0.125)*IT_1509;
    const ccomplex_t IT_1511 = IT_0993*IT_0994;
    const ccomplex_t IT_1512 = (-0.125)*IT_1511;
    const ccomplex_t IT_1513 = IT_1049*IT_1222;
    const ccomplex_t IT_1514 = (-0.125)*IT_1513;
    const ccomplex_t IT_1515 = IT_1010*IT_1214;
    const ccomplex_t IT_1516 = (-0.125)*IT_1515;
    const ccomplex_t IT_1517 = 4*IT_0865;
    const ccomplex_t IT_1518 = 4*IT_0501;
    const ccomplex_t IT_1519 = IT_1517 + IT_1518;
    const ccomplex_t IT_1520 = IT_1207*IT_1519;
    const ccomplex_t IT_1521 = (-0.125)*IT_1520;
    const ccomplex_t IT_1522 = (-4)*IT_0183;
    const ccomplex_t IT_1523 = (-4)*IT_0221;
    const ccomplex_t IT_1524 = IT_1522 + IT_1523;
    const ccomplex_t IT_1525 = IT_1233*IT_1524;
    const ccomplex_t IT_1526 = (-0.125)*IT_1525;
    const ccomplex_t IT_1527 = IT_0743*IT_1226;
    const ccomplex_t IT_1528 = (-0.125)*IT_1527;
    const ccomplex_t IT_1529 = (-8)*IT_0300;
    const ccomplex_t IT_1530 = (-0.125)*IT_1529;
    const ccomplex_t IT_1531 = IT_0330 + IT_1530;
    const ccomplex_t IT_1532 = IT_1217*IT_1531;
    const ccomplex_t IT_1533 = IT_0765*IT_0960 + IT_1143*IT_1200 + (-8)
      *IT_1501 + (-8)*IT_1503 + (-8)*IT_1505 + (-8)*IT_1507 + (-8)*IT_1510 + (-8
      )*IT_1512 + (-8)*IT_1514 + (-8)*IT_1516 + (-8)*IT_1521 + (-8)*IT_1526 + (
      -8)*IT_1528 + (-8)*IT_1532;
    const ccomplex_t IT_1534 = IT_0765*IT_0960 + IT_1171*(IT_1174 + IT_1179 +
       IT_1180) + IT_1143*IT_1239 + (-8)*IT_1501 + (-8)*IT_1503 + (-8)*IT_1505 +
       (-8)*IT_1507 + (-8)*IT_1510 + (-8)*IT_1512 + (-8)*IT_1514 + (-8)*IT_1516 
      + (-8)*IT_1521 + (-8)*IT_1526 + (-8)*IT_1528 + (-8)*IT_1532;
    const ccomplex_t IT_1535 = IT_0017*IT_0025*IT_0067 + IT_0076*(IT_0155
      *IT_0158 + IT_0067*(IT_0111 + IT_0154 + IT_0160) + (conj(IT_0017) + conj
      (IT_0183) + conj(IT_0221))*IT_0246 + (IT_0017 + IT_0183 + IT_0221)*IT_0249
      ) + IT_0263*(IT_0158*IT_0297 + IT_0067*IT_0299 + IT_0249*(IT_0300 +
       IT_0318 + IT_0329 + IT_0330) + IT_0246*(conj(IT_0300) + conj(IT_0318) +
       conj(IT_0329) + conj(IT_0330))) + IT_0343*(IT_0297*IT_0346 + IT_0299
      *IT_0349 + (conj(IT_0300) + conj(IT_0318) + conj(IT_0329) + conj(IT_0330))
      *IT_0352 + (IT_0300 + IT_0318 + IT_0329 + IT_0330)*IT_0355) + IT_0365*
      (IT_0296*IT_0392 + IT_0298*IT_0398 + IT_0399*IT_0421 + IT_0422*IT_0427) +
       4*IT_0433*(IT_0065*IT_0436 + conj(IT_0028)*IT_0439 + IT_0028*IT_0442 +
       0.25*conj(IT_0065)*IT_0445) + IT_0447*(IT_0249*IT_0452 + IT_0355*IT_0457 
      + IT_0465*IT_0466) + IT_0552 + (IT_0511 + IT_0560)*(IT_0562 + IT_0564) + 
      (IT_0534 + IT_0575)*(IT_0577 + IT_0578) + IT_0551*(IT_0579 + IT_0580 +
       IT_0581 + IT_0582 + IT_0583 + IT_0584 + IT_0585 + IT_0586) + (IT_0298
      *IT_0466 + IT_0579 + IT_0580 + IT_0581 + IT_0582 + IT_0583 + IT_0584 +
       IT_0585 + IT_0586)*IT_0594 + (IT_0521 + IT_0604)*(IT_0605 + IT_0606 +
       IT_0608 + IT_0610) + conj(IT_0447)*(IT_0246*IT_0452 + IT_0352*IT_0457 +
       IT_0465*IT_0611 + IT_0471*IT_0611 + IT_0614*IT_0615) + IT_0183*(IT_0466
      *IT_0604 + IT_0067*IT_0618 + IT_0625*IT_0626 + IT_0633*IT_0634) + IT_0666 
      + IT_0677*(IT_0296*(IT_0678 + IT_0679) + IT_0298*(IT_0680 + IT_0681) +
       IT_0399*(IT_0682 + IT_0683) + IT_0422*(IT_0684 + IT_0685)) + IT_0690
      *IT_0717 + (IT_0720 + IT_0723)*(IT_0725 + IT_0727) + 4*IT_0728*(0.25*conj
      (IT_0244)*IT_0437 + IT_0244*IT_0442 + 0.25*conj(IT_0222)*IT_0443 + IT_0028
      *IT_0730 + IT_0065*IT_0733) + 4*IT_0742*(conj(IT_0017)*IT_0390 + IT_0017
      *conj(IT_0390) + 0.25*IT_0652*IT_0743 + 0.25*IT_0543*IT_0744 + 0.25
      *IT_0644*IT_0745 + 0.25*IT_0648*IT_0746) + conj(IT_0017)*(IT_0025*IT_0158 
      + IT_0468*IT_0611 + IT_0748*IT_0749) + IT_0542*(IT_0561*IT_0644 + IT_0563
      *IT_0648 + IT_0611*IT_0652 + conj(IT_0390)*IT_0749 + IT_0390*IT_0752) +
       IT_0751*(IT_0543*IT_0563 + IT_0611*IT_0644 + IT_0466*IT_0648 + IT_0561
      *IT_0652 + conj(IT_0419)*IT_0749 + IT_0419*IT_0752) + conj(IT_0017)*
      (IT_0346*IT_0753 + IT_0352*IT_0754 + IT_0749*IT_0756) + IT_0492*(conj
      (IT_0419)*IT_0611 + IT_0644*IT_0749 + IT_0648*IT_0752 + IT_0652*IT_0757 +
       IT_0543*IT_0758) + conj(IT_0183)*(IT_0521*IT_0611 + IT_0604*IT_0611 +
       IT_0158*IT_0618 + IT_0615*IT_0625 + IT_0633*IT_0759 + IT_0645*IT_0761) +
       IT_0763*(IT_0615*IT_0644 + IT_0626*IT_0648 + IT_0652*IT_0765 + IT_0543
      *IT_0766 + conj(IT_0419)*IT_0767 + IT_0419*IT_0768) + IT_0777*(conj
      (IT_0390)*IT_0615 + IT_0390*IT_0626 + IT_0652*IT_0767 + IT_0543*IT_0768 +
       IT_0644*IT_0778 + IT_0648*IT_0779) + conj(IT_0501)*(IT_0511*IT_0611 +
       IT_0560*IT_0611 + IT_0352*IT_0754 + IT_0346*IT_0780 + IT_0615*IT_0787 +
       IT_0789*IT_0790) + conj(IT_0472)*(IT_0480*IT_0611 + IT_0803*IT_0806 +
       IT_0659*IT_0816 + IT_0615*IT_0818 + IT_0823*IT_0824 + IT_0835*IT_0836 +
       IT_0749*IT_0844 + IT_0790*IT_0852) + IT_0472*(IT_0660*IT_0816 + IT_0626
      *IT_0818 + IT_0752*IT_0844 + IT_0803*IT_0855 + IT_0823*IT_0856 + IT_0852
      *IT_0857) + IT_0872 + IT_0877*(IT_0644*IT_0824 + IT_0648*IT_0856 + IT_0652
      *IT_0878 + IT_0543*IT_0879 + conj(IT_0419)*IT_0880 + IT_0419*IT_0881) +
       conj(IT_0707)*IT_0934 + IT_0707*IT_0943 + (conj(IT_0390)*IT_0824 +
       IT_0390*IT_0856 + IT_0652*IT_0880 + IT_0543*IT_0881 + IT_0644*IT_0933 +
       IT_0648*IT_0942)*IT_0951 + conj(IT_0300)*(IT_1061 + IT_0391*IT_1067) +
       IT_1080 + IT_1087*IT_1092 + IT_1095*(conj(IT_0300)*IT_0381 + conj(IT_0330
      )*IT_0392 + IT_0300*IT_0396 + IT_0330*IT_0398 + IT_0296*IT_0413 + IT_0298
      *IT_0425 + IT_0420*IT_1096 + IT_0426*IT_1097) + (-4)*conj(IT_0576)*IT_1119
       + IT_0764*IT_1126 + IT_1132*(IT_0296*IT_1137 + IT_0298*IT_1142 + IT_0399
      *IT_1143 + IT_0422*IT_1144) + (IT_0399*IT_1137 + IT_0422*IT_1142 + IT_0296
      *IT_1143 + IT_0298*IT_1144)*IT_1152 + conj(IT_0293)*IT_1173 + conj(IT_0330
      )*(IT_1061 + IT_1164*IT_1174 + IT_1178*(IT_1179 + IT_1180)) + IT_1181
      *IT_1186 + (IT_0281 + IT_0292)*IT_1252 + 8*(s_12 + IT_1256 + IT_1257 +
       IT_1258)*IT_1293 + IT_1309 + 8*(s_12 + IT_0825 + IT_0826 + IT_1311)
      *IT_1314 + 8*IT_1318*IT_1319 + 6*conj(IT_0160)*IT_1339 + 6*IT_0160*IT_1340
       + (-4)*IT_0576*IT_1345 + 4*IT_0523*(IT_0442*IT_1099 + IT_0730*IT_1107 +
       IT_0768*IT_1342 + IT_0634*IT_1343 + IT_0660*IT_1344 + IT_0466*IT_1346) +
       4*conj(IT_0523)*(0.25*IT_0534*IT_0611 + IT_0439*IT_1099 + IT_1101*IT_1107
       + 0.25*IT_1111*IT_1341 + IT_0767*IT_1342 + IT_0759*IT_1343 + IT_0659
      *IT_1344 + IT_0611*IT_1346) + (-4)*IT_0330*IT_1365 + (-4)*IT_1096*IT_1374 
      + (-4)*IT_0329*IT_1396 + (conj(IT_0281) + conj(IT_0292))*IT_1406 + (-6)
      *conj(IT_0493)*IT_1420 + IT_1422 + 6*conj(IT_0865)*IT_1427 + 6*IT_0865
      *IT_1428 + (-6)*conj(IT_0764)*IT_1431 + (-6)*conj(IT_0221)*IT_1441 + (-6)
      *IT_0221*IT_1442 + IT_0017*(IT_0349*IT_0753 + IT_1443) + IT_0501*(IT_0466
      *IT_0560 + IT_0349*IT_0780 + IT_0626*IT_0787 + IT_1443) + conj(IT_0111)*
      (IT_0246*IT_0618 + IT_0653*IT_0761 + IT_1444 + IT_1445 + IT_1446 + IT_1448
       + IT_1450 + IT_1451) + conj(IT_0154)*(IT_0246*IT_0618 + IT_0655*IT_0761 +
       IT_0625*IT_0765 + IT_0633*IT_0907 + IT_1444 + IT_1445 + IT_1446 + IT_1448
       + IT_1450 + IT_1451) + IT_0111*(IT_1452 + IT_1453 + IT_1454 + IT_1455 +
       IT_1457 + IT_1459 + IT_1460) + IT_0154*(IT_0656*IT_0761 + IT_0625*IT_0766
       + IT_0633*IT_0940 + IT_1452 + IT_1453 + IT_1454 + IT_1455 + IT_1457 +
       IT_1459 + IT_1460) + IT_0293*(IT_0426*IT_1067 + IT_1467 + IT_1469 +
       IT_1470) + IT_1097*(IT_1467 + IT_1035*IT_1468 + IT_1469 + IT_1470) + conj
      (IT_0692)*(IT_0835*IT_1279 + IT_1471 + IT_1472 + IT_1473 + IT_1474 +
       IT_1475 + IT_1476 + IT_1477 + IT_1478 + IT_1479 + IT_1480 + IT_1481 +
       IT_1482 + IT_1483) + conj(IT_0696)*(IT_0835*IT_1287 + IT_1471 + IT_1472 +
       IT_1473 + IT_1474 + IT_1475 + IT_1476 + IT_1477 + IT_1478 + IT_1479 +
       IT_1480 + IT_1481 + IT_1482 + IT_1483) + conj(IT_0700)*(IT_0835*(IT_1278 
      + IT_1286 + IT_1291) + IT_1471 + IT_1472 + IT_1473 + IT_1474 + IT_1475 +
       IT_1476 + IT_1477 + IT_1478 + IT_1479 + IT_1480 + IT_1481 + IT_1482 +
       IT_1483) + IT_0692*(IT_0835*IT_1274 + IT_1484 + IT_1485 + IT_1486 +
       IT_1487 + IT_1488 + IT_1489 + IT_1490 + IT_1491 + IT_1492 + IT_1493 +
       IT_1494 + IT_1495 + IT_1496) + IT_0696*(IT_0835*IT_1275 + IT_1484 +
       IT_1485 + IT_1486 + IT_1487 + IT_1488 + IT_1489 + IT_1490 + IT_1491 +
       IT_1492 + IT_1493 + IT_1494 + IT_1495 + IT_1496) + IT_0700*(IT_0835
      *IT_1283 + IT_1484 + IT_1485 + IT_1486 + IT_1487 + IT_1488 + IT_1489 +
       IT_1490 + IT_1491 + IT_1492 + IT_1493 + IT_1494 + IT_1495 + IT_1496) +
       IT_0300*(IT_0941*IT_0972 + IT_0397*IT_1067 + IT_1048*IT_1356 + IT_1020
      *IT_1360 + IT_1009*IT_1363 + IT_1497 + IT_1498 + IT_1499) + (-8)*IT_0318
      *IT_1500 + conj(IT_0318)*IT_1533 + conj(IT_0329)*IT_1534;
    return create_ccomplex_return(IT_1535);
}

