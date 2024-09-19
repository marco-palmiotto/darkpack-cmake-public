#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sc_L_to_anti_c_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sc_L_to_anti_c_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sd_L = param->m_sd_L;
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
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = conj(N_B1)*e_em;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = 1.4142135623731*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = conj(N_W1)*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0005 + (-3)*IT_0010);
    const ccomplex_t IT_0012 = 0.166666666666667*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = sin(beta);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0007*IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0012*IT_0017;
    const ccomplex_t IT_0019 = pow(m_c, 2);
    const ccomplex_t IT_0020 = pow(m_sc_L, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0000*IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0024 = conj(N_B2)*e_em;
    const ccomplex_t IT_0025 = IT_0002*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W2)*e_em;
    const ccomplex_t IT_0028 = IT_0007*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0026 + (-3)*IT_0029);
    const ccomplex_t IT_0031 = 0.166666666666667*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0007*IT_0013*IT_0015;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0023*IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0038 = conj(N_B3)*e_em;
    const ccomplex_t IT_0039 = IT_0002*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = conj(N_W3)*e_em;
    const ccomplex_t IT_0042 = IT_0007*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0040 + (-3)*IT_0043);
    const ccomplex_t IT_0045 = 0.166666666666667*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0007*IT_0013*IT_0015;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0037*IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0052 = conj(N_B4)*e_em;
    const ccomplex_t IT_0053 = IT_0002*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = conj(N_W4)*e_em;
    const ccomplex_t IT_0056 = IT_0007*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + (-3)*IT_0057);
    const ccomplex_t IT_0059 = 0.166666666666667*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0007*IT_0013*IT_0015;
    const ccomplex_t IT_0061 = (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0051*IT_0062*IT_0063;
    const ccomplex_t IT_0065 = cos(beta);
    const ccomplex_t IT_0066 = cpow(IT_0065, -1);
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0007*IT_0013*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0017*IT_0068;
    const ccomplex_t IT_0070 = IT_0021*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = m_d*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0007*IT_0013*IT_0066;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0033*IT_0074;
    const ccomplex_t IT_0076 = IT_0035*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = m_d*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0007*IT_0013*IT_0066;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0047*IT_0080;
    const ccomplex_t IT_0082 = IT_0049*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = m_d*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0007*IT_0013*IT_0066;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0061*IT_0086;
    const ccomplex_t IT_0088 = IT_0063*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = m_d*IT_0089;
    const ccomplex_t IT_0091 = -IT_0022 + -IT_0036 + -IT_0050 + -IT_0064 +
       IT_0072 + IT_0078 + IT_0084 + IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *e_em*IT_0007*IT_0013*IT_0015*IT_0065;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = pow(m_d, 2);
    const ccomplex_t IT_0095 = (0 + _Complex_I*(-2.82842712474619))*V_cd*e_em
      *IT_0007*((m_W*IT_0014 + (-0.5)*IT_0013*IT_0015*IT_0019)*IT_0065 + -1./2
      *IT_0013*IT_0014*IT_0066*IT_0094);
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0093*IT_0096;
    const ccomplex_t IT_0098 = pow(m_sd_L, 2);
    const ccomplex_t IT_0099 = cpow((-2)*s_12 + -IT_0020 + -IT_0098 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0097*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*conj(V_cd)
      *e_em*IT_0007;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0020 + -IT_0098 + -reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*V_cd*e_em
      *IT_0007;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = m_c*IT_0108;
    const ccomplex_t IT_0110 = -IT_0107;
    const ccomplex_t IT_0111 = m_c*IT_0110;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0109 + IT_0112;
    const ccomplex_t IT_0114 = IT_0106*IT_0113;
    const ccomplex_t IT_0115 = pow(m_W, 2);
    const ccomplex_t IT_0116 = cpow((-2)*s_12 + -IT_0020 + -IT_0098 + IT_0115 
      + -reg_prop, -1);
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *e_em*IT_0007*IT_0013;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = cpow(IT_0065, 2);
    const ccomplex_t IT_0120 = cpow(IT_0014, 2);
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*V_cd*e_em
      *IT_0007*(IT_0013*(IT_0019 + -IT_0094) + m_W*(IT_0119 + -IT_0120));
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0118*IT_0122;
    const ccomplex_t IT_0124 = IT_0116*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = -IT_0101 + IT_0114 + -IT_0125;
    const ccomplex_t IT_0127 = N_B1*e_em;
    const ccomplex_t IT_0128 = IT_0002*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = N_W1*e_em;
    const ccomplex_t IT_0131 = IT_0007*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + 3*IT_0132);
    const ccomplex_t IT_0134 = 0.166666666666667*IT_0133;
    const ccomplex_t IT_0135 = IT_0068*IT_0134;
    const ccomplex_t IT_0136 = IT_0000*IT_0021*IT_0135;
    const ccomplex_t IT_0137 = N_B2*e_em;
    const ccomplex_t IT_0138 = IT_0002*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = N_W2*e_em;
    const ccomplex_t IT_0141 = IT_0007*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + 3*IT_0142);
    const ccomplex_t IT_0144 = 0.166666666666667*IT_0143;
    const ccomplex_t IT_0145 = IT_0074*IT_0144;
    const ccomplex_t IT_0146 = IT_0023*IT_0035*IT_0145;
    const ccomplex_t IT_0147 = N_B3*e_em;
    const ccomplex_t IT_0148 = IT_0002*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = N_W3*e_em;
    const ccomplex_t IT_0151 = IT_0007*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + 3*IT_0152);
    const ccomplex_t IT_0154 = 0.166666666666667*IT_0153;
    const ccomplex_t IT_0155 = IT_0080*IT_0154;
    const ccomplex_t IT_0156 = IT_0037*IT_0049*IT_0155;
    const ccomplex_t IT_0157 = N_B4*e_em;
    const ccomplex_t IT_0158 = IT_0002*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = N_W4*e_em;
    const ccomplex_t IT_0161 = IT_0007*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0159 + 3*IT_0162);
    const ccomplex_t IT_0164 = 0.166666666666667*IT_0163;
    const ccomplex_t IT_0165 = IT_0086*IT_0164;
    const ccomplex_t IT_0166 = IT_0051*IT_0063*IT_0165;
    const ccomplex_t IT_0167 = IT_0012*IT_0134;
    const ccomplex_t IT_0168 = IT_0021*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = m_d*IT_0169;
    const ccomplex_t IT_0171 = IT_0031*IT_0144;
    const ccomplex_t IT_0172 = IT_0035*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = m_d*IT_0173;
    const ccomplex_t IT_0175 = IT_0045*IT_0154;
    const ccomplex_t IT_0176 = IT_0049*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = m_d*IT_0177;
    const ccomplex_t IT_0179 = IT_0059*IT_0164;
    const ccomplex_t IT_0180 = IT_0063*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = m_d*IT_0181;
    const ccomplex_t IT_0183 = pow(m_sG, 2);
    const ccomplex_t IT_0184 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -IT_0183 +
       reg_prop, -1);
    const ccomplex_t IT_0185 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0186 = cpow(IT_0185, 2);
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = IT_0184*IT_0187;
    const ccomplex_t IT_0189 = m_d*IT_0188;
    const ccomplex_t IT_0190 = -IT_0136 + -IT_0146 + -IT_0156 + -IT_0166 +
       IT_0170 + IT_0174 + IT_0178 + IT_0182 + (-0.166666666666667)*IT_0189;
    const ccomplex_t IT_0191 = m_d*IT_0108;
    const ccomplex_t IT_0192 = m_d*IT_0110;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0106*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_cd)
      *e_em*IT_0007*IT_0013*IT_0014*IT_0066;
    const ccomplex_t IT_0197 = 0.5*IT_0196;
    const ccomplex_t IT_0198 = IT_0096*IT_0197;
    const ccomplex_t IT_0199 = IT_0099*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1.4142135623731)*m_d*conj(V_cd)
      *e_em*IT_0007*IT_0013;
    const ccomplex_t IT_0202 = (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = IT_0122*IT_0202;
    const ccomplex_t IT_0204 = IT_0116*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = 0.5*IT_0189 + -IT_0195 + -IT_0200 + -IT_0205;
    const ccomplex_t IT_0207 = 6*s_34;
    const ccomplex_t IT_0208 = 18*s_34;
    const ccomplex_t IT_0209 = s_34*IT_0098;
    const ccomplex_t IT_0210 = (-6)*IT_0209;
    const ccomplex_t IT_0211 = s_13*s_14;
    const ccomplex_t IT_0212 = 12*IT_0211;
    const ccomplex_t IT_0213 = IT_0210 + IT_0212;
    const ccomplex_t IT_0214 = IT_0106*IT_0110;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = -IT_0173 + -IT_0177 + -IT_0181;
    const ccomplex_t IT_0217 = -IT_0169 + 0.166666666666667*IT_0188;
    const ccomplex_t IT_0218 = IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = conj(IT_0216) + conj(IT_0217);
    const ccomplex_t IT_0220 = (-0.5)*IT_0188;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = (-18)*IT_0209;
    const ccomplex_t IT_0223 = 36*IT_0211;
    const ccomplex_t IT_0224 = IT_0222 + IT_0223;
    const ccomplex_t IT_0225 = -IT_0071 + -IT_0077 + -IT_0083 + -IT_0089;
    const ccomplex_t IT_0226 = m_c*m_d*IT_0098;
    const ccomplex_t IT_0227 = (-18)*conj(IT_0217);
    const ccomplex_t IT_0228 = (-6)*conj(IT_0215);
    const ccomplex_t IT_0229 = (-18)*conj(IT_0216);
    const ccomplex_t IT_0230 = 6*IT_0220 + IT_0227 + IT_0228 + IT_0229;
    const ccomplex_t IT_0231 = (-6)*IT_0220;
    const ccomplex_t IT_0232 = m_c*s_14;
    const ccomplex_t IT_0233 = (-18)*IT_0190;
    const ccomplex_t IT_0234 = (-6)*IT_0206;
    const ccomplex_t IT_0235 = IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = (-18)*conj(IT_0190);
    const ccomplex_t IT_0237 = (-6)*conj(IT_0206);
    const ccomplex_t IT_0238 = IT_0236 + IT_0237;
    const ccomplex_t IT_0239 = (-6)*IT_0126;
    const ccomplex_t IT_0240 = (-18)*conj(IT_0091);
    const ccomplex_t IT_0241 = (-6)*conj(IT_0126);
    const ccomplex_t IT_0242 = (-18)*conj(IT_0126);
    const ccomplex_t IT_0243 = (-0.166666666666667)*IT_0242;
    const ccomplex_t IT_0244 = (-18)*IT_0126*conj(IT_0215) + IT_0091*IT_0230 +
       conj(IT_0225)*IT_0235 + IT_0225*IT_0238 + IT_0219*IT_0239 + IT_0218*
      (IT_0240 + IT_0241) + (-6)*IT_0215*(conj(IT_0091) + IT_0243) + IT_0231*
      (conj(IT_0091) + (-3)*IT_0126 + IT_0243);
    const ccomplex_t IT_0245 = m_c*m_d;
    const ccomplex_t IT_0246 = (-18)*IT_0206;
    const ccomplex_t IT_0247 = (-0.166666666666667)*IT_0246;
    const ccomplex_t IT_0248 = m_d*s_13;
    const ccomplex_t IT_0249 = 6*IT_0206;
    const ccomplex_t IT_0250 = 18*IT_0190;
    const ccomplex_t IT_0251 = 6*conj(IT_0206);
    const ccomplex_t IT_0252 = 18*conj(IT_0190);
    const ccomplex_t IT_0253 = 6*conj(IT_0190);
    const ccomplex_t IT_0254 = 18*conj(IT_0206);
    const ccomplex_t IT_0255 = 6*(IT_0190 + 3*IT_0206)*conj(IT_0215) + (18
      *conj(IT_0091) + 6*conj(IT_0126))*IT_0225 + (18*IT_0091 + 6*IT_0126)*conj
      (IT_0225) + IT_0219*(IT_0249 + IT_0250) + IT_0218*(IT_0251 + IT_0252) +
       IT_0215*(IT_0253 + IT_0254) + IT_0231*(IT_0190 + IT_0247 + (
      -0.166666666666667)*IT_0253 + (-0.166666666666667)*IT_0254);
    const ccomplex_t IT_0256 = (conj(IT_0091)*IT_0126 + IT_0091*conj(IT_0126) 
      + conj(IT_0190)*IT_0206 + IT_0190*conj(IT_0206))*IT_0207 + (IT_0091*conj
      (IT_0091) + IT_0126*conj(IT_0126) + IT_0190*conj(IT_0190) + IT_0206*conj
      (IT_0206))*IT_0208 + IT_0213*(conj(IT_0215)*IT_0218 + IT_0215*IT_0219 + 
      (IT_0216 + -conj(IT_0216) + IT_0217 + -conj(IT_0217))*IT_0221) + IT_0224*(
      (IT_0216 + IT_0217)*IT_0219 + (IT_0215 + IT_0220)*(conj(IT_0215) + IT_0221
      ) + IT_0225*conj(IT_0225)) + IT_0226*(IT_0225*IT_0230 + (-6)*conj(IT_0225)
      *(IT_0215 + 3*IT_0216 + 3*IT_0217 + -1./6*IT_0231)) + IT_0232*IT_0244 + (
      -6)*IT_0245*(IT_0126*(conj(IT_0190) + 3*conj(IT_0206)) + (
      -0.166666666666667)*conj(IT_0091)*IT_0235 + (-0.166666666666667)*IT_0091
      *IT_0238 + conj(IT_0126)*(IT_0190 + IT_0247)) + IT_0248*IT_0255;
    return create_ccomplex_return(IT_0256);
}

