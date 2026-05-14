#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_sc_L_to_c_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_sc_L_to_c_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_c*N_u2*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_c, 2);
    const ccomplex_t IT_0008 = pow(m_sc_L, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = cos(alpha);
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = sin(alpha);
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0023 = IT_0015*IT_0022;
    const ccomplex_t IT_0024 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = IT_0017 + IT_0021 + IT_0023 + IT_0025;
    const ccomplex_t IT_0027 = IT_0004*IT_0012;
    const ccomplex_t IT_0028 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = IT_0004*IT_0018;
    const ccomplex_t IT_0031 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = -IT_0029 + -IT_0032 + -IT_0034 + -IT_0036;
    const ccomplex_t IT_0038 = IT_0026 + IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = m_N_2*IT_0040;
    const ccomplex_t IT_0042 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0043 = IT_0015*IT_0042;
    const ccomplex_t IT_0044 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0045 = IT_0019*IT_0044;
    const ccomplex_t IT_0046 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0047 = IT_0015*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0049 = IT_0019*IT_0048;
    const ccomplex_t IT_0050 = IT_0043 + IT_0045 + IT_0047 + IT_0049;
    const ccomplex_t IT_0051 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0052 = IT_0027*IT_0051;
    const ccomplex_t IT_0053 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0054 = IT_0030*IT_0053;
    const ccomplex_t IT_0055 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0056 = IT_0027*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0058 = IT_0030*IT_0057;
    const ccomplex_t IT_0059 = -IT_0052 + -IT_0054 + -IT_0056 + -IT_0058;
    const ccomplex_t IT_0060 = IT_0050 + IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = m_N_1*IT_0062;
    const ccomplex_t IT_0064 = IT_0041 + IT_0063;
    const ccomplex_t IT_0065 = IT_0011*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0000
      *IT_0002*IT_0004*IT_0018;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0014;
    const ccomplex_t IT_0069 = 0.666666666666667*IT_0068;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = pow(m_N_1, 2);
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0007 + IT_0071 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = cos(beta);
    const ccomplex_t IT_0076 = cpow(IT_0013, -2);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002
      *IT_0004*IT_0007*IT_0012 + (-0.5)*m_W*(IT_0001*IT_0012 + IT_0018*IT_0075)*
      (IT_0004 + -1./3*IT_0003*IT_0076));
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_13 + IT_0007 + IT_0071 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_12 + IT_0007 + -IT_0008 + -IT_0071 
      + -reg_prop, -1);
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_c*e_em*IT_0000*IT_0002
      *IT_0004*IT_0012;
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0080*IT_0087;
    const ccomplex_t IT_0089 = IT_0085*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = m_c*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0098 = IT_0015*IT_0097;
    const ccomplex_t IT_0099 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0100 = IT_0019*IT_0099;
    const ccomplex_t IT_0101 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0102 = IT_0015*IT_0101;
    const ccomplex_t IT_0103 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0104 = IT_0019*IT_0103;
    const ccomplex_t IT_0105 = IT_0098 + IT_0100 + IT_0102 + IT_0104;
    const ccomplex_t IT_0106 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0107 = IT_0027*IT_0106;
    const ccomplex_t IT_0108 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0109 = IT_0030*IT_0108;
    const ccomplex_t IT_0110 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0111 = IT_0027*IT_0110;
    const ccomplex_t IT_0112 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0113 = IT_0030*IT_0112;
    const ccomplex_t IT_0114 = -IT_0107 + -IT_0109 + -IT_0111 + -IT_0113;
    const ccomplex_t IT_0115 = IT_0105 + IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = m_N_3*IT_0117;
    const ccomplex_t IT_0119 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0120 = IT_0015*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0122 = IT_0019*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0124 = IT_0015*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0126 = IT_0019*IT_0125;
    const ccomplex_t IT_0127 = IT_0120 + IT_0122 + IT_0124 + IT_0126;
    const ccomplex_t IT_0128 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0129 = IT_0027*IT_0128;
    const ccomplex_t IT_0130 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0131 = IT_0030*IT_0130;
    const ccomplex_t IT_0132 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0133 = IT_0027*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0135 = IT_0030*IT_0134;
    const ccomplex_t IT_0136 = -IT_0129 + -IT_0131 + -IT_0133 + -IT_0135;
    const ccomplex_t IT_0137 = IT_0127 + IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = m_N_1*IT_0139;
    const ccomplex_t IT_0141 = IT_0118 + IT_0140;
    const ccomplex_t IT_0142 = IT_0096*IT_0141;
    const ccomplex_t IT_0143 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0144 = IT_0080*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0147 = IT_0015*IT_0146;
    const ccomplex_t IT_0148 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0149 = IT_0019*IT_0148;
    const ccomplex_t IT_0150 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0151 = IT_0027*IT_0150;
    const ccomplex_t IT_0152 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0153 = IT_0030*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0147 + IT_0149 + 
      -IT_0151 + -IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = m_N_1*IT_0155;
    const ccomplex_t IT_0157 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0158 = IT_0015*IT_0157;
    const ccomplex_t IT_0159 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0160 = IT_0019*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0162 = IT_0027*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0164 = IT_0030*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0158 + IT_0160 + 
      -IT_0162 + -IT_0164);
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = m_N_1*IT_0166;
    const ccomplex_t IT_0168 = IT_0156 + IT_0167;
    const ccomplex_t IT_0169 = IT_0145*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*m_c*N_u4*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0171 = (-0.5)*IT_0170;
    const ccomplex_t IT_0172 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0173 = IT_0171*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0176 = IT_0015*IT_0175;
    const ccomplex_t IT_0177 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0178 = IT_0019*IT_0177;
    const ccomplex_t IT_0179 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0180 = IT_0015*IT_0179;
    const ccomplex_t IT_0181 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0182 = IT_0019*IT_0181;
    const ccomplex_t IT_0183 = IT_0176 + IT_0178 + IT_0180 + IT_0182;
    const ccomplex_t IT_0184 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0185 = IT_0027*IT_0184;
    const ccomplex_t IT_0186 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0187 = IT_0030*IT_0186;
    const ccomplex_t IT_0188 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0189 = IT_0027*IT_0188;
    const ccomplex_t IT_0190 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0191 = IT_0030*IT_0190;
    const ccomplex_t IT_0192 = -IT_0185 + -IT_0187 + -IT_0189 + -IT_0191;
    const ccomplex_t IT_0193 = IT_0183 + IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = m_N_4*IT_0195;
    const ccomplex_t IT_0197 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0198 = IT_0015*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0200 = IT_0019*IT_0199;
    const ccomplex_t IT_0201 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0202 = IT_0015*IT_0201;
    const ccomplex_t IT_0203 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0204 = IT_0019*IT_0203;
    const ccomplex_t IT_0205 = IT_0198 + IT_0200 + IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0207 = IT_0027*IT_0206;
    const ccomplex_t IT_0208 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0209 = IT_0030*IT_0208;
    const ccomplex_t IT_0210 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0211 = IT_0027*IT_0210;
    const ccomplex_t IT_0212 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0213 = IT_0030*IT_0212;
    const ccomplex_t IT_0214 = -IT_0207 + -IT_0209 + -IT_0211 + -IT_0213;
    const ccomplex_t IT_0215 = IT_0205 + IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = (-0.5)*IT_0216;
    const ccomplex_t IT_0218 = m_N_1*IT_0217;
    const ccomplex_t IT_0219 = IT_0196 + IT_0218;
    const ccomplex_t IT_0220 = IT_0174*IT_0219;
    const ccomplex_t IT_0221 = IT_0065 + IT_0074 + IT_0084 + (-2)*IT_0091 +
       IT_0142 + IT_0169 + IT_0220;
    const ccomplex_t IT_0222 = conj(N_B1)*e_em;
    const ccomplex_t IT_0223 = IT_0014*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = conj(N_W1)*e_em;
    const ccomplex_t IT_0226 = IT_0004*IT_0225;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*(IT_0224 + 3*IT_0227);
    const ccomplex_t IT_0229 = 0.166666666666667*IT_0228;
    const ccomplex_t IT_0230 = IT_0087*IT_0229;
    const ccomplex_t IT_0231 = IT_0085*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = m_c*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0235 = 0.5*IT_0234;
    const ccomplex_t IT_0236 = IT_0067*IT_0235;
    const ccomplex_t IT_0237 = IT_0072*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = IT_0078*IT_0229;
    const ccomplex_t IT_0240 = IT_0082*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = m_N_1*IT_0040;
    const ccomplex_t IT_0243 = m_N_2*IT_0062;
    const ccomplex_t IT_0244 = IT_0242 + IT_0243;
    const ccomplex_t IT_0245 = conj(N_B2)*e_em;
    const ccomplex_t IT_0246 = IT_0014*IT_0245;
    const ccomplex_t IT_0247 = 1.4142135623731*IT_0246;
    const ccomplex_t IT_0248 = conj(N_W2)*e_em;
    const ccomplex_t IT_0249 = IT_0004*IT_0248;
    const ccomplex_t IT_0250 = 1.4142135623731*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*(IT_0247 + 3*IT_0250);
    const ccomplex_t IT_0252 = 0.166666666666667*IT_0251;
    const ccomplex_t IT_0253 = IT_0009*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0244*IT_0254;
    const ccomplex_t IT_0256 = IT_0143*IT_0229;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = IT_0168*IT_0257;
    const ccomplex_t IT_0259 = conj(N_B3)*e_em;
    const ccomplex_t IT_0260 = IT_0014*IT_0259;
    const ccomplex_t IT_0261 = 1.4142135623731*IT_0260;
    const ccomplex_t IT_0262 = conj(N_W3)*e_em;
    const ccomplex_t IT_0263 = IT_0004*IT_0262;
    const ccomplex_t IT_0264 = 1.4142135623731*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*(IT_0261 + 3*IT_0264);
    const ccomplex_t IT_0266 = 0.166666666666667*IT_0265;
    const ccomplex_t IT_0267 = IT_0094*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = m_N_1*IT_0117;
    const ccomplex_t IT_0270 = m_N_3*IT_0139;
    const ccomplex_t IT_0271 = IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = IT_0268*IT_0271;
    const ccomplex_t IT_0273 = conj(N_B4)*e_em;
    const ccomplex_t IT_0274 = IT_0014*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = conj(N_W4)*e_em;
    const ccomplex_t IT_0277 = IT_0004*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*(IT_0275 + 3*IT_0278);
    const ccomplex_t IT_0280 = 0.166666666666667*IT_0279;
    const ccomplex_t IT_0281 = IT_0172*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = m_N_1*IT_0195;
    const ccomplex_t IT_0284 = m_N_4*IT_0217;
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = IT_0282*IT_0285;
    const ccomplex_t IT_0287 = 2*IT_0233 + -IT_0238 + -IT_0241 + -IT_0255 + 
      -IT_0258 + -IT_0272 + -IT_0286;
    const ccomplex_t IT_0288 = 6*s_13;
    const ccomplex_t IT_0289 = m_c*m_N_1;
    const ccomplex_t IT_0290 = 6*IT_0221;
    const ccomplex_t IT_0291 = IT_0155*IT_0257;
    const ccomplex_t IT_0292 = IT_0040*IT_0254;
    const ccomplex_t IT_0293 = IT_0117*IT_0268;
    const ccomplex_t IT_0294 = IT_0195*IT_0282;
    const ccomplex_t IT_0295 = -IT_0090 + IT_0291 + IT_0292 + IT_0293 + IT_0294;
    const ccomplex_t IT_0296 = m_c*s_14;
    const ccomplex_t IT_0297 = IT_0096*IT_0139;
    const ccomplex_t IT_0298 = IT_0174*IT_0217;
    const ccomplex_t IT_0299 = IT_0145*IT_0166;
    const ccomplex_t IT_0300 = IT_0011*IT_0062;
    const ccomplex_t IT_0301 = IT_0232 + -IT_0297 + -IT_0298 + -IT_0299 + 
      -IT_0300;
    const ccomplex_t IT_0302 = s_34*m_N_1;
    const ccomplex_t IT_0303 = 6*conj(IT_0221);
    const ccomplex_t IT_0304 = 6*conj(IT_0287);
    const ccomplex_t IT_0305 = pow(m_h, 2);
    const ccomplex_t IT_0306 = m_c*m_N_1*IT_0305;
    const ccomplex_t IT_0307 = s_13*IT_0305;
    const ccomplex_t IT_0308 = s_14*s_34;
    const ccomplex_t IT_0309 = (-2)*IT_0308;
    const ccomplex_t IT_0310 = IT_0307 + IT_0309;
    const ccomplex_t IT_0311 = 6*IT_0287;
    const ccomplex_t IT_0312 = 0.166666666666667*IT_0302;
    const ccomplex_t IT_0313 = (IT_0221*conj(IT_0221) + IT_0287*conj(IT_0287))
      *IT_0288 + conj(IT_0287)*IT_0289*IT_0290 + (IT_0287*IT_0289 + IT_0295
      *IT_0296 + IT_0301*IT_0302)*IT_0303 + (IT_0296*IT_0301 + IT_0295*IT_0302)
      *IT_0304 + 6*conj(IT_0301)*(IT_0295*IT_0306 + -IT_0301*IT_0310 +
       0.166666666666667*IT_0296*IT_0311 + IT_0290*IT_0312) + 6*conj(IT_0295)*
      (0.166666666666667*IT_0290*IT_0296 + IT_0301*IT_0306 + -IT_0295*IT_0310 +
       IT_0311*IT_0312);
    return create_ccomplex_return(IT_0313);
}

