#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sd_L_to_d_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sd_L_to_d_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
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
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
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
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t IT_0000 = sin(alpha);
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0000 + m_d*e_em*mu_h*IT_0001*IT_0003*IT_0004*IT_0006);
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0010;
    const ccomplex_t IT_0012 = (-0.333333333333333)*IT_0011;
    const ccomplex_t IT_0013 = IT_0008*IT_0012;
    const ccomplex_t IT_0014 = pow(m_d, 2);
    const ccomplex_t IT_0015 = pow(m_N_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0014 + IT_0015 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0001
      *IT_0003*IT_0006;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0001*IT_0003*IT_0006;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = pow(m_sd_L, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + IT_0014 + -IT_0015 + -IT_0024 
      + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = m_d*IT_0027;
    const ccomplex_t IT_0029 = sin(beta);
    const ccomplex_t IT_0030 = cpow(IT_0009, -2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0006*IT_0014 + (-0.5)*m_W*(IT_0000*IT_0002 + IT_0004*IT_0029)*
      (IT_0006 + 1./3*IT_0005*IT_0030));
    const ccomplex_t IT_0032 = IT_0022*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0014 + IT_0015 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0004*IT_0010;
    const ccomplex_t IT_0037 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = IT_0000*IT_0010;
    const ccomplex_t IT_0040 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0043 = IT_0036*IT_0042;
    const ccomplex_t IT_0044 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0045 = IT_0039*IT_0044;
    const ccomplex_t IT_0046 = IT_0038 + IT_0041 + IT_0043 + IT_0045;
    const ccomplex_t IT_0047 = IT_0004*IT_0006;
    const ccomplex_t IT_0048 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = IT_0000*IT_0006;
    const ccomplex_t IT_0051 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0054 = IT_0047*IT_0053;
    const ccomplex_t IT_0055 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0056 = IT_0050*IT_0055;
    const ccomplex_t IT_0057 = -IT_0049 + -IT_0052 + -IT_0054 + -IT_0056;
    const ccomplex_t IT_0058 = IT_0046 + IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = m_N_3*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0063 = IT_0036*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0065 = IT_0039*IT_0064;
    const ccomplex_t IT_0066 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0067 = IT_0036*IT_0066;
    const ccomplex_t IT_0068 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0069 = IT_0039*IT_0068;
    const ccomplex_t IT_0070 = IT_0063 + IT_0065 + IT_0067 + IT_0069;
    const ccomplex_t IT_0071 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0072 = IT_0047*IT_0071;
    const ccomplex_t IT_0073 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0074 = IT_0050*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0076 = IT_0047*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0078 = IT_0050*IT_0077;
    const ccomplex_t IT_0079 = -IT_0072 + -IT_0074 + -IT_0076 + -IT_0078;
    const ccomplex_t IT_0080 = IT_0070 + IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = m_N_2*IT_0082;
    const ccomplex_t IT_0084 = IT_0061 + IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0001*IT_0003*IT_0006;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_23 + IT_0014 + IT_0024 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0084*IT_0089;
    const ccomplex_t IT_0091 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0092 = IT_0036*IT_0091;
    const ccomplex_t IT_0093 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0094 = IT_0039*IT_0093;
    const ccomplex_t IT_0095 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0096 = IT_0047*IT_0095;
    const ccomplex_t IT_0097 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0098 = IT_0050*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0092 + IT_0094 + 
      -IT_0096 + -IT_0098);
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = m_N_2*IT_0100;
    const ccomplex_t IT_0102 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0103 = IT_0036*IT_0102;
    const ccomplex_t IT_0104 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0105 = IT_0039*IT_0104;
    const ccomplex_t IT_0106 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0107 = IT_0047*IT_0106;
    const ccomplex_t IT_0108 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0109 = IT_0050*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + IT_0105 + 
      -IT_0107 + -IT_0109);
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = m_N_2*IT_0111;
    const ccomplex_t IT_0113 = IT_0101 + IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_23 + IT_0014 + IT_0024 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0115 = IT_0022*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0113*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0001*IT_0003*IT_0006;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_23 + IT_0014 + IT_0024 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0124 = IT_0036*IT_0123;
    const ccomplex_t IT_0125 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0126 = IT_0039*IT_0125;
    const ccomplex_t IT_0127 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0128 = IT_0036*IT_0127;
    const ccomplex_t IT_0129 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0130 = IT_0039*IT_0129;
    const ccomplex_t IT_0131 = IT_0124 + IT_0126 + IT_0128 + IT_0130;
    const ccomplex_t IT_0132 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0133 = IT_0047*IT_0132;
    const ccomplex_t IT_0134 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0135 = IT_0050*IT_0134;
    const ccomplex_t IT_0136 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0137 = IT_0047*IT_0136;
    const ccomplex_t IT_0138 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0139 = IT_0050*IT_0138;
    const ccomplex_t IT_0140 = -IT_0133 + -IT_0135 + -IT_0137 + -IT_0139;
    const ccomplex_t IT_0141 = IT_0131 + IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = m_N_1*IT_0143;
    const ccomplex_t IT_0145 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0146 = IT_0036*IT_0145;
    const ccomplex_t IT_0147 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0148 = IT_0039*IT_0147;
    const ccomplex_t IT_0149 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0150 = IT_0036*IT_0149;
    const ccomplex_t IT_0151 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0152 = IT_0039*IT_0151;
    const ccomplex_t IT_0153 = IT_0146 + IT_0148 + IT_0150 + IT_0152;
    const ccomplex_t IT_0154 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0155 = IT_0047*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0157 = IT_0050*IT_0156;
    const ccomplex_t IT_0158 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0159 = IT_0047*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0161 = IT_0050*IT_0160;
    const ccomplex_t IT_0162 = -IT_0155 + -IT_0157 + -IT_0159 + -IT_0161;
    const ccomplex_t IT_0163 = IT_0153 + IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = m_N_2*IT_0165;
    const ccomplex_t IT_0167 = IT_0144 + IT_0166;
    const ccomplex_t IT_0168 = IT_0122*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0001*IT_0003*IT_0006;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = cpow((-2)*s_23 + IT_0014 + IT_0024 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0172 = IT_0170*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0175 = IT_0036*IT_0174;
    const ccomplex_t IT_0176 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0177 = IT_0039*IT_0176;
    const ccomplex_t IT_0178 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0179 = IT_0036*IT_0178;
    const ccomplex_t IT_0180 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0181 = IT_0039*IT_0180;
    const ccomplex_t IT_0182 = IT_0175 + IT_0177 + IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0184 = IT_0047*IT_0183;
    const ccomplex_t IT_0185 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0186 = IT_0050*IT_0185;
    const ccomplex_t IT_0187 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0188 = IT_0047*IT_0187;
    const ccomplex_t IT_0189 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0190 = IT_0050*IT_0189;
    const ccomplex_t IT_0191 = -IT_0184 + -IT_0186 + -IT_0188 + -IT_0190;
    const ccomplex_t IT_0192 = IT_0182 + IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = m_N_4*IT_0194;
    const ccomplex_t IT_0196 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0197 = IT_0036*IT_0196;
    const ccomplex_t IT_0198 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0199 = IT_0039*IT_0198;
    const ccomplex_t IT_0200 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0201 = IT_0036*IT_0200;
    const ccomplex_t IT_0202 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0203 = IT_0039*IT_0202;
    const ccomplex_t IT_0204 = IT_0197 + IT_0199 + IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0206 = IT_0047*IT_0205;
    const ccomplex_t IT_0207 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0208 = IT_0050*IT_0207;
    const ccomplex_t IT_0209 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0210 = IT_0047*IT_0209;
    const ccomplex_t IT_0211 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0212 = IT_0050*IT_0211;
    const ccomplex_t IT_0213 = -IT_0206 + -IT_0208 + -IT_0210 + -IT_0212;
    const ccomplex_t IT_0214 = IT_0204 + IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = m_N_2*IT_0216;
    const ccomplex_t IT_0218 = IT_0195 + IT_0217;
    const ccomplex_t IT_0219 = IT_0173*IT_0218;
    const ccomplex_t IT_0220 = IT_0018 + (-2)*IT_0028 + IT_0035 + IT_0090 +
       IT_0117 + IT_0168 + IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0001*IT_0003*IT_0006;
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = IT_0008*IT_0222;
    const ccomplex_t IT_0224 = IT_0016*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = conj(N_B2)*e_em;
    const ccomplex_t IT_0227 = IT_0010*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = conj(N_W2)*e_em;
    const ccomplex_t IT_0230 = IT_0006*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0228 + (-3)*IT_0231);
    const ccomplex_t IT_0233 = 0.166666666666667*IT_0232;
    const ccomplex_t IT_0234 = IT_0020*IT_0233;
    const ccomplex_t IT_0235 = IT_0025*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = m_d*IT_0236;
    const ccomplex_t IT_0238 = m_N_2*IT_0060;
    const ccomplex_t IT_0239 = m_N_3*IT_0082;
    const ccomplex_t IT_0240 = IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = conj(N_B3)*e_em;
    const ccomplex_t IT_0242 = IT_0010*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = conj(N_W3)*e_em;
    const ccomplex_t IT_0245 = IT_0006*IT_0244;
    const ccomplex_t IT_0246 = 1.4142135623731*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*(IT_0243 + (-3)*IT_0246);
    const ccomplex_t IT_0248 = 0.166666666666667*IT_0247;
    const ccomplex_t IT_0249 = IT_0087*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = IT_0240*IT_0250;
    const ccomplex_t IT_0252 = IT_0114*IT_0233;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = IT_0113*IT_0253;
    const ccomplex_t IT_0255 = m_N_2*IT_0194;
    const ccomplex_t IT_0256 = m_N_4*IT_0216;
    const ccomplex_t IT_0257 = IT_0255 + IT_0256;
    const ccomplex_t IT_0258 = conj(N_B4)*e_em;
    const ccomplex_t IT_0259 = IT_0010*IT_0258;
    const ccomplex_t IT_0260 = 1.4142135623731*IT_0259;
    const ccomplex_t IT_0261 = conj(N_W4)*e_em;
    const ccomplex_t IT_0262 = IT_0006*IT_0261;
    const ccomplex_t IT_0263 = 1.4142135623731*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*(IT_0260 + (-3)*IT_0263);
    const ccomplex_t IT_0265 = 0.166666666666667*IT_0264;
    const ccomplex_t IT_0266 = IT_0171*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = IT_0257*IT_0267;
    const ccomplex_t IT_0269 = IT_0031*IT_0233;
    const ccomplex_t IT_0270 = IT_0033*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = m_N_2*IT_0143;
    const ccomplex_t IT_0273 = m_N_1*IT_0165;
    const ccomplex_t IT_0274 = IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = conj(N_B1)*e_em;
    const ccomplex_t IT_0276 = IT_0010*IT_0275;
    const ccomplex_t IT_0277 = 1.4142135623731*IT_0276;
    const ccomplex_t IT_0278 = conj(N_W1)*e_em;
    const ccomplex_t IT_0279 = IT_0006*IT_0278;
    const ccomplex_t IT_0280 = 1.4142135623731*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*(IT_0277 + (-3)*IT_0280);
    const ccomplex_t IT_0282 = 0.166666666666667*IT_0281;
    const ccomplex_t IT_0283 = IT_0120*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0274*IT_0284;
    const ccomplex_t IT_0286 = -IT_0225 + 2*IT_0237 + -IT_0251 + -IT_0254 + 
      -IT_0268 + -IT_0271 + -IT_0285;
    const ccomplex_t IT_0287 = 6*s_13;
    const ccomplex_t IT_0288 = m_d*m_N_2;
    const ccomplex_t IT_0289 = 6*IT_0220;
    const ccomplex_t IT_0290 = IT_0194*IT_0267;
    const ccomplex_t IT_0291 = IT_0060*IT_0250;
    const ccomplex_t IT_0292 = IT_0100*IT_0253;
    const ccomplex_t IT_0293 = IT_0143*IT_0284;
    const ccomplex_t IT_0294 = -IT_0027 + IT_0290 + IT_0291 + IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = m_d*s_14;
    const ccomplex_t IT_0296 = IT_0111*IT_0116;
    const ccomplex_t IT_0297 = IT_0173*IT_0216;
    const ccomplex_t IT_0298 = IT_0082*IT_0089;
    const ccomplex_t IT_0299 = IT_0122*IT_0165;
    const ccomplex_t IT_0300 = IT_0236 + -IT_0296 + -IT_0297 + -IT_0298 + 
      -IT_0299;
    const ccomplex_t IT_0301 = s_34*m_N_2;
    const ccomplex_t IT_0302 = 6*conj(IT_0220);
    const ccomplex_t IT_0303 = 6*conj(IT_0286);
    const ccomplex_t IT_0304 = pow(m_h, 2);
    const ccomplex_t IT_0305 = m_d*m_N_2*IT_0304;
    const ccomplex_t IT_0306 = s_13*IT_0304;
    const ccomplex_t IT_0307 = s_14*s_34;
    const ccomplex_t IT_0308 = (-2)*IT_0307;
    const ccomplex_t IT_0309 = IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = 6*IT_0286;
    const ccomplex_t IT_0311 = 0.166666666666667*IT_0301;
    const ccomplex_t IT_0312 = (IT_0220*conj(IT_0220) + IT_0286*conj(IT_0286))
      *IT_0287 + conj(IT_0286)*IT_0288*IT_0289 + (IT_0286*IT_0288 + IT_0294
      *IT_0295 + IT_0300*IT_0301)*IT_0302 + (IT_0295*IT_0300 + IT_0294*IT_0301)
      *IT_0303 + 6*conj(IT_0300)*(IT_0294*IT_0305 + -IT_0300*IT_0309 +
       0.166666666666667*IT_0295*IT_0310 + IT_0289*IT_0311) + 6*conj(IT_0294)*
      (0.166666666666667*IT_0289*IT_0295 + IT_0300*IT_0305 + -IT_0294*IT_0309 +
       IT_0310*IT_0311);
    return create_ccomplex_return(IT_0312);
}

