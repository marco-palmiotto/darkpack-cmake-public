#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_4_to_u_anti_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_4_to_u_anti_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = m_u*s_24*m_C_1;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = N_B4*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = N_W4*e_em;
    const ccomplex_t IT_0014 = IT_0005*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0012 + (-3)*IT_0015);
    const ccomplex_t IT_0017 = 0.166666666666667*IT_0016;
    const ccomplex_t IT_0018 = IT_0007*IT_0017;
    const ccomplex_t IT_0019 = pow(m_u, 2);
    const ccomplex_t IT_0020 = pow(m_C_1, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0019 + IT_0020 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0025 = IT_0005*IT_0024;
    const ccomplex_t IT_0026 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0027 = IT_0005*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + 0.5*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0005;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = pow(m_N_4, 2);
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0020 + -IT_0034 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cos(beta);
    const ccomplex_t IT_0039 = cpow(IT_0038, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d4)
      *e_em*IT_0001*IT_0005*IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = IT_0007*IT_0041;
    const ccomplex_t IT_0043 = IT_0021*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0001*IT_0005*IT_0039;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B4)*e_em;
    const ccomplex_t IT_0049 = IT_0009*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W4)*e_em;
    const ccomplex_t IT_0052 = IT_0005*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + 3*IT_0053);
    const ccomplex_t IT_0055 = 0.166666666666667*IT_0054;
    const ccomplex_t IT_0056 = IT_0047*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0019 + IT_0034 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = 6*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0005;
    const ccomplex_t IT_0062 = IT_0017*IT_0061;
    const ccomplex_t IT_0063 = IT_0021*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0001*IT_0005*IT_0039;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0047*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_13 + IT_0019 + IT_0020 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = -IT_0064 + -IT_0070;
    const ccomplex_t IT_0072 = 12*IT_0023;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0001*IT_0005;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0005*IT_0038;
    const ccomplex_t IT_0076 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (-2)*IT_0077;
    const ccomplex_t IT_0079 = IT_0009*IT_0038;
    const ccomplex_t IT_0080 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0084 = IT_0075*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = IT_0082 + IT_0085;
    const ccomplex_t IT_0087 = IT_0078 + IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0074*IT_0089;
    const ccomplex_t IT_0091 = pow(m_W, 2);
    const ccomplex_t IT_0092 = cpow((-2)*s_12 + -IT_0020 + -IT_0034 + IT_0091 
      + -reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0001*IT_0003*IT_0005*IT_0038;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0002*IT_0005;
    const ccomplex_t IT_0098 = IT_0076*IT_0097;
    const ccomplex_t IT_0099 = (-2)*IT_0098;
    const ccomplex_t IT_0100 = IT_0002*IT_0009;
    const ccomplex_t IT_0101 = IT_0080*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = IT_0083*IT_0097;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = IT_0102 + IT_0104;
    const ccomplex_t IT_0106 = IT_0099 + IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = IT_0096*IT_0108;
    const ccomplex_t IT_0110 = cpow((-2)*s_12 + -IT_0020 + -IT_0034 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = -IT_0094 + -IT_0112;
    const ccomplex_t IT_0114 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0115 = IT_0005*IT_0114;
    const ccomplex_t IT_0116 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0117 = IT_0005*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + (-0.5)*IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = IT_0032*IT_0120;
    const ccomplex_t IT_0122 = IT_0035*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = 12*IT_0123;
    const ccomplex_t IT_0125 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0126 = IT_0097*IT_0125;
    const ccomplex_t IT_0127 = 2*IT_0126;
    const ccomplex_t IT_0128 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0129 = IT_0100*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0132 = IT_0097*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = IT_0130 + IT_0133;
    const ccomplex_t IT_0135 = IT_0127 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0074*IT_0137;
    const ccomplex_t IT_0139 = IT_0092*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = IT_0075*IT_0125;
    const ccomplex_t IT_0142 = 2*IT_0141;
    const ccomplex_t IT_0143 = IT_0079*IT_0128;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = IT_0075*IT_0131;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = IT_0144 + IT_0146;
    const ccomplex_t IT_0148 = IT_0142 + IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0096*IT_0150;
    const ccomplex_t IT_0152 = IT_0110*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = IT_0140 + IT_0153;
    const ccomplex_t IT_0155 = s_13*s_24;
    const ccomplex_t IT_0156 = s_12*s_34;
    const ccomplex_t IT_0157 = s_14*s_23;
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = IT_0156 + IT_0158;
    const ccomplex_t IT_0160 = IT_0155 + IT_0159;
    const ccomplex_t IT_0161 = 6*IT_0023;
    const ccomplex_t IT_0162 = 12*IT_0037;
    const ccomplex_t IT_0163 = IT_0000*IT_0162;
    const ccomplex_t IT_0164 = 12*IT_0059;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = IT_0055*IT_0166;
    const ccomplex_t IT_0168 = IT_0057*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0171 = 0.5*IT_0170;
    const ccomplex_t IT_0172 = IT_0007*IT_0171;
    const ccomplex_t IT_0173 = cpow((-2)*s_23 + IT_0019 + IT_0034 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0174 = IT_0172*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = -IT_0169 + -IT_0175;
    const ccomplex_t IT_0177 = 24*IT_0123;
    const ccomplex_t IT_0178 = 24*IT_0037;
    const ccomplex_t IT_0179 = m_d*s_13*m_N_4;
    const ccomplex_t IT_0180 = 12*conj(IT_0059);
    const ccomplex_t IT_0181 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0001*IT_0005;
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = IT_0137*IT_0182;
    const ccomplex_t IT_0184 = IT_0092*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0001*IT_0002*IT_0005*IT_0039;
    const ccomplex_t IT_0187 = 0.5*IT_0186;
    const ccomplex_t IT_0188 = IT_0150*IT_0187;
    const ccomplex_t IT_0189 = IT_0110*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = IT_0185 + IT_0190;
    const ccomplex_t IT_0192 = IT_0089*IT_0182;
    const ccomplex_t IT_0193 = IT_0092*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = IT_0108*IT_0187;
    const ccomplex_t IT_0196 = IT_0110*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = -IT_0194 + -IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0009;
    const ccomplex_t IT_0200 = (-0.333333333333333)*IT_0199;
    const ccomplex_t IT_0201 = IT_0047*IT_0200;
    const ccomplex_t IT_0202 = IT_0068*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = IT_0041*IT_0061;
    const ccomplex_t IT_0205 = IT_0021*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = IT_0203 + IT_0206;
    const ccomplex_t IT_0208 = m_d*s_23*m_C_1;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0047*IT_0210;
    const ccomplex_t IT_0212 = IT_0057*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = (-12)*IT_0123;
    const ccomplex_t IT_0216 = (-12)*IT_0037;
    const ccomplex_t IT_0217 = m_u*s_14*m_N_4;
    const ccomplex_t IT_0218 = 12*conj(IT_0023);
    const ccomplex_t IT_0219 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0009;
    const ccomplex_t IT_0220 = 0.666666666666667*IT_0219;
    const ccomplex_t IT_0221 = IT_0007*IT_0220;
    const ccomplex_t IT_0222 = IT_0173*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = IT_0166*IT_0210;
    const ccomplex_t IT_0225 = IT_0057*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = IT_0223 + IT_0226;
    const ccomplex_t IT_0228 = m_d*m_u*s_12;
    const ccomplex_t IT_0229 = (-6)*IT_0023;
    const ccomplex_t IT_0230 = s_34*m_C_1*m_N_4;
    const ccomplex_t IT_0231 = 24*conj(IT_0037);
    const ccomplex_t IT_0232 = -IT_0155;
    const ccomplex_t IT_0233 = IT_0159 + IT_0232;
    const ccomplex_t IT_0234 = 6*conj(IT_0023);
    const ccomplex_t IT_0235 = 6*IT_0228;
    const ccomplex_t IT_0236 = 6*IT_0208;
    const ccomplex_t IT_0237 = 6*IT_0000;
    const ccomplex_t IT_0238 = 6*IT_0230;
    const ccomplex_t IT_0239 = 6*IT_0179;
    const ccomplex_t IT_0240 = 12*IT_0157;
    const ccomplex_t IT_0241 = m_d*m_u*m_C_1*m_N_4;
    const ccomplex_t IT_0242 = 12*IT_0241;
    const ccomplex_t IT_0243 = 12*IT_0217;
    const ccomplex_t IT_0244 = 12*conj(IT_0037);
    const ccomplex_t IT_0245 = 24*conj(IT_0123);
    const ccomplex_t IT_0246 = (-6)*IT_0208;
    const ccomplex_t IT_0247 = (-6)*IT_0217;
    const ccomplex_t IT_0248 = IT_0180*IT_0208 + IT_0217*IT_0234 + conj
      (IT_0045)*IT_0235 + conj(IT_0207)*IT_0236 + conj(IT_0113)*IT_0237 + conj
      (IT_0071)*IT_0238 + conj(IT_0191)*IT_0239 + conj(IT_0176)*IT_0240 + conj
      (IT_0214)*IT_0242 + conj(IT_0227)*IT_0243 + IT_0230*IT_0244 + IT_0157
      *IT_0245 + conj(IT_0198)*IT_0246 + conj(IT_0154)*IT_0247;
    const ccomplex_t IT_0249 = conj(IT_0154) + conj(IT_0227);
    const ccomplex_t IT_0250 = 6*conj(IT_0059);
    const ccomplex_t IT_0251 = 12*IT_0155;
    const ccomplex_t IT_0252 = 12*IT_0000;
    const ccomplex_t IT_0253 = 12*conj(IT_0123);
    const ccomplex_t IT_0254 = IT_0179*IT_0218 + conj(IT_0176)*IT_0235 + conj
      (IT_0198)*IT_0237 + conj(IT_0214)*IT_0238 + conj(IT_0071)*IT_0242 + conj
      (IT_0113)*IT_0246 + conj(IT_0191)*IT_0247 + IT_0239*IT_0249 + IT_0000
      *IT_0250 + conj(IT_0045)*IT_0251 + conj(IT_0207)*IT_0252 + IT_0228*IT_0253;
    const ccomplex_t IT_0255 = 12*IT_0156;
    const ccomplex_t IT_0256 = (-12)*IT_0228;
    const ccomplex_t IT_0257 = (-12)*IT_0230;
    const ccomplex_t IT_0258 = (-12)*conj(IT_0037);
    const ccomplex_t IT_0259 = (-6)*conj(IT_0023);
    const ccomplex_t IT_0260 = (-6)*IT_0228;
    const ccomplex_t IT_0261 = conj(IT_0176)*IT_0237 + conj(IT_0227)*IT_0238 +
       conj(IT_0214)*IT_0239 + conj(IT_0191)*IT_0242 + conj(IT_0045)*IT_0246 +
       conj(IT_0071)*IT_0247 + IT_0228*IT_0250 + IT_0000*IT_0253 + conj(IT_0113)
      *IT_0255 + conj(IT_0198)*IT_0256 + conj(IT_0154)*IT_0257 + IT_0217*IT_0258
       + IT_0230*IT_0259 + conj(IT_0207)*IT_0260;
    const ccomplex_t IT_0262 = 12*IT_0179;
    const ccomplex_t IT_0263 = IT_0000*IT_0218 + IT_0155*IT_0231 + conj
      (IT_0214)*IT_0235 + conj(IT_0176)*IT_0238 + conj(IT_0198)*IT_0239 + conj
      (IT_0045)*IT_0242 + conj(IT_0191)*IT_0246 + conj(IT_0113)*IT_0247 +
       IT_0237*IT_0249 + IT_0179*IT_0250 + conj(IT_0071)*IT_0251 + IT_0230
      *IT_0253 + conj(IT_0207)*IT_0262;
    const ccomplex_t IT_0264 = 6*IT_0217;
    const ccomplex_t IT_0265 = 12*IT_0208;
    const ccomplex_t IT_0266 = IT_0180*IT_0217 + IT_0208*IT_0234 + conj
      (IT_0071)*IT_0235 + conj(IT_0191)*IT_0237 + conj(IT_0045)*IT_0238 + conj
      (IT_0113)*IT_0239 + conj(IT_0214)*IT_0240 + conj(IT_0176)*IT_0242 +
       IT_0228*IT_0244 + conj(IT_0154)*IT_0246 + conj(IT_0198)*IT_0247 + conj
      (IT_0207)*IT_0264 + conj(IT_0227)*IT_0265;
    const ccomplex_t IT_0267 = (-6)*IT_0230;
    const ccomplex_t IT_0268 = conj(IT_0227)*IT_0235 + conj(IT_0214)*IT_0237 +
       conj(IT_0176)*IT_0239 + conj(IT_0113)*IT_0242 + conj(IT_0071)*IT_0246 +
       conj(IT_0045)*IT_0247 + IT_0230*IT_0250 + IT_0179*IT_0253 + conj(IT_0191)
      *IT_0255 + conj(IT_0154)*IT_0256 + conj(IT_0198)*IT_0257 + IT_0208*IT_0258
       + IT_0228*IT_0259 + conj(IT_0207)*IT_0267;
    const ccomplex_t IT_0269 = IT_0156 + IT_0157 + IT_0232;
    const ccomplex_t IT_0270 = (-6)*IT_0059;
    const ccomplex_t IT_0271 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0272 = conj(IT_0045)*IT_0178 + IT_0161*(conj(IT_0059) 
      + conj(IT_0198)) + IT_0072*conj(IT_0207) + IT_0177*conj(IT_0214) + IT_0207
      *IT_0218 + IT_0180*IT_0227 + IT_0164*conj(IT_0227) + IT_0045*IT_0231 + 
      (IT_0059 + IT_0198)*IT_0234 + IT_0214*IT_0245 + conj(IT_0154)*IT_0270 +
       IT_0154*IT_0271;
    const ccomplex_t IT_0273 = 6*IT_0160;
    const ccomplex_t IT_0274 = (-12)*conj(IT_0123);
    const ccomplex_t IT_0275 = (-6)*IT_0241;
    const ccomplex_t IT_0276 = conj(IT_0045)*IT_0237 + conj(IT_0071)*IT_0239 +
       conj(IT_0154)*IT_0242 + IT_0179*IT_0244 + conj(IT_0176)*IT_0246 + conj
      (IT_0214)*IT_0247 + conj(IT_0198)*IT_0255 + conj(IT_0113)*IT_0256 + conj
      (IT_0191)*IT_0257 + IT_0269*IT_0271 + conj(IT_0207)*IT_0273 + IT_0208
      *IT_0274 + conj(IT_0227)*IT_0275;
    const ccomplex_t IT_0277 = 6*IT_0241;
    const ccomplex_t IT_0278 = IT_0179*IT_0231 + conj(IT_0176)*IT_0236 + conj
      (IT_0207)*IT_0251 + conj(IT_0045)*IT_0252 + IT_0208*IT_0253 + conj(IT_0113
      )*IT_0260 + conj(IT_0071)*IT_0262 + conj(IT_0214)*IT_0264 + conj(IT_0191)
      *IT_0267 + IT_0233*IT_0271 + conj(IT_0198)*IT_0273 + IT_0249*IT_0277;
    const ccomplex_t IT_0279 = (-6)*IT_0269;
    const ccomplex_t IT_0280 = IT_0000*IT_0244;
    const ccomplex_t IT_0281 = conj(IT_0071)*IT_0237;
    const ccomplex_t IT_0282 = conj(IT_0045)*IT_0239;
    const ccomplex_t IT_0283 = conj(IT_0207)*IT_0277;
    const ccomplex_t IT_0284 = conj(IT_0191)*IT_0235 + conj(IT_0113)*IT_0238 +
       conj(IT_0227)*IT_0240 + conj(IT_0176)*IT_0243 + IT_0217*IT_0245 + IT_0233
      *IT_0259 + conj(IT_0214)*IT_0265 + conj(IT_0198)*IT_0275 + conj(IT_0154)
      *IT_0279 + IT_0280 + IT_0281 + IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = IT_0160*IT_0234 + conj(IT_0198)*IT_0242 + conj
      (IT_0214)*IT_0246 + conj(IT_0176)*IT_0247 + conj(IT_0154)*IT_0255 + conj
      (IT_0191)*IT_0256 + conj(IT_0113)*IT_0257 + IT_0217*IT_0274 + conj(IT_0227
      )*IT_0279 + IT_0280 + IT_0281 + IT_0282 + IT_0283;
    const ccomplex_t IT_0286 = IT_0000*(24*conj(IT_0023)*IT_0037 + 24*IT_0023
      *conj(IT_0037) + conj(IT_0045)*IT_0060 + conj(IT_0071)*IT_0072 + conj
      (IT_0113)*IT_0124) + conj(IT_0154)*(IT_0160*IT_0161 + IT_0163) + IT_0157*
      (48*IT_0123*conj(IT_0123) + conj(IT_0059)*IT_0164 + conj(IT_0176)*IT_0177)
       + IT_0155*(48*IT_0037*conj(IT_0037) + conj(IT_0023)*IT_0072 + conj
      (IT_0071)*IT_0178) + IT_0179*(IT_0060*conj(IT_0071) + conj(IT_0045)
      *IT_0072 + conj(IT_0037)*IT_0164 + IT_0037*IT_0180 + IT_0124*conj(IT_0191)
       + IT_0162*conj(IT_0198) + IT_0178*conj(IT_0207)) + IT_0208*(24*conj
      (IT_0059)*IT_0123 + 24*IT_0059*conj(IT_0123) + IT_0164*conj(IT_0176) +
       IT_0124*conj(IT_0207) + IT_0161*conj(IT_0214) + conj(IT_0198)*IT_0215 +
       conj(IT_0191)*IT_0216) + IT_0217*(IT_0072*conj(IT_0123) + IT_0161*conj
      (IT_0176) + IT_0164*conj(IT_0214) + conj(IT_0154)*IT_0215 + conj(IT_0113)
      *IT_0216 + IT_0123*IT_0218 + IT_0177*conj(IT_0227)) + IT_0228*(IT_0060
      *conj(IT_0113) + conj(IT_0045)*IT_0124 + IT_0162*conj(IT_0214) + conj
      (IT_0191)*IT_0229) + IT_0230*(conj(IT_0071)*IT_0124 + IT_0162*conj(IT_0176
      ) + conj(IT_0123)*IT_0178 + IT_0060*conj(IT_0191) + conj(IT_0113)*IT_0229 
      + IT_0123*IT_0231) + conj(IT_0227)*(IT_0163 + IT_0229*IT_0233) + IT_0176
      *IT_0248 + IT_0045*IT_0254 + IT_0113*IT_0261 + IT_0071*IT_0263 + IT_0214
      *IT_0266 + IT_0191*IT_0268 + (conj(IT_0207)*IT_0233 + conj(IT_0198)
      *IT_0269)*IT_0270 + IT_0241*IT_0272 + IT_0198*IT_0276 + IT_0207*IT_0278 +
       IT_0227*IT_0284 + IT_0154*IT_0285;
    return create_ccomplex_return(IT_0286);
}

