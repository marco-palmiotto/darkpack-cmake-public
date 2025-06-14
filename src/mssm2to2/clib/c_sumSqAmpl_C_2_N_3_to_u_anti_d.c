#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_u_anti_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_u_anti_d(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
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
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = s_13*s_24;
    const ccomplex_t IT_0001 = s_12*s_34;
    const ccomplex_t IT_0002 = s_14*s_23;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = IT_0001 + IT_0003;
    const ccomplex_t IT_0005 = IT_0000 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_u, 2);
    const ccomplex_t IT_0007 = pow(m_C_2, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0009*IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0017;
    const ccomplex_t IT_0019 = (-0.333333333333333)*IT_0018;
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = IT_0008*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0013;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d3)
      *e_em*IT_0009*IT_0011*IT_0013;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0022 + IT_0029;
    const ccomplex_t IT_0031 = pow(m_N_3, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + -IT_0007 + -IT_0031 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0033 = sin(beta);
    const ccomplex_t IT_0034 = IT_0013*IT_0033;
    const ccomplex_t IT_0035 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (-2)*IT_0036;
    const ccomplex_t IT_0038 = IT_0017*IT_0033;
    const ccomplex_t IT_0039 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = IT_0041 + IT_0044;
    const ccomplex_t IT_0046 = IT_0037 + IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0009*IT_0011*IT_0013*IT_0033;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = IT_0032*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = pow(m_W, 2);
    const ccomplex_t IT_0055 = cpow((-2)*s_12 + -IT_0007 + -IT_0031 + IT_0054 
      + -reg_prop, -1);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0009*IT_0013;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0010*IT_0013;
    const ccomplex_t IT_0059 = IT_0035*IT_0058;
    const ccomplex_t IT_0060 = (-2)*IT_0059;
    const ccomplex_t IT_0061 = IT_0010*IT_0017;
    const ccomplex_t IT_0062 = IT_0039*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = IT_0042*IT_0058;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = IT_0063 + IT_0065;
    const ccomplex_t IT_0067 = IT_0060 + IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0057*IT_0069;
    const ccomplex_t IT_0071 = IT_0055*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = -IT_0053 + -IT_0072;
    const ccomplex_t IT_0074 = 6*IT_0073;
    const ccomplex_t IT_0075 = 6*conj(IT_0073);
    const ccomplex_t IT_0076 = cpow(IT_0033, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0009*IT_0013*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = N_B3*e_em;
    const ccomplex_t IT_0080 = IT_0017*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = N_W3*e_em;
    const ccomplex_t IT_0083 = IT_0013*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + (-3)*IT_0084);
    const ccomplex_t IT_0086 = 0.166666666666667*IT_0085;
    const ccomplex_t IT_0087 = IT_0078*IT_0086;
    const ccomplex_t IT_0088 = IT_0027*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0009*IT_0013;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0093 = IT_0034*IT_0092;
    const ccomplex_t IT_0094 = 2*IT_0093;
    const ccomplex_t IT_0095 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0096 = IT_0038*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0099 = IT_0034*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = IT_0097 + IT_0100;
    const ccomplex_t IT_0102 = IT_0094 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0091*IT_0104;
    const ccomplex_t IT_0106 = IT_0055*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0009*IT_0010*IT_0013*IT_0076;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = IT_0058*IT_0092;
    const ccomplex_t IT_0111 = 2*IT_0110;
    const ccomplex_t IT_0112 = IT_0061*IT_0095;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = IT_0058*IT_0098;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = IT_0111 + IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0109*IT_0119;
    const ccomplex_t IT_0121 = IT_0032*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0107 + IT_0122;
    const ccomplex_t IT_0124 = 6*IT_0123;
    const ccomplex_t IT_0125 = 6*conj(IT_0123);
    const ccomplex_t IT_0126 = m_u*s_24*m_C_2;
    const ccomplex_t IT_0127 = IT_0048*IT_0109;
    const ccomplex_t IT_0128 = IT_0032*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0069*IT_0091;
    const ccomplex_t IT_0131 = IT_0055*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = -IT_0129 + -IT_0132;
    const ccomplex_t IT_0134 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0007 + -IT_0031 + -reg_prop, -1);
    const ccomplex_t IT_0135 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0136 = IT_0013*IT_0135;
    const ccomplex_t IT_0137 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0138 = IT_0013*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0136 + (-0.5)*IT_0139);
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0013;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = IT_0141*IT_0143;
    const ccomplex_t IT_0145 = IT_0134*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0025*IT_0078;
    const ccomplex_t IT_0148 = IT_0027*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = cpow((-2)*s_23 + IT_0006 + IT_0031 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0152 = conj(N_B3)*e_em;
    const ccomplex_t IT_0153 = IT_0017*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = conj(N_W3)*e_em;
    const ccomplex_t IT_0156 = IT_0013*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + 3*IT_0157);
    const ccomplex_t IT_0159 = 0.166666666666667*IT_0158;
    const ccomplex_t IT_0160 = IT_0015*IT_0159;
    const ccomplex_t IT_0161 = IT_0151*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = 0.166666666666667*IT_0074 + IT_0162;
    const ccomplex_t IT_0164 = 0.166666666666667*IT_0075;
    const ccomplex_t IT_0165 = conj(IT_0162) + IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0009*IT_0011*IT_0013;
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = IT_0015*IT_0167;
    const ccomplex_t IT_0169 = IT_0008*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = IT_0023*IT_0086;
    const ccomplex_t IT_0172 = IT_0027*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = -IT_0170 + -IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0017;
    const ccomplex_t IT_0176 = 0.666666666666667*IT_0175;
    const ccomplex_t IT_0177 = IT_0078*IT_0176;
    const ccomplex_t IT_0178 = cpow((-2)*s_23 + IT_0006 + IT_0031 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0179 = IT_0177*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0013;
    const ccomplex_t IT_0182 = -IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0009*IT_0013*IT_0076;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = IT_0182*IT_0184;
    const ccomplex_t IT_0186 = IT_0151*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = IT_0180 + IT_0187;
    const ccomplex_t IT_0189 = 6*IT_0188;
    const ccomplex_t IT_0190 = IT_0124 + IT_0189;
    const ccomplex_t IT_0191 = 6*conj(IT_0188);
    const ccomplex_t IT_0192 = IT_0125 + IT_0191;
    const ccomplex_t IT_0193 = IT_0015*IT_0184;
    const ccomplex_t IT_0194 = IT_0151*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0050*IT_0119;
    const ccomplex_t IT_0198 = IT_0032*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = IT_0057*IT_0104;
    const ccomplex_t IT_0201 = IT_0055*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0199 + IT_0202;
    const ccomplex_t IT_0204 = 6*IT_0203;
    const ccomplex_t IT_0205 = 6*conj(IT_0203);
    const ccomplex_t IT_0206 = IT_0159*IT_0182;
    const ccomplex_t IT_0207 = IT_0151*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0009*IT_0013*IT_0076;
    const ccomplex_t IT_0210 = 0.5*IT_0209;
    const ccomplex_t IT_0211 = IT_0078*IT_0210;
    const ccomplex_t IT_0212 = IT_0178*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = -IT_0208 + -IT_0213;
    const ccomplex_t IT_0215 = 6*IT_0133;
    const ccomplex_t IT_0216 = 6*conj(IT_0133);
    const ccomplex_t IT_0217 = 12*IT_0150;
    const ccomplex_t IT_0218 = 12*conj(IT_0150);
    const ccomplex_t IT_0219 = 12*IT_0133;
    const ccomplex_t IT_0220 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0221 = IT_0013*IT_0220;
    const ccomplex_t IT_0222 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0223 = IT_0013*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0221 + 0.5*IT_0224);
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = IT_0143*IT_0226;
    const ccomplex_t IT_0228 = IT_0134*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = 12*IT_0123;
    const ccomplex_t IT_0231 = 12*IT_0188;
    const ccomplex_t IT_0232 = 12*conj(IT_0123);
    const ccomplex_t IT_0233 = 12*conj(IT_0188);
    const ccomplex_t IT_0234 = 12*IT_0174;
    const ccomplex_t IT_0235 = 12*conj(IT_0174);
    const ccomplex_t IT_0236 = 12*conj(IT_0133)*IT_0146 + 6*conj(IT_0150)
      *IT_0163 + 6*IT_0150*IT_0165 + conj(IT_0174)*IT_0190 + IT_0174*IT_0192 +
       conj(IT_0196)*IT_0204 + IT_0196*IT_0205 + conj(IT_0214)*IT_0215 + IT_0214
      *IT_0216 + conj(IT_0030)*IT_0217 + IT_0030*IT_0218 + conj(IT_0146)*IT_0219
       + conj(IT_0229)*(24*IT_0089 + IT_0230 + IT_0231) + IT_0229*(24*conj
      (IT_0089) + IT_0232 + IT_0233) + conj(IT_0089)*IT_0234 + IT_0089*IT_0235;
    const ccomplex_t IT_0237 = 24*IT_0229;
    const ccomplex_t IT_0238 = IT_0234 + IT_0237;
    const ccomplex_t IT_0239 = 12*IT_0030;
    const ccomplex_t IT_0240 = 12*IT_0089;
    const ccomplex_t IT_0241 = 12*IT_0214;
    const ccomplex_t IT_0242 = 24*IT_0146;
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = 12*IT_0196;
    const ccomplex_t IT_0245 = 12*IT_0162;
    const ccomplex_t IT_0246 = m_d*m_u*s_12;
    const ccomplex_t IT_0247 = 6*IT_0196;
    const ccomplex_t IT_0248 = 6*conj(IT_0196);
    const ccomplex_t IT_0249 = 6*IT_0214;
    const ccomplex_t IT_0250 = 6*conj(IT_0214);
    const ccomplex_t IT_0251 = 12*conj(IT_0196);
    const ccomplex_t IT_0252 = (-12)*IT_0123;
    const ccomplex_t IT_0253 = (-12)*conj(IT_0123);
    const ccomplex_t IT_0254 = (-12)*IT_0133;
    const ccomplex_t IT_0255 = (-12)*conj(IT_0133);
    const ccomplex_t IT_0256 = (-6)*IT_0203;
    const ccomplex_t IT_0257 = (-6)*conj(IT_0203);
    const ccomplex_t IT_0258 = (-6)*IT_0133;
    const ccomplex_t IT_0259 = (-6)*conj(IT_0133);
    const ccomplex_t IT_0260 = conj(IT_0188)*IT_0204 + IT_0188*IT_0205 + conj
      (IT_0162)*IT_0215 + IT_0162*IT_0216 + conj(IT_0146)*IT_0217 + IT_0146
      *IT_0218 + conj(IT_0229)*IT_0244 + conj(IT_0174)*IT_0247 + IT_0174*IT_0248
       + conj(IT_0150)*IT_0249 + IT_0150*IT_0250 + IT_0229*IT_0251 + conj
      (IT_0203)*IT_0252 + IT_0203*IT_0253 + conj(IT_0073)*IT_0254 + IT_0073
      *IT_0255 + conj(IT_0089)*IT_0256 + IT_0089*IT_0257 + conj(IT_0030)*IT_0258
       + IT_0030*IT_0259;
    const ccomplex_t IT_0261 = s_34*m_C_2*m_N_3;
    const ccomplex_t IT_0262 = 6*IT_0174;
    const ccomplex_t IT_0263 = 6*conj(IT_0174);
    const ccomplex_t IT_0264 = 12*conj(IT_0214);
    const ccomplex_t IT_0265 = (-12)*IT_0203;
    const ccomplex_t IT_0266 = (-12)*conj(IT_0203);
    const ccomplex_t IT_0267 = conj(IT_0162)*IT_0204 + IT_0162*IT_0205 + conj
      (IT_0146)*IT_0234 + IT_0146*IT_0235 + conj(IT_0229)*IT_0243 + conj(IT_0150
      )*IT_0247 + IT_0150*IT_0248 + conj(IT_0133)*(IT_0189 + IT_0252) + IT_0133*
      (IT_0191 + IT_0253) + conj(IT_0030)*IT_0256 + IT_0030*IT_0257 + conj
      (IT_0089)*IT_0258 + IT_0089*IT_0259 + conj(IT_0214)*IT_0262 + IT_0214
      *IT_0263 + IT_0229*(24*conj(IT_0146) + IT_0264) + conj(IT_0073)*IT_0265 +
       IT_0073*IT_0266;
    const ccomplex_t IT_0268 = m_d*m_u*m_C_2*m_N_3;
    const ccomplex_t IT_0269 = (-6)*IT_0123;
    const ccomplex_t IT_0270 = IT_0231 + IT_0269;
    const ccomplex_t IT_0271 = (-6)*conj(IT_0123);
    const ccomplex_t IT_0272 = IT_0233 + IT_0271;
    const ccomplex_t IT_0273 = 12*IT_0203;
    const ccomplex_t IT_0274 = 12*conj(IT_0203);
    const ccomplex_t IT_0275 = 12*conj(IT_0030);
    const ccomplex_t IT_0276 = (-6)*IT_0188;
    const ccomplex_t IT_0277 = (-6)*conj(IT_0188);
    const ccomplex_t IT_0278 = conj(IT_0030)*IT_0190 + IT_0030*IT_0192 + 24
      *conj(IT_0146)*IT_0196 + 24*IT_0146*conj(IT_0196) + IT_0150*(24*conj
      (IT_0229) + IT_0235) + conj(IT_0150)*IT_0238 + 6*conj(IT_0089)*(IT_0163 +
       0.166666666666667*IT_0239) + conj(IT_0214)*IT_0244 + IT_0214*IT_0251 +
       conj(IT_0162)*IT_0270 + IT_0162*IT_0272 + conj(IT_0133)*IT_0273 + IT_0133
      *IT_0274 + 6*IT_0089*(IT_0165 + 0.166666666666667*IT_0275) + conj(IT_0073)
      *(IT_0230 + IT_0276) + IT_0073*(IT_0232 + IT_0277);
    const ccomplex_t IT_0279 = m_d*s_13*m_N_3;
    const ccomplex_t IT_0280 = 12*IT_0073;
    const ccomplex_t IT_0281 = conj(IT_0150)*IT_0190 + IT_0150*IT_0192 +
       IT_0205*IT_0214 + IT_0204*conj(IT_0214) + conj(IT_0089)*IT_0217 + IT_0089
      *IT_0218 + (24*conj(IT_0030) + 12*conj(IT_0073) + 12*conj(IT_0162))
      *IT_0229 + conj(IT_0030)*IT_0234 + IT_0030*IT_0235 + conj(IT_0133)*IT_0247
       + IT_0133*IT_0248 + (conj(IT_0073) + conj(IT_0162))*IT_0262 + (IT_0073 +
       IT_0162)*IT_0263 + conj(IT_0146)*IT_0273 + IT_0146*IT_0274 + conj(IT_0229
      )*(24*IT_0030 + IT_0245 + IT_0280);
    const ccomplex_t IT_0282 = -IT_0000;
    const ccomplex_t IT_0283 = IT_0004 + IT_0282;
    const ccomplex_t IT_0284 = (-0.166666666666667)*IT_0089;
    const ccomplex_t IT_0285 = (-0.166666666666667)*conj(IT_0089);
    const ccomplex_t IT_0286 = IT_0001 + IT_0002 + IT_0282;
    const ccomplex_t IT_0287 = (-6)*conj(IT_0073);
    const ccomplex_t IT_0288 = (-0.166666666666667)*IT_0162;
    const ccomplex_t IT_0289 = m_u*s_14*m_N_3;
    const ccomplex_t IT_0290 = (-0.166666666666667)*conj(IT_0150);
    const ccomplex_t IT_0291 = (-0.166666666666667)*IT_0214;
    const ccomplex_t IT_0292 = (-0.166666666666667)*IT_0030;
    const ccomplex_t IT_0293 = (-0.166666666666667)*conj(IT_0030);
    const ccomplex_t IT_0294 = (-0.166666666666667)*IT_0229;
    const ccomplex_t IT_0295 = (-0.166666666666667)*conj(IT_0229);
    const ccomplex_t IT_0296 = (-0.166666666666667)*conj(IT_0162);
    const ccomplex_t IT_0297 = conj(IT_0133)*IT_0174 + IT_0133*conj(IT_0174) +
       conj(IT_0073)*IT_0196 + IT_0073*conj(IT_0196) + (-4)*conj(IT_0146)*
      (IT_0188 + 0.0416666666666667*IT_0240 + 0.0416666666666667*IT_0252) + (-4)
      *IT_0146*(0.5*conj(IT_0089) + conj(IT_0188) + 0.0416666666666667*IT_0253) 
      + (-0.166666666666667)*IT_0150*IT_0257 + (-0.166666666666667)*conj(IT_0214
      )*IT_0270 + IT_0250*IT_0284 + IT_0249*IT_0285 + IT_0251*IT_0288 + IT_0256
      *IT_0290 + IT_0272*IT_0291 + IT_0248*IT_0292 + IT_0247*IT_0293 + IT_0255
      *IT_0294 + IT_0254*IT_0295 + IT_0244*IT_0296;
    const ccomplex_t IT_0298 = m_d*s_23*m_C_2;
    const ccomplex_t IT_0299 = IT_0073*conj(IT_0214) + 2*conj(IT_0146)*
      (IT_0073 + (-2)*IT_0162 + (-0.0833333333333333)*IT_0239) + (
      -0.166666666666667)*conj(IT_0174)*IT_0256 + (-0.166666666666667)*IT_0174
      *IT_0257 + (-0.166666666666667)*IT_0150*IT_0259 + (-0.166666666666667)
      *conj(IT_0196)*IT_0270 + (-0.166666666666667)*IT_0196*IT_0272 + 2*IT_0146*
      (conj(IT_0073) + (-2)*conj(IT_0162) + (-0.0833333333333333)*IT_0275) +
       IT_0248*IT_0284 + IT_0247*IT_0285 + IT_0264*IT_0288 + IT_0258*IT_0290 +
       IT_0287*IT_0291 + IT_0250*IT_0292 + IT_0249*IT_0293 + IT_0266*IT_0294 +
       IT_0265*IT_0295 + IT_0241*IT_0296;
    const ccomplex_t IT_0300 = IT_0005*(conj(IT_0030)*IT_0074 + IT_0030
      *IT_0075 + conj(IT_0089)*IT_0124 + IT_0089*IT_0125) + IT_0126*IT_0236 +
       IT_0000*(conj(IT_0150)*IT_0217 + (24*IT_0174 + 48*IT_0229)*conj(IT_0229) 
      + conj(IT_0174)*IT_0238 + conj(IT_0030)*IT_0239 + conj(IT_0089)*IT_0240) +
       IT_0002*(conj(IT_0146)*(48*IT_0146 + 24*IT_0214) + conj(IT_0188)*IT_0231 
      + conj(IT_0214)*IT_0243 + conj(IT_0196)*IT_0244 + conj(IT_0162)*IT_0245) +
       IT_0246*IT_0260 + IT_0261*IT_0267 + IT_0268*IT_0278 + IT_0001*(conj
      (IT_0133)*IT_0219 + conj(IT_0123)*IT_0230 + conj(IT_0203)*IT_0273 + conj
      (IT_0073)*IT_0280) + IT_0279*IT_0281 + (-6)*IT_0283*(conj(IT_0030)*IT_0162
       + IT_0030*conj(IT_0162) + IT_0277*IT_0284 + IT_0276*IT_0285) + (-6)
      *IT_0286*(IT_0073*conj(IT_0162) + (-0.166666666666667)*conj(IT_0188)
      *IT_0269 + (-0.166666666666667)*IT_0188*IT_0271 + IT_0287*IT_0288) + (-6)
      *IT_0289*IT_0297 + (-6)*IT_0298*IT_0299;
    return create_ccomplex_return(IT_0300);
}

