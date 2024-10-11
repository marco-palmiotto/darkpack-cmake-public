#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_cb = param->V_cb;
    const creal_t V_tb = param->V_tb;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_ub = param->V_ub;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(m_W, -4);
    const ccomplex_t IT_0003 = pow(m_W, 2);
    const ccomplex_t IT_0004 = pow(m_Z, 2);
    const ccomplex_t IT_0005 = (-0.5)*IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = pow(m_Z, -2);
    const ccomplex_t IT_0010 = 0.5*IT_0004;
    const ccomplex_t IT_0011 = cpow(IT_0010, 2);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (1 + 0.5*IT_0008)*(IT_0003 + -IT_0012);
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*IT_0014*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = (-2)*IT_0017;
    const ccomplex_t IT_0020 = -IT_0018 + IT_0019;
    const ccomplex_t IT_0021 = cpow(IT_0020, 2);
    const ccomplex_t IT_0022 = cpow(IT_0003 + (-0.5)*IT_0004, 2);
    const ccomplex_t IT_0023 = -IT_0003 + 0.5*IT_0001*IT_0004*(IT_0003 + (-0.5
      )*IT_0004) + IT_0001*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0023, 2);
    const ccomplex_t IT_0025 = IT_0009*IT_0024;
    const ccomplex_t IT_0026 = (1 + 0.5*IT_0008)*(IT_0003 + -IT_0025);
    const ccomplex_t IT_0027 = IT_0001*IT_0009;
    const ccomplex_t IT_0028 = IT_0023*IT_0027;
    const ccomplex_t IT_0029 = IT_0006*IT_0028;
    const ccomplex_t IT_0030 = IT_0001*IT_0023;
    const ccomplex_t IT_0031 = IT_0006*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = IT_0010 + IT_0032;
    const ccomplex_t IT_0034 = IT_0029*IT_0033;
    const ccomplex_t IT_0035 = 2*IT_0017;
    const ccomplex_t IT_0036 = IT_0001*IT_0009*IT_0024;
    const ccomplex_t IT_0037 = -IT_0008 + -IT_0036;
    const ccomplex_t IT_0038 = IT_0002*IT_0009;
    const ccomplex_t IT_0039 = IT_0023*IT_0038;
    const ccomplex_t IT_0040 = IT_0010*IT_0039;
    const ccomplex_t IT_0041 = IT_0006*IT_0040;
    const ccomplex_t IT_0042 = IT_0037 + IT_0041;
    const ccomplex_t IT_0043 = 1 + IT_0042;
    const ccomplex_t IT_0044 = IT_0006*IT_0043;
    const ccomplex_t IT_0045 = IT_0017 + -IT_0018;
    const ccomplex_t IT_0046 = IT_0009*IT_0023;
    const ccomplex_t IT_0047 = IT_0010*IT_0046;
    const ccomplex_t IT_0048 = (1 + 0.5*IT_0008)*(IT_0006 + -IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = IT_0001*IT_0007;
    const ccomplex_t IT_0052 = (1 + 0.5*IT_0036)*(IT_0003 + -IT_0051);
    const ccomplex_t IT_0053 = cpow(IT_0006, 3);
    const ccomplex_t IT_0054 = IT_0002*IT_0053;
    const ccomplex_t IT_0055 = IT_0001*IT_0009*IT_0011;
    const ccomplex_t IT_0056 = IT_0006*IT_0055;
    const ccomplex_t IT_0057 = IT_0006*IT_0036;
    const ccomplex_t IT_0058 = IT_0002*IT_0007*IT_0009;
    const ccomplex_t IT_0059 = IT_0023*IT_0058;
    const ccomplex_t IT_0060 = IT_0010*IT_0059;
    const ccomplex_t IT_0061 = -IT_0008 + -IT_0055;
    const ccomplex_t IT_0062 = IT_0041 + IT_0061;
    const ccomplex_t IT_0063 = 1 + IT_0062;
    const ccomplex_t IT_0064 = IT_0006*IT_0063;
    const ccomplex_t IT_0065 = (IT_0003 + -IT_0051)*(1 + 0.5*IT_0055);
    const ccomplex_t IT_0066 = IT_0010*IT_0027;
    const ccomplex_t IT_0067 = IT_0006*IT_0066;
    const ccomplex_t IT_0068 = IT_0001*IT_0010;
    const ccomplex_t IT_0069 = IT_0006*IT_0068;
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = IT_0023 + IT_0070;
    const ccomplex_t IT_0072 = IT_0067*IT_0071;
    const ccomplex_t IT_0073 = 0.666666666666667*IT_0013*IT_0021 +
       0.666666666666667*IT_0018*(IT_0018*IT_0026 + IT_0034*IT_0035 + -IT_0044
      *IT_0045 + 2*IT_0020*IT_0050) + 0.666666666666667*IT_0035*(IT_0035*IT_0052
       + -IT_0045*(IT_0006 + IT_0047 + -IT_0054 + -IT_0056 + -IT_0057 + IT_0060)
       + -IT_0020*IT_0064) + 0.666666666666667*IT_0045*(IT_0045*IT_0065 +
       IT_0020*IT_0072);
    const ccomplex_t IT_0074 = cabs(-IT_0003 + IT_0001*IT_0022);
    const ccomplex_t IT_0075 = cpow(IT_0074, 0.5);
    const ccomplex_t IT_0076 = ( ((m_W
) >= (m_W + m_Z
)) ? (0.125*IT_0000
      *IT_0001*IT_0073*IT_0075
) : (0
));
    const ccomplex_t IT_0077 = pow(m_c, 2);
    const ccomplex_t IT_0078 = pow(m_b, 2);
    const ccomplex_t IT_0079 = cpow(IT_0003 + IT_0077 + -IT_0078, 2);
    const ccomplex_t IT_0080 = -IT_0077 + 0.25*IT_0001*(IT_0003 + -IT_0077 +
       IT_0078)*(IT_0003 + IT_0077 + -IT_0078) + 0.25*IT_0001*IT_0079;
    const ccomplex_t IT_0081 = 0.5*IT_0003 + 0.5*IT_0077 + (-0.5)*IT_0078;
    const ccomplex_t IT_0082 = 0.5*IT_0003 + (-0.5)*IT_0077 + 0.5*IT_0078;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = IT_0003*IT_0080;
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0083 + IT_0085;
    const ccomplex_t IT_0087 = IT_0001*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0016;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = cpow(IT_0089, 2);
    const ccomplex_t IT_0091 = cabs(-IT_0077 + 0.25*IT_0001*IT_0079);
    const ccomplex_t IT_0092 = cpow(IT_0091, 0.5);
    const ccomplex_t IT_0093 = ( ((m_W
) >= (m_b + m_c
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0080 + IT_0087)*IT_0090*IT_0092
) : (0
));
    const ccomplex_t IT_0094 = pow(m_d, 2);
    const ccomplex_t IT_0095 = cpow(IT_0003 + IT_0077 + -IT_0094, 2);
    const ccomplex_t IT_0096 = -IT_0077 + 0.25*IT_0001*(IT_0003 + -IT_0077 +
       IT_0094)*(IT_0003 + IT_0077 + -IT_0094) + 0.25*IT_0001*IT_0095;
    const ccomplex_t IT_0097 = 0.5*IT_0003 + 0.5*IT_0077 + (-0.5)*IT_0094;
    const ccomplex_t IT_0098 = 0.5*IT_0003 + (-0.5)*IT_0077 + 0.5*IT_0094;
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = IT_0003*IT_0096;
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0099 + IT_0101;
    const ccomplex_t IT_0103 = IT_0001*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*V_cd*e_em
      *IT_0016;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = cabs(-IT_0077 + 0.25*IT_0001*IT_0095);
    const ccomplex_t IT_0107 = cpow(IT_0106, 0.5);
    const ccomplex_t IT_0108 = ( ((m_W
) >= (m_c + m_d
)) ? (0.5*IT_0000
      *IT_0001*(IT_0096 + IT_0103)*IT_0105*conj(IT_0105)*IT_0107
) : (0
));
    const ccomplex_t IT_0109 = pow(m_e, 2);
    const ccomplex_t IT_0110 = cpow(IT_0003 + -IT_0109, 2);
    const ccomplex_t IT_0111 = 0.25*IT_0001*(IT_0003 + IT_0109)*(IT_0003 + 
      -IT_0109) + 0.25*IT_0001*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0003 + (-0.5)*IT_0109;
    const ccomplex_t IT_0113 = 0.5*IT_0003 + 0.5*IT_0109;
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = IT_0003*IT_0111;
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = IT_0114 + IT_0116;
    const ccomplex_t IT_0118 = IT_0001*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0016;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = cpow(IT_0120, 2);
    const ccomplex_t IT_0122 = cabs(0.25*IT_0001*IT_0110);
    const ccomplex_t IT_0123 = cpow(IT_0122, 0.5);
    const ccomplex_t IT_0124 = ( ((m_W
) >= (m_e
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0111 + IT_0118)*IT_0121*IT_0123
) : (0
));
    const ccomplex_t IT_0125 = pow(m_h, 2);
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = IT_0003 + IT_0126;
    const ccomplex_t IT_0128 = cpow(IT_0127, 2);
    const ccomplex_t IT_0129 = IT_0002*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*m_W*e_em*IT_0016;
    const ccomplex_t IT_0131 = cpow(IT_0130, 2);
    const ccomplex_t IT_0132 = cpow(IT_0003 + (-0.5)*IT_0125, 2);
    const ccomplex_t IT_0133 = cabs(-IT_0003 + IT_0001*IT_0132);
    const ccomplex_t IT_0134 = cpow(IT_0133, 0.5);
    const ccomplex_t IT_0135 = ( ((m_W
) >= (m_W + m_h
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*(1 + 0.5*IT_0129)*IT_0131*IT_0134
) :
       (0
));
    const ccomplex_t IT_0136 = pow(m_s, 2);
    const ccomplex_t IT_0137 = cpow(IT_0003 + IT_0077 + -IT_0136, 2);
    const ccomplex_t IT_0138 = -IT_0077 + 0.25*IT_0001*(IT_0003 + -IT_0077 +
       IT_0136)*(IT_0003 + IT_0077 + -IT_0136) + 0.25*IT_0001*IT_0137;
    const ccomplex_t IT_0139 = 0.5*IT_0003 + 0.5*IT_0077 + (-0.5)*IT_0136;
    const ccomplex_t IT_0140 = 0.5*IT_0003 + (-0.5)*IT_0077 + 0.5*IT_0136;
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = IT_0003*IT_0138;
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = IT_0141 + IT_0143;
    const ccomplex_t IT_0145 = IT_0001*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0016;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = cabs(-IT_0077 + 0.25*IT_0001*IT_0137);
    const ccomplex_t IT_0149 = cpow(IT_0148, 0.5);
    const ccomplex_t IT_0150 = ( ((m_W
) >= (m_c + m_s
)) ? (0.5*IT_0000
      *IT_0001*(IT_0138 + IT_0145)*IT_0147*conj(IT_0147)*IT_0149
) : (0
));
    const ccomplex_t IT_0151 = pow(m_t, 2);
    const ccomplex_t IT_0152 = cpow(IT_0003 + -IT_0078 + IT_0151, 2);
    const ccomplex_t IT_0153 = -IT_0151 + 0.25*IT_0001*(IT_0003 + -IT_0078 +
       IT_0151)*(IT_0003 + IT_0078 + -IT_0151) + 0.25*IT_0001*IT_0152;
    const ccomplex_t IT_0154 = 0.5*IT_0003 + (-0.5)*IT_0078 + 0.5*IT_0151;
    const ccomplex_t IT_0155 = 0.5*IT_0003 + 0.5*IT_0078 + (-0.5)*IT_0151;
    const ccomplex_t IT_0156 = IT_0154*IT_0155;
    const ccomplex_t IT_0157 = IT_0003*IT_0153;
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = IT_0156 + IT_0158;
    const ccomplex_t IT_0160 = IT_0001*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0016;
    const ccomplex_t IT_0162 = 0.5*IT_0161;
    const ccomplex_t IT_0163 = cpow(IT_0162, 2);
    const ccomplex_t IT_0164 = cabs(-IT_0151 + 0.25*IT_0001*IT_0152);
    const ccomplex_t IT_0165 = cpow(IT_0164, 0.5);
    const ccomplex_t IT_0166 = ( ((m_W
) >= (m_b + m_t
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0153 + IT_0160)*IT_0163*IT_0165
) : (0
));
    const ccomplex_t IT_0167 = cpow(IT_0003 + -IT_0094 + IT_0151, 2);
    const ccomplex_t IT_0168 = -IT_0151 + 0.25*IT_0001*(IT_0003 + -IT_0094 +
       IT_0151)*(IT_0003 + IT_0094 + -IT_0151) + 0.25*IT_0001*IT_0167;
    const ccomplex_t IT_0169 = 0.5*IT_0003 + (-0.5)*IT_0094 + 0.5*IT_0151;
    const ccomplex_t IT_0170 = 0.5*IT_0003 + 0.5*IT_0094 + (-0.5)*IT_0151;
    const ccomplex_t IT_0171 = IT_0169*IT_0170;
    const ccomplex_t IT_0172 = IT_0003*IT_0168;
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = IT_0171 + IT_0173;
    const ccomplex_t IT_0175 = IT_0001*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0016;
    const ccomplex_t IT_0177 = 0.5*IT_0176;
    const ccomplex_t IT_0178 = cabs(-IT_0151 + 0.25*IT_0001*IT_0167);
    const ccomplex_t IT_0179 = cpow(IT_0178, 0.5);
    const ccomplex_t IT_0180 = ( ((m_W
) >= (m_d + m_t
)) ? (0.5*IT_0000
      *IT_0001*(IT_0168 + IT_0175)*IT_0177*conj(IT_0177)*IT_0179
) : (0
));
    const ccomplex_t IT_0181 = cpow(IT_0003 + -IT_0136 + IT_0151, 2);
    const ccomplex_t IT_0182 = -IT_0151 + 0.25*IT_0001*(IT_0003 + -IT_0136 +
       IT_0151)*(IT_0003 + IT_0136 + -IT_0151) + 0.25*IT_0001*IT_0181;
    const ccomplex_t IT_0183 = 0.5*IT_0003 + (-0.5)*IT_0136 + 0.5*IT_0151;
    const ccomplex_t IT_0184 = 0.5*IT_0003 + 0.5*IT_0136 + (-0.5)*IT_0151;
    const ccomplex_t IT_0185 = IT_0183*IT_0184;
    const ccomplex_t IT_0186 = IT_0003*IT_0182;
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = IT_0185 + IT_0187;
    const ccomplex_t IT_0189 = IT_0001*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0016;
    const ccomplex_t IT_0191 = 0.5*IT_0190;
    const ccomplex_t IT_0192 = cabs(-IT_0151 + 0.25*IT_0001*IT_0181);
    const ccomplex_t IT_0193 = cpow(IT_0192, 0.5);
    const ccomplex_t IT_0194 = ( ((m_W
) >= (m_s + m_t
)) ? (0.5*IT_0000
      *IT_0001*(IT_0182 + IT_0189)*IT_0191*conj(IT_0191)*IT_0193
) : (0
));
    const ccomplex_t IT_0195 = pow(m_u, 2);
    const ccomplex_t IT_0196 = cpow(IT_0003 + -IT_0078 + IT_0195, 2);
    const ccomplex_t IT_0197 = -IT_0195 + 0.25*IT_0001*(IT_0003 + -IT_0078 +
       IT_0195)*(IT_0003 + IT_0078 + -IT_0195) + 0.25*IT_0001*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0003 + (-0.5)*IT_0078 + 0.5*IT_0195;
    const ccomplex_t IT_0199 = 0.5*IT_0003 + 0.5*IT_0078 + (-0.5)*IT_0195;
    const ccomplex_t IT_0200 = IT_0198*IT_0199;
    const ccomplex_t IT_0201 = IT_0003*IT_0197;
    const ccomplex_t IT_0202 = (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = IT_0200 + IT_0202;
    const ccomplex_t IT_0204 = IT_0001*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1.4142135623731)*V_ub*e_em
      *IT_0016;
    const ccomplex_t IT_0206 = 0.5*IT_0205;
    const ccomplex_t IT_0207 = cabs(-IT_0195 + 0.25*IT_0001*IT_0196);
    const ccomplex_t IT_0208 = cpow(IT_0207, 0.5);
    const ccomplex_t IT_0209 = ( ((m_W
) >= (m_b + m_u
)) ? (0.5*IT_0000
      *IT_0001*(IT_0197 + IT_0204)*IT_0206*conj(IT_0206)*IT_0208
) : (0
));
    const ccomplex_t IT_0210 = cpow(IT_0003 + -IT_0094 + IT_0195, 2);
    const ccomplex_t IT_0211 = -IT_0195 + 0.25*IT_0001*(IT_0003 + -IT_0094 +
       IT_0195)*(IT_0003 + IT_0094 + -IT_0195) + 0.25*IT_0001*IT_0210;
    const ccomplex_t IT_0212 = 0.5*IT_0003 + (-0.5)*IT_0094 + 0.5*IT_0195;
    const ccomplex_t IT_0213 = 0.5*IT_0003 + 0.5*IT_0094 + (-0.5)*IT_0195;
    const ccomplex_t IT_0214 = IT_0212*IT_0213;
    const ccomplex_t IT_0215 = IT_0003*IT_0211;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = IT_0001*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0016;
    const ccomplex_t IT_0220 = 0.5*IT_0219;
    const ccomplex_t IT_0221 = cpow(IT_0220, 2);
    const ccomplex_t IT_0222 = cabs(-IT_0195 + 0.25*IT_0001*IT_0210);
    const ccomplex_t IT_0223 = cpow(IT_0222, 0.5);
    const ccomplex_t IT_0224 = ( ((m_W
) >= (m_d + m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0211 + IT_0218)*IT_0221*IT_0223
) : (0
));
    const ccomplex_t IT_0225 = cpow(IT_0003 + -IT_0136 + IT_0195, 2);
    const ccomplex_t IT_0226 = -IT_0195 + 0.25*IT_0001*(IT_0003 + -IT_0136 +
       IT_0195)*(IT_0003 + IT_0136 + -IT_0195) + 0.25*IT_0001*IT_0225;
    const ccomplex_t IT_0227 = 0.5*IT_0003 + (-0.5)*IT_0136 + 0.5*IT_0195;
    const ccomplex_t IT_0228 = 0.5*IT_0003 + 0.5*IT_0136 + (-0.5)*IT_0195;
    const ccomplex_t IT_0229 = IT_0227*IT_0228;
    const ccomplex_t IT_0230 = IT_0003*IT_0226;
    const ccomplex_t IT_0231 = (-0.5)*IT_0230;
    const ccomplex_t IT_0232 = IT_0229 + IT_0231;
    const ccomplex_t IT_0233 = IT_0001*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0016;
    const ccomplex_t IT_0235 = 0.5*IT_0234;
    const ccomplex_t IT_0236 = cpow(IT_0235, 2);
    const ccomplex_t IT_0237 = cabs(-IT_0195 + 0.25*IT_0001*IT_0225);
    const ccomplex_t IT_0238 = cpow(IT_0237, 0.5);
    const ccomplex_t IT_0239 = ( ((m_W
) >= (m_s + m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0226 + IT_0233)*IT_0236*IT_0238
) : (0
));
    const ccomplex_t IT_0240 = pow(m_mu, 2);
    const ccomplex_t IT_0241 = cpow(IT_0003 + -IT_0240, 2);
    const ccomplex_t IT_0242 = 0.25*IT_0001*(IT_0003 + IT_0240)*(IT_0003 + 
      -IT_0240) + 0.25*IT_0001*IT_0241;
    const ccomplex_t IT_0243 = 0.5*IT_0003 + (-0.5)*IT_0240;
    const ccomplex_t IT_0244 = 0.5*IT_0003 + 0.5*IT_0240;
    const ccomplex_t IT_0245 = IT_0243*IT_0244;
    const ccomplex_t IT_0246 = IT_0003*IT_0242;
    const ccomplex_t IT_0247 = (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = IT_0245 + IT_0247;
    const ccomplex_t IT_0249 = IT_0001*IT_0248;
    const ccomplex_t IT_0250 = cabs(0.25*IT_0001*IT_0241);
    const ccomplex_t IT_0251 = cpow(IT_0250, 0.5);
    const ccomplex_t IT_0252 = ( ((m_W
) >= (m_mu
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0121*(IT_0242 + IT_0249)*IT_0251
) : (0
));
    const ccomplex_t IT_0253 = pow(m_tau, 2);
    const ccomplex_t IT_0254 = cpow(IT_0003 + -IT_0253, 2);
    const ccomplex_t IT_0255 = 0.25*IT_0001*(IT_0003 + IT_0253)*(IT_0003 + 
      -IT_0253) + 0.25*IT_0001*IT_0254;
    const ccomplex_t IT_0256 = 0.5*IT_0003 + (-0.5)*IT_0253;
    const ccomplex_t IT_0257 = 0.5*IT_0003 + 0.5*IT_0253;
    const ccomplex_t IT_0258 = IT_0256*IT_0257;
    const ccomplex_t IT_0259 = IT_0003*IT_0255;
    const ccomplex_t IT_0260 = (-0.5)*IT_0259;
    const ccomplex_t IT_0261 = IT_0258 + IT_0260;
    const ccomplex_t IT_0262 = IT_0001*IT_0261;
    const ccomplex_t IT_0263 = cabs(0.25*IT_0001*IT_0254);
    const ccomplex_t IT_0264 = cpow(IT_0263, 0.5);
    const ccomplex_t IT_0265 = ( ((m_W
) >= (m_tau
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0121*(IT_0255 + IT_0262)*IT_0264
) : (0
));
    return create_ccomplex_return(IT_0076 + IT_0093 + IT_0108 + IT_0124 +
       IT_0135 + IT_0150 + IT_0166 + IT_0180 + IT_0194 + IT_0209 + IT_0224 +
       IT_0239 + IT_0252 + IT_0265);
}

