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
    const ccomplex_t IT_0002 = pow(m_W, 2);
    const ccomplex_t IT_0003 = pow(m_Z, 2);
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = pow(m_W, -4);
    const ccomplex_t IT_0007 = cpow(IT_0005, 3);
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0005, 2);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = pow(m_Z, -2);
    const ccomplex_t IT_0012 = 0.5*IT_0003;
    const ccomplex_t IT_0013 = cpow(IT_0012, 2);
    const ccomplex_t IT_0014 = IT_0001*IT_0011*IT_0013;
    const ccomplex_t IT_0015 = (IT_0002 + -IT_0010)*(1 + 0.5*IT_0014);
    const ccomplex_t IT_0016 = IT_0005*IT_0014;
    const ccomplex_t IT_0017 = cpow(IT_0002 + (-0.5)*IT_0003, 2);
    const ccomplex_t IT_0018 = -IT_0002 + 0.5*IT_0001*IT_0003*(IT_0002 + (-0.5
      )*IT_0003) + IT_0001*IT_0017;
    const ccomplex_t IT_0019 = cpow(IT_0018, 2);
    const ccomplex_t IT_0020 = IT_0001*IT_0011*IT_0019;
    const ccomplex_t IT_0021 = IT_0005*IT_0020;
    const ccomplex_t IT_0022 = (IT_0002 + -IT_0010)*(1 + 0.5*IT_0020);
    const ccomplex_t IT_0023 = IT_0011*IT_0018;
    const ccomplex_t IT_0024 = IT_0012*IT_0023;
    const ccomplex_t IT_0025 = IT_0006*IT_0009*IT_0011;
    const ccomplex_t IT_0026 = IT_0018*IT_0025;
    const ccomplex_t IT_0027 = IT_0012*IT_0026;
    const ccomplex_t IT_0028 = cos(theta_W);
    const ccomplex_t IT_0029 = sin(theta_W);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*IT_0028*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = IT_0031 + -IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0033, 2);
    const ccomplex_t IT_0035 = IT_0006*IT_0009;
    const ccomplex_t IT_0036 = IT_0011*IT_0013;
    const ccomplex_t IT_0037 = (1 + 0.5*IT_0035)*(IT_0002 + -IT_0036);
    const ccomplex_t IT_0038 = IT_0011*IT_0019;
    const ccomplex_t IT_0039 = (1 + 0.5*IT_0035)*(IT_0002 + -IT_0038);
    const ccomplex_t IT_0040 = (IT_0005 + -IT_0024)*(1 + 0.5*IT_0035);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0001*IT_0011;
    const ccomplex_t IT_0044 = IT_0012*IT_0043;
    const ccomplex_t IT_0045 = IT_0005*IT_0044;
    const ccomplex_t IT_0046 = IT_0001*IT_0012;
    const ccomplex_t IT_0047 = IT_0005*IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = IT_0018 + IT_0048;
    const ccomplex_t IT_0050 = IT_0045*IT_0049;
    const ccomplex_t IT_0051 = -IT_0020 + -IT_0035;
    const ccomplex_t IT_0052 = IT_0006*IT_0011;
    const ccomplex_t IT_0053 = IT_0018*IT_0052;
    const ccomplex_t IT_0054 = IT_0012*IT_0053;
    const ccomplex_t IT_0055 = IT_0005*IT_0054;
    const ccomplex_t IT_0056 = IT_0051 + IT_0055;
    const ccomplex_t IT_0057 = 1 + IT_0056;
    const ccomplex_t IT_0058 = IT_0005*IT_0057;
    const ccomplex_t IT_0059 = -IT_0014 + -IT_0035;
    const ccomplex_t IT_0060 = IT_0055 + IT_0059;
    const ccomplex_t IT_0061 = 1 + IT_0060;
    const ccomplex_t IT_0062 = IT_0005*IT_0061;
    const ccomplex_t IT_0063 = IT_0018*IT_0043;
    const ccomplex_t IT_0064 = IT_0005*IT_0063;
    const ccomplex_t IT_0065 = IT_0001*IT_0018;
    const ccomplex_t IT_0066 = IT_0005*IT_0065;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = IT_0012 + IT_0067;
    const ccomplex_t IT_0069 = IT_0064*IT_0068;
    const ccomplex_t IT_0070 = (-0.666666666666667)*(IT_0005 + -IT_0008 + 
      -IT_0015 + -IT_0016 + -IT_0021 + -IT_0022 + IT_0024 + IT_0027)*IT_0034 +
       0.666666666666667*IT_0032*(IT_0032*(IT_0037 + IT_0039 + 2*IT_0042) +
       IT_0033*(IT_0050 + -IT_0058 + -IT_0062 + IT_0069));
    const ccomplex_t IT_0071 = cabs(-IT_0002 + IT_0001*IT_0017);
    const ccomplex_t IT_0072 = cpow(IT_0071, 0.5);
    const ccomplex_t IT_0073 = ( ((m_W
) >= (m_W + m_Z
)) ? (0.125*IT_0000
      *IT_0001*IT_0070*IT_0072
) : (0
));
    const ccomplex_t IT_0074 = pow(m_c, 2);
    const ccomplex_t IT_0075 = pow(m_b, 2);
    const ccomplex_t IT_0076 = cpow(IT_0002 + IT_0074 + -IT_0075, 2);
    const ccomplex_t IT_0077 = -IT_0074 + 0.25*IT_0001*(IT_0002 + -IT_0074 +
       IT_0075)*(IT_0002 + IT_0074 + -IT_0075) + 0.25*IT_0001*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0002 + 0.5*IT_0074 + (-0.5)*IT_0075;
    const ccomplex_t IT_0079 = 0.5*IT_0002 + (-0.5)*IT_0074 + 0.5*IT_0075;
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = IT_0002*IT_0077;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0080 + IT_0082;
    const ccomplex_t IT_0084 = IT_0001*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0030;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = cpow(IT_0086, 2);
    const ccomplex_t IT_0088 = cabs(-IT_0074 + 0.25*IT_0001*IT_0076);
    const ccomplex_t IT_0089 = cpow(IT_0088, 0.5);
    const ccomplex_t IT_0090 = ( ((m_W
) >= (m_b + m_c
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0077 + IT_0084)*IT_0087*IT_0089
) : (0
));
    const ccomplex_t IT_0091 = pow(m_d, 2);
    const ccomplex_t IT_0092 = cpow(IT_0002 + IT_0074 + -IT_0091, 2);
    const ccomplex_t IT_0093 = -IT_0074 + 0.25*IT_0001*(IT_0002 + -IT_0074 +
       IT_0091)*(IT_0002 + IT_0074 + -IT_0091) + 0.25*IT_0001*IT_0092;
    const ccomplex_t IT_0094 = 0.5*IT_0002 + 0.5*IT_0074 + (-0.5)*IT_0091;
    const ccomplex_t IT_0095 = 0.5*IT_0002 + (-0.5)*IT_0074 + 0.5*IT_0091;
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = IT_0002*IT_0093;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0096 + IT_0098;
    const ccomplex_t IT_0100 = IT_0001*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*V_cd*e_em
      *IT_0030;
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = cabs(-IT_0074 + 0.25*IT_0001*IT_0092);
    const ccomplex_t IT_0104 = cpow(IT_0103, 0.5);
    const ccomplex_t IT_0105 = ( ((m_W
) >= (m_c + m_d
)) ? (0.5*IT_0000
      *IT_0001*(IT_0093 + IT_0100)*IT_0102*conj(IT_0102)*IT_0104
) : (0
));
    const ccomplex_t IT_0106 = pow(m_e, 2);
    const ccomplex_t IT_0107 = cpow(IT_0002 + -IT_0106, 2);
    const ccomplex_t IT_0108 = 0.25*IT_0001*(IT_0002 + IT_0106)*(IT_0002 + 
      -IT_0106) + 0.25*IT_0001*IT_0107;
    const ccomplex_t IT_0109 = 0.5*IT_0002 + (-0.5)*IT_0106;
    const ccomplex_t IT_0110 = 0.5*IT_0002 + 0.5*IT_0106;
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = IT_0002*IT_0108;
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0111 + IT_0113;
    const ccomplex_t IT_0115 = IT_0001*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0030;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = cpow(IT_0117, 2);
    const ccomplex_t IT_0119 = cabs(0.25*IT_0001*IT_0107);
    const ccomplex_t IT_0120 = cpow(IT_0119, 0.5);
    const ccomplex_t IT_0121 = ( ((m_W
) >= (m_e
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0108 + IT_0115)*IT_0118*IT_0120
) : (0
));
    const ccomplex_t IT_0122 = pow(m_h, 2);
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0002 + IT_0123;
    const ccomplex_t IT_0125 = cpow(IT_0124, 2);
    const ccomplex_t IT_0126 = IT_0006*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*m_W*e_em*IT_0030;
    const ccomplex_t IT_0128 = cpow(IT_0127, 2);
    const ccomplex_t IT_0129 = cpow(IT_0002 + (-0.5)*IT_0122, 2);
    const ccomplex_t IT_0130 = cabs(-IT_0002 + IT_0001*IT_0129);
    const ccomplex_t IT_0131 = cpow(IT_0130, 0.5);
    const ccomplex_t IT_0132 = ( ((m_W
) >= (m_W + m_h
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*(1 + 0.5*IT_0126)*IT_0128*IT_0131
) :
       (0
));
    const ccomplex_t IT_0133 = pow(m_s, 2);
    const ccomplex_t IT_0134 = cpow(IT_0002 + IT_0074 + -IT_0133, 2);
    const ccomplex_t IT_0135 = -IT_0074 + 0.25*IT_0001*(IT_0002 + -IT_0074 +
       IT_0133)*(IT_0002 + IT_0074 + -IT_0133) + 0.25*IT_0001*IT_0134;
    const ccomplex_t IT_0136 = 0.5*IT_0002 + 0.5*IT_0074 + (-0.5)*IT_0133;
    const ccomplex_t IT_0137 = 0.5*IT_0002 + (-0.5)*IT_0074 + 0.5*IT_0133;
    const ccomplex_t IT_0138 = IT_0136*IT_0137;
    const ccomplex_t IT_0139 = IT_0002*IT_0135;
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = IT_0138 + IT_0140;
    const ccomplex_t IT_0142 = IT_0001*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0030;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = cabs(-IT_0074 + 0.25*IT_0001*IT_0134);
    const ccomplex_t IT_0146 = cpow(IT_0145, 0.5);
    const ccomplex_t IT_0147 = ( ((m_W
) >= (m_c + m_s
)) ? (0.5*IT_0000
      *IT_0001*(IT_0135 + IT_0142)*IT_0144*conj(IT_0144)*IT_0146
) : (0
));
    const ccomplex_t IT_0148 = pow(m_t, 2);
    const ccomplex_t IT_0149 = cpow(IT_0002 + -IT_0075 + IT_0148, 2);
    const ccomplex_t IT_0150 = -IT_0148 + 0.25*IT_0001*(IT_0002 + -IT_0075 +
       IT_0148)*(IT_0002 + IT_0075 + -IT_0148) + 0.25*IT_0001*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0002 + (-0.5)*IT_0075 + 0.5*IT_0148;
    const ccomplex_t IT_0152 = 0.5*IT_0002 + 0.5*IT_0075 + (-0.5)*IT_0148;
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = IT_0002*IT_0150;
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = IT_0001*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0030;
    const ccomplex_t IT_0159 = 0.5*IT_0158;
    const ccomplex_t IT_0160 = cpow(IT_0159, 2);
    const ccomplex_t IT_0161 = cabs(-IT_0148 + 0.25*IT_0001*IT_0149);
    const ccomplex_t IT_0162 = cpow(IT_0161, 0.5);
    const ccomplex_t IT_0163 = ( ((m_W
) >= (m_b + m_t
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0150 + IT_0157)*IT_0160*IT_0162
) : (0
));
    const ccomplex_t IT_0164 = cpow(IT_0002 + -IT_0091 + IT_0148, 2);
    const ccomplex_t IT_0165 = -IT_0148 + 0.25*IT_0001*(IT_0002 + -IT_0091 +
       IT_0148)*(IT_0002 + IT_0091 + -IT_0148) + 0.25*IT_0001*IT_0164;
    const ccomplex_t IT_0166 = 0.5*IT_0002 + (-0.5)*IT_0091 + 0.5*IT_0148;
    const ccomplex_t IT_0167 = 0.5*IT_0002 + 0.5*IT_0091 + (-0.5)*IT_0148;
    const ccomplex_t IT_0168 = IT_0166*IT_0167;
    const ccomplex_t IT_0169 = IT_0002*IT_0165;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = IT_0168 + IT_0170;
    const ccomplex_t IT_0172 = IT_0001*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1.4142135623731)*V_td*e_em
      *IT_0030;
    const ccomplex_t IT_0174 = 0.5*IT_0173;
    const ccomplex_t IT_0175 = cabs(-IT_0148 + 0.25*IT_0001*IT_0164);
    const ccomplex_t IT_0176 = cpow(IT_0175, 0.5);
    const ccomplex_t IT_0177 = ( ((m_W
) >= (m_d + m_t
)) ? (0.5*IT_0000
      *IT_0001*(IT_0165 + IT_0172)*IT_0174*conj(IT_0174)*IT_0176
) : (0
));
    const ccomplex_t IT_0178 = cpow(IT_0002 + -IT_0133 + IT_0148, 2);
    const ccomplex_t IT_0179 = -IT_0148 + 0.25*IT_0001*(IT_0002 + -IT_0133 +
       IT_0148)*(IT_0002 + IT_0133 + -IT_0148) + 0.25*IT_0001*IT_0178;
    const ccomplex_t IT_0180 = 0.5*IT_0002 + (-0.5)*IT_0133 + 0.5*IT_0148;
    const ccomplex_t IT_0181 = 0.5*IT_0002 + 0.5*IT_0133 + (-0.5)*IT_0148;
    const ccomplex_t IT_0182 = IT_0180*IT_0181;
    const ccomplex_t IT_0183 = IT_0002*IT_0179;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = IT_0182 + IT_0184;
    const ccomplex_t IT_0186 = IT_0001*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1.4142135623731)*V_ts*e_em
      *IT_0030;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = cabs(-IT_0148 + 0.25*IT_0001*IT_0178);
    const ccomplex_t IT_0190 = cpow(IT_0189, 0.5);
    const ccomplex_t IT_0191 = ( ((m_W
) >= (m_s + m_t
)) ? (0.5*IT_0000
      *IT_0001*(IT_0179 + IT_0186)*IT_0188*conj(IT_0188)*IT_0190
) : (0
));
    const ccomplex_t IT_0192 = pow(m_u, 2);
    const ccomplex_t IT_0193 = cpow(IT_0002 + -IT_0075 + IT_0192, 2);
    const ccomplex_t IT_0194 = -IT_0192 + 0.25*IT_0001*(IT_0002 + -IT_0075 +
       IT_0192)*(IT_0002 + IT_0075 + -IT_0192) + 0.25*IT_0001*IT_0193;
    const ccomplex_t IT_0195 = 0.5*IT_0002 + (-0.5)*IT_0075 + 0.5*IT_0192;
    const ccomplex_t IT_0196 = 0.5*IT_0002 + 0.5*IT_0075 + (-0.5)*IT_0192;
    const ccomplex_t IT_0197 = IT_0195*IT_0196;
    const ccomplex_t IT_0198 = IT_0002*IT_0194;
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0197 + IT_0199;
    const ccomplex_t IT_0201 = IT_0001*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1.4142135623731)*V_ub*e_em
      *IT_0030;
    const ccomplex_t IT_0203 = 0.5*IT_0202;
    const ccomplex_t IT_0204 = cabs(-IT_0192 + 0.25*IT_0001*IT_0193);
    const ccomplex_t IT_0205 = cpow(IT_0204, 0.5);
    const ccomplex_t IT_0206 = ( ((m_W
) >= (m_b + m_u
)) ? (0.5*IT_0000
      *IT_0001*(IT_0194 + IT_0201)*IT_0203*conj(IT_0203)*IT_0205
) : (0
));
    const ccomplex_t IT_0207 = cpow(IT_0002 + -IT_0091 + IT_0192, 2);
    const ccomplex_t IT_0208 = -IT_0192 + 0.25*IT_0001*(IT_0002 + -IT_0091 +
       IT_0192)*(IT_0002 + IT_0091 + -IT_0192) + 0.25*IT_0001*IT_0207;
    const ccomplex_t IT_0209 = 0.5*IT_0002 + (-0.5)*IT_0091 + 0.5*IT_0192;
    const ccomplex_t IT_0210 = 0.5*IT_0002 + 0.5*IT_0091 + (-0.5)*IT_0192;
    const ccomplex_t IT_0211 = IT_0209*IT_0210;
    const ccomplex_t IT_0212 = IT_0002*IT_0208;
    const ccomplex_t IT_0213 = (-0.5)*IT_0212;
    const ccomplex_t IT_0214 = IT_0211 + IT_0213;
    const ccomplex_t IT_0215 = IT_0001*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0030;
    const ccomplex_t IT_0217 = 0.5*IT_0216;
    const ccomplex_t IT_0218 = cpow(IT_0217, 2);
    const ccomplex_t IT_0219 = cabs(-IT_0192 + 0.25*IT_0001*IT_0207);
    const ccomplex_t IT_0220 = cpow(IT_0219, 0.5);
    const ccomplex_t IT_0221 = ( ((m_W
) >= (m_d + m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0208 + IT_0215)*IT_0218*IT_0220
) : (0
));
    const ccomplex_t IT_0222 = cpow(IT_0002 + -IT_0133 + IT_0192, 2);
    const ccomplex_t IT_0223 = -IT_0192 + 0.25*IT_0001*(IT_0002 + -IT_0133 +
       IT_0192)*(IT_0002 + IT_0133 + -IT_0192) + 0.25*IT_0001*IT_0222;
    const ccomplex_t IT_0224 = 0.5*IT_0002 + (-0.5)*IT_0133 + 0.5*IT_0192;
    const ccomplex_t IT_0225 = 0.5*IT_0002 + 0.5*IT_0133 + (-0.5)*IT_0192;
    const ccomplex_t IT_0226 = IT_0224*IT_0225;
    const ccomplex_t IT_0227 = IT_0002*IT_0223;
    const ccomplex_t IT_0228 = (-0.5)*IT_0227;
    const ccomplex_t IT_0229 = IT_0226 + IT_0228;
    const ccomplex_t IT_0230 = IT_0001*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0030;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = cpow(IT_0232, 2);
    const ccomplex_t IT_0234 = cabs(-IT_0192 + 0.25*IT_0001*IT_0222);
    const ccomplex_t IT_0235 = cpow(IT_0234, 0.5);
    const ccomplex_t IT_0236 = ( ((m_W
) >= (m_s + m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0223 + IT_0230)*IT_0233*IT_0235
) : (0
));
    const ccomplex_t IT_0237 = pow(m_mu, 2);
    const ccomplex_t IT_0238 = cpow(IT_0002 + -IT_0237, 2);
    const ccomplex_t IT_0239 = 0.25*IT_0001*(IT_0002 + IT_0237)*(IT_0002 + 
      -IT_0237) + 0.25*IT_0001*IT_0238;
    const ccomplex_t IT_0240 = 0.5*IT_0002 + (-0.5)*IT_0237;
    const ccomplex_t IT_0241 = 0.5*IT_0002 + 0.5*IT_0237;
    const ccomplex_t IT_0242 = IT_0240*IT_0241;
    const ccomplex_t IT_0243 = IT_0002*IT_0239;
    const ccomplex_t IT_0244 = (-0.5)*IT_0243;
    const ccomplex_t IT_0245 = IT_0242 + IT_0244;
    const ccomplex_t IT_0246 = IT_0001*IT_0245;
    const ccomplex_t IT_0247 = cabs(0.25*IT_0001*IT_0238);
    const ccomplex_t IT_0248 = cpow(IT_0247, 0.5);
    const ccomplex_t IT_0249 = ( ((m_W
) >= (m_mu
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0118*(IT_0239 + IT_0246)*IT_0248
) : (0
));
    const ccomplex_t IT_0250 = pow(m_tau, 2);
    const ccomplex_t IT_0251 = cpow(IT_0002 + -IT_0250, 2);
    const ccomplex_t IT_0252 = 0.25*IT_0001*(IT_0002 + IT_0250)*(IT_0002 + 
      -IT_0250) + 0.25*IT_0001*IT_0251;
    const ccomplex_t IT_0253 = 0.5*IT_0002 + (-0.5)*IT_0250;
    const ccomplex_t IT_0254 = 0.5*IT_0002 + 0.5*IT_0250;
    const ccomplex_t IT_0255 = IT_0253*IT_0254;
    const ccomplex_t IT_0256 = IT_0002*IT_0252;
    const ccomplex_t IT_0257 = (-0.5)*IT_0256;
    const ccomplex_t IT_0258 = IT_0255 + IT_0257;
    const ccomplex_t IT_0259 = IT_0001*IT_0258;
    const ccomplex_t IT_0260 = cabs(0.25*IT_0001*IT_0251);
    const ccomplex_t IT_0261 = cpow(IT_0260, 0.5);
    const ccomplex_t IT_0262 = ( ((m_W
) >= (m_tau
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0118*(IT_0252 + IT_0259)*IT_0261
) : (0
));
    return create_ccomplex_return(IT_0073 + IT_0090 + IT_0105 + IT_0121 +
       IT_0132 + IT_0147 + IT_0163 + IT_0177 + IT_0191 + IT_0206 + IT_0221 +
       IT_0236 + IT_0249 + IT_0262);
}

