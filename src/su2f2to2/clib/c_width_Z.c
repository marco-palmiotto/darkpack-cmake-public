#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_Z(
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
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = e_em*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0004, -1);
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0011, 2);
    const ccomplex_t IT_0013 = pow(m_Z, 2);
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    const ccomplex_t IT_0016 = IT_0013*IT_0014;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0015 + IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = cabs(0.25*IT_0013);
    const ccomplex_t IT_0021 = cpow(IT_0020, 0.5);
    const ccomplex_t IT_0022 = ( ((m_Z
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0012*(IT_0014 + IT_0019)*IT_0021
) : (0
));
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*IT_0002*IT_0007;
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -2);
    const ccomplex_t IT_0026 = IT_0001*IT_0015*IT_0025;
    const ccomplex_t IT_0027 = IT_0015*IT_0025;
    const ccomplex_t IT_0028 = (1 + 0.5*IT_0026)*(IT_0013 + -IT_0027);
    const ccomplex_t IT_0029 = pow(m_W, 2);
    const ccomplex_t IT_0030 = 0.5*IT_0013 + -IT_0029;
    const ccomplex_t IT_0031 = cpow(IT_0030, 2);
    const ccomplex_t IT_0032 = IT_0025*IT_0031;
    const ccomplex_t IT_0033 = (1 + 0.5*IT_0026)*(IT_0029 + -IT_0032);
    const ccomplex_t IT_0034 = IT_0025*IT_0030;
    const ccomplex_t IT_0035 = IT_0014*IT_0034;
    const ccomplex_t IT_0036 = (1 + 0.5*IT_0026)*(IT_0014 + -IT_0035);
    const ccomplex_t IT_0037 = 2*IT_0023;
    const ccomplex_t IT_0038 = pow(m_W, -4);
    const ccomplex_t IT_0039 = IT_0014*IT_0038;
    const ccomplex_t IT_0040 = IT_0014*IT_0039;
    const ccomplex_t IT_0041 = IT_0001*IT_0014;
    const ccomplex_t IT_0042 = IT_0014*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0030 + IT_0043;
    const ccomplex_t IT_0045 = IT_0040*IT_0044;
    const ccomplex_t IT_0046 = (-2)*IT_0023;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = IT_0001*IT_0025;
    const ccomplex_t IT_0049 = IT_0030*IT_0048;
    const ccomplex_t IT_0050 = IT_0014*IT_0049;
    const ccomplex_t IT_0051 = -IT_0035;
    const ccomplex_t IT_0052 = IT_0014 + IT_0051;
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = IT_0001*IT_0015;
    const ccomplex_t IT_0055 = IT_0031*IT_0038;
    const ccomplex_t IT_0056 = (IT_0029 + -IT_0054)*(1 + 0.5*IT_0055);
    const ccomplex_t IT_0057 = (-2)*IT_0026;
    const ccomplex_t IT_0058 = IT_0001*IT_0038;
    const ccomplex_t IT_0059 = IT_0030*IT_0058;
    const ccomplex_t IT_0060 = IT_0014*IT_0059;
    const ccomplex_t IT_0061 = IT_0014*IT_0060;
    const ccomplex_t IT_0062 = IT_0057 + IT_0061;
    const ccomplex_t IT_0063 = 1 + IT_0062;
    const ccomplex_t IT_0064 = IT_0014*IT_0063;
    const ccomplex_t IT_0065 = (-0.666666666666667)*IT_0024;
    const ccomplex_t IT_0066 = cpow(IT_0014, 3);
    const ccomplex_t IT_0067 = IT_0001*IT_0025*IT_0066;
    const ccomplex_t IT_0068 = IT_0014*IT_0026;
    const ccomplex_t IT_0069 = IT_0014*IT_0055;
    const ccomplex_t IT_0070 = IT_0001*IT_0038*IT_0066;
    const ccomplex_t IT_0071 = IT_0030*IT_0070;
    const ccomplex_t IT_0072 = -IT_0026 + -IT_0055;
    const ccomplex_t IT_0073 = IT_0061 + IT_0072;
    const ccomplex_t IT_0074 = 1 + IT_0073;
    const ccomplex_t IT_0075 = IT_0014*IT_0074;
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0024*(IT_0024*(IT_0028 +
       IT_0033 + 2*IT_0036) + IT_0037*IT_0045 + IT_0047*IT_0053) +
       0.666666666666667*IT_0047*(IT_0047*IT_0056 + 1.5*IT_0064*IT_0065) +
       0.666666666666667*IT_0037*(IT_0028*IT_0037 + -IT_0047*(IT_0014 + IT_0035 
      + -IT_0067 + -IT_0068 + -IT_0069 + IT_0071) + 1.5*IT_0065*IT_0075);
    const ccomplex_t IT_0077 = cabs(0.25*IT_0013 + -IT_0029);
    const ccomplex_t IT_0078 = cpow(IT_0077, 0.5);
    const ccomplex_t IT_0079 = ( ((m_Z
) >= (2*m_W
)) ? (0.125*IT_0000*IT_0001
      *IT_0076*IT_0078
) : (0
));
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0081 = 0.333333333333333*IT_0080;
    const ccomplex_t IT_0082 = cpow(IT_0081, 2);
    const ccomplex_t IT_0083 = pow(m_b, 2);
    const ccomplex_t IT_0084 = 0.5*IT_0013 + -IT_0083;
    const ccomplex_t IT_0085 = IT_0013*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0015 + IT_0086;
    const ccomplex_t IT_0088 = IT_0001*IT_0087;
    const ccomplex_t IT_0089 = IT_0084 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0006 + 3*IT_0009);
    const ccomplex_t IT_0091 = (-0.166666666666667)*IT_0090;
    const ccomplex_t IT_0092 = cabs(0.25*IT_0013 + -IT_0083);
    const ccomplex_t IT_0093 = cpow(IT_0092, 0.5);
    const ccomplex_t IT_0094 = ( ((m_Z
) >= (2*m_b
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0082*IT_0089 + (-0.25)*IT_0091*((-12)*IT_0081*IT_0083 + (-4)
      *IT_0089*IT_0091))*IT_0093
) : (0
));
    const ccomplex_t IT_0095 = (-0.666666666666667)*IT_0080;
    const ccomplex_t IT_0096 = cpow(IT_0095, 2);
    const ccomplex_t IT_0097 = pow(m_c, 2);
    const ccomplex_t IT_0098 = 0.5*IT_0013 + -IT_0097;
    const ccomplex_t IT_0099 = IT_0013*IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0015 + IT_0100;
    const ccomplex_t IT_0102 = IT_0001*IT_0101;
    const ccomplex_t IT_0103 = IT_0098 + IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0006 + (-3)*IT_0009);
    const ccomplex_t IT_0105 = (-0.166666666666667)*IT_0104;
    const ccomplex_t IT_0106 = cabs(0.25*IT_0013 + -IT_0097);
    const ccomplex_t IT_0107 = cpow(IT_0106, 0.5);
    const ccomplex_t IT_0108 = ( ((m_Z
) >= (2*m_c
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0096*IT_0103 + (-0.25)*IT_0105*((-12)*IT_0095*IT_0097 + (-4)
      *IT_0103*IT_0105))*IT_0107
) : (0
));
    const ccomplex_t IT_0109 = pow(m_d, 2);
    const ccomplex_t IT_0110 = 0.5*IT_0013 + -IT_0109;
    const ccomplex_t IT_0111 = IT_0013*IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = IT_0015 + IT_0112;
    const ccomplex_t IT_0114 = IT_0001*IT_0113;
    const ccomplex_t IT_0115 = IT_0110 + IT_0114;
    const ccomplex_t IT_0116 = cabs(0.25*IT_0013 + -IT_0109);
    const ccomplex_t IT_0117 = cpow(IT_0116, 0.5);
    const ccomplex_t IT_0118 = ( ((m_Z
) >= (2*m_d
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0082*IT_0115 + (-0.25)*IT_0091*((-12)*IT_0081*IT_0109 + (-4)
      *IT_0091*IT_0115))*IT_0117
) : (0
));
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0006 + -IT_0009);
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = cpow(IT_0120, 2);
    const ccomplex_t IT_0122 = pow(m_e, 2);
    const ccomplex_t IT_0123 = 0.5*IT_0013 + -IT_0122;
    const ccomplex_t IT_0124 = IT_0013*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0015 + IT_0125;
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = IT_0123 + IT_0127;
    const ccomplex_t IT_0129 = cabs(0.25*IT_0013 + -IT_0122);
    const ccomplex_t IT_0130 = cpow(IT_0129, 0.5);
    const ccomplex_t IT_0131 = ( ((m_Z
) >= (2*m_e
)) ? (0.125*IT_0000*IT_0001
      *((-1.33333333333333)*IT_0121*IT_0128 + (-4)*IT_0080*(IT_0120*IT_0122 +
       0.333333333333333*IT_0080*IT_0128))*IT_0130
) : (0
));
    const ccomplex_t IT_0132 = m_W*e_em;
    const ccomplex_t IT_0133 = cpow(IT_0002, -2);
    const ccomplex_t IT_0134 = cpow(IT_0004, 3);
    const ccomplex_t IT_0135 = IT_0133*IT_0134;
    const ccomplex_t IT_0136 = cpow(IT_0002, 2);
    const ccomplex_t IT_0137 = IT_0007*IT_0136;
    const ccomplex_t IT_0138 = IT_0135 + IT_0137;
    const ccomplex_t IT_0139 = 2*IT_0004;
    const ccomplex_t IT_0140 = IT_0138 + IT_0139;
    const ccomplex_t IT_0141 = IT_0132*IT_0140;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = 2*IT_0143;
    const ccomplex_t IT_0145 = cpow(IT_0144, 2);
    const ccomplex_t IT_0146 = pow(m_Z, -4);
    const ccomplex_t IT_0147 = pow(m_h, 2);
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0013 + IT_0148;
    const ccomplex_t IT_0150 = cpow(IT_0149, 2);
    const ccomplex_t IT_0151 = IT_0146*IT_0150;
    const ccomplex_t IT_0152 = pow(m_h, 4);
    const ccomplex_t IT_0153 = cabs(-IT_0147 + 0.25*IT_0001*IT_0152);
    const ccomplex_t IT_0154 = cpow(IT_0153, 0.5);
    const ccomplex_t IT_0155 = ( ((m_Z
) >= (m_Z + m_h
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0145*(1 + 0.5*IT_0151)*IT_0154
) :
       (0
));
    const ccomplex_t IT_0156 = pow(m_s, 2);
    const ccomplex_t IT_0157 = 0.5*IT_0013 + -IT_0156;
    const ccomplex_t IT_0158 = IT_0013*IT_0157;
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0015 + IT_0159;
    const ccomplex_t IT_0161 = IT_0001*IT_0160;
    const ccomplex_t IT_0162 = IT_0157 + IT_0161;
    const ccomplex_t IT_0163 = cabs(0.25*IT_0013 + -IT_0156);
    const ccomplex_t IT_0164 = cpow(IT_0163, 0.5);
    const ccomplex_t IT_0165 = ( ((m_Z
) >= (2*m_s
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0082*IT_0162 + (-0.25)*IT_0091*((-12)*IT_0081*IT_0156 + (-4)
      *IT_0091*IT_0162))*IT_0164
) : (0
));
    const ccomplex_t IT_0166 = pow(m_t, 2);
    const ccomplex_t IT_0167 = 0.5*IT_0013 + -IT_0166;
    const ccomplex_t IT_0168 = IT_0013*IT_0167;
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0015 + IT_0169;
    const ccomplex_t IT_0171 = IT_0001*IT_0170;
    const ccomplex_t IT_0172 = IT_0167 + IT_0171;
    const ccomplex_t IT_0173 = cabs(0.25*IT_0013 + -IT_0166);
    const ccomplex_t IT_0174 = cpow(IT_0173, 0.5);
    const ccomplex_t IT_0175 = ( ((m_Z
) >= (2*m_t
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0096*IT_0172 + (-0.25)*IT_0105*((-12)*IT_0095*IT_0166 + (-4)
      *IT_0105*IT_0172))*IT_0174
) : (0
));
    const ccomplex_t IT_0176 = pow(m_u, 2);
    const ccomplex_t IT_0177 = 0.5*IT_0013 + -IT_0176;
    const ccomplex_t IT_0178 = IT_0013*IT_0177;
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = IT_0015 + IT_0179;
    const ccomplex_t IT_0181 = IT_0001*IT_0180;
    const ccomplex_t IT_0182 = IT_0177 + IT_0181;
    const ccomplex_t IT_0183 = cabs(0.25*IT_0013 + -IT_0176);
    const ccomplex_t IT_0184 = cpow(IT_0183, 0.5);
    const ccomplex_t IT_0185 = ( ((m_Z
) >= (2*m_u
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0096*IT_0182 + (-0.25)*IT_0105*((-12)*IT_0095*IT_0176 + (-4)
      *IT_0105*IT_0182))*IT_0184
) : (0
));
    const ccomplex_t IT_0186 = pow(m_mu, 2);
    const ccomplex_t IT_0187 = 0.5*IT_0013 + -IT_0186;
    const ccomplex_t IT_0188 = IT_0013*IT_0187;
    const ccomplex_t IT_0189 = (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = IT_0015 + IT_0189;
    const ccomplex_t IT_0191 = IT_0001*IT_0190;
    const ccomplex_t IT_0192 = IT_0187 + IT_0191;
    const ccomplex_t IT_0193 = cabs(0.25*IT_0013 + -IT_0186);
    const ccomplex_t IT_0194 = cpow(IT_0193, 0.5);
    const ccomplex_t IT_0195 = ( ((m_Z
) >= (2*m_mu
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0121*IT_0192 + (-4)*IT_0080*(IT_0120
      *IT_0186 + 0.333333333333333*IT_0080*IT_0192))*IT_0194
) : (0
));
    const ccomplex_t IT_0196 = pow(m_tau, 2);
    const ccomplex_t IT_0197 = 0.5*IT_0013 + -IT_0196;
    const ccomplex_t IT_0198 = IT_0013*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0015 + IT_0199;
    const ccomplex_t IT_0201 = IT_0001*IT_0200;
    const ccomplex_t IT_0202 = IT_0197 + IT_0201;
    const ccomplex_t IT_0203 = cabs(0.25*IT_0013 + -IT_0196);
    const ccomplex_t IT_0204 = cpow(IT_0203, 0.5);
    const ccomplex_t IT_0205 = ( ((m_Z
) >= (2*m_tau
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0121*IT_0202 + (-4)*IT_0080*(IT_0120
      *IT_0196 + 0.333333333333333*IT_0080*IT_0202))*IT_0204
) : (0
));
    return create_ccomplex_return(3*IT_0022 + IT_0079 + IT_0094 + IT_0108 +
       IT_0118 + IT_0131 + IT_0155 + IT_0165 + IT_0175 + IT_0185 + IT_0195 +
       IT_0205);
}

