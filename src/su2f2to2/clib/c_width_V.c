#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_V.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_V(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t g_f = param->g_f;
    const creal_t m_V = param->m_V;
    const creal_t m_e = param->m_e;
    const creal_t m_mu = param->m_mu;
    const creal_t m_phi = param->m_phi;
    const creal_t m_tau = param->m_tau;
    const creal_t m_chi_dm_1 = param->m_chi_dm_1;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const creal_t theta_12_L = param->theta_12_L;
    const creal_t theta_12_R = param->theta_12_R;
    const creal_t theta_13_L = param->theta_13_L;
    const creal_t theta_13_R = param->theta_13_R;
    const creal_t theta_23_L = param->theta_23_L;
    const creal_t theta_23_R = param->theta_23_R;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_V, -2);
    const ccomplex_t IT_0002 = pow(m_V, 2);
    const ccomplex_t IT_0003 = 0.5*IT_0002;
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = IT_0002*IT_0003;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = IT_0004 + IT_0006;
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = cos(theta_23_L);
    const ccomplex_t IT_0010 = sin(theta_13_L);
    const ccomplex_t IT_0011 = cos(theta_13_L);
    const ccomplex_t IT_0012 = sin(theta_12_L);
    const ccomplex_t IT_0013 = cos(theta_12_L);
    const ccomplex_t IT_0014 = sin(theta_23_L);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0009
      *IT_0010*(IT_0011*IT_0012 + -IT_0010*IT_0013*IT_0014);
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0016, 2);
    const ccomplex_t IT_0018 = cabs(0.25*IT_0002);
    const ccomplex_t IT_0019 = cpow(IT_0018, 0.5);
    const ccomplex_t IT_0020 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0017*IT_0019
) : (0
));
    const ccomplex_t IT_0021 = cpow(IT_0009, 2);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0010
      *IT_0013*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0023, 2);
    const ccomplex_t IT_0025 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0019*IT_0024
) : (0
));
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0014*
      (IT_0011*IT_0012 + -IT_0010*IT_0013*IT_0014);
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0027, 2);
    const ccomplex_t IT_0029 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0019*IT_0028
) : (0
));
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0009
      *IT_0010*(IT_0010*IT_0012 + IT_0011*IT_0013*IT_0014);
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = cpow(IT_0031, 2);
    const ccomplex_t IT_0033 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0019*IT_0032
) : (0
));
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0009
      *IT_0011*(IT_0011*IT_0012 + -IT_0010*IT_0013*IT_0014);
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0035, 2);
    const ccomplex_t IT_0037 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0019*IT_0036
) : (0
));
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0009
      *IT_0013*IT_0014;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0039, 2);
    const ccomplex_t IT_0041 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0019*IT_0040
) : (0
));
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0014*
      (IT_0010*IT_0012 + IT_0011*IT_0013*IT_0014);
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = cpow(IT_0043, 2);
    const ccomplex_t IT_0045 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0019*IT_0044
) : (0
));
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0011
      *IT_0013*IT_0021;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0047, 2);
    const ccomplex_t IT_0049 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0019*IT_0048
) : (0
));
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0009
      *IT_0011*(IT_0010*IT_0012 + IT_0011*IT_0013*IT_0014);
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = cpow(IT_0051, 2);
    const ccomplex_t IT_0053 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0019*IT_0052
) : (0
));
    const ccomplex_t IT_0054 = IT_0001*IT_0004;
    const ccomplex_t IT_0055 = pow(m_V, -4);
    const ccomplex_t IT_0056 = IT_0004*IT_0055;
    const ccomplex_t IT_0057 = (IT_0002 + -IT_0054)*(1 + 0.5*IT_0056);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = -IT_0058 + IT_0059;
    const ccomplex_t IT_0061 = cpow(IT_0060, 2);
    const ccomplex_t IT_0062 = (-2)*IT_0058;
    const ccomplex_t IT_0063 = (-0.5)*IT_0002;
    const ccomplex_t IT_0064 = cpow(IT_0063, 2);
    const ccomplex_t IT_0065 = IT_0055*IT_0064;
    const ccomplex_t IT_0066 = (IT_0002 + -IT_0054)*(1 + 0.5*IT_0065);
    const ccomplex_t IT_0067 = (-2)*IT_0056;
    const ccomplex_t IT_0068 = pow(m_V, -6);
    const ccomplex_t IT_0069 = IT_0063*IT_0068;
    const ccomplex_t IT_0070 = IT_0003*IT_0069;
    const ccomplex_t IT_0071 = IT_0003*IT_0070;
    const ccomplex_t IT_0072 = IT_0067 + IT_0071;
    const ccomplex_t IT_0073 = 1 + IT_0072;
    const ccomplex_t IT_0074 = IT_0003*IT_0073;
    const ccomplex_t IT_0075 = 2*IT_0058;
    const ccomplex_t IT_0076 = -IT_0058 + IT_0075;
    const ccomplex_t IT_0077 = cpow(IT_0003, 3);
    const ccomplex_t IT_0078 = IT_0055*IT_0077;
    const ccomplex_t IT_0079 = IT_0003*IT_0056;
    const ccomplex_t IT_0080 = IT_0003*IT_0065;
    const ccomplex_t IT_0081 = IT_0001*IT_0063;
    const ccomplex_t IT_0082 = IT_0003*IT_0081;
    const ccomplex_t IT_0083 = IT_0068*IT_0077;
    const ccomplex_t IT_0084 = IT_0063*IT_0083;
    const ccomplex_t IT_0085 = IT_0001*IT_0064;
    const ccomplex_t IT_0086 = (1 + 0.5*IT_0056)*(IT_0002 + -IT_0085);
    const ccomplex_t IT_0087 = IT_0055*IT_0063;
    const ccomplex_t IT_0088 = IT_0003*IT_0087;
    const ccomplex_t IT_0089 = -IT_0082;
    const ccomplex_t IT_0090 = IT_0003 + IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = (1 + 0.5*IT_0056)*(IT_0003 + -IT_0082);
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = -IT_0056 + -IT_0065;
    const ccomplex_t IT_0096 = IT_0071 + IT_0095;
    const ccomplex_t IT_0097 = 1 + IT_0096;
    const ccomplex_t IT_0098 = IT_0003*IT_0097;
    const ccomplex_t IT_0099 = IT_0003*IT_0055;
    const ccomplex_t IT_0100 = IT_0003*IT_0099;
    const ccomplex_t IT_0101 = IT_0001*IT_0003;
    const ccomplex_t IT_0102 = IT_0003*IT_0101;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0063 + IT_0103;
    const ccomplex_t IT_0105 = IT_0100*IT_0104;
    const ccomplex_t IT_0106 = 0.666666666666667*IT_0057*IT_0061 +
       0.666666666666667*IT_0062*(IT_0062*IT_0066 + -IT_0074*IT_0076 + -IT_0060*
      (IT_0003 + -IT_0078 + -IT_0079 + -IT_0080 + IT_0082 + IT_0084)) +
       0.666666666666667*IT_0058*(IT_0058*IT_0086 + IT_0062*IT_0091 + 2*IT_0076
      *IT_0094 + -IT_0060*IT_0098) + 0.666666666666667*IT_0076*(IT_0057*IT_0076 
      + IT_0060*IT_0105);
    const ccomplex_t IT_0107 = cabs((-0.75)*IT_0002);
    const ccomplex_t IT_0108 = cpow(IT_0107, 0.5);
    const ccomplex_t IT_0109 = ( ((m_V
) >= (2*m_V
)) ? (0.125*IT_0000*IT_0001
      *IT_0106*IT_0108
) : (0
));
    const ccomplex_t IT_0110 = cos(theta_23_R);
    const ccomplex_t IT_0111 = sin(theta_13_R);
    const ccomplex_t IT_0112 = cos(theta_13_R);
    const ccomplex_t IT_0113 = sin(theta_12_R);
    const ccomplex_t IT_0114 = cos(theta_12_R);
    const ccomplex_t IT_0115 = sin(theta_23_R);
    const ccomplex_t IT_0116 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0110
      *IT_0111*(IT_0112*IT_0113 + -IT_0111*IT_0114*IT_0115);
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = cpow(IT_0117, 2);
    const ccomplex_t IT_0119 = pow(m_e, 2);
    const ccomplex_t IT_0120 = 0.5*IT_0002 + -IT_0119;
    const ccomplex_t IT_0121 = IT_0002*IT_0120;
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = IT_0004 + IT_0122;
    const ccomplex_t IT_0124 = IT_0001*IT_0123;
    const ccomplex_t IT_0125 = IT_0120 + IT_0124;
    const ccomplex_t IT_0126 = (-1.33333333333333)*IT_0125;
    const ccomplex_t IT_0127 = cabs(0.25*IT_0002 + -IT_0119);
    const ccomplex_t IT_0128 = cpow(IT_0127, 0.5);
    const ccomplex_t IT_0129 = ( ((m_V
) >= (2*m_e
)) ? (0.125*IT_0000*IT_0001
      *(IT_0118*IT_0126 + (-4)*IT_0016*(IT_0117*IT_0119 + (-0.25)*IT_0016
      *IT_0126))*IT_0128
) : (0
));
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0110
      *IT_0114*IT_0115;
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = cpow(IT_0131, 2);
    const ccomplex_t IT_0133 = pow(m_mu, 2);
    const ccomplex_t IT_0134 = 0.5*IT_0002 + -IT_0133;
    const ccomplex_t IT_0135 = IT_0002*IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = IT_0004 + IT_0136;
    const ccomplex_t IT_0138 = IT_0001*IT_0137;
    const ccomplex_t IT_0139 = IT_0134 + IT_0138;
    const ccomplex_t IT_0140 = (-1.33333333333333)*IT_0139;
    const ccomplex_t IT_0141 = cabs(0.25*IT_0002 + -IT_0133);
    const ccomplex_t IT_0142 = cpow(IT_0141, 0.5);
    const ccomplex_t IT_0143 = ( ((m_V
) >= (2*m_mu
)) ? (0.125*IT_0000
      *IT_0001*(IT_0132*IT_0140 + (-4)*IT_0039*(IT_0131*IT_0133 + (-0.25)
      *IT_0039*IT_0140))*IT_0142
) : (0
));
    const ccomplex_t IT_0144 = cpow(IT_0110, 2);
    const ccomplex_t IT_0145 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0111
      *IT_0114*IT_0144;
    const ccomplex_t IT_0146 = (-0.5)*IT_0145;
    const ccomplex_t IT_0147 = cpow(IT_0146, 2);
    const ccomplex_t IT_0148 = cpow(IT_0002 + IT_0119 + -IT_0133, 2);
    const ccomplex_t IT_0149 = -IT_0119 + 0.25*IT_0001*(IT_0002 + -IT_0119 +
       IT_0133)*(IT_0002 + IT_0119 + -IT_0133) + 0.25*IT_0001*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0002 + 0.5*IT_0119 + (-0.5)*IT_0133;
    const ccomplex_t IT_0151 = 0.5*IT_0002 + (-0.5)*IT_0119 + 0.5*IT_0133;
    const ccomplex_t IT_0152 = IT_0150*IT_0151;
    const ccomplex_t IT_0153 = IT_0002*IT_0149;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = IT_0152 + IT_0154;
    const ccomplex_t IT_0156 = IT_0001*IT_0155;
    const ccomplex_t IT_0157 = IT_0149 + IT_0156;
    const ccomplex_t IT_0158 = (-1.33333333333333)*IT_0157;
    const ccomplex_t IT_0159 = m_e*m_mu;
    const ccomplex_t IT_0160 = cabs(-IT_0119 + 0.25*IT_0001*IT_0148);
    const ccomplex_t IT_0161 = cpow(IT_0160, 0.5);
    const ccomplex_t IT_0162 = ( ((m_V
) >= (m_e + m_mu
)) ? (0.125*IT_0000
      *IT_0001*(IT_0147*IT_0158 + (-4)*IT_0023*((-0.25)*IT_0023*IT_0158 +
       IT_0146*IT_0159))*IT_0161
) : (0
));
    const ccomplex_t IT_0163 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0115*
      (IT_0112*IT_0113 + -IT_0111*IT_0114*IT_0115);
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = cpow(IT_0164, 2);
    const ccomplex_t IT_0166 = ( ((m_V
) >= (m_e + m_mu
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0161*(IT_0027*((-0.25)*IT_0027*IT_0158 + IT_0159*IT_0164) + 
      -1./4*IT_0158*IT_0165)
) : (0
));
    const ccomplex_t IT_0167 = pow(m_phi, 2);
    const ccomplex_t IT_0168 = (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = IT_0002 + IT_0168;
    const ccomplex_t IT_0170 = cpow(IT_0169, 2);
    const ccomplex_t IT_0171 = IT_0055*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*g_f*m_V;
    const ccomplex_t IT_0173 = cpow(IT_0172, 2);
    const ccomplex_t IT_0174 = cpow(IT_0002 + (-0.5)*IT_0167, 2);
    const ccomplex_t IT_0175 = cabs(-IT_0002 + IT_0001*IT_0174);
    const ccomplex_t IT_0176 = cpow(IT_0175, 0.5);
    const ccomplex_t IT_0177 = ( ((m_V
) >= (m_V + m_phi
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*(1 + 0.5*IT_0171)*IT_0173*IT_0176
) :
       (0
));
    const ccomplex_t IT_0178 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0110
      *IT_0112*(IT_0111*IT_0113 + IT_0112*IT_0114*IT_0115);
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = cpow(IT_0179, 2);
    const ccomplex_t IT_0181 = pow(m_tau, 2);
    const ccomplex_t IT_0182 = 0.5*IT_0002 + -IT_0181;
    const ccomplex_t IT_0183 = IT_0002*IT_0182;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = IT_0004 + IT_0184;
    const ccomplex_t IT_0186 = IT_0001*IT_0185;
    const ccomplex_t IT_0187 = IT_0182 + IT_0186;
    const ccomplex_t IT_0188 = (-1.33333333333333)*IT_0187;
    const ccomplex_t IT_0189 = cabs(0.25*IT_0002 + -IT_0181);
    const ccomplex_t IT_0190 = cpow(IT_0189, 0.5);
    const ccomplex_t IT_0191 = ( ((m_V
) >= (2*m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0180*IT_0188 + (-4)*IT_0051*(IT_0179*IT_0181 + (-0.25)
      *IT_0051*IT_0188))*IT_0190
) : (0
));
    const ccomplex_t IT_0192 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0110
      *IT_0111*(IT_0111*IT_0113 + IT_0112*IT_0114*IT_0115);
    const ccomplex_t IT_0193 = 0.5*IT_0192;
    const ccomplex_t IT_0194 = cpow(IT_0193, 2);
    const ccomplex_t IT_0195 = cpow(IT_0002 + IT_0119 + -IT_0181, 2);
    const ccomplex_t IT_0196 = -IT_0119 + 0.25*IT_0001*(IT_0002 + -IT_0119 +
       IT_0181)*(IT_0002 + IT_0119 + -IT_0181) + 0.25*IT_0001*IT_0195;
    const ccomplex_t IT_0197 = 0.5*IT_0002 + 0.5*IT_0119 + (-0.5)*IT_0181;
    const ccomplex_t IT_0198 = 0.5*IT_0002 + (-0.5)*IT_0119 + 0.5*IT_0181;
    const ccomplex_t IT_0199 = IT_0197*IT_0198;
    const ccomplex_t IT_0200 = IT_0002*IT_0196;
    const ccomplex_t IT_0201 = (-0.5)*IT_0200;
    const ccomplex_t IT_0202 = IT_0199 + IT_0201;
    const ccomplex_t IT_0203 = IT_0001*IT_0202;
    const ccomplex_t IT_0204 = IT_0196 + IT_0203;
    const ccomplex_t IT_0205 = (-1.33333333333333)*IT_0204;
    const ccomplex_t IT_0206 = m_e*m_tau;
    const ccomplex_t IT_0207 = cabs(-IT_0119 + 0.25*IT_0001*IT_0195);
    const ccomplex_t IT_0208 = cpow(IT_0207, 0.5);
    const ccomplex_t IT_0209 = ( ((m_V
) >= (m_e + m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0194*IT_0205 + (-4)*IT_0031*((-0.25)*IT_0031*IT_0205 +
       IT_0193*IT_0206))*IT_0208
) : (0
));
    const ccomplex_t IT_0210 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0110
      *IT_0112*(IT_0112*IT_0113 + -IT_0111*IT_0114*IT_0115);
    const ccomplex_t IT_0211 = (-0.5)*IT_0210;
    const ccomplex_t IT_0212 = cpow(IT_0211, 2);
    const ccomplex_t IT_0213 = ( ((m_V
) >= (m_e + m_tau
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0208*(IT_0035*((-0.25)*IT_0035*IT_0205 + IT_0206*IT_0211) + 
      -1./4*IT_0205*IT_0212)
) : (0
));
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0115*
      (IT_0111*IT_0113 + IT_0112*IT_0114*IT_0115);
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = cpow(IT_0215, 2);
    const ccomplex_t IT_0217 = cpow(IT_0002 + IT_0133 + -IT_0181, 2);
    const ccomplex_t IT_0218 = -IT_0133 + 0.25*IT_0001*(IT_0002 + -IT_0133 +
       IT_0181)*(IT_0002 + IT_0133 + -IT_0181) + 0.25*IT_0001*IT_0217;
    const ccomplex_t IT_0219 = 0.5*IT_0002 + 0.5*IT_0133 + (-0.5)*IT_0181;
    const ccomplex_t IT_0220 = 0.5*IT_0002 + (-0.5)*IT_0133 + 0.5*IT_0181;
    const ccomplex_t IT_0221 = IT_0219*IT_0220;
    const ccomplex_t IT_0222 = IT_0002*IT_0218;
    const ccomplex_t IT_0223 = (-0.5)*IT_0222;
    const ccomplex_t IT_0224 = IT_0221 + IT_0223;
    const ccomplex_t IT_0225 = IT_0001*IT_0224;
    const ccomplex_t IT_0226 = IT_0218 + IT_0225;
    const ccomplex_t IT_0227 = (-1.33333333333333)*IT_0226;
    const ccomplex_t IT_0228 = m_mu*m_tau;
    const ccomplex_t IT_0229 = cabs(-IT_0133 + 0.25*IT_0001*IT_0217);
    const ccomplex_t IT_0230 = cpow(IT_0229, 0.5);
    const ccomplex_t IT_0231 = ( ((m_V
) >= (m_mu + m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0216*IT_0227 + (-4)*IT_0043*((-0.25)*IT_0043*IT_0227 +
       IT_0215*IT_0228))*IT_0230
) : (0
));
    const ccomplex_t IT_0232 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0112
      *IT_0114*IT_0144;
    const ccomplex_t IT_0233 = (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = cpow(IT_0233, 2);
    const ccomplex_t IT_0235 = ( ((m_V
) >= (m_mu + m_tau
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0230*(IT_0047*((-0.25)*IT_0047*IT_0227 + IT_0228*IT_0233) + 
      -1./4*IT_0227*IT_0234)
) : (0
));
    const ccomplex_t IT_0236 = m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0237 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0238 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0239 = cpow(IT_0002 + IT_0237 + -IT_0238, 2);
    const ccomplex_t IT_0240 = -IT_0237 + 0.25*IT_0001*(IT_0002 + -IT_0237 +
       IT_0238)*(IT_0002 + IT_0237 + -IT_0238) + 0.25*IT_0001*IT_0239;
    const ccomplex_t IT_0241 = 0.5*IT_0002 + 0.5*IT_0237 + (-0.5)*IT_0238;
    const ccomplex_t IT_0242 = 0.5*IT_0002 + (-0.5)*IT_0237 + 0.5*IT_0238;
    const ccomplex_t IT_0243 = IT_0241*IT_0242;
    const ccomplex_t IT_0244 = IT_0002*IT_0240;
    const ccomplex_t IT_0245 = (-0.5)*IT_0244;
    const ccomplex_t IT_0246 = IT_0243 + IT_0245;
    const ccomplex_t IT_0247 = IT_0001*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1.4142135623731)*g_f;
    const ccomplex_t IT_0249 = 0.5*IT_0248;
    const ccomplex_t IT_0250 = cpow(IT_0249, 2);
    const ccomplex_t IT_0251 = cabs(-IT_0237 + 0.25*IT_0001*IT_0239);
    const ccomplex_t IT_0252 = cpow(IT_0251, 0.5);
    const ccomplex_t IT_0253 = ( ((m_V
) >= (m_chi_dm_1 + m_chi_dm_2
)) ? ((
      -0.5)*IT_0000*IT_0001*(IT_0236 + 0.666666666666667*IT_0240 +
       0.666666666666667*IT_0247)*IT_0250*IT_0252
) : (0
));
    const ccomplex_t IT_0254 = IT_0020 + IT_0025 + IT_0029 + IT_0033 + IT_0037
       + IT_0041 + IT_0045 + IT_0049 + IT_0053 + IT_0109 + IT_0129 + IT_0143 +
       IT_0162 + IT_0166 + IT_0177 + IT_0191 + IT_0209 + IT_0213 + IT_0231 +
       IT_0235 + IT_0253;
    return create_ccomplex_return(IT_0254);
}

