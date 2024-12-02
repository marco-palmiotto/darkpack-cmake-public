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
    const ccomplex_t IT_0002 = cos(theta_23_L);
    const ccomplex_t IT_0003 = sin(theta_13_L);
    const ccomplex_t IT_0004 = cos(theta_13_L);
    const ccomplex_t IT_0005 = sin(theta_12_L);
    const ccomplex_t IT_0006 = cos(theta_12_L);
    const ccomplex_t IT_0007 = sin(theta_23_L);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0002
      *IT_0003*(IT_0004*IT_0005 + -IT_0003*IT_0006*IT_0007);
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = pow(m_V, 2);
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0012, 2);
    const ccomplex_t IT_0014 = IT_0011*IT_0012;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = IT_0013 + IT_0015;
    const ccomplex_t IT_0017 = IT_0001*IT_0016;
    const ccomplex_t IT_0018 = cabs(0.25*IT_0011);
    const ccomplex_t IT_0019 = cpow(IT_0018, 0.5);
    const ccomplex_t IT_0020 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0010*(IT_0012 + IT_0017)*IT_0019
) : (0
));
    const ccomplex_t IT_0021 = cpow(IT_0002, 2);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0003
      *IT_0006*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = cpow(IT_0023, 2);
    const ccomplex_t IT_0025 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0012 + IT_0017)*IT_0019*IT_0024
) : (0
));
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0007*
      (IT_0004*IT_0005 + -IT_0003*IT_0006*IT_0007);
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0027, 2);
    const ccomplex_t IT_0029 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0012 + IT_0017)*IT_0019*IT_0028
) : (0
));
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0002
      *IT_0003*(IT_0003*IT_0005 + IT_0004*IT_0006*IT_0007);
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = cpow(IT_0031, 2);
    const ccomplex_t IT_0033 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0012 + IT_0017)*IT_0019*IT_0032
) : (0
));
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0002
      *IT_0004*(IT_0004*IT_0005 + -IT_0003*IT_0006*IT_0007);
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0035, 2);
    const ccomplex_t IT_0037 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0012 + IT_0017)*IT_0019*IT_0036
) : (0
));
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0002
      *IT_0006*IT_0007;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0039, 2);
    const ccomplex_t IT_0041 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0012 + IT_0017)*IT_0019*IT_0040
) : (0
));
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0007*
      (IT_0003*IT_0005 + IT_0004*IT_0006*IT_0007);
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = cpow(IT_0043, 2);
    const ccomplex_t IT_0045 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0012 + IT_0017)*IT_0019*IT_0044
) : (0
));
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0004
      *IT_0006*IT_0021;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = cpow(IT_0047, 2);
    const ccomplex_t IT_0049 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0012 + IT_0017)*IT_0019*IT_0048
) : (0
));
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0002
      *IT_0004*(IT_0003*IT_0005 + IT_0004*IT_0006*IT_0007);
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = cpow(IT_0051, 2);
    const ccomplex_t IT_0053 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0012 + IT_0017)*IT_0019*IT_0052
) : (0
));
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0055 = cpow(IT_0054, 2);
    const ccomplex_t IT_0056 = pow(m_V, -4);
    const ccomplex_t IT_0057 = IT_0013*IT_0056;
    const ccomplex_t IT_0058 = (-0.5)*IT_0011;
    const ccomplex_t IT_0059 = cpow(IT_0058, 2);
    const ccomplex_t IT_0060 = IT_0001*IT_0059;
    const ccomplex_t IT_0061 = (1 + 0.5*IT_0057)*(IT_0011 + -IT_0060);
    const ccomplex_t IT_0062 = (-2)*IT_0054;
    const ccomplex_t IT_0063 = IT_0001*IT_0013;
    const ccomplex_t IT_0064 = IT_0056*IT_0059;
    const ccomplex_t IT_0065 = (IT_0011 + -IT_0063)*(1 + 0.5*IT_0064);
    const ccomplex_t IT_0066 = IT_0056*IT_0058;
    const ccomplex_t IT_0067 = IT_0012*IT_0066;
    const ccomplex_t IT_0068 = IT_0001*IT_0058;
    const ccomplex_t IT_0069 = IT_0012*IT_0068;
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = IT_0012 + IT_0070;
    const ccomplex_t IT_0072 = IT_0067*IT_0071;
    const ccomplex_t IT_0073 = -IT_0054;
    const ccomplex_t IT_0074 = -IT_0054 + IT_0073;
    const ccomplex_t IT_0075 = (1 + 0.5*IT_0057)*(IT_0011 + -IT_0063);
    const ccomplex_t IT_0076 = -IT_0057 + -IT_0064;
    const ccomplex_t IT_0077 = pow(m_V, -6);
    const ccomplex_t IT_0078 = IT_0058*IT_0077;
    const ccomplex_t IT_0079 = IT_0012*IT_0078;
    const ccomplex_t IT_0080 = IT_0012*IT_0079;
    const ccomplex_t IT_0081 = IT_0076 + IT_0080;
    const ccomplex_t IT_0082 = 1 + IT_0081;
    const ccomplex_t IT_0083 = IT_0012*IT_0082;
    const ccomplex_t IT_0084 = cpow(IT_0012, 3);
    const ccomplex_t IT_0085 = IT_0056*IT_0084;
    const ccomplex_t IT_0086 = IT_0012*IT_0057;
    const ccomplex_t IT_0087 = IT_0012*IT_0064;
    const ccomplex_t IT_0088 = IT_0077*IT_0084;
    const ccomplex_t IT_0089 = IT_0058*IT_0088;
    const ccomplex_t IT_0090 = (-0.666666666666667)*IT_0062;
    const ccomplex_t IT_0091 = 2*IT_0054;
    const ccomplex_t IT_0092 = -IT_0054 + IT_0091;
    const ccomplex_t IT_0093 = (1 + 0.5*IT_0057)*(IT_0012 + -IT_0069);
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0012*IT_0056;
    const ccomplex_t IT_0097 = IT_0012*IT_0096;
    const ccomplex_t IT_0098 = IT_0001*IT_0012;
    const ccomplex_t IT_0099 = IT_0012*IT_0098;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = IT_0058 + IT_0100;
    const ccomplex_t IT_0102 = IT_0097*IT_0101;
    const ccomplex_t IT_0103 = (-2)*IT_0057;
    const ccomplex_t IT_0104 = IT_0080 + IT_0103;
    const ccomplex_t IT_0105 = 1 + IT_0104;
    const ccomplex_t IT_0106 = IT_0012*IT_0105;
    const ccomplex_t IT_0107 = 0.666666666666667*IT_0055*IT_0061 +
       0.666666666666667*IT_0062*(IT_0062*IT_0065 + IT_0054*IT_0072) +
       0.666666666666667*IT_0074*(IT_0074*IT_0075 + -IT_0054*IT_0083 + 1.5*
      (IT_0012 + IT_0069 + -IT_0085 + -IT_0086 + -IT_0087 + IT_0089)*IT_0090) +
       0.666666666666667*IT_0092*(IT_0075*IT_0092 + 2*IT_0054*IT_0095 + IT_0074
      *IT_0102 + 1.5*IT_0090*IT_0106);
    const ccomplex_t IT_0108 = cabs((-0.75)*IT_0011);
    const ccomplex_t IT_0109 = cpow(IT_0108, 0.5);
    const ccomplex_t IT_0110 = ( ((m_V
) >= (2*m_V
)) ? (0.125*IT_0000*IT_0001
      *IT_0107*IT_0109
) : (0
));
    const ccomplex_t IT_0111 = cos(theta_23_R);
    const ccomplex_t IT_0112 = sin(theta_13_R);
    const ccomplex_t IT_0113 = cos(theta_13_R);
    const ccomplex_t IT_0114 = sin(theta_12_R);
    const ccomplex_t IT_0115 = cos(theta_12_R);
    const ccomplex_t IT_0116 = sin(theta_23_R);
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0111
      *IT_0112*(IT_0113*IT_0114 + -IT_0112*IT_0115*IT_0116);
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = cpow(IT_0118, 2);
    const ccomplex_t IT_0120 = pow(m_e, 2);
    const ccomplex_t IT_0121 = 0.5*IT_0011 + -IT_0120;
    const ccomplex_t IT_0122 = IT_0011*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0013 + IT_0123;
    const ccomplex_t IT_0125 = IT_0001*IT_0124;
    const ccomplex_t IT_0126 = IT_0121 + IT_0125;
    const ccomplex_t IT_0127 = cabs(0.25*IT_0011 + -IT_0120);
    const ccomplex_t IT_0128 = cpow(IT_0127, 0.5);
    const ccomplex_t IT_0129 = ( ((m_V
) >= (2*m_e
)) ? (0.125*IT_0000*IT_0001
      *((-1.33333333333333)*IT_0119*IT_0126 + (-4)*IT_0009*(IT_0118*IT_0120 +
       0.333333333333333*IT_0009*IT_0126))*IT_0128
) : (0
));
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0111
      *IT_0115*IT_0116;
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = cpow(IT_0131, 2);
    const ccomplex_t IT_0133 = pow(m_mu, 2);
    const ccomplex_t IT_0134 = 0.5*IT_0011 + -IT_0133;
    const ccomplex_t IT_0135 = IT_0011*IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = IT_0013 + IT_0136;
    const ccomplex_t IT_0138 = IT_0001*IT_0137;
    const ccomplex_t IT_0139 = IT_0134 + IT_0138;
    const ccomplex_t IT_0140 = cabs(0.25*IT_0011 + -IT_0133);
    const ccomplex_t IT_0141 = cpow(IT_0140, 0.5);
    const ccomplex_t IT_0142 = ( ((m_V
) >= (2*m_mu
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0132*IT_0139 + (-4)*IT_0039*(IT_0131
      *IT_0133 + 0.333333333333333*IT_0039*IT_0139))*IT_0141
) : (0
));
    const ccomplex_t IT_0143 = cpow(IT_0111, 2);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0112
      *IT_0115*IT_0143;
    const ccomplex_t IT_0145 = (-0.5)*IT_0144;
    const ccomplex_t IT_0146 = cpow(IT_0145, 2);
    const ccomplex_t IT_0147 = cpow(IT_0011 + IT_0120 + -IT_0133, 2);
    const ccomplex_t IT_0148 = -IT_0120 + 0.25*IT_0001*(IT_0011 + -IT_0120 +
       IT_0133)*(IT_0011 + IT_0120 + -IT_0133) + 0.25*IT_0001*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0011 + 0.5*IT_0120 + (-0.5)*IT_0133;
    const ccomplex_t IT_0150 = 0.5*IT_0011 + (-0.5)*IT_0120 + 0.5*IT_0133;
    const ccomplex_t IT_0151 = IT_0149*IT_0150;
    const ccomplex_t IT_0152 = IT_0011*IT_0148;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = IT_0151 + IT_0153;
    const ccomplex_t IT_0155 = IT_0001*IT_0154;
    const ccomplex_t IT_0156 = IT_0148 + IT_0155;
    const ccomplex_t IT_0157 = m_e*m_mu;
    const ccomplex_t IT_0158 = cabs(-IT_0120 + 0.25*IT_0001*IT_0147);
    const ccomplex_t IT_0159 = cpow(IT_0158, 0.5);
    const ccomplex_t IT_0160 = ( ((m_V
) >= (m_e + m_mu
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0146*IT_0156 + (-4)*IT_0023*
      (0.333333333333333*IT_0023*IT_0156 + IT_0145*IT_0157))*IT_0159
) : (0
));
    const ccomplex_t IT_0161 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0116*
      (IT_0113*IT_0114 + -IT_0112*IT_0115*IT_0116);
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = cpow(IT_0162, 2);
    const ccomplex_t IT_0164 = ( ((m_V
) >= (m_e + m_mu
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0159*(IT_0027*(0.333333333333333*IT_0027*IT_0156 + IT_0157
      *IT_0162) + 0.333333333333333*IT_0156*IT_0163)
) : (0
));
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*g_f*m_V;
    const ccomplex_t IT_0166 = cpow(IT_0165, 2);
    const ccomplex_t IT_0167 = pow(m_phi, 2);
    const ccomplex_t IT_0168 = (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = IT_0011 + IT_0168;
    const ccomplex_t IT_0170 = cpow(IT_0169, 2);
    const ccomplex_t IT_0171 = IT_0056*IT_0170;
    const ccomplex_t IT_0172 = cpow(IT_0011 + (-0.5)*IT_0167, 2);
    const ccomplex_t IT_0173 = cabs(-IT_0011 + IT_0001*IT_0172);
    const ccomplex_t IT_0174 = cpow(IT_0173, 0.5);
    const ccomplex_t IT_0175 = ( ((m_V
) >= (m_V + m_phi
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0166*(1 + 0.5*IT_0171)*IT_0174
) :
       (0
));
    const ccomplex_t IT_0176 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0111
      *IT_0113*(IT_0112*IT_0114 + IT_0113*IT_0115*IT_0116);
    const ccomplex_t IT_0177 = 0.5*IT_0176;
    const ccomplex_t IT_0178 = cpow(IT_0177, 2);
    const ccomplex_t IT_0179 = pow(m_tau, 2);
    const ccomplex_t IT_0180 = 0.5*IT_0011 + -IT_0179;
    const ccomplex_t IT_0181 = IT_0011*IT_0180;
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = IT_0013 + IT_0182;
    const ccomplex_t IT_0184 = IT_0001*IT_0183;
    const ccomplex_t IT_0185 = IT_0180 + IT_0184;
    const ccomplex_t IT_0186 = cabs(0.25*IT_0011 + -IT_0179);
    const ccomplex_t IT_0187 = cpow(IT_0186, 0.5);
    const ccomplex_t IT_0188 = ( ((m_V
) >= (2*m_tau
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0178*IT_0185 + (-4)*IT_0051*(IT_0177
      *IT_0179 + 0.333333333333333*IT_0051*IT_0185))*IT_0187
) : (0
));
    const ccomplex_t IT_0189 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0111
      *IT_0112*(IT_0112*IT_0114 + IT_0113*IT_0115*IT_0116);
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = cpow(IT_0190, 2);
    const ccomplex_t IT_0192 = cpow(IT_0011 + IT_0120 + -IT_0179, 2);
    const ccomplex_t IT_0193 = -IT_0120 + 0.25*IT_0001*(IT_0011 + -IT_0120 +
       IT_0179)*(IT_0011 + IT_0120 + -IT_0179) + 0.25*IT_0001*IT_0192;
    const ccomplex_t IT_0194 = 0.5*IT_0011 + 0.5*IT_0120 + (-0.5)*IT_0179;
    const ccomplex_t IT_0195 = 0.5*IT_0011 + (-0.5)*IT_0120 + 0.5*IT_0179;
    const ccomplex_t IT_0196 = IT_0194*IT_0195;
    const ccomplex_t IT_0197 = IT_0011*IT_0193;
    const ccomplex_t IT_0198 = (-0.5)*IT_0197;
    const ccomplex_t IT_0199 = IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = IT_0001*IT_0199;
    const ccomplex_t IT_0201 = IT_0193 + IT_0200;
    const ccomplex_t IT_0202 = m_e*m_tau;
    const ccomplex_t IT_0203 = cabs(-IT_0120 + 0.25*IT_0001*IT_0192);
    const ccomplex_t IT_0204 = cpow(IT_0203, 0.5);
    const ccomplex_t IT_0205 = ( ((m_V
) >= (m_e + m_tau
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0191*IT_0201 + (-4)*IT_0031*
      (0.333333333333333*IT_0031*IT_0201 + IT_0190*IT_0202))*IT_0204
) : (0
));
    const ccomplex_t IT_0206 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0111
      *IT_0113*(IT_0113*IT_0114 + -IT_0112*IT_0115*IT_0116);
    const ccomplex_t IT_0207 = (-0.5)*IT_0206;
    const ccomplex_t IT_0208 = cpow(IT_0207, 2);
    const ccomplex_t IT_0209 = ( ((m_V
) >= (m_e + m_tau
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0204*(IT_0035*(0.333333333333333*IT_0035*IT_0201 + IT_0202
      *IT_0207) + 0.333333333333333*IT_0201*IT_0208)
) : (0
));
    const ccomplex_t IT_0210 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0116*
      (IT_0112*IT_0114 + IT_0113*IT_0115*IT_0116);
    const ccomplex_t IT_0211 = 0.5*IT_0210;
    const ccomplex_t IT_0212 = cpow(IT_0211, 2);
    const ccomplex_t IT_0213 = cpow(IT_0011 + IT_0133 + -IT_0179, 2);
    const ccomplex_t IT_0214 = -IT_0133 + 0.25*IT_0001*(IT_0011 + -IT_0133 +
       IT_0179)*(IT_0011 + IT_0133 + -IT_0179) + 0.25*IT_0001*IT_0213;
    const ccomplex_t IT_0215 = 0.5*IT_0011 + 0.5*IT_0133 + (-0.5)*IT_0179;
    const ccomplex_t IT_0216 = 0.5*IT_0011 + (-0.5)*IT_0133 + 0.5*IT_0179;
    const ccomplex_t IT_0217 = IT_0215*IT_0216;
    const ccomplex_t IT_0218 = IT_0011*IT_0214;
    const ccomplex_t IT_0219 = (-0.5)*IT_0218;
    const ccomplex_t IT_0220 = IT_0217 + IT_0219;
    const ccomplex_t IT_0221 = IT_0001*IT_0220;
    const ccomplex_t IT_0222 = IT_0214 + IT_0221;
    const ccomplex_t IT_0223 = m_mu*m_tau;
    const ccomplex_t IT_0224 = cabs(-IT_0133 + 0.25*IT_0001*IT_0213);
    const ccomplex_t IT_0225 = cpow(IT_0224, 0.5);
    const ccomplex_t IT_0226 = ( ((m_V
) >= (m_mu + m_tau
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0212*IT_0222 + (-4)*IT_0043*
      (0.333333333333333*IT_0043*IT_0222 + IT_0211*IT_0223))*IT_0225
) : (0
));
    const ccomplex_t IT_0227 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0113
      *IT_0115*IT_0143;
    const ccomplex_t IT_0228 = (-0.5)*IT_0227;
    const ccomplex_t IT_0229 = cpow(IT_0228, 2);
    const ccomplex_t IT_0230 = ( ((m_V
) >= (m_mu + m_tau
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0225*(IT_0047*(0.333333333333333*IT_0047*IT_0222 + IT_0223
      *IT_0228) + 0.333333333333333*IT_0222*IT_0229)
) : (0
));
    const ccomplex_t IT_0231 = (0 + _Complex_I*1.4142135623731)*g_f;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = cpow(IT_0232, 2);
    const ccomplex_t IT_0234 = m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0235 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0236 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0237 = cpow(IT_0011 + IT_0235 + -IT_0236, 2);
    const ccomplex_t IT_0238 = -IT_0235 + 0.25*IT_0001*(IT_0011 + -IT_0235 +
       IT_0236)*(IT_0011 + IT_0235 + -IT_0236) + 0.25*IT_0001*IT_0237;
    const ccomplex_t IT_0239 = 0.5*IT_0011 + 0.5*IT_0235 + (-0.5)*IT_0236;
    const ccomplex_t IT_0240 = 0.5*IT_0011 + (-0.5)*IT_0235 + 0.5*IT_0236;
    const ccomplex_t IT_0241 = IT_0239*IT_0240;
    const ccomplex_t IT_0242 = IT_0011*IT_0238;
    const ccomplex_t IT_0243 = (-0.5)*IT_0242;
    const ccomplex_t IT_0244 = IT_0241 + IT_0243;
    const ccomplex_t IT_0245 = IT_0001*IT_0244;
    const ccomplex_t IT_0246 = cabs(-IT_0235 + 0.25*IT_0001*IT_0237);
    const ccomplex_t IT_0247 = cpow(IT_0246, 0.5);
    const ccomplex_t IT_0248 = ( ((m_V
) >= (m_chi_dm_1 + m_chi_dm_2
)) ? ((
      -0.5)*IT_0000*IT_0001*IT_0233*(IT_0234 + 0.666666666666667*IT_0238 +
       0.666666666666667*IT_0245)*IT_0247
) : (0
));
    const ccomplex_t IT_0249 = IT_0020 + IT_0025 + IT_0029 + IT_0033 + IT_0037
       + IT_0041 + IT_0045 + IT_0049 + IT_0053 + IT_0110 + IT_0129 + IT_0142 +
       IT_0160 + IT_0164 + IT_0175 + IT_0188 + IT_0205 + IT_0209 + IT_0226 +
       IT_0230 + IT_0248;
    return create_ccomplex_return(IT_0249);
}

