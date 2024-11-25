#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_V_3.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_V_3(
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
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = sin(theta_13_L);
    const ccomplex_t IT_0012 = cpow(IT_0011, 2);
    const ccomplex_t IT_0013 = cos(theta_13_L);
    const ccomplex_t IT_0014 = sin(theta_12_L);
    const ccomplex_t IT_0015 = cos(theta_12_L);
    const ccomplex_t IT_0016 = sin(theta_23_L);
    const ccomplex_t IT_0017 = cpow(IT_0013*IT_0014 + -IT_0011*IT_0015
      *IT_0016, 2);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*g_f*(IT_0010*IT_0012 + 
      -IT_0017);
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0019, 2);
    const ccomplex_t IT_0021 = cabs(0.25*IT_0002);
    const ccomplex_t IT_0022 = cpow(IT_0021, 0.5);
    const ccomplex_t IT_0023 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0020*IT_0022
) : (0
));
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*g_f*IT_0009*(IT_0011*IT_0016
       + -IT_0015*(IT_0013*IT_0014 + -IT_0011*IT_0015*IT_0016));
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0025, 2);
    const ccomplex_t IT_0027 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0022*IT_0026
) : (0
));
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*g_f*(IT_0010*IT_0011*IT_0013
       + (IT_0013*IT_0014 + -IT_0011*IT_0015*IT_0016)*(IT_0011*IT_0014 + IT_0013
      *IT_0015*IT_0016));
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0029, 2);
    const ccomplex_t IT_0031 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0022*IT_0030
) : (0
));
    const ccomplex_t IT_0032 = cpow(IT_0015, 2);
    const ccomplex_t IT_0033 = cpow(IT_0016, 2);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*g_f*(IT_0010*IT_0032 + 
      -IT_0033);
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0035, 2);
    const ccomplex_t IT_0037 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0022*IT_0036
) : (0
));
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*g_f*IT_0009*(IT_0013*IT_0016
       + IT_0015*(IT_0011*IT_0014 + IT_0013*IT_0015*IT_0016));
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0039, 2);
    const ccomplex_t IT_0041 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0022*IT_0040
) : (0
));
    const ccomplex_t IT_0042 = cpow(IT_0013, 2);
    const ccomplex_t IT_0043 = cpow(IT_0011*IT_0014 + IT_0013*IT_0015*IT_0016,
       2);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*g_f*(IT_0010*IT_0042 + 
      -IT_0043);
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = cpow(IT_0045, 2);
    const ccomplex_t IT_0047 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0022*IT_0046
) : (0
));
    const ccomplex_t IT_0048 = IT_0001*IT_0004;
    const ccomplex_t IT_0049 = pow(m_V, -4);
    const ccomplex_t IT_0050 = (-0.5)*IT_0002;
    const ccomplex_t IT_0051 = cpow(IT_0050, 2);
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = (IT_0002 + -IT_0048)*(1 + 0.5*IT_0052);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0054 + -IT_0055;
    const ccomplex_t IT_0057 = cpow(IT_0056, 2);
    const ccomplex_t IT_0058 = IT_0004*IT_0049;
    const ccomplex_t IT_0059 = (IT_0002 + -IT_0048)*(1 + 0.5*IT_0058);
    const ccomplex_t IT_0060 = IT_0003*IT_0049;
    const ccomplex_t IT_0061 = IT_0003*IT_0060;
    const ccomplex_t IT_0062 = IT_0001*IT_0003;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0050 + IT_0064;
    const ccomplex_t IT_0066 = IT_0061*IT_0065;
    const ccomplex_t IT_0067 = (-2)*IT_0054;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0001*IT_0050;
    const ccomplex_t IT_0070 = IT_0003*IT_0069;
    const ccomplex_t IT_0071 = (1 + 0.5*IT_0058)*(IT_0003 + -IT_0070);
    const ccomplex_t IT_0072 = 2*IT_0054;
    const ccomplex_t IT_0073 = IT_0054 + -IT_0072;
    const ccomplex_t IT_0074 = (-2)*IT_0058;
    const ccomplex_t IT_0075 = pow(m_V, -6);
    const ccomplex_t IT_0076 = IT_0050*IT_0075;
    const ccomplex_t IT_0077 = IT_0003*IT_0076;
    const ccomplex_t IT_0078 = IT_0003*IT_0077;
    const ccomplex_t IT_0079 = IT_0074 + IT_0078;
    const ccomplex_t IT_0080 = 1 + IT_0079;
    const ccomplex_t IT_0081 = IT_0003*IT_0080;
    const ccomplex_t IT_0082 = IT_0001*IT_0051;
    const ccomplex_t IT_0083 = (1 + 0.5*IT_0058)*(IT_0002 + -IT_0082);
    const ccomplex_t IT_0084 = IT_0049*IT_0050;
    const ccomplex_t IT_0085 = IT_0003*IT_0084;
    const ccomplex_t IT_0086 = -IT_0070;
    const ccomplex_t IT_0087 = IT_0003 + IT_0086;
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = -IT_0052 + -IT_0058;
    const ccomplex_t IT_0090 = IT_0078 + IT_0089;
    const ccomplex_t IT_0091 = 1 + IT_0090;
    const ccomplex_t IT_0092 = IT_0003*IT_0091;
    const ccomplex_t IT_0093 = cpow(IT_0003, 3);
    const ccomplex_t IT_0094 = IT_0049*IT_0093;
    const ccomplex_t IT_0095 = IT_0003*IT_0052;
    const ccomplex_t IT_0096 = IT_0003*IT_0058;
    const ccomplex_t IT_0097 = IT_0075*IT_0093;
    const ccomplex_t IT_0098 = IT_0050*IT_0097;
    const ccomplex_t IT_0099 = 0.666666666666667*IT_0053*IT_0057 +
       0.666666666666667*IT_0055*(IT_0055*IT_0059 + IT_0066*IT_0068 + 2*IT_0071
      *IT_0073 + -IT_0056*IT_0081) + 0.666666666666667*IT_0073*(IT_0073*IT_0083 
      + IT_0056*IT_0088) + 0.666666666666667*IT_0068*(IT_0059*IT_0068 + -IT_0073
      *IT_0092 + -IT_0056*(IT_0003 + IT_0070 + -IT_0094 + -IT_0095 + -IT_0096 +
       IT_0098));
    const ccomplex_t IT_0100 = cabs((-0.75)*IT_0002);
    const ccomplex_t IT_0101 = cpow(IT_0100, 0.5);
    const ccomplex_t IT_0102 = ( ((m_V
) >= (2*m_V
)) ? (0.125*IT_0000*IT_0001
      *IT_0099*IT_0101
) : (0
));
    const ccomplex_t IT_0103 = cos(theta_23_R);
    const ccomplex_t IT_0104 = cpow(IT_0103, 2);
    const ccomplex_t IT_0105 = sin(theta_13_R);
    const ccomplex_t IT_0106 = cpow(IT_0105, 2);
    const ccomplex_t IT_0107 = cos(theta_13_R);
    const ccomplex_t IT_0108 = sin(theta_12_R);
    const ccomplex_t IT_0109 = cos(theta_12_R);
    const ccomplex_t IT_0110 = sin(theta_23_R);
    const ccomplex_t IT_0111 = cpow(IT_0107*IT_0108 + -IT_0105*IT_0109
      *IT_0110, 2);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*g_f*(IT_0104*IT_0106 + 
      -IT_0111);
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = cpow(IT_0113, 2);
    const ccomplex_t IT_0115 = pow(m_e, 2);
    const ccomplex_t IT_0116 = 0.5*IT_0002 + -IT_0115;
    const ccomplex_t IT_0117 = IT_0002*IT_0116;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = IT_0004 + IT_0118;
    const ccomplex_t IT_0120 = IT_0001*IT_0119;
    const ccomplex_t IT_0121 = IT_0116 + IT_0120;
    const ccomplex_t IT_0122 = (-1.33333333333333)*IT_0121;
    const ccomplex_t IT_0123 = cabs(0.25*IT_0002 + -IT_0115);
    const ccomplex_t IT_0124 = cpow(IT_0123, 0.5);
    const ccomplex_t IT_0125 = ( ((m_V
) >= (2*m_e
)) ? (0.125*IT_0000*IT_0001
      *(IT_0114*IT_0122 + (-4)*IT_0019*(IT_0113*IT_0115 + (-0.25)*IT_0019
      *IT_0122))*IT_0124
) : (0
));
    const ccomplex_t IT_0126 = cpow(IT_0109, 2);
    const ccomplex_t IT_0127 = cpow(IT_0110, 2);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*g_f*(IT_0104*IT_0126 + 
      -IT_0127);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = cpow(IT_0129, 2);
    const ccomplex_t IT_0131 = pow(m_mu, 2);
    const ccomplex_t IT_0132 = 0.5*IT_0002 + -IT_0131;
    const ccomplex_t IT_0133 = IT_0002*IT_0132;
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = IT_0004 + IT_0134;
    const ccomplex_t IT_0136 = IT_0001*IT_0135;
    const ccomplex_t IT_0137 = IT_0132 + IT_0136;
    const ccomplex_t IT_0138 = (-1.33333333333333)*IT_0137;
    const ccomplex_t IT_0139 = cabs(0.25*IT_0002 + -IT_0131);
    const ccomplex_t IT_0140 = cpow(IT_0139, 0.5);
    const ccomplex_t IT_0141 = ( ((m_V
) >= (2*m_mu
)) ? (0.125*IT_0000
      *IT_0001*(IT_0130*IT_0138 + (-4)*IT_0035*(IT_0129*IT_0131 + (-0.25)
      *IT_0035*IT_0138))*IT_0140
) : (0
));
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*g_f*IT_0103*(IT_0105*IT_0110
       + -IT_0109*(IT_0107*IT_0108 + -IT_0105*IT_0109*IT_0110));
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = cpow(IT_0143, 2);
    const ccomplex_t IT_0145 = cpow(IT_0002 + IT_0115 + -IT_0131, 2);
    const ccomplex_t IT_0146 = -IT_0115 + 0.25*IT_0001*(IT_0002 + -IT_0115 +
       IT_0131)*(IT_0002 + IT_0115 + -IT_0131) + 0.25*IT_0001*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0002 + 0.5*IT_0115 + (-0.5)*IT_0131;
    const ccomplex_t IT_0148 = 0.5*IT_0002 + (-0.5)*IT_0115 + 0.5*IT_0131;
    const ccomplex_t IT_0149 = IT_0147*IT_0148;
    const ccomplex_t IT_0150 = IT_0002*IT_0146;
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = IT_0149 + IT_0151;
    const ccomplex_t IT_0153 = IT_0001*IT_0152;
    const ccomplex_t IT_0154 = IT_0146 + IT_0153;
    const ccomplex_t IT_0155 = (-1.33333333333333)*IT_0154;
    const ccomplex_t IT_0156 = m_e*m_mu;
    const ccomplex_t IT_0157 = cabs(-IT_0115 + 0.25*IT_0001*IT_0145);
    const ccomplex_t IT_0158 = cpow(IT_0157, 0.5);
    const ccomplex_t IT_0159 = ( ((m_V
) >= (m_e + m_mu
)) ? (0.125*IT_0000
      *IT_0001*(IT_0144*IT_0155 + (-4)*IT_0025*((-0.25)*IT_0025*IT_0155 +
       IT_0143*IT_0156))*IT_0158
) : (0
));
    const ccomplex_t IT_0160 = pow(m_phi, 2);
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = IT_0002 + IT_0161;
    const ccomplex_t IT_0163 = cpow(IT_0162, 2);
    const ccomplex_t IT_0164 = IT_0049*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*g_f*m_V;
    const ccomplex_t IT_0166 = cpow(IT_0165, 2);
    const ccomplex_t IT_0167 = cpow(IT_0002 + (-0.5)*IT_0160, 2);
    const ccomplex_t IT_0168 = cabs(-IT_0002 + IT_0001*IT_0167);
    const ccomplex_t IT_0169 = cpow(IT_0168, 0.5);
    const ccomplex_t IT_0170 = ( ((m_V
) >= (m_V + m_phi
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*(1 + 0.5*IT_0164)*IT_0166*IT_0169
) :
       (0
));
    const ccomplex_t IT_0171 = cpow(IT_0107, 2);
    const ccomplex_t IT_0172 = cpow(IT_0105*IT_0108 + IT_0107*IT_0109*IT_0110,
       2);
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*g_f*(IT_0104*IT_0171 + 
      -IT_0172);
    const ccomplex_t IT_0174 = 0.5*IT_0173;
    const ccomplex_t IT_0175 = cpow(IT_0174, 2);
    const ccomplex_t IT_0176 = pow(m_tau, 2);
    const ccomplex_t IT_0177 = 0.5*IT_0002 + -IT_0176;
    const ccomplex_t IT_0178 = IT_0002*IT_0177;
    const ccomplex_t IT_0179 = (-0.5)*IT_0178;
    const ccomplex_t IT_0180 = IT_0004 + IT_0179;
    const ccomplex_t IT_0181 = IT_0001*IT_0180;
    const ccomplex_t IT_0182 = IT_0177 + IT_0181;
    const ccomplex_t IT_0183 = (-1.33333333333333)*IT_0182;
    const ccomplex_t IT_0184 = cabs(0.25*IT_0002 + -IT_0176);
    const ccomplex_t IT_0185 = cpow(IT_0184, 0.5);
    const ccomplex_t IT_0186 = ( ((m_V
) >= (2*m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0175*IT_0183 + (-4)*IT_0045*(IT_0174*IT_0176 + (-0.25)
      *IT_0045*IT_0183))*IT_0185
) : (0
));
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*g_f*(IT_0104*IT_0105*IT_0107
       + (IT_0107*IT_0108 + -IT_0105*IT_0109*IT_0110)*(IT_0105*IT_0108 + IT_0107
      *IT_0109*IT_0110));
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = cpow(IT_0188, 2);
    const ccomplex_t IT_0190 = cpow(IT_0002 + IT_0115 + -IT_0176, 2);
    const ccomplex_t IT_0191 = -IT_0115 + 0.25*IT_0001*(IT_0002 + -IT_0115 +
       IT_0176)*(IT_0002 + IT_0115 + -IT_0176) + 0.25*IT_0001*IT_0190;
    const ccomplex_t IT_0192 = 0.5*IT_0002 + 0.5*IT_0115 + (-0.5)*IT_0176;
    const ccomplex_t IT_0193 = 0.5*IT_0002 + (-0.5)*IT_0115 + 0.5*IT_0176;
    const ccomplex_t IT_0194 = IT_0192*IT_0193;
    const ccomplex_t IT_0195 = IT_0002*IT_0191;
    const ccomplex_t IT_0196 = (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = IT_0001*IT_0197;
    const ccomplex_t IT_0199 = IT_0191 + IT_0198;
    const ccomplex_t IT_0200 = (-1.33333333333333)*IT_0199;
    const ccomplex_t IT_0201 = m_e*m_tau;
    const ccomplex_t IT_0202 = cabs(-IT_0115 + 0.25*IT_0001*IT_0190);
    const ccomplex_t IT_0203 = cpow(IT_0202, 0.5);
    const ccomplex_t IT_0204 = ( ((m_V
) >= (m_e + m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0189*IT_0200 + (-4)*IT_0029*((-0.25)*IT_0029*IT_0200 +
       IT_0188*IT_0201))*IT_0203
) : (0
));
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*g_f*IT_0103*(IT_0107*IT_0110
       + IT_0109*(IT_0105*IT_0108 + IT_0107*IT_0109*IT_0110));
    const ccomplex_t IT_0206 = 0.5*IT_0205;
    const ccomplex_t IT_0207 = cpow(IT_0206, 2);
    const ccomplex_t IT_0208 = cpow(IT_0002 + IT_0131 + -IT_0176, 2);
    const ccomplex_t IT_0209 = -IT_0131 + 0.25*IT_0001*(IT_0002 + -IT_0131 +
       IT_0176)*(IT_0002 + IT_0131 + -IT_0176) + 0.25*IT_0001*IT_0208;
    const ccomplex_t IT_0210 = 0.5*IT_0002 + 0.5*IT_0131 + (-0.5)*IT_0176;
    const ccomplex_t IT_0211 = 0.5*IT_0002 + (-0.5)*IT_0131 + 0.5*IT_0176;
    const ccomplex_t IT_0212 = IT_0210*IT_0211;
    const ccomplex_t IT_0213 = IT_0002*IT_0209;
    const ccomplex_t IT_0214 = (-0.5)*IT_0213;
    const ccomplex_t IT_0215 = IT_0212 + IT_0214;
    const ccomplex_t IT_0216 = IT_0001*IT_0215;
    const ccomplex_t IT_0217 = IT_0209 + IT_0216;
    const ccomplex_t IT_0218 = (-1.33333333333333)*IT_0217;
    const ccomplex_t IT_0219 = m_mu*m_tau;
    const ccomplex_t IT_0220 = cabs(-IT_0131 + 0.25*IT_0001*IT_0208);
    const ccomplex_t IT_0221 = cpow(IT_0220, 0.5);
    const ccomplex_t IT_0222 = ( ((m_V
) >= (m_mu + m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0207*IT_0218 + (-4)*IT_0039*((-0.25)*IT_0039*IT_0218 +
       IT_0206*IT_0219))*IT_0221
) : (0
));
    const ccomplex_t IT_0223 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0224 = 0.5*IT_0002 + -IT_0223;
    const ccomplex_t IT_0225 = IT_0002*IT_0224;
    const ccomplex_t IT_0226 = (-0.5)*IT_0225;
    const ccomplex_t IT_0227 = IT_0004 + IT_0226;
    const ccomplex_t IT_0228 = IT_0001*IT_0227;
    const ccomplex_t IT_0229 = 0.5*IT_0054;
    const ccomplex_t IT_0230 = cpow(IT_0229, 2);
    const ccomplex_t IT_0231 = cabs(0.25*IT_0002 + -IT_0223);
    const ccomplex_t IT_0232 = cpow(IT_0231, 0.5);
    const ccomplex_t IT_0233 = ( ((m_V
) >= (2*m_chi_dm_1
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0223 + 0.666666666666667*IT_0224 + 0.666666666666667*IT_0228)
      *IT_0230*IT_0232
) : (0
));
    const ccomplex_t IT_0234 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0235 = 0.5*IT_0002 + -IT_0234;
    const ccomplex_t IT_0236 = IT_0002*IT_0235;
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0004 + IT_0237;
    const ccomplex_t IT_0239 = IT_0001*IT_0238;
    const ccomplex_t IT_0240 = (-0.5)*IT_0054;
    const ccomplex_t IT_0241 = cpow(IT_0240, 2);
    const ccomplex_t IT_0242 = cabs(0.25*IT_0002 + -IT_0234);
    const ccomplex_t IT_0243 = cpow(IT_0242, 0.5);
    const ccomplex_t IT_0244 = ( ((m_V
) >= (2*m_chi_dm_2
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0234 + 0.666666666666667*IT_0235 + 0.666666666666667*IT_0239)
      *IT_0241*IT_0243
) : (0
));
    const ccomplex_t IT_0245 = IT_0023 + 2*IT_0027 + 2*IT_0031 + IT_0037 + 2
      *IT_0041 + IT_0047 + IT_0102 + IT_0125 + IT_0141 + 2*IT_0159 + IT_0170 +
       IT_0186 + 2*IT_0204 + 2*IT_0222 + IT_0233 + IT_0244;
    return create_ccomplex_return(IT_0245);
}

