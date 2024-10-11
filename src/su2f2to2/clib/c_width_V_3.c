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
    const ccomplex_t IT_0048 = pow(m_V, -4);
    const ccomplex_t IT_0049 = IT_0004*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0002;
    const ccomplex_t IT_0051 = cpow(IT_0050, 2);
    const ccomplex_t IT_0052 = IT_0001*IT_0051;
    const ccomplex_t IT_0053 = (1 + 0.5*IT_0049)*(IT_0002 + -IT_0052);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0055 = 2*IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = cpow(IT_0056, 2);
    const ccomplex_t IT_0058 = -IT_0054;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = IT_0001*IT_0004;
    const ccomplex_t IT_0061 = IT_0048*IT_0051;
    const ccomplex_t IT_0062 = (IT_0002 + -IT_0060)*(1 + 0.5*IT_0061);
    const ccomplex_t IT_0063 = pow(m_V, -6);
    const ccomplex_t IT_0064 = IT_0050*IT_0063;
    const ccomplex_t IT_0065 = IT_0003*IT_0064;
    const ccomplex_t IT_0066 = IT_0003*IT_0065;
    const ccomplex_t IT_0067 = -IT_0049 + -IT_0061;
    const ccomplex_t IT_0068 = IT_0066 + IT_0067;
    const ccomplex_t IT_0069 = 1 + IT_0068;
    const ccomplex_t IT_0070 = IT_0050*IT_0069;
    const ccomplex_t IT_0071 = (-2)*IT_0054;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0048*IT_0050;
    const ccomplex_t IT_0074 = IT_0003*IT_0073;
    const ccomplex_t IT_0075 = IT_0001*IT_0050;
    const ccomplex_t IT_0076 = IT_0003*IT_0075;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0003 + IT_0077;
    const ccomplex_t IT_0079 = IT_0074*IT_0078;
    const ccomplex_t IT_0080 = IT_0001*IT_0003;
    const ccomplex_t IT_0081 = IT_0003*IT_0080;
    const ccomplex_t IT_0082 = cpow(IT_0050, 3);
    const ccomplex_t IT_0083 = IT_0048*IT_0082;
    const ccomplex_t IT_0084 = IT_0049*IT_0050;
    const ccomplex_t IT_0085 = IT_0051*IT_0063;
    const ccomplex_t IT_0086 = IT_0003*IT_0085;
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = (1 + 0.5*IT_0061)*(IT_0050 + -IT_0081);
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = 0.666666666666667*IT_0053*IT_0057 +
       0.666666666666667*IT_0059*(IT_0059*IT_0062 + -IT_0070*IT_0072 + IT_0056
      *IT_0079) + 0.666666666666667*IT_0072*(IT_0053*IT_0072 + -IT_0056*(IT_0050
       + IT_0081 + -IT_0083 + (-2)*IT_0084 + IT_0087)) + 0.666666666666667
      *IT_0058*(IT_0058*IT_0062 + -IT_0056*IT_0070 + IT_0072*IT_0079 + (-2)
      *IT_0059*IT_0089);
    const ccomplex_t IT_0091 = cabs((-0.75)*IT_0002);
    const ccomplex_t IT_0092 = cpow(IT_0091, 0.5);
    const ccomplex_t IT_0093 = ( ((m_V
) >= (2*m_V
)) ? (0.125*IT_0000*IT_0001
      *IT_0090*IT_0092
) : (0
));
    const ccomplex_t IT_0094 = cos(theta_23_R);
    const ccomplex_t IT_0095 = cpow(IT_0094, 2);
    const ccomplex_t IT_0096 = sin(theta_13_R);
    const ccomplex_t IT_0097 = cpow(IT_0096, 2);
    const ccomplex_t IT_0098 = cos(theta_13_R);
    const ccomplex_t IT_0099 = sin(theta_12_R);
    const ccomplex_t IT_0100 = cos(theta_12_R);
    const ccomplex_t IT_0101 = sin(theta_23_R);
    const ccomplex_t IT_0102 = cpow(IT_0098*IT_0099 + -IT_0096*IT_0100
      *IT_0101, 2);
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*g_f*(IT_0095*IT_0097 + 
      -IT_0102);
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = cpow(IT_0104, 2);
    const ccomplex_t IT_0106 = pow(m_e, 2);
    const ccomplex_t IT_0107 = 0.5*IT_0002 + -IT_0106;
    const ccomplex_t IT_0108 = IT_0002*IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = IT_0004 + IT_0109;
    const ccomplex_t IT_0111 = IT_0001*IT_0110;
    const ccomplex_t IT_0112 = IT_0107 + IT_0111;
    const ccomplex_t IT_0113 = (-1.33333333333333)*IT_0112;
    const ccomplex_t IT_0114 = cabs(0.25*IT_0002 + -IT_0106);
    const ccomplex_t IT_0115 = cpow(IT_0114, 0.5);
    const ccomplex_t IT_0116 = ( ((m_V
) >= (2*m_e
)) ? (0.125*IT_0000*IT_0001
      *(IT_0105*IT_0113 + (-4)*IT_0019*(IT_0104*IT_0106 + (-0.25)*IT_0019
      *IT_0113))*IT_0115
) : (0
));
    const ccomplex_t IT_0117 = cpow(IT_0100, 2);
    const ccomplex_t IT_0118 = cpow(IT_0101, 2);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*g_f*(IT_0095*IT_0117 + 
      -IT_0118);
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = cpow(IT_0120, 2);
    const ccomplex_t IT_0122 = pow(m_mu, 2);
    const ccomplex_t IT_0123 = 0.5*IT_0002 + -IT_0122;
    const ccomplex_t IT_0124 = IT_0002*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0004 + IT_0125;
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = IT_0123 + IT_0127;
    const ccomplex_t IT_0129 = (-1.33333333333333)*IT_0128;
    const ccomplex_t IT_0130 = cabs(0.25*IT_0002 + -IT_0122);
    const ccomplex_t IT_0131 = cpow(IT_0130, 0.5);
    const ccomplex_t IT_0132 = ( ((m_V
) >= (2*m_mu
)) ? (0.125*IT_0000
      *IT_0001*(IT_0121*IT_0129 + (-4)*IT_0035*(IT_0120*IT_0122 + (-0.25)
      *IT_0035*IT_0129))*IT_0131
) : (0
));
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*g_f*IT_0094*(IT_0096*IT_0101
       + -IT_0100*(IT_0098*IT_0099 + -IT_0096*IT_0100*IT_0101));
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = cpow(IT_0134, 2);
    const ccomplex_t IT_0136 = cpow(IT_0002 + IT_0106 + -IT_0122, 2);
    const ccomplex_t IT_0137 = -IT_0106 + 0.25*IT_0001*(IT_0002 + -IT_0106 +
       IT_0122)*(IT_0002 + IT_0106 + -IT_0122) + 0.25*IT_0001*IT_0136;
    const ccomplex_t IT_0138 = 0.5*IT_0002 + 0.5*IT_0106 + (-0.5)*IT_0122;
    const ccomplex_t IT_0139 = 0.5*IT_0002 + (-0.5)*IT_0106 + 0.5*IT_0122;
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = IT_0002*IT_0137;
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = IT_0140 + IT_0142;
    const ccomplex_t IT_0144 = IT_0001*IT_0143;
    const ccomplex_t IT_0145 = IT_0137 + IT_0144;
    const ccomplex_t IT_0146 = (-1.33333333333333)*IT_0145;
    const ccomplex_t IT_0147 = m_e*m_mu;
    const ccomplex_t IT_0148 = cabs(-IT_0106 + 0.25*IT_0001*IT_0136);
    const ccomplex_t IT_0149 = cpow(IT_0148, 0.5);
    const ccomplex_t IT_0150 = ( ((m_V
) >= (m_e + m_mu
)) ? (0.125*IT_0000
      *IT_0001*(IT_0135*IT_0146 + (-4)*IT_0025*((-0.25)*IT_0025*IT_0146 +
       IT_0134*IT_0147))*IT_0149
) : (0
));
    const ccomplex_t IT_0151 = pow(m_phi, 2);
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0002 + IT_0152;
    const ccomplex_t IT_0154 = cpow(IT_0153, 2);
    const ccomplex_t IT_0155 = IT_0048*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*g_f*m_V;
    const ccomplex_t IT_0157 = cpow(IT_0156, 2);
    const ccomplex_t IT_0158 = cpow(IT_0002 + (-0.5)*IT_0151, 2);
    const ccomplex_t IT_0159 = cabs(-IT_0002 + IT_0001*IT_0158);
    const ccomplex_t IT_0160 = cpow(IT_0159, 0.5);
    const ccomplex_t IT_0161 = ( ((m_V
) >= (m_V + m_phi
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*(1 + 0.5*IT_0155)*IT_0157*IT_0160
) :
       (0
));
    const ccomplex_t IT_0162 = cpow(IT_0098, 2);
    const ccomplex_t IT_0163 = cpow(IT_0096*IT_0099 + IT_0098*IT_0100*IT_0101,
       2);
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*g_f*(IT_0095*IT_0162 + 
      -IT_0163);
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = cpow(IT_0165, 2);
    const ccomplex_t IT_0167 = pow(m_tau, 2);
    const ccomplex_t IT_0168 = 0.5*IT_0002 + -IT_0167;
    const ccomplex_t IT_0169 = IT_0002*IT_0168;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = IT_0004 + IT_0170;
    const ccomplex_t IT_0172 = IT_0001*IT_0171;
    const ccomplex_t IT_0173 = IT_0168 + IT_0172;
    const ccomplex_t IT_0174 = (-1.33333333333333)*IT_0173;
    const ccomplex_t IT_0175 = cabs(0.25*IT_0002 + -IT_0167);
    const ccomplex_t IT_0176 = cpow(IT_0175, 0.5);
    const ccomplex_t IT_0177 = ( ((m_V
) >= (2*m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0166*IT_0174 + (-4)*IT_0045*(IT_0165*IT_0167 + (-0.25)
      *IT_0045*IT_0174))*IT_0176
) : (0
));
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*g_f*(IT_0095*IT_0096*IT_0098
       + (IT_0098*IT_0099 + -IT_0096*IT_0100*IT_0101)*(IT_0096*IT_0099 + IT_0098
      *IT_0100*IT_0101));
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = cpow(IT_0179, 2);
    const ccomplex_t IT_0181 = cpow(IT_0002 + IT_0106 + -IT_0167, 2);
    const ccomplex_t IT_0182 = -IT_0106 + 0.25*IT_0001*(IT_0002 + -IT_0106 +
       IT_0167)*(IT_0002 + IT_0106 + -IT_0167) + 0.25*IT_0001*IT_0181;
    const ccomplex_t IT_0183 = 0.5*IT_0002 + 0.5*IT_0106 + (-0.5)*IT_0167;
    const ccomplex_t IT_0184 = 0.5*IT_0002 + (-0.5)*IT_0106 + 0.5*IT_0167;
    const ccomplex_t IT_0185 = IT_0183*IT_0184;
    const ccomplex_t IT_0186 = IT_0002*IT_0182;
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = IT_0185 + IT_0187;
    const ccomplex_t IT_0189 = IT_0001*IT_0188;
    const ccomplex_t IT_0190 = IT_0182 + IT_0189;
    const ccomplex_t IT_0191 = (-1.33333333333333)*IT_0190;
    const ccomplex_t IT_0192 = m_e*m_tau;
    const ccomplex_t IT_0193 = cabs(-IT_0106 + 0.25*IT_0001*IT_0181);
    const ccomplex_t IT_0194 = cpow(IT_0193, 0.5);
    const ccomplex_t IT_0195 = ( ((m_V
) >= (m_e + m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0180*IT_0191 + (-4)*IT_0029*((-0.25)*IT_0029*IT_0191 +
       IT_0179*IT_0192))*IT_0194
) : (0
));
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*g_f*IT_0094*(IT_0098*IT_0101
       + IT_0100*(IT_0096*IT_0099 + IT_0098*IT_0100*IT_0101));
    const ccomplex_t IT_0197 = 0.5*IT_0196;
    const ccomplex_t IT_0198 = cpow(IT_0197, 2);
    const ccomplex_t IT_0199 = cpow(IT_0002 + IT_0122 + -IT_0167, 2);
    const ccomplex_t IT_0200 = -IT_0122 + 0.25*IT_0001*(IT_0002 + -IT_0122 +
       IT_0167)*(IT_0002 + IT_0122 + -IT_0167) + 0.25*IT_0001*IT_0199;
    const ccomplex_t IT_0201 = 0.5*IT_0002 + 0.5*IT_0122 + (-0.5)*IT_0167;
    const ccomplex_t IT_0202 = 0.5*IT_0002 + (-0.5)*IT_0122 + 0.5*IT_0167;
    const ccomplex_t IT_0203 = IT_0201*IT_0202;
    const ccomplex_t IT_0204 = IT_0002*IT_0200;
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = IT_0001*IT_0206;
    const ccomplex_t IT_0208 = IT_0200 + IT_0207;
    const ccomplex_t IT_0209 = (-1.33333333333333)*IT_0208;
    const ccomplex_t IT_0210 = m_mu*m_tau;
    const ccomplex_t IT_0211 = cabs(-IT_0122 + 0.25*IT_0001*IT_0199);
    const ccomplex_t IT_0212 = cpow(IT_0211, 0.5);
    const ccomplex_t IT_0213 = ( ((m_V
) >= (m_mu + m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0198*IT_0209 + (-4)*IT_0039*((-0.25)*IT_0039*IT_0209 +
       IT_0197*IT_0210))*IT_0212
) : (0
));
    const ccomplex_t IT_0214 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0215 = 0.5*IT_0002 + -IT_0214;
    const ccomplex_t IT_0216 = IT_0002*IT_0215;
    const ccomplex_t IT_0217 = (-0.5)*IT_0216;
    const ccomplex_t IT_0218 = IT_0004 + IT_0217;
    const ccomplex_t IT_0219 = IT_0001*IT_0218;
    const ccomplex_t IT_0220 = 0.5*IT_0054;
    const ccomplex_t IT_0221 = cpow(IT_0220, 2);
    const ccomplex_t IT_0222 = cabs(0.25*IT_0002 + -IT_0214);
    const ccomplex_t IT_0223 = cpow(IT_0222, 0.5);
    const ccomplex_t IT_0224 = ( ((m_V
) >= (2*m_chi_dm_1
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0214 + 0.666666666666667*IT_0215 + 0.666666666666667*IT_0219)
      *IT_0221*IT_0223
) : (0
));
    const ccomplex_t IT_0225 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0226 = 0.5*IT_0002 + -IT_0225;
    const ccomplex_t IT_0227 = IT_0002*IT_0226;
    const ccomplex_t IT_0228 = (-0.5)*IT_0227;
    const ccomplex_t IT_0229 = IT_0004 + IT_0228;
    const ccomplex_t IT_0230 = IT_0001*IT_0229;
    const ccomplex_t IT_0231 = (-0.5)*IT_0054;
    const ccomplex_t IT_0232 = cpow(IT_0231, 2);
    const ccomplex_t IT_0233 = cabs(0.25*IT_0002 + -IT_0225);
    const ccomplex_t IT_0234 = cpow(IT_0233, 0.5);
    const ccomplex_t IT_0235 = ( ((m_V
) >= (2*m_chi_dm_2
)) ? ((-0.5)*IT_0000
      *IT_0001*(IT_0225 + 0.666666666666667*IT_0226 + 0.666666666666667*IT_0230)
      *IT_0232*IT_0234
) : (0
));
    const ccomplex_t IT_0236 = IT_0023 + 2*IT_0027 + 2*IT_0031 + IT_0037 + 2
      *IT_0041 + IT_0047 + IT_0093 + IT_0116 + IT_0132 + 2*IT_0150 + IT_0161 +
       IT_0177 + 2*IT_0195 + 2*IT_0213 + IT_0224 + IT_0235;
    return create_ccomplex_return(IT_0236);
}

