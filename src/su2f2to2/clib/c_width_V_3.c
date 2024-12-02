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
    const ccomplex_t IT_0002 = cos(theta_23_L);
    const ccomplex_t IT_0003 = cpow(IT_0002, 2);
    const ccomplex_t IT_0004 = sin(theta_13_L);
    const ccomplex_t IT_0005 = cpow(IT_0004, 2);
    const ccomplex_t IT_0006 = cos(theta_13_L);
    const ccomplex_t IT_0007 = sin(theta_12_L);
    const ccomplex_t IT_0008 = cos(theta_12_L);
    const ccomplex_t IT_0009 = sin(theta_23_L);
    const ccomplex_t IT_0010 = cpow(IT_0006*IT_0007 + -IT_0004*IT_0008
      *IT_0009, 2);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*g_f*(IT_0003*IT_0005 + 
      -IT_0010);
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0012, 2);
    const ccomplex_t IT_0014 = pow(m_V, 2);
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = cpow(IT_0015, 2);
    const ccomplex_t IT_0017 = IT_0014*IT_0015;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0016 + IT_0018;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = cabs(0.25*IT_0014);
    const ccomplex_t IT_0022 = cpow(IT_0021, 0.5);
    const ccomplex_t IT_0023 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*IT_0013*(IT_0015 + IT_0020)*IT_0022
) : (0
));
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*g_f*IT_0002*(IT_0004*IT_0009
       + -IT_0008*(IT_0006*IT_0007 + -IT_0004*IT_0008*IT_0009));
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0025, 2);
    const ccomplex_t IT_0027 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0015 + IT_0020)*IT_0022*IT_0026
) : (0
));
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*g_f*(IT_0003*IT_0004*IT_0006
       + (IT_0006*IT_0007 + -IT_0004*IT_0008*IT_0009)*(IT_0004*IT_0007 + IT_0006
      *IT_0008*IT_0009));
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0029, 2);
    const ccomplex_t IT_0031 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0015 + IT_0020)*IT_0022*IT_0030
) : (0
));
    const ccomplex_t IT_0032 = cpow(IT_0008, 2);
    const ccomplex_t IT_0033 = cpow(IT_0009, 2);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*g_f*(IT_0003*IT_0032 + 
      -IT_0033);
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0035, 2);
    const ccomplex_t IT_0037 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0015 + IT_0020)*IT_0022*IT_0036
) : (0
));
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*g_f*IT_0002*(IT_0006*IT_0009
       + IT_0008*(IT_0004*IT_0007 + IT_0006*IT_0008*IT_0009));
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0039, 2);
    const ccomplex_t IT_0041 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0015 + IT_0020)*IT_0022*IT_0040
) : (0
));
    const ccomplex_t IT_0042 = cpow(IT_0006, 2);
    const ccomplex_t IT_0043 = cpow(IT_0004*IT_0007 + IT_0006*IT_0008*IT_0009,
       2);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*g_f*(IT_0003*IT_0042 + 
      -IT_0043);
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = cpow(IT_0045, 2);
    const ccomplex_t IT_0047 = ( ((m_V
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0015 + IT_0020)*IT_0022*IT_0046
) : (0
));
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0049 = (-2)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = cpow(IT_0050, 2);
    const ccomplex_t IT_0052 = pow(m_V, -4);
    const ccomplex_t IT_0053 = IT_0016*IT_0052;
    const ccomplex_t IT_0054 = IT_0001*IT_0016;
    const ccomplex_t IT_0055 = (1 + 0.5*IT_0053)*(IT_0014 + -IT_0054);
    const ccomplex_t IT_0056 = -IT_0048;
    const ccomplex_t IT_0057 = 2*IT_0048;
    const ccomplex_t IT_0058 = IT_0048 + -IT_0057;
    const ccomplex_t IT_0059 = (-0.5)*IT_0014;
    const ccomplex_t IT_0060 = IT_0001*IT_0059;
    const ccomplex_t IT_0061 = IT_0015*IT_0060;
    const ccomplex_t IT_0062 = (1 + 0.5*IT_0053)*(IT_0015 + -IT_0061);
    const ccomplex_t IT_0063 = IT_0015*IT_0052;
    const ccomplex_t IT_0064 = IT_0015*IT_0063;
    const ccomplex_t IT_0065 = IT_0001*IT_0015;
    const ccomplex_t IT_0066 = IT_0015*IT_0065;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = IT_0059 + IT_0067;
    const ccomplex_t IT_0069 = IT_0064*IT_0068;
    const ccomplex_t IT_0070 = IT_0048 + -IT_0056;
    const ccomplex_t IT_0071 = (-2)*IT_0053;
    const ccomplex_t IT_0072 = pow(m_V, -6);
    const ccomplex_t IT_0073 = IT_0059*IT_0072;
    const ccomplex_t IT_0074 = IT_0015*IT_0073;
    const ccomplex_t IT_0075 = IT_0015*IT_0074;
    const ccomplex_t IT_0076 = IT_0071 + IT_0075;
    const ccomplex_t IT_0077 = 1 + IT_0076;
    const ccomplex_t IT_0078 = IT_0015*IT_0077;
    const ccomplex_t IT_0079 = cpow(IT_0059, 2);
    const ccomplex_t IT_0080 = IT_0001*IT_0079;
    const ccomplex_t IT_0081 = (1 + 0.5*IT_0053)*(IT_0014 + -IT_0080);
    const ccomplex_t IT_0082 = IT_0052*IT_0059;
    const ccomplex_t IT_0083 = IT_0015*IT_0082;
    const ccomplex_t IT_0084 = -IT_0061;
    const ccomplex_t IT_0085 = IT_0015 + IT_0084;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = IT_0052*IT_0079;
    const ccomplex_t IT_0088 = -IT_0053 + -IT_0087;
    const ccomplex_t IT_0089 = IT_0075 + IT_0088;
    const ccomplex_t IT_0090 = 1 + IT_0089;
    const ccomplex_t IT_0091 = IT_0015*IT_0090;
    const ccomplex_t IT_0092 = (-0.666666666666667)*IT_0050;
    const ccomplex_t IT_0093 = (IT_0014 + -IT_0054)*(1 + 0.5*IT_0087);
    const ccomplex_t IT_0094 = cpow(IT_0015, 3);
    const ccomplex_t IT_0095 = IT_0052*IT_0094;
    const ccomplex_t IT_0096 = IT_0015*IT_0053;
    const ccomplex_t IT_0097 = IT_0015*IT_0087;
    const ccomplex_t IT_0098 = IT_0072*IT_0094;
    const ccomplex_t IT_0099 = IT_0059*IT_0098;
    const ccomplex_t IT_0100 = 0.666666666666667*IT_0051*IT_0055 +
       0.666666666666667*IT_0056*(IT_0055*IT_0056 + 2*IT_0058*IT_0062 + IT_0050
      *IT_0069 + -IT_0070*IT_0078) + 0.666666666666667*IT_0058*(IT_0058*IT_0081 
      + IT_0070*IT_0086 + 1.5*IT_0091*IT_0092) + 0.666666666666667*IT_0070*
      (IT_0070*IT_0093 + 1.5*IT_0092*(IT_0015 + IT_0061 + -IT_0095 + -IT_0096 + 
      -IT_0097 + IT_0099));
    const ccomplex_t IT_0101 = cabs((-0.75)*IT_0014);
    const ccomplex_t IT_0102 = cpow(IT_0101, 0.5);
    const ccomplex_t IT_0103 = ( ((m_V
) >= (2*m_V
)) ? (0.125*IT_0000*IT_0001
      *IT_0100*IT_0102
) : (0
));
    const ccomplex_t IT_0104 = cos(theta_23_R);
    const ccomplex_t IT_0105 = cpow(IT_0104, 2);
    const ccomplex_t IT_0106 = sin(theta_13_R);
    const ccomplex_t IT_0107 = cpow(IT_0106, 2);
    const ccomplex_t IT_0108 = cos(theta_13_R);
    const ccomplex_t IT_0109 = sin(theta_12_R);
    const ccomplex_t IT_0110 = cos(theta_12_R);
    const ccomplex_t IT_0111 = sin(theta_23_R);
    const ccomplex_t IT_0112 = cpow(IT_0108*IT_0109 + -IT_0106*IT_0110
      *IT_0111, 2);
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*g_f*(IT_0105*IT_0107 + 
      -IT_0112);
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = cpow(IT_0114, 2);
    const ccomplex_t IT_0116 = pow(m_e, 2);
    const ccomplex_t IT_0117 = 0.5*IT_0014 + -IT_0116;
    const ccomplex_t IT_0118 = IT_0014*IT_0117;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = IT_0016 + IT_0119;
    const ccomplex_t IT_0121 = IT_0001*IT_0120;
    const ccomplex_t IT_0122 = IT_0117 + IT_0121;
    const ccomplex_t IT_0123 = cabs(0.25*IT_0014 + -IT_0116);
    const ccomplex_t IT_0124 = cpow(IT_0123, 0.5);
    const ccomplex_t IT_0125 = ( ((m_V
) >= (2*m_e
)) ? (0.125*IT_0000*IT_0001
      *((-1.33333333333333)*IT_0115*IT_0122 + (-4)*IT_0012*(IT_0114*IT_0116 +
       0.333333333333333*IT_0012*IT_0122))*IT_0124
) : (0
));
    const ccomplex_t IT_0126 = cpow(IT_0110, 2);
    const ccomplex_t IT_0127 = cpow(IT_0111, 2);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*g_f*(IT_0105*IT_0126 + 
      -IT_0127);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = cpow(IT_0129, 2);
    const ccomplex_t IT_0131 = pow(m_mu, 2);
    const ccomplex_t IT_0132 = 0.5*IT_0014 + -IT_0131;
    const ccomplex_t IT_0133 = IT_0014*IT_0132;
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = IT_0016 + IT_0134;
    const ccomplex_t IT_0136 = IT_0001*IT_0135;
    const ccomplex_t IT_0137 = IT_0132 + IT_0136;
    const ccomplex_t IT_0138 = cabs(0.25*IT_0014 + -IT_0131);
    const ccomplex_t IT_0139 = cpow(IT_0138, 0.5);
    const ccomplex_t IT_0140 = ( ((m_V
) >= (2*m_mu
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0130*IT_0137 + (-4)*IT_0035*(IT_0129
      *IT_0131 + 0.333333333333333*IT_0035*IT_0137))*IT_0139
) : (0
));
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*g_f*IT_0104*(IT_0106*IT_0111
       + -IT_0110*(IT_0108*IT_0109 + -IT_0106*IT_0110*IT_0111));
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = cpow(IT_0142, 2);
    const ccomplex_t IT_0144 = cpow(IT_0014 + IT_0116 + -IT_0131, 2);
    const ccomplex_t IT_0145 = -IT_0116 + 0.25*IT_0001*(IT_0014 + -IT_0116 +
       IT_0131)*(IT_0014 + IT_0116 + -IT_0131) + 0.25*IT_0001*IT_0144;
    const ccomplex_t IT_0146 = 0.5*IT_0014 + 0.5*IT_0116 + (-0.5)*IT_0131;
    const ccomplex_t IT_0147 = 0.5*IT_0014 + (-0.5)*IT_0116 + 0.5*IT_0131;
    const ccomplex_t IT_0148 = IT_0146*IT_0147;
    const ccomplex_t IT_0149 = IT_0014*IT_0145;
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = IT_0148 + IT_0150;
    const ccomplex_t IT_0152 = IT_0001*IT_0151;
    const ccomplex_t IT_0153 = IT_0145 + IT_0152;
    const ccomplex_t IT_0154 = m_e*m_mu;
    const ccomplex_t IT_0155 = cabs(-IT_0116 + 0.25*IT_0001*IT_0144);
    const ccomplex_t IT_0156 = cpow(IT_0155, 0.5);
    const ccomplex_t IT_0157 = ( ((m_V
) >= (m_e + m_mu
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0143*IT_0153 + (-4)*IT_0025*
      (0.333333333333333*IT_0025*IT_0153 + IT_0142*IT_0154))*IT_0156
) : (0
));
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*g_f*m_V;
    const ccomplex_t IT_0159 = cpow(IT_0158, 2);
    const ccomplex_t IT_0160 = pow(m_phi, 2);
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = IT_0014 + IT_0161;
    const ccomplex_t IT_0163 = cpow(IT_0162, 2);
    const ccomplex_t IT_0164 = IT_0052*IT_0163;
    const ccomplex_t IT_0165 = cpow(IT_0014 + (-0.5)*IT_0160, 2);
    const ccomplex_t IT_0166 = cabs(-IT_0014 + IT_0001*IT_0165);
    const ccomplex_t IT_0167 = cpow(IT_0166, 0.5);
    const ccomplex_t IT_0168 = ( ((m_V
) >= (m_V + m_phi
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*IT_0159*(1 + 0.5*IT_0164)*IT_0167
) :
       (0
));
    const ccomplex_t IT_0169 = cpow(IT_0108, 2);
    const ccomplex_t IT_0170 = cpow(IT_0106*IT_0109 + IT_0108*IT_0110*IT_0111,
       2);
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*g_f*(IT_0105*IT_0169 + 
      -IT_0170);
    const ccomplex_t IT_0172 = 0.5*IT_0171;
    const ccomplex_t IT_0173 = cpow(IT_0172, 2);
    const ccomplex_t IT_0174 = pow(m_tau, 2);
    const ccomplex_t IT_0175 = 0.5*IT_0014 + -IT_0174;
    const ccomplex_t IT_0176 = IT_0014*IT_0175;
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = IT_0016 + IT_0177;
    const ccomplex_t IT_0179 = IT_0001*IT_0178;
    const ccomplex_t IT_0180 = IT_0175 + IT_0179;
    const ccomplex_t IT_0181 = cabs(0.25*IT_0014 + -IT_0174);
    const ccomplex_t IT_0182 = cpow(IT_0181, 0.5);
    const ccomplex_t IT_0183 = ( ((m_V
) >= (2*m_tau
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0173*IT_0180 + (-4)*IT_0045*(IT_0172
      *IT_0174 + 0.333333333333333*IT_0045*IT_0180))*IT_0182
) : (0
));
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*g_f*(IT_0105*IT_0106*IT_0108
       + (IT_0108*IT_0109 + -IT_0106*IT_0110*IT_0111)*(IT_0106*IT_0109 + IT_0108
      *IT_0110*IT_0111));
    const ccomplex_t IT_0185 = 0.5*IT_0184;
    const ccomplex_t IT_0186 = cpow(IT_0185, 2);
    const ccomplex_t IT_0187 = cpow(IT_0014 + IT_0116 + -IT_0174, 2);
    const ccomplex_t IT_0188 = -IT_0116 + 0.25*IT_0001*(IT_0014 + -IT_0116 +
       IT_0174)*(IT_0014 + IT_0116 + -IT_0174) + 0.25*IT_0001*IT_0187;
    const ccomplex_t IT_0189 = 0.5*IT_0014 + 0.5*IT_0116 + (-0.5)*IT_0174;
    const ccomplex_t IT_0190 = 0.5*IT_0014 + (-0.5)*IT_0116 + 0.5*IT_0174;
    const ccomplex_t IT_0191 = IT_0189*IT_0190;
    const ccomplex_t IT_0192 = IT_0014*IT_0188;
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0001*IT_0194;
    const ccomplex_t IT_0196 = IT_0188 + IT_0195;
    const ccomplex_t IT_0197 = m_e*m_tau;
    const ccomplex_t IT_0198 = cabs(-IT_0116 + 0.25*IT_0001*IT_0187);
    const ccomplex_t IT_0199 = cpow(IT_0198, 0.5);
    const ccomplex_t IT_0200 = ( ((m_V
) >= (m_e + m_tau
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0186*IT_0196 + (-4)*IT_0029*
      (0.333333333333333*IT_0029*IT_0196 + IT_0185*IT_0197))*IT_0199
) : (0
));
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*g_f*IT_0104*(IT_0108*IT_0111
       + IT_0110*(IT_0106*IT_0109 + IT_0108*IT_0110*IT_0111));
    const ccomplex_t IT_0202 = 0.5*IT_0201;
    const ccomplex_t IT_0203 = cpow(IT_0202, 2);
    const ccomplex_t IT_0204 = cpow(IT_0014 + IT_0131 + -IT_0174, 2);
    const ccomplex_t IT_0205 = -IT_0131 + 0.25*IT_0001*(IT_0014 + -IT_0131 +
       IT_0174)*(IT_0014 + IT_0131 + -IT_0174) + 0.25*IT_0001*IT_0204;
    const ccomplex_t IT_0206 = 0.5*IT_0014 + 0.5*IT_0131 + (-0.5)*IT_0174;
    const ccomplex_t IT_0207 = 0.5*IT_0014 + (-0.5)*IT_0131 + 0.5*IT_0174;
    const ccomplex_t IT_0208 = IT_0206*IT_0207;
    const ccomplex_t IT_0209 = IT_0014*IT_0205;
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0208 + IT_0210;
    const ccomplex_t IT_0212 = IT_0001*IT_0211;
    const ccomplex_t IT_0213 = IT_0205 + IT_0212;
    const ccomplex_t IT_0214 = m_mu*m_tau;
    const ccomplex_t IT_0215 = cabs(-IT_0131 + 0.25*IT_0001*IT_0204);
    const ccomplex_t IT_0216 = cpow(IT_0215, 0.5);
    const ccomplex_t IT_0217 = ( ((m_V
) >= (m_mu + m_tau
)) ? (0.125*IT_0000
      *IT_0001*((-1.33333333333333)*IT_0203*IT_0213 + (-4)*IT_0039*
      (0.333333333333333*IT_0039*IT_0213 + IT_0202*IT_0214))*IT_0216
) : (0
));
    const ccomplex_t IT_0218 = 0.5*IT_0048;
    const ccomplex_t IT_0219 = cpow(IT_0218, 2);
    const ccomplex_t IT_0220 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0221 = 0.5*IT_0014 + -IT_0220;
    const ccomplex_t IT_0222 = IT_0014*IT_0221;
    const ccomplex_t IT_0223 = (-0.5)*IT_0222;
    const ccomplex_t IT_0224 = IT_0016 + IT_0223;
    const ccomplex_t IT_0225 = IT_0001*IT_0224;
    const ccomplex_t IT_0226 = cabs(0.25*IT_0014 + -IT_0220);
    const ccomplex_t IT_0227 = cpow(IT_0226, 0.5);
    const ccomplex_t IT_0228 = ( ((m_V
) >= (2*m_chi_dm_1
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0219*(IT_0220 + 0.666666666666667*IT_0221 + 0.666666666666667
      *IT_0225)*IT_0227
) : (0
));
    const ccomplex_t IT_0229 = (-0.5)*IT_0048;
    const ccomplex_t IT_0230 = cpow(IT_0229, 2);
    const ccomplex_t IT_0231 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0232 = 0.5*IT_0014 + -IT_0231;
    const ccomplex_t IT_0233 = IT_0014*IT_0232;
    const ccomplex_t IT_0234 = (-0.5)*IT_0233;
    const ccomplex_t IT_0235 = IT_0016 + IT_0234;
    const ccomplex_t IT_0236 = IT_0001*IT_0235;
    const ccomplex_t IT_0237 = cabs(0.25*IT_0014 + -IT_0231);
    const ccomplex_t IT_0238 = cpow(IT_0237, 0.5);
    const ccomplex_t IT_0239 = ( ((m_V
) >= (2*m_chi_dm_2
)) ? ((-0.5)*IT_0000
      *IT_0001*IT_0230*(IT_0231 + 0.666666666666667*IT_0232 + 0.666666666666667
      *IT_0236)*IT_0238
) : (0
));
    const ccomplex_t IT_0240 = IT_0023 + 2*IT_0027 + 2*IT_0031 + IT_0037 + 2
      *IT_0041 + IT_0047 + IT_0103 + IT_0125 + IT_0140 + 2*IT_0157 + IT_0168 +
       IT_0183 + 2*IT_0200 + 2*IT_0217 + IT_0228 + IT_0239;
    return create_ccomplex_return(IT_0240);
}

