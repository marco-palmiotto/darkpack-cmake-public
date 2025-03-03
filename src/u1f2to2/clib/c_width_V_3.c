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
    const creal_t m_e = param->m_e;
    const creal_t m_mu = param->m_mu;
    const creal_t m_V_3 = param->m_V_3;
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
    const ccomplex_t IT_0001 = pow(m_V_3, -2);
    const ccomplex_t IT_0002 = pow(m_V_3, 2);
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
    const ccomplex_t IT_0023 = ( ((m_V_3
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0020*IT_0022
) : (0
));
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*g_f*IT_0009*(IT_0011*IT_0016
       + -IT_0015*(IT_0013*IT_0014 + -IT_0011*IT_0015*IT_0016));
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0025, 2);
    const ccomplex_t IT_0027 = ( ((m_V_3
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
    const ccomplex_t IT_0031 = ( ((m_V_3
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
    const ccomplex_t IT_0037 = ( ((m_V_3
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0022*IT_0036
) : (0
));
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*g_f*IT_0009*(IT_0013*IT_0016
       + IT_0015*(IT_0011*IT_0014 + IT_0013*IT_0015*IT_0016));
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0039, 2);
    const ccomplex_t IT_0041 = ( ((m_V_3
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
    const ccomplex_t IT_0047 = ( ((m_V_3
) >= (0
)) ? ((-0.166666666666667)
      *IT_0000*IT_0001*(IT_0003 + IT_0008)*IT_0022*IT_0046
) : (0
));
    const ccomplex_t IT_0048 = cos(theta_23_R);
    const ccomplex_t IT_0049 = cpow(IT_0048, 2);
    const ccomplex_t IT_0050 = sin(theta_13_R);
    const ccomplex_t IT_0051 = cpow(IT_0050, 2);
    const ccomplex_t IT_0052 = cos(theta_13_R);
    const ccomplex_t IT_0053 = sin(theta_12_R);
    const ccomplex_t IT_0054 = cos(theta_12_R);
    const ccomplex_t IT_0055 = sin(theta_23_R);
    const ccomplex_t IT_0056 = cpow(IT_0052*IT_0053 + -IT_0050*IT_0054
      *IT_0055, 2);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*g_f*(IT_0049*IT_0051 + 
      -IT_0056);
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = cpow(IT_0058, 2);
    const ccomplex_t IT_0060 = pow(m_e, 2);
    const ccomplex_t IT_0061 = 0.5*IT_0002 + -IT_0060;
    const ccomplex_t IT_0062 = IT_0002*IT_0061;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0004 + IT_0063;
    const ccomplex_t IT_0065 = IT_0001*IT_0064;
    const ccomplex_t IT_0066 = IT_0061 + IT_0065;
    const ccomplex_t IT_0067 = (-1.33333333333333)*IT_0066;
    const ccomplex_t IT_0068 = cabs(0.25*IT_0002 + -IT_0060);
    const ccomplex_t IT_0069 = cpow(IT_0068, 0.5);
    const ccomplex_t IT_0070 = ( ((m_V_3
) >= (2*m_e
)) ? (0.125*IT_0000
      *IT_0001*(IT_0059*IT_0067 + (-4)*IT_0019*(IT_0058*IT_0060 + (-0.25)
      *IT_0019*IT_0067))*IT_0069
) : (0
));
    const ccomplex_t IT_0071 = cpow(IT_0054, 2);
    const ccomplex_t IT_0072 = cpow(IT_0055, 2);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*g_f*(IT_0049*IT_0071 + 
      -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = cpow(IT_0074, 2);
    const ccomplex_t IT_0076 = pow(m_mu, 2);
    const ccomplex_t IT_0077 = 0.5*IT_0002 + -IT_0076;
    const ccomplex_t IT_0078 = IT_0002*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0004 + IT_0079;
    const ccomplex_t IT_0081 = IT_0001*IT_0080;
    const ccomplex_t IT_0082 = IT_0077 + IT_0081;
    const ccomplex_t IT_0083 = (-1.33333333333333)*IT_0082;
    const ccomplex_t IT_0084 = cabs(0.25*IT_0002 + -IT_0076);
    const ccomplex_t IT_0085 = cpow(IT_0084, 0.5);
    const ccomplex_t IT_0086 = ( ((m_V_3
) >= (2*m_mu
)) ? (0.125*IT_0000
      *IT_0001*(IT_0075*IT_0083 + (-4)*IT_0035*(IT_0074*IT_0076 + (-0.25)
      *IT_0035*IT_0083))*IT_0085
) : (0
));
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*g_f*IT_0048*(IT_0050*IT_0055
       + -IT_0054*(IT_0052*IT_0053 + -IT_0050*IT_0054*IT_0055));
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = cpow(IT_0088, 2);
    const ccomplex_t IT_0090 = cpow(IT_0002 + IT_0060 + -IT_0076, 2);
    const ccomplex_t IT_0091 = -IT_0060 + 0.25*IT_0001*(IT_0002 + -IT_0060 +
       IT_0076)*(IT_0002 + IT_0060 + -IT_0076) + 0.25*IT_0001*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0002 + 0.5*IT_0060 + (-0.5)*IT_0076;
    const ccomplex_t IT_0093 = 0.5*IT_0002 + (-0.5)*IT_0060 + 0.5*IT_0076;
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = IT_0002*IT_0091;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = IT_0094 + IT_0096;
    const ccomplex_t IT_0098 = IT_0001*IT_0097;
    const ccomplex_t IT_0099 = IT_0091 + IT_0098;
    const ccomplex_t IT_0100 = (-1.33333333333333)*IT_0099;
    const ccomplex_t IT_0101 = m_e*m_mu;
    const ccomplex_t IT_0102 = cabs(-IT_0060 + 0.25*IT_0001*IT_0090);
    const ccomplex_t IT_0103 = cpow(IT_0102, 0.5);
    const ccomplex_t IT_0104 = ( ((m_V_3
) >= (m_e + m_mu
)) ? (0.125*IT_0000
      *IT_0001*(IT_0089*IT_0100 + (-4)*IT_0025*((-0.25)*IT_0025*IT_0100 +
       IT_0088*IT_0101))*IT_0103
) : (0
));
    const ccomplex_t IT_0105 = pow(m_V_3, -4);
    const ccomplex_t IT_0106 = pow(m_phi, 2);
    const ccomplex_t IT_0107 = (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = IT_0002 + IT_0107;
    const ccomplex_t IT_0109 = cpow(IT_0108, 2);
    const ccomplex_t IT_0110 = IT_0105*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*g_f*m_V_3;
    const ccomplex_t IT_0112 = cpow(IT_0111, 2);
    const ccomplex_t IT_0113 = cpow(IT_0002 + (-0.5)*IT_0106, 2);
    const ccomplex_t IT_0114 = cabs(-IT_0002 + IT_0001*IT_0113);
    const ccomplex_t IT_0115 = cpow(IT_0114, 0.5);
    const ccomplex_t IT_0116 = ( ((m_V_3
) >= (m_V_3 + m_phi
)) ? ((
      -0.0833333333333333)*IT_0000*IT_0001*(1 + 0.5*IT_0110)*IT_0112*IT_0115
) :
       (0
));
    const ccomplex_t IT_0117 = cpow(IT_0052, 2);
    const ccomplex_t IT_0118 = cpow(IT_0050*IT_0053 + IT_0052*IT_0054*IT_0055,
       2);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*g_f*(IT_0049*IT_0117 + 
      -IT_0118);
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = cpow(IT_0120, 2);
    const ccomplex_t IT_0122 = pow(m_tau, 2);
    const ccomplex_t IT_0123 = 0.5*IT_0002 + -IT_0122;
    const ccomplex_t IT_0124 = IT_0002*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0004 + IT_0125;
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = IT_0123 + IT_0127;
    const ccomplex_t IT_0129 = (-1.33333333333333)*IT_0128;
    const ccomplex_t IT_0130 = cabs(0.25*IT_0002 + -IT_0122);
    const ccomplex_t IT_0131 = cpow(IT_0130, 0.5);
    const ccomplex_t IT_0132 = ( ((m_V_3
) >= (2*m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0121*IT_0129 + (-4)*IT_0045*(IT_0120*IT_0122 + (-0.25)
      *IT_0045*IT_0129))*IT_0131
) : (0
));
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*g_f*(IT_0049*IT_0050*IT_0052
       + (IT_0052*IT_0053 + -IT_0050*IT_0054*IT_0055)*(IT_0050*IT_0053 + IT_0052
      *IT_0054*IT_0055));
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = cpow(IT_0134, 2);
    const ccomplex_t IT_0136 = cpow(IT_0002 + IT_0060 + -IT_0122, 2);
    const ccomplex_t IT_0137 = -IT_0060 + 0.25*IT_0001*(IT_0002 + -IT_0060 +
       IT_0122)*(IT_0002 + IT_0060 + -IT_0122) + 0.25*IT_0001*IT_0136;
    const ccomplex_t IT_0138 = 0.5*IT_0002 + 0.5*IT_0060 + (-0.5)*IT_0122;
    const ccomplex_t IT_0139 = 0.5*IT_0002 + (-0.5)*IT_0060 + 0.5*IT_0122;
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = IT_0002*IT_0137;
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = IT_0140 + IT_0142;
    const ccomplex_t IT_0144 = IT_0001*IT_0143;
    const ccomplex_t IT_0145 = IT_0137 + IT_0144;
    const ccomplex_t IT_0146 = (-1.33333333333333)*IT_0145;
    const ccomplex_t IT_0147 = m_e*m_tau;
    const ccomplex_t IT_0148 = cabs(-IT_0060 + 0.25*IT_0001*IT_0136);
    const ccomplex_t IT_0149 = cpow(IT_0148, 0.5);
    const ccomplex_t IT_0150 = ( ((m_V_3
) >= (m_e + m_tau
)) ? (0.125*IT_0000
      *IT_0001*(IT_0135*IT_0146 + (-4)*IT_0029*((-0.25)*IT_0029*IT_0146 +
       IT_0134*IT_0147))*IT_0149
) : (0
));
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*g_f*IT_0048*(IT_0052*IT_0055
       + IT_0054*(IT_0050*IT_0053 + IT_0052*IT_0054*IT_0055));
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = cpow(IT_0152, 2);
    const ccomplex_t IT_0154 = cpow(IT_0002 + IT_0076 + -IT_0122, 2);
    const ccomplex_t IT_0155 = -IT_0076 + 0.25*IT_0001*(IT_0002 + -IT_0076 +
       IT_0122)*(IT_0002 + IT_0076 + -IT_0122) + 0.25*IT_0001*IT_0154;
    const ccomplex_t IT_0156 = 0.5*IT_0002 + 0.5*IT_0076 + (-0.5)*IT_0122;
    const ccomplex_t IT_0157 = 0.5*IT_0002 + (-0.5)*IT_0076 + 0.5*IT_0122;
    const ccomplex_t IT_0158 = IT_0156*IT_0157;
    const ccomplex_t IT_0159 = IT_0002*IT_0155;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = IT_0158 + IT_0160;
    const ccomplex_t IT_0162 = IT_0001*IT_0161;
    const ccomplex_t IT_0163 = IT_0155 + IT_0162;
    const ccomplex_t IT_0164 = (-1.33333333333333)*IT_0163;
    const ccomplex_t IT_0165 = m_mu*m_tau;
    const ccomplex_t IT_0166 = cabs(-IT_0076 + 0.25*IT_0001*IT_0154);
    const ccomplex_t IT_0167 = cpow(IT_0166, 0.5);
    const ccomplex_t IT_0168 = ( ((m_V_3
) >= (m_mu + m_tau
)) ? (0.125
      *IT_0000*IT_0001*(IT_0153*IT_0164 + (-4)*IT_0039*((-0.25)*IT_0039*IT_0164 
      + IT_0152*IT_0165))*IT_0167
) : (0
));
    const ccomplex_t IT_0169 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0170 = 0.5*IT_0002 + -IT_0169;
    const ccomplex_t IT_0171 = IT_0002*IT_0170;
    const ccomplex_t IT_0172 = (-0.5)*IT_0171;
    const ccomplex_t IT_0173 = IT_0004 + IT_0172;
    const ccomplex_t IT_0174 = IT_0001*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0176 = 0.5*IT_0175;
    const ccomplex_t IT_0177 = cpow(IT_0176, 2);
    const ccomplex_t IT_0178 = cabs(0.25*IT_0002 + -IT_0169);
    const ccomplex_t IT_0179 = cpow(IT_0178, 0.5);
    const ccomplex_t IT_0180 = ( ((m_V_3
) >= (2*m_chi_dm_1
)) ? ((-0.5)
      *IT_0000*IT_0001*(IT_0169 + 0.666666666666667*IT_0170 + 0.666666666666667
      *IT_0174)*IT_0177*IT_0179
) : (0
));
    const ccomplex_t IT_0181 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0182 = 0.5*IT_0002 + -IT_0181;
    const ccomplex_t IT_0183 = IT_0002*IT_0182;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = IT_0004 + IT_0184;
    const ccomplex_t IT_0186 = IT_0001*IT_0185;
    const ccomplex_t IT_0187 = (-0.5)*IT_0175;
    const ccomplex_t IT_0188 = cpow(IT_0187, 2);
    const ccomplex_t IT_0189 = cabs(0.25*IT_0002 + -IT_0181);
    const ccomplex_t IT_0190 = cpow(IT_0189, 0.5);
    const ccomplex_t IT_0191 = ( ((m_V_3
) >= (2*m_chi_dm_2
)) ? ((-0.5)
      *IT_0000*IT_0001*(IT_0181 + 0.666666666666667*IT_0182 + 0.666666666666667
      *IT_0186)*IT_0188*IT_0190
) : (0
));
    const ccomplex_t IT_0192 = IT_0023 + 2*IT_0027 + 2*IT_0031 + IT_0037 + 2
      *IT_0041 + IT_0047 + IT_0070 + IT_0086 + 2*IT_0104 + IT_0116 + IT_0132 + 2
      *IT_0150 + 2*IT_0168 + IT_0180 + IT_0191;
    return create_ccomplex_return(IT_0192);
}

