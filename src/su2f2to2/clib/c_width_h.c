#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_h(
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
    const ccomplex_t IT_0001 = pow(m_h, -2);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = pow(m_h, 2);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-1.5)*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = cabs((-0.75)*IT_0003);
    const ccomplex_t IT_0010 = cpow(IT_0009, 0.5);
    const ccomplex_t IT_0011 = ( ((m_h
) >= (2*m_h
)) ? ((-0.0625)*IT_0000
      *IT_0001*IT_0008*IT_0010
) : (0
));
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*m_W*e_em*IT_0005;
    const ccomplex_t IT_0013 = cpow(IT_0012, 2);
    const ccomplex_t IT_0014 = pow(m_W, -4);
    const ccomplex_t IT_0015 = pow(m_W, 2);
    const ccomplex_t IT_0016 = 0.5*IT_0003 + -IT_0015;
    const ccomplex_t IT_0017 = cpow(IT_0016, 2);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = cabs(0.25*IT_0003 + -IT_0015);
    const ccomplex_t IT_0020 = cpow(IT_0019, 0.5);
    const ccomplex_t IT_0021 = ( ((m_h
) >= (2*m_W
)) ? ((-0.25)*IT_0000
      *IT_0001*IT_0013*(1 + 0.5*IT_0018)*IT_0020
) : (0
));
    const ccomplex_t IT_0022 = m_W*e_em;
    const ccomplex_t IT_0023 = cos(theta_W);
    const ccomplex_t IT_0024 = cpow(IT_0023, -2);
    const ccomplex_t IT_0025 = cpow(IT_0004, 3);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0023, 2);
    const ccomplex_t IT_0028 = IT_0005*IT_0027;
    const ccomplex_t IT_0029 = IT_0026 + IT_0028;
    const ccomplex_t IT_0030 = 2*IT_0004;
    const ccomplex_t IT_0031 = IT_0029 + IT_0030;
    const ccomplex_t IT_0032 = IT_0022*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = 2*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0035, 2);
    const ccomplex_t IT_0037 = pow(m_Z, -4);
    const ccomplex_t IT_0038 = pow(m_Z, 2);
    const ccomplex_t IT_0039 = 0.5*IT_0003 + -IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0039, 2);
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = cabs(0.25*IT_0003 + -IT_0038);
    const ccomplex_t IT_0043 = cpow(IT_0042, 0.5);
    const ccomplex_t IT_0044 = ( ((m_h
) >= (2*m_Z
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0036*(1 + 0.5*IT_0041)*IT_0043
) : (0
));
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*m_b*e_em*IT_0002*IT_0005;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = cpow(IT_0046, 2);
    const ccomplex_t IT_0048 = pow(m_b, 2);
    const ccomplex_t IT_0049 = 0.5*IT_0003 + -IT_0048;
    const ccomplex_t IT_0050 = cabs(0.25*IT_0003 + -IT_0048);
    const ccomplex_t IT_0051 = cpow(IT_0050, 0.5);
    const ccomplex_t IT_0052 = ( ((m_h
) >= (2*m_b
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0047*((-12)*IT_0048 + 12*IT_0049)*IT_0051
) : (0
));
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_c*e_em*IT_0002*IT_0005;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = cpow(IT_0054, 2);
    const ccomplex_t IT_0056 = pow(m_c, 2);
    const ccomplex_t IT_0057 = 0.5*IT_0003 + -IT_0056;
    const ccomplex_t IT_0058 = cabs(0.25*IT_0003 + -IT_0056);
    const ccomplex_t IT_0059 = cpow(IT_0058, 0.5);
    const ccomplex_t IT_0060 = ( ((m_h
) >= (2*m_c
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0055*((-12)*IT_0056 + 12*IT_0057)*IT_0059
) : (0
));
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*m_d*e_em*IT_0002*IT_0005;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = cpow(IT_0062, 2);
    const ccomplex_t IT_0064 = pow(m_d, 2);
    const ccomplex_t IT_0065 = 0.5*IT_0003 + -IT_0064;
    const ccomplex_t IT_0066 = cabs(0.25*IT_0003 + -IT_0064);
    const ccomplex_t IT_0067 = cpow(IT_0066, 0.5);
    const ccomplex_t IT_0068 = ( ((m_h
) >= (2*m_d
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0063*((-12)*IT_0064 + 12*IT_0065)*IT_0067
) : (0
));
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*m_e*e_em*IT_0002*IT_0005;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = cpow(IT_0070, 2);
    const ccomplex_t IT_0072 = pow(m_e, 2);
    const ccomplex_t IT_0073 = 0.5*IT_0003 + -IT_0072;
    const ccomplex_t IT_0074 = cabs(0.25*IT_0003 + -IT_0072);
    const ccomplex_t IT_0075 = cpow(IT_0074, 0.5);
    const ccomplex_t IT_0076 = ( ((m_h
) >= (2*m_e
)) ? (0.5*IT_0000*IT_0001
      *IT_0071*(IT_0072 + -IT_0073)*IT_0075
) : (0
));
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_s*e_em*IT_0002*IT_0005;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = cpow(IT_0078, 2);
    const ccomplex_t IT_0080 = pow(m_s, 2);
    const ccomplex_t IT_0081 = 0.5*IT_0003 + -IT_0080;
    const ccomplex_t IT_0082 = cabs(0.25*IT_0003 + -IT_0080);
    const ccomplex_t IT_0083 = cpow(IT_0082, 0.5);
    const ccomplex_t IT_0084 = ( ((m_h
) >= (2*m_s
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0079*((-12)*IT_0080 + 12*IT_0081)*IT_0083
) : (0
));
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_t*e_em*IT_0002*IT_0005;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = cpow(IT_0086, 2);
    const ccomplex_t IT_0088 = pow(m_t, 2);
    const ccomplex_t IT_0089 = 0.5*IT_0003 + -IT_0088;
    const ccomplex_t IT_0090 = cabs(0.25*IT_0003 + -IT_0088);
    const ccomplex_t IT_0091 = cpow(IT_0090, 0.5);
    const ccomplex_t IT_0092 = ( ((m_h
) >= (2*m_t
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0087*((-12)*IT_0088 + 12*IT_0089)*IT_0091
) : (0
));
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_u*e_em*IT_0002*IT_0005;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = cpow(IT_0094, 2);
    const ccomplex_t IT_0096 = pow(m_u, 2);
    const ccomplex_t IT_0097 = 0.5*IT_0003 + -IT_0096;
    const ccomplex_t IT_0098 = cabs(0.25*IT_0003 + -IT_0096);
    const ccomplex_t IT_0099 = cpow(IT_0098, 0.5);
    const ccomplex_t IT_0100 = ( ((m_h
) >= (2*m_u
)) ? ((-0.125)*IT_0000
      *IT_0001*IT_0095*((-12)*IT_0096 + 12*IT_0097)*IT_0099
) : (0
));
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*e_em*m_mu*IT_0002*IT_0005;
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = cpow(IT_0102, 2);
    const ccomplex_t IT_0104 = pow(m_mu, 2);
    const ccomplex_t IT_0105 = 0.5*IT_0003 + -IT_0104;
    const ccomplex_t IT_0106 = cabs(0.25*IT_0003 + -IT_0104);
    const ccomplex_t IT_0107 = cpow(IT_0106, 0.5);
    const ccomplex_t IT_0108 = ( ((m_h
) >= (2*m_mu
)) ? (0.5*IT_0000*IT_0001
      *IT_0103*(IT_0104 + -IT_0105)*IT_0107
) : (0
));
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*e_em*m_tau*IT_0002*IT_0005;
    const ccomplex_t IT_0110 = (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = cpow(IT_0110, 2);
    const ccomplex_t IT_0112 = pow(m_tau, 2);
    const ccomplex_t IT_0113 = 0.5*IT_0003 + -IT_0112;
    const ccomplex_t IT_0114 = cabs(0.25*IT_0003 + -IT_0112);
    const ccomplex_t IT_0115 = cpow(IT_0114, 0.5);
    const ccomplex_t IT_0116 = ( ((m_h
) >= (2*m_tau
)) ? (0.5*IT_0000*IT_0001
      *IT_0111*(IT_0112 + -IT_0113)*IT_0115
) : (0
));
    return create_ccomplex_return(IT_0011 + IT_0021 + IT_0044 + IT_0052 +
       IT_0060 + IT_0068 + IT_0076 + IT_0084 + IT_0092 + IT_0100 + IT_0108 +
       IT_0116);
}

