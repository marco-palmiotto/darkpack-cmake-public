#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t V_tb = param->V_tb;
    const creal_t e_em = param->e_em;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_t, -2);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0006 = (-0.666666666666667)*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = pow(m_t, 2);
    const ccomplex_t IT_0009 = pow(m_Z, 2);
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0008 + IT_0010;
    const ccomplex_t IT_0012 = pow(m_Z, -2);
    const ccomplex_t IT_0013 = IT_0009*IT_0011;
    const ccomplex_t IT_0014 = 0.5*IT_0009;
    const ccomplex_t IT_0015 = pow(m_Z, 4);
    const ccomplex_t IT_0016 = -IT_0009 + 0.5*IT_0001*IT_0009*(IT_0008 + -1./2
      *IT_0009) + 0.25*IT_0001*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (-2)*IT_0017;
    const ccomplex_t IT_0019 = IT_0013 + IT_0018;
    const ccomplex_t IT_0020 = IT_0012*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0011 + IT_0021;
    const ccomplex_t IT_0023 = (-2)*IT_0022;
    const ccomplex_t IT_0024 = IT_0003*IT_0004;
    const ccomplex_t IT_0025 = e_em*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0004, -1);
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = e_em*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + (-3)*IT_0028);
    const ccomplex_t IT_0030 = (-0.166666666666667)*IT_0029;
    const ccomplex_t IT_0031 = cabs(-IT_0009 + 0.25*IT_0001*IT_0015);
    const ccomplex_t IT_0032 = cpow(IT_0031, 0.5);
    const ccomplex_t IT_0033 = ( ((m_t
) >= (m_Z + m_t
)) ? (0.125*IT_0000
      *IT_0001*(IT_0007*IT_0023 + 6*IT_0030*(IT_0006*IT_0008 + 0.166666666666667
      *IT_0023*IT_0030))*IT_0032
) : (0
));
    const ccomplex_t IT_0034 = pow(m_h, 2);
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0008 + IT_0035;
    const ccomplex_t IT_0037 = pow(m_W, -1);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_t*e_em*IT_0026*IT_0037;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0039, 2);
    const ccomplex_t IT_0041 = pow(m_h, 4);
    const ccomplex_t IT_0042 = cabs(-IT_0034 + 0.25*IT_0001*IT_0041);
    const ccomplex_t IT_0043 = cpow(IT_0042, 0.5);
    const ccomplex_t IT_0044 = ( ((m_t
) >= (m_h + m_t
)) ? ((-0.25)*IT_0000
      *IT_0001*(IT_0008 + IT_0036)*IT_0040*IT_0043
) : (0
));
    const ccomplex_t IT_0045 = pow(m_W, 2);
    const ccomplex_t IT_0046 = pow(m_b, 2);
    const ccomplex_t IT_0047 = 0.5*IT_0008 + (-0.5)*IT_0045 + 0.5*IT_0046;
    const ccomplex_t IT_0048 = pow(m_W, -2);
    const ccomplex_t IT_0049 = IT_0045*IT_0047;
    const ccomplex_t IT_0050 = 0.5*IT_0008 + 0.5*IT_0045 + (-0.5)*IT_0046;
    const ccomplex_t IT_0051 = cpow(IT_0008 + IT_0045 + -IT_0046, 2);
    const ccomplex_t IT_0052 = -IT_0045 + 0.25*IT_0001*(IT_0008 + -IT_0045 +
       IT_0046)*(IT_0008 + IT_0045 + -IT_0046) + 0.25*IT_0001*IT_0051;
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = (-2)*IT_0053;
    const ccomplex_t IT_0055 = IT_0049 + IT_0054;
    const ccomplex_t IT_0056 = IT_0048*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0026;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = cpow(IT_0058, 2);
    const ccomplex_t IT_0060 = cabs(-IT_0045 + 0.25*IT_0001*IT_0051);
    const ccomplex_t IT_0061 = cpow(IT_0060, 0.5);
    const ccomplex_t IT_0062 = ( ((m_t
) >= (m_W + m_b
)) ? ((-0.25)*IT_0000
      *IT_0001*(IT_0047 + (-0.5)*IT_0056)*IT_0059*IT_0061
) : (0
));
    const ccomplex_t IT_0063 = pow(m_d, 2);
    const ccomplex_t IT_0064 = 0.5*IT_0008 + (-0.5)*IT_0045 + 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0045*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0008 + 0.5*IT_0045 + (-0.5)*IT_0063;
    const ccomplex_t IT_0067 = cpow(IT_0008 + IT_0045 + -IT_0063, 2);
    const ccomplex_t IT_0068 = -IT_0045 + 0.25*IT_0001*(IT_0008 + -IT_0045 +
       IT_0063)*(IT_0008 + IT_0045 + -IT_0063) + 0.25*IT_0001*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = (-2)*IT_0069;
    const ccomplex_t IT_0071 = IT_0065 + IT_0070;
    const ccomplex_t IT_0072 = IT_0048*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0026;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = cabs(-IT_0045 + 0.25*IT_0001*IT_0067);
    const ccomplex_t IT_0076 = cpow(IT_0075, 0.5);
    const ccomplex_t IT_0077 = ( ((m_t
) >= (m_W + m_d
)) ? (0.25*IT_0000
      *IT_0001*(IT_0064 + (-0.5)*IT_0072)*IT_0074*conj(IT_0074)*IT_0076
) : (0
));
    const ccomplex_t IT_0078 = pow(m_s, 2);
    const ccomplex_t IT_0079 = 0.5*IT_0008 + (-0.5)*IT_0045 + 0.5*IT_0078;
    const ccomplex_t IT_0080 = IT_0045*IT_0079;
    const ccomplex_t IT_0081 = 0.5*IT_0008 + 0.5*IT_0045 + (-0.5)*IT_0078;
    const ccomplex_t IT_0082 = cpow(IT_0008 + IT_0045 + -IT_0078, 2);
    const ccomplex_t IT_0083 = -IT_0045 + 0.25*IT_0001*(IT_0008 + -IT_0045 +
       IT_0078)*(IT_0008 + IT_0045 + -IT_0078) + 0.25*IT_0001*IT_0082;
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = IT_0080 + IT_0085;
    const ccomplex_t IT_0087 = IT_0048*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*conj(V_ts)
      *e_em*IT_0026;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = cabs(-IT_0045 + 0.25*IT_0001*IT_0082);
    const ccomplex_t IT_0091 = cpow(IT_0090, 0.5);
    const ccomplex_t IT_0092 = ( ((m_t
) >= (m_W + m_s
)) ? (0.25*IT_0000
      *IT_0001*(IT_0079 + (-0.5)*IT_0087)*IT_0089*conj(IT_0089)*IT_0091
) : (0
));
    return create_ccomplex_return(IT_0033 + IT_0044 + IT_0062 + IT_0077 +
       IT_0092);
}

