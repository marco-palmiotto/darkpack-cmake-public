#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_width_phi.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_width_phi(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t pi = param->pi;
    const creal_t g_d = param->g_d;
    const creal_t g_l = param->g_l;
    const creal_t g_u = param->g_u;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_e = param->m_e;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t g_chi = param->g_chi;
    const creal_t m_chi = param->m_chi;
    const creal_t m_phi = param->m_phi;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const ccomplex_t IT_0000 = cpow(pi, -1);
    const ccomplex_t IT_0001 = pow(m_phi, -2);
    const ccomplex_t IT_0002 = pow(m_b, 2);
    const ccomplex_t IT_0003 = pow(m_phi, 2);
    const ccomplex_t IT_0004 = -IT_0002 + 0.5*IT_0003;
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*g_d*m_b*e_em*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0009, 2);
    const ccomplex_t IT_0011 = cabs(-IT_0002 + 0.25*IT_0003);
    const ccomplex_t IT_0012 = cpow(IT_0011, 0.5);
    const ccomplex_t IT_0013 = ( ((m_phi
) >= (2*m_b
)) ? ((-0.125)*IT_0000
      *IT_0001*((-12)*IT_0002 + 12*IT_0004)*IT_0010*IT_0012
) : (0
));
    const ccomplex_t IT_0014 = pow(m_c, 2);
    const ccomplex_t IT_0015 = 0.5*IT_0003 + -IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*g_u*m_c*e_em*IT_0005*IT_0007;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0017, 2);
    const ccomplex_t IT_0019 = cabs(0.25*IT_0003 + -IT_0014);
    const ccomplex_t IT_0020 = cpow(IT_0019, 0.5);
    const ccomplex_t IT_0021 = ( ((m_phi
) >= (2*m_c
)) ? ((-0.125)*IT_0000
      *IT_0001*((-12)*IT_0014 + 12*IT_0015)*IT_0018*IT_0020
) : (0
));
    const ccomplex_t IT_0022 = pow(m_d, 2);
    const ccomplex_t IT_0023 = 0.5*IT_0003 + -IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*g_d*m_d*e_em*IT_0005*IT_0007;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0025, 2);
    const ccomplex_t IT_0027 = cabs(0.25*IT_0003 + -IT_0022);
    const ccomplex_t IT_0028 = cpow(IT_0027, 0.5);
    const ccomplex_t IT_0029 = ( ((m_phi
) >= (2*m_d
)) ? ((-0.125)*IT_0000
      *IT_0001*((-12)*IT_0022 + 12*IT_0023)*IT_0026*IT_0028
) : (0
));
    const ccomplex_t IT_0030 = pow(m_e, 2);
    const ccomplex_t IT_0031 = 0.5*IT_0003 + -IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*g_l*m_e*e_em*IT_0005*IT_0007;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = cpow(IT_0033, 2);
    const ccomplex_t IT_0035 = cabs(0.25*IT_0003 + -IT_0030);
    const ccomplex_t IT_0036 = cpow(IT_0035, 0.5);
    const ccomplex_t IT_0037 = ( ((m_phi
) >= (2*m_e
)) ? (0.5*IT_0000*IT_0001
      *(IT_0030 + -IT_0031)*IT_0034*IT_0036
) : (0
));
    const ccomplex_t IT_0038 = pow(m_s, 2);
    const ccomplex_t IT_0039 = 0.5*IT_0003 + -IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*g_d*m_s*e_em*IT_0005*IT_0007;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = cpow(IT_0041, 2);
    const ccomplex_t IT_0043 = cabs(0.25*IT_0003 + -IT_0038);
    const ccomplex_t IT_0044 = cpow(IT_0043, 0.5);
    const ccomplex_t IT_0045 = ( ((m_phi
) >= (2*m_s
)) ? ((-0.125)*IT_0000
      *IT_0001*((-12)*IT_0038 + 12*IT_0039)*IT_0042*IT_0044
) : (0
));
    const ccomplex_t IT_0046 = pow(m_t, 2);
    const ccomplex_t IT_0047 = 0.5*IT_0003 + -IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*g_u*m_t*e_em*IT_0005*IT_0007;
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = cpow(IT_0049, 2);
    const ccomplex_t IT_0051 = cabs(0.25*IT_0003 + -IT_0046);
    const ccomplex_t IT_0052 = cpow(IT_0051, 0.5);
    const ccomplex_t IT_0053 = ( ((m_phi
) >= (2*m_t
)) ? ((-0.125)*IT_0000
      *IT_0001*((-12)*IT_0046 + 12*IT_0047)*IT_0050*IT_0052
) : (0
));
    const ccomplex_t IT_0054 = pow(m_u, 2);
    const ccomplex_t IT_0055 = 0.5*IT_0003 + -IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*g_u*m_u*e_em*IT_0005*IT_0007;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = cpow(IT_0057, 2);
    const ccomplex_t IT_0059 = cabs(0.25*IT_0003 + -IT_0054);
    const ccomplex_t IT_0060 = cpow(IT_0059, 0.5);
    const ccomplex_t IT_0061 = ( ((m_phi
) >= (2*m_u
)) ? ((-0.125)*IT_0000
      *IT_0001*((-12)*IT_0054 + 12*IT_0055)*IT_0058*IT_0060
) : (0
));
    const ccomplex_t IT_0062 = pow(m_mu, 2);
    const ccomplex_t IT_0063 = 0.5*IT_0003 + -IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*g_l*e_em*m_mu*IT_0005*IT_0007;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = cpow(IT_0065, 2);
    const ccomplex_t IT_0067 = cabs(0.25*IT_0003 + -IT_0062);
    const ccomplex_t IT_0068 = cpow(IT_0067, 0.5);
    const ccomplex_t IT_0069 = ( ((m_phi
) >= (2*m_mu
)) ? (0.5*IT_0000
      *IT_0001*(IT_0062 + -IT_0063)*IT_0066*IT_0068
) : (0
));
    const ccomplex_t IT_0070 = pow(m_chi, 2);
    const ccomplex_t IT_0071 = 0.5*IT_0003 + -IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*g_chi;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = cpow(IT_0073, 2);
    const ccomplex_t IT_0075 = cabs(0.25*IT_0003 + -IT_0070);
    const ccomplex_t IT_0076 = cpow(IT_0075, 0.5);
    const ccomplex_t IT_0077 = ( ((m_phi
) >= (2*m_chi
)) ? (0.5*IT_0000
      *IT_0001*(IT_0070 + -IT_0071)*IT_0074*IT_0076
) : (0
));
    const ccomplex_t IT_0078 = pow(m_tau, 2);
    const ccomplex_t IT_0079 = 0.5*IT_0003 + -IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*g_l*e_em*m_tau*IT_0005
      *IT_0007;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = cpow(IT_0081, 2);
    const ccomplex_t IT_0083 = cabs(0.25*IT_0003 + -IT_0078);
    const ccomplex_t IT_0084 = cpow(IT_0083, 0.5);
    const ccomplex_t IT_0085 = ( ((m_phi
) >= (2*m_tau
)) ? (0.5*IT_0000
      *IT_0001*(IT_0078 + -IT_0079)*IT_0082*IT_0084
) : (0
));
    return create_ccomplex_return(IT_0013 + IT_0021 + IT_0029 + IT_0037 +
       IT_0045 + IT_0053 + IT_0061 + IT_0069 + IT_0077 + IT_0085);
}

