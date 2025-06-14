#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_su_R_to_u_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_su_R_to_u_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_Z = param->m_Z;
    const creal_t m_u = param->m_u;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.33333333333333*IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0006 = pow(m_u, 2);
    const ccomplex_t IT_0007 = pow(m_sG, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0006 + IT_0007 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = IT_0004*IT_0010;
    const ccomplex_t IT_0012 = pow(m_Z, -2);
    const ccomplex_t IT_0013 = pow(s_14, 2);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0007 + IT_0015;
    const ccomplex_t IT_0017 = s_13*IT_0016;
    const ccomplex_t IT_0018 = pow(m_su_R, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + IT_0006 + -IT_0007 + -IT_0018 
      + -reg_prop, -1);
    const ccomplex_t IT_0020 = (-0.666666666666667)*IT_0003;
    const ccomplex_t IT_0021 = IT_0005*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = -IT_0011 + (-2)*IT_0023;
    const ccomplex_t IT_0025 = s_14*s_34*IT_0012;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = s_13 + IT_0026;
    const ccomplex_t IT_0028 = s_13*IT_0027;
    const ccomplex_t IT_0029 = IT_0001*IT_0002;
    const ccomplex_t IT_0030 = e_em*IT_0029;
    const ccomplex_t IT_0031 = cpow(IT_0002, -1);
    const ccomplex_t IT_0032 = IT_0000*IT_0031;
    const ccomplex_t IT_0033 = e_em*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + (-3)*IT_0033);
    const ccomplex_t IT_0035 = (-0.166666666666667)*IT_0034;
    const ccomplex_t IT_0036 = IT_0005*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0038 = IT_0019*IT_0036*IT_0037;
    const ccomplex_t IT_0039 = m_u*IT_0023;
    const ccomplex_t IT_0040 = -IT_0038 + IT_0039;
    const ccomplex_t IT_0041 = pow(m_Z, 2);
    const ccomplex_t IT_0042 = s_13*IT_0041;
    const ccomplex_t IT_0043 = s_14*s_34;
    const ccomplex_t IT_0044 = (-2)*IT_0043;
    const ccomplex_t IT_0045 = IT_0042 + IT_0044;
    const ccomplex_t IT_0046 = IT_0012*IT_0045;
    const ccomplex_t IT_0047 = m_u*IT_0016;
    const ccomplex_t IT_0048 = (-8)*conj(IT_0011);
    const ccomplex_t IT_0049 = 0.125*IT_0048;
    const ccomplex_t IT_0050 = IT_0011 + IT_0049;
    const ccomplex_t IT_0051 = m_u*IT_0027;
    const ccomplex_t IT_0052 = 8*IT_0024;
    const ccomplex_t IT_0053 = (-8)*conj(IT_0024);
    const ccomplex_t IT_0054 = IT_0052 + IT_0053;
    const ccomplex_t IT_0055 = -IT_0023;
    const ccomplex_t IT_0056 = pow(m_Z, 4);
    const ccomplex_t IT_0057 = s_13*IT_0056;
    const ccomplex_t IT_0058 = IT_0012*IT_0057;
    const ccomplex_t IT_0059 = s_14*IT_0041;
    const ccomplex_t IT_0060 = IT_0012*IT_0059;
    const ccomplex_t IT_0061 = (-0.25)*IT_0060;
    const ccomplex_t IT_0062 = s_14 + IT_0061;
    const ccomplex_t IT_0063 = m_u*IT_0062;
    const ccomplex_t IT_0064 = (-4)*IT_0063;
    const ccomplex_t IT_0065 = s_13*s_14;
    const ccomplex_t IT_0066 = s_34*IT_0007;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = IT_0065 + IT_0067;
    const ccomplex_t IT_0069 = s_14*IT_0006;
    const ccomplex_t IT_0070 = s_13*s_34;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0069 + IT_0071;
    const ccomplex_t IT_0073 = pow(s_34, 2);
    const ccomplex_t IT_0074 = IT_0012*IT_0073;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0006 + IT_0075;
    const ccomplex_t IT_0077 = s_13*IT_0076;
    const ccomplex_t IT_0078 = (-8)*IT_0011*(conj(IT_0011)*IT_0017 + conj
      (IT_0024)*IT_0028) + -IT_0040*(IT_0040*(16*s_13 + (-8)*IT_0046) + (-8)
      *IT_0047*IT_0050 + -IT_0051*IT_0054) + -IT_0055*(IT_0055*((-16)*IT_0042 +
       32*IT_0043 + 8*IT_0058) + 16*IT_0040*IT_0064 + (-8)*IT_0050*IT_0068 + 
      -IT_0054*IT_0072) + (-8)*IT_0024*((-0.125)*IT_0028*IT_0048 + conj(IT_0024)
      *IT_0077);
    return create_ccomplex_return(IT_0078);
}

