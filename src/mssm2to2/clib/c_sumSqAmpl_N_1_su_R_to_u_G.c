#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_su_R_to_u_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_su_R_to_u_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_N_1, 2);
    const ccomplex_t IT_0002 = pow(m_su_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0006;
    const ccomplex_t IT_0008 = 0.666666666666667*IT_0007;
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = (-2)*IT_0004;
    const ccomplex_t IT_0013 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (-2)*IT_0011 + -IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = sin(beta);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = sin(theta_W);
    const ccomplex_t IT_0022 = cpow(IT_0021, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0018*IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0013*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = IT_0012*IT_0026;
    const ccomplex_t IT_0028 = IT_0004*IT_0024;
    const ccomplex_t IT_0029 = IT_0003*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0027 + 2*IT_0030;
    const ccomplex_t IT_0032 = (-8)*IT_0000;
    const ccomplex_t IT_0033 = -IT_0004;
    const ccomplex_t IT_0034 = IT_0015*IT_0033;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0033;
    const ccomplex_t IT_0037 = (-8)*s_13;
    const ccomplex_t IT_0038 = s_34*IT_0037;
    const ccomplex_t IT_0039 = m_u*s_14*m_N_1;
    const ccomplex_t IT_0040 = -IT_0011;
    const ccomplex_t IT_0041 = s_34*IT_0001;
    const ccomplex_t IT_0042 = 8*IT_0041;
    const ccomplex_t IT_0043 = s_13*s_14;
    const ccomplex_t IT_0044 = (-16)*IT_0043;
    const ccomplex_t IT_0045 = IT_0042 + IT_0044;
    const ccomplex_t IT_0046 = s_14*IT_0037;
    const ccomplex_t IT_0047 = IT_0001*IT_0037;
    const ccomplex_t IT_0048 = pow(s_13, 2);
    const ccomplex_t IT_0049 = (-8)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0027;
    const ccomplex_t IT_0051 = m_u*m_N_1;
    const ccomplex_t IT_0052 = (-8)*s_14;
    const ccomplex_t IT_0053 = (-8)*s_34;
    const ccomplex_t IT_0054 = (-8)*IT_0001;
    const ccomplex_t IT_0055 = (conj(IT_0017)*IT_0031 + IT_0017*conj(IT_0031))
      *IT_0032 + IT_0037*(conj(IT_0016)*IT_0031 + IT_0016*conj(IT_0031) + conj
      (IT_0017)*IT_0050 + IT_0017*conj(IT_0050)) + (conj(IT_0016)*IT_0036 +
       IT_0016*conj(IT_0036) + conj(IT_0035)*IT_0050 + IT_0035*conj(IT_0050))
      *IT_0052 + (conj(IT_0031)*IT_0035 + IT_0031*conj(IT_0035) + conj(IT_0017)
      *IT_0036 + IT_0017*conj(IT_0036))*IT_0053 + (conj(IT_0016)*IT_0050 +
       IT_0016*conj(IT_0050))*IT_0054;
    const ccomplex_t IT_0056 = m_u*s_34*m_N_1;
    const ccomplex_t IT_0057 = (-8)*conj(IT_0017);
    const ccomplex_t IT_0058 = s_14*IT_0000;
    const ccomplex_t IT_0059 = (-8)*conj(IT_0031);
    const ccomplex_t IT_0060 = s_14*s_34;
    const ccomplex_t IT_0061 = (-8)*IT_0031;
    const ccomplex_t IT_0062 = (-0.125)*IT_0058;
    const ccomplex_t IT_0063 = (-8)*IT_0017;
    const ccomplex_t IT_0064 = (-0.125)*IT_0056;
    const ccomplex_t IT_0065 = s_13*(IT_0017*conj(IT_0017) + IT_0031*conj
      (IT_0031))*IT_0032 + (conj(IT_0017)*IT_0035 + IT_0017*conj(IT_0035) + conj
      (IT_0031)*IT_0036 + IT_0031*conj(IT_0036))*IT_0038 + (-8)*conj(IT_0016)*
      (IT_0030*IT_0039 + (-0.125)*IT_0040*IT_0045 + (-0.125)*IT_0035*IT_0046 + (
      -0.125)*IT_0016*IT_0047 + (-0.125)*IT_0017*IT_0049) + IT_0016*(conj
      (IT_0035)*IT_0046 + conj(IT_0017)*IT_0049) + (conj(IT_0036)*IT_0046 + conj
      (IT_0031)*IT_0049)*IT_0050 + (-8)*(IT_0039*IT_0040 + (-0.125)*IT_0030
      *IT_0045 + (-0.125)*IT_0036*IT_0046 + (-0.125)*IT_0031*IT_0049 + (-0.125)
      *IT_0047*IT_0050)*conj(IT_0050) + IT_0051*IT_0055 + IT_0040*(IT_0057
      *IT_0058 + IT_0056*IT_0059) + IT_0030*(IT_0056*IT_0057 + IT_0058*IT_0059) 
      + (-8)*conj(IT_0040)*((-0.125)*IT_0016*IT_0045 + IT_0039*IT_0050 + (-4)
      *IT_0040*IT_0060 + IT_0062*IT_0063 + IT_0061*IT_0064) + (-8)*conj(IT_0030)
      *(IT_0016*IT_0039 + (-0.125)*IT_0045*IT_0050 + (-4)*IT_0030*IT_0060 +
       IT_0061*IT_0062 + IT_0063*IT_0064);
    return create_ccomplex_return(IT_0065);
}

