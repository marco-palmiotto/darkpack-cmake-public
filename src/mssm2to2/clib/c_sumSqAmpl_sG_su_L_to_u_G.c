#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_su_L_to_u_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_su_L_to_u_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_u = param->m_u;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0003 = pow(m_u, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0002*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = -g_s;
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = pow(m_su_L, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0011 +
       reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = pow(s_13, 2);
    const ccomplex_t IT_0018 = 1.33333333333333*IT_0017;
    const ccomplex_t IT_0019 = -IT_0000;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0022 = s_34*IT_0021;
    const ccomplex_t IT_0023 = -IT_0008 + (0 + _Complex_I*-1)*IT_0014;
    const ccomplex_t IT_0024 = (-10.6666666666667)*IT_0017;
    const ccomplex_t IT_0025 = 1.33333333333333*s_13;
    const ccomplex_t IT_0026 = (-1.33333333333333)*s_13;
    const ccomplex_t IT_0027 = (-1.33333333333333)*IT_0017;
    const ccomplex_t IT_0028 = IT_0000*IT_0002;
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + IT_0003 + -IT_0004 + -IT_0011 
      + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = 2*IT_0031;
    const ccomplex_t IT_0033 = s_34*IT_0004;
    const ccomplex_t IT_0034 = s_13*s_14;
    const ccomplex_t IT_0035 = s_14*s_34;
    const ccomplex_t IT_0036 = (-2)*IT_0034;
    const ccomplex_t IT_0037 = IT_0033 + IT_0036;
    const ccomplex_t IT_0038 = s_14*IT_0003;
    const ccomplex_t IT_0039 = 1.33333333333333*conj(IT_0008);
    const ccomplex_t IT_0040 = (-1.33333333333333)*IT_0008;
    const ccomplex_t IT_0041 = IT_0016*((-10.6666666666667)*IT_0033 +
       21.3333333333333*IT_0034) + IT_0015*(10.6666666666667*IT_0033 + (
      -21.3333333333333)*IT_0034) + 21.3333333333333*IT_0031*IT_0035 +
       1.33333333333333*(IT_0015 + (-0.5)*IT_0023 + 0.5*conj(IT_0023))*IT_0037 +
       (-0.5)*IT_0038*(21.3333333333333*IT_0032 + IT_0039 + IT_0040);
    const ccomplex_t IT_0042 = 0.5*IT_0015;
    const ccomplex_t IT_0043 = s_13*s_34;
    const ccomplex_t IT_0044 = (-21.3333333333333)*IT_0043;
    const ccomplex_t IT_0045 = 10.6666666666667*IT_0038;
    const ccomplex_t IT_0046 = IT_0044 + IT_0045;
    const ccomplex_t IT_0047 = -IT_0008;
    const ccomplex_t IT_0048 = conj(IT_0008) + IT_0047;
    const ccomplex_t IT_0049 = 21.3333333333333*IT_0043;
    const ccomplex_t IT_0050 = (-10.6666666666667)*IT_0038;
    const ccomplex_t IT_0051 = IT_0049 + IT_0050;
    const ccomplex_t IT_0052 = (-2)*IT_0032;
    const ccomplex_t IT_0053 = (-2)*IT_0043;
    const ccomplex_t IT_0054 = IT_0038 + IT_0053;
    const ccomplex_t IT_0055 = ((-48)*IT_0015 + 24*IT_0016 + 12*IT_0023 + (-12
      )*conj(IT_0023))*IT_0033 + IT_0035*((-96)*IT_0031 + (-96)*IT_0042) +
       IT_0046*IT_0048 + IT_0051*IT_0052 + 2.66666666666667*(IT_0032 + 0.375
      *IT_0039 + 0.375*IT_0040)*IT_0054;
    const ccomplex_t IT_0056 = (-24)*IT_0016;
    const ccomplex_t IT_0057 = IT_0015*(s_14*((-12)*IT_0020 + 12*conj(IT_0020)
      ) + s_13*((-12)*IT_0008 + 12*conj(IT_0008) + 24*IT_0032) + IT_0004*(24
      *IT_0015 + (-12)*IT_0023 + 12*conj(IT_0023) + IT_0056));
    const ccomplex_t IT_0058 = (-10.6666666666667)*IT_0003*IT_0008*conj
      (IT_0008) + IT_0004*((-10.6666666666667)*IT_0023*conj(IT_0023) + IT_0016*
      (10.6666666666667*IT_0016 + 1.33333333333333*IT_0023 + (-1.33333333333333)
      *conj(IT_0023))) + IT_0032*((-21.3333333333333)*s_13*IT_0016 + IT_0003*((
      -1.33333333333333)*IT_0008 + 1.33333333333333*conj(IT_0008) +
       10.6666666666667*IT_0032)) + IT_0057;
    const ccomplex_t IT_0059 = (-0.5)*IT_0015;
    const ccomplex_t IT_0060 = IT_0048*IT_0051 + IT_0046*IT_0052 + (
      -2.66666666666667)*((-0.5)*IT_0008 + 0.5*conj(IT_0008) + IT_0032)*IT_0054 
      + IT_0033*(48*IT_0015 + (-12)*IT_0023 + 12*conj(IT_0023) + IT_0056) +
       IT_0035*(96*IT_0031 + 192*IT_0042 + (-96)*IT_0059);
    const ccomplex_t IT_0061 = IT_0008*(IT_0016*IT_0018 + conj(IT_0020)
      *IT_0022 + conj(IT_0023)*IT_0024) + s_14*(IT_0021*(conj(IT_0020)*IT_0023 +
       IT_0020*conj(IT_0023)) + IT_0016*(IT_0020*IT_0025 + conj(IT_0020)*IT_0026
      )) + conj(IT_0008)*(IT_0020*IT_0022 + IT_0023*IT_0024 + IT_0016*IT_0027) +
       (IT_0018*conj(IT_0023) + s_34*(conj(IT_0020)*IT_0025 + IT_0020*IT_0026) +
       IT_0023*IT_0027)*IT_0032 + (-2)*IT_0031*IT_0041 + IT_0042*IT_0055 + s_13
      *IT_0058 + IT_0059*IT_0060;
    return create_ccomplex_return(IT_0061);
}

