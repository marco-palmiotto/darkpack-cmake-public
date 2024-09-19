#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_se_L_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_se_L_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0000*IT_0005*(IT_0007 + 
      -IT_0009);
    const ccomplex_t IT_0011 = 0.25*IT_0010;
    const ccomplex_t IT_0012 = pow(m_se_L, 2);
    const ccomplex_t IT_0013 = s_12 + IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = 2*IT_0014;
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = s_13 + s_23;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0018;
    const ccomplex_t IT_0022 = cpow(IT_0017, 2);
    const ccomplex_t IT_0023 = s_23*IT_0022;
    const ccomplex_t IT_0024 = (-0.5)*IT_0020 + (-0.5)*IT_0021 + 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0016 + IT_0024;
    const ccomplex_t IT_0026 = cpow(s_12 + IT_0012 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0001, -1);
    const ccomplex_t IT_0030 = cpow(IT_0003, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0029 + 
      -IT_0001*IT_0030);
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*s_13 + (-0.5)*s_23;
    const ccomplex_t IT_0035 = IT_0013 + IT_0034;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0031;
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = IT_0013*IT_0038;
    const ccomplex_t IT_0040 = -IT_0031;
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = IT_0013*IT_0041;
    const ccomplex_t IT_0043 = IT_0031*IT_0040;
    const ccomplex_t IT_0044 = -s_23;
    const ccomplex_t IT_0045 = IT_0013 + IT_0044;
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = -IT_0039 + 0.5*IT_0042 + (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0036 + IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = pow(m_W, -1);
    const ccomplex_t IT_0053 = pow(m_e, 2);
    const ccomplex_t IT_0054 = cpow(IT_0006, -1);
    const ccomplex_t IT_0055 = sin(alpha);
    const ccomplex_t IT_0056 = cos(alpha);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*e_em*(IT_0030*IT_0052
      *IT_0053*IT_0054*IT_0055 + 0.5*m_W*(IT_0002*IT_0003 + -IT_0030)*(IT_0006
      *IT_0055 + IT_0008*IT_0056));
    const ccomplex_t IT_0058 = IT_0002*IT_0003;
    const ccomplex_t IT_0059 = -IT_0030 + IT_0058;
    const ccomplex_t IT_0060 = cpow(IT_0008, 3);
    const ccomplex_t IT_0061 = 3*IT_0030 + IT_0058;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_W*e_em*(IT_0056*IT_0059
      *IT_0060 + IT_0006*(IT_0009*IT_0055*IT_0061 + -IT_0006*(IT_0006*IT_0055
      *IT_0059 + IT_0008*IT_0056*IT_0061)));
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = IT_0057*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*e_em*(IT_0030*IT_0052
      *IT_0053*IT_0054*IT_0056 + (-0.5)*m_W*(IT_0002*IT_0003 + -IT_0030)*
      (IT_0008*IT_0055 + -IT_0006*IT_0056));
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = -IT_0030;
    const ccomplex_t IT_0071 = IT_0058 + IT_0070;
    const ccomplex_t IT_0072 = -IT_0006;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_W*e_em*(IT_0055*IT_0060
      *IT_0071 + IT_0072*(IT_0009*IT_0056*IT_0061 + -(IT_0008*IT_0055*IT_0061 + 
      -IT_0006*IT_0056*IT_0071)*IT_0072));
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0069*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + (-2)*IT_0012 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0008
      *IT_0030*(m_W*IT_0006 + (-0.5)*IT_0052*IT_0053*IT_0054);
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = cpow(IT_0080, 2);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = pow(m_Hp, 2);
    const ccomplex_t IT_0084 = cpow((-2)*s_13 + IT_0012 + IT_0083 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = IT_0011 + IT_0028 + 2*IT_0051 + -IT_0067 + 
      -IT_0078 + IT_0085;
    return create_ccomplex_return(IT_0086*conj(IT_0086));
}

