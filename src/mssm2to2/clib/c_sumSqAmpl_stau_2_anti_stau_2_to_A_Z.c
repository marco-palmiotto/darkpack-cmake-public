#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_A_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t theta_W = param->theta_W;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0001, -2);
    const ccomplex_t IT_0005 = cpow(IT_0003, -2);
    const ccomplex_t IT_0006 = tan(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = cpow(1 + IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0000*IT_0002*IT_0003
      *IT_0008*(0.5*(IT_0004 + -IT_0005)*U_stau_01*conj(U_stau_01) + IT_0004
      *U_stau_11*conj(U_stau_11));
    const ccomplex_t IT_0010 = (-2)*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0010, 2);
    const ccomplex_t IT_0012 = pow(m_Z, -2);
    const ccomplex_t IT_0013 = s_14*s_24*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = s_12 + IT_0014;
    const ccomplex_t IT_0016 = cpow(s_23 + (-0.5)*reg_prop + (0 + _Complex_I*(
      -0.5))*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0019 = 2*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0003, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0002*IT_0003 +
       -IT_0001*IT_0020)*U_stau_01*conj(U_stau_01) + IT_0002*IT_0003*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = cpow(s_13 + (-0.5)*reg_prop + (0 + _Complex_I*(
      -0.5))*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (-2)*IT_0021;
    const ccomplex_t IT_0029 = IT_0019*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*(IT_0023 + -IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0023*IT_0027;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = pow(m_stau_2, 2);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = IT_0010 + IT_0035;
    const ccomplex_t IT_0037 = s_12*IT_0025;
    const ccomplex_t IT_0038 = IT_0018*IT_0022;
    const ccomplex_t IT_0039 = IT_0017*IT_0038;
    const ccomplex_t IT_0040 = -IT_0038;
    const ccomplex_t IT_0041 = IT_0023 + IT_0040;
    const ccomplex_t IT_0042 = IT_0027*IT_0041;
    const ccomplex_t IT_0043 = (-0.5)*IT_0039 + 0.5*IT_0042;
    const ccomplex_t IT_0044 = s_13*IT_0043;
    const ccomplex_t IT_0045 = IT_0036 + IT_0037 + IT_0044;
    const ccomplex_t IT_0046 = -IT_0010;
    const ccomplex_t IT_0047 = s_12*conj(IT_0025);
    const ccomplex_t IT_0048 = s_13*conj(IT_0043);
    const ccomplex_t IT_0049 = conj(IT_0033)*IT_0034 + IT_0046 + IT_0047 +
       IT_0048;
    const ccomplex_t IT_0050 = IT_0017*(IT_0023 + -IT_0029);
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = s_12*IT_0033;
    const ccomplex_t IT_0053 = s_23*IT_0043;
    const ccomplex_t IT_0054 = IT_0025*IT_0034 + IT_0052 + IT_0053;
    const ccomplex_t IT_0055 = s_12*conj(IT_0033);
    const ccomplex_t IT_0056 = conj(IT_0025)*IT_0034;
    const ccomplex_t IT_0057 = IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = s_23*conj(IT_0043);
    const ccomplex_t IT_0059 = IT_0057 + IT_0058;
    const ccomplex_t IT_0060 = IT_0018*IT_0028;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0038 + IT_0061;
    const ccomplex_t IT_0063 = IT_0017*IT_0062;
    const ccomplex_t IT_0064 = IT_0023 + IT_0060;
    const ccomplex_t IT_0065 = -IT_0029 + -IT_0038;
    const ccomplex_t IT_0066 = IT_0064 + IT_0065;
    const ccomplex_t IT_0067 = IT_0027*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0063 + 0.5*IT_0067;
    const ccomplex_t IT_0069 = s_13*IT_0033;
    const ccomplex_t IT_0070 = s_23*IT_0025;
    const ccomplex_t IT_0071 = IT_0069 + IT_0070;
    const ccomplex_t IT_0072 = s_13*conj(IT_0033);
    const ccomplex_t IT_0073 = s_23*conj(IT_0025);
    const ccomplex_t IT_0074 = IT_0072 + IT_0073;
    const ccomplex_t IT_0075 = s_14*s_34*IT_0012;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = s_13 + IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0032;
    const ccomplex_t IT_0079 = (-0.5)*IT_0050;
    const ccomplex_t IT_0080 = (-0.5)*IT_0042 + 0.5*IT_0063;
    const ccomplex_t IT_0081 = pow(s_14, 2);
    const ccomplex_t IT_0082 = IT_0012*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = IT_0034 + IT_0083;
    const ccomplex_t IT_0085 = pow(s_24, 2);
    const ccomplex_t IT_0086 = IT_0012*IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0034 + IT_0087;
    const ccomplex_t IT_0089 = s_13*IT_0031;
    const ccomplex_t IT_0090 = s_23*IT_0051;
    const ccomplex_t IT_0091 = IT_0089 + IT_0090;
    const ccomplex_t IT_0092 = s_13*conj(IT_0031);
    const ccomplex_t IT_0093 = s_23*conj(IT_0051);
    const ccomplex_t IT_0094 = IT_0092 + IT_0093;
    const ccomplex_t IT_0095 = s_12*IT_0031;
    const ccomplex_t IT_0096 = s_12*conj(IT_0031);
    const ccomplex_t IT_0097 = s_24*s_34*IT_0012;
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = s_23 + IT_0098;
    const ccomplex_t IT_0100 = IT_0046*IT_0068 + IT_0010*conj(IT_0068) + (conj
      (IT_0031)*IT_0034 + s_12*conj(IT_0051) + s_13*conj(IT_0068))*IT_0078 + 
      (IT_0031*IT_0034 + s_12*IT_0051 + s_13*IT_0068)*conj(IT_0078) + conj
      (IT_0080)*IT_0091 + IT_0080*IT_0094 + conj(IT_0079)*(IT_0010 + IT_0034
      *IT_0051 + s_23*IT_0068 + IT_0095) + IT_0079*(IT_0046 + IT_0034*conj
      (IT_0051) + s_23*conj(IT_0068) + IT_0096);
    const ccomplex_t IT_0101 = pow(s_34, 2);
    const ccomplex_t IT_0102 = IT_0012*IT_0101;
    const ccomplex_t IT_0103 = -s_13;
    const ccomplex_t IT_0104 = -s_23;
    const ccomplex_t IT_0105 = -conj(IT_0080);
    const ccomplex_t IT_0106 = IT_0080*(IT_0010 + conj(IT_0078)*IT_0103 + conj
      (IT_0079)*IT_0104) + -IT_0046*IT_0105 + -IT_0078*(IT_0034*conj(IT_0078) +
       s_12*conj(IT_0079) + IT_0103*IT_0105) + -IT_0079*(s_12*conj(IT_0078) +
       IT_0034*conj(IT_0079) + IT_0104*IT_0105);
    const ccomplex_t IT_0107 = (-3)*IT_0011 + IT_0015*(IT_0010*conj(IT_0025) +
       conj(IT_0031)*IT_0045 + IT_0025*IT_0046 + IT_0031*IT_0049 + conj(IT_0051)
      *IT_0054 + IT_0051*IT_0059 + conj(IT_0068)*IT_0071 + IT_0068*IT_0074) +
       IT_0077*(IT_0010*conj(IT_0043) + IT_0043*IT_0046 + IT_0049*IT_0078 +
       IT_0045*conj(IT_0078) + IT_0059*IT_0079 + IT_0054*conj(IT_0079) + IT_0074
      *IT_0080 + IT_0071*conj(IT_0080)) + (conj(IT_0033)*IT_0036 + IT_0033
      *IT_0046 + conj(IT_0025)*IT_0052 + IT_0025*IT_0057 + conj(IT_0043)*IT_0071
       + IT_0043*IT_0074)*IT_0084 + IT_0088*(IT_0031*conj(IT_0031)*IT_0034 +
       conj(IT_0068)*IT_0091 + IT_0068*IT_0094 + conj(IT_0051)*(IT_0010 +
       IT_0034*IT_0051 + IT_0095) + IT_0051*(IT_0046 + IT_0096)) + IT_0099
      *IT_0100 + IT_0102*IT_0106;
    return create_ccomplex_return(IT_0107);
}

