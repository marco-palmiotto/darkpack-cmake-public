#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_snu_tau_to_A_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_snu_tau_to_A_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = cpow(s_13 + (-0.5)*reg_prop + (0 + _Complex_I*(
      -0.5))*m_stau_1*Gamma_tau1, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*IT_0000;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *U_stau_00;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = 2*IT_0002;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = IT_0001*(IT_0007 + -IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = pow(m_stau_1, 2);
    const ccomplex_t IT_0013 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = (-2)*IT_0002;
    const ccomplex_t IT_0017 = IT_0006*IT_0016;
    const ccomplex_t IT_0018 = IT_0007 + IT_0017;
    const ccomplex_t IT_0019 = -IT_0005;
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = IT_0018 + IT_0021;
    const ccomplex_t IT_0023 = IT_0015*IT_0022;
    const ccomplex_t IT_0024 = IT_0008*IT_0019;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0009 + IT_0025;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = IT_0023 + (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0002;
    const ccomplex_t IT_0030 = IT_0006*IT_0029;
    const ccomplex_t IT_0031 = IT_0020 + IT_0030;
    const ccomplex_t IT_0032 = IT_0019*IT_0029;
    const ccomplex_t IT_0033 = -IT_0017 + -IT_0032;
    const ccomplex_t IT_0034 = IT_0031 + IT_0033;
    const ccomplex_t IT_0035 = IT_0015*IT_0034;
    const ccomplex_t IT_0036 = IT_0002*IT_0019;
    const ccomplex_t IT_0037 = IT_0009 + IT_0036;
    const ccomplex_t IT_0038 = -IT_0007 + -IT_0024;
    const ccomplex_t IT_0039 = IT_0037 + IT_0038;
    const ccomplex_t IT_0040 = IT_0001*IT_0039;
    const ccomplex_t IT_0041 = IT_0035 + 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0001*IT_0009;
    const ccomplex_t IT_0043 = -IT_0009;
    const ccomplex_t IT_0044 = IT_0007 + IT_0043;
    const ccomplex_t IT_0045 = IT_0015*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0042 + -IT_0045;
    const ccomplex_t IT_0047 = pow(e_em, 2);
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*IT_0004
      *IT_0047*U_stau_00;
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = s_13 + s_23;
    const ccomplex_t IT_0051 = IT_0030*IT_0050;
    const ccomplex_t IT_0052 = s_23*IT_0032;
    const ccomplex_t IT_0053 = -s_23;
    const ccomplex_t IT_0054 = s_12 + IT_0013;
    const ccomplex_t IT_0055 = IT_0053 + IT_0054;
    const ccomplex_t IT_0056 = IT_0036*IT_0055;
    const ccomplex_t IT_0057 = (-0.5)*s_13 + (-0.5)*s_23 + 0.5*IT_0012 + 0.5
      *IT_0013;
    const ccomplex_t IT_0058 = s_12 + IT_0057;
    const ccomplex_t IT_0059 = IT_0007*IT_0058;
    const ccomplex_t IT_0060 = -IT_0052 + -IT_0056 + 2*IT_0059;
    const ccomplex_t IT_0061 = IT_0051 + IT_0060;
    const ccomplex_t IT_0062 = IT_0015*IT_0061;
    const ccomplex_t IT_0063 = cos(beta);
    const ccomplex_t IT_0064 = pow(m_W, -1);
    const ccomplex_t IT_0065 = cpow(IT_0063, -1);
    const ccomplex_t IT_0066 = sin(beta);
    const ccomplex_t IT_0067 = pow(m_tau, 2);
    const ccomplex_t IT_0068 = cpow(IT_0063, 2);
    const ccomplex_t IT_0069 = cpow(IT_0066, 2);
    const ccomplex_t IT_0070 = A_tau*IT_0063*U_stau_10 + 0.707106781186548
      *e_em*IT_0004*((IT_0064*IT_0067 + -m_W*(IT_0068 + -IT_0069))*U_stau_00 + 
      -mu_h*m_tau*IT_0064*IT_0065*IT_0066*U_stau_10);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_W*e_em;
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = pow(m_W, 2);
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + IT_0075 
      + -reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0049 + IT_0062 + -IT_0078;
    const ccomplex_t IT_0080 = -IT_0015*(IT_0007 + -IT_0009 + IT_0024 + 
      -IT_0036);
    const ccomplex_t IT_0081 = -IT_0015*(IT_0007 + IT_0017 + -IT_0020 + 
      -IT_0036);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = IT_0012*(conj(IT_0028)*IT_0046 + IT_0028*conj
      (IT_0046)) + s_13*(conj(IT_0011)*IT_0028 + IT_0011*conj(IT_0028) + conj
      (IT_0041)*IT_0046 + IT_0041*conj(IT_0046)) + conj(IT_0079)*(IT_0041 +
       IT_0080) + IT_0079*(conj(IT_0041) + conj(IT_0080)) + s_23*(conj(IT_0041)
      *IT_0080 + IT_0041*conj(IT_0080) + conj(IT_0011)*IT_0082 + IT_0011*conj
      (IT_0082)) + s_12*(conj(IT_0028)*IT_0080 + IT_0028*conj(IT_0080) + conj
      (IT_0046)*IT_0082 + IT_0046*conj(IT_0082)) + IT_0013*(conj(IT_0080)
      *IT_0082 + IT_0080*conj(IT_0082));
    const ccomplex_t IT_0084 = pow(m_W, -2);
    const ccomplex_t IT_0085 = s_24*s_34*IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = s_23 + IT_0086;
    const ccomplex_t IT_0088 = IT_0015*IT_0018;
    const ccomplex_t IT_0089 = (-0.5)*IT_0042 + IT_0088;
    const ccomplex_t IT_0090 = pow(s_14, 2);
    const ccomplex_t IT_0091 = IT_0084*IT_0090;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = IT_0012 + IT_0092;
    const ccomplex_t IT_0094 = pow(s_24, 2);
    const ccomplex_t IT_0095 = IT_0084*IT_0094;
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = IT_0013 + IT_0096;
    const ccomplex_t IT_0098 = -IT_0015*(IT_0007 + IT_0017 + -IT_0036);
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = s_12*conj(IT_0089);
    const ccomplex_t IT_0101 = IT_0013*conj(IT_0099);
    const ccomplex_t IT_0102 = IT_0100 + IT_0101;
    const ccomplex_t IT_0103 = -IT_0017;
    const ccomplex_t IT_0104 = IT_0030 + IT_0103;
    const ccomplex_t IT_0105 = IT_0015*IT_0104;
    const ccomplex_t IT_0106 = IT_0001*IT_0044;
    const ccomplex_t IT_0107 = IT_0105 + (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = s_13*IT_0089;
    const ccomplex_t IT_0109 = s_23*IT_0099;
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = s_13*conj(IT_0089);
    const ccomplex_t IT_0112 = s_23*conj(IT_0099);
    const ccomplex_t IT_0113 = IT_0111 + IT_0112;
    const ccomplex_t IT_0114 = s_12*IT_0089;
    const ccomplex_t IT_0115 = -IT_0013;
    const ccomplex_t IT_0116 = IT_0095 + IT_0115;
    const ccomplex_t IT_0117 = -s_13;
    const ccomplex_t IT_0118 = -s_12;
    const ccomplex_t IT_0119 = -conj(IT_0041);
    const ccomplex_t IT_0120 = IT_0041*(IT_0053*conj(IT_0082) + conj(IT_0028)
      *IT_0117) + -IT_0082*(IT_0013*conj(IT_0082) + -conj(IT_0028)*IT_0118 +
       IT_0053*IT_0119) + -IT_0028*(IT_0012*conj(IT_0028) + -conj(IT_0082)
      *IT_0118 + IT_0117*IT_0119);
    const ccomplex_t IT_0121 = pow(s_34, 2);
    const ccomplex_t IT_0122 = IT_0084*IT_0121;
    const ccomplex_t IT_0123 = -conj(IT_0046);
    const ccomplex_t IT_0124 = -conj(IT_0080);
    const ccomplex_t IT_0125 = conj(IT_0011)*(IT_0079 + -IT_0053*IT_0080 + 
      -IT_0046*IT_0117) + IT_0080*(IT_0013*conj(IT_0080) + IT_0118*IT_0123) +
       IT_0011*(conj(IT_0079) + IT_0117*IT_0123 + IT_0053*IT_0124) + IT_0046*
      (IT_0012*conj(IT_0046) + IT_0118*IT_0124);
    const ccomplex_t IT_0126 = s_14*s_24*IT_0084;
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = s_12 + IT_0127;
    const ccomplex_t IT_0129 = s_12*IT_0099;
    const ccomplex_t IT_0130 = s_13*IT_0107;
    const ccomplex_t IT_0131 = IT_0012*IT_0089;
    const ccomplex_t IT_0132 = IT_0079 + IT_0129 + IT_0130 + IT_0131;
    const ccomplex_t IT_0133 = s_23*IT_0107;
    const ccomplex_t IT_0134 = IT_0013*IT_0099 + IT_0114 + IT_0133;
    const ccomplex_t IT_0135 = s_23*conj(IT_0107);
    const ccomplex_t IT_0136 = IT_0102 + IT_0135;
    const ccomplex_t IT_0137 = s_12*conj(IT_0099);
    const ccomplex_t IT_0138 = s_13*conj(IT_0107);
    const ccomplex_t IT_0139 = IT_0012*conj(IT_0089);
    const ccomplex_t IT_0140 = s_14*s_34*IT_0084;
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = s_13 + IT_0141;
    const ccomplex_t IT_0143 = IT_0083*IT_0087 + conj(IT_0079)*(IT_0089
      *IT_0093 + IT_0082*IT_0097) + 3*IT_0079*(conj(IT_0079) + 0.333333333333333
      *conj(IT_0089)*IT_0093 + 0.333333333333333*conj(IT_0082)*IT_0097) +
       IT_0093*(IT_0012*IT_0089*conj(IT_0089) + IT_0099*IT_0102 + conj(IT_0107)
      *IT_0110 + IT_0107*IT_0113 + conj(IT_0099)*IT_0114) + IT_0116*IT_0120 + 
      -IT_0122*IT_0125 + IT_0128*(conj(IT_0079)*IT_0099 + IT_0079*conj(IT_0099) 
      + conj(IT_0041)*IT_0110 + IT_0041*IT_0113 + conj(IT_0028)*IT_0132 + conj
      (IT_0082)*IT_0134 + IT_0082*IT_0136 + IT_0028*(conj(IT_0079) + IT_0137 +
       IT_0138 + IT_0139)) + (conj(IT_0079)*IT_0107 + IT_0079*conj(IT_0107) +
       conj(IT_0011)*IT_0110 + IT_0011*IT_0113 + conj(IT_0046)*IT_0132 + conj
      (IT_0080)*IT_0134 + IT_0080*IT_0136 + IT_0046*(conj(IT_0079) + IT_0137 +
       IT_0138 + IT_0139))*IT_0142;
    return create_ccomplex_return(IT_0143);
}

