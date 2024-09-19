#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_2_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_2_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = -IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0000, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0000
      *IT_0002 + -IT_0001*IT_0005)*U_stau_00*conj(U_stau_01) + IT_0001*IT_0005
      *U_stau_10*conj(U_stau_11));
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 2*IT_0003;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = pow(m_stau_1, 2);
    const ccomplex_t IT_0012 = pow(m_stau_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (IT_0008 + IT_0010)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, -2);
    const ccomplex_t IT_0018 = pow(s_13, 2);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = pow(s_14, 2);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (IT_0011 + -IT_0019)*(IT_0011 + -IT_0021);
    const ccomplex_t IT_0023 = (-2)*IT_0006;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = IT_0009*IT_0023;
    const ccomplex_t IT_0026 = -IT_0014*(IT_0008 + IT_0010 + -IT_0024 + 
      -IT_0025);
    const ccomplex_t IT_0027 = s_13*s_23*IT_0017;
    const ccomplex_t IT_0028 = s_14*s_24*IT_0017;
    const ccomplex_t IT_0029 = (s_12 + -IT_0027)*(s_12 + -IT_0028);
    const ccomplex_t IT_0030 = IT_0003*IT_0007;
    const ccomplex_t IT_0031 = IT_0003*IT_0023;
    const ccomplex_t IT_0032 = IT_0014*(IT_0030 + -IT_0031);
    const ccomplex_t IT_0033 = s_14*s_34*IT_0017;
    const ccomplex_t IT_0034 = (s_12 + -IT_0027)*(s_13 + -IT_0033);
    const ccomplex_t IT_0035 = pow(e_em, 2);
    const ccomplex_t IT_0036 = cpow(IT_0001, -2);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0035*IT_0036*U_stau_00
      *conj(U_stau_01);
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = (-2)*IT_0003;
    const ccomplex_t IT_0040 = IT_0007*IT_0039;
    const ccomplex_t IT_0041 = s_13 + s_23;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = s_12 + IT_0012;
    const ccomplex_t IT_0044 = IT_0024*IT_0043;
    const ccomplex_t IT_0045 = IT_0042 + IT_0044;
    const ccomplex_t IT_0046 = IT_0023*IT_0039;
    const ccomplex_t IT_0047 = s_23*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0011 + 0.5*IT_0012;
    const ccomplex_t IT_0049 = s_12 + IT_0048;
    const ccomplex_t IT_0050 = IT_0008*IT_0049;
    const ccomplex_t IT_0051 = -IT_0047 + (-2)*IT_0050;
    const ccomplex_t IT_0052 = IT_0045 + IT_0051;
    const ccomplex_t IT_0053 = IT_0014*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0055 = cos(alpha);
    const ccomplex_t IT_0056 = sin(beta);
    const ccomplex_t IT_0057 = cos(beta);
    const ccomplex_t IT_0058 = sin(alpha);
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_W*e_em*IT_0002*(IT_0055
      *IT_0056 + -IT_0057*IT_0058);
    const ccomplex_t IT_0060 = pow(m_W, -1);
    const ccomplex_t IT_0061 = mu_h*m_tau*IT_0060;
    const ccomplex_t IT_0062 = cpow(IT_0057, -1);
    const ccomplex_t IT_0063 = IT_0057*IT_0058;
    const ccomplex_t IT_0064 = IT_0055*IT_0056;
    const ccomplex_t IT_0065 = IT_0063 + IT_0064;
    const ccomplex_t IT_0066 = cpow(IT_0000, -2);
    const ccomplex_t IT_0067 = pow(m_tau, 2);
    const ccomplex_t IT_0068 = 2*U_stau_00*(((-0.5)*m_W*IT_0065*(IT_0002 + 
      -IT_0001*IT_0066) + IT_0002*IT_0058*IT_0060*IT_0062*IT_0067)*conj
      (U_stau_01) + 1./2*IT_0002*IT_0055*IT_0061*IT_0062*conj(U_stau_11)) +
       U_stau_10*(IT_0002*IT_0055*IT_0061*IT_0062*conj(U_stau_01) + (-2)*(m_W
      *IT_0001*IT_0065*IT_0066 + -IT_0002*IT_0058*IT_0060*IT_0062*IT_0067)*conj
      (U_stau_11));
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(e_em*IT_0068 +
       1.4142135623731*A_tau*IT_0058*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0059*IT_0070;
    const ccomplex_t IT_0072 = IT_0054*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_W*e_em*IT_0002*(IT_0055
      *IT_0057 + IT_0056*IT_0058);
    const ccomplex_t IT_0076 = -IT_0056*IT_0058;
    const ccomplex_t IT_0077 = IT_0055*IT_0057;
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = (-2)*U_stau_00*((IT_0002*IT_0055*IT_0060
      *IT_0062*IT_0067 + (-0.5)*m_W*(IT_0002 + -IT_0001*IT_0066)*IT_0078)*conj
      (U_stau_01) + -1./2*IT_0002*IT_0058*IT_0061*IT_0062*conj(U_stau_11)) +
       U_stau_10*(IT_0002*IT_0058*IT_0061*IT_0062*conj(U_stau_01) + (-2)*
      (IT_0002*IT_0055*IT_0060*IT_0062*IT_0067 + -m_W*IT_0001*IT_0066*IT_0078)
      *conj(U_stau_11));
    const ccomplex_t IT_0080 = (0 + _Complex_I*-1)*(e_em*IT_0079 + (
      -1.4142135623731)*A_tau*IT_0055*(conj(U_stau_01)*U_stau_10 + U_stau_00
      *conj(U_stau_11)));
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = IT_0075*IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0038 + IT_0053 + -IT_0073 + -IT_0084;
    const ccomplex_t IT_0086 = pow(m_W, -4);
    const ccomplex_t IT_0087 = s_13*s_14*s_34*IT_0086;
    const ccomplex_t IT_0088 = -IT_0021;
    const ccomplex_t IT_0089 = IT_0011 + -IT_0019 + IT_0087 + IT_0088;
    const ccomplex_t IT_0090 = pow(s_34, 2);
    const ccomplex_t IT_0091 = IT_0086*IT_0090;
    const ccomplex_t IT_0092 = IT_0014*IT_0030;
    const ccomplex_t IT_0093 = s_34*IT_0017;
    const ccomplex_t IT_0094 = s_13*s_34*IT_0017;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = s_14 + IT_0095;
    const ccomplex_t IT_0097 = IT_0093*IT_0096;
    const ccomplex_t IT_0098 = s_23*s_34*IT_0017;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = s_24 + IT_0099;
    const ccomplex_t IT_0101 = IT_0093*IT_0100;
    const ccomplex_t IT_0102 = IT_0008 + IT_0010;
    const ccomplex_t IT_0103 = -IT_0025;
    const ccomplex_t IT_0104 = IT_0102 + IT_0103;
    const ccomplex_t IT_0105 = IT_0014*IT_0104;
    const ccomplex_t IT_0106 = pow(m_W, 2);
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0011 + IT_0106 + reg_prop +
       -m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0002
      *conj(U_stau_01);
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0002
      *U_stau_00;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0110*IT_0112;
    const ccomplex_t IT_0114 = IT_0108*IT_0113;
    const ccomplex_t IT_0115 = -IT_0105 + -IT_0114;
    const ccomplex_t IT_0116 = s_13*s_24*s_34*IT_0086;
    const ccomplex_t IT_0117 = -IT_0028;
    const ccomplex_t IT_0118 = -IT_0027;
    const ccomplex_t IT_0119 = s_12 + IT_0117 + IT_0118;
    const ccomplex_t IT_0120 = IT_0116 + IT_0119;
    const ccomplex_t IT_0121 = -IT_0014*(IT_0008 + IT_0010 + -IT_0024);
    const ccomplex_t IT_0122 = s_14*s_23*s_34*IT_0086;
    const ccomplex_t IT_0123 = IT_0119 + IT_0122;
    const ccomplex_t IT_0124 = s_23*s_24*s_34*IT_0086;
    const ccomplex_t IT_0125 = pow(s_24, 2);
    const ccomplex_t IT_0126 = IT_0017*IT_0125;
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = pow(s_23, 2);
    const ccomplex_t IT_0129 = IT_0017*IT_0128;
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0012 + IT_0124 + IT_0127 + IT_0130;
    const ccomplex_t IT_0132 = 0.5*conj(IT_0016)*IT_0089 + conj(IT_0085)*(1 +
       0.5*IT_0091) + (-0.5)*conj(IT_0092)*IT_0097 + (-0.5)*conj(IT_0032)
      *IT_0101 + 0.5*conj(IT_0115)*IT_0120 + 0.5*conj(IT_0121)*IT_0123 + 0.5
      *conj(IT_0026)*IT_0131;
    const ccomplex_t IT_0133 = (IT_0011 + -IT_0019)*(s_12 + -IT_0028);
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = (IT_0011 + -IT_0019)*(s_13 + -IT_0033);
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = (IT_0011 + -IT_0021)*(s_12 + -IT_0027);
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = -IT_0016;
    const ccomplex_t IT_0140 = (IT_0011 + -IT_0019)*(IT_0012 + -IT_0126);
    const ccomplex_t IT_0141 = s_24*s_34*IT_0017;
    const ccomplex_t IT_0142 = (s_12 + -IT_0027)*(s_23 + -IT_0141);
    const ccomplex_t IT_0143 = -conj(IT_0016);
    const ccomplex_t IT_0144 = (IT_0011 + -IT_0021)*(IT_0012 + -IT_0129);
    const ccomplex_t IT_0145 = (IT_0011 + -IT_0019)*(s_23 + -IT_0141);
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = (s_12 + -IT_0027)*(IT_0012 + -IT_0126);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = -IT_0115;
    const ccomplex_t IT_0150 = IT_0017*IT_0090;
    const ccomplex_t IT_0151 = (IT_0011 + -IT_0019)*(IT_0106 + -IT_0150);
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = (s_12 + -IT_0027)*(IT_0106 + -IT_0150);
    const ccomplex_t IT_0154 = -conj(IT_0085);
    const ccomplex_t IT_0155 = -conj(IT_0115);
    const ccomplex_t IT_0156 = (s_12 + -IT_0028)*(IT_0012 + -IT_0129);
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = (s_13 + -IT_0033)*(IT_0012 + -IT_0129);
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = -IT_0121;
    const ccomplex_t IT_0161 = (IT_0012 + -IT_0126)*(IT_0012 + -IT_0129);
    const ccomplex_t IT_0162 = -conj(IT_0121);
    const ccomplex_t IT_0163 = (IT_0012 + -IT_0129)*(s_23 + -IT_0141);
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = -conj(IT_0032);
    const ccomplex_t IT_0166 = (IT_0012 + -IT_0129)*(IT_0106 + -IT_0150);
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0016*(conj(IT_0016)*IT_0022 + conj(IT_0026)
      *IT_0029 + conj(IT_0032)*IT_0034 + conj(IT_0085)*IT_0089) + 2*IT_0085
      *IT_0132 + (conj(IT_0115)*IT_0134 + conj(IT_0092)*IT_0136 + conj(IT_0121)
      *IT_0138)*IT_0139 + IT_0115*(conj(IT_0085)*IT_0120 + IT_0029*conj(IT_0121)
       + conj(IT_0115)*IT_0140 + conj(IT_0032)*IT_0142 + IT_0134*IT_0143) +
       IT_0121*(IT_0034*conj(IT_0092) + IT_0029*conj(IT_0115) + conj(IT_0085)
      *IT_0123 + IT_0138*IT_0143 + conj(IT_0121)*IT_0144) + (conj(IT_0092)
      *IT_0146 + conj(IT_0026)*IT_0148)*IT_0149 + IT_0092*(IT_0034*conj(IT_0121)
       + conj(IT_0026)*IT_0142 + IT_0136*IT_0143 + -conj(IT_0092)*IT_0152 + conj
      (IT_0032)*IT_0153 + IT_0097*IT_0154 + IT_0146*IT_0155) + (conj(IT_0026)
      *IT_0157 + conj(IT_0032)*IT_0159)*IT_0160 + IT_0026*(conj(IT_0016)*IT_0029
       + conj(IT_0085)*IT_0131 + conj(IT_0092)*IT_0142 + IT_0148*IT_0155 + conj
      (IT_0026)*IT_0161 + IT_0157*IT_0162 + IT_0164*IT_0165) + IT_0032*(conj
      (IT_0016)*IT_0034 + conj(IT_0115)*IT_0142 + conj(IT_0092)*IT_0153 +
       IT_0101*IT_0154 + IT_0159*IT_0162 + -conj(IT_0026)*IT_0164 + IT_0165
      *IT_0167);
    return create_ccomplex_return(IT_0168);
}

