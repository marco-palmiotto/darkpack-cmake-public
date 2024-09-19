#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_snu_tau_to_anti_tau_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_snu_tau_to_anti_tau_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = pow(m_C_1, 2);
    const ccomplex_t IT_0012 = pow(m_tau, 2);
    const ccomplex_t IT_0013 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0014 = cpow(s_12 + 0.5*IT_0011 + (-0.5)*IT_0012 + 0.5
      *IT_0013 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_tau*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = cos(alpha);
    const ccomplex_t IT_0021 = V_u1*U_Wm2;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = U_d2*V_Wp1;
    const ccomplex_t IT_0024 = IT_0003*IT_0023;
    const ccomplex_t IT_0025 = IT_0005*(IT_0022 + -IT_0024);
    const ccomplex_t IT_0026 = 1.4142135623731*e_em*IT_0025;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0019*IT_0028;
    const ccomplex_t IT_0030 = cpow(s_23 + (-0.5)*IT_0012 + (-0.5)*IT_0013 +
       0.5*m_C_2*(m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = V_u1*U_Wm1;
    const ccomplex_t IT_0036 = IT_0020*IT_0035;
    const ccomplex_t IT_0037 = U_d1*V_Wp1;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = IT_0005*(IT_0036 + -IT_0038);
    const ccomplex_t IT_0040 = 1.4142135623731*e_em*IT_0039;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0034*IT_0042;
    const ccomplex_t IT_0044 = cpow(s_23 + (-0.5)*IT_0012 + (-0.5)*IT_0013 +
       0.5*m_C_1*(m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0032 + IT_0046;
    const ccomplex_t IT_0048 = m_C_1*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = s_34*m_C_1;
    const ccomplex_t IT_0051 = 2*IT_0050;
    const ccomplex_t IT_0052 = conj(IT_0049)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0005;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0056 = IT_0020*IT_0055;
    const ccomplex_t IT_0057 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0058 = IT_0003*IT_0057;
    const ccomplex_t IT_0059 = IT_0005*(IT_0056 + -IT_0058);
    const ccomplex_t IT_0060 = 1.4142135623731*e_em*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0054*IT_0062;
    const ccomplex_t IT_0064 = IT_0030*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0067 = IT_0020*IT_0066;
    const ccomplex_t IT_0068 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0069 = IT_0003*IT_0068;
    const ccomplex_t IT_0070 = IT_0005*(IT_0067 + -IT_0069);
    const ccomplex_t IT_0071 = 1.4142135623731*e_em*IT_0070;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0009*IT_0073;
    const ccomplex_t IT_0075 = IT_0044*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0065 + IT_0076;
    const ccomplex_t IT_0078 = m_C_1*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = s_14*m_tau;
    const ccomplex_t IT_0081 = 2*IT_0080;
    const ccomplex_t IT_0082 = conj(IT_0079)*IT_0081;
    const ccomplex_t IT_0083 = IT_0007*IT_0034;
    const ccomplex_t IT_0084 = IT_0014*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = pow(m_h, 2);
    const ccomplex_t IT_0088 = m_C_1*m_tau*IT_0087;
    const ccomplex_t IT_0089 = 2*IT_0088;
    const ccomplex_t IT_0090 = conj(IT_0086)*IT_0089;
    const ccomplex_t IT_0091 = cpow(s_13 + (-0.5)*IT_0011 + (-0.5)*IT_0012 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0092 = cos(theta_W);
    const ccomplex_t IT_0093 = cpow(IT_0092, -2);
    const ccomplex_t IT_0094 = IT_0004*IT_0093;
    const ccomplex_t IT_0095 = IT_0005 + IT_0094;
    const ccomplex_t IT_0096 = sin(beta);
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_W*e_em*IT_0095*(IT_0001
      *IT_0003 + IT_0020*IT_0096);
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0034*IT_0098;
    const ccomplex_t IT_0100 = IT_0091*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0103 = IT_0019*IT_0062;
    const ccomplex_t IT_0104 = IT_0030*IT_0102*IT_0103;
    const ccomplex_t IT_0105 = m_tau*IT_0085;
    const ccomplex_t IT_0106 = IT_0034*IT_0073;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0108 = IT_0044*IT_0106*IT_0107;
    const ccomplex_t IT_0109 = 0.5*IT_0101 + (-0.5)*IT_0104 + -IT_0105 + (-0.5
      )*IT_0108;
    const ccomplex_t IT_0110 = (-2)*IT_0050;
    const ccomplex_t IT_0111 = conj(IT_0109)*IT_0110;
    const ccomplex_t IT_0112 = IT_0028*IT_0054;
    const ccomplex_t IT_0113 = IT_0030*IT_0102*IT_0112;
    const ccomplex_t IT_0114 = IT_0009*IT_0098;
    const ccomplex_t IT_0115 = IT_0091*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0009*IT_0042;
    const ccomplex_t IT_0118 = IT_0044*IT_0107*IT_0117;
    const ccomplex_t IT_0119 = m_tau*IT_0016;
    const ccomplex_t IT_0120 = 0.5*IT_0113 + (-0.5)*IT_0116 + 0.5*IT_0118 +
       IT_0119;
    const ccomplex_t IT_0121 = (-2)*IT_0080;
    const ccomplex_t IT_0122 = conj(IT_0120)*IT_0121;
    const ccomplex_t IT_0123 = 0.5*IT_0065 + 0.5*IT_0076;
    const ccomplex_t IT_0124 = IT_0089*conj(IT_0123);
    const ccomplex_t IT_0125 = IT_0052 + IT_0082 + IT_0090 + IT_0111 + IT_0122
       + IT_0124;
    const ccomplex_t IT_0126 = 2*s_13;
    const ccomplex_t IT_0127 = m_C_1*m_tau;
    const ccomplex_t IT_0128 = 2*IT_0127;
    const ccomplex_t IT_0129 = conj(IT_0086) + conj(IT_0123);
    const ccomplex_t IT_0130 = (-0.5)*IT_0032 + (-0.5)*IT_0046;
    const ccomplex_t IT_0131 = conj(IT_0017) + conj(IT_0130);
    const ccomplex_t IT_0132 = (-2)*s_13;
    const ccomplex_t IT_0133 = (-2)*IT_0127;
    const ccomplex_t IT_0134 = IT_0049*IT_0081;
    const ccomplex_t IT_0135 = IT_0110*IT_0120;
    const ccomplex_t IT_0136 = IT_0051*conj(IT_0079);
    const ccomplex_t IT_0137 = conj(IT_0049)*IT_0081;
    const ccomplex_t IT_0138 = IT_0110*conj(IT_0120);
    const ccomplex_t IT_0139 = conj(IT_0109)*IT_0121;
    const ccomplex_t IT_0140 = s_13*IT_0087;
    const ccomplex_t IT_0141 = s_14*s_34;
    const ccomplex_t IT_0142 = (-2)*IT_0141;
    const ccomplex_t IT_0143 = IT_0140 + IT_0142;
    const ccomplex_t IT_0144 = (-2)*IT_0143;
    const ccomplex_t IT_0145 = IT_0049*IT_0051;
    const ccomplex_t IT_0146 = IT_0120*IT_0121;
    const ccomplex_t IT_0147 = IT_0017*IT_0144;
    const ccomplex_t IT_0148 = IT_0017*IT_0125 + IT_0120*(conj(IT_0120)
      *IT_0126 + conj(IT_0109)*IT_0128 + conj(IT_0079)*IT_0132 + conj(IT_0049)
      *IT_0133) + IT_0109*(conj(IT_0109)*IT_0126 + conj(IT_0120)*IT_0128 +
       IT_0121*IT_0129 + IT_0110*IT_0131 + conj(IT_0049)*IT_0132 + conj(IT_0079)
      *IT_0133) + IT_0079*(conj(IT_0079)*IT_0126 + conj(IT_0049)*IT_0128 +
       IT_0051*IT_0129 + IT_0081*IT_0131 + conj(IT_0120)*IT_0132 + conj(IT_0109)
      *IT_0133) + IT_0049*(conj(IT_0049)*IT_0126 + conj(IT_0079)*IT_0128 + conj
      (IT_0109)*IT_0132 + conj(IT_0120)*IT_0133) + IT_0129*(IT_0134 + IT_0135) +
       (IT_0086 + IT_0123)*(IT_0089*IT_0131 + IT_0136 + IT_0137 + IT_0138 +
       IT_0139 + IT_0129*IT_0144) + IT_0130*(IT_0125 + IT_0131*IT_0144) +
       IT_0131*(IT_0145 + IT_0146 + IT_0147);
    return create_ccomplex_return(IT_0148);
}

