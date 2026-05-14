#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_00) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0005 = (-0.333333333333333)*IT_0004;
    const ccomplex_t IT_0006 = m_d*IT_0005;
    const ccomplex_t IT_0007 = (-0.666666666666667)*IT_0004;
    const ccomplex_t IT_0008 = m_d*IT_0007;
    const ccomplex_t IT_0009 = pow(m_sb_1, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0012 = 0.333333333333333*IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0001*IT_0002;
    const ccomplex_t IT_0016 = e_em*IT_0015;
    const ccomplex_t IT_0017 = IT_0000*IT_0003;
    const ccomplex_t IT_0018 = e_em*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0016 + 3*IT_0018);
    const ccomplex_t IT_0020 = (-0.166666666666667)*IT_0019;
    const ccomplex_t IT_0021 = IT_0010*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (IT_0006 + -IT_0008)*(IT_0014 + -IT_0022);
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = cos(alpha);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_d*e_em*IT_0003*IT_0025
      *IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = cpow(V_tb, 2);
    const ccomplex_t IT_0032 = m_b*mu_h*IT_0025;
    const ccomplex_t IT_0033 = sin(alpha);
    const ccomplex_t IT_0034 = sin(beta);
    const ccomplex_t IT_0035 = -IT_0033*IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0028;
    const ccomplex_t IT_0037 = IT_0035 + IT_0036;
    const ccomplex_t IT_0038 = cpow(IT_0000, -2);
    const ccomplex_t IT_0039 = pow(m_b, 2);
    const ccomplex_t IT_0040 = (-2)*U_sb_00*(((-0.5)*m_W*IT_0037*(IT_0003 + 1.
      /3*IT_0002*IT_0038) + IT_0003*IT_0025*IT_0027*IT_0028*IT_0039)*conj
      (U_sb_00) + -1./2*IT_0003*IT_0027*IT_0032*IT_0033*conj(U_sb_10)) + U_sb_10
      *(IT_0003*IT_0027*IT_0032*IT_0033*conj(U_sb_00) + 0.666666666666667*(m_W
      *IT_0002*IT_0037*IT_0038 + (-3)*IT_0003*IT_0025*IT_0027*IT_0028*IT_0039)
      *conj(U_sb_10));
    const ccomplex_t IT_0041 = (0 + _Complex_I*-1)*(e_em*IT_0040 + (
      -1.4142135623731)*A_b*IT_0028*IT_0031*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0030*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + (-2)*IT_0009 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*m_d*e_em*IT_0003*IT_0025
      *IT_0027*IT_0033;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0026*IT_0033;
    const ccomplex_t IT_0050 = IT_0028*IT_0034;
    const ccomplex_t IT_0051 = IT_0049 + IT_0050;
    const ccomplex_t IT_0052 = 2*U_sb_00*((IT_0003*IT_0025*IT_0027*IT_0033
      *IT_0039 + (-0.5)*m_W*(IT_0003 + 1./3*IT_0002*IT_0038)*IT_0051)*conj
      (U_sb_00) + 1./2*IT_0003*IT_0027*IT_0028*IT_0032*conj(U_sb_10)) + U_sb_10*
      (IT_0003*IT_0027*IT_0028*IT_0032*conj(U_sb_00) + 2*(IT_0003*IT_0025
      *IT_0027*IT_0033*IT_0039 + -1./3*m_W*IT_0002*IT_0038*IT_0051)*conj(U_sb_10
      ));
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(e_em*IT_0052 +
       1.4142135623731*A_b*IT_0031*IT_0033*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0048*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = -IT_0046 + -IT_0058;
    const ccomplex_t IT_0060 = -IT_0024;
    const ccomplex_t IT_0061 = pow(m_d, 2);
    const ccomplex_t IT_0062 = (-18)*IT_0061;
    const ccomplex_t IT_0063 = 18*s_34;
    const ccomplex_t IT_0064 = IT_0062 + IT_0063;
    const ccomplex_t IT_0065 = m_d*e_em*IT_0003*IT_0025*IT_0027*IT_0034;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = m_b*e_em*mu_h*IT_0025*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0068 = IT_0003*IT_0067;
    const ccomplex_t IT_0069 = m_b*e_em*mu_h*IT_0025*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0070 = IT_0003*IT_0069;
    const ccomplex_t IT_0071 = A_b*IT_0031*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0072 = IT_0034*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = A_b*IT_0031*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0075 = IT_0034*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0068 + 0.5*IT_0070 + (-0.5)*IT_0073 +
       0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0066*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + (-2)*IT_0009 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = m_d*e_em*IT_0003*IT_0025;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0003*IT_0027*IT_0034;
    const ccomplex_t IT_0085 = IT_0067*IT_0084;
    const ccomplex_t IT_0086 = IT_0069*IT_0084;
    const ccomplex_t IT_0087 = IT_0026*IT_0071;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = IT_0026*IT_0074;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0085 + 0.5*IT_0086 + 0.5*IT_0088 + (
      -0.5)*IT_0090;
    const ccomplex_t IT_0092 = IT_0083*IT_0091;
    const ccomplex_t IT_0093 = pow(m_Z, 2);
    const ccomplex_t IT_0094 = cpow((-2)*s_12 + (-2)*IT_0009 + IT_0093 + 
      -reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = -IT_0081 + -IT_0096;
    const ccomplex_t IT_0098 = 18*IT_0061 + IT_0063;
    const ccomplex_t IT_0099 = 36*s_34;
    const ccomplex_t IT_0100 = (-18)*s_34 + IT_0062;
    const ccomplex_t IT_0101 = cpow(s_12 + IT_0009 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (-2)*IT_0102;
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0108 = 0.666666666666667*IT_0107;
    const ccomplex_t IT_0109 = (-0.333333333333333)*IT_0107;
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = IT_0101*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = (-0.0833333333333333)*IT_0106 + 0.5*IT_0112;
    const ccomplex_t IT_0114 = -s_14;
    const ccomplex_t IT_0115 = s_13 + IT_0114;
    const ccomplex_t IT_0116 = m_d*IT_0115;
    const ccomplex_t IT_0117 = 18*IT_0116;
    const ccomplex_t IT_0118 = (-18)*IT_0024;
    const ccomplex_t IT_0119 = 0.25*IT_0106;
    const ccomplex_t IT_0120 = (-3)*IT_0113*(conj(IT_0024) + (-2)*IT_0059 + 2
      *conj(IT_0059) + -IT_0060 + 0.0555555555555556*IT_0118) + (IT_0024 + -conj
      (IT_0024) + 2*IT_0059 + (-2)*conj(IT_0059) + IT_0060 + -conj(IT_0060))
      *IT_0119;
    const ccomplex_t IT_0121 = (-18)*conj(IT_0024);
    const ccomplex_t IT_0122 = cpow(IT_0113, 2);
    const ccomplex_t IT_0123 = s_34*IT_0009;
    const ccomplex_t IT_0124 = (-36)*IT_0123;
    const ccomplex_t IT_0125 = IT_0009*IT_0061;
    const ccomplex_t IT_0126 = (-36)*IT_0125;
    const ccomplex_t IT_0127 = s_13*s_14;
    const ccomplex_t IT_0128 = 72*IT_0127;
    const ccomplex_t IT_0129 = IT_0124 + IT_0126 + IT_0128;
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0007*IT_0014;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = m_d*s_13;
    const ccomplex_t IT_0134 = (-18)*IT_0123;
    const ccomplex_t IT_0135 = 36*IT_0127;
    const ccomplex_t IT_0136 = IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = IT_0007*IT_0022;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = (-18)*IT_0125;
    const ccomplex_t IT_0140 = s_13 + s_14;
    const ccomplex_t IT_0141 = m_d*IT_0140;
    const ccomplex_t IT_0142 = 18*IT_0141;
    const ccomplex_t IT_0143 = conj(IT_0059)*IT_0117;
    const ccomplex_t IT_0144 = m_d*s_14;
    const ccomplex_t IT_0145 = (-18)*IT_0144;
    const ccomplex_t IT_0146 = IT_0134 + IT_0135 + IT_0139;
    const ccomplex_t IT_0147 = -IT_0113*IT_0146;
    const ccomplex_t IT_0148 = (-6)*IT_0123;
    const ccomplex_t IT_0149 = (-6)*IT_0125;
    const ccomplex_t IT_0150 = 12*IT_0127;
    const ccomplex_t IT_0151 = IT_0148 + IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = -IT_0119*IT_0151;
    const ccomplex_t IT_0153 = 18*IT_0133;
    const ccomplex_t IT_0154 = (-18)*IT_0141;
    const ccomplex_t IT_0155 = IT_0113*IT_0146;
    const ccomplex_t IT_0156 = IT_0119*IT_0151;
    const ccomplex_t IT_0157 = IT_0059*IT_0117;
    const ccomplex_t IT_0158 = (conj(IT_0024)*IT_0059 + conj(IT_0059)*(IT_0024
       + IT_0060))*IT_0064 + conj(IT_0024)*(IT_0024*IT_0063 + IT_0097*IT_0098) +
       IT_0059*conj(IT_0059)*((-36)*IT_0061 + IT_0099) + conj(IT_0097)*(IT_0024
      *IT_0098 + IT_0097*(36*IT_0061 + IT_0099) + IT_0060*IT_0100) + conj
      (IT_0060)*(IT_0060*IT_0063 + IT_0059*IT_0064 + IT_0097*IT_0100 + IT_0113
      *IT_0117 + IT_0061*IT_0118) + (-6)*IT_0116*IT_0120 + IT_0060*IT_0061
      *IT_0121 + IT_0122*IT_0130 + (-2)*IT_0119*(IT_0113*((-12)*IT_0123 + (-12)
      *IT_0125 + 24*IT_0127) + (-0.5)*IT_0119*IT_0130) + IT_0132*(18*conj
      (IT_0024)*IT_0133 + conj(IT_0132)*IT_0136 + conj(IT_0138)*IT_0139 + conj
      (IT_0097)*IT_0142 + IT_0143 + conj(IT_0060)*IT_0145 + IT_0147 + IT_0152) +
       IT_0138*(IT_0136*conj(IT_0138) + conj(IT_0132)*IT_0139 + IT_0143 +
       IT_0121*IT_0144 + IT_0147 + IT_0152 + conj(IT_0060)*IT_0153 + conj
      (IT_0097)*IT_0154) + conj(IT_0132)*(18*IT_0024*IT_0133 + IT_0097*IT_0142 +
       IT_0060*IT_0145 + IT_0155 + IT_0156 + IT_0157) + conj(IT_0138)*(IT_0118
      *IT_0144 + IT_0060*IT_0153 + IT_0097*IT_0154 + IT_0155 + IT_0156 + IT_0157);
    return create_ccomplex_return(IT_0158);
}

