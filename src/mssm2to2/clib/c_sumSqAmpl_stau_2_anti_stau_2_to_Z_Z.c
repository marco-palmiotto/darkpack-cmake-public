#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_Z_Z(
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
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = U_stau_01*conj(U_stau_01);
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0001, -4);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = U_stau_11*conj(U_stau_11);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0006, -2);
    const ccomplex_t IT_0013 = IT_0003*IT_0012;
    const ccomplex_t IT_0014 = tan(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, 2);
    const ccomplex_t IT_0016 = cpow(1 + IT_0015, -1);
    const ccomplex_t IT_0017 = (IT_0004 + (-0.5)*IT_0009 + (-2)*IT_0011 + (
      -0.5)*IT_0013)*IT_0016;
    const ccomplex_t IT_0018 = IT_0000*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = m_W*e_em;
    const ccomplex_t IT_0022 = cos(alpha);
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = IT_0005*IT_0006*IT_0022*IT_0023;
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = sin(alpha);
    const ccomplex_t IT_0027 = IT_0005*IT_0006*IT_0025*IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0006, -3);
    const ccomplex_t IT_0029 = IT_0022*IT_0023*IT_0028;
    const ccomplex_t IT_0030 = IT_0025*IT_0026*IT_0028;
    const ccomplex_t IT_0031 = cpow(IT_0006, -1);
    const ccomplex_t IT_0032 = IT_0002*IT_0022*IT_0023*IT_0031;
    const ccomplex_t IT_0033 = IT_0002*IT_0025*IT_0026*IT_0031;
    const ccomplex_t IT_0034 = cpow(IT_0002 + IT_0012, -1);
    const ccomplex_t IT_0035 = (IT_0024 + -IT_0027 + IT_0029 + -IT_0030 + 2
      *IT_0032 + (-2)*IT_0033)*IT_0034;
    const ccomplex_t IT_0036 = IT_0021*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = pow(m_W, -1);
    const ccomplex_t IT_0039 = mu_h*m_tau*IT_0038;
    const ccomplex_t IT_0040 = cpow(IT_0025, -1);
    const ccomplex_t IT_0041 = IT_0025*IT_0026;
    const ccomplex_t IT_0042 = IT_0022*IT_0023;
    const ccomplex_t IT_0043 = IT_0041 + IT_0042;
    const ccomplex_t IT_0044 = pow(m_tau, 2);
    const ccomplex_t IT_0045 = 2*U_stau_01*((0.5*m_W*(IT_0002*IT_0006 + 
      -IT_0031)*IT_0043 + IT_0026*IT_0031*IT_0038*IT_0040*IT_0044)*conj
      (U_stau_01) + 1./2*IT_0022*IT_0031*IT_0039*IT_0040*conj(U_stau_11)) +
       U_stau_11*(IT_0022*IT_0031*IT_0039*IT_0040*conj(U_stau_01) + (-2)*(m_W
      *IT_0002*IT_0006*IT_0043 + -IT_0026*IT_0031*IT_0038*IT_0040*IT_0044)*conj
      (U_stau_11));
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(e_em*IT_0045 +
       1.4142135623731*A_tau*IT_0026*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0037*IT_0047;
    const ccomplex_t IT_0049 = pow(m_stau_2, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0049 + -reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0005*IT_0006*IT_0022*IT_0025;
    const ccomplex_t IT_0054 = IT_0005*IT_0006*IT_0023*IT_0026;
    const ccomplex_t IT_0055 = IT_0022*IT_0025*IT_0028;
    const ccomplex_t IT_0056 = IT_0023*IT_0026*IT_0028;
    const ccomplex_t IT_0057 = IT_0002*IT_0022*IT_0025*IT_0031;
    const ccomplex_t IT_0058 = IT_0002*IT_0023*IT_0026*IT_0031;
    const ccomplex_t IT_0059 = IT_0034*(IT_0053 + IT_0054 + IT_0055 + IT_0056 
      + 2*IT_0057 + 2*IT_0058);
    const ccomplex_t IT_0060 = IT_0021*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = -IT_0023*IT_0026;
    const ccomplex_t IT_0063 = IT_0022*IT_0025;
    const ccomplex_t IT_0064 = IT_0062 + IT_0063;
    const ccomplex_t IT_0065 = (-2)*U_stau_01*((IT_0022*IT_0031*IT_0038
      *IT_0040*IT_0044 + 0.5*m_W*(IT_0002*IT_0006 + -IT_0031)*IT_0064)*conj
      (U_stau_01) + -1./2*IT_0026*IT_0031*IT_0039*IT_0040*conj(U_stau_11)) +
       U_stau_11*(IT_0026*IT_0031*IT_0039*IT_0040*conj(U_stau_01) + (-2)*
      (IT_0022*IT_0031*IT_0038*IT_0040*IT_0044 + -m_W*IT_0002*IT_0006*IT_0064)
      *conj(U_stau_11));
    const ccomplex_t IT_0066 = (0 + _Complex_I*-1)*(e_em*IT_0065 + (
      -1.4142135623731)*A_tau*IT_0022*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0061*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + (-2)*IT_0049 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0020 + -IT_0052 + -IT_0071;
    const ccomplex_t IT_0073 = pow(m_Z, -4);
    const ccomplex_t IT_0074 = pow(s_34, 2);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = pow(m_Z, 2);
    const ccomplex_t IT_0077 = cpow((-2)*s_13 + IT_0049 + IT_0076 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = cpow(IT_0001, -1);
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0031 + -IT_0006*IT_0079)*U_stau_00*conj(U_stau_01) + IT_0006*IT_0079
      *U_stau_10*conj(U_stau_11));
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0031 + -IT_0006*IT_0079)*conj(U_stau_00)*U_stau_01 + IT_0006*IT_0079
      *conj(U_stau_10)*U_stau_11);
    const ccomplex_t IT_0083 = (-2)*IT_0082;
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = -IT_0080;
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = IT_0078*IT_0088;
    const ccomplex_t IT_0090 = cpow((-2)*s_13 + IT_0076 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0031 + -IT_0006*IT_0079)*U_stau_01*conj(U_stau_01) + IT_0006*IT_0079
      *U_stau_11*conj(U_stau_11));
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = (-2)*IT_0092;
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = cpow(IT_0094, 2);
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = IT_0095 + IT_0097;
    const ccomplex_t IT_0099 = IT_0091*IT_0098;
    const ccomplex_t IT_0100 = -IT_0089 + IT_0099;
    const ccomplex_t IT_0101 = s_13*s_24*s_34*IT_0073;
    const ccomplex_t IT_0102 = pow(m_Z, -2);
    const ccomplex_t IT_0103 = s_14*s_24*IT_0102;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = s_13*s_23*IT_0102;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = s_12 + IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = IT_0101 + IT_0107;
    const ccomplex_t IT_0109 = -IT_0082;
    const ccomplex_t IT_0110 = IT_0081*IT_0109;
    const ccomplex_t IT_0111 = cpow((-2)*s_23 + IT_0049 + IT_0076 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0110*IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_23 + IT_0076 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = IT_0095*IT_0115;
    const ccomplex_t IT_0117 = -IT_0113 + -IT_0116;
    const ccomplex_t IT_0118 = s_14*s_23*s_34*IT_0073;
    const ccomplex_t IT_0119 = IT_0107 + IT_0118;
    const ccomplex_t IT_0120 = -IT_0084;
    const ccomplex_t IT_0121 = IT_0110 + IT_0120;
    const ccomplex_t IT_0122 = IT_0112*IT_0121;
    const ccomplex_t IT_0123 = IT_0098*IT_0115;
    const ccomplex_t IT_0124 = -IT_0122 + -IT_0123;
    const ccomplex_t IT_0125 = s_23*s_24*s_34*IT_0073;
    const ccomplex_t IT_0126 = pow(s_24, 2);
    const ccomplex_t IT_0127 = IT_0102*IT_0126;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = pow(s_23, 2);
    const ccomplex_t IT_0130 = IT_0102*IT_0129;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = IT_0049 + IT_0125 + IT_0128 + IT_0131;
    const ccomplex_t IT_0133 = pow(s_13, 2);
    const ccomplex_t IT_0134 = IT_0102*IT_0133;
    const ccomplex_t IT_0135 = -(s_12 + -IT_0103)*(IT_0049 + -IT_0134);
    const ccomplex_t IT_0136 = IT_0078*IT_0086;
    const ccomplex_t IT_0137 = IT_0091*IT_0095;
    const ccomplex_t IT_0138 = -IT_0136 + -IT_0137;
    const ccomplex_t IT_0139 = s_14*s_34*IT_0102;
    const ccomplex_t IT_0140 = -(IT_0049 + -IT_0134)*(s_13 + -IT_0139);
    const ccomplex_t IT_0141 = pow(s_14, 2);
    const ccomplex_t IT_0142 = IT_0102*IT_0141;
    const ccomplex_t IT_0143 = -(s_12 + -IT_0105)*(IT_0049 + -IT_0142);
    const ccomplex_t IT_0144 = IT_0136 + IT_0137;
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = (s_12 + -IT_0103)*(s_12 + -IT_0105);
    const ccomplex_t IT_0147 = -(IT_0049 + -IT_0127)*(IT_0049 + -IT_0134);
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = IT_0122 + IT_0123;
    const ccomplex_t IT_0150 = s_24*s_34*IT_0102;
    const ccomplex_t IT_0151 = (s_12 + -IT_0105)*(s_23 + -IT_0150);
    const ccomplex_t IT_0152 = -conj(IT_0144);
    const ccomplex_t IT_0153 = (s_12 + -IT_0105)*(s_13 + -IT_0139);
    const ccomplex_t IT_0154 = -(IT_0049 + -IT_0130)*(IT_0049 + -IT_0142);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = -(IT_0049 + -IT_0134)*(s_23 + -IT_0150);
    const ccomplex_t IT_0157 = (s_12 + -IT_0105)*(IT_0049 + -IT_0127);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = -IT_0100;
    const ccomplex_t IT_0160 = s_34*IT_0102;
    const ccomplex_t IT_0161 = s_13*s_34*IT_0102;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = s_14 + IT_0162;
    const ccomplex_t IT_0164 = IT_0160*IT_0163;
    const ccomplex_t IT_0165 = s_23*s_34*IT_0102;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = s_24 + IT_0166;
    const ccomplex_t IT_0168 = IT_0160*IT_0167;
    const ccomplex_t IT_0169 = -IT_0049;
    const ccomplex_t IT_0170 = s_13*s_14*s_34*IT_0073;
    const ccomplex_t IT_0171 = -IT_0170;
    const ccomplex_t IT_0172 = IT_0134 + IT_0142 + IT_0169 + IT_0171;
    const ccomplex_t IT_0173 = -IT_0072;
    const ccomplex_t IT_0174 = IT_0074*IT_0102;
    const ccomplex_t IT_0175 = -(IT_0049 + -IT_0134)*(IT_0076 + -IT_0174);
    const ccomplex_t IT_0176 = (s_12 + -IT_0105)*(IT_0076 + -IT_0174);
    const ccomplex_t IT_0177 = -conj(IT_0100);
    const ccomplex_t IT_0178 = -conj(IT_0072);
    const ccomplex_t IT_0179 = (IT_0134 + IT_0169)*(IT_0142 + IT_0169);
    const ccomplex_t IT_0180 = (s_12 + -IT_0103)*(IT_0049 + -IT_0130);
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = (IT_0049 + -IT_0130)*(s_13 + -IT_0139);
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = -IT_0117;
    const ccomplex_t IT_0185 = (IT_0049 + -IT_0127)*(IT_0049 + -IT_0130);
    const ccomplex_t IT_0186 = -conj(IT_0117);
    const ccomplex_t IT_0187 = (IT_0049 + -IT_0130)*(s_23 + -IT_0150);
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = -conj(IT_0149);
    const ccomplex_t IT_0190 = (IT_0049 + -IT_0130)*(IT_0076 + -IT_0174);
    const ccomplex_t IT_0191 = -IT_0190;
    const ccomplex_t IT_0192 = 2*IT_0072*(conj(IT_0072)*(1 + 0.5*IT_0075) +
       0.5*conj(IT_0100)*IT_0108 + 0.5*conj(IT_0117)*IT_0119 + 0.5*conj(IT_0124)
      *IT_0132) + (conj(IT_0100)*IT_0135 + conj(IT_0138)*IT_0140 + conj(IT_0117)
      *IT_0143)*IT_0145 + IT_0100*(conj(IT_0072)*IT_0108 + conj(IT_0117)*IT_0146
       + conj(IT_0100)*IT_0148 + conj(IT_0149)*IT_0151 + IT_0135*IT_0152) +
       IT_0117*(conj(IT_0072)*IT_0119 + conj(IT_0100)*IT_0146 + IT_0143*IT_0152 
      + conj(IT_0138)*IT_0153 + conj(IT_0117)*IT_0155) + (conj(IT_0138)*IT_0156 
      + conj(IT_0124)*IT_0158)*IT_0159 + (conj(IT_0138)*IT_0164 + conj(IT_0149)
      *IT_0168 + conj(IT_0144)*IT_0172)*IT_0173 + IT_0138*(conj(IT_0124)*IT_0151
       + IT_0140*IT_0152 + conj(IT_0117)*IT_0153 + -conj(IT_0138)*IT_0175 + conj
      (IT_0149)*IT_0176 + IT_0156*IT_0177 + IT_0164*IT_0178) + IT_0144*(conj
      (IT_0124)*IT_0146 + conj(IT_0149)*IT_0153 + IT_0172*IT_0178 + conj(IT_0144
      )*IT_0179) + (conj(IT_0124)*IT_0181 + conj(IT_0149)*IT_0183)*IT_0184 +
       IT_0124*(conj(IT_0072)*IT_0132 + conj(IT_0144)*IT_0146 + conj(IT_0138)
      *IT_0151 + IT_0158*IT_0177 + conj(IT_0124)*IT_0185 + IT_0181*IT_0186 +
       IT_0188*IT_0189) + IT_0149*(conj(IT_0100)*IT_0151 + conj(IT_0144)*IT_0153
       + conj(IT_0138)*IT_0176 + IT_0168*IT_0178 + IT_0183*IT_0186 + -conj
      (IT_0124)*IT_0188 + IT_0189*IT_0191);
    return create_ccomplex_return(IT_0192);
}

