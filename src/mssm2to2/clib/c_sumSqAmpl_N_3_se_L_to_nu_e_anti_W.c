#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_se_L_to_nu_e_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_se_L_to_nu_e_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -2);
    const ccomplex_t IT_0001 = pow(m_W, 2);
    const ccomplex_t IT_0002 = s_13*IT_0001;
    const ccomplex_t IT_0003 = s_14*s_34;
    const ccomplex_t IT_0004 = (-2)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = s_13 + IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0010;
    const ccomplex_t IT_0012 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + 0.5*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = pow(m_se_L, 2);
    const ccomplex_t IT_0021 = cpow(s_23 + (-0.5)*IT_0020 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = m_N_3*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0010;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = cos(beta);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0010*IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0027*IT_0032;
    const ccomplex_t IT_0034 = pow(m_e, 2);
    const ccomplex_t IT_0035 = pow(m_N_3, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + -IT_0020 + IT_0034 + -IT_0035 
      + -reg_prop, -1);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_e;
    const ccomplex_t IT_0038 = IT_0033*IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0010;
    const ccomplex_t IT_0040 = N_W3*e_em*conj(V_Wp1);
    const ccomplex_t IT_0041 = IT_0010*IT_0040;
    const ccomplex_t IT_0042 = N_u3*conj(V_u1)*e_em;
    const ccomplex_t IT_0043 = IT_0010*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0041 + (-0.5)*IT_0044);
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = IT_0039*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0049 = cpow(s_23 + (-0.5)*IT_0020 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0047*IT_0048*IT_0049;
    const ccomplex_t IT_0051 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0052 = IT_0010*IT_0051;
    const ccomplex_t IT_0053 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0054 = IT_0010*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0052 + (-0.5)*IT_0055);
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0011*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0060 = IT_0021*IT_0058*IT_0059;
    const ccomplex_t IT_0061 = conj(N_W3)*e_em*U_Wm1;
    const ccomplex_t IT_0062 = IT_0010*IT_0061;
    const ccomplex_t IT_0063 = conj(N_d3)*U_d1*e_em;
    const ccomplex_t IT_0064 = IT_0010*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + 0.5*IT_0065);
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = IT_0039*IT_0067;
    const ccomplex_t IT_0069 = IT_0049*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = m_N_3*IT_0070;
    const ccomplex_t IT_0072 = -IT_0038 + 0.5*IT_0050 + 0.5*IT_0060 + (-0.5)
      *IT_0071;
    const ccomplex_t IT_0073 = IT_0025 + IT_0072;
    const ccomplex_t IT_0074 = conj(IT_0025) + conj(IT_0072);
    const ccomplex_t IT_0075 = pow(m_W, 4);
    const ccomplex_t IT_0076 = s_13*IT_0075;
    const ccomplex_t IT_0077 = IT_0000*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0002 + IT_0004 + IT_0078;
    const ccomplex_t IT_0080 = cos(theta_W);
    const ccomplex_t IT_0081 = cpow(IT_0080, -1);
    const ccomplex_t IT_0082 = conj(N_B3)*e_em;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = conj(N_W3)*e_em;
    const ccomplex_t IT_0086 = IT_0010*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0084 + IT_0087);
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0027*IT_0089;
    const ccomplex_t IT_0091 = IT_0036*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = 0.5*IT_0023;
    const ccomplex_t IT_0094 = (-0.5)*IT_0070;
    const ccomplex_t IT_0095 = IT_0093 + -IT_0094;
    const ccomplex_t IT_0096 = -conj(IT_0094);
    const ccomplex_t IT_0097 = conj(IT_0093) + IT_0096;
    const ccomplex_t IT_0098 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = s_13 + IT_0099;
    const ccomplex_t IT_0101 = s_13*IT_0100;
    const ccomplex_t IT_0102 = -IT_0026;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0084 + -IT_0087);
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = cpow(s_13 + (-0.5)*IT_0035 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0102*IT_0107;
    const ccomplex_t IT_0109 = 2*IT_0092 + (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = -IT_0023;
    const ccomplex_t IT_0111 = 2*IT_0110;
    const ccomplex_t IT_0112 = IT_0070 + 0.5*IT_0108;
    const ccomplex_t IT_0113 = (-2)*IT_0112;
    const ccomplex_t IT_0114 = IT_0111 + IT_0113;
    const ccomplex_t IT_0115 = 2*conj(IT_0110);
    const ccomplex_t IT_0116 = (-2)*conj(IT_0112);
    const ccomplex_t IT_0117 = IT_0115 + IT_0116;
    const ccomplex_t IT_0118 = m_N_3*IT_0100;
    const ccomplex_t IT_0119 = s_34*IT_0001;
    const ccomplex_t IT_0120 = IT_0000*IT_0119;
    const ccomplex_t IT_0121 = (-0.25)*IT_0120;
    const ccomplex_t IT_0122 = s_34 + IT_0121;
    const ccomplex_t IT_0123 = m_N_3*IT_0122;
    const ccomplex_t IT_0124 = (-4)*IT_0123;
    const ccomplex_t IT_0125 = 2*IT_0025;
    const ccomplex_t IT_0126 = 2*IT_0072;
    const ccomplex_t IT_0127 = IT_0125 + IT_0126;
    const ccomplex_t IT_0128 = 2*conj(IT_0025);
    const ccomplex_t IT_0129 = 2*conj(IT_0072);
    const ccomplex_t IT_0130 = IT_0128 + IT_0129;
    const ccomplex_t IT_0131 = pow(s_34, 2);
    const ccomplex_t IT_0132 = s_13*IT_0000*IT_0131;
    const ccomplex_t IT_0133 = 2*IT_0109;
    const ccomplex_t IT_0134 = m_N_3*IT_0000*IT_0131;
    const ccomplex_t IT_0135 = 2*conj(IT_0109);
    const ccomplex_t IT_0136 = s_13*s_34;
    const ccomplex_t IT_0137 = conj(IT_0092) + conj(IT_0093);
    const ccomplex_t IT_0138 = IT_0092 + IT_0093;
    const ccomplex_t IT_0139 = (-0.25)*IT_0077;
    const ccomplex_t IT_0140 = IT_0003 + IT_0139;
    const ccomplex_t IT_0141 = 0.5*IT_0120;
    const ccomplex_t IT_0142 = s_34 + IT_0141;
    const ccomplex_t IT_0143 = m_N_3*IT_0142;
    const ccomplex_t IT_0144 = 2*IT_0143;
    const ccomplex_t IT_0145 = 0.125*IT_0144;
    const ccomplex_t IT_0146 = pow(s_14, 2);
    const ccomplex_t IT_0147 = IT_0000*IT_0146;
    const ccomplex_t IT_0148 = -IT_0035;
    const ccomplex_t IT_0149 = IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = s_13*IT_0149;
    const ccomplex_t IT_0151 = (-2)*conj(IT_0110);
    const ccomplex_t IT_0152 = 2*IT_0112;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = s_13*s_14;
    const ccomplex_t IT_0155 = s_34*IT_0035;
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = 4*IT_0008*IT_0073*IT_0074 + (-4)*IT_0079*
      (IT_0092*conj(IT_0092) + IT_0095*IT_0097) + IT_0101*(conj(IT_0109)*IT_0114
       + IT_0109*IT_0117) + (IT_0074*IT_0114 + IT_0073*IT_0117)*IT_0118 + (-2)
      *IT_0124*(IT_0074*IT_0093 + IT_0073*conj(IT_0093) + (-0.5)*conj(IT_0094)
      *IT_0127 + (-0.5)*IT_0094*IT_0130) + conj(IT_0109)*IT_0132*IT_0133 +
       IT_0134*(IT_0074*IT_0133 + IT_0073*IT_0135) + (-2)*IT_0136*(conj(IT_0094)
      *IT_0109 + IT_0094*conj(IT_0109) + (-0.5)*IT_0133*IT_0137 + (-0.5)*IT_0135
      *IT_0138) + 8*conj(IT_0092)*(IT_0095*IT_0140 + IT_0127*IT_0145) + 8
      *IT_0092*(IT_0097*IT_0140 + IT_0130*IT_0145) + IT_0150*(conj(IT_0110)
      *IT_0111 + IT_0112*IT_0151 + (-2)*conj(IT_0112)*(IT_0110 + IT_0153)) + 
      (IT_0114*IT_0137 + IT_0117*IT_0138 + 2*IT_0094*(conj(IT_0112) + 0.5
      *IT_0151) + 2*IT_0096*(IT_0110 + IT_0153))*IT_0157;
    return create_ccomplex_return(IT_0158);
}

