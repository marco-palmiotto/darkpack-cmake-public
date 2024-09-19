#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_ss_L_to_mu_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_ss_L_to_mu_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0008 = (-0.5)*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_mu, 2);
    const ccomplex_t IT_0011 = pow(m_smu_R, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = -IT_0014 + -IT_0022 + -IT_0030 + -IT_0038;
    const ccomplex_t IT_0040 = cos(theta_W);
    const ccomplex_t IT_0041 = cpow(IT_0040, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = conj(N_B1)*e_em;
    const ccomplex_t IT_0045 = IT_0041*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W1)*e_em;
    const ccomplex_t IT_0048 = IT_0004*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + (-3)*IT_0049);
    const ccomplex_t IT_0051 = 0.166666666666667*IT_0050;
    const ccomplex_t IT_0052 = IT_0043*IT_0051;
    const ccomplex_t IT_0053 = IT_0012*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B2)*e_em;
    const ccomplex_t IT_0056 = IT_0041*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W2)*e_em;
    const ccomplex_t IT_0059 = IT_0004*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + (-3)*IT_0060);
    const ccomplex_t IT_0062 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0041;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = IT_0020*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = conj(N_B3)*e_em;
    const ccomplex_t IT_0069 = IT_0041*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = conj(N_W3)*e_em;
    const ccomplex_t IT_0072 = IT_0004*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + (-3)*IT_0073);
    const ccomplex_t IT_0075 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0041;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0075*IT_0077;
    const ccomplex_t IT_0079 = IT_0028*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = conj(N_B4)*e_em;
    const ccomplex_t IT_0082 = IT_0041*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W4)*e_em;
    const ccomplex_t IT_0085 = IT_0004*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + (-3)*IT_0086);
    const ccomplex_t IT_0088 = 0.166666666666667*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0041;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = IT_0036*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = -IT_0054 + -IT_0067 + -IT_0080 + -IT_0093;
    const ccomplex_t IT_0095 = pow(m_ss_L, 2);
    const ccomplex_t IT_0096 = m_s*m_mu*IT_0095;
    const ccomplex_t IT_0097 = s_34*IT_0095;
    const ccomplex_t IT_0098 = (-6)*IT_0097;
    const ccomplex_t IT_0099 = s_23*s_24;
    const ccomplex_t IT_0100 = 12*IT_0099;
    const ccomplex_t IT_0101 = IT_0098 + IT_0100;
    const ccomplex_t IT_0102 = m_s*s_23;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0104 = IT_0016*IT_0062;
    const ccomplex_t IT_0105 = IT_0020*IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0107 = IT_0024*IT_0075;
    const ccomplex_t IT_0108 = IT_0028*IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0110 = IT_0032*IT_0088;
    const ccomplex_t IT_0111 = IT_0036*IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0113 = IT_0006*IT_0051;
    const ccomplex_t IT_0114 = IT_0012*IT_0112*IT_0113;
    const ccomplex_t IT_0115 = -IT_0105 + -IT_0108 + -IT_0111 + -IT_0114;
    const ccomplex_t IT_0116 = IT_0014 + IT_0022 + IT_0030 + IT_0038;
    const ccomplex_t IT_0117 = m_s*IT_0116;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = IT_0115 + -IT_0118;
    const ccomplex_t IT_0120 = -conj(IT_0118);
    const ccomplex_t IT_0121 = conj(IT_0115) + IT_0120;
    const ccomplex_t IT_0122 = IT_0008*IT_0043;
    const ccomplex_t IT_0123 = IT_0012*IT_0112*IT_0122;
    const ccomplex_t IT_0124 = IT_0018*IT_0064;
    const ccomplex_t IT_0125 = IT_0020*IT_0103*IT_0124;
    const ccomplex_t IT_0126 = IT_0026*IT_0077;
    const ccomplex_t IT_0127 = IT_0028*IT_0106*IT_0126;
    const ccomplex_t IT_0128 = IT_0034*IT_0090;
    const ccomplex_t IT_0129 = IT_0036*IT_0109*IT_0128;
    const ccomplex_t IT_0130 = -IT_0123 + -IT_0125 + -IT_0127 + -IT_0129;
    const ccomplex_t IT_0131 = 6*IT_0130;
    const ccomplex_t IT_0132 = IT_0054 + IT_0067 + IT_0080 + IT_0093;
    const ccomplex_t IT_0133 = m_s*IT_0132;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = (-6)*IT_0134;
    const ccomplex_t IT_0136 = IT_0131 + IT_0135;
    const ccomplex_t IT_0137 = 6*conj(IT_0130);
    const ccomplex_t IT_0138 = (-6)*conj(IT_0134);
    const ccomplex_t IT_0139 = IT_0137 + IT_0138;
    const ccomplex_t IT_0140 = m_s*m_mu;
    const ccomplex_t IT_0141 = 6*IT_0134;
    const ccomplex_t IT_0142 = (-6)*IT_0130;
    const ccomplex_t IT_0143 = IT_0141 + IT_0142;
    const ccomplex_t IT_0144 = 6*conj(IT_0134);
    const ccomplex_t IT_0145 = (-6)*conj(IT_0130);
    const ccomplex_t IT_0146 = IT_0144 + IT_0145;
    const ccomplex_t IT_0147 = m_mu*s_24;
    const ccomplex_t IT_0148 = 6*s_34;
    const ccomplex_t IT_0149 = (-6)*s_34;
    const ccomplex_t IT_0150 = (-6)*conj(IT_0039)*(IT_0094*IT_0096 + (
      -0.166666666666667)*IT_0039*IT_0101) + (-6)*conj(IT_0094)*(IT_0039*IT_0096
       + (-0.166666666666667)*IT_0094*IT_0101) + 6*IT_0102*(conj(IT_0039)
      *IT_0119 + IT_0039*IT_0121 + 0.166666666666667*conj(IT_0094)*IT_0136 +
       0.166666666666667*IT_0094*IT_0139) + IT_0140*(conj(IT_0118)*IT_0136 +
       IT_0118*IT_0139 + conj(IT_0115)*IT_0143 + IT_0115*IT_0146) + (-6)*(conj
      (IT_0094)*IT_0119 + IT_0094*IT_0121 + (-0.166666666666667)*conj(IT_0039)
      *IT_0143 + (-0.166666666666667)*IT_0039*IT_0146)*IT_0147 + (IT_0115*conj
      (IT_0115) + IT_0118*conj(IT_0118) + IT_0130*conj(IT_0130) + IT_0134*conj
      (IT_0134))*IT_0148 + (conj(IT_0115)*IT_0118 + IT_0115*conj(IT_0118) + conj
      (IT_0130)*IT_0134 + IT_0130*conj(IT_0134))*IT_0149;
    return create_ccomplex_return(IT_0150);
}

