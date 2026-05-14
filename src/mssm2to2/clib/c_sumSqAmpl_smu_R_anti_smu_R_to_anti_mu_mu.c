#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_smu_R_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_smu_R_to_anti_mu_mu(
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
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = e_em*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + -IT_0007);
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = pow(m_smu_R, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = m_mu*IT_0015;
    const ccomplex_t IT_0017 = (-2)*IT_0014;
    const ccomplex_t IT_0018 = m_mu*IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = IT_0016 + IT_0019;
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = IT_0011*IT_0014;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0026 = pow(m_mu, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_23 + IT_0010 + IT_0026 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = cos(beta);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0005*IT_0028*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0001;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0027*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0001;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0005*IT_0028*IT_0030;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0010 + IT_0026 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0037*IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0001;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0005*IT_0028*IT_0030;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0010 + IT_0026 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0045*IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0001;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0005*IT_0028*IT_0030;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0010 + IT_0026 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0053*IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0001;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = IT_0034*IT_0062;
    const ccomplex_t IT_0064 = IT_0027*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = m_mu*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0001;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0039*IT_0068;
    const ccomplex_t IT_0070 = IT_0043*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = m_mu*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0001;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = IT_0047*IT_0074;
    const ccomplex_t IT_0076 = IT_0051*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = m_mu*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0001;
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = IT_0055*IT_0080;
    const ccomplex_t IT_0082 = IT_0059*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = m_mu*IT_0083;
    const ccomplex_t IT_0085 = IT_0021 + -IT_0024 + -IT_0036 + -IT_0044 + 
      -IT_0052 + -IT_0060 + IT_0066 + IT_0072 + IT_0078 + IT_0084;
    const ccomplex_t IT_0086 = IT_0017*IT_0023;
    const ccomplex_t IT_0087 = -IT_0065 + -IT_0071 + -IT_0077 + -IT_0083 + 
      -IT_0086;
    const ccomplex_t IT_0088 = m_mu*s_14;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*N_d2*e_em*m_mu
      *IT_0005*IT_0028*IT_0030;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = IT_0068*IT_0090;
    const ccomplex_t IT_0092 = IT_0037*IT_0043*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*N_d3*e_em*m_mu
      *IT_0005*IT_0028*IT_0030;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0074*IT_0094;
    const ccomplex_t IT_0096 = IT_0045*IT_0051*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*N_d4*e_em*m_mu
      *IT_0005*IT_0028*IT_0030;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0080*IT_0098;
    const ccomplex_t IT_0100 = IT_0053*IT_0059*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*N_d1*e_em*m_mu
      *IT_0005*IT_0028*IT_0030;
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0032*IT_0102;
    const ccomplex_t IT_0104 = IT_0027*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = m_mu*IT_0105;
    const ccomplex_t IT_0107 = IT_0041*IT_0090;
    const ccomplex_t IT_0108 = IT_0043*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = m_mu*IT_0109;
    const ccomplex_t IT_0111 = IT_0049*IT_0094;
    const ccomplex_t IT_0112 = IT_0051*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = m_mu*IT_0113;
    const ccomplex_t IT_0115 = IT_0057*IT_0098;
    const ccomplex_t IT_0116 = IT_0059*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = m_mu*IT_0117;
    const ccomplex_t IT_0119 = IT_0062*IT_0102;
    const ccomplex_t IT_0120 = IT_0025*IT_0027*IT_0119;
    const ccomplex_t IT_0121 = -IT_0021 + IT_0024 + -IT_0092 + -IT_0096 + 
      -IT_0100 + IT_0106 + IT_0110 + IT_0114 + IT_0118 + -IT_0120;
    const ccomplex_t IT_0122 = (-2)*conj(IT_0121);
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = 2*IT_0123;
    const ccomplex_t IT_0126 = IT_0124*IT_0125;
    const ccomplex_t IT_0127 = cpow(s_12 + IT_0010 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0126*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0013*IT_0017;
    const ccomplex_t IT_0132 = -IT_0105 + -IT_0109 + -IT_0113 + -IT_0117 + 
      -IT_0131;
    const ccomplex_t IT_0133 = s_13*s_14;
    const ccomplex_t IT_0134 = s_34*IT_0010;
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = IT_0010*IT_0026;
    const ccomplex_t IT_0137 = (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = IT_0133 + IT_0135 + IT_0137;
    const ccomplex_t IT_0139 = -s_14;
    const ccomplex_t IT_0140 = s_13 + IT_0139;
    const ccomplex_t IT_0141 = m_mu*IT_0140;
    const ccomplex_t IT_0142 = sin(alpha);
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*e_em*m_mu*IT_0005*IT_0028
      *IT_0030*IT_0142;
    const ccomplex_t IT_0144 = 0.5*IT_0143;
    const ccomplex_t IT_0145 = cpow(IT_0000, -2);
    const ccomplex_t IT_0146 = cos(alpha);
    const ccomplex_t IT_0147 = sin(beta);
    const ccomplex_t IT_0148 = (0 + _Complex_I*-1)*e_em*(IT_0005*IT_0026
      *IT_0028*IT_0030*IT_0142 + -m_W*IT_0002*IT_0145*(IT_0029*IT_0142 + IT_0146
      *IT_0147));
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = IT_0144*IT_0149;
    const ccomplex_t IT_0151 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0152 = IT_0150*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*e_em*m_mu*IT_0005*IT_0028
      *IT_0030*IT_0146;
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*-1)*e_em*(IT_0005*IT_0026
      *IT_0028*IT_0030*IT_0146 + -m_W*IT_0002*IT_0145*(IT_0029*IT_0146 + 
      -IT_0142*IT_0147));
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = cpow((-2)*s_12 + (-2)*IT_0010 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = -IT_0153 + -IT_0160;
    const ccomplex_t IT_0162 = (-2)*IT_0121;
    const ccomplex_t IT_0163 = (-2)*IT_0085;
    const ccomplex_t IT_0164 = 2*conj(IT_0161);
    const ccomplex_t IT_0165 = (0.5*IT_0087 + (-0.5)*conj(IT_0087) + IT_0130 +
       0.5*IT_0132 + (-0.5)*conj(IT_0132))*IT_0138 + (-0.25)*IT_0141*(conj
      (IT_0085) + conj(IT_0121) + (-2)*IT_0161 + 0.5*IT_0162 + 0.5*IT_0163 +
       IT_0164);
    const ccomplex_t IT_0166 = -IT_0026;
    const ccomplex_t IT_0167 = s_34 + IT_0166;
    const ccomplex_t IT_0168 = m_mu*s_13;
    const ccomplex_t IT_0169 = IT_0133 + IT_0135;
    const ccomplex_t IT_0170 = -IT_0132;
    const ccomplex_t IT_0171 = 0.5*IT_0088;
    const ccomplex_t IT_0172 = 2*s_34;
    const ccomplex_t IT_0173 = IT_0161*IT_0167;
    const ccomplex_t IT_0174 = (IT_0026*IT_0085 + IT_0087*IT_0088)*IT_0122 + (
      -8)*IT_0130*IT_0165 + IT_0164*((IT_0087 + IT_0132)*IT_0141 + (IT_0085 +
       IT_0121 + 2*IT_0161)*IT_0167) + 2*conj(IT_0087)*(IT_0141*IT_0161 +
       IT_0085*IT_0168 + 2*IT_0087*IT_0169 + IT_0136*IT_0170 + IT_0162*IT_0171) 
      + (-2)*conj(IT_0132)*(IT_0087*IT_0136 + -IT_0141*IT_0161 + -IT_0121
      *IT_0168 + (-2)*IT_0132*IT_0169 + -IT_0163*IT_0171) + 2*conj(IT_0085)*(0.5
      *IT_0026*IT_0162 + IT_0087*IT_0168 + IT_0088*IT_0170 + 0.5*IT_0085*IT_0172
       + IT_0173) + 2*conj(IT_0121)*(IT_0132*IT_0168 + 0.5*IT_0121*IT_0172 +
       IT_0173);
    return create_ccomplex_return(IT_0174);
}

