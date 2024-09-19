#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_snu_mu_to_mu_anti_nu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_snu_mu_to_mu_anti_nu_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_snu_mu = param->m_snu_mu;
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
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0000*IT_0002;
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0000*IT_0002*IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, 2);
    const ccomplex_t IT_0012 = pow(m_smu_R, 2);
    const ccomplex_t IT_0013 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + IT_0011 + -IT_0012 + -IT_0013 
      + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0000*IT_0002*IT_0006*IT_0007;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0000*IT_0002;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = cos(theta_W);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = N_B1*e_em;
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = N_W1*e_em;
    const ccomplex_t IT_0033 = IT_0002*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + -IT_0034);
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0028*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0039 = pow(m_mu, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0012 + IT_0039 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0037*IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0026;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = N_B2*e_em;
    const ccomplex_t IT_0045 = IT_0026*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = N_W2*e_em;
    const ccomplex_t IT_0048 = IT_0002*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + -IT_0049);
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0043*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0054 = cpow((-2)*s_13 + IT_0012 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0052*IT_0053*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0026;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = N_B3*e_em;
    const ccomplex_t IT_0059 = IT_0026*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = N_W3*e_em;
    const ccomplex_t IT_0062 = IT_0002*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0060 + -IT_0063);
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0057*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0068 = cpow((-2)*s_13 + IT_0012 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0066*IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0026;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = N_B4*e_em;
    const ccomplex_t IT_0073 = IT_0026*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = N_W4*e_em;
    const ccomplex_t IT_0076 = IT_0002*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + -IT_0077);
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0071*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0082 = cpow((-2)*s_13 + IT_0012 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0080*IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0006;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0036*IT_0085;
    const ccomplex_t IT_0087 = IT_0040*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = m_mu*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0006;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = IT_0051*IT_0091;
    const ccomplex_t IT_0093 = IT_0054*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = m_mu*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0000*IT_0002*IT_0006;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = IT_0065*IT_0097;
    const ccomplex_t IT_0099 = IT_0068*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = m_mu*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0000*IT_0002*IT_0006;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0079*IT_0103;
    const ccomplex_t IT_0105 = IT_0082*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = m_mu*IT_0106;
    const ccomplex_t IT_0108 = -IT_0016 + -IT_0024 + -IT_0041 + -IT_0055 + 
      -IT_0069 + -IT_0083 + IT_0089 + IT_0095 + IT_0101 + IT_0107;
    const ccomplex_t IT_0109 = -IT_0088 + -IT_0094 + -IT_0100 + -IT_0106;
    const ccomplex_t IT_0110 = m_mu*s_14;
    const ccomplex_t IT_0111 = s_13*s_14;
    const ccomplex_t IT_0112 = s_34*IT_0012;
    return create_ccomplex_return(2*IT_0108*(s_34*conj(IT_0108) + conj(IT_0109
      )*IT_0110) + 2*IT_0109*(conj(IT_0108)*IT_0110 + 2*conj(IT_0109)*(IT_0111 +
       (-0.5)*IT_0112)));
}

