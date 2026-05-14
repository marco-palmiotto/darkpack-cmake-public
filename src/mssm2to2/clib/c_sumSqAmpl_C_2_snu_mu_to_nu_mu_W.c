#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_snu_mu_to_nu_mu_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_snu_mu_to_nu_mu_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ml = param->Gamma_ml;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0005 = pow(m_C_2, 2);
    const ccomplex_t IT_0006 = cpow(s_13 + (-0.5)*IT_0005 + 0.5*m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -2);
    const ccomplex_t IT_0012 = pow(s_34, 2);
    const ccomplex_t IT_0013 = s_13*IT_0011*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = conj(N_B2)*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = conj(N_W2)*e_em;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0018 + -IT_0021);
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + (-0.5)*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0023*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0033 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0034 = cpow(s_23 + (-0.5)*IT_0033 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0031*IT_0032*IT_0034;
    const ccomplex_t IT_0036 = conj(N_B3)*e_em;
    const ccomplex_t IT_0037 = IT_0015*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W3)*e_em;
    const ccomplex_t IT_0040 = IT_0001*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0045 = IT_0001*IT_0044;
    const ccomplex_t IT_0046 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0047 = IT_0001*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + (-0.5)*IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = IT_0043*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0053 = cpow(s_23 + (-0.5)*IT_0033 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0051*IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0056 = cpow(s_23 + (-0.5)*IT_0033 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0057 = conj(N_B1)*e_em;
    const ccomplex_t IT_0058 = IT_0015*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W1)*e_em;
    const ccomplex_t IT_0061 = IT_0001*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + -IT_0062);
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0066 = IT_0001*IT_0065;
    const ccomplex_t IT_0067 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0068 = IT_0001*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + (-0.5)*IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0064*IT_0071;
    const ccomplex_t IT_0073 = IT_0055*IT_0056*IT_0072;
    const ccomplex_t IT_0074 = conj(N_B4)*e_em;
    const ccomplex_t IT_0075 = IT_0015*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W4)*e_em;
    const ccomplex_t IT_0078 = IT_0001*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0076 + -IT_0079);
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0083 = IT_0001*IT_0082;
    const ccomplex_t IT_0084 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0085 = IT_0001*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + (-0.5)*IT_0086);
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0081*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0091 = cpow(s_23 + (-0.5)*IT_0033 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0089*IT_0090*IT_0091;
    const ccomplex_t IT_0093 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0094 = IT_0001*IT_0093;
    const ccomplex_t IT_0095 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0096 = IT_0001*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + 0.5*IT_0097);
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0064*IT_0099;
    const ccomplex_t IT_0101 = IT_0056*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = m_C_2*IT_0102;
    const ccomplex_t IT_0104 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0105 = IT_0001*IT_0104;
    const ccomplex_t IT_0106 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0107 = IT_0001*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + 0.5*IT_0108);
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0023*IT_0110;
    const ccomplex_t IT_0112 = IT_0034*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = m_C_2*IT_0113;
    const ccomplex_t IT_0115 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0116 = IT_0001*IT_0115;
    const ccomplex_t IT_0117 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0118 = IT_0001*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*(IT_0116 + 0.5*IT_0119);
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = IT_0043*IT_0121;
    const ccomplex_t IT_0123 = IT_0053*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = m_C_2*IT_0124;
    const ccomplex_t IT_0126 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0127 = IT_0001*IT_0126;
    const ccomplex_t IT_0128 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0129 = IT_0001*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 0.5*IT_0130);
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0081*IT_0132;
    const ccomplex_t IT_0134 = IT_0091*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = m_C_2*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0035 + 0.5*IT_0054 + 0.5*IT_0073 + 0.5
      *IT_0092 + (-0.5)*IT_0103 + (-0.5)*IT_0114 + (-0.5)*IT_0125 + (-0.5)
      *IT_0136;
    const ccomplex_t IT_0138 = m_C_2*IT_0011*IT_0012;
    const ccomplex_t IT_0139 = (-0.5)*IT_0102 + (-0.5)*IT_0113 + (-0.5)
      *IT_0124 + (-0.5)*IT_0135;
    const ccomplex_t IT_0140 = s_13*s_34;
    const ccomplex_t IT_0141 = s_14*s_34*IT_0011;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = s_13 + IT_0142;
    const ccomplex_t IT_0144 = s_13*IT_0143;
    const ccomplex_t IT_0145 = 0.5*IT_0009 + IT_0102 + IT_0113 + IT_0124 +
       IT_0135;
    const ccomplex_t IT_0146 = (-2)*conj(IT_0145);
    const ccomplex_t IT_0147 = pow(m_W, 2);
    const ccomplex_t IT_0148 = s_13*IT_0147;
    const ccomplex_t IT_0149 = s_14*s_34;
    const ccomplex_t IT_0150 = (-2)*IT_0149;
    const ccomplex_t IT_0151 = IT_0148 + IT_0150;
    const ccomplex_t IT_0152 = IT_0011*IT_0151;
    const ccomplex_t IT_0153 = s_34*IT_0147;
    const ccomplex_t IT_0154 = IT_0011*IT_0153;
    const ccomplex_t IT_0155 = (-0.25)*IT_0154;
    const ccomplex_t IT_0156 = s_34 + IT_0155;
    const ccomplex_t IT_0157 = m_C_2*IT_0156;
    const ccomplex_t IT_0158 = (-4)*IT_0157;
    const ccomplex_t IT_0159 = m_C_2*IT_0143;
    const ccomplex_t IT_0160 = pow(s_14, 2);
    const ccomplex_t IT_0161 = IT_0011*IT_0160;
    const ccomplex_t IT_0162 = -IT_0005;
    const ccomplex_t IT_0163 = IT_0161 + IT_0162;
    const ccomplex_t IT_0164 = s_13*IT_0163;
    const ccomplex_t IT_0165 = s_13*s_14;
    const ccomplex_t IT_0166 = s_34*IT_0005;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = -conj(IT_0137);
    const ccomplex_t IT_0170 = pow(m_W, 4);
    const ccomplex_t IT_0171 = s_13*IT_0170;
    const ccomplex_t IT_0172 = IT_0011*IT_0171;
    const ccomplex_t IT_0173 = 2*IT_0010*(conj(IT_0010)*IT_0013 + conj(IT_0137
      )*IT_0138 + -conj(IT_0139)*IT_0140 + 0.5*IT_0144*IT_0146) + 2*IT_0137*
      (conj(IT_0010)*IT_0138 + 2*conj(IT_0137)*(s_13 + (-0.5)*IT_0152) + conj
      (IT_0139)*IT_0158 + 0.5*IT_0146*IT_0159) + (-2)*IT_0145*(conj(IT_0010)
      *IT_0144 + -conj(IT_0145)*IT_0164 + -conj(IT_0139)*IT_0168 + -IT_0159
      *IT_0169) + (-2)*IT_0139*(conj(IT_0010)*IT_0140 + -conj(IT_0145)*IT_0168 +
       IT_0158*IT_0169 + 2*conj(IT_0139)*(IT_0148 + (-2)*IT_0149 + (-0.5)
      *IT_0172));
    return create_ccomplex_return(IT_0173);
}

