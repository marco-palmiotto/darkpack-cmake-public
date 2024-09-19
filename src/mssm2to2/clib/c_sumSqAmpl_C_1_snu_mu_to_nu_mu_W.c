#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_snu_mu_to_nu_mu_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_snu_mu_to_nu_mu_W(
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
    const creal_t m_C_1 = param->m_C_1;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B3)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W3)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0013 = IT_0006*IT_0012;
    const ccomplex_t IT_0014 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + 0.5*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0021 = cpow(s_23 + (-0.5)*IT_0020 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = conj(N_B4)*e_em;
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W4)*e_em;
    const ccomplex_t IT_0028 = IT_0006*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0026 + -IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0033 = IT_0006*IT_0032;
    const ccomplex_t IT_0034 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0035 = IT_0006*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + 0.5*IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = IT_0031*IT_0038;
    const ccomplex_t IT_0040 = cpow(s_23 + (-0.5)*IT_0020 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = conj(N_B1)*e_em;
    const ccomplex_t IT_0044 = IT_0001*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = conj(N_W1)*e_em;
    const ccomplex_t IT_0047 = IT_0006*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + -IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0052 = IT_0006*IT_0051;
    const ccomplex_t IT_0053 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0054 = IT_0006*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0052 + 0.5*IT_0055);
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0050*IT_0057;
    const ccomplex_t IT_0059 = cpow(s_23 + (-0.5)*IT_0020 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B2)*e_em;
    const ccomplex_t IT_0063 = IT_0001*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = conj(N_W2)*e_em;
    const ccomplex_t IT_0066 = IT_0006*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + -IT_0067);
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0071 = IT_0006*IT_0070;
    const ccomplex_t IT_0072 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0073 = IT_0006*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0071 + 0.5*IT_0074);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0069*IT_0076;
    const ccomplex_t IT_0078 = cpow(s_23 + (-0.5)*IT_0020 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0023 + (-0.5)*IT_0042 + (-0.5)
      *IT_0061 + (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = s_13*s_34;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0006;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0006;
    const ccomplex_t IT_0086 = pow(m_C_1, 2);
    const ccomplex_t IT_0087 = cpow(s_13 + (-0.5)*IT_0086 + 0.5*m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0084*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = pow(m_W, 2);
    const ccomplex_t IT_0093 = s_13*IT_0092;
    const ccomplex_t IT_0094 = s_14*s_34;
    const ccomplex_t IT_0095 = pow(m_W, -2);
    const ccomplex_t IT_0096 = pow(m_W, 4);
    const ccomplex_t IT_0097 = s_13*IT_0096;
    const ccomplex_t IT_0098 = IT_0095*IT_0097;
    const ccomplex_t IT_0099 = IT_0023 + IT_0042 + IT_0061 + IT_0080 + 0.5
      *IT_0090;
    const ccomplex_t IT_0100 = s_13*s_14;
    const ccomplex_t IT_0101 = s_34*IT_0086;
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0105 = IT_0006*IT_0104;
    const ccomplex_t IT_0106 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0107 = IT_0006*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + (-0.5)*IT_0108);
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0069*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0113 = IT_0078*IT_0111*IT_0112;
    const ccomplex_t IT_0114 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0115 = IT_0006*IT_0114;
    const ccomplex_t IT_0116 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0117 = IT_0006*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + (-0.5)*IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = IT_0011*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0123 = IT_0021*IT_0121*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0125 = IT_0006*IT_0124;
    const ccomplex_t IT_0126 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0127 = IT_0006*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + (-0.5)*IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0031*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0133 = IT_0040*IT_0131*IT_0132;
    const ccomplex_t IT_0134 = m_C_1*IT_0023;
    const ccomplex_t IT_0135 = m_C_1*IT_0042;
    const ccomplex_t IT_0136 = m_C_1*IT_0061;
    const ccomplex_t IT_0137 = m_C_1*IT_0080;
    const ccomplex_t IT_0138 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0139 = IT_0006*IT_0138;
    const ccomplex_t IT_0140 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0141 = IT_0006*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + (-0.5)*IT_0142);
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0050*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0147 = IT_0059*IT_0145*IT_0146;
    const ccomplex_t IT_0148 = 0.5*IT_0113 + 0.5*IT_0123 + 0.5*IT_0133 + (-0.5
      )*IT_0134 + (-0.5)*IT_0135 + (-0.5)*IT_0136 + (-0.5)*IT_0137 + 0.5*IT_0147;
    const ccomplex_t IT_0149 = s_34*IT_0092;
    const ccomplex_t IT_0150 = IT_0095*IT_0149;
    const ccomplex_t IT_0151 = (-0.25)*IT_0150;
    const ccomplex_t IT_0152 = s_34 + IT_0151;
    const ccomplex_t IT_0153 = m_C_1*IT_0152;
    const ccomplex_t IT_0154 = (-4)*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0154;
    const ccomplex_t IT_0156 = (-2)*IT_0094;
    const ccomplex_t IT_0157 = IT_0093 + IT_0156;
    const ccomplex_t IT_0158 = IT_0095*IT_0157;
    const ccomplex_t IT_0159 = pow(s_34, 2);
    const ccomplex_t IT_0160 = m_C_1*IT_0095*IT_0159;
    const ccomplex_t IT_0161 = 2*IT_0160;
    const ccomplex_t IT_0162 = s_14*s_34*IT_0095;
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = s_13 + IT_0163;
    const ccomplex_t IT_0165 = m_C_1*IT_0164;
    const ccomplex_t IT_0166 = (-2)*IT_0165;
    const ccomplex_t IT_0167 = s_13*IT_0095*IT_0159;
    const ccomplex_t IT_0168 = s_13*IT_0164;
    const ccomplex_t IT_0169 = (-2)*IT_0168;
    const ccomplex_t IT_0170 = 0.5*conj(IT_0148);
    const ccomplex_t IT_0171 = pow(s_14, 2);
    const ccomplex_t IT_0172 = IT_0095*IT_0171;
    const ccomplex_t IT_0173 = -IT_0086;
    const ccomplex_t IT_0174 = IT_0172 + IT_0173;
    const ccomplex_t IT_0175 = s_13*IT_0174;
    const ccomplex_t IT_0176 = (-2)*IT_0081*(IT_0082*conj(IT_0091) + 2*conj
      (IT_0081)*(IT_0093 + (-2)*IT_0094 + (-0.5)*IT_0098) + -conj(IT_0099)
      *IT_0103 + (-0.5)*conj(IT_0148)*IT_0155) + 4*IT_0148*(0.25*conj(IT_0081)
      *IT_0155 + conj(IT_0148)*(s_13 + (-0.5)*IT_0158) + 0.25*conj(IT_0091)
      *IT_0161 + 0.25*conj(IT_0099)*IT_0166) + (-2)*IT_0091*(conj(IT_0081)
      *IT_0082 + -conj(IT_0091)*IT_0167 + (-0.5)*conj(IT_0099)*IT_0169 + 
      -IT_0161*IT_0170) + 2*IT_0099*(conj(IT_0081)*IT_0103 + 0.5*conj(IT_0091)
      *IT_0169 + IT_0166*IT_0170 + conj(IT_0099)*IT_0175);
    return create_ccomplex_return(IT_0176);
}

