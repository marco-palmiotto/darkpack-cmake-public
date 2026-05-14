#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_stau_1_to_nu_tau_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_stau_1_to_nu_tau_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
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
    const ccomplex_t IT_0011 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0012 + 0.5*IT_0015);
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0019 = IT_0010*IT_0018;
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = IT_0010*IT_0021;
    const ccomplex_t IT_0023 = pow(m_W, -1);
    const ccomplex_t IT_0024 = conj(U_d2)*e_em*m_tau*IT_0023*U_stau_10;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0019 + (-0.5)*IT_0026);
    const ccomplex_t IT_0028 = IT_0017*IT_0027;
    const ccomplex_t IT_0029 = pow(m_stau_1, 2);
    const ccomplex_t IT_0030 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = m_N_1*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0036 = IT_0010*IT_0035;
    const ccomplex_t IT_0037 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0038 = IT_0010*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + (-0.5)*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0043 = IT_0010*IT_0042;
    const ccomplex_t IT_0044 = conj(U_d1)*e_em*m_tau*IT_0023*U_stau_10;
    const ccomplex_t IT_0045 = IT_0022*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + (-0.5)*IT_0046);
    const ccomplex_t IT_0048 = IT_0041*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0050 = cpow(s_23 + (-0.5)*IT_0029 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0048*IT_0049*IT_0050;
    const ccomplex_t IT_0052 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0053 = IT_0010*IT_0052;
    const ccomplex_t IT_0054 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0055 = IT_0010*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + (-0.5)*IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = IT_0027*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0061 = IT_0030*IT_0059*IT_0060;
    const ccomplex_t IT_0062 = pow(m_N_1, 2);
    const ccomplex_t IT_0063 = pow(m_tau, 2);
    const ccomplex_t IT_0064 = cpow(s_12 + 0.5*IT_0029 + 0.5*IT_0062 + (-0.5)
      *IT_0063 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0010;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = cos(theta_W);
    const ccomplex_t IT_0069 = cpow(IT_0068, -1);
    const ccomplex_t IT_0070 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = N_d1*e_em*m_tau*IT_0023*U_stau_00;
    const ccomplex_t IT_0074 = IT_0022*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + 0.5*IT_0075);
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0067*IT_0077;
    const ccomplex_t IT_0079 = IT_0064*IT_0065*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0081 = IT_0010*IT_0080;
    const ccomplex_t IT_0082 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0083 = IT_0010*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + 0.5*IT_0084);
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = IT_0047*IT_0086;
    const ccomplex_t IT_0088 = IT_0050*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = m_N_1*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0051 + 0.5*IT_0061 + 0.5*IT_0079 + (-0.5
      )*IT_0090;
    const ccomplex_t IT_0092 = IT_0034 + IT_0091;
    const ccomplex_t IT_0093 = conj(IT_0034) + conj(IT_0091);
    const ccomplex_t IT_0094 = pow(m_W, 4);
    const ccomplex_t IT_0095 = s_13*IT_0094;
    const ccomplex_t IT_0096 = IT_0000*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = IT_0002 + IT_0004 + IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0032;
    const ccomplex_t IT_0100 = (-0.5)*IT_0089;
    const ccomplex_t IT_0101 = IT_0099 + -IT_0100;
    const ccomplex_t IT_0102 = -conj(IT_0100);
    const ccomplex_t IT_0103 = conj(IT_0099) + IT_0102;
    const ccomplex_t IT_0104 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = s_13 + IT_0105;
    const ccomplex_t IT_0107 = s_13*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0010
      *U_stau_00;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = conj(N_B1)*e_em;
    const ccomplex_t IT_0111 = IT_0069*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W1)*e_em;
    const ccomplex_t IT_0114 = IT_0010*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + -IT_0115);
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = cpow(s_13 + (-0.5)*IT_0062 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0119 = IT_0117*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = IT_0109*IT_0120;
    const ccomplex_t IT_0122 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0123 = IT_0069*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0126 = IT_0010*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = conj(N_d1)*e_em*m_tau*IT_0023*U_stau_10;
    const ccomplex_t IT_0129 = IT_0022*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0124 + IT_0127 + 
      -IT_0130);
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = IT_0067*IT_0132;
    const ccomplex_t IT_0134 = IT_0064*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0121 + -IT_0135;
    const ccomplex_t IT_0137 = -IT_0032;
    const ccomplex_t IT_0138 = 2*IT_0137;
    const ccomplex_t IT_0139 = IT_0089 + 0.5*IT_0121;
    const ccomplex_t IT_0140 = (-2)*IT_0139;
    const ccomplex_t IT_0141 = IT_0138 + IT_0140;
    const ccomplex_t IT_0142 = 2*conj(IT_0137);
    const ccomplex_t IT_0143 = (-2)*conj(IT_0139);
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = m_N_1*IT_0106;
    const ccomplex_t IT_0146 = (-0.5)*IT_0135;
    const ccomplex_t IT_0147 = (-0.25)*IT_0096;
    const ccomplex_t IT_0148 = IT_0003 + IT_0147;
    const ccomplex_t IT_0149 = s_34*IT_0001;
    const ccomplex_t IT_0150 = IT_0000*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = s_34 + IT_0151;
    const ccomplex_t IT_0153 = m_N_1*IT_0152;
    const ccomplex_t IT_0154 = 2*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0034;
    const ccomplex_t IT_0156 = 2*IT_0091;
    const ccomplex_t IT_0157 = IT_0155 + IT_0156;
    const ccomplex_t IT_0158 = (-0.25)*IT_0150;
    const ccomplex_t IT_0159 = s_34 + IT_0158;
    const ccomplex_t IT_0160 = m_N_1*IT_0159;
    const ccomplex_t IT_0161 = (-4)*IT_0160;
    const ccomplex_t IT_0162 = 2*conj(IT_0034);
    const ccomplex_t IT_0163 = 2*conj(IT_0091);
    const ccomplex_t IT_0164 = IT_0162 + IT_0163;
    const ccomplex_t IT_0165 = pow(s_34, 2);
    const ccomplex_t IT_0166 = s_13*IT_0000*IT_0165;
    const ccomplex_t IT_0167 = 2*IT_0136;
    const ccomplex_t IT_0168 = m_N_1*IT_0000*IT_0165;
    const ccomplex_t IT_0169 = 2*conj(IT_0136);
    const ccomplex_t IT_0170 = s_13*s_34;
    const ccomplex_t IT_0171 = conj(IT_0099) + conj(IT_0146);
    const ccomplex_t IT_0172 = IT_0099 + IT_0146;
    const ccomplex_t IT_0173 = pow(s_14, 2);
    const ccomplex_t IT_0174 = IT_0000*IT_0173;
    const ccomplex_t IT_0175 = -IT_0062;
    const ccomplex_t IT_0176 = IT_0174 + IT_0175;
    const ccomplex_t IT_0177 = s_13*IT_0176;
    const ccomplex_t IT_0178 = (-2)*conj(IT_0137);
    const ccomplex_t IT_0179 = 2*IT_0139;
    const ccomplex_t IT_0180 = (-0.5)*IT_0179;
    const ccomplex_t IT_0181 = s_13*s_14;
    const ccomplex_t IT_0182 = s_34*IT_0062;
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = IT_0181 + IT_0183;
    const ccomplex_t IT_0185 = 4*IT_0008*IT_0092*IT_0093 + (-4)*IT_0098
      *IT_0101*IT_0103 + IT_0107*(conj(IT_0136)*IT_0141 + IT_0136*IT_0144) + 
      (IT_0093*IT_0141 + IT_0092*IT_0144)*IT_0145 + 8*conj(IT_0146)*(IT_0101
      *IT_0148 + 0.125*IT_0154*IT_0157) + (-2)*IT_0161*(IT_0093*IT_0099 +
       IT_0092*conj(IT_0099) + (-0.5)*conj(IT_0100)*IT_0157 + (-0.5)*IT_0100
      *IT_0164) + (-4)*IT_0146*(IT_0098*conj(IT_0146) + (-2)*IT_0103*IT_0148 + (
      -0.25)*IT_0154*IT_0164) + conj(IT_0136)*IT_0166*IT_0167 + IT_0168*(IT_0093
      *IT_0167 + IT_0092*IT_0169) + (-2)*IT_0170*(conj(IT_0100)*IT_0136 +
       IT_0100*conj(IT_0136) + (-0.5)*IT_0167*IT_0171 + (-0.5)*IT_0169*IT_0172) 
      + IT_0177*(conj(IT_0137)*IT_0138 + IT_0139*IT_0178 + (-2)*conj(IT_0139)*
      (IT_0137 + IT_0180)) + (IT_0141*IT_0171 + IT_0144*IT_0172 + 2*IT_0100*
      (conj(IT_0139) + 0.5*IT_0178) + 2*IT_0102*(IT_0137 + IT_0180))*IT_0184;
    return create_ccomplex_return(IT_0185);
}

