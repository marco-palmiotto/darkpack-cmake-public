#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_stau_2_to_nu_tau_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_stau_2_to_nu_tau_anti_W(
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
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + 0.5*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = cos(beta);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = pow(m_W, -1);
    const ccomplex_t IT_0015 = conj(U_d2)*e_em*m_tau*IT_0014*U_stau_11;
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0010 + (-0.5)*IT_0017);
    const ccomplex_t IT_0019 = IT_0008*IT_0018;
    const ccomplex_t IT_0020 = pow(m_stau_2, 2);
    const ccomplex_t IT_0021 = cpow(s_23 + (-0.5)*IT_0020 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W4)*e_em*U_Wm1;
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = conj(N_d4)*U_d1*e_em;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + 0.5*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = conj(U_d1)*e_em*m_tau*IT_0014*U_stau_11;
    const ccomplex_t IT_0034 = IT_0013*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0032 + (-0.5)*IT_0035);
    const ccomplex_t IT_0037 = IT_0030*IT_0036;
    const ccomplex_t IT_0038 = cpow(s_23 + (-0.5)*IT_0020 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0023 + IT_0040;
    const ccomplex_t IT_0042 = m_N_4*IT_0041;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = pow(m_N_4, 2);
    const ccomplex_t IT_0045 = pow(m_tau, 2);
    const ccomplex_t IT_0046 = cpow(s_12 + 0.5*IT_0020 + 0.5*IT_0044 + (-0.5)
      *IT_0045 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = cos(theta_W);
    const ccomplex_t IT_0051 = cpow(IT_0050, -1);
    const ccomplex_t IT_0052 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = N_d4*e_em*m_tau*IT_0014*U_stau_01;
    const ccomplex_t IT_0056 = IT_0013*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + 0.5*IT_0057);
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = IT_0049*IT_0059;
    const ccomplex_t IT_0061 = IT_0046*IT_0047*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0063 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0064 = IT_0001*IT_0063;
    const ccomplex_t IT_0065 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0066 = IT_0001*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + (-0.5)*IT_0067);
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = IT_0018*IT_0069;
    const ccomplex_t IT_0071 = IT_0021*IT_0062*IT_0070;
    const ccomplex_t IT_0072 = N_W4*e_em*conj(V_Wp1);
    const ccomplex_t IT_0073 = IT_0001*IT_0072;
    const ccomplex_t IT_0074 = N_u4*conj(V_u1)*e_em;
    const ccomplex_t IT_0075 = IT_0001*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + (-0.5)*IT_0076);
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = IT_0036*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0081 = IT_0038*IT_0079*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0061 + 0.5*IT_0071 + 0.5*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0084 = IT_0051*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0087 = IT_0001*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = conj(N_d4)*e_em*m_tau*IT_0014*U_stau_11;
    const ccomplex_t IT_0090 = IT_0013*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + IT_0088 + 
      -IT_0091);
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = IT_0049*IT_0093;
    const ccomplex_t IT_0095 = IT_0046*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = pow(m_W, -2);
    const ccomplex_t IT_0099 = pow(m_W, 2);
    const ccomplex_t IT_0100 = s_34*IT_0099;
    const ccomplex_t IT_0101 = IT_0098*IT_0100;
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = s_34 + IT_0102;
    const ccomplex_t IT_0104 = m_N_4*IT_0103;
    const ccomplex_t IT_0105 = 2*IT_0104;
    const ccomplex_t IT_0106 = pow(s_34, 2);
    const ccomplex_t IT_0107 = m_N_4*IT_0098*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *U_stau_01;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = conj(N_B4)*e_em;
    const ccomplex_t IT_0111 = IT_0051*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W4)*e_em;
    const ccomplex_t IT_0114 = IT_0001*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + -IT_0115);
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = cpow(s_13 + (-0.5)*IT_0044 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0119 = IT_0117*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = IT_0109*IT_0120;
    const ccomplex_t IT_0122 = -IT_0096 + (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = 2*IT_0122;
    const ccomplex_t IT_0124 = (-0.5)*IT_0107*IT_0123;
    const ccomplex_t IT_0125 = s_13*IT_0099;
    const ccomplex_t IT_0126 = s_14*s_34;
    const ccomplex_t IT_0127 = (-2)*IT_0126;
    const ccomplex_t IT_0128 = IT_0125 + IT_0127;
    const ccomplex_t IT_0129 = IT_0098*IT_0128;
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = s_13 + IT_0130;
    const ccomplex_t IT_0132 = IT_0043 + IT_0082;
    const ccomplex_t IT_0133 = (-2)*IT_0131*IT_0132;
    const ccomplex_t IT_0134 = 0.5*IT_0023 + 0.5*IT_0040;
    const ccomplex_t IT_0135 = (-0.25)*IT_0101;
    const ccomplex_t IT_0136 = s_34 + IT_0135;
    const ccomplex_t IT_0137 = m_N_4*IT_0136;
    const ccomplex_t IT_0138 = (-4)*IT_0137;
    const ccomplex_t IT_0139 = IT_0134*IT_0138;
    const ccomplex_t IT_0140 = -IT_0097*IT_0105 + IT_0124 + IT_0133 + IT_0139;
    const ccomplex_t IT_0141 = s_14*s_34*IT_0098;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = s_13 + IT_0142;
    const ccomplex_t IT_0144 = s_13*IT_0143;
    const ccomplex_t IT_0145 = -IT_0023 + -IT_0040;
    const ccomplex_t IT_0146 = 2*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0121;
    const ccomplex_t IT_0148 = (-2)*IT_0147;
    const ccomplex_t IT_0149 = IT_0146 + IT_0148;
    const ccomplex_t IT_0150 = 2*conj(IT_0145);
    const ccomplex_t IT_0151 = (-2)*conj(IT_0147);
    const ccomplex_t IT_0152 = IT_0150 + IT_0151;
    const ccomplex_t IT_0153 = m_N_4*IT_0143;
    const ccomplex_t IT_0154 = s_13*s_14;
    const ccomplex_t IT_0155 = s_34*IT_0044;
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = IT_0097 + IT_0134;
    const ccomplex_t IT_0159 = pow(s_14, 2);
    const ccomplex_t IT_0160 = IT_0098*IT_0159;
    const ccomplex_t IT_0161 = -IT_0044;
    const ccomplex_t IT_0162 = IT_0160 + IT_0161;
    const ccomplex_t IT_0163 = s_13*IT_0162;
    const ccomplex_t IT_0164 = s_13*IT_0098*IT_0106;
    const ccomplex_t IT_0165 = pow(m_W, 4);
    const ccomplex_t IT_0166 = s_13*IT_0165;
    const ccomplex_t IT_0167 = IT_0098*IT_0166;
    const ccomplex_t IT_0168 = IT_0125 + IT_0127 + (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = (-0.25)*IT_0167;
    const ccomplex_t IT_0170 = IT_0126 + IT_0169;
    const ccomplex_t IT_0171 = s_13*s_34;
    const ccomplex_t IT_0172 = 2*conj(IT_0122);
    const ccomplex_t IT_0173 = (-2)*(conj(IT_0043) + conj(IT_0082))*IT_0140 +
       IT_0144*(conj(IT_0122)*IT_0149 + IT_0122*IT_0152) + ((conj(IT_0043) +
       conj(IT_0082))*IT_0149 + IT_0132*IT_0152)*IT_0153 + IT_0157*((conj
      (IT_0097) + conj(IT_0134))*IT_0149 + IT_0152*IT_0158) + (-2)*((-0.5)*conj
      (IT_0145)*IT_0146 + IT_0145*conj(IT_0147) + IT_0147*(conj(IT_0145) + -conj
      (IT_0147)))*IT_0163 + conj(IT_0122)*IT_0123*IT_0164 + (-4)*conj(IT_0134)*
      (0.5*IT_0132*IT_0138 + IT_0134*IT_0168 + (-2)*IT_0097*IT_0170 + (-0.25)
      *IT_0123*IT_0171) + (-4)*conj(IT_0097)*((-0.5)*IT_0105*IT_0132 + IT_0097
      *IT_0168 + (-2)*IT_0134*IT_0170 + (-0.25)*IT_0123*IT_0171) + (IT_0107
      *IT_0132 + IT_0158*IT_0171)*IT_0172;
    return create_ccomplex_return(IT_0173);
}

