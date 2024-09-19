#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_snu_tau_to_nu_tau_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_snu_tau_to_nu_tau_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
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
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
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
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *conj(U_stau_00);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = conj(U_d2)*e_em*m_tau*IT_0009*U_stau_10;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0005 + (-0.5)*IT_0012);
    const ccomplex_t IT_0014 = pow(m_C_2, 2);
    const ccomplex_t IT_0015 = cpow(s_13 + (-0.5)*IT_0014 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0003*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *conj(U_stau_01);
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = conj(U_d2)*e_em*m_tau*IT_0009*U_stau_11;
    const ccomplex_t IT_0024 = IT_0008*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + (-0.5)*IT_0025);
    const ccomplex_t IT_0027 = cpow(s_13 + (-0.5)*IT_0014 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0020*IT_0029;
    const ccomplex_t IT_0031 = (-0.5)*IT_0018 + (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = pow(m_W, -2);
    const ccomplex_t IT_0033 = pow(s_34, 2);
    const ccomplex_t IT_0034 = s_13*IT_0032*IT_0033;
    const ccomplex_t IT_0035 = cos(theta_W);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = conj(N_B2)*e_em;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W2)*e_em;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + -IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0046 = IT_0001*IT_0045;
    const ccomplex_t IT_0047 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0048 = IT_0001*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + 0.5*IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0044*IT_0051;
    const ccomplex_t IT_0053 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0054 = cpow(s_23 + (-0.5)*IT_0053 + 0.5*m_N_2*(m_N_2 +
       (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = m_C_2*IT_0056;
    const ccomplex_t IT_0058 = conj(N_B4)*e_em;
    const ccomplex_t IT_0059 = IT_0036*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = conj(N_W4)*e_em;
    const ccomplex_t IT_0062 = IT_0001*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0060 + -IT_0063);
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0067 = IT_0001*IT_0066;
    const ccomplex_t IT_0068 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0069 = IT_0001*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0067 + 0.5*IT_0070);
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0065*IT_0072;
    const ccomplex_t IT_0074 = cpow(s_23 + (-0.5)*IT_0053 + 0.5*m_N_4*(m_N_4 +
       (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = m_C_2*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B3)*e_em;
    const ccomplex_t IT_0079 = IT_0036*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = conj(N_W3)*e_em;
    const ccomplex_t IT_0082 = IT_0001*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0087 = IT_0001*IT_0086;
    const ccomplex_t IT_0088 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0089 = IT_0001*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + (-0.5)*IT_0090);
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = IT_0085*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0095 = cpow(s_23 + (-0.5)*IT_0053 + 0.5*m_N_3*(m_N_3 +
       (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0093*IT_0094*IT_0095;
    const ccomplex_t IT_0097 = cpow(s_23 + (-0.5)*IT_0053 + 0.5*m_N_1*(m_N_1 +
       (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0098 = conj(N_B1)*e_em;
    const ccomplex_t IT_0099 = IT_0036*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W1)*e_em;
    const ccomplex_t IT_0102 = IT_0001*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + -IT_0103);
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0107 = IT_0001*IT_0106;
    const ccomplex_t IT_0108 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0109 = IT_0001*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + 0.5*IT_0110);
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0105*IT_0112;
    const ccomplex_t IT_0114 = IT_0097*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = m_C_2*IT_0115;
    const ccomplex_t IT_0117 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0118 = IT_0001*IT_0117;
    const ccomplex_t IT_0119 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0120 = IT_0001*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0118 + 0.5*IT_0121);
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = IT_0085*IT_0123;
    const ccomplex_t IT_0125 = IT_0095*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = m_C_2*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0129 = IT_0001*IT_0128;
    const ccomplex_t IT_0130 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0131 = IT_0001*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + (-0.5)*IT_0132);
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = IT_0105*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0137 = IT_0097*IT_0135*IT_0136;
    const ccomplex_t IT_0138 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0139 = IT_0001*IT_0138;
    const ccomplex_t IT_0140 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0141 = IT_0001*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + (-0.5)*IT_0142);
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0065*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0147 = IT_0074*IT_0145*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0149 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0150 = IT_0001*IT_0149;
    const ccomplex_t IT_0151 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0152 = IT_0001*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0150 + (-0.5)*IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = IT_0044*IT_0155;
    const ccomplex_t IT_0157 = IT_0054*IT_0148*IT_0156;
    const ccomplex_t IT_0158 = (-0.5)*IT_0057 + (-0.5)*IT_0077 + 0.5*IT_0096 +
       (-0.5)*IT_0116 + (-0.5)*IT_0127 + 0.5*IT_0137 + 0.5*IT_0147 + 0.5*IT_0157;
    const ccomplex_t IT_0159 = m_C_2*IT_0032*IT_0033;
    const ccomplex_t IT_0160 = (-0.5)*IT_0056 + (-0.5)*IT_0076 + (-0.5)
      *IT_0115 + (-0.5)*IT_0126;
    const ccomplex_t IT_0161 = s_13*s_34;
    const ccomplex_t IT_0162 = s_14*s_34*IT_0032;
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = s_13 + IT_0163;
    const ccomplex_t IT_0165 = s_13*IT_0164;
    const ccomplex_t IT_0166 = 0.5*IT_0018 + 0.5*IT_0030 + IT_0056 + IT_0076 +
       IT_0115 + IT_0126;
    const ccomplex_t IT_0167 = (-2)*conj(IT_0166);
    const ccomplex_t IT_0168 = pow(m_W, 2);
    const ccomplex_t IT_0169 = s_13*IT_0168;
    const ccomplex_t IT_0170 = s_14*s_34;
    const ccomplex_t IT_0171 = (-2)*IT_0170;
    const ccomplex_t IT_0172 = IT_0169 + IT_0171;
    const ccomplex_t IT_0173 = IT_0032*IT_0172;
    const ccomplex_t IT_0174 = s_34*IT_0168;
    const ccomplex_t IT_0175 = IT_0032*IT_0174;
    const ccomplex_t IT_0176 = (-0.25)*IT_0175;
    const ccomplex_t IT_0177 = s_34 + IT_0176;
    const ccomplex_t IT_0178 = m_C_2*IT_0177;
    const ccomplex_t IT_0179 = (-4)*IT_0178;
    const ccomplex_t IT_0180 = m_C_2*IT_0164;
    const ccomplex_t IT_0181 = pow(s_14, 2);
    const ccomplex_t IT_0182 = IT_0032*IT_0181;
    const ccomplex_t IT_0183 = -IT_0014;
    const ccomplex_t IT_0184 = IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = s_13*IT_0184;
    const ccomplex_t IT_0186 = s_13*s_14;
    const ccomplex_t IT_0187 = s_34*IT_0014;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = IT_0186 + IT_0188;
    const ccomplex_t IT_0190 = -conj(IT_0158);
    const ccomplex_t IT_0191 = pow(m_W, 4);
    const ccomplex_t IT_0192 = s_13*IT_0191;
    const ccomplex_t IT_0193 = IT_0032*IT_0192;
    const ccomplex_t IT_0194 = 2*IT_0031*(conj(IT_0031)*IT_0034 + conj(IT_0158
      )*IT_0159 + -conj(IT_0160)*IT_0161 + 0.5*IT_0165*IT_0167) + 2*IT_0158*
      (conj(IT_0031)*IT_0159 + 2*conj(IT_0158)*(s_13 + (-0.5)*IT_0173) + conj
      (IT_0160)*IT_0179 + 0.5*IT_0167*IT_0180) + (-2)*IT_0166*(conj(IT_0031)
      *IT_0165 + -conj(IT_0166)*IT_0185 + -conj(IT_0160)*IT_0189 + -IT_0180
      *IT_0190) + (-2)*IT_0160*(conj(IT_0031)*IT_0161 + -conj(IT_0166)*IT_0189 +
       IT_0179*IT_0190 + 2*conj(IT_0160)*(IT_0169 + (-2)*IT_0170 + (-0.5)
      *IT_0193));
    return create_ccomplex_return(IT_0194);
}

