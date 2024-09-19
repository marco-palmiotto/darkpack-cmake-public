#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_ss_L_to_c_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_ss_L_to_c_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0001 = pow(m_c, 2);
    const ccomplex_t IT_0002 = pow(m_ss_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = IT_0005*IT_0008;
    const ccomplex_t IT_0012 = N_d1*e_em*U_Wm1;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = (-2)*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0005*IT_0016;
    const ccomplex_t IT_0018 = N_B1*U_d1*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = N_W1*U_d1*e_em;
    const ccomplex_t IT_0022 = IT_0011*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = IT_0020 + IT_0023;
    const ccomplex_t IT_0025 = IT_0014 + IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0010*IT_0027;
    const ccomplex_t IT_0029 = IT_0000*IT_0003*IT_0028;
    const ccomplex_t IT_0030 = pow(m_s, 2);
    const ccomplex_t IT_0031 = pow(m_N_1, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + -IT_0002 + IT_0030 + -IT_0031 
      + -reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0034 = cos(beta);
    const ccomplex_t IT_0035 = cpow(IT_0034, -1);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0004*IT_0008*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0004*IT_0006*IT_0008*IT_0034;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = IT_0032*IT_0033*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0008;
    const ccomplex_t IT_0043 = IT_0027*IT_0042;
    const ccomplex_t IT_0044 = IT_0003*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = m_c*IT_0045;
    const ccomplex_t IT_0047 = conj(N_B1)*e_em;
    const ccomplex_t IT_0048 = IT_0016*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = conj(N_W1)*e_em;
    const ccomplex_t IT_0051 = IT_0008*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-3)*IT_0052);
    const ccomplex_t IT_0054 = 0.166666666666667*IT_0053;
    const ccomplex_t IT_0055 = IT_0039*IT_0054;
    const ccomplex_t IT_0056 = IT_0032*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = m_N_1*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_13 + IT_0001 + IT_0031 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0060 = (0 + _Complex_I*(-2.82842712474619))*V_cs*e_em
      *IT_0008*((m_W*IT_0005 + (-0.5)*IT_0001*IT_0004*IT_0006)*IT_0034 + -1./2
      *IT_0004*IT_0005*IT_0030*IT_0035);
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*m_c*N_u1*e_em
      *IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = IT_0059*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_13 + IT_0001 + IT_0031 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *mu_h*IT_0004*IT_0008;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0016;
    const ccomplex_t IT_0071 = 0.666666666666667*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = IT_0067*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0076 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = N_d1*e_em*U_Wm2;
    const ccomplex_t IT_0080 = IT_0011*IT_0079;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = N_B1*U_d2*e_em;
    const ccomplex_t IT_0083 = IT_0017*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = N_W1*U_d2*e_em;
    const ccomplex_t IT_0086 = IT_0011*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = IT_0081 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = IT_0078*IT_0091;
    const ccomplex_t IT_0093 = IT_0075*IT_0076*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0008;
    const ccomplex_t IT_0095 = IT_0091*IT_0094;
    const ccomplex_t IT_0096 = IT_0076*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = m_c*IT_0097;
    const ccomplex_t IT_0099 = IT_0029 + -IT_0041 + -IT_0046 + IT_0058 +
       IT_0066 + IT_0074 + IT_0093 + -IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0004*IT_0005*IT_0008*IT_0035;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0054*IT_0101;
    const ccomplex_t IT_0103 = IT_0032*IT_0033*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0049 + 3*IT_0052);
    const ccomplex_t IT_0105 = 0.166666666666667*IT_0104;
    const ccomplex_t IT_0106 = IT_0061*IT_0105;
    const ccomplex_t IT_0107 = IT_0059*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0069*IT_0110;
    const ccomplex_t IT_0112 = IT_0067*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0008*IT_0034;
    const ccomplex_t IT_0115 = conj(N_u1)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = 2*IT_0116;
    const ccomplex_t IT_0118 = IT_0016*IT_0034;
    const ccomplex_t IT_0119 = conj(N_B1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = conj(N_W1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0123 = IT_0114*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = IT_0121 + IT_0124;
    const ccomplex_t IT_0126 = IT_0117 + IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0042*IT_0128;
    const ccomplex_t IT_0130 = IT_0000*IT_0003*IT_0129;
    const ccomplex_t IT_0131 = conj(N_u1)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0132 = IT_0114*IT_0131;
    const ccomplex_t IT_0133 = 2*IT_0132;
    const ccomplex_t IT_0134 = conj(N_B1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0135 = IT_0118*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = conj(N_W1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0138 = IT_0114*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = IT_0136 + IT_0139;
    const ccomplex_t IT_0141 = IT_0133 + IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = IT_0094*IT_0143;
    const ccomplex_t IT_0145 = IT_0075*IT_0076*IT_0144;
    const ccomplex_t IT_0146 = IT_0010*IT_0128;
    const ccomplex_t IT_0147 = IT_0003*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = m_c*IT_0148;
    const ccomplex_t IT_0150 = IT_0078*IT_0143;
    const ccomplex_t IT_0151 = IT_0076*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = m_c*IT_0152;
    const ccomplex_t IT_0154 = IT_0037*IT_0101;
    const ccomplex_t IT_0155 = IT_0032*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = m_N_1*IT_0156;
    const ccomplex_t IT_0158 = IT_0103 + -IT_0108 + -IT_0113 + IT_0130 +
       IT_0145 + -IT_0149 + -IT_0153 + -IT_0157;
    const ccomplex_t IT_0159 = 6*s_13;
    const ccomplex_t IT_0160 = m_c*m_N_1;
    const ccomplex_t IT_0161 = 6*IT_0099;
    const ccomplex_t IT_0162 = IT_0045 + IT_0097 + -IT_0156;
    const ccomplex_t IT_0163 = m_c*s_12;
    const ccomplex_t IT_0164 = IT_0057 + IT_0148 + IT_0152;
    const ccomplex_t IT_0165 = s_23*m_N_1;
    const ccomplex_t IT_0166 = 6*conj(IT_0099);
    const ccomplex_t IT_0167 = 6*conj(IT_0158);
    const ccomplex_t IT_0168 = m_c*m_N_1*IT_0002;
    const ccomplex_t IT_0169 = s_13*IT_0002;
    const ccomplex_t IT_0170 = (-6)*IT_0169;
    const ccomplex_t IT_0171 = s_12*s_23;
    const ccomplex_t IT_0172 = 12*IT_0171;
    const ccomplex_t IT_0173 = IT_0170 + IT_0172;
    const ccomplex_t IT_0174 = 6*IT_0158;
    const ccomplex_t IT_0175 = 0.166666666666667*IT_0163;
    const ccomplex_t IT_0176 = 0.166666666666667*IT_0165;
    const ccomplex_t IT_0177 = (IT_0099*conj(IT_0099) + IT_0158*conj(IT_0158))
      *IT_0159 + conj(IT_0158)*IT_0160*IT_0161 + (IT_0158*IT_0160 + IT_0162
      *IT_0163 + IT_0164*IT_0165)*IT_0166 + (IT_0163*IT_0164 + IT_0162*IT_0165)
      *IT_0167 + 6*conj(IT_0164)*(IT_0162*IT_0168 + 0.166666666666667*IT_0164
      *IT_0173 + IT_0174*IT_0175 + IT_0161*IT_0176) + 6*conj(IT_0162)*(IT_0164
      *IT_0168 + 0.166666666666667*IT_0162*IT_0173 + IT_0161*IT_0175 + IT_0174
      *IT_0176);
    return create_ccomplex_return(IT_0177);
}

