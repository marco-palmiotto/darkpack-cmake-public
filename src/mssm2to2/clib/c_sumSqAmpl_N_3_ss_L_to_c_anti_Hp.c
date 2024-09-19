#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_ss_L_to_c_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_ss_L_to_c_anti_Hp(
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
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_N_3, 2);
    const ccomplex_t IT_0002 = pow(m_ss_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0005*IT_0007*IT_0009;
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0005*IT_0006*IT_0009*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = IT_0011*IT_0015;
    const ccomplex_t IT_0017 = IT_0003*IT_0004*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0005*IT_0009*IT_0013;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0009*IT_0012;
    const ccomplex_t IT_0021 = N_d3*e_em*U_Wm1;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (-2)*IT_0022;
    const ccomplex_t IT_0024 = cos(theta_W);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = IT_0012*IT_0025;
    const ccomplex_t IT_0027 = N_B3*U_d1*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = N_W3*U_d1*e_em;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = IT_0029 + IT_0032;
    const ccomplex_t IT_0034 = IT_0023 + IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0019*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0039 = pow(m_c, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0002 + IT_0039 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0037*IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0005*IT_0009*IT_0013;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = N_W3*U_d2*e_em;
    const ccomplex_t IT_0045 = IT_0020*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = N_B3*U_d2*e_em;
    const ccomplex_t IT_0048 = IT_0026*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = IT_0046 + IT_0049;
    const ccomplex_t IT_0051 = N_d3*e_em*U_Wm2;
    const ccomplex_t IT_0052 = IT_0020*IT_0051;
    const ccomplex_t IT_0053 = (-2)*IT_0052;
    const ccomplex_t IT_0054 = IT_0050 + IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0043*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0002 + IT_0039 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0057*IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0009;
    const ccomplex_t IT_0062 = IT_0036*IT_0061;
    const ccomplex_t IT_0063 = IT_0040*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = m_c*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0009;
    const ccomplex_t IT_0067 = IT_0056*IT_0066;
    const ccomplex_t IT_0068 = IT_0059*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = m_c*IT_0069;
    const ccomplex_t IT_0071 = conj(N_B3)*e_em;
    const ccomplex_t IT_0072 = IT_0025*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = conj(N_W3)*e_em;
    const ccomplex_t IT_0075 = IT_0009*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0073 + (-3)*IT_0076);
    const ccomplex_t IT_0078 = 0.166666666666667*IT_0077;
    const ccomplex_t IT_0079 = IT_0015*IT_0078;
    const ccomplex_t IT_0080 = IT_0003*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = m_N_3*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*V_cs*e_em
      *IT_0009*(IT_0000*IT_0005*IT_0007*IT_0012 + (-2)*IT_0006*(m_W*IT_0012 + (
      -0.5)*IT_0005*IT_0013*IT_0039));
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0005*IT_0009*IT_0013;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0001 + IT_0039 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *mu_h*IT_0005*IT_0009;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0025;
    const ccomplex_t IT_0094 = 0.666666666666667*IT_0093;
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_13 + IT_0001 + IT_0039 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = -IT_0017 + IT_0041 + IT_0060 + -IT_0065 + 
      -IT_0070 + IT_0082 + IT_0090 + IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0073 + 3*IT_0076);
    const ccomplex_t IT_0101 = 0.166666666666667*IT_0100;
    const ccomplex_t IT_0102 = IT_0084*IT_0101;
    const ccomplex_t IT_0103 = IT_0088*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0005*IT_0009*IT_0013;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0092*IT_0106;
    const ccomplex_t IT_0108 = IT_0096*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0006*IT_0009;
    const ccomplex_t IT_0111 = conj(N_u3)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = 2*IT_0112;
    const ccomplex_t IT_0114 = IT_0006*IT_0025;
    const ccomplex_t IT_0115 = conj(N_B3)*conj(V_u1)*e_em;
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W3)*conj(V_u1)*e_em;
    const ccomplex_t IT_0119 = IT_0110*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = IT_0117 + IT_0120;
    const ccomplex_t IT_0122 = IT_0113 + IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = IT_0061*IT_0124;
    const ccomplex_t IT_0126 = IT_0038*IT_0040*IT_0125;
    const ccomplex_t IT_0127 = conj(N_u3)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0128 = IT_0110*IT_0127;
    const ccomplex_t IT_0129 = 2*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B3)*conj(V_u2)*e_em;
    const ccomplex_t IT_0131 = IT_0114*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = conj(N_W3)*conj(V_u2)*e_em;
    const ccomplex_t IT_0134 = IT_0110*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = IT_0132 + IT_0135;
    const ccomplex_t IT_0137 = IT_0129 + IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = IT_0066*IT_0139;
    const ccomplex_t IT_0141 = IT_0058*IT_0059*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0005*IT_0007*IT_0009*IT_0012;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = IT_0078*IT_0143;
    const ccomplex_t IT_0145 = IT_0003*IT_0004*IT_0144;
    const ccomplex_t IT_0146 = IT_0019*IT_0124;
    const ccomplex_t IT_0147 = IT_0040*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = m_c*IT_0148;
    const ccomplex_t IT_0150 = IT_0043*IT_0139;
    const ccomplex_t IT_0151 = IT_0059*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = m_c*IT_0152;
    const ccomplex_t IT_0154 = IT_0011*IT_0143;
    const ccomplex_t IT_0155 = IT_0003*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = m_N_3*IT_0156;
    const ccomplex_t IT_0158 = -IT_0104 + -IT_0109 + IT_0126 + IT_0141 +
       IT_0145 + -IT_0149 + -IT_0153 + -IT_0157;
    const ccomplex_t IT_0159 = 6*s_13;
    const ccomplex_t IT_0160 = m_c*m_N_3;
    const ccomplex_t IT_0161 = 6*IT_0099;
    const ccomplex_t IT_0162 = IT_0064 + IT_0069 + -IT_0156;
    const ccomplex_t IT_0163 = m_c*s_12;
    const ccomplex_t IT_0164 = IT_0081 + IT_0148 + IT_0152;
    const ccomplex_t IT_0165 = s_23*m_N_3;
    const ccomplex_t IT_0166 = 6*conj(IT_0099);
    const ccomplex_t IT_0167 = 6*conj(IT_0158);
    const ccomplex_t IT_0168 = m_c*m_N_3*IT_0002;
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

