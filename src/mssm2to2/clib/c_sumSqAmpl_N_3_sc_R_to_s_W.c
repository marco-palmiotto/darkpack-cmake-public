#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sc_R_to_s_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sc_R_to_s_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u1)*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0008 + (-0.5)*IT_0011);
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = pow(m_s, 2);
    const ccomplex_t IT_0016 = pow(m_sc_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u2)*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0025 = IT_0004*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + (-0.5)*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0015 + IT_0016 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0019 + IT_0032;
    const ccomplex_t IT_0034 = m_N_3*IT_0033;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = pow(m_W, -2);
    const ccomplex_t IT_0037 = s_14*s_34*IT_0036;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = s_13 + IT_0038;
    const ccomplex_t IT_0040 = m_N_3*IT_0039;
    const ccomplex_t IT_0041 = pow(m_c, 2);
    const ccomplex_t IT_0042 = pow(m_N_3, 2);
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + -IT_0016 + IT_0041 + -IT_0042 
      + -reg_prop, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0004;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = IT_0043*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = 2*IT_0050;
    const ccomplex_t IT_0052 = s_13*IT_0039;
    const ccomplex_t IT_0053 = s_13*s_14;
    const ccomplex_t IT_0054 = s_34*IT_0042;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0053 + IT_0055;
    const ccomplex_t IT_0057 = IT_0033 + IT_0050;
    const ccomplex_t IT_0058 = (-2)*IT_0019 + (-2)*IT_0032;
    const ccomplex_t IT_0059 = 6*conj(IT_0058);
    const ccomplex_t IT_0060 = m_s*IT_0050;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = pow(s_14, 2);
    const ccomplex_t IT_0063 = IT_0036*IT_0062;
    const ccomplex_t IT_0064 = -IT_0042;
    const ccomplex_t IT_0065 = IT_0063 + IT_0064;
    const ccomplex_t IT_0066 = m_s*IT_0065;
    const ccomplex_t IT_0067 = s_13*IT_0065;
    const ccomplex_t IT_0068 = 6*IT_0058;
    const ccomplex_t IT_0069 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0070 = IT_0004*IT_0069;
    const ccomplex_t IT_0071 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0072 = IT_0004*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0070 + 0.5*IT_0073);
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0006*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0078 = IT_0017*IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0080 = cos(theta_W);
    const ccomplex_t IT_0081 = cpow(IT_0080, -1);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0081;
    const ccomplex_t IT_0083 = 0.666666666666667*IT_0082;
    const ccomplex_t IT_0084 = IT_0047*IT_0083;
    const ccomplex_t IT_0085 = IT_0043*IT_0079*IT_0084;
    const ccomplex_t IT_0086 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0087 = IT_0004*IT_0086;
    const ccomplex_t IT_0088 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0089 = IT_0004*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + 0.5*IT_0090);
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = IT_0021*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0095 = IT_0030*IT_0093*IT_0094;
    const ccomplex_t IT_0096 = IT_0078 + -IT_0085 + IT_0095;
    const ccomplex_t IT_0097 = 6*IT_0096;
    const ccomplex_t IT_0098 = pow(m_W, 2);
    const ccomplex_t IT_0099 = s_14*IT_0098;
    const ccomplex_t IT_0100 = IT_0036*IT_0099;
    const ccomplex_t IT_0101 = (-0.25)*IT_0100;
    const ccomplex_t IT_0102 = s_14 + IT_0101;
    const ccomplex_t IT_0103 = m_s*IT_0102;
    const ccomplex_t IT_0104 = (-4)*IT_0103;
    const ccomplex_t IT_0105 = pow(m_W, 4);
    const ccomplex_t IT_0106 = s_13*IT_0105;
    const ccomplex_t IT_0107 = IT_0036*IT_0106;
    const ccomplex_t IT_0108 = 6*IT_0107;
    const ccomplex_t IT_0109 = s_14*s_34;
    const ccomplex_t IT_0110 = 24*IT_0109;
    const ccomplex_t IT_0111 = s_13*IT_0098;
    const ccomplex_t IT_0112 = (-12)*IT_0111;
    const ccomplex_t IT_0113 = IT_0108 + IT_0110 + IT_0112;
    const ccomplex_t IT_0114 = (-6)*IT_0107;
    const ccomplex_t IT_0115 = IT_0110 + IT_0114;
    const ccomplex_t IT_0116 = s_34*IT_0098;
    const ccomplex_t IT_0117 = IT_0036*IT_0116;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = s_34 + IT_0118;
    const ccomplex_t IT_0120 = m_N_3*IT_0119;
    const ccomplex_t IT_0121 = 2*IT_0120;
    const ccomplex_t IT_0122 = s_14*IT_0015;
    const ccomplex_t IT_0123 = s_13*s_34;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = (-6)*IT_0051;
    const ccomplex_t IT_0127 = pow(s_34, 2);
    const ccomplex_t IT_0128 = IT_0036*IT_0127;
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = IT_0015 + IT_0129;
    const ccomplex_t IT_0131 = m_N_3*IT_0130;
    const ccomplex_t IT_0132 = m_s*IT_0039;
    const ccomplex_t IT_0133 = (-6)*conj(IT_0051);
    const ccomplex_t IT_0134 = s_13*IT_0130;
    const ccomplex_t IT_0135 = (-6)*IT_0096;
    const ccomplex_t IT_0136 = 0.5*IT_0100;
    const ccomplex_t IT_0137 = s_14 + IT_0136;
    const ccomplex_t IT_0138 = m_s*IT_0137;
    const ccomplex_t IT_0139 = 2*IT_0138;
    const ccomplex_t IT_0140 = (-0.25)*IT_0117;
    const ccomplex_t IT_0141 = s_34 + IT_0140;
    const ccomplex_t IT_0142 = m_N_3*IT_0141;
    const ccomplex_t IT_0143 = (-4)*IT_0142;
    const ccomplex_t IT_0144 = (-6)*conj(IT_0096);
    const ccomplex_t IT_0145 = m_s*m_N_3;
    const ccomplex_t IT_0146 = (-2)*IT_0109;
    const ccomplex_t IT_0147 = IT_0111 + IT_0146;
    const ccomplex_t IT_0148 = IT_0036*IT_0147;
    const ccomplex_t IT_0149 = (-6)*IT_0148;
    const ccomplex_t IT_0150 = 12*s_13;
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = 6*IT_0050;
    const ccomplex_t IT_0153 = (-6)*IT_0033;
    const ccomplex_t IT_0154 = (-18)*IT_0061*IT_0145;
    const ccomplex_t IT_0155 = IT_0121*IT_0152;
    const ccomplex_t IT_0156 = (IT_0035 + IT_0096)*IT_0151;
    const ccomplex_t IT_0157 = (IT_0035*IT_0040 + IT_0051*IT_0052 + IT_0056
      *IT_0057)*IT_0059 + (-6)*conj(IT_0058)*(IT_0061*IT_0066 + (
      -0.166666666666667)*IT_0067*IT_0068 + (-0.166666666666667)*IT_0040*IT_0097
      ) + 6*conj(IT_0050)*(0.166666666666667*IT_0056*IT_0068 + IT_0061*IT_0104 +
       0.166666666666667*IT_0050*IT_0113 + 0.166666666666667*IT_0033*IT_0115 + 
      (IT_0035 + 0.166666666666667*IT_0097)*IT_0121 + 0.166666666666667*IT_0125
      *IT_0126) + (IT_0057*IT_0125 + IT_0035*IT_0131 + IT_0061*IT_0132)*IT_0133 
      + conj(IT_0051)*(IT_0052*IT_0068 + IT_0126*IT_0134 + IT_0131*IT_0135) + (
      -6)*conj(IT_0033)*((-0.166666666666667)*IT_0056*IT_0068 + (
      -0.166666666666667)*IT_0033*IT_0113 + (-0.166666666666667)*IT_0050*IT_0115
       + (-0.166666666666667)*IT_0125*IT_0126 + IT_0061*IT_0139 + (IT_0035 + (
      -0.166666666666667)*IT_0135)*IT_0143) + (IT_0051*IT_0131 + IT_0033*IT_0143
      )*IT_0144 + (-6)*conj(IT_0061)*(IT_0058*IT_0066 + -1./6*IT_0126*IT_0132 + 
      -1./6*((-18)*IT_0035 + (-18)*IT_0096)*IT_0145 + -1./6*IT_0061*IT_0151 + 
      -1./6*IT_0104*IT_0152 + -1./6*IT_0139*IT_0153) + conj(IT_0035)*(IT_0040
      *IT_0068 + IT_0126*IT_0131 + IT_0143*IT_0153 + IT_0154 + IT_0155 + IT_0156
      ) + 6*conj(IT_0096)*(IT_0040*IT_0058 + 0.166666666666667*IT_0154 +
       0.166666666666667*IT_0155 + 0.166666666666667*IT_0156);
    return create_ccomplex_return(IT_0157);
}

