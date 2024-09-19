#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_R_to_u_anti_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_R_to_u_anti_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0000*IT_0004*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0013 = pow(m_u, 2);
    const ccomplex_t IT_0014 = pow(m_su_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0011*IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0000*IT_0004*IT_0008;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0000*IT_0004*IT_0008;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0029*IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0000*IT_0004*IT_0008;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0037*IT_0038*IT_0039;
    const ccomplex_t IT_0041 = cos(theta_W);
    const ccomplex_t IT_0042 = cpow(IT_0041, -1);
    const ccomplex_t IT_0043 = conj(N_B1)*e_em;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = conj(N_W1)*e_em;
    const ccomplex_t IT_0047 = IT_0004*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + 3*IT_0048);
    const ccomplex_t IT_0050 = 0.166666666666667*IT_0049;
    const ccomplex_t IT_0051 = IT_0006*IT_0050;
    const ccomplex_t IT_0052 = IT_0015*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = m_u*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B2)*e_em;
    const ccomplex_t IT_0056 = IT_0042*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W2)*e_em;
    const ccomplex_t IT_0059 = IT_0004*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + 3*IT_0060);
    const ccomplex_t IT_0062 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0063 = IT_0018*IT_0062;
    const ccomplex_t IT_0064 = IT_0023*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = m_u*IT_0065;
    const ccomplex_t IT_0067 = conj(N_B3)*e_em;
    const ccomplex_t IT_0068 = IT_0042*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = conj(N_W3)*e_em;
    const ccomplex_t IT_0071 = IT_0004*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 3*IT_0072);
    const ccomplex_t IT_0074 = 0.166666666666667*IT_0073;
    const ccomplex_t IT_0075 = IT_0026*IT_0074;
    const ccomplex_t IT_0076 = IT_0031*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = m_u*IT_0077;
    const ccomplex_t IT_0079 = conj(N_B4)*e_em;
    const ccomplex_t IT_0080 = IT_0042*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W4)*e_em;
    const ccomplex_t IT_0083 = IT_0004*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0081 + 3*IT_0084);
    const ccomplex_t IT_0086 = 0.166666666666667*IT_0085;
    const ccomplex_t IT_0087 = IT_0034*IT_0086;
    const ccomplex_t IT_0088 = IT_0039*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = m_u*IT_0089;
    const ccomplex_t IT_0091 = -IT_0016 + -IT_0024 + -IT_0032 + -IT_0040 +
       IT_0054 + IT_0066 + IT_0078 + IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0000*IT_0004;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0000*IT_0002*IT_0004*IT_0007;
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = pow(m_W, 2);
    const ccomplex_t IT_0098 = pow(m_sd_R, 2);
    const ccomplex_t IT_0099 = cpow((-2)*s_12 + -IT_0014 + IT_0097 + -IT_0098 
      + -reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0096*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0000*IT_0001*IT_0004*IT_0008;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0000*IT_0004;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_12 + -IT_0014 + -IT_0098 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = -IT_0101 + -IT_0109;
    const ccomplex_t IT_0111 = -IT_0053 + -IT_0065 + -IT_0077 + -IT_0089;
    const ccomplex_t IT_0112 = pow(m_sG, 2);
    const ccomplex_t IT_0113 = cpow((-2)*s_13 + IT_0013 + IT_0014 + -IT_0112 +
       reg_prop, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0115 = cpow(IT_0114, 2);
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*m_sG*IT_0115;
    const ccomplex_t IT_0117 = IT_0113*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0042;
    const ccomplex_t IT_0119 = (-0.333333333333333)*IT_0118;
    const ccomplex_t IT_0120 = IT_0050*IT_0119;
    const ccomplex_t IT_0121 = IT_0012*IT_0015*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0042;
    const ccomplex_t IT_0123 = (-0.333333333333333)*IT_0122;
    const ccomplex_t IT_0124 = IT_0062*IT_0123;
    const ccomplex_t IT_0125 = IT_0022*IT_0023*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0042;
    const ccomplex_t IT_0127 = (-0.333333333333333)*IT_0126;
    const ccomplex_t IT_0128 = IT_0074*IT_0127;
    const ccomplex_t IT_0129 = IT_0030*IT_0031*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0042;
    const ccomplex_t IT_0131 = (-0.333333333333333)*IT_0130;
    const ccomplex_t IT_0132 = IT_0086*IT_0131;
    const ccomplex_t IT_0133 = IT_0038*IT_0039*IT_0132;
    const ccomplex_t IT_0134 = IT_0010*IT_0119;
    const ccomplex_t IT_0135 = IT_0015*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = m_u*IT_0136;
    const ccomplex_t IT_0138 = IT_0020*IT_0123;
    const ccomplex_t IT_0139 = IT_0023*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_u*IT_0140;
    const ccomplex_t IT_0142 = IT_0028*IT_0127;
    const ccomplex_t IT_0143 = IT_0031*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = m_u*IT_0144;
    const ccomplex_t IT_0146 = IT_0036*IT_0131;
    const ccomplex_t IT_0147 = IT_0039*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = m_u*IT_0148;
    const ccomplex_t IT_0150 = 0.166666666666667*IT_0117 + -IT_0121 + -IT_0125
       + -IT_0129 + -IT_0133 + IT_0137 + IT_0141 + IT_0145 + IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0000*IT_0004;
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0095*IT_0152;
    const ccomplex_t IT_0154 = IT_0099*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0000*IT_0002*IT_0004*IT_0007;
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = IT_0105*IT_0157;
    const ccomplex_t IT_0159 = IT_0107*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (-0.5)*IT_0117 + -IT_0155 + -IT_0160;
    const ccomplex_t IT_0162 = -IT_0136 + -IT_0140 + -IT_0144 + -IT_0148;
    const ccomplex_t IT_0163 = m_u*s_14;
    const ccomplex_t IT_0164 = ((18*conj(IT_0091) + 6*conj(IT_0110))*IT_0111 +
       (18*IT_0091 + 6*IT_0110)*conj(IT_0111) + (18*conj(IT_0150) + 6*conj
      (IT_0161))*IT_0162 + (18*IT_0150 + 6*IT_0161)*conj(IT_0162))*IT_0163;
    const ccomplex_t IT_0165 = m_d*m_u*IT_0014;
    const ccomplex_t IT_0166 = s_34*IT_0014;
    const ccomplex_t IT_0167 = (-18)*IT_0166;
    const ccomplex_t IT_0168 = s_13*s_14;
    const ccomplex_t IT_0169 = 36*IT_0168;
    const ccomplex_t IT_0170 = IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = 6*s_34;
    const ccomplex_t IT_0172 = 18*s_34;
    const ccomplex_t IT_0173 = m_d*m_u;
    const ccomplex_t IT_0174 = (-6)*conj(IT_0161);
    const ccomplex_t IT_0175 = (-18)*conj(IT_0150);
    const ccomplex_t IT_0176 = IT_0174 + IT_0175;
    const ccomplex_t IT_0177 = (-18)*conj(IT_0091);
    const ccomplex_t IT_0178 = (-6)*conj(IT_0110);
    const ccomplex_t IT_0179 = IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = (-6)*IT_0110;
    const ccomplex_t IT_0181 = m_d*s_13;
    const ccomplex_t IT_0182 = IT_0164 + conj(IT_0111)*((-18)*IT_0162*IT_0165 
      + IT_0111*IT_0170) + conj(IT_0162)*((-18)*IT_0111*IT_0165 + IT_0162
      *IT_0170) + (conj(IT_0091)*IT_0110 + IT_0091*conj(IT_0110) + conj(IT_0150)
      *IT_0161 + IT_0150*conj(IT_0161))*IT_0171 + (IT_0091*conj(IT_0091) +
       IT_0110*conj(IT_0110) + IT_0150*conj(IT_0150) + IT_0161*conj(IT_0161))
      *IT_0172 + (-6)*IT_0173*((conj(IT_0091) + 3*conj(IT_0110))*IT_0161 + 3
      *IT_0110*conj(IT_0161) + (-0.166666666666667)*IT_0091*IT_0176 + (
      -0.166666666666667)*IT_0150*IT_0179 + (-0.166666666666667)*conj(IT_0150)
      *IT_0180) + (conj(IT_0111)*((-18)*IT_0150 + (-6)*IT_0161) + IT_0111
      *IT_0176 + IT_0162*IT_0179 + conj(IT_0162)*((-18)*IT_0091 + IT_0180))
      *IT_0181;
    return create_ccomplex_return(IT_0182);
}

