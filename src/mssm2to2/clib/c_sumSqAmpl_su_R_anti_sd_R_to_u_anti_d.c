#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_sd_R_to_u_anti_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_sd_R_to_u_anti_d(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_R = param->m_su_R;
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
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_su_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_d*N_d1*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0003*IT_0007*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = IT_0002*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0003*IT_0007*IT_0011;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_d*N_d2*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0003*IT_0007*IT_0011;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_d*N_d3*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0003*IT_0007*IT_0011;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_d*N_d4*e_em
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0016 + -IT_0024 + -IT_0032 + -IT_0040;
    const ccomplex_t IT_0042 = cos(theta_W);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0043;
    const ccomplex_t IT_0045 = (-0.333333333333333)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0043;
    const ccomplex_t IT_0047 = 0.666666666666667*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = IT_0030*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0043;
    const ccomplex_t IT_0052 = 0.666666666666667*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0043;
    const ccomplex_t IT_0054 = (-0.333333333333333)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = IT_0002*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0043;
    const ccomplex_t IT_0059 = (-0.333333333333333)*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0043;
    const ccomplex_t IT_0061 = 0.666666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = IT_0038*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0043;
    const ccomplex_t IT_0066 = 0.666666666666667*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0043;
    const ccomplex_t IT_0068 = (-0.333333333333333)*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = IT_0022*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0050 + -IT_0057 + -IT_0064 + -IT_0071;
    const ccomplex_t IT_0073 = m_d*m_u*IT_0001;
    const ccomplex_t IT_0074 = s_34*IT_0001;
    const ccomplex_t IT_0075 = (-18)*IT_0074;
    const ccomplex_t IT_0076 = s_13*s_14;
    const ccomplex_t IT_0077 = 36*IT_0076;
    const ccomplex_t IT_0078 = IT_0075 + IT_0077;
    const ccomplex_t IT_0079 = m_u*IT_0032;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0081 = IT_0020*IT_0066;
    const ccomplex_t IT_0082 = IT_0022*IT_0080*IT_0081;
    const ccomplex_t IT_0083 = IT_0028*IT_0047;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0085 = IT_0030*IT_0083*IT_0084;
    const ccomplex_t IT_0086 = IT_0036*IT_0061;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0088 = IT_0038*IT_0086*IT_0087;
    const ccomplex_t IT_0089 = m_u*IT_0016;
    const ccomplex_t IT_0090 = m_u*IT_0024;
    const ccomplex_t IT_0091 = m_u*IT_0040;
    const ccomplex_t IT_0092 = IT_0009*IT_0052;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0094 = IT_0002*IT_0092*IT_0093;
    const ccomplex_t IT_0095 = IT_0079 + -IT_0082 + -IT_0085 + -IT_0088 +
       IT_0089 + IT_0090 + IT_0091 + -IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0003*IT_0004*IT_0007*IT_0011;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1.4142135623731)*m_d*m_u*V_ud
      *e_em*IT_0003*IT_0005*IT_0007*IT_0011;
    const ccomplex_t IT_0099 = 0.5*IT_0098;
    const ccomplex_t IT_0100 = IT_0097*IT_0099;
    const ccomplex_t IT_0101 = pow(m_sd_R, 2);
    const ccomplex_t IT_0102 = cpow((-2)*s_12 + -IT_0001 + -IT_0101 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0100*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0003*IT_0005*IT_0007*IT_0010;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0099*IT_0107;
    const ccomplex_t IT_0109 = IT_0102*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = IT_0026*IT_0045;
    const ccomplex_t IT_0113 = IT_0030*IT_0084*IT_0112;
    const ccomplex_t IT_0114 = IT_0034*IT_0059;
    const ccomplex_t IT_0115 = IT_0038*IT_0087*IT_0114;
    const ccomplex_t IT_0116 = IT_0018*IT_0068;
    const ccomplex_t IT_0117 = IT_0022*IT_0080*IT_0116;
    const ccomplex_t IT_0118 = pow(m_sG, 2);
    const ccomplex_t IT_0119 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0118 +
       reg_prop, -1);
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0121 = cpow(IT_0120, 2);
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0119*IT_0122;
    const ccomplex_t IT_0124 = m_u*IT_0123;
    const ccomplex_t IT_0125 = m_u*IT_0071;
    const ccomplex_t IT_0126 = IT_0013*IT_0054;
    const ccomplex_t IT_0127 = IT_0002*IT_0093*IT_0126;
    const ccomplex_t IT_0128 = m_u*IT_0057;
    const ccomplex_t IT_0129 = m_u*IT_0050;
    const ccomplex_t IT_0130 = m_u*IT_0064;
    const ccomplex_t IT_0131 = -IT_0113 + -IT_0115 + -IT_0117 + (
      -0.166666666666667)*IT_0124 + IT_0125 + -IT_0127 + IT_0128 + IT_0129 +
       IT_0130;
    const ccomplex_t IT_0132 = 6*s_34;
    const ccomplex_t IT_0133 = 18*s_34;
    const ccomplex_t IT_0134 = (-18)*conj(IT_0041);
    const ccomplex_t IT_0135 = m_u*s_14;
    const ccomplex_t IT_0136 = (-0.5)*IT_0123;
    const ccomplex_t IT_0137 = 0.5*IT_0124;
    const ccomplex_t IT_0138 = 6*IT_0105;
    const ccomplex_t IT_0139 = 18*IT_0095;
    const ccomplex_t IT_0140 = IT_0138 + IT_0139;
    const ccomplex_t IT_0141 = 6*IT_0137;
    const ccomplex_t IT_0142 = 6*conj(IT_0111);
    const ccomplex_t IT_0143 = 18*conj(IT_0131);
    const ccomplex_t IT_0144 = (-6)*IT_0137;
    const ccomplex_t IT_0145 = 0.166666666666667*IT_0123;
    const ccomplex_t IT_0146 = (-18)*IT_0131;
    const ccomplex_t IT_0147 = (-6)*IT_0111;
    const ccomplex_t IT_0148 = IT_0041*(18*conj(IT_0095) + 6*conj(IT_0105)) +
       IT_0136*((-18)*IT_0111 + 18*conj(IT_0111) + (-6)*IT_0131 + 6*conj(IT_0131
      ) + (-36)*IT_0137) + conj(IT_0041)*IT_0140 + 6*conj(IT_0072)*(IT_0111 + 3
      *IT_0131 + 0.166666666666667*IT_0141) + IT_0072*(IT_0142 + IT_0143 +
       IT_0144) + IT_0145*((-12)*IT_0137 + IT_0142 + IT_0143 + IT_0146 + IT_0147);
    const ccomplex_t IT_0149 = -IT_0072;
    const ccomplex_t IT_0150 = (-6)*IT_0074;
    const ccomplex_t IT_0151 = 12*IT_0076;
    const ccomplex_t IT_0152 = IT_0150 + IT_0151;
    const ccomplex_t IT_0153 = (-0.166666666666667)*IT_0136;
    const ccomplex_t IT_0154 = m_d*s_13;
    const ccomplex_t IT_0155 = (-18)*IT_0095;
    const ccomplex_t IT_0156 = (-6)*IT_0105;
    const ccomplex_t IT_0157 = IT_0155 + IT_0156;
    const ccomplex_t IT_0158 = (-18)*conj(IT_0095);
    const ccomplex_t IT_0159 = (-6)*conj(IT_0105);
    const ccomplex_t IT_0160 = IT_0158 + IT_0159;
    const ccomplex_t IT_0161 = 6*IT_0095;
    const ccomplex_t IT_0162 = 18*IT_0105;
    const ccomplex_t IT_0163 = (-18)*conj(IT_0105);
    const ccomplex_t IT_0164 = (-6)*conj(IT_0095);
    const ccomplex_t IT_0165 = IT_0161 + IT_0162 + IT_0163 + IT_0164;
    const ccomplex_t IT_0166 = (-6)*IT_0041*(conj(IT_0111) + 3*conj(IT_0131) +
       (-0.166666666666667)*IT_0141) + conj(IT_0041)*(IT_0144 + IT_0146 +
       IT_0147) + conj(IT_0072)*IT_0157 + IT_0145*(IT_0140 + IT_0158 + IT_0159) 
      + IT_0072*IT_0160 + (-6)*IT_0153*IT_0165;
    const ccomplex_t IT_0167 = m_d*m_u;
    const ccomplex_t IT_0168 = (-0.166666666666667)*conj(IT_0131);
    const ccomplex_t IT_0169 = (-18)*s_34;
    const ccomplex_t IT_0170 = IT_0041*((-18)*conj(IT_0072)*IT_0073 + conj
      (IT_0041)*IT_0078) + (conj(IT_0095)*IT_0105 + IT_0095*conj(IT_0105) + conj
      (IT_0111)*IT_0131 + IT_0111*conj(IT_0131))*IT_0132 + (IT_0095*conj(IT_0095
      ) + IT_0105*conj(IT_0105) + IT_0111*conj(IT_0111) + IT_0131*conj(IT_0131))
      *IT_0133 + IT_0072*(conj(IT_0072)*IT_0078 + IT_0073*IT_0134) + IT_0135
      *IT_0148 + IT_0145*(IT_0073*(18*IT_0041 + IT_0134) + IT_0078*(conj(IT_0072
      ) + -IT_0145 + IT_0149)) + 6*IT_0136*((IT_0041 + -conj(IT_0041))*IT_0073 +
       (-0.333333333333333)*((-0.5)*conj(IT_0072) + IT_0145 + (-0.5)*IT_0149)
      *IT_0152 + IT_0078*IT_0153) + IT_0154*IT_0166 + (-6)*IT_0167*((IT_0095 + 3
      *IT_0105)*conj(IT_0111) + (-0.166666666666667)*IT_0131*IT_0160 + (
      -0.166666666666667)*IT_0111*(IT_0163 + IT_0164) + (-0.166666666666667)
      *IT_0137*IT_0165 + IT_0157*IT_0168) + (-6)*IT_0137*(s_34*IT_0131 + (
      -0.166666666666667)*conj(IT_0111)*IT_0133 + IT_0132*IT_0168 + (
      -0.166666666666667)*(IT_0111 + IT_0137)*IT_0169);
    return create_ccomplex_return(IT_0170);
}

