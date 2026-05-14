#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sG_to_t_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sG_to_t_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0000*IT_0002*IT_0004*conj(U_st_00);
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = pow(m_t, 2);
    const ccomplex_t IT_0010 = pow(m_sG, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0000*IT_0002*IT_0004*conj(U_st_01);
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0009 + IT_0010 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0013 + IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_10);
    const ccomplex_t IT_0023 = V_tb*e_em*conj(U_Wm1)*U_sb_00;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = IT_0002*IT_0004;
    const ccomplex_t IT_0026 = m_b*conj(U_d1)*V_tb*e_em*IT_0000*U_sb_10;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0024 + (-0.5)*IT_0028);
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = pow(m_C_1, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0009 + IT_0031 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_11);
    const ccomplex_t IT_0036 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0037 = IT_0004*IT_0036;
    const ccomplex_t IT_0038 = m_b*conj(U_d1)*V_tb*e_em*IT_0000*U_sb_11;
    const ccomplex_t IT_0039 = IT_0025*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + (-0.5)*IT_0040);
    const ccomplex_t IT_0042 = IT_0035*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0009 + IT_0031 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0034 + IT_0045;
    const ccomplex_t IT_0047 = s_13*s_24;
    const ccomplex_t IT_0048 = s_14*s_23;
    const ccomplex_t IT_0049 = s_12*s_34;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = IT_0047 + IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = 8*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_00);
    const ccomplex_t IT_0054 = IT_0029*IT_0053;
    const ccomplex_t IT_0055 = IT_0032*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_sb_01);
    const ccomplex_t IT_0058 = IT_0041*IT_0057;
    const ccomplex_t IT_0059 = IT_0043*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = -IT_0056 + -IT_0060;
    const ccomplex_t IT_0062 = m_b*m_sG*s_13;
    const ccomplex_t IT_0063 = 8*IT_0062;
    const ccomplex_t IT_0064 = sin(beta);
    const ccomplex_t IT_0065 = cpow(IT_0064, -1);
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0000*IT_0004*IT_0065*U_sb_01;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0035*IT_0067;
    const ccomplex_t IT_0069 = IT_0043*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0000*IT_0004*IT_0065*U_sb_00;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0022*IT_0072;
    const ccomplex_t IT_0074 = IT_0032*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = -IT_0070 + -IT_0075;
    const ccomplex_t IT_0077 = m_t*s_24*m_C_1;
    const ccomplex_t IT_0078 = 8*IT_0077;
    const ccomplex_t IT_0079 = IT_0053*IT_0072;
    const ccomplex_t IT_0080 = IT_0032*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = IT_0057*IT_0067;
    const ccomplex_t IT_0083 = IT_0043*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0081 + IT_0084;
    const ccomplex_t IT_0086 = m_b*m_t*m_sG*m_C_1;
    const ccomplex_t IT_0087 = 8*IT_0086;
    const ccomplex_t IT_0088 = 16*IT_0048;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0090 = IT_0006*IT_0089;
    const ccomplex_t IT_0091 = IT_0011*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0094 = IT_0015*IT_0093;
    const ccomplex_t IT_0095 = IT_0018*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = -IT_0092 + -IT_0096;
    const ccomplex_t IT_0098 = m_t*m_sG*s_14;
    const ccomplex_t IT_0099 = 16*IT_0098;
    const ccomplex_t IT_0100 = IT_0004*IT_0065;
    const ccomplex_t IT_0101 = m_t*V_tb*V_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0107 = IT_0004*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0105 + IT_0110;
    const ccomplex_t IT_0112 = IT_0007*IT_0111;
    const ccomplex_t IT_0113 = IT_0011*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = m_t*V_tb*V_u1*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0116 = IT_0100*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0121 = IT_0004*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0119 + IT_0124;
    const ccomplex_t IT_0126 = IT_0016*IT_0125;
    const ccomplex_t IT_0127 = IT_0018*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0114 + IT_0128;
    const ccomplex_t IT_0130 = m_b*s_23*m_C_1;
    const ccomplex_t IT_0131 = 16*IT_0130;
    const ccomplex_t IT_0132 = IT_0089*IT_0111;
    const ccomplex_t IT_0133 = IT_0011*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0093*IT_0125;
    const ccomplex_t IT_0136 = IT_0018*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = -IT_0134 + -IT_0137;
    const ccomplex_t IT_0139 = 16*IT_0086;
    const ccomplex_t IT_0140 = 8*IT_0098;
    const ccomplex_t IT_0141 = 8*IT_0130;
    const ccomplex_t IT_0142 = m_b*m_t*s_12;
    const ccomplex_t IT_0143 = 8*IT_0142;
    const ccomplex_t IT_0144 = m_sG*s_34*m_C_1;
    const ccomplex_t IT_0145 = 8*IT_0144;
    const ccomplex_t IT_0146 = 16*IT_0047;
    const ccomplex_t IT_0147 = 16*IT_0062;
    const ccomplex_t IT_0148 = 16*IT_0077;
    const ccomplex_t IT_0149 = conj(IT_0138)*(IT_0063*IT_0076 + IT_0061
      *IT_0078 + IT_0052*IT_0085 + IT_0046*IT_0087 + IT_0099*IT_0129 + IT_0097
      *IT_0131 + IT_0088*IT_0138 + IT_0021*IT_0139) + conj(IT_0021)*(IT_0046
      *IT_0052 + IT_0061*IT_0063 + IT_0076*IT_0078 + IT_0085*IT_0087 + IT_0021
      *IT_0088 + IT_0097*IT_0099 + IT_0129*IT_0131 + IT_0138*IT_0139) + conj
      (IT_0129)*(IT_0088*IT_0129 + IT_0021*IT_0131 + IT_0099*IT_0138 + IT_0097
      *IT_0139 + IT_0085*IT_0140 + IT_0046*IT_0141 + IT_0076*IT_0143 + IT_0061
      *IT_0145) + conj(IT_0097)*(IT_0088*IT_0097 + IT_0021*IT_0099 + IT_0131
      *IT_0138 + IT_0129*IT_0139 + IT_0046*IT_0140 + IT_0085*IT_0141 + IT_0061
      *IT_0143 + IT_0076*IT_0145) + conj(IT_0076)*(IT_0021*IT_0078 + IT_0063
      *IT_0138 + IT_0061*IT_0139 + IT_0129*IT_0143 + IT_0097*IT_0145 + IT_0076
      *IT_0146 + IT_0085*IT_0147 + IT_0046*IT_0148) + conj(IT_0085)*(IT_0021
      *IT_0087 + IT_0052*IT_0138 + IT_0046*IT_0139 + IT_0129*IT_0140 + IT_0097
      *IT_0141 + IT_0085*IT_0146 + IT_0076*IT_0147 + IT_0061*IT_0148) + conj
      (IT_0046)*(IT_0021*IT_0052 + IT_0087*IT_0138 + IT_0085*IT_0139 + IT_0097
      *IT_0140 + IT_0129*IT_0141 + IT_0046*IT_0146 + IT_0061*IT_0147 + IT_0076
      *IT_0148) + conj(IT_0061)*(IT_0021*IT_0063 + IT_0078*IT_0138 + IT_0076
      *IT_0139 + IT_0097*IT_0143 + IT_0129*IT_0145 + IT_0061*IT_0146 + IT_0046
      *IT_0147 + IT_0085*IT_0148);
    return create_ccomplex_return(IT_0149);
}

