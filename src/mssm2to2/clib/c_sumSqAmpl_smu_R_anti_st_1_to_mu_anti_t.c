#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_st_1_to_mu_anti_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_st_1_to_mu_anti_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0004*IT_0016;
    const ccomplex_t IT_0018 = m_t*N_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0011 + 3*IT_0014 + 3
      *IT_0020);
    const ccomplex_t IT_0022 = 0.166666666666667*IT_0021;
    const ccomplex_t IT_0023 = IT_0006*IT_0022;
    const ccomplex_t IT_0024 = pow(m_mu, 2);
    const ccomplex_t IT_0025 = pow(m_smu_R, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_13 + IT_0024 + IT_0025 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = N_B2*e_em*conj(U_st_00);
    const ccomplex_t IT_0032 = IT_0008*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = N_W2*e_em*conj(U_st_00);
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = m_t*N_u2*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0038 = IT_0017*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0033 + 3*IT_0036 + 3
      *IT_0039);
    const ccomplex_t IT_0041 = 0.166666666666667*IT_0040;
    const ccomplex_t IT_0042 = IT_0030*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0024 + IT_0025 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0049 = IT_0008*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0052 = IT_0004*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = m_t*N_u4*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0055 = IT_0017*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0050 + 3*IT_0053 + 3
      *IT_0056);
    const ccomplex_t IT_0058 = 0.166666666666667*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0024 + IT_0025 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0024 + IT_0025 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0064 = N_B3*e_em*conj(U_st_00);
    const ccomplex_t IT_0065 = IT_0008*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_W3*e_em*conj(U_st_00);
    const ccomplex_t IT_0068 = IT_0004*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = m_t*N_u3*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0071 = IT_0017*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0066 + 3*IT_0069 + 3
      *IT_0072);
    const ccomplex_t IT_0074 = 0.166666666666667*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = IT_0063*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = -IT_0028 + -IT_0045 + -IT_0062 + -IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0008;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = IT_0022*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0085 = IT_0026*IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0008;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0041*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0090 = IT_0043*IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0008;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = IT_0074*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0095 = IT_0063*IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0008;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = IT_0058*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0100 = IT_0060*IT_0098*IT_0099;
    const ccomplex_t IT_0101 = m_mu*IT_0028;
    const ccomplex_t IT_0102 = m_mu*IT_0045;
    const ccomplex_t IT_0103 = m_mu*IT_0062;
    const ccomplex_t IT_0104 = m_mu*IT_0079;
    const ccomplex_t IT_0105 = -IT_0085 + -IT_0090 + -IT_0095 + -IT_0100 +
       IT_0101 + IT_0102 + IT_0103 + IT_0104;
    const ccomplex_t IT_0106 = m_mu*s_14;
    const ccomplex_t IT_0107 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0108 = IT_0008*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = m_t*conj(N_u1)*e_em*IT_0000*conj(U_st_00);
    const ccomplex_t IT_0111 = IT_0017*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + (-0.75)*IT_0112);
    const ccomplex_t IT_0114 = 0.666666666666667*IT_0113;
    const ccomplex_t IT_0115 = IT_0006*IT_0114;
    const ccomplex_t IT_0116 = IT_0026*IT_0084*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B2)*e_em*conj(U_st_10);
    const ccomplex_t IT_0118 = IT_0008*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = m_t*conj(N_u2)*e_em*IT_0000*conj(U_st_00);
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0119 + (-0.75)*IT_0122);
    const ccomplex_t IT_0124 = 0.666666666666667*IT_0123;
    const ccomplex_t IT_0125 = IT_0030*IT_0124;
    const ccomplex_t IT_0126 = IT_0043*IT_0089*IT_0125;
    const ccomplex_t IT_0127 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0128 = IT_0008*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = m_t*conj(N_u4)*e_em*IT_0000*conj(U_st_00);
    const ccomplex_t IT_0131 = IT_0017*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + (-0.75)*IT_0132);
    const ccomplex_t IT_0134 = 0.666666666666667*IT_0133;
    const ccomplex_t IT_0135 = IT_0047*IT_0134;
    const ccomplex_t IT_0136 = IT_0060*IT_0099*IT_0135;
    const ccomplex_t IT_0137 = IT_0082*IT_0114;
    const ccomplex_t IT_0138 = IT_0026*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = m_mu*IT_0139;
    const ccomplex_t IT_0141 = IT_0087*IT_0124;
    const ccomplex_t IT_0142 = IT_0043*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = m_mu*IT_0143;
    const ccomplex_t IT_0145 = conj(N_B3)*e_em*conj(U_st_10);
    const ccomplex_t IT_0146 = IT_0008*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = m_t*conj(N_u3)*e_em*IT_0000*conj(U_st_00);
    const ccomplex_t IT_0149 = IT_0017*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + (-0.75)*IT_0150);
    const ccomplex_t IT_0152 = 0.666666666666667*IT_0151;
    const ccomplex_t IT_0153 = IT_0092*IT_0152;
    const ccomplex_t IT_0154 = IT_0063*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = m_mu*IT_0155;
    const ccomplex_t IT_0157 = IT_0097*IT_0134;
    const ccomplex_t IT_0158 = IT_0060*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = m_mu*IT_0159;
    const ccomplex_t IT_0161 = IT_0076*IT_0152;
    const ccomplex_t IT_0162 = IT_0063*IT_0094*IT_0161;
    const ccomplex_t IT_0163 = -IT_0116 + -IT_0126 + -IT_0136 + IT_0140 +
       IT_0144 + IT_0156 + IT_0160 + -IT_0162;
    const ccomplex_t IT_0164 = m_t*s_13;
    const ccomplex_t IT_0165 = -IT_0139 + -IT_0143 + -IT_0155 + -IT_0159;
    const ccomplex_t IT_0166 = m_t*m_mu*IT_0025;
    const ccomplex_t IT_0167 = s_34*IT_0025;
    const ccomplex_t IT_0168 = (-6)*IT_0167;
    const ccomplex_t IT_0169 = s_13*s_14;
    const ccomplex_t IT_0170 = 12*IT_0169;
    const ccomplex_t IT_0171 = IT_0168 + IT_0170;
    const ccomplex_t IT_0172 = m_t*m_mu;
    const ccomplex_t IT_0173 = 6*s_34;
    const ccomplex_t IT_0174 = -conj(IT_0080);
    const ccomplex_t IT_0175 = (-6)*conj(IT_0105);
    const ccomplex_t IT_0176 = 6*IT_0080*(conj(IT_0105)*IT_0106 + -conj
      (IT_0163)*IT_0164 + -conj(IT_0165)*IT_0166 + 0.166666666666667*conj
      (IT_0080)*IT_0171) + (-6)*IT_0105*(IT_0164*conj(IT_0165) + conj(IT_0163)
      *IT_0172 + (-0.166666666666667)*conj(IT_0105)*IT_0173 + IT_0106*IT_0174) +
       6*IT_0165*(IT_0106*conj(IT_0163) + 0.166666666666667*conj(IT_0165)
      *IT_0171 + IT_0166*IT_0174 + 0.166666666666667*IT_0164*IT_0175) + 6
      *IT_0163*(IT_0106*conj(IT_0165) + 0.166666666666667*conj(IT_0163)*IT_0173 
      + IT_0164*IT_0174 + 0.166666666666667*IT_0172*IT_0175);
    return create_ccomplex_return(IT_0176);
}

