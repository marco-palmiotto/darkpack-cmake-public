#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_sb_1_to_mu_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_sb_1_to_mu_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_mu, 2);
    const ccomplex_t IT_0001 = pow(m_smu_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*conj(N_d1)
      *e_em*m_mu*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = IT_0005*IT_0007;
    const ccomplex_t IT_0011 = m_b*N_d1*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0013 + -IT_0016 +
       0.333333333333333*IT_0021);
    const ccomplex_t IT_0023 = 3*IT_0022;
    const ccomplex_t IT_0024 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0025 = IT_0009*IT_0024;
    const ccomplex_t IT_0026 = IT_0002*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*conj(N_d2)
      *e_em*m_mu*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = m_b*N_d2*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0032 = IT_0010*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0035 = IT_0007*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0038 = IT_0018*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0033 + -IT_0036 +
       0.333333333333333*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = (-3)*IT_0041;
    const ccomplex_t IT_0043 = 0.166666666666667*IT_0042;
    const ccomplex_t IT_0044 = IT_0030*IT_0043;
    const ccomplex_t IT_0045 = IT_0028*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*conj(N_d4)
      *e_em*m_mu*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0051 = IT_0018*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0054 = IT_0007*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_b*N_d4*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0057 = IT_0010*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0052 + (-3)*IT_0055 + 3
      *IT_0058);
    const ccomplex_t IT_0060 = 0.166666666666667*IT_0059;
    const ccomplex_t IT_0061 = IT_0049*IT_0060;
    const ccomplex_t IT_0062 = IT_0047*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0065 = m_b*N_d3*e_em*IT_0003*conj(U_sb_10);
    const ccomplex_t IT_0066 = IT_0010*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0069 = IT_0018*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0072 = IT_0007*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0067 + 0.333333333333333
      *IT_0070 + -IT_0073);
    const ccomplex_t IT_0075 = 3*IT_0074;
    const ccomplex_t IT_0076 = 0.166666666666667*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1.4142135623731)*conj(N_d3)
      *e_em*m_mu*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = IT_0064*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0027 + -IT_0046 + -IT_0063 + -IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0018;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0076*IT_0085;
    const ccomplex_t IT_0087 = IT_0064*IT_0083*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0018;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = IT_0060*IT_0090;
    const ccomplex_t IT_0092 = IT_0047*IT_0088*IT_0091;
    const ccomplex_t IT_0093 = m_mu*IT_0027;
    const ccomplex_t IT_0094 = m_mu*IT_0046;
    const ccomplex_t IT_0095 = m_mu*IT_0063;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0018;
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = IT_0024*IT_0098;
    const ccomplex_t IT_0100 = IT_0002*IT_0096*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0018;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = IT_0043*IT_0103;
    const ccomplex_t IT_0105 = IT_0028*IT_0101*IT_0104;
    const ccomplex_t IT_0106 = m_mu*IT_0081;
    const ccomplex_t IT_0107 = -IT_0087 + -IT_0092 + IT_0093 + IT_0094 +
       IT_0095 + -IT_0100 + -IT_0105 + IT_0106;
    const ccomplex_t IT_0108 = m_mu*s_14;
    const ccomplex_t IT_0109 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0110 = IT_0018*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = m_b*conj(N_d3)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0113 = IT_0010*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + 1.5*IT_0114);
    const ccomplex_t IT_0116 = (-0.333333333333333)*IT_0115;
    const ccomplex_t IT_0117 = IT_0078*IT_0116;
    const ccomplex_t IT_0118 = IT_0064*IT_0083*IT_0117;
    const ccomplex_t IT_0119 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0120 = IT_0018*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = m_b*conj(N_d1)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0123 = IT_0010*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 1.5*IT_0124);
    const ccomplex_t IT_0126 = (-0.333333333333333)*IT_0125;
    const ccomplex_t IT_0127 = IT_0009*IT_0126;
    const ccomplex_t IT_0128 = IT_0002*IT_0096*IT_0127;
    const ccomplex_t IT_0129 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0130 = IT_0018*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = m_b*conj(N_d2)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0133 = IT_0010*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + 1.5*IT_0134);
    const ccomplex_t IT_0136 = (-0.333333333333333)*IT_0135;
    const ccomplex_t IT_0137 = IT_0030*IT_0136;
    const ccomplex_t IT_0138 = IT_0028*IT_0101*IT_0137;
    const ccomplex_t IT_0139 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0140 = IT_0018*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = m_b*conj(N_d4)*e_em*IT_0003*conj(U_sb_00);
    const ccomplex_t IT_0143 = IT_0010*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0141 + 1.5*IT_0144);
    const ccomplex_t IT_0146 = (-0.333333333333333)*IT_0145;
    const ccomplex_t IT_0147 = IT_0049*IT_0146;
    const ccomplex_t IT_0148 = IT_0047*IT_0088*IT_0147;
    const ccomplex_t IT_0149 = IT_0098*IT_0126;
    const ccomplex_t IT_0150 = IT_0002*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = m_mu*IT_0151;
    const ccomplex_t IT_0153 = IT_0103*IT_0136;
    const ccomplex_t IT_0154 = IT_0028*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = m_mu*IT_0155;
    const ccomplex_t IT_0157 = IT_0085*IT_0116;
    const ccomplex_t IT_0158 = IT_0064*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = m_mu*IT_0159;
    const ccomplex_t IT_0161 = IT_0090*IT_0146;
    const ccomplex_t IT_0162 = IT_0047*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = m_mu*IT_0163;
    const ccomplex_t IT_0165 = -IT_0118 + -IT_0128 + -IT_0138 + -IT_0148 +
       IT_0152 + IT_0156 + IT_0160 + IT_0164;
    const ccomplex_t IT_0166 = m_b*s_13;
    const ccomplex_t IT_0167 = -IT_0151 + -IT_0155 + -IT_0159 + -IT_0163;
    const ccomplex_t IT_0168 = m_b*m_mu*IT_0001;
    const ccomplex_t IT_0169 = s_34*IT_0001;
    const ccomplex_t IT_0170 = (-6)*IT_0169;
    const ccomplex_t IT_0171 = s_13*s_14;
    const ccomplex_t IT_0172 = 12*IT_0171;
    const ccomplex_t IT_0173 = IT_0170 + IT_0172;
    const ccomplex_t IT_0174 = m_b*m_mu;
    const ccomplex_t IT_0175 = 6*s_34;
    const ccomplex_t IT_0176 = -conj(IT_0082);
    const ccomplex_t IT_0177 = (-6)*conj(IT_0107);
    const ccomplex_t IT_0178 = 6*IT_0082*(conj(IT_0107)*IT_0108 + -conj
      (IT_0165)*IT_0166 + -conj(IT_0167)*IT_0168 + 0.166666666666667*conj
      (IT_0082)*IT_0173) + (-6)*IT_0107*(IT_0166*conj(IT_0167) + conj(IT_0165)
      *IT_0174 + (-0.166666666666667)*conj(IT_0107)*IT_0175 + IT_0108*IT_0176) +
       6*IT_0167*(IT_0108*conj(IT_0165) + 0.166666666666667*conj(IT_0167)
      *IT_0173 + IT_0168*IT_0176 + 0.166666666666667*IT_0166*IT_0177) + 6
      *IT_0165*(IT_0108*conj(IT_0167) + 0.166666666666667*conj(IT_0165)*IT_0175 
      + IT_0166*IT_0176 + 0.166666666666667*IT_0174*IT_0177);
    return create_ccomplex_return(IT_0178);
}

