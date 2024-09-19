#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_su_R_to_tau_anti_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_su_R_to_tau_anti_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
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
    const creal_t m_stau_2 = param->m_stau_2;
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
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*IT_0010;
    const ccomplex_t IT_0012 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0004*IT_0016;
    const ccomplex_t IT_0018 = conj(N_d1)*e_em*m_tau*IT_0000*U_stau_11;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0011 + IT_0014 + 
      -IT_0020);
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0006*IT_0022;
    const ccomplex_t IT_0024 = pow(m_tau, 2);
    const ccomplex_t IT_0025 = pow(m_stau_2, 2);
    const ccomplex_t IT_0026 = cpow(s_13 + (-0.5)*IT_0024 + (-0.5)*IT_0025 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0032 = IT_0008*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_d2)*e_em*m_tau*IT_0000*U_stau_11;
    const ccomplex_t IT_0038 = IT_0017*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0033 + IT_0036 + 
      -IT_0039);
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = IT_0030*IT_0041;
    const ccomplex_t IT_0043 = cpow(s_13 + (-0.5)*IT_0024 + (-0.5)*IT_0025 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0049 = IT_0008*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0052 = IT_0004*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_d3)*e_em*m_tau*IT_0000*U_stau_11;
    const ccomplex_t IT_0055 = IT_0017*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0050 + IT_0053 + 
      -IT_0056);
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0047*IT_0058;
    const ccomplex_t IT_0060 = cpow(s_13 + (-0.5)*IT_0024 + (-0.5)*IT_0025 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0066 = IT_0008*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0069 = IT_0004*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = conj(N_d4)*e_em*m_tau*IT_0000*U_stau_11;
    const ccomplex_t IT_0072 = IT_0017*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0067 + IT_0070 + 
      -IT_0073);
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = IT_0064*IT_0075;
    const ccomplex_t IT_0077 = cpow(s_13 + (-0.5)*IT_0024 + (-0.5)*IT_0025 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = 0.5*IT_0028 + 0.5*IT_0045 + 0.5*IT_0062 + 0.5
      *IT_0079;
    const ccomplex_t IT_0081 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0082 = IT_0008*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = N_d2*e_em*m_tau*IT_0000*U_stau_01;
    const ccomplex_t IT_0085 = IT_0017*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + 0.5*IT_0086);
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0030*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0091 = IT_0043*IT_0089*IT_0090;
    const ccomplex_t IT_0092 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0093 = IT_0008*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = N_d3*e_em*m_tau*IT_0000*U_stau_01;
    const ccomplex_t IT_0096 = IT_0017*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + 0.5*IT_0097);
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0047*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0102 = IT_0060*IT_0100*IT_0101;
    const ccomplex_t IT_0103 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0104 = IT_0008*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = N_d4*e_em*m_tau*IT_0000*U_stau_01;
    const ccomplex_t IT_0107 = IT_0017*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + 0.5*IT_0108);
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0064*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0113 = IT_0077*IT_0111*IT_0112;
    const ccomplex_t IT_0114 = m_tau*IT_0028;
    const ccomplex_t IT_0115 = m_tau*IT_0045;
    const ccomplex_t IT_0116 = m_tau*IT_0062;
    const ccomplex_t IT_0117 = m_tau*IT_0079;
    const ccomplex_t IT_0118 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0119 = IT_0008*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = N_d1*e_em*m_tau*IT_0000*U_stau_01;
    const ccomplex_t IT_0122 = IT_0017*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + 0.5*IT_0123);
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = IT_0006*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0128 = IT_0026*IT_0126*IT_0127;
    const ccomplex_t IT_0129 = 0.5*IT_0091 + 0.5*IT_0102 + 0.5*IT_0113 + (-0.5
      )*IT_0114 + (-0.5)*IT_0115 + (-0.5)*IT_0116 + (-0.5)*IT_0117 + 0.5*IT_0128;
    const ccomplex_t IT_0130 = s_14*m_tau;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0008;
    const ccomplex_t IT_0132 = 0.666666666666667*IT_0131;
    const ccomplex_t IT_0133 = IT_0088*IT_0132;
    const ccomplex_t IT_0134 = IT_0043*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = m_tau*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0008;
    const ccomplex_t IT_0138 = 0.666666666666667*IT_0137;
    const ccomplex_t IT_0139 = IT_0099*IT_0138;
    const ccomplex_t IT_0140 = IT_0060*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = m_tau*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0008;
    const ccomplex_t IT_0144 = 0.666666666666667*IT_0143;
    const ccomplex_t IT_0145 = IT_0110*IT_0144;
    const ccomplex_t IT_0146 = IT_0077*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = m_tau*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0008;
    const ccomplex_t IT_0150 = 0.666666666666667*IT_0149;
    const ccomplex_t IT_0151 = IT_0022*IT_0150;
    const ccomplex_t IT_0152 = IT_0026*IT_0127*IT_0151;
    const ccomplex_t IT_0153 = IT_0041*IT_0132;
    const ccomplex_t IT_0154 = IT_0043*IT_0090*IT_0153;
    const ccomplex_t IT_0155 = IT_0058*IT_0138;
    const ccomplex_t IT_0156 = IT_0060*IT_0101*IT_0155;
    const ccomplex_t IT_0157 = IT_0075*IT_0144;
    const ccomplex_t IT_0158 = IT_0077*IT_0112*IT_0157;
    const ccomplex_t IT_0159 = IT_0125*IT_0150;
    const ccomplex_t IT_0160 = IT_0026*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = m_tau*IT_0161;
    const ccomplex_t IT_0163 = (-0.5)*IT_0136 + (-0.5)*IT_0142 + (-0.5)
      *IT_0148 + 0.5*IT_0152 + 0.5*IT_0154 + 0.5*IT_0156 + 0.5*IT_0158 + (-0.5)
      *IT_0162;
    const ccomplex_t IT_0164 = m_u*s_13;
    const ccomplex_t IT_0165 = 0.5*IT_0135 + 0.5*IT_0141 + 0.5*IT_0147 + 0.5
      *IT_0161;
    const ccomplex_t IT_0166 = m_u*m_tau*IT_0025;
    const ccomplex_t IT_0167 = s_34*IT_0025;
    const ccomplex_t IT_0168 = (-6)*IT_0167;
    const ccomplex_t IT_0169 = s_13*s_14;
    const ccomplex_t IT_0170 = 12*IT_0169;
    const ccomplex_t IT_0171 = IT_0168 + IT_0170;
    const ccomplex_t IT_0172 = 6*s_34;
    const ccomplex_t IT_0173 = m_u*m_tau;
    const ccomplex_t IT_0174 = (-6)*conj(IT_0129);
    const ccomplex_t IT_0175 = (-6)*IT_0129;
    const ccomplex_t IT_0176 = -IT_0080;
    const ccomplex_t IT_0177 = 6*conj(IT_0080)*(IT_0129*IT_0130 + -IT_0163
      *IT_0164 + -IT_0165*IT_0166 + 0.166666666666667*IT_0080*IT_0171) + 6*conj
      (IT_0129)*(IT_0080*IT_0130 + 0.166666666666667*IT_0129*IT_0172) + (IT_0164
      *IT_0165 + IT_0163*IT_0173)*IT_0174 + 6*conj(IT_0163)*(IT_0130*IT_0165 +
       0.166666666666667*IT_0163*IT_0172 + 0.166666666666667*IT_0173*IT_0175 +
       IT_0164*IT_0176) + 6*conj(IT_0165)*(IT_0130*IT_0163 + 0.166666666666667
      *IT_0165*IT_0171 + 0.166666666666667*IT_0164*IT_0175 + IT_0166*IT_0176);
    return create_ccomplex_return(IT_0177);
}

