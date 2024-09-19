#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_snu_tau_to_anti_A0_nu_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_snu_tau_to_anti_A0_nu_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
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
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0011;
    const ccomplex_t IT_0016 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = -IT_0005 + IT_0009 + IT_0014 + -IT_0017;
    const ccomplex_t IT_0019 = m_N_3*IT_0018;
    const ccomplex_t IT_0020 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0021 = IT_0003*IT_0020;
    const ccomplex_t IT_0022 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0023 = IT_0007*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0025 = IT_0012*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0027 = IT_0015*IT_0026;
    const ccomplex_t IT_0028 = IT_0021 + -IT_0023 + -IT_0025 + IT_0027;
    const ccomplex_t IT_0029 = m_N_3*IT_0028;
    const ccomplex_t IT_0030 = IT_0019 + IT_0029;
    const ccomplex_t IT_0031 = conj(N_B3)*e_em;
    const ccomplex_t IT_0032 = IT_0002*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W3)*e_em;
    const ccomplex_t IT_0035 = IT_0011*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + -IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = pow(m_A0, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0039 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0030*IT_0042;
    const ccomplex_t IT_0044 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0045 = IT_0003*IT_0044;
    const ccomplex_t IT_0046 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0047 = IT_0007*IT_0046;
    const ccomplex_t IT_0048 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0049 = IT_0012*IT_0048;
    const ccomplex_t IT_0050 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0051 = IT_0015*IT_0050;
    const ccomplex_t IT_0052 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0053 = IT_0003*IT_0052;
    const ccomplex_t IT_0054 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0055 = IT_0007*IT_0054;
    const ccomplex_t IT_0056 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0057 = IT_0012*IT_0056;
    const ccomplex_t IT_0058 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0059 = IT_0015*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0045 + 0.5*IT_0047 + 0.5*IT_0049 + (
      -0.5)*IT_0051 + (-0.5)*IT_0053 + 0.5*IT_0055 + 0.5*IT_0057 + (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = m_N_3*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0065 = IT_0007*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0067 = IT_0012*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0069 = IT_0015*IT_0068;
    const ccomplex_t IT_0070 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0071 = IT_0003*IT_0070;
    const ccomplex_t IT_0072 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0073 = IT_0007*IT_0072;
    const ccomplex_t IT_0074 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0075 = IT_0012*IT_0074;
    const ccomplex_t IT_0076 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0077 = IT_0015*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0063 + (-0.5)*IT_0065 + (-0.5)*IT_0067 +
       0.5*IT_0069 + 0.5*IT_0071 + (-0.5)*IT_0073 + (-0.5)*IT_0075 + 0.5*IT_0077;
    const ccomplex_t IT_0079 = m_N_1*IT_0078;
    const ccomplex_t IT_0080 = IT_0061 + IT_0079;
    const ccomplex_t IT_0081 = conj(N_B1)*e_em;
    const ccomplex_t IT_0082 = IT_0002*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W1)*e_em;
    const ccomplex_t IT_0085 = IT_0011*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + -IT_0086);
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = pow(m_N_3, 2);
    const ccomplex_t IT_0090 = cpow((-2)*s_13 + IT_0039 + IT_0089 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0080*IT_0092;
    const ccomplex_t IT_0094 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0095 = IT_0003*IT_0094;
    const ccomplex_t IT_0096 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0097 = IT_0007*IT_0096;
    const ccomplex_t IT_0098 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0099 = IT_0012*IT_0098;
    const ccomplex_t IT_0100 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0101 = IT_0015*IT_0100;
    const ccomplex_t IT_0102 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0103 = IT_0003*IT_0102;
    const ccomplex_t IT_0104 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0105 = IT_0007*IT_0104;
    const ccomplex_t IT_0106 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0107 = IT_0012*IT_0106;
    const ccomplex_t IT_0108 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0109 = IT_0015*IT_0108;
    const ccomplex_t IT_0110 = (-0.5)*IT_0095 + 0.5*IT_0097 + 0.5*IT_0099 + (
      -0.5)*IT_0101 + (-0.5)*IT_0103 + 0.5*IT_0105 + 0.5*IT_0107 + (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = m_N_3*IT_0110;
    const ccomplex_t IT_0112 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0113 = IT_0003*IT_0112;
    const ccomplex_t IT_0114 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0115 = IT_0007*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0117 = IT_0012*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0119 = IT_0015*IT_0118;
    const ccomplex_t IT_0120 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0121 = IT_0003*IT_0120;
    const ccomplex_t IT_0122 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0123 = IT_0007*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0125 = IT_0012*IT_0124;
    const ccomplex_t IT_0126 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0127 = IT_0015*IT_0126;
    const ccomplex_t IT_0128 = 0.5*IT_0113 + (-0.5)*IT_0115 + (-0.5)*IT_0117 +
       0.5*IT_0119 + 0.5*IT_0121 + (-0.5)*IT_0123 + (-0.5)*IT_0125 + 0.5*IT_0127;
    const ccomplex_t IT_0129 = m_N_2*IT_0128;
    const ccomplex_t IT_0130 = IT_0111 + IT_0129;
    const ccomplex_t IT_0131 = conj(N_B2)*e_em;
    const ccomplex_t IT_0132 = IT_0002*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W2)*e_em;
    const ccomplex_t IT_0135 = IT_0011*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + -IT_0136);
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = cpow((-2)*s_13 + IT_0039 + IT_0089 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = IT_0130*IT_0141;
    const ccomplex_t IT_0143 = conj(N_B4)*e_em;
    const ccomplex_t IT_0144 = IT_0002*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W4)*e_em;
    const ccomplex_t IT_0147 = IT_0011*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + -IT_0148);
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = cpow((-2)*s_13 + IT_0039 + IT_0089 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0152 = IT_0150*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0155 = IT_0003*IT_0154;
    const ccomplex_t IT_0156 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0157 = IT_0007*IT_0156;
    const ccomplex_t IT_0158 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0159 = IT_0012*IT_0158;
    const ccomplex_t IT_0160 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0161 = IT_0015*IT_0160;
    const ccomplex_t IT_0162 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0163 = IT_0003*IT_0162;
    const ccomplex_t IT_0164 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0165 = IT_0007*IT_0164;
    const ccomplex_t IT_0166 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0167 = IT_0012*IT_0166;
    const ccomplex_t IT_0168 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0169 = IT_0015*IT_0168;
    const ccomplex_t IT_0170 = (-0.5)*IT_0155 + 0.5*IT_0157 + 0.5*IT_0159 + (
      -0.5)*IT_0161 + (-0.5)*IT_0163 + 0.5*IT_0165 + 0.5*IT_0167 + (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = m_N_3*IT_0170;
    const ccomplex_t IT_0172 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0173 = IT_0003*IT_0172;
    const ccomplex_t IT_0174 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0175 = IT_0007*IT_0174;
    const ccomplex_t IT_0176 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0177 = IT_0012*IT_0176;
    const ccomplex_t IT_0178 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0179 = IT_0015*IT_0178;
    const ccomplex_t IT_0180 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0181 = IT_0003*IT_0180;
    const ccomplex_t IT_0182 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0183 = IT_0007*IT_0182;
    const ccomplex_t IT_0184 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0185 = IT_0012*IT_0184;
    const ccomplex_t IT_0186 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0187 = IT_0015*IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0173 + (-0.5)*IT_0175 + (-0.5)*IT_0177 +
       0.5*IT_0179 + 0.5*IT_0181 + (-0.5)*IT_0183 + (-0.5)*IT_0185 + 0.5*IT_0187;
    const ccomplex_t IT_0189 = m_N_4*IT_0188;
    const ccomplex_t IT_0190 = IT_0171 + IT_0189;
    const ccomplex_t IT_0191 = IT_0153*IT_0190;
    const ccomplex_t IT_0192 = -IT_0043 + -IT_0093 + -IT_0142 + -IT_0191;
    const ccomplex_t IT_0193 = IT_0060*IT_0092;
    const ccomplex_t IT_0194 = IT_0153*IT_0170;
    const ccomplex_t IT_0195 = IT_0018*IT_0042;
    const ccomplex_t IT_0196 = IT_0110*IT_0141;
    const ccomplex_t IT_0197 = IT_0193 + IT_0194 + IT_0195 + IT_0196;
    const ccomplex_t IT_0198 = s_34*m_N_3;
    const ccomplex_t IT_0199 = s_14*IT_0039;
    const ccomplex_t IT_0200 = s_13*s_34;
    return create_ccomplex_return(2*IT_0192*(s_14*conj(IT_0192) + conj(IT_0197
      )*IT_0198) + 2*IT_0197*(conj(IT_0192)*IT_0198 + -conj(IT_0197)*(IT_0199 + 
      (-2)*IT_0200)));
}

