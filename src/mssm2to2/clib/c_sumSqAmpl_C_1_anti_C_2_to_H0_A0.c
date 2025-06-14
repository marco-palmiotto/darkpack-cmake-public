#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_H0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_H0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_H0, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = U_d1*V_Wp1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = V_u1*U_Wm1;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (IT_0005 + IT_0008)*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = V_u1*U_Wm2;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = sin(beta);
    const ccomplex_t IT_0019 = U_d2*V_Wp1;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = IT_0010*(IT_0017 + IT_0020);
    const ccomplex_t IT_0022 = 1.4142135623731*e_em*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = IT_0014*IT_0023;
    const ccomplex_t IT_0025 = IT_0001*IT_0002*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0028 = U_d2*V_Wp2;
    const ccomplex_t IT_0029 = IT_0003*IT_0028;
    const ccomplex_t IT_0030 = V_u2*U_Wm2;
    const ccomplex_t IT_0031 = IT_0006*IT_0030;
    const ccomplex_t IT_0032 = IT_0010*(IT_0029 + IT_0031);
    const ccomplex_t IT_0033 = 1.4142135623731*e_em*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0023*IT_0035;
    const ccomplex_t IT_0037 = IT_0026*IT_0027*IT_0036;
    const ccomplex_t IT_0038 = IT_0007*IT_0015;
    const ccomplex_t IT_0039 = IT_0004*IT_0018;
    const ccomplex_t IT_0040 = IT_0010*(IT_0038 + IT_0039);
    const ccomplex_t IT_0041 = 1.4142135623731*e_em*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0003*IT_0019;
    const ccomplex_t IT_0044 = IT_0006*IT_0016;
    const ccomplex_t IT_0045 = IT_0010*(IT_0043 + IT_0044);
    const ccomplex_t IT_0046 = 1.4142135623731*e_em*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0042*IT_0048;
    const ccomplex_t IT_0050 = pow(m_C_2, 2);
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0000 + IT_0050 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0002*IT_0049*IT_0051;
    const ccomplex_t IT_0053 = IT_0015*IT_0030;
    const ccomplex_t IT_0054 = IT_0018*IT_0028;
    const ccomplex_t IT_0055 = IT_0010*(IT_0053 + IT_0054);
    const ccomplex_t IT_0056 = 1.4142135623731*e_em*IT_0055;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = IT_0048*IT_0057;
    const ccomplex_t IT_0059 = pow(m_C_1, 2);
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0000 + IT_0059 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0027*IT_0058*IT_0060;
    const ccomplex_t IT_0062 = cos(theta_W);
    const ccomplex_t IT_0063 = cpow(IT_0062, -1);
    const ccomplex_t IT_0064 = -e_em*IT_0010*(IT_0006*IT_0015 + -IT_0003
      *IT_0018)*IT_0063;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = m_C_1*IT_0065;
    const ccomplex_t IT_0067 = IT_0009*IT_0063;
    const ccomplex_t IT_0068 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = IT_0010*IT_0062;
    const ccomplex_t IT_0071 = IT_0068*IT_0070;
    const ccomplex_t IT_0072 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0069 + -IT_0071 + (-2)
      *IT_0073);
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0050 + -IT_0059 + -reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0066*IT_0078;
    const ccomplex_t IT_0080 = m_C_2*IT_0065;
    const ccomplex_t IT_0081 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0082 = IT_0067*IT_0081;
    const ccomplex_t IT_0083 = IT_0070*IT_0081;
    const ccomplex_t IT_0084 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0085 = IT_0070*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + -IT_0083 + (-2)
      *IT_0085);
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0076*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0080*IT_0089;
    const ccomplex_t IT_0091 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0092 = IT_0003*IT_0091;
    const ccomplex_t IT_0093 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0094 = IT_0006*IT_0093;
    const ccomplex_t IT_0095 = IT_0010*(IT_0092 + IT_0094);
    const ccomplex_t IT_0096 = 1.4142135623731*e_em*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = IT_0023*IT_0098;
    const ccomplex_t IT_0100 = IT_0026*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = m_C_2*IT_0101;
    const ccomplex_t IT_0103 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0104 = IT_0003*IT_0103;
    const ccomplex_t IT_0105 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0106 = IT_0006*IT_0105;
    const ccomplex_t IT_0107 = IT_0010*(IT_0104 + IT_0106);
    const ccomplex_t IT_0108 = 1.4142135623731*e_em*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0057*IT_0110;
    const ccomplex_t IT_0112 = IT_0060*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = m_C_1*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_12 + -IT_0050 + -IT_0059 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0116 = cpow(IT_0062, -2);
    const ccomplex_t IT_0117 = IT_0009*IT_0116;
    const ccomplex_t IT_0118 = IT_0010 + IT_0117;
    const ccomplex_t IT_0119 = cpow(IT_0018, 3);
    const ccomplex_t IT_0120 = -IT_0015;
    const ccomplex_t IT_0121 = cpow(IT_0018, 2);
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*m_W*e_em*IT_0118*(IT_0006
      *IT_0119 + IT_0120*((IT_0003*IT_0015 + -IT_0006*IT_0018)*IT_0120 + IT_0003
      *IT_0121));
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = IT_0023*IT_0123;
    const ccomplex_t IT_0125 = IT_0115*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = pow(m_Z, 2);
    const ccomplex_t IT_0128 = cpow((-2)*s_12 + -IT_0050 + -IT_0059 + IT_0127 
      + -reg_prop, -1);
    const ccomplex_t IT_0129 = IT_0015*IT_0019;
    const ccomplex_t IT_0130 = IT_0016*IT_0018;
    const ccomplex_t IT_0131 = IT_0010*(IT_0129 + -IT_0130);
    const ccomplex_t IT_0132 = 1.4142135623731*e_em*IT_0131;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*m_W*e_em*IT_0015*IT_0018*
      (IT_0003*IT_0015 + -IT_0006*IT_0018)*IT_0118;
    const ccomplex_t IT_0135 = IT_0133*IT_0134;
    const ccomplex_t IT_0136 = IT_0128*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0139 = IT_0003*IT_0138;
    const ccomplex_t IT_0140 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0141 = IT_0006*IT_0140;
    const ccomplex_t IT_0142 = IT_0010*(IT_0139 + IT_0141);
    const ccomplex_t IT_0143 = 1.4142135623731*e_em*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = 0.5*IT_0144;
    const ccomplex_t IT_0146 = IT_0023*IT_0145;
    const ccomplex_t IT_0147 = IT_0001*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = m_C_1*IT_0148;
    const ccomplex_t IT_0150 = IT_0042*IT_0110;
    const ccomplex_t IT_0151 = IT_0051*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = m_C_2*IT_0152;
    const ccomplex_t IT_0154 = IT_0025 + IT_0037 + IT_0052 + IT_0061 + 
      -IT_0079 + -IT_0090 + IT_0102 + IT_0114 + IT_0126 + IT_0137 + IT_0149 +
       IT_0153;
    const ccomplex_t IT_0155 = IT_0064*IT_0078;
    const ccomplex_t IT_0156 = IT_0015*IT_0140;
    const ccomplex_t IT_0157 = IT_0018*IT_0138;
    const ccomplex_t IT_0158 = IT_0010*(IT_0156 + IT_0157);
    const ccomplex_t IT_0159 = 1.4142135623731*e_em*IT_0158;
    const ccomplex_t IT_0160 = (-0.5)*IT_0159;
    const ccomplex_t IT_0161 = IT_0048*IT_0160;
    const ccomplex_t IT_0162 = IT_0051*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = IT_0015*IT_0105;
    const ccomplex_t IT_0165 = IT_0018*IT_0103;
    const ccomplex_t IT_0166 = IT_0010*(IT_0164 + IT_0165);
    const ccomplex_t IT_0167 = 1.4142135623731*e_em*IT_0166;
    const ccomplex_t IT_0168 = (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = IT_0035*IT_0168;
    const ccomplex_t IT_0170 = IT_0026*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = -IT_0113 + -IT_0148 + IT_0155 + IT_0163 +
       IT_0171;
    const ccomplex_t IT_0173 = s_23*m_C_1;
    const ccomplex_t IT_0174 = 2*IT_0173;
    const ccomplex_t IT_0175 = IT_0110*IT_0160;
    const ccomplex_t IT_0176 = IT_0002*IT_0051*IT_0175;
    const ccomplex_t IT_0177 = IT_0123*IT_0168;
    const ccomplex_t IT_0178 = IT_0115*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = IT_0015*IT_0103;
    const ccomplex_t IT_0181 = IT_0018*IT_0105;
    const ccomplex_t IT_0182 = IT_0010*(IT_0180 + -IT_0181);
    const ccomplex_t IT_0183 = 1.4142135623731*e_em*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = IT_0134*IT_0184;
    const ccomplex_t IT_0186 = IT_0128*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = IT_0078*IT_0080;
    const ccomplex_t IT_0189 = IT_0066*IT_0089;
    const ccomplex_t IT_0190 = IT_0145*IT_0168;
    const ccomplex_t IT_0191 = IT_0001*IT_0002*IT_0190;
    const ccomplex_t IT_0192 = IT_0015*IT_0093;
    const ccomplex_t IT_0193 = IT_0018*IT_0091;
    const ccomplex_t IT_0194 = IT_0010*(IT_0192 + IT_0193);
    const ccomplex_t IT_0195 = 1.4142135623731*e_em*IT_0194;
    const ccomplex_t IT_0196 = (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = IT_0110*IT_0196;
    const ccomplex_t IT_0198 = IT_0027*IT_0060*IT_0197;
    const ccomplex_t IT_0199 = m_C_2*IT_0171;
    const ccomplex_t IT_0200 = IT_0098*IT_0168;
    const ccomplex_t IT_0201 = IT_0026*IT_0027*IT_0200;
    const ccomplex_t IT_0202 = m_C_2*IT_0163;
    const ccomplex_t IT_0203 = IT_0014*IT_0168;
    const ccomplex_t IT_0204 = IT_0001*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = m_C_1*IT_0205;
    const ccomplex_t IT_0207 = IT_0048*IT_0196;
    const ccomplex_t IT_0208 = IT_0060*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = m_C_1*IT_0209;
    const ccomplex_t IT_0211 = IT_0176 + IT_0179 + IT_0187 + IT_0188 + IT_0189
       + IT_0191 + IT_0198 + IT_0199 + IT_0201 + IT_0202 + IT_0206 + IT_0210;
    const ccomplex_t IT_0212 = m_C_1*m_C_2;
    const ccomplex_t IT_0213 = (-2)*IT_0212;
    const ccomplex_t IT_0214 = IT_0064*IT_0089;
    const ccomplex_t IT_0215 = IT_0101 + IT_0152 + -IT_0205 + -IT_0209 + 
      -IT_0214;
    const ccomplex_t IT_0216 = s_13*m_C_2;
    const ccomplex_t IT_0217 = (-2)*IT_0216;
    const ccomplex_t IT_0218 = s_12*IT_0000;
    const ccomplex_t IT_0219 = s_13*s_23;
    const ccomplex_t IT_0220 = (-2)*IT_0219;
    const ccomplex_t IT_0221 = IT_0218 + IT_0220;
    const ccomplex_t IT_0222 = (-2)*IT_0221;
    const ccomplex_t IT_0223 = m_C_1*m_C_2*IT_0000;
    const ccomplex_t IT_0224 = (-2)*IT_0223;
    const ccomplex_t IT_0225 = 2*IT_0211*(s_12*conj(IT_0211) + 0.5*conj
      (IT_0154)*IT_0213 + 0.5*IT_0174*conj(IT_0215) + 0.5*conj(IT_0172)*IT_0217)
       + 2*IT_0154*(s_12*conj(IT_0154) + 0.5*conj(IT_0172)*IT_0174 + 0.5*conj
      (IT_0211)*IT_0213 + 0.5*conj(IT_0215)*IT_0217) + IT_0215*(IT_0174*conj
      (IT_0211) + conj(IT_0154)*IT_0217 + conj(IT_0215)*IT_0222 + conj(IT_0172)
      *IT_0224) + IT_0172*(conj(IT_0154)*IT_0174 + conj(IT_0211)*IT_0217 + conj
      (IT_0172)*IT_0222 + conj(IT_0215)*IT_0224);
    return create_ccomplex_return(IT_0225);
}

