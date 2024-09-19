#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_H0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_H0_A0(
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + -IT_0004
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = IT_0001*IT_0006;
    const ccomplex_t IT_0010 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0000*IT_0007;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0011 + -IT_0013 + (-2)
      *IT_0015);
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0008*IT_0021;
    const ccomplex_t IT_0023 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0024 = IT_0002*IT_0023;
    const ccomplex_t IT_0025 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0026 = IT_0005*IT_0025;
    const ccomplex_t IT_0027 = IT_0007*(IT_0024 + IT_0026);
    const ccomplex_t IT_0028 = 1.4142135623731*e_em*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = V_u2*U_Wm1;
    const ccomplex_t IT_0032 = IT_0004*IT_0031;
    const ccomplex_t IT_0033 = U_d1*V_Wp2;
    const ccomplex_t IT_0034 = IT_0003*IT_0033;
    const ccomplex_t IT_0035 = IT_0007*(IT_0032 + IT_0034);
    const ccomplex_t IT_0036 = 1.4142135623731*e_em*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0030*IT_0037;
    const ccomplex_t IT_0039 = pow(m_H0, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0018 + IT_0039 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0039 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0044 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0045 = IT_0004*IT_0044;
    const ccomplex_t IT_0046 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0047 = IT_0003*IT_0046;
    const ccomplex_t IT_0048 = IT_0007*(IT_0045 + IT_0047);
    const ccomplex_t IT_0049 = 1.4142135623731*e_em*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = U_d1*V_Wp1;
    const ccomplex_t IT_0052 = IT_0002*IT_0051;
    const ccomplex_t IT_0053 = V_u1*U_Wm1;
    const ccomplex_t IT_0054 = IT_0005*IT_0053;
    const ccomplex_t IT_0055 = IT_0007*(IT_0052 + IT_0054);
    const ccomplex_t IT_0056 = 1.4142135623731*e_em*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0058;
    const ccomplex_t IT_0060 = IT_0043*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0004*IT_0025;
    const ccomplex_t IT_0063 = IT_0003*IT_0023;
    const ccomplex_t IT_0064 = IT_0007*(IT_0062 + IT_0063);
    const ccomplex_t IT_0065 = 1.4142135623731*e_em*IT_0064;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0002*IT_0033;
    const ccomplex_t IT_0068 = IT_0005*IT_0031;
    const ccomplex_t IT_0069 = IT_0007*(IT_0067 + IT_0068);
    const ccomplex_t IT_0070 = 1.4142135623731*e_em*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0066*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_23 + IT_0018 + IT_0039 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_13 + IT_0039 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0002*IT_0046;
    const ccomplex_t IT_0079 = IT_0005*IT_0044;
    const ccomplex_t IT_0080 = IT_0007*(IT_0078 + IT_0079);
    const ccomplex_t IT_0081 = 1.4142135623731*e_em*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0004*IT_0053;
    const ccomplex_t IT_0085 = IT_0003*IT_0051;
    const ccomplex_t IT_0086 = IT_0007*(IT_0084 + IT_0085);
    const ccomplex_t IT_0087 = 1.4142135623731*e_em*IT_0086;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0083*IT_0088;
    const ccomplex_t IT_0090 = IT_0077*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0022 + -IT_0042 + IT_0061 + IT_0076 + 
      -IT_0091;
    const ccomplex_t IT_0093 = 0.5*IT_0008;
    const ccomplex_t IT_0094 = m_C_1*IT_0093;
    const ccomplex_t IT_0095 = IT_0021*IT_0094;
    const ccomplex_t IT_0096 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0097 = IT_0009*IT_0096;
    const ccomplex_t IT_0098 = IT_0012*IT_0096;
    const ccomplex_t IT_0099 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0100 = IT_0012*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0097 + -IT_0098 + (-2)
      *IT_0100);
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = IT_0019*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0094*IT_0104;
    const ccomplex_t IT_0106 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0107 = IT_0002*IT_0106;
    const ccomplex_t IT_0108 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0109 = IT_0005*IT_0108;
    const ccomplex_t IT_0110 = IT_0007*(IT_0107 + IT_0109);
    const ccomplex_t IT_0111 = 1.4142135623731*e_em*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0066*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0116 = IT_0074*IT_0114*IT_0115;
    const ccomplex_t IT_0117 = IT_0004*IT_0108;
    const ccomplex_t IT_0118 = IT_0003*IT_0106;
    const ccomplex_t IT_0119 = IT_0007*(IT_0117 + IT_0118);
    const ccomplex_t IT_0120 = 1.4142135623731*e_em*IT_0119;
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = IT_0030*IT_0121;
    const ccomplex_t IT_0123 = IT_0040*IT_0115*IT_0122;
    const ccomplex_t IT_0124 = cpow((-2)*s_12 + (-2)*IT_0018 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0125 = cpow(IT_0000, -2);
    const ccomplex_t IT_0126 = IT_0006*IT_0125;
    const ccomplex_t IT_0127 = IT_0007 + IT_0126;
    const ccomplex_t IT_0128 = cpow(IT_0003, 3);
    const ccomplex_t IT_0129 = -IT_0004;
    const ccomplex_t IT_0130 = cpow(IT_0003, 2);
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*m_W*e_em*IT_0127*(IT_0005
      *IT_0128 + IT_0129*((IT_0002*IT_0004 + -IT_0003*IT_0005)*IT_0129 + IT_0002
      *IT_0130));
    const ccomplex_t IT_0132 = 0.5*IT_0131;
    const ccomplex_t IT_0133 = IT_0050*IT_0132;
    const ccomplex_t IT_0134 = IT_0124*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0050*IT_0083;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0138 = IT_0043*IT_0136*IT_0137;
    const ccomplex_t IT_0139 = IT_0077*IT_0136*IT_0137;
    const ccomplex_t IT_0140 = IT_0059*IT_0077;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = m_C_1*IT_0141;
    const ccomplex_t IT_0143 = m_C_1*IT_0061;
    const ccomplex_t IT_0144 = U_d2*V_Wp1;
    const ccomplex_t IT_0145 = IT_0002*IT_0144;
    const ccomplex_t IT_0146 = V_u1*U_Wm2;
    const ccomplex_t IT_0147 = IT_0005*IT_0146;
    const ccomplex_t IT_0148 = IT_0007*(IT_0145 + IT_0147);
    const ccomplex_t IT_0149 = 1.4142135623731*e_em*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0121*IT_0151;
    const ccomplex_t IT_0153 = IT_0040*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = m_C_1*IT_0154;
    const ccomplex_t IT_0156 = m_C_1*IT_0076;
    const ccomplex_t IT_0157 = IT_0004*IT_0046;
    const ccomplex_t IT_0158 = IT_0003*IT_0044;
    const ccomplex_t IT_0159 = IT_0007*(IT_0157 + -IT_0158);
    const ccomplex_t IT_0160 = 1.4142135623731*e_em*IT_0159;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*m_W*e_em*IT_0003*IT_0004*
      (IT_0002*IT_0004 + -IT_0003*IT_0005)*IT_0127;
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = pow(m_Z, 2);
    const ccomplex_t IT_0165 = cpow((-2)*s_12 + (-2)*IT_0018 + IT_0164 + 
      -reg_prop, -1);
    const ccomplex_t IT_0166 = IT_0163*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = IT_0095 + IT_0105 + IT_0116 + IT_0123 + IT_0135
       + IT_0138 + IT_0139 + IT_0142 + IT_0143 + IT_0155 + IT_0156 + IT_0167;
    const ccomplex_t IT_0169 = s_13*m_C_1;
    const ccomplex_t IT_0170 = IT_0058*IT_0088;
    const ccomplex_t IT_0171 = IT_0077*IT_0137*IT_0170;
    const ccomplex_t IT_0172 = IT_0043*IT_0137*IT_0170;
    const ccomplex_t IT_0173 = IT_0004*IT_0146;
    const ccomplex_t IT_0174 = IT_0003*IT_0144;
    const ccomplex_t IT_0175 = IT_0007*(IT_0173 + IT_0174);
    const ccomplex_t IT_0176 = 1.4142135623731*e_em*IT_0175;
    const ccomplex_t IT_0177 = 0.5*IT_0176;
    const ccomplex_t IT_0178 = IT_0072*IT_0177;
    const ccomplex_t IT_0179 = IT_0074*IT_0115*IT_0178;
    const ccomplex_t IT_0180 = IT_0037*IT_0151;
    const ccomplex_t IT_0181 = IT_0040*IT_0115*IT_0180;
    const ccomplex_t IT_0182 = IT_0113*IT_0177;
    const ccomplex_t IT_0183 = IT_0074*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = m_C_1*IT_0184;
    const ccomplex_t IT_0186 = m_C_1*IT_0042;
    const ccomplex_t IT_0187 = IT_0088*IT_0132;
    const ccomplex_t IT_0188 = IT_0124*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = IT_0004*IT_0051;
    const ccomplex_t IT_0191 = IT_0003*IT_0053;
    const ccomplex_t IT_0192 = IT_0007*(IT_0190 + -IT_0191);
    const ccomplex_t IT_0193 = 1.4142135623731*e_em*IT_0192;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = IT_0162*IT_0194;
    const ccomplex_t IT_0196 = IT_0165*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = IT_0043*IT_0089;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = m_C_1*IT_0199;
    const ccomplex_t IT_0201 = m_C_1*IT_0091;
    const ccomplex_t IT_0202 = -IT_0095 + -IT_0105 + IT_0171 + IT_0172 +
       IT_0179 + IT_0181 + IT_0185 + IT_0186 + IT_0189 + IT_0197 + IT_0200 +
       IT_0201;
    const ccomplex_t IT_0203 = s_23*m_C_1;
    const ccomplex_t IT_0204 = IT_0008*IT_0104;
    const ccomplex_t IT_0205 = -IT_0141 + -IT_0154 + IT_0184 + IT_0199 + 
      -IT_0204;
    const ccomplex_t IT_0206 = IT_0018*IT_0039;
    const ccomplex_t IT_0207 = s_12*IT_0039;
    const ccomplex_t IT_0208 = s_13*s_23;
    const ccomplex_t IT_0209 = (-2)*IT_0208;
    const ccomplex_t IT_0210 = IT_0207 + IT_0209;
    const ccomplex_t IT_0211 = 2*s_12;
    const ccomplex_t IT_0212 = (-2)*conj(IT_0202);
    const ccomplex_t IT_0213 = -conj(IT_0092);
    const ccomplex_t IT_0214 = (-2)*IT_0092*(conj(IT_0168)*IT_0169 + -conj
      (IT_0202)*IT_0203 + conj(IT_0205)*IT_0206 + conj(IT_0092)*IT_0210) + (-2)
      *IT_0168*(conj(IT_0092)*IT_0169 + -IT_0203*conj(IT_0205) + (-0.5)*conj
      (IT_0168)*IT_0211 + (-0.5)*IT_0018*IT_0212) + (-2)*IT_0202*(IT_0018*conj
      (IT_0168) + IT_0169*conj(IT_0205) + (-0.5)*conj(IT_0202)*IT_0211 + IT_0203
      *IT_0213) + 2*IT_0205*(conj(IT_0168)*IT_0203 + -conj(IT_0205)*IT_0210 +
       0.5*IT_0169*IT_0212 + IT_0206*IT_0213);
    return create_ccomplex_return(IT_0214);
}

