#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_A_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0004, -1);
    const ccomplex_t IT_0009 = IT_0002*IT_0008;
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0007 + -IT_0010 + (-2)
      *IT_0012);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = IT_0000*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0019 = IT_0005*IT_0018;
    const ccomplex_t IT_0020 = IT_0009*IT_0018;
    const ccomplex_t IT_0021 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0022 = IT_0009*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + -IT_0020 + (-2)
      *IT_0022);
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = pow(m_C_2, 2);
    const ccomplex_t IT_0030 = pow(m_Z, -2);
    const ccomplex_t IT_0031 = s_14*s_24*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = s_12 + IT_0032;
    const ccomplex_t IT_0034 = IT_0029*IT_0033;
    const ccomplex_t IT_0035 = pow(s_24, 2);
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = IT_0029 + IT_0037;
    const ccomplex_t IT_0039 = IT_0029*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = IT_0000*IT_0025;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = pow(s_14, 2);
    const ccomplex_t IT_0045 = IT_0030*IT_0044;
    const ccomplex_t IT_0046 = -IT_0029;
    const ccomplex_t IT_0047 = IT_0045 + IT_0046;
    const ccomplex_t IT_0048 = s_12*IT_0047;
    const ccomplex_t IT_0049 = (-0.5)*IT_0028;
    const ccomplex_t IT_0050 = s_12*s_14;
    const ccomplex_t IT_0051 = s_24*IT_0029;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = IT_0050 + IT_0052;
    const ccomplex_t IT_0054 = s_12*IT_0038;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = (-0.5)*IT_0017;
    const ccomplex_t IT_0057 = s_12*s_24;
    const ccomplex_t IT_0058 = s_14*IT_0029;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = IT_0057 + IT_0059;
    const ccomplex_t IT_0061 = 4*IT_0017;
    const ccomplex_t IT_0062 = m_C_2*IT_0017;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = m_C_2*IT_0042;
    const ccomplex_t IT_0065 = 0.5*IT_0064;
    const ccomplex_t IT_0066 = conj(IT_0063) + conj(IT_0065);
    const ccomplex_t IT_0067 = m_C_2*IT_0033;
    const ccomplex_t IT_0068 = IT_0061*IT_0067;
    const ccomplex_t IT_0069 = 8*IT_0028;
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = pow(m_Z, 2);
    const ccomplex_t IT_0072 = s_12*IT_0071;
    const ccomplex_t IT_0073 = s_14*s_24;
    const ccomplex_t IT_0074 = (-2)*IT_0073;
    const ccomplex_t IT_0075 = IT_0072 + IT_0074;
    const ccomplex_t IT_0076 = IT_0030*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = s_12 + IT_0077;
    const ccomplex_t IT_0079 = 8*IT_0078;
    const ccomplex_t IT_0080 = IT_0063*IT_0079;
    const ccomplex_t IT_0081 = IT_0015*IT_0026;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = m_C_2*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = m_C_2*IT_0028;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = conj(IT_0084) + conj(IT_0086);
    const ccomplex_t IT_0088 = m_C_2*IT_0038;
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = IT_0061*IT_0089;
    const ccomplex_t IT_0091 = m_C_2*IT_0047;
    const ccomplex_t IT_0092 = IT_0069*IT_0091;
    const ccomplex_t IT_0093 = (-0.5)*IT_0064;
    const ccomplex_t IT_0094 = 8*IT_0033;
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = IT_0079*IT_0084;
    const ccomplex_t IT_0097 = 12*IT_0029;
    const ccomplex_t IT_0098 = IT_0063*IT_0097;
    const ccomplex_t IT_0099 = IT_0065*IT_0097;
    const ccomplex_t IT_0100 = (-0.5)*IT_0083;
    const ccomplex_t IT_0101 = (-24)*IT_0029;
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0085;
    const ccomplex_t IT_0104 = conj(IT_0100) + conj(IT_0103);
    const ccomplex_t IT_0105 = (-8)*IT_0017;
    const ccomplex_t IT_0106 = IT_0067*IT_0105;
    const ccomplex_t IT_0107 = s_12*IT_0033;
    const ccomplex_t IT_0108 = 0.5*IT_0042;
    const ccomplex_t IT_0109 = (-4)*IT_0028;
    const ccomplex_t IT_0110 = (-0.5)*IT_0062;
    const ccomplex_t IT_0111 = conj(IT_0093) + conj(IT_0110);
    const ccomplex_t IT_0112 = IT_0079*IT_0093;
    const ccomplex_t IT_0113 = IT_0097*IT_0100;
    const ccomplex_t IT_0114 = IT_0063*IT_0101;
    const ccomplex_t IT_0115 = IT_0065*IT_0101;
    const ccomplex_t IT_0116 = IT_0089*IT_0105;
    const ccomplex_t IT_0117 = IT_0091*IT_0109;
    const ccomplex_t IT_0118 = 0.5*IT_0082;
    const ccomplex_t IT_0119 = (-4)*IT_0017;
    const ccomplex_t IT_0120 = -IT_0082;
    const ccomplex_t IT_0121 = IT_0029*IT_0047;
    const ccomplex_t IT_0122 = 4*IT_0067;
    const ccomplex_t IT_0123 = 4*IT_0060;
    const ccomplex_t IT_0124 = 4*IT_0089;
    const ccomplex_t IT_0125 = (-8)*IT_0067;
    const ccomplex_t IT_0126 = (-8)*IT_0034;
    const ccomplex_t IT_0127 = (-8)*IT_0089;
    const ccomplex_t IT_0128 = (-4)*conj(IT_0028);
    const ccomplex_t IT_0129 = (-4)*IT_0060;
    const ccomplex_t IT_0130 = conj(IT_0017)*IT_0040 + 0.5*conj(IT_0043)
      *IT_0055 + 0.125*IT_0111*IT_0122 + 0.125*conj(IT_0108)*IT_0123 + 0.125
      *IT_0104*IT_0124 + 0.125*IT_0087*IT_0125 + 0.125*conj(IT_0120)*IT_0126 +
       0.125*IT_0066*IT_0127 + 0.125*IT_0107*IT_0128 + 0.125*conj(IT_0049)
      *IT_0129;
    const ccomplex_t IT_0131 = 4*conj(IT_0017);
    const ccomplex_t IT_0132 = s_24*m_C_2;
    const ccomplex_t IT_0133 = 12*IT_0132;
    const ccomplex_t IT_0134 = s_14*m_C_2;
    const ccomplex_t IT_0135 = 12*IT_0134;
    const ccomplex_t IT_0136 = (-24)*IT_0132;
    const ccomplex_t IT_0137 = IT_0029*IT_0071;
    const ccomplex_t IT_0138 = (-24)*IT_0137;
    const ccomplex_t IT_0139 = (-12)*IT_0072;
    const ccomplex_t IT_0140 = (-4)*IT_0053;
    const ccomplex_t IT_0141 = (-4)*IT_0073;
    const ccomplex_t IT_0142 = IT_0072 + IT_0141;
    const ccomplex_t IT_0143 = (-4)*IT_0142;
    const ccomplex_t IT_0144 = 4*IT_0053;
    const ccomplex_t IT_0145 = (-3)*IT_0137;
    const ccomplex_t IT_0146 = 8*IT_0145;
    const ccomplex_t IT_0147 = s_14*IT_0071;
    const ccomplex_t IT_0148 = IT_0030*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = s_14 + IT_0149;
    const ccomplex_t IT_0151 = m_C_2*IT_0150;
    const ccomplex_t IT_0152 = 2*IT_0151;
    const ccomplex_t IT_0153 = 8*IT_0152;
    const ccomplex_t IT_0154 = (-12)*IT_0132;
    const ccomplex_t IT_0155 = (-12)*IT_0134;
    const ccomplex_t IT_0156 = (-4)*conj(IT_0017);
    const ccomplex_t IT_0157 = 8*conj(IT_0028);
    const ccomplex_t IT_0158 = (-4)*IT_0091;
    const ccomplex_t IT_0159 = 8*IT_0091;
    const ccomplex_t IT_0160 = 0.25*IT_0104;
    const ccomplex_t IT_0161 = 0.25*conj(IT_0043);
    const ccomplex_t IT_0162 = 8*IT_0067;
    const ccomplex_t IT_0163 = 0.25*IT_0111;
    const ccomplex_t IT_0164 = (-4)*IT_0067;
    const ccomplex_t IT_0165 = 0.25*IT_0087;
    const ccomplex_t IT_0166 = 0.25*conj(IT_0118);
    const ccomplex_t IT_0167 = IT_0048*conj(IT_0120) + 0.25*conj(IT_0056)
      *IT_0140 + 0.25*IT_0107*IT_0156 + 0.25*IT_0121*IT_0157 + 0.25*IT_0066
      *IT_0158 + IT_0159*IT_0160 + IT_0126*IT_0161 + IT_0162*IT_0163 + IT_0164
      *IT_0165 + IT_0144*IT_0166;
    const ccomplex_t IT_0168 = IT_0067*IT_0131;
    const ccomplex_t IT_0169 = IT_0094*conj(IT_0103);
    const ccomplex_t IT_0170 = IT_0094*conj(IT_0100);
    const ccomplex_t IT_0171 = IT_0067*IT_0157;
    const ccomplex_t IT_0172 = IT_0168 + IT_0169 + IT_0170 + IT_0171;
    const ccomplex_t IT_0173 = conj(IT_0120)*IT_0158;
    const ccomplex_t IT_0174 = conj(IT_0043)*IT_0127;
    const ccomplex_t IT_0175 = conj(IT_0056)*IT_0135;
    const ccomplex_t IT_0176 = conj(IT_0118)*IT_0155;
    const ccomplex_t IT_0177 = conj(IT_0108)*IT_0136;
    const ccomplex_t IT_0178 = IT_0089*IT_0131;
    const ccomplex_t IT_0179 = IT_0091*IT_0157;
    const ccomplex_t IT_0180 = conj(IT_0093)*IT_0094;
    const ccomplex_t IT_0181 = conj(IT_0063)*IT_0097;
    const ccomplex_t IT_0182 = conj(IT_0065)*IT_0097;
    const ccomplex_t IT_0183 = IT_0101*conj(IT_0103);
    const ccomplex_t IT_0184 = conj(IT_0100)*IT_0101;
    const ccomplex_t IT_0185 = conj(IT_0043)*IT_0125;
    const ccomplex_t IT_0186 = conj(IT_0120)*IT_0164;
    const ccomplex_t IT_0187 = IT_0094*conj(IT_0110);
    const ccomplex_t IT_0188 = conj(IT_0056)*IT_0133;
    const ccomplex_t IT_0189 = conj(IT_0118)*IT_0154;
    const ccomplex_t IT_0190 = conj(IT_0049)*IT_0153;
    const ccomplex_t IT_0191 = conj(IT_0120)*IT_0159;
    const ccomplex_t IT_0192 = (-8)*conj(IT_0017);
    const ccomplex_t IT_0193 = IT_0067*IT_0192;
    const ccomplex_t IT_0194 = IT_0067*IT_0128;
    const ccomplex_t IT_0195 = conj(IT_0043)*IT_0124;
    const ccomplex_t IT_0196 = conj(IT_0118)*IT_0153;
    const ccomplex_t IT_0197 = conj(IT_0108)*IT_0133;
    const ccomplex_t IT_0198 = conj(IT_0049)*IT_0154;
    const ccomplex_t IT_0199 = conj(IT_0043)*IT_0122;
    const ccomplex_t IT_0200 = conj(IT_0120)*IT_0162;
    const ccomplex_t IT_0201 = IT_0097*conj(IT_0103);
    const ccomplex_t IT_0202 = IT_0097*conj(IT_0100);
    const ccomplex_t IT_0203 = conj(IT_0063)*IT_0101;
    const ccomplex_t IT_0204 = conj(IT_0065)*IT_0101;
    const ccomplex_t IT_0205 = IT_0089*IT_0192;
    const ccomplex_t IT_0206 = IT_0091*IT_0128;
    const ccomplex_t IT_0207 = conj(IT_0056)*IT_0136;
    const ccomplex_t IT_0208 = conj(IT_0108)*IT_0135;
    const ccomplex_t IT_0209 = conj(IT_0049)*IT_0155;
    const ccomplex_t IT_0210 = (-8)*IT_0017*(conj(IT_0028)*IT_0034 + -IT_0040
      *conj(IT_0043)) + (-8)*IT_0028*(conj(IT_0017)*IT_0034 + -1./2*conj(IT_0028
      )*IT_0048 + -1./2*conj(IT_0049)*IT_0053) + (conj(IT_0017)*IT_0055 + conj
      (IT_0056)*IT_0060)*IT_0061 + IT_0066*(IT_0068 + IT_0070 + IT_0080) +
       IT_0087*(IT_0090 + IT_0092 + IT_0095 + IT_0096 + IT_0098 + IT_0099 +
       IT_0102) + IT_0104*IT_0106 + (IT_0067*IT_0104 + conj(IT_0043)*IT_0107 +
       IT_0053*conj(IT_0108))*IT_0109 + IT_0111*(IT_0112 + IT_0113 + IT_0114 +
       IT_0115 + IT_0116 + IT_0117) + IT_0060*conj(IT_0118)*IT_0119 + conj
      (IT_0120)*(IT_0107*IT_0119 + IT_0069*IT_0121) + 8*IT_0043*IT_0130 +
       IT_0056*(IT_0060*IT_0131 + IT_0087*IT_0133 + IT_0066*IT_0135 + IT_0111
      *IT_0136 + conj(IT_0108)*IT_0138 + conj(IT_0118)*IT_0139 + conj(IT_0120)
      *IT_0140 + conj(IT_0056)*IT_0143) + IT_0108*(conj(IT_0043)*IT_0123 +
       IT_0053*IT_0128 + IT_0104*IT_0133 + IT_0111*IT_0135 + IT_0066*IT_0136 +
       conj(IT_0056)*IT_0138 + conj(IT_0049)*IT_0139 + conj(IT_0108)*IT_0143) +
       IT_0118*(conj(IT_0056)*IT_0139 + conj(IT_0118)*IT_0143 + conj(IT_0120)
      *IT_0144 + conj(IT_0049)*IT_0146 + IT_0104*IT_0153 + IT_0087*IT_0154 +
       IT_0066*IT_0155 + IT_0060*IT_0156) + 4*IT_0049*(conj(IT_0028)*IT_0053 +
       0.25*conj(IT_0108)*IT_0139 + 0.25*conj(IT_0049)*IT_0143 + IT_0154*IT_0160
       + IT_0129*IT_0161 + IT_0155*IT_0163 + IT_0153*IT_0165 + IT_0146*IT_0166) 
      + 4*IT_0120*IT_0167 + IT_0063*(IT_0172 + IT_0173 + IT_0174 + IT_0175 +
       IT_0176 + IT_0177) + IT_0065*(IT_0066*IT_0079 + IT_0172 + IT_0173 +
       IT_0174 + IT_0175 + IT_0176 + IT_0177) + IT_0084*(IT_0178 + IT_0179 +
       IT_0180 + IT_0181 + IT_0182 + IT_0183 + IT_0184 + IT_0185 + IT_0186 +
       IT_0187 + IT_0188 + IT_0189 + IT_0190) + IT_0086*(IT_0079*IT_0087 +
       IT_0178 + IT_0179 + IT_0180 + IT_0181 + IT_0182 + IT_0183 + IT_0184 +
       IT_0185 + IT_0186 + IT_0187 + IT_0188 + IT_0189 + IT_0190) + IT_0100*
      (IT_0066*IT_0094 + IT_0079*IT_0104 + IT_0191 + IT_0193 + IT_0194 + IT_0195
       + IT_0196 + IT_0197 + IT_0198) + IT_0103*(IT_0066*IT_0094 + IT_0087
      *IT_0101 + IT_0079*IT_0104 + IT_0097*IT_0111 + IT_0191 + IT_0193 + IT_0194
       + IT_0195 + IT_0196 + IT_0197 + IT_0198) + IT_0093*(IT_0199 + IT_0200 +
       IT_0201 + IT_0202 + IT_0203 + IT_0204 + IT_0205 + IT_0206 + IT_0207 +
       IT_0208 + IT_0209) + IT_0110*(IT_0087*IT_0094 + IT_0079*IT_0111 + IT_0199
       + IT_0200 + IT_0201 + IT_0202 + IT_0203 + IT_0204 + IT_0205 + IT_0206 +
       IT_0207 + IT_0208 + IT_0209);
    return create_ccomplex_return(IT_0210);
}

