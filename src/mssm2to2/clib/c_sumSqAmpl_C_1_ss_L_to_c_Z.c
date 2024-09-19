#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_ss_L_to_c_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_ss_L_to_c_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = IT_0000*IT_0004;
    const ccomplex_t IT_0006 = e_em*IT_0005;
    const ccomplex_t IT_0007 = IT_0001*IT_0003;
    const ccomplex_t IT_0008 = e_em*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0006 + (-3)*IT_0008);
    const ccomplex_t IT_0010 = (-0.166666666666667)*IT_0009;
    const ccomplex_t IT_0011 = IT_0002*IT_0010;
    const ccomplex_t IT_0012 = pow(m_c, 2);
    const ccomplex_t IT_0013 = pow(m_C_1, 2);
    const ccomplex_t IT_0014 = pow(m_ss_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + IT_0012 + -IT_0013 + -IT_0014 
      + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0011*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0019 = IT_0005*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0018;
    const ccomplex_t IT_0021 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0022 = IT_0007*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + -IT_0020 + (-2)
      *IT_0022);
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0002*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0012 + IT_0014 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*3)*e_em*(IT_0001*IT_0003 + 1./3
      *IT_0000*IT_0004);
    const ccomplex_t IT_0030 = 0.333333333333333*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0002*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*V_cs*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0036 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0037 = IT_0005*IT_0036;
    const ccomplex_t IT_0038 = IT_0007*IT_0036;
    const ccomplex_t IT_0039 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0040 = IT_0007*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0038 + (-2)
      *IT_0040);
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0035*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0012 + IT_0014 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = 2*IT_0028 + -IT_0034 + 2*IT_0046;
    const ccomplex_t IT_0048 = s_13*s_14;
    const ccomplex_t IT_0049 = s_34*IT_0013;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = 2*IT_0017 + IT_0034;
    const ccomplex_t IT_0053 = s_14*IT_0012;
    const ccomplex_t IT_0054 = s_13*s_34;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0053 + IT_0055;
    const ccomplex_t IT_0057 = pow(m_Z, -2);
    const ccomplex_t IT_0058 = pow(m_Z, 4);
    const ccomplex_t IT_0059 = s_13*IT_0058;
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = 6*IT_0060;
    const ccomplex_t IT_0062 = s_14*s_34;
    const ccomplex_t IT_0063 = 24*IT_0062;
    const ccomplex_t IT_0064 = pow(m_Z, 2);
    const ccomplex_t IT_0065 = s_13*IT_0064;
    const ccomplex_t IT_0066 = (-12)*IT_0065;
    const ccomplex_t IT_0067 = IT_0061 + IT_0063 + IT_0066;
    const ccomplex_t IT_0068 = -IT_0028 + -IT_0046;
    const ccomplex_t IT_0069 = (-24)*IT_0062;
    const ccomplex_t IT_0070 = IT_0061 + IT_0069;
    const ccomplex_t IT_0071 = pow(m_W, -1);
    const ccomplex_t IT_0072 = sin(beta);
    const ccomplex_t IT_0073 = cpow(IT_0072, -1);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u2
      *e_em*IT_0001*IT_0071*IT_0073;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0077 = IT_0005*IT_0076;
    const ccomplex_t IT_0078 = IT_0007*IT_0076;
    const ccomplex_t IT_0079 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0080 = IT_0007*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + -IT_0078 + (-2)
      *IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0075*IT_0082;
    const ccomplex_t IT_0084 = IT_0044*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0001*IT_0071*IT_0073;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0090 = IT_0005*IT_0089;
    const ccomplex_t IT_0091 = IT_0007*IT_0089;
    const ccomplex_t IT_0092 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0093 = IT_0007*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + -IT_0091 + (-2)
      *IT_0093);
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0088*IT_0095;
    const ccomplex_t IT_0097 = IT_0026*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = m_c*m_C_1;
    const ccomplex_t IT_0101 = IT_0064*IT_0100;
    const ccomplex_t IT_0102 = 3*IT_0101;
    const ccomplex_t IT_0103 = 6*IT_0017;
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = IT_0067*IT_0086;
    const ccomplex_t IT_0106 = IT_0104 + IT_0105;
    const ccomplex_t IT_0107 = (-3)*IT_0101;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*e_em*IT_0000*IT_0004;
    const ccomplex_t IT_0109 = (-0.666666666666667)*IT_0108;
    const ccomplex_t IT_0110 = IT_0088*IT_0109;
    const ccomplex_t IT_0111 = IT_0015*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = s_34*IT_0064;
    const ccomplex_t IT_0115 = IT_0057*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = s_34 + IT_0116;
    const ccomplex_t IT_0118 = m_C_1*IT_0117;
    const ccomplex_t IT_0119 = 2*IT_0118;
    const ccomplex_t IT_0120 = IT_0010*IT_0088;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0122 = IT_0015*IT_0120*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0124 = IT_0002*IT_0095;
    const ccomplex_t IT_0125 = IT_0026*IT_0123*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0127 = IT_0035*IT_0082;
    const ccomplex_t IT_0128 = IT_0044*IT_0126*IT_0127;
    const ccomplex_t IT_0129 = m_C_1*IT_0028;
    const ccomplex_t IT_0130 = m_c*IT_0112;
    const ccomplex_t IT_0131 = m_C_1*IT_0046;
    const ccomplex_t IT_0132 = -IT_0122 + -IT_0125 + -IT_0128 + -IT_0129 +
       IT_0130 + -IT_0131;
    const ccomplex_t IT_0133 = s_14*IT_0064;
    const ccomplex_t IT_0134 = IT_0057*IT_0133;
    const ccomplex_t IT_0135 = (-0.25)*IT_0134;
    const ccomplex_t IT_0136 = s_14 + IT_0135;
    const ccomplex_t IT_0137 = m_c*IT_0136;
    const ccomplex_t IT_0138 = (-4)*IT_0137;
    const ccomplex_t IT_0139 = IT_0024*IT_0088;
    const ccomplex_t IT_0140 = IT_0026*IT_0123*IT_0139;
    const ccomplex_t IT_0141 = m_c*IT_0017;
    const ccomplex_t IT_0142 = IT_0002*IT_0109;
    const ccomplex_t IT_0143 = IT_0015*IT_0121*IT_0142;
    const ccomplex_t IT_0144 = IT_0042*IT_0075;
    const ccomplex_t IT_0145 = IT_0044*IT_0126*IT_0144;
    const ccomplex_t IT_0146 = m_C_1*IT_0098;
    const ccomplex_t IT_0147 = m_C_1*IT_0085;
    const ccomplex_t IT_0148 = -IT_0140 + -IT_0141 + IT_0143 + -IT_0145 + 
      -IT_0146 + -IT_0147;
    const ccomplex_t IT_0149 = 6*IT_0051;
    const ccomplex_t IT_0150 = 6*IT_0056;
    const ccomplex_t IT_0151 = conj(IT_0086) + conj(IT_0099);
    const ccomplex_t IT_0152 = 6*IT_0107;
    const ccomplex_t IT_0153 = 6*IT_0102;
    const ccomplex_t IT_0154 = 0.5*IT_0134;
    const ccomplex_t IT_0155 = s_14 + IT_0154;
    const ccomplex_t IT_0156 = m_c*IT_0155;
    const ccomplex_t IT_0157 = 2*IT_0156;
    const ccomplex_t IT_0158 = 6*IT_0157;
    const ccomplex_t IT_0159 = (-0.25)*IT_0115;
    const ccomplex_t IT_0160 = s_34 + IT_0159;
    const ccomplex_t IT_0161 = m_C_1*IT_0160;
    const ccomplex_t IT_0162 = (-4)*IT_0161;
    const ccomplex_t IT_0163 = 6*IT_0162;
    const ccomplex_t IT_0164 = 6*conj(IT_0017);
    const ccomplex_t IT_0165 = 6*IT_0119;
    const ccomplex_t IT_0166 = 6*IT_0138;
    const ccomplex_t IT_0167 = IT_0031*IT_0088;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0030*IT_0168;
    const ccomplex_t IT_0170 = 2*IT_0085 + 2*IT_0098 + IT_0169;
    const ccomplex_t IT_0171 = (-6)*IT_0051;
    const ccomplex_t IT_0172 = (-2)*IT_0112 + -IT_0169;
    const ccomplex_t IT_0173 = (-6)*IT_0056;
    const ccomplex_t IT_0174 = pow(s_14, 2);
    const ccomplex_t IT_0175 = IT_0057*IT_0174;
    const ccomplex_t IT_0176 = -IT_0013;
    const ccomplex_t IT_0177 = IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = s_13*IT_0177;
    const ccomplex_t IT_0179 = 6*IT_0178;
    const ccomplex_t IT_0180 = m_c*IT_0177;
    const ccomplex_t IT_0181 = 6*IT_0180;
    const ccomplex_t IT_0182 = IT_0100*IT_0177;
    const ccomplex_t IT_0183 = 6*IT_0182;
    const ccomplex_t IT_0184 = (-6)*conj(IT_0017);
    const ccomplex_t IT_0185 = s_14*s_34*IT_0057;
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = s_13 + IT_0186;
    const ccomplex_t IT_0188 = s_13*IT_0187;
    const ccomplex_t IT_0189 = (-6)*IT_0188;
    const ccomplex_t IT_0190 = IT_0100*IT_0187;
    const ccomplex_t IT_0191 = (-6)*IT_0190;
    const ccomplex_t IT_0192 = m_C_1*IT_0187;
    const ccomplex_t IT_0193 = (-6)*IT_0192;
    const ccomplex_t IT_0194 = 12*s_13;
    const ccomplex_t IT_0195 = (-2)*IT_0062;
    const ccomplex_t IT_0196 = IT_0065 + IT_0195;
    const ccomplex_t IT_0197 = IT_0057*IT_0196;
    const ccomplex_t IT_0198 = (-6)*IT_0197;
    const ccomplex_t IT_0199 = IT_0194 + IT_0198;
    const ccomplex_t IT_0200 = (-18)*IT_0100;
    const ccomplex_t IT_0201 = m_c*IT_0187;
    const ccomplex_t IT_0202 = (-6)*IT_0201;
    const ccomplex_t IT_0203 = pow(s_34, 2);
    const ccomplex_t IT_0204 = IT_0057*IT_0203;
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0012 + IT_0205;
    const ccomplex_t IT_0207 = m_C_1*IT_0206;
    const ccomplex_t IT_0208 = (-6)*IT_0207;
    const ccomplex_t IT_0209 = IT_0148*(conj(IT_0068)*IT_0158 + IT_0151
      *IT_0163 + IT_0138*IT_0164 + conj(IT_0113)*IT_0165 + conj(IT_0047)*IT_0181
       + conj(IT_0170)*IT_0193 + conj(IT_0148)*IT_0199 + conj(IT_0132)*IT_0200 +
       conj(IT_0052)*IT_0202 + conj(IT_0172)*IT_0208);
    const ccomplex_t IT_0210 = IT_0132*(IT_0151*IT_0158 + conj(IT_0068)
      *IT_0163 + IT_0119*IT_0164 + conj(IT_0113)*IT_0166 + conj(IT_0170)*IT_0181
       + conj(IT_0047)*IT_0193 + conj(IT_0132)*IT_0199 + conj(IT_0148)*IT_0200 +
       conj(IT_0172)*IT_0202 + conj(IT_0052)*IT_0208);
    const ccomplex_t IT_0211 = s_13*IT_0206;
    const ccomplex_t IT_0212 = (-6)*IT_0211;
    const ccomplex_t IT_0213 = IT_0100*IT_0206;
    const ccomplex_t IT_0214 = (-6)*IT_0213;
    const ccomplex_t IT_0215 = IT_0102*IT_0164;
    const ccomplex_t IT_0216 = IT_0149*conj(IT_0170);
    const ccomplex_t IT_0217 = IT_0150*conj(IT_0172);
    const ccomplex_t IT_0218 = conj(IT_0132)*IT_0158;
    const ccomplex_t IT_0219 = conj(IT_0148)*IT_0163;
    const ccomplex_t IT_0220 = conj(IT_0068)*IT_0152;
    const ccomplex_t IT_0221 = (-6)*IT_0017*(conj(IT_0047)*IT_0051 + conj
      (IT_0052)*IT_0056 + (-0.166666666666667)*conj(IT_0017)*IT_0067 + (
      -0.166666666666667)*conj(IT_0068)*IT_0070) + (conj(IT_0086) + conj(IT_0099
      ))*IT_0106 + IT_0103*(IT_0107*conj(IT_0113) + IT_0119*conj(IT_0132) +
       IT_0138*conj(IT_0148)) + IT_0068*(IT_0067*conj(IT_0068) + conj(IT_0017)
      *IT_0070 + conj(IT_0047)*IT_0149 + conj(IT_0052)*IT_0150 + IT_0151*IT_0152
       + conj(IT_0113)*IT_0153 + conj(IT_0148)*IT_0158 + conj(IT_0132)*IT_0163) 
      + IT_0113*(IT_0067*conj(IT_0113) + IT_0070*IT_0151 + conj(IT_0068)*IT_0153
       + IT_0107*IT_0164 + conj(IT_0148)*IT_0165 + conj(IT_0132)*IT_0166 + conj
      (IT_0170)*IT_0171 + conj(IT_0172)*IT_0173) + IT_0047*(conj(IT_0068)
      *IT_0149 + conj(IT_0047)*IT_0179 + conj(IT_0148)*IT_0181 + conj(IT_0170)
      *IT_0183 + IT_0051*IT_0184 + conj(IT_0052)*IT_0189 + conj(IT_0172)*IT_0191
       + conj(IT_0132)*IT_0193) + IT_0170*(IT_0149*IT_0151 + conj(IT_0113)
      *IT_0171 + conj(IT_0170)*IT_0179 + conj(IT_0132)*IT_0181 + conj(IT_0047)
      *IT_0183 + conj(IT_0172)*IT_0189 + conj(IT_0052)*IT_0191 + conj(IT_0148)
      *IT_0193) + IT_0209 + IT_0210 + IT_0172*(IT_0150*IT_0151 + conj(IT_0113)
      *IT_0173 + conj(IT_0170)*IT_0189 + conj(IT_0047)*IT_0191 + conj(IT_0132)
      *IT_0202 + conj(IT_0148)*IT_0208 + conj(IT_0172)*IT_0212 + conj(IT_0052)
      *IT_0214) + IT_0052*(conj(IT_0068)*IT_0150 + IT_0056*IT_0184 + conj
      (IT_0047)*IT_0189 + conj(IT_0170)*IT_0191 + conj(IT_0148)*IT_0202 + conj
      (IT_0132)*IT_0208 + conj(IT_0052)*IT_0212 + conj(IT_0172)*IT_0214) +
       IT_0086*(IT_0070*conj(IT_0113) + IT_0215 + IT_0216 + IT_0217 + IT_0218 +
       IT_0219 + IT_0220) + IT_0099*(IT_0070*conj(IT_0113) + IT_0067*IT_0151 +
       IT_0215 + IT_0216 + IT_0217 + IT_0218 + IT_0219 + IT_0220);
    return create_ccomplex_return(IT_0221);
}

