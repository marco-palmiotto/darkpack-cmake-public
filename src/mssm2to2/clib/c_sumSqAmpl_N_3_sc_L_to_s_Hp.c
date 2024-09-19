#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_sc_L_to_s_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_sc_L_to_s_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_sc_L, 2);
    const ccomplex_t IT_0001 = s_13*IT_0000;
    const ccomplex_t IT_0002 = (-6)*IT_0001;
    const ccomplex_t IT_0003 = s_12*s_23;
    const ccomplex_t IT_0004 = 12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp1)
      *IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 2*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0009*IT_0015;
    const ccomplex_t IT_0017 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0021 = IT_0010*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = IT_0019 + IT_0022;
    const ccomplex_t IT_0024 = IT_0013 + IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = pow(m_s, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0000 + IT_0028 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp2)
      *IT_0007;
    const ccomplex_t IT_0033 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0034 = IT_0010*IT_0033;
    const ccomplex_t IT_0035 = 2*IT_0034;
    const ccomplex_t IT_0036 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0037 = IT_0016*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0040 = IT_0010*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = IT_0038 + IT_0041;
    const ccomplex_t IT_0043 = IT_0035 + IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0032*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0000 + IT_0028 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0031 + -IT_0049;
    const ccomplex_t IT_0051 = pow(m_c, 2);
    const ccomplex_t IT_0052 = pow(m_N_3, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_12 + -IT_0000 + IT_0051 + -IT_0052 
      + -reg_prop, -1);
    const ccomplex_t IT_0054 = pow(m_W, -1);
    const ccomplex_t IT_0055 = sin(beta);
    const ccomplex_t IT_0056 = cpow(IT_0055, -1);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0007*IT_0009*IT_0054*IT_0056;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_c*N_u3*e_em
      *IT_0007*IT_0054*IT_0056;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0058*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0050 + IT_0064;
    const ccomplex_t IT_0066 = conj(IT_0050) + conj(IT_0064);
    const ccomplex_t IT_0067 = IT_0007*IT_0055;
    const ccomplex_t IT_0068 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (-2)*IT_0069;
    const ccomplex_t IT_0071 = IT_0015*IT_0055;
    const ccomplex_t IT_0072 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0076 = IT_0067*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = IT_0074 + IT_0077;
    const ccomplex_t IT_0079 = IT_0070 + IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = cpow(IT_0009, -1);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_cs)*e_em*IT_0007*IT_0054*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0081*IT_0084;
    const ccomplex_t IT_0086 = IT_0047*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0089 = IT_0067*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0092 = IT_0071*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = IT_0090 + IT_0093;
    const ccomplex_t IT_0095 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0096 = IT_0067*IT_0095;
    const ccomplex_t IT_0097 = (-2)*IT_0096;
    const ccomplex_t IT_0098 = IT_0094 + IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_cs)*e_em*IT_0007*IT_0054*IT_0082;
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0100*IT_0102;
    const ccomplex_t IT_0104 = IT_0029*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = -IT_0087 + -IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0007*IT_0054*IT_0055*IT_0082;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B3)*e_em;
    const ccomplex_t IT_0110 = IT_0015*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W3)*e_em;
    const ccomplex_t IT_0113 = IT_0007*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + 3*IT_0114);
    const ccomplex_t IT_0116 = 0.166666666666667*IT_0115;
    const ccomplex_t IT_0117 = IT_0108*IT_0116;
    const ccomplex_t IT_0118 = IT_0053*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = IT_0106 + IT_0119;
    const ccomplex_t IT_0121 = conj(IT_0106) + conj(IT_0119);
    const ccomplex_t IT_0122 = m_s*m_N_3;
    const ccomplex_t IT_0123 = IT_0031 + IT_0049;
    const ccomplex_t IT_0124 = m_s*IT_0123;
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*(-2.82842712474619))*conj(V_cs)
      *e_em*IT_0007*(IT_0009*(m_W*IT_0055 + (-0.5)*IT_0051*IT_0054*IT_0056) + 
      -1./2*IT_0028*IT_0054*IT_0055*IT_0082);
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0007*IT_0054*IT_0082;
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = IT_0127*IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_13 + IT_0028 + IT_0052 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0132 = IT_0130*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0135 = IT_0026*IT_0102;
    const ccomplex_t IT_0136 = IT_0029*IT_0134*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0138 = IT_0060*IT_0108;
    const ccomplex_t IT_0139 = IT_0053*IT_0137*IT_0138;
    const ccomplex_t IT_0140 = m_N_3*IT_0119;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0007*IT_0054;
    const ccomplex_t IT_0142 = 0.5*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0015;
    const ccomplex_t IT_0144 = (-0.333333333333333)*IT_0143;
    const ccomplex_t IT_0145 = IT_0142*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_13 + IT_0028 + IT_0052 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0147 = IT_0145*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0150 = IT_0045*IT_0084;
    const ccomplex_t IT_0151 = IT_0047*IT_0149*IT_0150;
    const ccomplex_t IT_0152 = IT_0133 + -IT_0136 + -IT_0139 + IT_0140 +
       IT_0148 + -IT_0151;
    const ccomplex_t IT_0153 = IT_0087 + IT_0105;
    const ccomplex_t IT_0154 = m_s*IT_0153;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = IT_0032*IT_0081;
    const ccomplex_t IT_0157 = IT_0047*IT_0149*IT_0156;
    const ccomplex_t IT_0158 = m_N_3*IT_0063;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0111 + (-3)*IT_0114);
    const ccomplex_t IT_0160 = 0.166666666666667*IT_0159;
    const ccomplex_t IT_0161 = IT_0127*IT_0160;
    const ccomplex_t IT_0162 = IT_0131*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = IT_0058*IT_0116;
    const ccomplex_t IT_0165 = IT_0053*IT_0137*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0007*IT_0054*IT_0082;
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = IT_0142*IT_0167;
    const ccomplex_t IT_0169 = IT_0146*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = IT_0008*IT_0100;
    const ccomplex_t IT_0172 = IT_0029*IT_0134*IT_0171;
    const ccomplex_t IT_0173 = -IT_0157 + -IT_0158 + -IT_0163 + IT_0165 + 
      -IT_0170 + -IT_0172;
    const ccomplex_t IT_0174 = 6*IT_0152;
    const ccomplex_t IT_0175 = (-6)*IT_0125;
    const ccomplex_t IT_0176 = IT_0174 + IT_0175;
    const ccomplex_t IT_0177 = 6*conj(IT_0152);
    const ccomplex_t IT_0178 = (-6)*conj(IT_0125);
    const ccomplex_t IT_0179 = IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = s_23*m_N_3;
    const ccomplex_t IT_0181 = 6*IT_0173;
    const ccomplex_t IT_0182 = (-6)*IT_0155;
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = 6*conj(IT_0173);
    const ccomplex_t IT_0185 = (-6)*conj(IT_0155);
    const ccomplex_t IT_0186 = IT_0184 + IT_0185;
    const ccomplex_t IT_0187 = m_s*s_12;
    const ccomplex_t IT_0188 = 6*s_13;
    const ccomplex_t IT_0189 = m_s*m_N_3*IT_0000;
    const ccomplex_t IT_0190 = 6*IT_0064;
    const ccomplex_t IT_0191 = 6*conj(IT_0064);
    const ccomplex_t IT_0192 = (-6)*s_13;
    const ccomplex_t IT_0193 = IT_0005*(IT_0065*IT_0066 + IT_0120*IT_0121) + 6
      *IT_0122*((conj(IT_0125) + -conj(IT_0152))*IT_0155 + (IT_0125 + -IT_0152)
      *conj(IT_0155) + 0.166666666666667*conj(IT_0173)*IT_0176 +
       0.166666666666667*IT_0173*IT_0179) + IT_0180*(IT_0121*IT_0176 + IT_0120
      *IT_0179 + IT_0066*IT_0183 + IT_0065*IT_0186) + (IT_0066*IT_0176 + IT_0065
      *IT_0179 + IT_0121*IT_0183 + IT_0120*IT_0186)*IT_0187 + (IT_0125*conj
      (IT_0125) + IT_0152*conj(IT_0152) + IT_0155*conj(IT_0155) + IT_0173*conj
      (IT_0173))*IT_0188 + 6*IT_0189*(IT_0121*(IT_0050 + 0.166666666666667
      *IT_0190) + IT_0120*(conj(IT_0050) + 0.166666666666667*IT_0191)) + (conj
      (IT_0125)*IT_0152 + IT_0125*conj(IT_0152) + conj(IT_0155)*IT_0173 +
       IT_0155*conj(IT_0173))*IT_0192;
    return create_ccomplex_return(IT_0193);
}

