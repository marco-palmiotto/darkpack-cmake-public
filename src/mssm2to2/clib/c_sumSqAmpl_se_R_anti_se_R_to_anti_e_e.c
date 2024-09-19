#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_se_R_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_se_R_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = pow(m_e, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0003 + -IT_0010*IT_0011) + -IT_0000*IT_0002*IT_0003*IT_0005
      *IT_0012);
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = pow(m_se_R, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + (-2)*IT_0015 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0011;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0003*IT_0010 + IT_0001*IT_0011) + -IT_0000*IT_0002*IT_0005*IT_0011
      *IT_0012);
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0018 + -IT_0026;
    const ccomplex_t IT_0028 = cpow(IT_0008, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_e*N_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0012 + IT_0015 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0033*IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0028;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_e*N_d2*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0012 + IT_0015 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0041*IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0028;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_e*N_d3*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0012 + IT_0015 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0049*IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0028;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0012 + IT_0015 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0057*IT_0058*IT_0059;
    const ccomplex_t IT_0061 = IT_0004*IT_0028;
    const ccomplex_t IT_0062 = e_em*IT_0061;
    const ccomplex_t IT_0063 = IT_0005*IT_0008;
    const ccomplex_t IT_0064 = e_em*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0062 + -IT_0064);
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0028;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = m_e*IT_0071;
    const ccomplex_t IT_0073 = (-2)*IT_0070;
    const ccomplex_t IT_0074 = m_e*IT_0073;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0072 + IT_0075;
    const ccomplex_t IT_0077 = IT_0069*IT_0076;
    const ccomplex_t IT_0078 = IT_0067*IT_0070;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0076*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = IT_0032*IT_0082;
    const ccomplex_t IT_0084 = IT_0035*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = m_e*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d2)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = IT_0040*IT_0088;
    const ccomplex_t IT_0090 = IT_0043*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = m_e*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d3)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0048*IT_0094;
    const ccomplex_t IT_0096 = IT_0051*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = m_e*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1.4142135623731)*m_e*conj(N_d4)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0056*IT_0100;
    const ccomplex_t IT_0102 = IT_0059*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = m_e*IT_0103;
    const ccomplex_t IT_0105 = -IT_0036 + -IT_0044 + -IT_0052 + -IT_0060 + 
      -IT_0077 + IT_0080 + IT_0086 + IT_0092 + IT_0098 + IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0028;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = IT_0088*IT_0107;
    const ccomplex_t IT_0109 = IT_0042*IT_0043*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0028;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = IT_0094*IT_0111;
    const ccomplex_t IT_0113 = IT_0050*IT_0051*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0028;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = IT_0100*IT_0115;
    const ccomplex_t IT_0117 = IT_0058*IT_0059*IT_0116;
    const ccomplex_t IT_0118 = IT_0038*IT_0107;
    const ccomplex_t IT_0119 = IT_0043*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = m_e*IT_0120;
    const ccomplex_t IT_0122 = IT_0046*IT_0111;
    const ccomplex_t IT_0123 = IT_0051*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = m_e*IT_0124;
    const ccomplex_t IT_0126 = IT_0054*IT_0115;
    const ccomplex_t IT_0127 = IT_0059*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = m_e*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0028;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = IT_0082*IT_0131;
    const ccomplex_t IT_0133 = IT_0034*IT_0035*IT_0132;
    const ccomplex_t IT_0134 = IT_0030*IT_0131;
    const ccomplex_t IT_0135 = IT_0035*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = m_e*IT_0136;
    const ccomplex_t IT_0138 = IT_0077 + -IT_0080 + -IT_0109 + -IT_0113 + 
      -IT_0117 + IT_0121 + IT_0125 + IT_0129 + -IT_0133 + IT_0137;
    const ccomplex_t IT_0139 = -IT_0012;
    const ccomplex_t IT_0140 = s_34 + IT_0139;
    const ccomplex_t IT_0141 = IT_0069*IT_0073;
    const ccomplex_t IT_0142 = -IT_0085 + -IT_0091 + -IT_0097 + -IT_0103 + 
      -IT_0141;
    const ccomplex_t IT_0143 = m_e*s_13;
    const ccomplex_t IT_0144 = m_e*s_14;
    const ccomplex_t IT_0145 = IT_0073*IT_0079;
    const ccomplex_t IT_0146 = -IT_0120 + -IT_0124 + -IT_0128 + -IT_0136 + 
      -IT_0145;
    const ccomplex_t IT_0147 = IT_0012*IT_0015;
    const ccomplex_t IT_0148 = s_13*s_14;
    const ccomplex_t IT_0149 = s_34*IT_0015;
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = IT_0148 + IT_0150;
    const ccomplex_t IT_0152 = 2*s_34;
    const ccomplex_t IT_0153 = (-0.5)*IT_0147;
    const ccomplex_t IT_0154 = IT_0148 + IT_0150 + IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = 2*IT_0155;
    const ccomplex_t IT_0158 = IT_0156*IT_0157;
    const ccomplex_t IT_0159 = cpow(s_12 + IT_0015 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0160 = IT_0158*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = 0.5*IT_0161;
    const ccomplex_t IT_0163 = cpow(IT_0162, 2);
    const ccomplex_t IT_0164 = conj(IT_0142) + conj(IT_0146);
    const ccomplex_t IT_0165 = 4*IT_0162;
    const ccomplex_t IT_0166 = IT_0142 + IT_0146;
    const ccomplex_t IT_0167 = (-4)*IT_0162;
    const ccomplex_t IT_0168 = (-2)*IT_0105;
    const ccomplex_t IT_0169 = (-2)*conj(IT_0105);
    const ccomplex_t IT_0170 = (-2)*IT_0138;
    const ccomplex_t IT_0171 = -s_14;
    const ccomplex_t IT_0172 = s_13 + IT_0171;
    const ccomplex_t IT_0173 = m_e*IT_0172;
    const ccomplex_t IT_0174 = 2*(conj(IT_0027)*(IT_0105 + IT_0138) + 2
      *IT_0027*(conj(IT_0027) + 0.5*conj(IT_0105) + 0.5*conj(IT_0138)))*IT_0140 
      + 2*IT_0142*(conj(IT_0105)*IT_0143 + -conj(IT_0138)*IT_0144 + -conj
      (IT_0146)*IT_0147 + 2*conj(IT_0142)*IT_0151) + 2*IT_0105*(conj(IT_0142)
      *IT_0143 + 0.5*conj(IT_0105)*IT_0152) + (-8)*IT_0154*(IT_0163 + (-0.125)
      *IT_0164*IT_0165 + (-0.125)*IT_0166*IT_0167) + (IT_0012*conj(IT_0138) +
       IT_0144*conj(IT_0146))*IT_0168 + 2*IT_0138*(IT_0143*conj(IT_0146) + 0.5
      *conj(IT_0138)*IT_0152 + 0.5*IT_0012*IT_0169) + 2*IT_0146*(conj(IT_0138)
      *IT_0143 + -conj(IT_0142)*IT_0147 + 2*conj(IT_0146)*IT_0151 + 0.5*IT_0144
      *IT_0169) + conj(IT_0142)*IT_0144*IT_0170 + 2*(IT_0027*IT_0164 + conj
      (IT_0027)*IT_0166 + (-2)*IT_0162*(IT_0027 + -conj(IT_0027) + (-0.5)*conj
      (IT_0105) + (-0.5)*conj(IT_0138) + (-0.25)*IT_0168 + (-0.25)*IT_0170))
      *IT_0173;
    return create_ccomplex_return(IT_0174);
}

