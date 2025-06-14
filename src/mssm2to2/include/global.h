#ifndef CSL_LIB_GLOBAL
#define CSL_LIB_GLOBAL
#include "params.h"
#include "common.h"

namespace mssm2to2 {

void updateSpectrum(param_t &params);

struct SpectrumInput;
struct SpectrumOutput;

SpectrumOutput updateDiagonalization(SpectrumInput const&);

void updateDiagonalization(param_t &params);

////////////////////////////////////////////////////
// Here are the parameters to set before calling    
// updateDiagonalization()                          
////////////////////////////////////////////////////
struct SpectrumInput {
    complex_t A_b;
    complex_t A_t;
    complex_t A_tau;
    complex_t M_1;
    complex_t M_2;
    complex_t M_q1L;
    complex_t M_q3L;
    complex_t M_qbR;
    complex_t M_qdR;
    complex_t M_qtR;
    complex_t M_tauL;
    complex_t M_tauR;
    complex_t V_cb;
    complex_t V_tb;
    complex_t beta;
    complex_t e_em;
    complex_t m_W;
    complex_t m_Z;
    complex_t m_b;
    complex_t m_t;
    complex_t m_tau;
    complex_t mu_h;
    complex_t theta_W;
};

////////////////////////////////////////////////////
// Here are the masses and mixings                 
// result of the diagonalization                   
////////////////////////////////////////////////////
struct SpectrumOutput {
    real_t m_N_1;
    real_t m_N_2;
    real_t m_N_3;
    real_t m_N_4;
    real_t m_C_1;
    real_t m_C_2;
    real_t m_st_1;
    real_t m_st_2;
    real_t m_sb_1;
    real_t m_sb_2;
    real_t m_stau_1;
    real_t m_stau_2;

    complex_t N_B1;
    complex_t N_B2;
    complex_t N_B3;
    complex_t N_B4;
    complex_t N_W1;
    complex_t N_W2;
    complex_t N_W3;
    complex_t N_W4;
    complex_t N_d1;
    complex_t N_d2;
    complex_t N_d3;
    complex_t N_d4;
    complex_t N_u1;
    complex_t N_u2;
    complex_t N_u3;
    complex_t N_u4;
    complex_t U_Wm1;
    complex_t U_Wm2;
    complex_t U_d1;
    complex_t U_d2;
    complex_t U_sb_00;
    complex_t U_sb_01;
    complex_t U_sb_10;
    complex_t U_sb_11;
    complex_t U_st_00;
    complex_t U_st_01;
    complex_t U_st_10;
    complex_t U_st_11;
    complex_t U_stau_00;
    complex_t U_stau_01;
    complex_t U_stau_10;
    complex_t U_stau_11;
    complex_t V_Wp1;
    complex_t V_Wp2;
    complex_t V_u1;
    complex_t V_u2;
};

////////////////////////////////////////////////////
// Here is a generic function to read results      
// of the diagonalization in a corresponding struct
////////////////////////////////////////////////////

template<class Type>
void readDiagonalizationInputs(
        SpectrumInput &diagData,
        Type    const &input
        )
{
    diagData.A_b = input.A_b;
    diagData.A_t = input.A_t;
    diagData.A_tau = input.A_tau;
    diagData.M_1 = input.M_1;
    diagData.M_2 = input.M_2;
    diagData.M_q1L = input.M_q1L;
    diagData.M_q3L = input.M_q3L;
    diagData.M_qbR = input.M_qbR;
    diagData.M_qdR = input.M_qdR;
    diagData.M_qtR = input.M_qtR;
    diagData.M_tauL = input.M_tauL;
    diagData.M_tauR = input.M_tauR;
    diagData.V_cb = input.V_cb;
    diagData.V_tb = input.V_tb;
    diagData.beta = input.beta;
    diagData.e_em = input.e_em;
    diagData.m_W = input.m_W;
    diagData.m_Z = input.m_Z;
    diagData.m_b = input.m_b;
    diagData.m_t = input.m_t;
    diagData.m_tau = input.m_tau;
    diagData.mu_h = input.mu_h;
    diagData.theta_W = input.theta_W;
}

template<class Type>
void readDiagonalizationOutputs(
        SpectrumOutput const &diagData,
        Type                 &output
        )
{
    output.m_N_1 = diagData.m_N_1;
    output.m_N_2 = diagData.m_N_2;
    output.m_N_3 = diagData.m_N_3;
    output.m_N_4 = diagData.m_N_4;
    output.m_C_1 = diagData.m_C_1;
    output.m_C_2 = diagData.m_C_2;
    output.m_st_1 = diagData.m_st_1;
    output.m_st_2 = diagData.m_st_2;
    output.m_sb_1 = diagData.m_sb_1;
    output.m_sb_2 = diagData.m_sb_2;
    output.m_stau_1 = diagData.m_stau_1;
    output.m_stau_2 = diagData.m_stau_2;
    output.N_B1 = diagData.N_B1;
    output.N_B2 = diagData.N_B2;
    output.N_B3 = diagData.N_B3;
    output.N_B4 = diagData.N_B4;
    output.N_W1 = diagData.N_W1;
    output.N_W2 = diagData.N_W2;
    output.N_W3 = diagData.N_W3;
    output.N_W4 = diagData.N_W4;
    output.N_d1 = diagData.N_d1;
    output.N_d2 = diagData.N_d2;
    output.N_d3 = diagData.N_d3;
    output.N_d4 = diagData.N_d4;
    output.N_u1 = diagData.N_u1;
    output.N_u2 = diagData.N_u2;
    output.N_u3 = diagData.N_u3;
    output.N_u4 = diagData.N_u4;
    output.U_Wm1 = diagData.U_Wm1;
    output.U_Wm2 = diagData.U_Wm2;
    output.U_d1 = diagData.U_d1;
    output.U_d2 = diagData.U_d2;
    output.U_sb_00 = diagData.U_sb_00;
    output.U_sb_01 = diagData.U_sb_01;
    output.U_sb_10 = diagData.U_sb_10;
    output.U_sb_11 = diagData.U_sb_11;
    output.U_st_00 = diagData.U_st_00;
    output.U_st_01 = diagData.U_st_01;
    output.U_st_10 = diagData.U_st_10;
    output.U_st_11 = diagData.U_st_11;
    output.U_stau_00 = diagData.U_stau_00;
    output.U_stau_01 = diagData.U_stau_01;
    output.U_stau_10 = diagData.U_stau_10;
    output.U_stau_11 = diagData.U_stau_11;
    output.V_Wp1 = diagData.V_Wp1;
    output.V_Wp2 = diagData.V_Wp2;
    output.V_u1 = diagData.V_u1;
    output.V_u2 = diagData.V_u2;
}
void updateMassExpressions(param_t &params);


} // End of namespace mssm2to2

#endif
