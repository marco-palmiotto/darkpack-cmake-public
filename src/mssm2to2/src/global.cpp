#include "global.h"
#include "libcomplexop.h"
#include "libdiagonalization.h"
#include "mssm2to2.h"

namespace mssm2to2
{


  void updateSpectrum(param_t& params)
  {
    updateDiagonalization(params);
    updateMassExpressions(params);
  }

  void updateDiagonalization(param_t& params)
  {
    SpectrumInput inputs;
    readDiagonalizationInputs(inputs, params);
    SpectrumOutput outputs = updateDiagonalization(inputs);
    readDiagonalizationOutputs(outputs, params);
  }

  SpectrumOutput updateDiagonalization(SpectrumInput const& inputs)
  {
    auto const& A_b = inputs.A_b;
    auto const& A_t = inputs.A_t;
    auto const& A_tau = inputs.A_tau;
    auto const& M_1 = inputs.M_1;
    auto const& M_2 = inputs.M_2;
    auto const& M_q1L = inputs.M_q1L;
    auto const& M_q3L = inputs.M_q3L;
    auto const& M_qbR = inputs.M_qbR;
    auto const& M_qdR = inputs.M_qdR;
    auto const& M_qtR = inputs.M_qtR;
    auto const& M_tauL = inputs.M_tauL;
    auto const& M_tauR = inputs.M_tauR;
    auto const& V_cb = inputs.V_cb;
    auto const& V_tb = inputs.V_tb;
    auto const& beta = inputs.beta;
    auto const& e_em = inputs.e_em;
    auto const& m_W = inputs.m_W;
    auto const& m_Z = inputs.m_Z;
    auto const& m_b = inputs.m_b;
    auto const& m_t = inputs.m_t;
    auto const& m_tau = inputs.m_tau;
    auto const& mu_h = inputs.mu_h;
    auto const& theta_W = inputs.theta_W;

    SpectrumOutput outputs;

    Diagonalizer::applyDiagonalization(
        {
            M_1,
            0,
            (-2) * m_Z * std::cos(beta) * std::sin(theta_W),
            2 * m_Z * std::sin(beta) * std::sin(theta_W),
            0,
            M_2,
            2 * m_W * std::cos(beta),
            (-2) * m_W * std::sin(beta),
            0,
            0,
            0,
            (-2) * mu_h,
            0,
            0,
            0,
            0,
        },
        {
            &outputs.N_B1,
            &outputs.N_B2,
            &outputs.N_B3,
            &outputs.N_B4,
            &outputs.N_W1,
            &outputs.N_W2,
            &outputs.N_W3,
            &outputs.N_W4,
            &outputs.N_d1,
            &outputs.N_d2,
            &outputs.N_d3,
            &outputs.N_d4,
            &outputs.N_u1,
            &outputs.N_u2,
            &outputs.N_u3,
            &outputs.N_u4,
        },
        {
            &outputs.m_N_1,
            &outputs.m_N_2,
            &outputs.m_N_3,
            &outputs.m_N_4,
        });

    Diagonalizer::applyBiDiagonalization(
        {
            M_2,
            1.41421 * m_W * std::sin(beta),
            1.41421 * m_W * std::cos(beta),
            mu_h,
        },
        {
            &outputs.V_Wp1,
            &outputs.V_Wp2,
            &outputs.V_u1,
            &outputs.V_u2,
        },
        {
            &outputs.U_Wm1,
            &outputs.U_Wm2,
            &outputs.U_d1,
            &outputs.U_d2,
        },
        {
            &outputs.m_C_1,
            &outputs.m_C_2,
        });

    Diagonalizer::applyDiagonalization(
        {
            std::pow(m_t, 2) + std::pow(M_q3L, 2) + (-0.5) * std::pow(m_W, 2) * std::pow(std::sin(beta), 2) +
                (-0.5) * std::pow(m_W, 2) * (1 + -std::pow(std::sin(beta), 2)) +
                (-0.166667) * std::pow(m_Z, 2) * std::pow(std::sin(beta), 2) * std::pow(std::sin(theta_W), 2) +
                0.166667 * std::pow(m_Z, 2) * (1 + -std::pow(std::sin(beta), 2)) * std::pow(std::sin(theta_W), 2),
            (-2) * m_t * mu_h * std::cos(beta) / std::sin(beta) +
                2.82843 * A_t * m_W * std::sin(beta) * std::sin(theta_W) / e_em,
            0,
            std::pow(m_t, 2) + std::pow(M_qtR, 2) +
                (-0.666667) * std::pow(m_Z, 2) * std::pow(std::sin(beta), 2) * std::pow(std::sin(theta_W), 2) +
                0.666667 * std::pow(m_Z, 2) * (1 + -std::pow(std::sin(beta), 2)) * std::pow(std::sin(theta_W), 2),
        },
        {
            &outputs.U_st_00,
            &outputs.U_st_01,
            &outputs.U_st_10,
            &outputs.U_st_11,
        },
        {
            &outputs.m_st_1,
            &outputs.m_st_2,
        });
    if (0 > outputs.m_st_1)
    {
      std::cerr << "Warning: negative squared mass for " << "m_st_1" << ".\n";
    }
    outputs.m_st_1 = std::sqrt(outputs.m_st_1);
    if (0 > outputs.m_st_2)
    {
      std::cerr << "Warning: negative squared mass for " << "m_st_2" << ".\n";
    }
    outputs.m_st_2 = std::sqrt(outputs.m_st_2);

    Diagonalizer::applyDiagonalization(
        {
            1.45924e-05 * std::pow(M_q1L, 2) + std::pow(V_cb, 2) * std::pow(M_q1L, 2) +
                V_tb * std::conj(V_tb) * std::pow(M_q3L, 2) + (-0.5) * std::pow(m_W, 2) * std::pow(std::sin(beta), 2) +
                std::pow(m_b, 2) * std::pow(std::cos(beta), -2) * (1 + -std::pow(std::sin(beta), 2)) +
                (-0.166667) * std::pow(m_Z, 2) * std::pow(std::sin(beta), 2) * std::pow(std::sin(theta_W), 2),
            (-2) * m_b * mu_h * std::sin(beta) / std::cos(beta) +
                2.82843 * A_b * m_W * V_tb * std::conj(V_tb) * std::cos(beta) * std::sin(theta_W) / e_em,
            0,
            V_tb * std::conj(V_tb) * std::pow(M_qbR, 2) + 1.45924e-05 * std::pow(M_qdR, 2) +
                std::pow(V_cb, 2) * std::pow(M_qdR, 2) +
                std::pow(m_b, 2) * std::pow(std::cos(beta), -2) * (1 + -std::pow(std::sin(beta), 2)) +
                0.333333 * std::pow(m_Z, 2) * std::pow(std::sin(beta), 2) * std::pow(std::sin(theta_W), 2) +
                (-0.333333) * std::pow(m_Z, 2) * (1 + -std::pow(std::sin(beta), 2)) * std::pow(std::sin(theta_W), 2),
        },
        {
            &outputs.U_sb_00,
            &outputs.U_sb_01,
            &outputs.U_sb_10,
            &outputs.U_sb_11,
        },
        {
            &outputs.m_sb_1,
            &outputs.m_sb_2,
        });
    if (0 > outputs.m_sb_1)
    {
      std::cerr << "Warning: negative squared mass for " << "m_sb_1" << ".\n";
    }
    outputs.m_sb_1 = std::sqrt(outputs.m_sb_1);
    if (0 > outputs.m_sb_2)
    {
      std::cerr << "Warning: negative squared mass for " << "m_sb_2" << ".\n";
    }
    outputs.m_sb_2 = std::sqrt(outputs.m_sb_2);

    Diagonalizer::applyDiagonalization(
        {
            std::pow(M_tauL, 2) + (-0.5) * std::pow(m_W, 2) * std::pow(std::sin(beta), 2) +
                (-0.5) * std::pow(m_W, 2) * (1 + -std::pow(std::sin(beta), 2)) +
                std::pow(m_tau, 2) * std::pow(std::cos(beta), -2) * (1 + -std::pow(std::sin(beta), 2)) +
                0.5 * std::pow(m_Z, 2) * std::pow(std::sin(beta), 2) * std::pow(std::sin(theta_W), 2) +
                (-0.5) * std::pow(m_Z, 2) * (1 + -std::pow(std::sin(beta), 2)) * std::pow(std::sin(theta_W), 2),
            (-2) * mu_h * m_tau * std::sin(beta) / std::cos(beta) +
                2.82843 * m_W * A_tau * std::cos(beta) * std::sin(theta_W) / e_em,
            0,
            std::pow(M_tauR, 2) +
                std::pow(m_tau, 2) * std::pow(std::cos(beta), -2) * (1 + -std::pow(std::sin(beta), 2)) +
                std::pow(m_Z, 2) * std::pow(std::sin(beta), 2) * std::pow(std::sin(theta_W), 2) +
                -std::pow(m_Z, 2) * (1 + -std::pow(std::sin(beta), 2)) * std::pow(std::sin(theta_W), 2),
        },
        {
            &outputs.U_stau_00,
            &outputs.U_stau_01,
            &outputs.U_stau_10,
            &outputs.U_stau_11,
        },
        {
            &outputs.m_stau_1,
            &outputs.m_stau_2,
        });
    if (0 > outputs.m_stau_1)
    {
      std::cerr << "Warning: negative squared mass for " << "m_stau_1" << ".\n";
    }
    outputs.m_stau_1 = std::sqrt(outputs.m_stau_1);
    if (0 > outputs.m_stau_2)
    {
      std::cerr << "Warning: negative squared mass for " << "m_stau_2" << ".\n";
    }
    outputs.m_stau_2 = std::sqrt(outputs.m_stau_2);

    return outputs;
  }

  void updateMassExpressions(param_t& params)
  {
    params.m_Hp = m_Hp(params).real();
    params.m_h = m_h(params).real();
    params.m_H0 = m_H0(params).real();
    params.m_sc_L = m_sc_L(params).real();
    params.m_su_L = m_su_L(params).real();
    params.m_sc_R = m_sc_R(params).real();
    params.m_su_R = m_su_R(params).real();
    params.m_sd_R = m_sd_R(params).real();
    params.m_ss_R = m_ss_R(params).real();
    params.m_sd_L = m_sd_L(params).real();
    params.m_ss_L = m_ss_L(params).real();
    params.m_se_R = m_se_R(params).real();
    params.m_smu_R = m_smu_R(params).real();
    params.m_se_L = m_se_L(params).real();
    params.m_smu_L = m_smu_L(params).real();
    params.m_snu_e = m_snu_e(params).real();
    params.m_snu_tau = m_snu_tau(params).real();
  }

} // End of namespace mssm2to2
