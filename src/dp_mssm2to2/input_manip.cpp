#include "input_manip.hpp"

namespace mssm2to2::input_manip
{

  void FixTrilinears(param_t& input)
  {
#ifdef DEBUG
    std::cout << "Calling input_manip::FixTrilinears\n";
#endif
    // Those are defined at the m_b(m_b) scale
    const double prefactor = input.e_em / (std::sqrt(2) * std::cos(input.beta) * std::sin(input.theta_W) * input.m_W);

    input.A_tau = prefactor * input.m_tau * input.A_tau;
    input.A_t = input.e_em / (std::sqrt(2) * std::sin(input.beta) * std::sin(input.theta_W) * input.m_W) * input.m_t *
                input.A_t;
    input.A_b = prefactor *
                std::abs(input.m_d * input.V_td * std::conj(input.V_td.get()) +
                         input.m_s * input.V_ts * std::conj(input.V_ts.get()) +
                         input.m_b * input.V_tb * std::conj(input.V_tb.get())) *
                input.A_b;
  }

  void CKMIdentity(param_t& input)
  {
    std::cout << "Calling input_manip::CKMIdentity\n";
    input.V_ud = 1.0;
    input.V_us = 0.0;
    input.V_ub_mod = 0.0;
    input.delta_wolf = 0.0;

    input.V_cd = 0.0;
    input.V_cs = 1.0;
    input.V_cb = 0.0;

    input.V_td = 0.0;
    input.V_ts = 0.0;
    input.V_tb = 1.0;
    return;
  }

  bool CheckUnitarity(param_t& input)
  {
    // Defining the CKM matrix

    complex_t V[3][3];

    V[0][0] = input.V_ud;
    V[0][1] = input.V_us;
    V[0][2] = std::polar(double(input.V_ub_mod), double(input.delta_wolf));

    V[1][0] = input.V_cd;
    V[1][1] = input.V_cs;
    V[1][2] = input.V_cb;

    V[2][0] = input.V_td;
    V[2][1] = input.V_ts;
    V[2][2] = input.V_tb;

    bool unitary = true;
    complex_t sum1, sum2;

    for (auto j = 0; j < 3; j++)
    {
      for (auto k = 0; k < 3; k++)
      {
        sum1 = 0.0;
        sum2 = 0.0;
        for (auto i = 0; i < 3; i++)
        {
          sum1 += V[i][j] * std::conj(V[i][k]);
          sum2 += V[j][i] * std::conj(V[k][i]);
        }

        std::cout << " j = " << j << "  k = " << k << " sum1 = " << sum1 << " sum2 = " << sum2 << std::endl;


        if (k == j && (std::abs(sum1 - 1.0) > 1.0e-12 || std::abs(sum2 - 1.0) > 1.0e-12))
          unitary = false;
        else if (std::abs(sum1) > 1.0e-12 || std::abs(sum2) > 1.0e-12)
          unitary = false;

        if (!unitary)
          break;
      }
    }

    return unitary;
  }

  bool ImposeUnitarity(param_t& input)
  {
    int i = 0;
    if (!CheckUnitarity(input))
    {
      do
      {
        // I recreate the matrix using the Wolfenstein parametrization
        constexpr const double lambda = 0.22650;
        constexpr const double A = 0.790;
        constexpr const double rhobar = 0.141;
        constexpr const double etabar = 0.357;
        constexpr const complex_t im_u{0.0, 1.0};

        const double rho = rhobar / (1 - lambda * lambda / 2);
        const double eta = etabar / (1 - lambda * lambda / 2);

        input.V_ud = 1 - std::pow(lambda, 2) / 2;
        input.V_us = lambda;
        const complex_t Vub = std::abs(A * std::pow(lambda, 3) * (rho - im_u * eta));
        input.V_ub_mod = std::abs(Vub);

        input.V_cd = -lambda;
        input.V_cs = -lambda * lambda / 2;
        input.V_cb = A * lambda * lambda;

        input.V_td = A * std::pow(lambda, 3) * (1 - rho - im_u * eta);
        input.V_ts = -A * lambda * lambda;
        input.V_tb = 1;

        input.delta_wolf = std::atan(etabar / rhobar); // Angle between eta_bar and rho_bar

        i++;

        std::cout << "CKM matrix is reconstructed" << std::endl;
      } while (!CheckUnitarity(input) && i < 4);
    }

    if (!CheckUnitarity(input))
      exit(2);
    else
      return true;
  }

  void imposePDG2016Values(mssm2to2::Param_t& input)
  {
    input.mbmb = pdg2016Value::mb_mb;

    input.m_el_pole = pdg2016Value::m_el_pole;
    input.m_mu_pole = pdg2016Value::m_mu_pole;
    input.m_d2GeV = pdg2016Value::m_d2GeV;
    input.m_u2GeV = pdg2016Value::m_u2GeV;
    input.m_s2GeV = pdg2016Value::m_s2GeV;
    input.m_c_m_c = pdg2016Value::m_c_m_c;

    input.e_em = std::sqrt(4. * M_PI / pdg2016Value::alpha_em_inv);

    input.m_e = pdg2016Value::m_el_pole;
    input.m_mu = pdg2016Value::m_mu_pole;
    //       input.m_tau = pdg2016Value::m_tau_pole;

    input.m_u = pdg2016Value::m_u2GeV;
    input.m_d = pdg2016Value::m_d2GeV;

    input.m_s = pdg2016Value::m_s2GeV;
    input.m_c = pdg2016Value::m_c_m_c;
    input.m_b = pdg2016Value::mb_mb;

    input.theta_W = pdg2016Value::theta_W;
    input.delta_wolf = pdg2016Value::delta_wolf;
    input.V_ud = pdg2016Value::V_ud;
    input.V_us = pdg2016Value::V_us;
    input.V_ub_mod = pdg2016Value::V_ub_mod;
    input.V_cd = pdg2016Value::V_cd;
    input.V_cs = pdg2016Value::V_cs;
    input.V_cb = pdg2016Value::V_cb;
    input.V_td = pdg2016Value::V_td;
    input.V_ts = pdg2016Value::V_ts;
    input.V_tb = pdg2016Value::V_tb;

#ifdef WIDTH_Z
    input.Gamma_Z = pdg2016Value::Gamma_Z;
#endif
#ifdef WIDTH_W
    input.Gamma_W = pdg2016Value::Gamma_W;
#endif
  }

} // namespace mssm2to2::input_manip