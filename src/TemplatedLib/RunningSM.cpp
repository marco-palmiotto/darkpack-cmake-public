#include "RunningSM.hpp"
#include <algorithm>
// Needed for TestRunning masses
#include <fstream>
#include <iomanip>

// #define NO_RUN

namespace __SPEC_LIB_NAME__
{

  // Mathematical macros
  constexpr const double ZETA3 = 1.2020569031595942855;
  constexpr const double NEGLIGIBLE_ENERGY_DIFFERENCE = 1.0e-3;

  template <typename T> inline T SQUARE(const T x) { return x * x; }

  void RunningSM::init(const Param_t& input)
  {
    setSanitNames();
    int i;
    runlightquarks = false;
    runcharm = true;
    higgsloops = true;

    // Initialising class fixed parametrs
    alphas_MZ = input.alpha_str_Mz;
    mass_Z = input.m_Z_pole;

    m_u2GeV = input.m_u2GeV;
    m_d2GeV = input.m_d2GeV;
    m_s2GeV = input.m_s2GeV;
    mc_mc = input.m_c_m_c;

    mb_mb = input.mbmb;

    mass_top_pole = input.mtop_pole;

    QuarkMass[UP] = input.m_u2GeV;
    QuarkMass[DOWN] = input.m_d2GeV;
    QuarkMass[STRANGE] = input.m_s2GeV;
    for (i = STRANGE; i >= UP; i--)
      QuarkMassQ[i] = 2.;

    QuarkMass[CHARM] = input.m_c_m_c;
    QuarkMassQ[CHARM] = input.m_c_m_c;

    QuarkMass[BEAUTY] = mb_mb;
    QuarkMassQ[BEAUTY] = mb_mb;

    alphas = alphas_MZ;
    Q_alphas = mass_Z;

    Lambda3 = Lambda4 = Lambda5 = Lambda6 = 0.;
    alphasMZ_Lambda3 = alphasMZ_Lambda4 = alphasMZ_Lambda5 = alphasMZ_Lambda6 = alphas_MZ;

    QuarkMass[TOP] = GetMtopMtop();
    QuarkMassQ[TOP] = QuarkMass[TOP];

    mass_b_pole = GetMbPole();
    mass_c_pole = GetMcPole();

#ifdef DEBUG
    std::cout << "Running parameters assigned, setting LastQuarkMass, and calling HandleParamRunning\n";
#endif
    if (input.Running_scale < 0.)
    {
      for (i = TOP; i >= UP; i--)
      {
        LastQuarkMass[i] = QuarkMass[i];
        LastQuarkMassQ[i] = QuarkMassQ[i];
      }
      Param_t input_copy(input);
      HandleParamRunning(input_copy, mass_Z);
    }
    else
    {
      LastQuarkMass[UP] = input.m_u;
      LastQuarkMass[DOWN] = input.m_d;
      LastQuarkMass[STRANGE] = input.m_s;
      LastQuarkMass[CHARM] = input.m_c;
      LastQuarkMass[BEAUTY] = input.m_b;
      LastQuarkMass[TOP] = input.m_t;
      for (i = TOP; i >= UP; i--)
        LastQuarkMassQ[i] = input.Running_scale;
    }
#ifdef DEBUG
    std::cout << "Running class correctly built\n";
#endif
    return;
  }

  // Constructor of the class (automatic parameters)
  RunningSM::RunningSM()
  {
    setSanitNames();
    runlightquarks = false;
    runcharm = true;
    higgsloops = true;

    // Initialising class fixed parametrs
    alphas_MZ = pdgValue::alpha_str_Mz;
    mass_Z = pdgValue::m_Z_pole;

    m_u2GeV = pdgValue::m_u2GeV;
    m_d2GeV = pdgValue::m_d2GeV;
    m_s2GeV = pdgValue::m_s2GeV;
    mc_mc = pdgValue::m_c_m_c;

    mb_mb = pdgValue::mb_mb;

    mass_top_pole = pdgValue::m_top_pole;

    QuarkMass[UP] = pdgValue::m_u2GeV;
    QuarkMass[DOWN] = pdgValue::m_d2GeV;
    QuarkMass[STRANGE] = pdgValue::m_s2GeV;
    for (short int i = STRANGE; i >= UP; i--)
      QuarkMassQ[i] = 2.;

    QuarkMass[CHARM] = pdgValue::m_c_m_c;
    QuarkMassQ[CHARM] = pdgValue::m_c_m_c;

    QuarkMass[BEAUTY] = mb_mb;
    QuarkMassQ[BEAUTY] = mb_mb;

    alphas = alphas_MZ;
    Q_alphas = mass_Z;

    Lambda3 = Lambda4 = Lambda5 = Lambda6 = 0.;
    alphasMZ_Lambda3 = alphasMZ_Lambda4 = alphasMZ_Lambda5 = alphasMZ_Lambda6 = alphas_MZ;


#ifdef TRY_CONSTEXPR
    QuarkMass[TOP] = pdgValue::m_top_m_top_computed;
#else
    QuarkMass[TOP] = GetMtopMtop();
#endif
    QuarkMassQ[TOP] = QuarkMass[TOP];

    for (short int i = TOP; i >= UP; i--)
    {
      LastQuarkMass[i] = QuarkMass[i];
      LastQuarkMassQ[i] = QuarkMassQ[i];
    }

#ifdef TRY_CONSTEXPR
    mass_b_pole = pdgValue::mb_pole;
    mass_c_pole = pdgValue::m_c_pole_computed;
#else
    mass_b_pole = GetMbPole();
    mass_c_pole = GetMcPole();
#endif
  }

  bool RunningSM::operator==(const RunningSM& other) const
  {
    if (this->alphas_MZ != other.alphas_MZ)
      return false;
    if (this->mass_Z != other.mass_Z)
      return false;
    if (this->m_d2GeV != other.m_d2GeV)
      return false;
    if (this->m_u2GeV != other.m_u2GeV)
      return false;
    if (this->m_s2GeV != other.m_s2GeV)
      return false;
    if (this->mc_mc != other.mc_mc)
      return false;
    if (this->mb_mb != other.mb_mb)
      return false;
    if (this->mass_top_pole != other.mass_top_pole)
      return false;
    if (this->mass_b_pole != other.mass_b_pole)
      return false;
    if (this->mass_c_pole != other.mass_c_pole)
      return false;
    for (int i = 0; i < 6; i++)
    {
      if (this->QuarkMass[i] != other.QuarkMass[i])
        return false;
      if (this->QuarkMassQ[i] != other.QuarkMassQ[i])
        return false;
      if (this->LastQuarkMass[i] != other.LastQuarkMass[i])
        return false;
      if (this->LastQuarkMassQ[i] != other.LastQuarkMassQ[i])
        return false;
    }
    if (this->alphas != other.alphas)
      return false;
    if (this->Q_alphas != other.Q_alphas)
      return false;
    if (this->Lambda3 != other.Lambda3)
      return false;
    if (this->alphasMZ_Lambda3 != other.alphasMZ_Lambda3)
      return false;
    if (this->Lambda4 != other.Lambda4)
      return false;
    if (this->alphasMZ_Lambda4 != other.alphasMZ_Lambda4)
      return false;
    if (this->Lambda5 != other.Lambda5)
      return false;
    if (this->alphasMZ_Lambda5 != other.alphasMZ_Lambda5)
      return false;
    if (this->Lambda6 != other.Lambda6)
      return false;
    if (this->alphasMZ_Lambda6 != other.alphasMZ_Lambda6)
      return false;
    return true;
  }

  void RunningSM::print(std::ostream& out) const
  {
    out << (alphas_MZ);
    out << (mass_Z);
    out << (m_d2GeV);
    out << (m_u2GeV);
    out << (m_s2GeV);
    out << (mc_mc);
    out << (mb_mb);
    out << (mass_top_pole);
    out << (mass_b_pole);
    out << (mass_c_pole);

    for (int i = 0; i < 6; i++)
    {
      out << QuarkMass[i];
      out << (QuarkMassQ[i]);
      out << (LastQuarkMass[i]);
      out << (LastQuarkMassQ[i]);
    }

    out << (alphas);
    out << (Q_alphas);
    out << (Lambda3);
    out << (alphasMZ_Lambda3);
    out << (Lambda4);
    out << (alphasMZ_Lambda4);
    out << (Lambda5);
    out << (alphasMZ_Lambda5);
    out << (Lambda6);
    out << (alphasMZ_Lambda6);
  }

  // computes the QCD coupling constant alphas at the energy Q, fixing at mtop and mbot the matching scales between
  // the flavours
  real_t RunningSM::AlphaStrong(real_t Q, real_t mtop)
  {
#ifdef DEBUG
    std::cout << "RunningSM::AlphaStrong: Called AlphaStrong(Q = " << Q << " , mtop = " << mtop << ")\n";
#endif
    if (std::fabs(Q - Q_alphas) < NEGLIGIBLE_ENERGY_DIFFERENCE)
      return alphas;

#ifdef DEBUG
    std::cout << "RunningSM::AlphaStrong:  Doing the computation \n";
#endif

    const real_t mbot = mb_mb;
    if (mtop < 0.0)
      mtop = mass_top_pole;
    real_t alphas_running;

    real_t Lambda_min, Lambda_max, Lambda_avg, alphas_min, alphas_avg;

    int nf = 5;
    real_t beta0 = 11. - 2. / 3. * nf;
    real_t beta0_sq = beta0 * beta0;
    real_t beta1 = 51. - 19. / 3. * nf;
    real_t beta1_sq = beta1 * beta1;
    real_t beta2 = 2857. - 5033. * nf / 9. + 325. / 27. * nf * nf;

    if (Lambda5 == -1.)
      return -1.;

    if (Lambda5 == 0. || std::fabs(1. - alphasMZ_Lambda5 / alphas_MZ) > 1.e-5)
    {
      Lambda_min = 1.e-3;
      Lambda_max = 1.;
      alphas_min = 0.;

      while ((std::fabs(1. - alphas_min / alphas_MZ) >= 1.e-4) && (std::fabs(1. - Lambda_min / Lambda_max) > 1.e-5))
      {
        real_t log_MZ_Lambda_min = std::log(std::pow(mass_Z / Lambda_min, 2.));

        alphas_min =
            4. * M_PI / beta0 / log_MZ_Lambda_min *
            (1. - 2. * beta1 / beta0_sq * std::log(log_MZ_Lambda_min) / log_MZ_Lambda_min +
             4. * beta1_sq / std::pow(beta0_sq * log_MZ_Lambda_min, 2.) *
                 (std::pow(std::log(log_MZ_Lambda_min) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

        Lambda_avg = (Lambda_min + Lambda_max) / 2.;
        real_t log_MZ_Lambda_avg = std::log(std::pow(mass_Z / Lambda_avg, 2.));

        alphas_avg =
            4. * M_PI / beta0 / log_MZ_Lambda_avg *
            (1. - 2. * beta1 / beta0_sq * std::log(log_MZ_Lambda_avg) / log_MZ_Lambda_avg +
             4. * beta1_sq / std::pow(beta0_sq * log_MZ_Lambda_avg, 2.) *
                 (std::pow(std::log(log_MZ_Lambda_avg) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

        if ((alphas_MZ >= alphas_min) && (alphas_MZ <= alphas_avg))
          Lambda_max = Lambda_avg;
        else
          Lambda_min = Lambda_avg;
      } // End of the while

      Lambda5 = Lambda_min;
      alphasMZ_Lambda5 = alphas_MZ;

      if (std::fabs(1. - Lambda_min / Lambda_max) <= 1.e-5)
      {
        Lambda5 = -1.;
        return -1.;
      }
    } // End of the 1st if

    /*In the end of this if cycle Lmabda5L (local) is or equal to the one in the original class, or
     *equal to something determined in the if. In every case, once the if statement contentent is
     *executed, the two values coincide*/

    if ((Q <= mtop) && (Q >= mbot))
    /* 5 active flavours */
    {
      alphas_running =
          4. * M_PI / beta0 / std::log(std::pow(Q / Lambda5, 2.)) *
          (1. -
           2. * beta1 / beta0_sq * std::log(std::log(std::pow(Q / Lambda5, 2.))) / std::log(std::pow(Q / Lambda5, 2.)) +
           4. * beta1_sq / std::pow(beta0_sq * std::log(std::pow(Q / Lambda5, 2.)), 2.) *
               (std::pow(std::log(std::log(std::pow(Q / Lambda5, 2.))) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq -
                5. / 4.));
      alphas = alphas_running;
      Q_alphas = Q;
      return alphas_running;
    }
    else if ((Q > mtop))
    /* 6 active flavours */
    {
      alphas_running = 4. * M_PI / beta0 / std::log(std::pow(mtop / Lambda5, 2.)) *
                       (1. -
                        2. * beta1 / beta0_sq * std::log(std::log(std::pow(mtop / Lambda5, 2.))) /
                            std::log(std::pow(mtop / Lambda5, 2.)) +
                        4. * beta1_sq / std::pow(beta0_sq * std::log(std::pow(mtop / Lambda5, 2.)), 2.) *
                            (std::pow(std::log(std::log(std::pow(mtop / Lambda5, 2.))) - 1. / 2., 2.) +
                             beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

      nf = 6;
      beta0 = 11. - 2. / 3. * nf;
      beta0_sq = beta0 * beta0;
      beta1 = 51. - 19. / 3. * nf;
      beta1_sq = beta1 * beta1;
      beta2 = 2857. - 5033. * nf / 9. + 325. / 27. * nf * nf;

      if (Lambda6 == 0. || std::fabs(1. - alphasMZ_Lambda6 / alphas_MZ) > 1.e-5)
      {
        Lambda_min = 1.e-3;
        Lambda_max = 1.;
        alphas_min = 0.;

        while ((std::fabs(1. - alphas_min / alphas_running) >= 1.e-4) &&
               (std::fabs(1. - Lambda_min / Lambda_max) > 1.e-5))
        {
          real_t log_mtop_Lambda_min = std::log(std::pow(mtop / Lambda_min, 2.));

          alphas_min =
              4. * M_PI / beta0 / log_mtop_Lambda_min *
              (1. - 2. * beta1 / beta0_sq * std::log(log_mtop_Lambda_min) / log_mtop_Lambda_min +
               4. * beta1_sq / std::pow(beta0_sq * log_mtop_Lambda_min, 2.) *
                   (std::pow(std::log(log_mtop_Lambda_min) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

          Lambda_avg = (Lambda_min + Lambda_max) / 2.;
          real_t log_mtop_Lambda_avg = std::log(std::pow(mtop / Lambda_avg, 2.));

          alphas_avg =
              4. * M_PI / beta0 / log_mtop_Lambda_avg *
              (1. - 2. * beta1 / beta0_sq * std::log(log_mtop_Lambda_avg) / log_mtop_Lambda_avg +
               4. * beta1_sq / std::pow(beta0_sq * log_mtop_Lambda_avg, 2.) *
                   (std::pow(std::log(log_mtop_Lambda_avg) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

          if ((alphas_running >= alphas_min) && (alphas_running <= alphas_avg))
            Lambda_max = Lambda_avg;
          else
            Lambda_min = Lambda_avg;
        }

        Lambda6 = Lambda_min;
        alphasMZ_Lambda6 = alphas_MZ;

        if (std::fabs(1. - Lambda_min / Lambda_max) <= 1.e-5)
        {
          Lambda5 = -1.;
          return -1.;
        }
      }
      /*In the end of this if cycle Lmabda6L (local) is or equal to the one in the original class, or
       *equal to something determined in the if. In every case, once the if statement contentent is
       *executed, the two values coincide*/

      alphas_running =
          4. * M_PI / beta0 / std::log(std::pow(Q / Lambda6, 2.)) *
          (1. -
           2. * beta1 / beta0_sq * std::log(std::log(std::pow(Q / Lambda6, 2.))) / std::log(std::pow(Q / Lambda6, 2.)) +
           4. * beta1_sq / std::pow(beta0_sq * std::log(std::pow(Q / Lambda6, 2.)), 2.) *
               (std::pow(std::log(std::log(std::pow(Q / Lambda6, 2.))) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq -
                5. / 4.));

      alphas = alphas_running;
      Q_alphas = Q;
      return alphas_running;
    }
    else
    /* 4 active flavours */
    {
      alphas_running = 4. * M_PI / beta0 / std::log(std::pow(mbot / Lambda5, 2.)) *
                       (1. -
                        2. * beta1 / beta0_sq * std::log(std::log(std::pow(mbot / Lambda5, 2.))) /
                            std::log(std::pow(mbot / Lambda5, 2.)) +
                        4. * beta1_sq / std::pow(beta0_sq * std::log(std::pow(mbot / Lambda5, 2.)), 2.) *
                            (std::pow(std::log(std::log(std::pow(mbot / Lambda5, 2.))) - 1. / 2., 2.) +
                             beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

      nf = 4;
      beta0 = 11. - 2. / 3. * nf;
      beta0_sq = beta0 * beta0;
      beta1 = 51. - 19. / 3. * nf;
      beta1_sq = beta1 * beta1;
      beta2 = 2857. - 5033. * nf / 9. + 325. / 27. * nf * nf;

      if (Lambda4 == 0. || std::fabs(1. - alphasMZ_Lambda4 / alphas_MZ) > 1.e-5)
      {
        Lambda_min = 1.e-3;
        Lambda_max = 1.;
        alphas_min = 0.;

        while ((std::fabs(1. - alphas_min / alphas_running) >= 1.e-4) &&
               (std::fabs(1. - Lambda_min / Lambda_max) > 1.e-5))
        {
          real_t log_mbot_Lambda_min = std::log(std::pow(mbot / Lambda_min, 2.));

          alphas_min =
              4. * M_PI / beta0 / log_mbot_Lambda_min *
              (1. - 2. * beta1 / beta0_sq * std::log(log_mbot_Lambda_min) / log_mbot_Lambda_min +
               4. * beta1_sq / std::pow(beta0_sq * log_mbot_Lambda_min, 2.) *
                   (std::pow(std::log(log_mbot_Lambda_min) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

          Lambda_avg = (Lambda_min + Lambda_max) / 2.;
          real_t log_mbot_Lambda_moy = std::log(std::pow(mbot / Lambda_avg, 2.));

          alphas_avg =
              4. * M_PI / beta0 / log_mbot_Lambda_moy *
              (1. - 2. * beta1 / beta0_sq * std::log(log_mbot_Lambda_moy) / log_mbot_Lambda_moy +
               4. * beta1_sq / std::pow(beta0_sq * log_mbot_Lambda_moy, 2.) *
                   (std::pow(std::log(log_mbot_Lambda_moy) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

          if ((alphas_running >= alphas_min) && (alphas_running <= alphas_avg))
            Lambda_max = Lambda_avg;
          else
            Lambda_min = Lambda_avg;
        }

        Lambda4 = Lambda_min;
        alphasMZ_Lambda4 = alphas_MZ;

        if (std::fabs(1. - Lambda_min / Lambda_max) <= 1.e-5)
        {
          Lambda5 = -1.;
          return -1.;
        }
      }

      if (Q > QuarkMass[CHARM])
      {

        real_t log_Q_lambda4 = std::log(std::pow(Q / Lambda4, 2.));
        alphas_running =
            4. * M_PI / beta0 / log_Q_lambda4 *
            (1. - 2. * beta1 / beta0_sq * std::log(log_Q_lambda4) / log_Q_lambda4 +
             4. * beta1_sq / std::pow(beta0_sq * log_Q_lambda4, 2.) *
                 (std::pow(std::log(log_Q_lambda4) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

        alphas = alphas_running;
        Q_alphas = Q;
        return alphas_running;
      }
      else
      {
        alphas_running = 4. * M_PI / beta0 / std::log(std::pow(QuarkMass[CHARM] / Lambda4, 2.)) *
                         (1. -
                          2. * beta1 / beta0_sq * std::log(std::log(std::pow(QuarkMass[CHARM] / Lambda4, 2.))) /
                              std::log(std::pow(QuarkMass[CHARM] / Lambda4, 2.)) +
                          4. * beta1_sq / std::pow(beta0_sq * std::log(std::pow(QuarkMass[CHARM] / Lambda4, 2.)), 2.) *
                              (std::pow(std::log(std::log(std::pow(QuarkMass[CHARM] / Lambda4, 2.))) - 1. / 2., 2.) +
                               beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

        nf = 3;
        beta0 = 11. - 2. / 3. * nf;
        beta0_sq = beta0 * beta0;
        beta1 = 51. - 19. / 3. * nf;
        beta1_sq = beta1 * beta1;
        beta2 = 2857. - 5033. * nf / 9. + 325. / 27. * nf * nf;

        if (Lambda3 == 0. || std::fabs(1. - alphasMZ_Lambda3 / alphas_MZ) > 1.e-5)
        {
          Lambda_min = 1.e-3;
          Lambda_max = 1.;
          alphas_min = 0.;

          while ((std::fabs(1. - alphas_min / alphas_running) >= 1.e-4) &&
                 (std::fabs(1. - Lambda_min / Lambda_max) > 1.e-5))
          {
            real_t log_mcha_Lambda_min = std::log(std::pow(QuarkMass[CHARM] / Lambda_min, 2.));

            alphas_min =
                4. * M_PI / beta0 / log_mcha_Lambda_min *
                (1. - 2. * beta1 / beta0_sq * std::log(log_mcha_Lambda_min) / log_mcha_Lambda_min +
                 4. * beta1_sq / std::pow(beta0_sq * log_mcha_Lambda_min, 2.) *
                     (std::pow(std::log(log_mcha_Lambda_min) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

            Lambda_avg = (Lambda_min + Lambda_max) / 2.;
            real_t log_mcha_Lambda_moy = std::log(std::pow(QuarkMass[CHARM] / Lambda_avg, 2.));

            alphas_avg =
                4. * M_PI / beta0 / log_mcha_Lambda_moy *
                (1. - 2. * beta1 / beta0_sq * std::log(log_mcha_Lambda_moy) / log_mcha_Lambda_moy +
                 4. * beta1_sq / std::pow(beta0_sq * log_mcha_Lambda_moy, 2.) *
                     (std::pow(std::log(log_mcha_Lambda_moy) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

            if ((alphas_running >= alphas_min) && (alphas_running <= alphas_avg))
              Lambda_max = Lambda_avg;
            else
              Lambda_min = Lambda_avg;
          }
          Lambda3 = Lambda_min;
          alphasMZ_Lambda3 = alphas_MZ;

          if (std::fabs(1. - Lambda_min / Lambda_max) <= 1.e-5)
          {
            Lambda5 = -1.;
            return -1.;
          }
        }

        if (Q < Lambda3)
          Q = Lambda3 + 0.3; /* see discussion in 1604.08082, p125-127 - SuperIso choice: alpha_s(0) ~ 1 */

        real_t log_Q_lambda3 = std::log(std::pow(Q / Lambda3, 2.));
        alphas_running =
            4. * M_PI / beta0 / log_Q_lambda3 *
            (1. - 2. * beta1 / beta0_sq * std::log(log_Q_lambda3) / log_Q_lambda3 +
             4. * beta1_sq / std::pow(beta0_sq * log_Q_lambda3, 2.) *
                 (std::pow(std::log(log_Q_lambda3) - 1. / 2., 2.) + beta2 * beta0 / 8. / beta1_sq - 5. / 4.));

        alphas = alphas_running;
        Q_alphas = Q;
        return alphas_running;
      }
    }
  }

  real_t RunningSM::GetMtopMtop()
  /* computes the top mass mt(mt) */
  {
    real_t alphas_mtop = AlphaStrong(mass_top_pole);

    real_t mtop = mass_top_pole / (1. + alphas_mtop / M_PI *
                                            (4. / 3. + alphas_mtop / M_PI *
                                                           (307. / 32. + M_PI * M_PI / 3. + M_PI * M_PI / 9. * log(2.) -
                                                            1. / 6. * ZETA3 - 71. / 144. * 5.)));

    alphas_mtop = AlphaStrong(mtop, mtop);

    return mass_top_pole / (1. + alphas_mtop / M_PI *
                                     (4. / 3. + alphas_mtop / M_PI *
                                                    (307. / 32. + M_PI * M_PI / 3. + M_PI * M_PI / 9. * std::log(2.) -
                                                     1. / 6. * ZETA3 - 71. / 144. * 5.))); /* 0906.5273 */
  }
  real_t RunningSM::RunQuarkMass(real_t quark_mass, const real_t& Qinit, const real_t& Qfin)
  /* computes the running quark mass at the energy Qfin, from a given running quark mass quark_mass at energy Qinit */
  {
    const real_t mtop = mass_top_pole;
    const real_t mbot = mb_mb;
    real_t alphas_Qinit, alphas_Qfin, running_mass;
    real_t beta0, beta1, beta2, gamma0, gamma1, gamma2;
    int nf;

    alphas_Qinit = AlphaStrong(Qinit);

    real_t R_Qinit, R_Qfin;

    if (Qinit <= mbot)
    /* 4 active flavours at Qinit */
    {
      nf = 4;

      beta0 = 11. - 2. / 3. * nf;
      beta1 = 51. - 19. / 3. * nf;
      beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

      gamma0 = 2.;
      gamma1 = 101. / 12. - 5. / 18. * nf;
      gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

      R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                 1. / 2. *
                     (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                      beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                      beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                     std::pow(alphas_Qinit / M_PI, 2.));

      if (Qfin <= mbot)
      {
        /* if Qinit and Qfin are in the same range, just calculate R_Qfin */
        /* 4 active flavours at Qinit and Qfin */
        nf = 4;

        alphas_Qfin = AlphaStrong(Qfin);

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;
      }
      else if (Qfin <= mtop)
      {
        /* if Qinit and Qfin are NOT in the same range, evolve first the running mass towards the range
         * limit(s), then from the limit(s) towards Qfin */
        /* 4 active flavours at Qinit, 5 active flavours at Qfin */
        nf = 4;

        alphas_Qfin = AlphaStrong(mbot); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;


        alphas_Qinit = alphas_Qfin;

        nf = 5;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                  (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                   1. / 2. *
                       (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                        beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                        beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                       std::pow(alphas_Qinit / M_PI, 2.));


        nf = 5;

        alphas_Qfin = AlphaStrong(Qfin); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * running_mass;
      }
      else
      /* 4 active flavours at Qinit, 6 active flavours at Qfin */
      {
        nf = 4;

        alphas_Qfin = AlphaStrong(mbot); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;

        alphas_Qinit = alphas_Qfin;

        nf = 5;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                  (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                   1. / 2. *
                       (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                        beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                        beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                       std::pow(alphas_Qinit / M_PI, 2.));

        nf = 5;

        alphas_Qfin = AlphaStrong(mtop);

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * running_mass;

        alphas_Qinit = alphas_Qfin;

        nf = 6;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                  (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                   1. / 2. *
                       (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                        beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                        beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                       std::pow(alphas_Qinit / M_PI, 2.));

        nf = 6;

        alphas_Qfin = AlphaStrong(Qfin); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * running_mass;
      }
    }
    else if (Qinit <= mtop)
    /* 5 active flavours at Qinit */
    {
      nf = 5;

      beta0 = 11. - 2. / 3. * nf;
      beta1 = 51. - 19. / 3. * nf;
      beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

      gamma0 = 2.;
      gamma1 = 101. / 12. - 5. / 18. * nf;
      gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

      R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                 1. / 2. *
                     (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                      beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                      beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                     std::pow(alphas_Qinit / M_PI, 2.));

      if (Qfin <= mbot)
      /* 5 active flavours at Qinit, 4 active flavours at Qfin */
      {
        nf = 5;

        alphas_Qfin = AlphaStrong(mbot); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;

        alphas_Qinit = alphas_Qfin;

        nf = 4;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                  (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                   1. / 2. *
                       (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                        beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                        beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                       std::pow(alphas_Qinit / M_PI, 2.));

        nf = 4;

        alphas_Qfin = AlphaStrong(Qfin);

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * running_mass;
      }
      else if (Qfin <= mtop)
      /* 5 active flavours at Qinit and Qfin */
      {
        nf = 5;

        alphas_Qfin = AlphaStrong(Qfin); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;
      }
      else
      /* 5 active flavours at Qinit, 6 active flavours at Qfin */
      {
        nf = 5;

        alphas_Qfin = AlphaStrong(mtop); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;

        alphas_Qinit = alphas_Qfin;

        nf = 6;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                  (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                   1. / 2. *
                       (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                        beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                        beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                       std::pow(alphas_Qinit / M_PI, 2.));

        nf = 6;

        alphas_Qfin = AlphaStrong(Qfin); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * running_mass;
      }
    }
    else
    /* 6 active flavours at Qinit */
    {
      nf = 6;

      beta0 = 11. - 2. / 3. * nf;
      beta1 = 51. - 19. / 3. * nf;
      beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

      gamma0 = 2.;
      gamma1 = 101. / 12. - 5. / 18. * nf;
      gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

      R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                 1. / 2. *
                     (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                      beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                      beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                     std::pow(alphas_Qinit / M_PI, 2.));

      if (Qfin >= mtop)
      /* 6 active flavours at Qinit and Qfin */
      {
        nf = 6;

        alphas_Qfin = AlphaStrong(Qfin); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;
      }
      else if (Qfin >= mbot)
      /* 6 active flavours at Qinit, 5 active flavours at Qfin */
      {
        nf = 6;

        alphas_Qfin = AlphaStrong(mtop); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;

        alphas_Qinit = alphas_Qfin;

        nf = 5;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                  (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                   1. / 2. *
                       (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                        beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                        beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                       std::pow(alphas_Qinit / M_PI, 2.));

        nf = 5;

        alphas_Qfin = AlphaStrong(Qfin);

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * running_mass;
      }
      else
      /* 6 active flavours at Qinit, 4 active flavours at Qfin */
      {
        nf = 6;

        alphas_Qfin = AlphaStrong(mtop);

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * quark_mass;

        alphas_Qinit = alphas_Qfin;


        nf = 5;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                  (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                   1. / 2. *
                       (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                        beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                        beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                       std::pow(alphas_Qinit / M_PI, 2.));

        alphas_Qfin = AlphaStrong(mbot);

        nf = 5;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * running_mass;

        alphas_Qinit = alphas_Qfin;


        nf = 4;

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qinit = std::pow(beta0 / 2. / M_PI * alphas_Qinit, 2. * gamma0 / beta0) *
                  (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qinit / M_PI +
                   1. / 2. *
                       (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                        beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                        beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                       std::pow(alphas_Qinit / M_PI, 2.));


        nf = 4;

        alphas_Qfin = AlphaStrong(Qfin); // Check!

        beta0 = 11. - 2. / 3. * nf;
        beta1 = 51. - 19. / 3. * nf;
        beta2 = 2857. - 5033. / 9. * nf + 325. / 27. * nf * nf;

        gamma0 = 2.;
        gamma1 = 101. / 12. - 5. / 18. * nf;
        gamma2 = 1. / 32. * (1249. - (2216. / 27. + 160. / 3. * ZETA3) * nf - 140. / 81. * nf * nf);

        R_Qfin = std::pow(beta0 / 2. / M_PI * alphas_Qfin, 2. * gamma0 / beta0) *
                 (1. + (2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0) * alphas_Qfin / M_PI +
                  1. / 2. *
                      (std::pow(2. * gamma1 / beta0 - beta1 * gamma0 / beta0 / beta0, 2.) + 2. * gamma2 / beta0 -
                       beta1 * gamma1 / beta0 / beta0 - beta2 * gamma0 / 16. / beta0 / beta0 +
                       beta1 * beta1 * gamma0 / 2. / beta0 / beta0 / beta0) *
                      std::pow(alphas_Qfin / M_PI, 2.));

        running_mass = R_Qfin / R_Qinit * running_mass;
      }
    }

    return running_mass;
  }

  real_t RunningSM::GetMcPole()
  /* computes the c pole mass */
  {
    real_t alphas_mc = AlphaStrong(mc_mc);

    return mc_mc *
           (1 +
            alphas_mc / M_PI *
                (4. / 3. + alphas_mc / M_PI *
                               ((13.4434 - 1.0414 * 3. + 1.0414 * 4. / 3. * ((m_u2GeV + m_d2GeV + m_s2GeV) / mc_mc)))));
  }

  real_t RunningSM::GetMbPole()
  /* computes the b pole mass */
  {
#ifdef DEBUG
    std::cout << "RunningSM::GetMbPole: Called GetMbPole\n";
#endif

    const real_t alphas_mb = AlphaStrong(mb_mb);

#ifdef DEBUG
    std::cout << "RunningSM::GetMbPole: Evaluated AlphaStrong(mb_mb)\n";
#endif

    return mb_mb * (1. + alphas_mb / M_PI *
                             (4. / 3. + alphas_mb / M_PI *
                                            ((13.4434 - 1.0414 * 4. +
                                              1.0414 * 4. / 3. * ((m_u2GeV + m_d2GeV + m_s2GeV + mc_mc) / mb_mb))
                                             /* +alphas_mb/M_PI*(190.595-4.*(26.655-4.*0.6527)) */
                                             )));
  }


  real_t RunningSM::GetMcPole1Loop()
  /* computes the c pole mass at 1 loop */
  {
    const real_t alphas_mc = AlphaStrong(mc_mc);

    return mc_mc * (1 + alphas_mc / M_PI * 4. / 3.);
  }

  real_t RunningSM::GetMbPole1Loop()
  /* computes the b pole mass at 1 loop */
  {
    const real_t alphas_mb = AlphaStrong(mb_mb);

    return mb_mb * (1. + alphas_mb / M_PI * 4. / 3.);
  }

  real_t RunningSM::GetMcPole3Loops()
  /* computes the c pole mass */
  {
    const real_t alphas_mc = AlphaStrong(mc_mc);

    return mc_mc * (1 + alphas_mc / M_PI *
                            (4. / 3. +
                             alphas_mc / M_PI *
                                 ((13.4434 - 1.0414 * 3. + 1.0414 * 4. / 3. * ((m_u2GeV + m_d2GeV + m_s2GeV) / mc_mc)) +
                                  alphas_mc / M_PI * (190.595 - 3. * (26.655 - 3. * 0.6527)))));
  }

  real_t RunningSM::GetMb1S()
  /* computes the 1S b mass */
  {
    const real_t mb = GetMbPole();
    const real_t mu = mb_mb / 2.;
    const real_t as = AlphaStrong(mu);

    return mb * (1. - 2. / 9. * pow(as, 2.));
  }

  real_t RunningSM::GetMcMcFromPole(const real_t& mcpole, const int loop)
  /* computes the running c mass mc(mc) */
  {
    real_t mcmc_min = mcpole / 2.;
    real_t mcmc_max = mcpole;
    real_t mcmc_med, mcpole_med;

    real_t mcmcsav = QuarkMass[CHARM];

    QuarkMass[CHARM] = mcmc_min;

    do
    {
      mcmc_med = (mcmc_min + mcmc_max) / 2.;
      QuarkMass[CHARM] = mcmc_med;

      switch (loop)
      {
      case 2:
        mcpole_med = GetMcPole();
        break;

      case 1:
        mcpole_med = GetMcPole1Loop();
        break;

      default:
        mcpole_med = GetMcPole3Loops();
        break;
      }

      if (mcpole_med < mcpole)
      {
        mcmc_min = mcmc_med;
      }
      else
      {
        mcmc_max = mcmc_med;
      }
    } while (std::fabs(1. - mcpole_med / mcpole) > 1.e-4);

    QuarkMass[CHARM] = mcmcsav;
    QuarkMassQ[CHARM] = QuarkMass[CHARM];

    return mcmc_med;
  }

  real_t RunningSM::GetMbMbFromPole(const real_t& mbpole, const int loop)
  /* computes the running c mass mb(mb) */
  {
    real_t mbmb_min = mbpole / 2.;
    real_t mbmb_max = mbpole;
    real_t mbmb_med, mbpole_med;

    real_t mbmbsav = QuarkMass[BEAUTY];

    if (loop == 2)
    {
      QuarkMass[BEAUTY] = mbmb_min;
      QuarkMass[BEAUTY] = mbmb_max;
    }
    else
    {
      QuarkMass[BEAUTY] = mbmb_min;
      QuarkMass[BEAUTY] = mbmb_max;
    }

    do
    {
      if (loop == 2)
      {
        mbmb_med = (mbmb_min + mbmb_max) / 2.;
        QuarkMass[BEAUTY] = mbmb_med;
        mbpole_med = GetMbPole();
      }
      else
      {
        mbmb_med = (mbmb_min + mbmb_max) / 2.;
        QuarkMass[BEAUTY] = mbmb_med;
        mbpole_med = GetMbPole1Loop();
      }

      if (mbpole_med < mbpole)
      {
        mbmb_min = mbmb_med;
      }
      else
      {
        mbmb_max = mbmb_med;
      }
    } while (std::fabs(1. - mbpole_med / mbpole) > 1.e-4);

    QuarkMass[BEAUTY] = mbmbsav;
    QuarkMassQ[BEAUTY] = QuarkMass[BEAUTY];
    return mbmb_med;
  }

  real_t RunningSM::GetQuarkMass(enum ParticlesList part, const real_t& Qf)
  {
    if (part > RunningSM::TOP || part < RunningSM::UP)
    {
      std::cerr << "Error! GetQuarkMass has been called to return a non-quark mass\n";
      exit(2);
    }
#ifdef DEBUG
    std::cout << "I run all the masses at once, if they're not already ran\n";
#endif

    // Anonymous function to calculate the Higgs corrections
    auto loophiggs = [&](const unsigned short int nf)
    {
      if (!higgsloops)
        return 1.;
      alphas = AlphaStrong(Qf);
      Q_alphas = Qf;
      const real_t a = alphas / M_PI;
      return std::sqrt(1. + 5.67 * a + (35.94 - 1.36 * nf) * a * a + (164.14 - nf * (25.77 - 0.259 * nf)) * a * a * a);
    };

    // if (std::abs(Qf - LastQuarkMassQ[BEAUTY]) > NEGLIGIBLE_ENERGY_DIFFERENCE)
    // {
    //   std::cout << "The higgs corrections for the beauty quark at scale " << Qf << " is: " << loophiggs(5) << "\n";
    // }

    // if (std::abs(Qf - LastQuarkMassQ[TOP]) > NEGLIGIBLE_ENERGY_DIFFERENCE)
    // {
    //   std::cout << "The higgs corrections for the t quark at scale " << Qf << " is: " << loophiggs(6) << "\n";
    // }

    // I run all the masses at once, if they're not already ran
    if (std::fabs(Qf - LastQuarkMassQ[part]) > NEGLIGIBLE_ENERGY_DIFFERENCE)
    {
      // Difference between the current scale Qf and the scale of the initial quark mass
      const real_t Qfmq = std::abs(Qf - QuarkMassQ[part]);

      // Difference between the current scale Qf and the scale of the last quark mass
      const real_t QfQ = std::abs(Qf - LastQuarkMassQ[part]);

      // Taking as input parameters the ones with the least difference
      const real_t mq_in = (Qfmq < QfQ) ? QuarkMass[part] : LastQuarkMass[part];
      const real_t Q_in = (Qfmq < QfQ) ? QuarkMassQ[part] : LastQuarkMassQ[part];

      LastQuarkMass[part] = RunQuarkMass(mq_in, Q_in, Qf);
      LastQuarkMassQ[part] = Qf;
    }
#ifdef DEBUG
    std::cout << "Setting returning values according to lines 272-280 of omega.c\n";
#endif

    // Setting returning values according to lines 272-280 of omega.c
    // the multiplication is the (2.11) from the micromega 1.3 manual
    //  **  it takes into account Higgs corrections
    // the return value has finally always to be less than the pole mass
    if (part == BEAUTY)
    {
      const real_t mb = LastQuarkMass[BEAUTY] * loophiggs(5);
      return std::min(mb, mass_b_pole.get());
    }

    if (part == TOP)
    {
      const real_t mt = LastQuarkMass[TOP] * loophiggs(6);
      return std::min(mt, mass_top_pole.get());
    }

    // If the quark is nor top nor beauty, we return the mass with no Higgs corrections
    return LastQuarkMass[part];
  }

  void RunningSM::HandleParamRunning(Param_t& input, const real_t& Q)
  {
#ifdef DEBUG
    std::cout << "Entering HandleParamRunning\n";
#endif

    // Checking if top quark is massive: if it is not, the user is not interested
    // in the running, so we take g_strong at the mass of the Z
    if (!(mass_top_pole > 0.))
    {
      input.m_t = 0.;
      input.m_b = 0.;
      input.g_s = std::sqrt(4 * M_PI * input.alpha_str_Mz);
      return;
    }

    // Checking if beauty quark is massive: if it is not, the user is not interested
    // in the running of it, so we take g_strong at the mass of the Z, and the
    // mass of the top is left unchanged

    if (!(mb_mb > 0.))
    {
      //     input.m_t=mass_top_pole.get();
      input.m_b = 0.;
      input.g_s = std::sqrt(4 * M_PI * input.alpha_str_Mz);
      return;
    }

#ifdef NO_RUN
    input.m_b = mb_mb;
    input.m_t = mass_top_pole;
    input.g_s = std::sqrt(4 * M_PI * input.alpha_str_Mz);
    return;
#endif

    // Handling running
// Factor needed for the trilinear couplings renormalisation in the MSSM
#ifdef MODEL_MSSM
  #define FACABREF                                                                                                     \
    std::abs(input.m_d* input.V_td* std::conj(input.V_td.get()) +                                                      \
             input.m_s * input.V_ts * std::conj(input.V_ts.get()) +                                                    \
             input.m_b * input.V_tb * std::conj(input.V_tb.get()))
    const real_t Atref = input.A_t; // m_t is the factor
    const real_t Abref = input.A_b; // std::abs(input.m_d*input.V_td*std::conj(input.V_td.get()) +
                                    // input.m_s*input.V_ts*std::conj(input.V_ts.get()) +
                                    // input.m_b*input.V_tb*std::conj(input.V_tb.get()))
    const real_t mtref = input.m_t;
    const real_t facabref = FACABREF;
#endif

#ifdef DEBUG
    std::cout << "Computing g_s\n";
#endif

    // If everything is defined, we

    // 1. Compute g_s

    input.g_s = std::sqrt(4 * M_PI * AlphaStrong(Q, -1));

    // 2. Run quark's masses: from lightest to heaviest
    // COMMENTING ONLY FOR DEBUG
    if (runlightquarks)
    {
      input.m_u = GetQuarkMass(RunningSM::UP, Q);
      input.m_d = GetQuarkMass(RunningSM::DOWN, Q);
      input.m_s = GetQuarkMass(RunningSM::STRANGE, Q);
      input.masses_vector[corr::u] = input.m_u;
      input.masses_vector[corr::d] = input.m_d;
      input.masses_vector[corr::u] = input.m_s;
    }

    if (runcharm)
    {
      input.m_c = GetQuarkMass(RunningSM::CHARM, Q);
      input.masses_vector[corr::c] = input.m_c;
    }

#ifdef DEBUG
    std::cout << "Running heavy quark masses\n";
#endif

    input.m_b = GetQuarkMass(RunningSM::BEAUTY, Q);
    input.m_t = GetQuarkMass(RunningSM::TOP, Q);
    input.masses_vector[corr::b] = input.m_b;
    input.masses_vector[corr::t] = input.m_t;

    input.Running_scale = Q;

#ifdef MODEL_MSSM
    // Updating trilinear couplings struct
    input.A_t = Atref / mtref * input.m_t;
    input.A_b = Abref / facabref * FACABREF;
#endif
#undef FACABREF
    return;
  }


  void RunningSM::TestRunningMasses(Param_t& input)
  {
    const real_t Emin = 0.01, Emax = 10.e+4, deltae = 1.0;
    const int precision = 10;
    std::string pathplots = static_cast<std::string>(OUTPATH) + static_cast<std::string>("running_MARTY");
    char filename[7][1000];
    char tempstring[1000];

    // Creating the directory, if needed, with a system call
    sprintf(tempstring, "[ -d \"%s\" ] || mkdir -p \"%s\"", pathplots.c_str(), pathplots.c_str());
    [[maybe_unused]] int exit_code = std::system(tempstring);

    sprintf(filename[0], "%s/alphaS_Q.dat", pathplots.c_str());

    std::fstream file[7];

    file[0].open(filename[0], std::ios::out);

    for (int i = 1; i < 7; i++)
    {
      sprintf(filename[i], "%s/quark_%d.dat", pathplots.c_str(), i);
      file[i].open(filename[i], std::ios::out);
      file[i] << std::setprecision(precision) << std::scientific;
    }

    csl::InitSanitizer<real_t>* values[7] = {&input.g_s, &input.m_u, &input.m_d, &input.m_s,
                                             &input.m_c, &input.m_b, &input.m_t};

    for (real_t Q = Emin; Q < Emax; Q += deltae)
    {
      this->HandleParamRunning(input, Q);
      file[0] << Q << "\t" << (*values[0]) * (*values[0]) / 4. / M_PI << std::endl;

      for (int i = 1; i < 7; i++)
      {
        file[i] << Q << "\t" << *values[i] << std::endl;
      }
    }
    std::cout << "Closing the files\n";
    for (int i = 0; i < 7; i++)
      file[i].close();
  }
} // end of namespace __SPEC_LIB_NAME__
