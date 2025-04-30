#include "boltzmann.hpp"
#include "macros.hpp"

#define DEBUG
#ifdef DEBUG
  #include <fstream>
  #include <iomanip>
#endif

template <typename T> static inline T SQUARE(const T x) { return x * x; }

namespace __SPEC_LIB_NAME__
{

  void BoltzmannSolver::reset_prevs()
  {
#ifdef DEBUG
    std::cout << "Resetting prev values\n";
#endif
    T_prev_boltzphi = 0.;
    rhophi_prev_boltzphi = 0.;
    T_prev_boltzphi_log = 0.;
    rhophi_prev_boltzphi_log = 0.;
  }

  real_t BoltzmannSolver::Yeq(const real_t& T)
  {
    real_t Yeq = 0.;
    const real_t x = getMassLBSM() / T;

    // Adding BSM particles' contribution
    for (const auto& part : corr::bsm_particles)
    {
      Yeq += corr::part_hel_dof[part] * SQUARE(getMass(part)) * advmath::K2(x * getMass(part) / getMassLBSM());
    }

    Yeq *= 45. * x * x / (4. * std::pow(pi, 4.) * getgeff(T)) / SQUARE(getMassLBSM());

    return Yeq;
  }

  real_t BoltzmannSolver::Sigmatildestar_init(const real_t& T, const real_t& delta)
  {
    /* Computes Sigma tilde* */

    /* Notations:
     *  s is the entropy density:
     *    srad is for radiation
     *    sdark is for entropy injection
     *    stildedark = sdark/srad
     *  rho is the energy density:
     *    rhorad is for radiation
     *    rhodark is for entropy injection
     *    rhotildedark = sdark/srad
     *  n is the number density
     */

    const real_t rho_phi = rhot_phi_Tmax * pi * pi / 15. * std::pow(T, 4.);

    const real_t Y = (1. + delta) * Yeq(T);

    const real_t heff = getheff(T);
    const real_t geff = getgeff(T);

    const real_t s_rad = 2. * SQUARE(pi) / 45. * heff * std::pow(T, 3);
    const real_t rho_rad = SQUARE(pi) / 30. * geff * std::pow(T, 4);

    const real_t s_dark = dark_entropy(T);
    const real_t s_tilde_dark = s_dark / s_rad;

    const real_t n_tilde_dm = Y * s_rad * (1. + s_tilde_dark) * getMassLBSM() / rho_rad;

    const real_t rho_dark = (phi_model.get() != 0 && rho_phi != 0.) ? rho_phi : dark_density(T);
    const real_t rho_tilde_dark = rho_dark / rho_rad;

    const real_t sigma_rad = (phi_model.get() != 0 && rho_phi != 0.) ? Gamma_phi * rho_phi / T : entropy_Sigmarad(T);

    const real_t H_tilde = std::sqrt(1. + rho_tilde_dark + n_tilde_dm);
    // H_tilde = H / std::sqrt( 8 pi / 3 / Mplank^2 ) / rho_rad

    const real_t numerator = 45. * std::sqrt(5) * Mplanck * sigma_rad;
    const real_t denominator = 4 * std::pow(pi, 3.5) * heff * std::sqrt(geff) * std::pow(T, 5) * H_tilde;
    return numerator / denominator;
  }

  real_t BoltzmannSolver::dYeq_dT(const real_t& T) { return (Yeq(T * 1.001) - Yeq(T * 0.999)) / 0.002 / T; }

  void BoltzmannSolver::setTfo(real_t delta)
  {
#ifdef DEBUG
    std::cout << "Called setTfo\n";
#endif
    if (this->full_comput)
      return;

    real_t Tmin, Tmax_local, T, Y, Yphi, rho_phi, srad;
    real_t testTfo = -1.;
    real_t testTfo_prev;
    int compteur = 0;

    Tmax_local = this->Tmax;
#ifdef DEBUG
    std::cout << " Tmax=" << this->Tmax << '\n';
#endif
    Tmin = 1.e-13;

    while (Yeq(Tmin) == 0. && Tmin < Tmax_local)
      Tmin *= 2.;
    if (Tmin >= Tmax_local)
    {
      this->Tfo = -1.;
      return;
    }

    int signSigmatildestar;

    if (1. - Sigmatildestar_init(Tmin, delta) > 0.)
      signSigmatildestar = 1;
    else
      signSigmatildestar = -1;

    compteur = 0;
    testTfo = -1.;

    real_t Tmax_local0, Tmin0;
    Tmax_local0 = Tmax_local;
    Tmin0 = Tmin;
    do
    {
      testTfo_prev = testTfo;
      T = std::exp((std::log(Tmax_local0) + std::log(Tmin0)) / 2.);

      testTfo = 1. - Sigmatildestar_init(T, delta);

      if (std::fabs(1. - testTfo / testTfo_prev) > 1.e-2)
      {
        if (testTfo * signSigmatildestar < 0.)
          Tmax_local0 = T;
        else
          Tmin0 = T;
      }

      compteur++;
    } while (std::fabs(1. - testTfo / testTfo_prev) > 1.e-2 && T > Tmin && compteur < 100);

    if (T < Tmin)
    {
      this->Tfo = -2.;
      return;
    }

    if (compteur < 100)
    {
      Tmax_local = Tmin0;

      if ((1. - Sigmatildestar_init(Tmax_local, delta)) * signSigmatildestar <= 0.)
        Tmax_local -= Tmax_local0 - Tmin0;
    }

    compteur = 0;
    testTfo = -1.;
    /* We want Tfo such that testTfo<0 at higher temperatures and >0 at lower temperatures */
    do
    {
      testTfo_prev = testTfo;
      T = (Tmax_local + Tmin) / 2.;

      Y = (1. + delta) * Yeq(T);

      rho_phi = this->rhot_phi_Tmax * pi * pi / 15. * std::pow(T, 4.);
      srad = 2. * pi * pi / 45. * getheff(T) * std::pow(T, 3.); // MODIFIED: exchanged with down line
      Yphi = rho_phi / srad;                                    // MODIFIED: exchanged with upper line

      testTfo = (1. + delta) * T * dYeq_dT(T) - boltzright(T, getAverageSigmav(T), Y, Yphi);

      if (testTfo < 0.)
        Tmax_local = T;
      else
        Tmin = T;

      compteur++;
    } while (std::fabs(1. - testTfo / testTfo_prev) > 1.e-2 && compteur < 100);

    if (compteur == 100 || T < 3.e-13 || T > 9.e3)
      this->Tfo = -3.;
    else
    {
#ifdef DEBUG
      std::cout << "Setting Tmax = Tfo = " << this->Tfo << std::endl;
#endif
      this->Tfo = this->Tmax = T;
    }

    if (phi_model.get() != 0)
    {
      this->rhot_phi0 = this->rhot_phi_Tmax;
    }

#ifdef DEBUG
    std::cout << "Tfo has been computed. Tfo = " << this->Tfo << std::endl;
    if (Tfo < 0.)
      full_comput = 1;
#endif
    return;
  }

  real_t BoltzmannSolver::boltzright(const real_t& T, const real_t& sv, const real_t& Y, const real_t& Yphi)
  {
    const real_t heffT = getheff(T);
    const real_t geffT = getgeff(T);
    const real_t srad = 2. * pi * pi / 45. * heffT * std::pow(T, 3.);
    const real_t rho_phi = Yphi * srad;
    const real_t rhorad = pi * pi / 30. * geffT * std::pow(T, 4.);
    const real_t darkentropytilde = dark_entropy(T) / srad;
    const real_t dmdensitytilde = Y * srad * (1. + darkentropytilde) * this->relicmass / rhorad;

    //   real_t darkdensitytilde;
    //   if(this->phi_model&&Yphi!=0.) darkdensitytilde=rho_phi/rhorad; else
    //   darkdensitytilde=dark_density(T)/rhorad;
    const real_t darkdensitytilde = (phi_model.get() == 1 && Yphi != 0. ? rho_phi / rhorad : dark_density(T) / rhorad);

    const real_t Htilde =
        std::sqrt(1. + darkdensitytilde + dmdensitytilde); /*Htilde = H / sqrt(8 pi / 3 M_P^2) / rho_rad */

    //   real_t Sigmarad;
    //   if(this->phi_model&&Yphi!=0.) Sigmarad=this->Gamma_phi*rho_phi/T; else Sigmarad=entropy_Sigmarad(T);
    const real_t Sigmarad = (phi_model.get() == 1 && Yphi != 0. ? this->Gamma_phi * rho_phi / T : entropy_Sigmarad(T));

    const real_t Sigmatildestar_local = 45. * std::sqrt(5.) / 4. / std::pow(pi, 3.5) * Mplanck / heffT /
                                        std::sqrt(geffT) / std::pow(T, 5.) / Htilde * Sigmarad;
    const real_t YeqT = Yeq(T);

    //   real_t nontherm;
    //   if(this->phi_model&&Yphi!=0.) nontherm=this->eta_phi*this->Gamma_phi*rho_phi; else nontherm=nonthermal(T);
    const real_t nontherm =
        (phi_model.get() == 1 && Yphi != 0. ? this->eta_phi * this->Gamma_phi * rho_phi : nonthermal(T));

    return std::sqrt(pi / 45.) * Mplanck * T * getsqrtgstar(T) * (1. + darkentropytilde) / Htilde *
           (sv * (Y * Y - YeqT * YeqT) +
            (Y * (Sigmarad + dark_entropy_Sigmad(T)) - nontherm) / std::pow(srad * (1. + darkentropytilde), 2.)) /
           (1. - Sigmatildestar_local);
  }

  /*--------------------------------------------------------------*/

  real_t BoltzmannSolver::boltzright_phi(const real_t& T, const real_t& Y, const real_t& Yphi)
  /* Computes the right-handed term of the decaying scalar_local field Klein-Gordon equation */
  {
    if (phi_model.get() == 0 || Yphi == 0.)
      return 0.;

    // old static definition
    // old static definition

    const real_t heffT = getheff(T);
    const real_t srad = 2. * pi * pi / 45. * heffT * std::pow(T, 3.);
    const real_t rho_phi = Yphi * srad;

    rhot_phi0 = (T_prev_boltzphi > Tinit * K_to_eV) && (T < Tinit * K_to_eV)
                    ? (rhophi_prev_boltzphi + (std::log(this->Tinit * K_to_eV) - std::log(T_prev_boltzphi)) *
                                                  (Yphi * srad - rhophi_prev_boltzphi) /
                                                  (std::log(T) - std::log(T_prev_boltzphi))) /
                          (pi * pi / 15. * std::pow(T, 4.))
                    : rhot_phi0.get();

    T_prev_boltzphi = T;

    //   rho_phi=rhophi_prev_boltzphi=Yphi*srad;
    rhophi_prev_boltzphi = rho_phi;

    const real_t geffT = getgeff(T);
    const real_t rhorad = pi * pi / 30. * geffT * std::pow(T, 4.);
    const real_t dmdensitytilde = Y * srad * this->relicmass / rhorad;
    const real_t darkdensitytilde = rho_phi / rhorad;

    const real_t Htilde =
        std::sqrt(1. + darkdensitytilde + dmdensitytilde); /*Htilde = H / std::sqrt(8 pi / 3 M_P^2) / rho_rad */

    const real_t Sigmarad = Gamma_phi * rho_phi / T;

    const real_t Sigmatildestar_local = 45. * std::sqrt(5.) / 4. / std::pow(pi, 3.5) * Mplanck / heffT /
                                        std::sqrt(geffT) / std::pow(T, 5.) / Htilde * Sigmarad;

    return std::sqrt(pi / 45.) * Mplanck * T * getsqrtgstar(T) / Htilde * Gamma_phi * Yphi * (1. + Yphi / T) / srad /
           (1. - Sigmatildestar_local);
  }

  /*--------------------------------------------------------------*/

  real_t BoltzmannSolver::boltzright_log(const real_t& T, const real_t& sv, const real_t& Y, const real_t& Yphi)
  /* Computes the right-handed term of the Boltzmann equation, with logarithmic derivative */
  {
    const real_t boltz = boltzright(T, sv, Y, Yphi);
    //   if(Y==0.)
    //   {
    //     if(boltz==0.) return boltzright(T,sv,1.e-100,Yphi)/1.e-100;
    //     else return -1.e30;
    //   }
    //   else return boltz/Y;
    return Y == 0. ? (boltz == 0. ? boltzright(T, sv, 1.e-100, Yphi) / 1.e-100 : -1.0e30) : boltz / Y;
  }

  /*--------------------------------------------------------------*/

  real_t BoltzmannSolver::boltzright_phi_log(const real_t& T, const real_t& Y, const real_t& Yphi)
  /* Computes the right-handed term of the decaying scalar_local field Klein-Gordon equation, with logarithmic
   * derivative */
  {
    if (phi_model.get() == 0 || Yphi == 0.)
      return 0.;

    // old static definition
    // old static definition

    const real_t heffT = getheff(T);
    const real_t srad = 2. * pi * pi / 45. * heffT * std::pow(T, 3.);
    const real_t rho_phi = Yphi * srad;

    rhot_phi0 = (T_prev_boltzphi > Tinit * K_to_eV) && (T < Tinit * K_to_eV)
                    ? (rhophi_prev_boltzphi + (std::log(Tinit * K_to_eV) - std::log(T_prev_boltzphi)) *
                                                  (Yphi * srad - rhophi_prev_boltzphi) /
                                                  (std::log(T) - std::log(T_prev_boltzphi))) /
                          (pi * pi / 15. * std::pow(T, 4.))
                    : rhot_phi0.get();

    T_prev_boltzphi = T;
    //   rho_phi=rhophi_prev_boltzphi=Yphi*srad;
    rhophi_prev_boltzphi = rho_phi;

    const real_t geffT = getgeff(T);
    const real_t rhorad = pi * pi / 30. * geffT * std::pow(T, 4.);
    const real_t dmdensitytilde = Y * srad * relicmass / rhorad;
    const real_t darkdensitytilde = rho_phi / rhorad;

    const real_t Htilde =
        std::sqrt(1. + darkdensitytilde + dmdensitytilde); /*Htilde = H / std::sqrt(8 pi / 3 M_P^2) / rho_rad */

    const real_t Sigmarad = Gamma_phi * rho_phi / T;

    const real_t Sigmatildestar_local = 45. * std::sqrt(5.) / 4. / std::pow(pi, 3.5) * Mplanck / heffT /
                                        std::sqrt(geffT) / std::pow(T, 5.) / Htilde * Sigmarad;

    return std::sqrt(pi / 45.) * Mplanck * T * getsqrtgstar(T) / Htilde * Gamma_phi * (1. + Yphi / T) / srad /
           (1. - Sigmatildestar_local);
  }
  /*--------------------------------------------------------------*/

  real_t BoltzmannSolver::relic_density_linsolver()
  /* Computes the relic density, using standard derivatives */
  {
    int counter = 0;
#ifdef DEBUG
    std::cout << "Entering BoltzmannSolver::relic_density_linsolver()\n";
    std::cout << std::setprecision(5) << std::scientific;
    // int scalar_local;
#endif
    real_t sv;
    real_t T, Y, Yb, Ytemp1, Ytemp2, x, Y_eq;
    csl::InitSanitizer<real_t> k1, k2, k3, k4, k5, k6, k7, k8;
    real_t Yphi, Yphib, Yphitemp1, Yphitemp2;
    csl::InitSanitizer<real_t> krhophi1, krhophi2, krhophi3, krhophi4, krhophi5, krhophi6, krhophi7, krhophi8;
    real_t rho_phi = 0.;
    int testY = 1, testYphi = 1, testdx = 1;
    real_t testdiv;
    //   testY=testYphi=testdx=1;

    if (phi_model.get() == 0)
      this->rhot_phi_Tmax = this->eta_phi = 0.;

    // Relicmass assignment
    constexpr const real_t Tfin = 2.3482260e-13;
    real_t delta = 1.5;
    setTfo(delta);

    if (this->Tmax < Tfin)
      return -1.;

    rho_phi = this->rhot_phi_Tmax * pi * pi / 15. * std::pow(this->Tmax, 4.);

#ifdef DEBUG
    if (this->full_comput == 0)
      std::cout << "Tfo=" << this->Tfo << '\n';
    //   FILE *sig,*phi0;
    //   sig=fopen("sigmav.out","w");
    std::ofstream fsig("sigmav.out"), fphi0;
    if (!fsig)
    {
      std::cerr << "Impossible to open sigmav.out\n";
      exit(1);
    }
    if (phi_model != 0)
    {
      fphi0.open("rhophi0.out");
      if (!fphi0)
      {
        std::cerr << "Impossible to open rhophi0.out\n";
        exit(1);
      }
    }
    fsig << std::setprecision(5) << std::scientific;
#endif

    //   if(this->full_comput==0) adapt_error=adapt_error_phi=1.e-3; else adapt_error=adapt_error_phi=1.e-2;
    real_t adapt_error_phi = (this->full_comput == 0 ? 1.e-3 : 1.e-2);
    real_t adapt_error = adapt_error_phi;

    T = this->Tmax;
    x = std::log(T);

    //   if(phi_model.get()==0) dx=(std::log(Tfin)-std::log(this->Tmax))*adapt_error*10.;
    //   else dx=(std::log(Tfin)-std::log(this->Tmax))*std::min(adapt_error,adapt_error_phi)*10.;
    real_t dx =
        (phi_model.get() == 0 ? (std::log(Tfin) - std::log(this->Tmax)) * adapt_error * 10.
                              : (std::log(Tfin) - std::log(this->Tmax)) * std::min(adapt_error, adapt_error_phi) * 10.);
#ifdef DEBUG
    std::cout << "Tfin=" << Tfin << ", Tmax=" << Tmax.get() << ", adapt_error=" << adapt_error << "dx=" << dx << '\n';
#endif

    if (this->full_comput)
      Y = Yeq(T);
    else
      Y = Yeq(T) * (1. + delta);
    Yphi = rho_phi / (2. * pi * pi / 45. * getheff(T) * std::pow(T, 3.));

    if (Yphi < 1.e-100)
      Yphi = 0.;
    if (Y < 1.e-100)
      Y = 0.;

    do
    {
      ++counter;
      if (this->full_comput == 0)
      {
        if (counter > 100000)
          return -1.;
      }
      else
      {
        if (counter > 1000000)
          return -1.;
      }

      if (Y < 1.e-30)
        return 0.;

#ifdef DEBUG
      // if (counter % 10000 == 0)
      //       printf("T=%.5e    Y=%.5e  Yeq=%.5e  Yphi=%.5e
      //       dx=%.3e\n",std::exp(x),Y,Yeq(std::exp(x)),Yphi,dx);
      std::cout << "T=" << std::scientific << std::setprecision(5) << std::exp(x) << "  Y=" << Y
                << "  Yeq=" << Yeq(std::exp(x)) << "  Yphi=" << Yphi << "  dx=" << std::scientific
                << std::setprecision(3) << dx << std::endl;
#endif

      if (this->full_comput == 0)
        if (testdx == 0)
        {
          if ((adapt_error > 0.25) || (adapt_error > 0.5))
          {
#ifdef DEBUG
            //         printf("Warning: adaptative step too small and accepted error too large,
            //         abandoning...\n");
            std::cerr << "Warning: adaptative step too small and accepted error too large, abandoning...\n";
#endif
            return -1.;
          }
          else
          {
            if (testY)
            {
#ifdef DEBUG
              //           printf("Warning: adaptative step too small (%.3e), increasing the accepted
              //           error (%.3e).\n",std::fabs(dx),adapt_error);
              std::cout << "Warning: adaptive step too small (" << std::scientific << std::setprecision(3)
                        << std::fabs(dx) << "), increasing the accepted error (" << std::scientific
                        << std::setprecision(3) << adapt_error << ")." << std::endl;
#endif
              adapt_error *= 2.;
            }
            if (testYphi)
            {
#ifdef DEBUG
              //           printf("Warning: adaptative step too small (%.3e) for phi, increasing the
              //           accepted error (%.3e).\n",std::fabs(dx),adapt_error_phi);
              std::cerr << "Warning: adaptive step too small (" << std::scientific << std::setprecision(3)
                        << std::fabs(dx) << ") for phi, increasing the accepted error (" << std::scientific
                        << std::setprecision(3) << adapt_error_phi << ")." << std::endl;
#endif
              adapt_error_phi *= 2.;
            }
            dx *= 10.;
          }
        }

      if (std::isnan(Y) || std::isinf(Y) || std::isnan(Yphi) || std::isinf(Yphi))
      {
        if (this->full_comput == 0)
        {
          delta *= 1.1;
          if (delta > 10.)
            return -1.;

          setTfo(delta);

          if (this->Tmax < Tfin)
            return -1.;

          rho_phi = this->rhot_phi_Tmax * pi * pi / 15. * std::pow(this->Tmax, 4.);

#ifdef DEBUG
          //         fclose(sig);
          //         printf("T_fo=%.5e\n",this->Tfo);
          //         sig=fopen("sigmav.out","w");

          //         fsig.close();
          //         fsig.open("sigmav.out");
          //         fsig << std::setprecision(5) << std::scientific;
          std::cout << this->Tfo;
#endif

          adapt_error = adapt_error_phi = 1.e-3;

          T = this->Tmax;
          x = std::log(T);
          dx = (std::log(Tfin) - std::log(this->Tmax)) * adapt_error * 10.;

          Y = Yeq(T) * (1. + delta);
          Yphi = rho_phi / (2. * pi * pi / 45. * getheff(T) * std::pow(T, 3.));
        }
        else
          return -1.;
      }

#ifdef DEBUG
      //     fprintf(sig,"%.5e\t %.5e\n",this->relicmass/T, sv);
      fsig << this->relicmass / T << '\t' << sv << '\n';
      if (this->phi_model != 0)
      {
        //       fprintf(phi0,"%.5e\t%.5e\t%.5e\t%.5e\t%.5e \t%.5e\n",T,this->relicmass/T,
        //               Yphi*(2.*pi*pi/45.*getheff(T)*std::pow(T,3.)),Y*(2.*pi*pi/45.*heff(T)*std::pow(T,3.)*this->relicmass),(pi*pi/30.*getgeff(T)*std::pow(T,4.)),2.*pi*pi/45.*heff(T)*std::pow(T,3.));
        fphi0 << std::scientific << std::setprecision(5) << T << "\t" << relicmass / T << "\t"
              << Yphi * (2. * pi * pi / 45. * getheff(T) * std::pow(T, 3.)) << "\t"
              << Y * (2. * pi * pi / 45. * heff(T) * std::pow(T, 3.) * relicmass) << "\t"
              << (pi * pi / 30. * getgeff(T) * std::pow(T, 4.)) << "\t"
              << 2. * pi * pi / 45. * heff(T) * std::pow(T, 3.) << "\n";
      }
#endif

      if (Yphi < 1.e-100)
        Yphi = 0.;
      if (Y < 1.e-100)
        Y = 0.;

      testY = testYphi = testdx = 1;
      /* testdx=1 (OK) if dx/x is not too small
         testY=1 if precision on Y is bad, i.e. dx too large
         testY=1 if precision on Yphi is bad, i.e. dx too large
      */

      while ((testY || testYphi) && testdx)
      /* if testY=1 or testphi=1, the precision is too bad, so the next step is recomputed with a smaller dx
         if testdx=0, the step is too small, the loop has to be broken, and the accepted precision decreased */
      {
        T = std::exp(x);
        sv = getAverageSigmav(T);
        k1 = boltzright(T, sv, Y, Yphi);
        krhophi1 = boltzright_phi(T, Y, Yphi);

        if (k1 < -1.e5 * Y)
          return 0.;

        if (Y + k1 * dx / 2. < 0.)
          dx *= std::fabs(Y / k1) / 2.2;

        testdiv = 0;
        do
        {
#ifdef DEBUG
          if (testdiv)
            std::cout << "Avoiding negative densities...\n";
#endif
          testdiv = 0;
          T = std::exp(x + dx / 2.);
          sv = getAverageSigmav(T);
          k2 = boltzright(T, sv, Y + k1 * dx / 2., Yphi + krhophi1 * dx / 2.);
          krhophi2 = boltzright_phi(T, Y + k1 * dx / 2., Yphi + krhophi1 * dx / 2.);

          if (Y + k2 * dx / 2. < 0.)
          {
            testdiv = 1;
            dx *= std::fabs(Y / k2) / 2.2;
          }

          if (!testdiv)
          {
            k3 = boltzright(T, sv, Y + k2 * dx / 2., Yphi + krhophi2 * dx / 2.);
            krhophi3 = boltzright_phi(T, Y + k2 * dx / 2., Yphi + krhophi2 * dx / 2.);

            if (Y + k3 * dx < 0.)
            {
              testdiv = 1;
              dx *= std::fabs(Y / k3) / 1.1;
            }
          }

          if (!testdiv)
          {
            T = std::exp(x + dx);
            sv = getAverageSigmav(T);
            k4 = boltzright(T, sv, Y + k3 * dx, Yphi + krhophi3 * dx);
            krhophi4 = boltzright_phi(T, Y + k3 * dx, Yphi + krhophi3 * dx);

            if (Y + dx / 6. * (k1 + 2. * k2 + 2. * k3 + k4) < 0.)
            {
              testdiv = 1;
              dx *= std::fabs(Y / (k1 + 2. * k2 + 2. * k3 + k4)) / 6.6;
            }
          }

          testdx = std::fabs(dx / x) > 1.e-8;
        } while (testdiv && testdx);

        Ytemp1 = Y + dx / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        Yphitemp1 = Yphi + dx / 6. * (krhophi1 + 2. * krhophi2 + 2. * krhophi3 + krhophi4);

        /*-------------*/

        dx = dx / 2.;

        T = std::exp(x + dx / 2.);
        sv = getAverageSigmav(T);
        k2 = boltzright(T, sv, Y + k1 * dx / 2., Yphi + krhophi1 * dx / 2.);
        krhophi2 = boltzright_phi(T, Y + k1 * dx / 2., Yphi + krhophi1 * dx / 2.);
        k3 = boltzright(T, sv, Y + k2 * dx / 2., Yphi + krhophi2 * dx / 2.);
        krhophi3 = boltzright_phi(T, Y + k2 * dx / 2., Yphi + krhophi2 * dx / 2.);

        T = std::exp(x + dx);
        sv = getAverageSigmav(T);
        k4 = boltzright(T, sv, Y + k3 * dx, Yphi + krhophi3 * dx);
        krhophi4 = boltzright_phi(T, Y + k3 * dx, Yphi + krhophi3 * dx);

        Yb = Y + dx / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        Yphib = Yphi + dx / 6. * (krhophi1 + 2. * krhophi2 + 2. * krhophi3 + krhophi4);

        k5 = boltzright(T, sv, Yb, Yphib);
        krhophi5 = boltzright_phi(T, Yb, Yphib);

        T = std::exp(x + 1.5 * dx);
        sv = getAverageSigmav(T);
        k6 = boltzright(T, sv, Yb + k5 * dx / 2., Yphib + krhophi5 * dx / 2.);
        krhophi6 = boltzright_phi(T, Yb + k5 * dx / 2., Yphib + krhophi5 * dx / 2.);
        k7 = boltzright(T, sv, Yb + k6 * dx / 2., Yphib + krhophi6 * dx / 2.);
        krhophi7 = boltzright_phi(T, Yb + k6 * dx / 2., Yphib + krhophi6 * dx / 2.);

        T = std::exp(x + 2. * dx);
        sv = getAverageSigmav(T);
        k8 = boltzright(T, sv, Yb + k7 * dx, Yphib + krhophi7 * dx);
        krhophi8 = boltzright_phi(T, Yb + k7 * dx, Yphib + krhophi7 * dx);

        Ytemp2 = Yb + dx / 6. * (k5 + 2. * k6 + 2. * k7 + k8);
        Yphitemp2 = Yphib + dx / 6. * (krhophi5 + 2. * krhophi6 + 2. * krhophi7 + krhophi8);

        testY = (std::isnan(Ytemp1) || (Ytemp1 < 0.) || std::isnan(Ytemp2) || (Ytemp2 < 0.) ||
                 (std::fabs((Ytemp2 - Ytemp1) / (Ytemp2)) > adapt_error));

        testYphi = (phi_model.get()) * (std::isnan(Yphitemp1) || (Yphitemp1 < 0.) || std::isnan(Yphitemp2) ||
                                        (Yphitemp2 < 0.) ||
                                        (std::fabs((Yphitemp2 - Yphitemp1) / (Yphitemp2)) > adapt_error_phi));

        if (this->full_comput)
        {
          Y_eq = Yeq(std::exp(x + 2. * dx));
          if (Y_eq > 0.)
            if (((!std::isnan(Ytemp2) && !std::isinf(Ytemp2)) &&
                 ((Ytemp2 < Y_eq) || (std::fabs(1. - Ytemp2 / Y_eq) < 0.01))) ||
                ((std::isnan(Ytemp2) || std::isinf(Ytemp2)) &&
                 ((Y + 2. * dx * k1 < Y_eq) || (std::fabs(1. - (Y + 2. * dx * k1) / Y_eq) < 0.001))))
            {
              Ytemp2 = Ytemp1 = Y_eq;
              testY = 0;
            }
        }

        testdx = std::fabs(dx / x) > 1.e-8;

        if (testY || testYphi)
        {
          if ((phi_model.get()) && (Yphitemp2 - Yphitemp1 != 0.))
            dx *= 2. * 0.9 *
                  std::min(std::max(std::min(std::fabs(Ytemp2 * adapt_error / (Ytemp2 - Ytemp1)) *
                                                 !std::isnan(Ytemp2 / (Ytemp2 - Ytemp1)),
                                             std::fabs(Yphitemp2 * adapt_error_phi / (Yphitemp2 - Yphitemp1)) *
                                                 !std::isnan(Yphitemp2 / (Yphitemp2 - Yphitemp1))),
                                    0.3),
                           1.);
          else
            dx *= 2. * 0.9 *
                  std::min(std::max(std::fabs(Ytemp2 * adapt_error / (Ytemp2 - Ytemp1)) *
                                        !std::isnan(Ytemp2 / (Ytemp2 - Ytemp1)),
                                    0.3),
                           1.);
        }
      }

      if (this->full_comput)
      {
        if ((Y < (1. + delta) * Yeq(std::exp(x))) && (Yb > (1. + delta) * Yeq(std::exp(x + dx))))
        {
          const real_t y1 = Y - (1. + delta) * Yeq(std::exp(x));
          const real_t x1 = x;
          const real_t y2 = Yb - Yeq(std::exp(x + dx));
          const real_t x2 = x + dx;
          const real_t y0 = 0.;
          const real_t x0 = (y0 - y1) * (x2 - x1) / (y2 - y1) + x1;

          this->Tfo = std::exp(x0);
        }
        else if ((Yb < (1. + delta) * Yeq(std::exp(x + dx))) && (Ytemp2 > (1. + delta) * Yeq(std::exp(x + 2. * dx))))
        {
          const real_t y1 = Yb - Yeq(std::exp(x + dx));
          const real_t x1 = x + dx;
          const real_t y2 = Ytemp2 - Yeq(std::exp(x + 2. * dx));
          const real_t x2 = x + 2. * dx;
          const real_t y0 = 0.;
          const real_t x0 = (y0 - y1) * (x2 - x1) / (y2 - y1) + x1;

          this->Tfo = std::exp(x0);
        }
      }

      Y = Ytemp2;
      Yphi = Yphitemp2;
      x = x + 2. * dx;

      if (Yphi < 1e-10 * Y)
        Yphi = 0.;

      if ((phi_model.get()) && (Yphitemp2 - Yphitemp1 != 0.))
        dx *= 2. * 0.9 *
              std::max(1., std::min(std::min(std::fabs(Ytemp2 * adapt_error / (Ytemp2 - Ytemp1)) *
                                                 !std::isnan(Ytemp2 / (Ytemp2 - Ytemp1)),
                                             std::fabs(Yphitemp2 * adapt_error_phi / (Yphitemp2 - Yphitemp1)) *
                                                 !std::isnan(Yphitemp2 / (Yphitemp2 - Yphitemp1))),
                                    2.));
      else
        dx *= 2. * 0.9 *
              std::max(1., std::min(std::fabs(Ytemp2 * adapt_error / (Ytemp2 - Ytemp1)) *
                                        !std::isnan(Ytemp2 / (Ytemp2 - Ytemp1)),
                                    2.));

      if (dx > x - std::log(Tfin))
        dx = x - std::log(Tfin);
    } while (T > Tfin);

#ifdef DEBUG
    //   fclose(sig);
    //     fclose(phi0);
    //   printf("%d iterations where needed for the integration of the Boltzmann equation.\n\n",counter);

    fsig.close();
    if (fphi0)
      fphi0.close();
    std::cout << counter << " iterations were needed for the integration of the Boltzmann equation.\n\n";
#endif

    if ((this->mgravitino != 0.) && (this->mgravitino < this->relicmass))
      return 2.755e8 * this->mgravitino * Y;

    return 2.472e8 * this->relicmass * Y;
  }

  /*--------------------------------------------------------------*/

  real_t BoltzmannSolver::relic_density_logsolver()
  /* Computes the relic density, using logarithmic derivatives */
  {
    int counter = 0;
#ifdef DEBUG
    std::cout << "Entering BoltzmannSolver::relic_density_logsolver()\n";
    std::cout << std::scientific << std::setprecision(5);
    // real_t boltz;
    // real_t Yb;
    // real_t Yphib;
    // real_t phi,phib,phitemp1,phitemp2,kphi1,kphi2,kphi3,kphi4,kphi5,kphi6,kphi7,kphi8;
    // real_t logphi,logphib,logphitemp1,logphitemp2;
    // int scalar_local;
    // real_t Ytemp1,Ytemp2;
#endif
    real_t sv, Y_eq;
    real_t k1, k2, k3, k4, k5, k6, k7, k8;
    real_t Yphitemp1, Yphitemp2, krhophi1, krhophi2, krhophi3, krhophi4, krhophi5, krhophi6, krhophi7, krhophi8;
    real_t logYb, logYtemp1, logYtemp2, logYphib, logYphitemp1, logYphitemp2;

    int testY = 1, testYphi = 1, testdx = 1; // Can be converted into booleans

    if (phi_model.get() == 0)
      this->rhot_phi_Tmax = this->eta_phi = 0.;

    constexpr const real_t Tfin = 2.3482260e-13;

    real_t delta = 1.5;
    setTfo(delta);

    if (this->Tmax < Tfin)
      return -1.;

    real_t rho_phi = (phi_model.get() == 1 ? this->rhot_phi_Tmax * pi * pi / 15. * std::pow(this->Tmax, 4.) : 0.);
#ifdef DEBUG
    std::cout << "rho_phi = " << rho_phi << '\n';
    if (this->full_comput == 0)
      std::cout << "T_fo=" << this->Tfo << '\n';

    std::cout << "T_max=" << this->Tmax << '\n';
    //   FILE *sig,*phi0;
    //   sig=fopen("sigmav.out","w");
    //   if(this->phi_model !=0 ) phi0=fopen("rhophi.out","w");
    std::ofstream fsig("sigmav.out"), fphi0;
    if (!fsig)
    {
      std::cerr << "Impossible to open sigmav.out\n";
      exit(1);
    }
    if (phi_model.get() != 0)
    {
      fphi0.open("rhophi0.out");
      if (!fphi0)
      {
        std::cerr << "Impossible to open rhophi0.out\n";
        exit(1);
      }
    }
    fsig << std::setprecision(5) << std::scientific;
    std::ofstream fYt("Yt.dat");
    if (!fYt)
    {
      std::cerr << "Impossible to open Yt.dat\n";
      exit(1);
    }
    fYt << std::setprecision(5) << std::scientific;
#endif

    real_t adapt_error = 0.001;
    real_t adapt_error_phi = 1.e-4;

#ifdef DEBUG
    std::cout << "Tfin=" << Tfin << ", Tmax=" << Tmax.get() << ", adapt_error=" << adapt_error << '\n';
#endif

    real_t T = this->Tmax;
    real_t x = std::log(T);
    real_t dx = phi_model.get() == 0
                    ? (std::log(Tfin) - std::log(this->Tmax)) * adapt_error * 10.
                    : (std::log(Tfin) - std::log(this->Tmax)) * std::min(adapt_error, adapt_error_phi) * 10.;

    //   if(phi_model.get()==0) dx=(std::log(Tfin)-std::log(this->Tmax))*adapt_error*10.;
    //   else dx=(std::log(Tfin)-std::log(this->Tmax))*std::min(adapt_error,adapt_error_phi)*10.;

    real_t YeqT = Yeq(T);
#ifdef DEBUG
    std::cout << "dx=" << dx << " Y_eq=" << YeqT << '\n';
#endif

    real_t Y = (this->full_comput) ? YeqT : YeqT * (1. + delta);
    real_t Yphi =
        (phi_model.get() == 1 && rho_phi != 0.) ? rho_phi / (2. * pi * pi / 45. * getheff(T) * std::pow(T, 3.)) : 0.;
    real_t logY = (Y == 0.) ? -1.e100 : std::log(Y);
    real_t logYphi = (Yphi == 0.) ? -1.e100 : std::log(Yphi);

//   if(this->full_comput) Y=Yeq(T); else Y=Yeq(T)*(1.+delta);
//   if(phi_model.get() !=0 && rho_phi!=0.) Yphi=rho_phi/(2.*pi*pi/45.*getheff(T)*std::pow(T,3.)); else Yphi=0.;
//   if(Y==0.) logY=-1.e100; else logY=std::log(Y);
//   if(Yphi==0.) logYphi=-1.e100; else logYphi=std::log(Yphi);
#ifdef DEBUG
    std::cout << "Y=" << Y << ", Yphi=" << Yphi << ", logY=" << logY << ", logYphi=" << logYphi << '\n';
#endif

    do
    {
      constexpr const real_t Y_inf_limit(1.e-30);
      ++counter;
      if (this->full_comput == 0)
      {
        if (counter > 10000)
          return -1.;
      }
      else
      {
        if (counter > 1000000)
          return -1.;
      }

      if (Y < Y_inf_limit)
      {
#ifdef DEBUG
        std::cout << "Y < " << Y_inf_limit << ", returning 0.\n";
#endif
        return 0.;
      }

#ifdef DEBUG
      std::cout << "counter=" << counter << ", T=" << std::exp(x) << ", Yeq=" << Yeq(std::exp(x))
                << ", Tfo=" << this->Tfo << ", Yphi=" << Yphi << ", dx=" << dx << ", s=" << rho_phi / Yphi << '\n';
      //     fprintf(fptr, "%.4e\t%.4e\n", 1./T, Y);
      fYt << 1. / T << '\t' << Y << '\n';
#endif

      if (this->full_comput == 0)
        if (testdx == 0)
        {
          if ((adapt_error > 0.25) || (adapt_error > 0.5))
          {
#ifdef DEBUG
            std::cout << "Warning: adaptive step too small and accepted error too large, abandoning...\n";
#endif
            return -1.;
          }
          else
          {
            if (testY)
            {
#ifdef DEBUG
              //           printf("Warning: adaptive step too small (%.3e), increasing the accepted error
              //           (%.3e).\n",std::fabs(dx),adapt_error);
              std::cout << "Warning: adaptive step too small (" << std::scientific << std::setprecision(3)
                        << std::fabs(dx) << "), increasing the accepted error (" << std::scientific
                        << std::setprecision(3) << adapt_error << ")." << std::endl;
#endif
              adapt_error *= 2.;
            }
            if (testYphi)
            {
#ifdef DEBUG
              //           printf("Warning: adaptive step too small (%.3e) for phi, increasing the
              //           accepted error (%.3e).\n",std::fabs(dx),adapt_error_phi);
              std::cout << "Warning: adaptive step too small (" << std::scientific << std::setprecision(3)
                        << std::fabs(dx) << ") for phi, increasing the accepted error (" << std::scientific
                        << std::setprecision(3) << adapt_error_phi << ")." << std::endl;
#endif
              adapt_error_phi *= 2.;
            }
            dx *= 10.;
          }
        }

      if (std::isnan(logY) || std::isinf(logY) || std::isnan(logYphi) || std::isinf(logYphi))
      {
        if (this->full_comput == 0)
        {
          delta *= 1.1;
          if (delta > 10.)
            return -1.;

          setTfo(delta);

          if (this->Tmax < Tfin)
            return -1.;

          //         if(this->phi_model) rho_phi=this->rhot_phi_Tmax*pi*pi/15.*std::pow(this->Tmax,4.); else
          //         rho_phi=0.;
          rho_phi = (phi_model.get() == 1) ? this->rhot_phi_Tmax * pi * pi / 15. * std::pow(this->Tmax, 4.) : 0.;

#ifdef DEBUG
          //         fclose(sig);
          //         sig=fopen("sigmav.out","w");]
          //         fsig.close();
          //         fsig.open("sigmav.out");
          std::cout << this->Tfo;
#endif

          adapt_error = adapt_error_phi = 1.e-3;
          T = this->Tmax;
          x = std::log(T);
          dx = (std::log(Tfin) - std::log(this->Tmax)) * adapt_error * 10.;

          Y = Yeq(T) * (1. + delta);
          //         if(phi_model.get() !=0 && rho_phi!=0.)
          //         Yphi=rho_phi/(2.*pi*pi/45.*getheff(T)*std::pow(T,3.)); else Yphi=0.;
          Yphi = (phi_model.get() == 1 && rho_phi != 0.) ? rho_phi / (2. * pi * pi / 45. * getheff(T) * std::pow(T, 3.))
                                                         : 0.;
        }
        else
          return -1.;
      }

#ifdef DEBUG
      //     fprintf(sig,"%.5e\t %.5e\n",this->relicmass/T, sv);
      //     if(this->phi_model) fprintf(phi0,"%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t
      //     %.10e\n",T,this->relicmass/T,Yphi*(2.*pi*pi/45.*getheff(T)*std::pow(T,3.)),Y*(2.*pi*pi/45.*heff(T)*std::pow(T,3.)),(pi*pi/30.*getgeff(T)*std::pow(T,4.)),2.*pi*pi/45.*heff(T)*std::pow(T,3.));
      fsig << std::setprecision(5) << std::scientific << this->relicmass / T << '\t' << sv << '\n';
      if (phi_model.get() != 0)
        fphi0 << std::scientific << std::setprecision(10) << T << "\t" << relicmass / T << "\t"
              << Yphi * (2. * pi * pi / 45. * getheff(T) * std::pow(T, 3.)) << "\t"
              << Y * (2. * pi * pi / 45. * heff(T) * std::pow(T, 3.) * relicmass) << "\t"
              << (pi * pi / 30. * getgeff(T) * std::pow(T, 4.)) << "\t"
              << 2. * pi * pi / 45. * heff(T) * std::pow(T, 3.) << "\n";
#endif

      //     if(Y==0.) logY=-1.e100; else logY=std::log(Y);
      //     if(Yphi==0.) logYphi=-1.e100; else logYphi=std::log(Yphi);
      logY = (Y == 0.) ? -1.e100 : std::log(Y);
      logYphi = (Yphi == 0.) ? -1.e100 : std::log(Yphi);

      testY = testYphi = testdx = 1;
      /* testdx=1 (OK) if dx/x is not too small
         testY=1 if precision on Y is bad, i.e. dx too large
         testY=1 if precision on Yphi is bad, i.e. dx too large
      */

      while ((testY || testYphi) && testdx)
      /* if testY=1 or testphi=1, the precision is too bad, so the next step is recomputed with a smaller dx
         if testdx=0, the step is too small, the loop has to be broken, and the accepted precision decreased */
      {
        T = std::exp(x);
        krhophi1 = boltzright_phi_log(T, std::exp(logY), std::exp(logYphi));
        sv = getAverageSigmav(T);
        k1 = boltzright_log(T, sv, std::exp(logY), std::exp(logYphi));

        T = std::exp(x + dx / 2.);
        krhophi2 = boltzright_phi_log(T, std::exp(logY + k1 * dx / 2.), std::exp(logYphi + krhophi1 * dx / 2.));
        sv = getAverageSigmav(T);
        k2 = boltzright_log(T, sv, std::exp(logY + k1 * dx / 2.), std::exp(logYphi + krhophi1 * dx / 2.));

        krhophi3 = boltzright_phi_log(T, std::exp(logY + k2 * dx / 2.), std::exp(logYphi + krhophi2 * dx / 2.));
        k3 = boltzright_log(T, sv, std::exp(logY + k2 * dx / 2.), std::exp(logYphi + krhophi2 * dx / 2.));

        T = std::exp(x + dx);
        krhophi4 = boltzright_phi_log(T, std::exp(logY + k3 * dx), std::exp(logYphi + krhophi3 * dx));
        sv = getAverageSigmav(T);
        k4 = boltzright_log(T, sv, std::exp(logY + k3 * dx), std::exp(logYphi + krhophi3 * dx));

        logYtemp1 = logY + dx / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        logYphitemp1 = logYphi + dx / 6. * (krhophi1 + 2. * krhophi2 + 2. * krhophi3 + krhophi4);

        // Ytemp1=std::exp(logYtemp1);
        Yphitemp1 = std::exp(logYphitemp1);

        /*-------------*/

        dx = dx / 2.;

        T = std::exp(x + dx / 2.);
        krhophi2 = boltzright_phi_log(T, std::exp(logY + k1 * dx / 2.), std::exp(logYphi + krhophi1 * dx / 2.));
        sv = getAverageSigmav(T);
        k2 = boltzright_log(T, sv, std::exp(logY + k1 * dx / 2.), std::exp(logYphi + krhophi1 * dx / 2.));
        krhophi3 = boltzright_phi_log(T, std::exp(logY + k2 * dx / 2.), std::exp(logYphi + krhophi2 * dx / 2.));
        k3 = boltzright_log(T, sv, std::exp(logY + k2 * dx / 2.), std::exp(logYphi + krhophi2 * dx / 2.));

        T = std::exp(x + dx);
        krhophi4 = boltzright_phi_log(T, std::exp(logY + k3 * dx), std::exp(logYphi + krhophi3 * dx));
        sv = getAverageSigmav(T);
        k4 = boltzright_log(T, sv, std::exp(logY + k3 * dx), std::exp(logYphi + krhophi3 * dx));

        logYb = logY + dx / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        logYphib = logYphi + dx / 6. * (krhophi1 + 2. * krhophi2 + 2. * krhophi3 + krhophi4);

        krhophi5 = boltzright_phi_log(T, std::exp(logYb), std::exp(logYphib));
        k5 = boltzright_log(T, sv, std::exp(logYb), std::exp(logYphib));

        T = std::exp(x + 1.5 * dx);
        krhophi6 = boltzright_phi_log(T, std::exp(logYb + k5 * dx / 2.), std::exp(logYphib + krhophi5 * dx / 2.));
        sv = getAverageSigmav(T);
        k6 = boltzright_log(T, sv, std::exp(logYb + k5 * dx / 2.), std::exp(logYphib + krhophi5 * dx / 2.));
        krhophi7 = boltzright_phi_log(T, std::exp(logYb + k6 * dx / 2.), std::exp(logYphib + krhophi6 * dx / 2.));
        k7 = boltzright_log(T, sv, std::exp(logYb + k6 * dx / 2.), std::exp(logYphib + krhophi6 * dx / 2.));

        T = std::exp(x + 2. * dx);
        krhophi8 = boltzright_phi_log(T, std::exp(logYb + k7 * dx), std::exp(logYphib + krhophi7 * dx));
        sv = getAverageSigmav(T);
        k8 = boltzright_log(T, sv, std::exp(logYb + k7 * dx), std::exp(logYphib + krhophi7 * dx));

        logYtemp2 = logYb + dx / 6. * (k5 + 2. * k6 + 2. * k7 + k8);
        logYphitemp2 = logYphib + dx / 6. * (krhophi5 + 2. * krhophi6 + 2. * krhophi7 + krhophi8);

        // Ytemp2=std::exp(logYtemp2);
        Yphitemp2 = std::exp(logYphitemp2);

        testY = (std::isinf(logYtemp1) || std::isnan(logYtemp1) || std::isinf(logYtemp2) || std::isnan(logYtemp2) ||
                 (logYtemp2 < -1.e5 && std::fabs(1. + logYtemp2 / 1.e100) > 1.e-5) ||
                 ((logYtemp2 != 0.) && (std::fabs((logYtemp2 - logYtemp1) / (logYtemp2)) > adapt_error)));

        testYphi =
            (phi_model.get()) *
            (std::isnan(logYphitemp1) || std::isnan(logYphitemp1) || std::isnan(logYphitemp2) ||
             std::isnan(logYphitemp2) || (logYphitemp2 < -1.e5 && std::fabs(1. + logYphitemp2 / 1.e100) > 1.e-5) ||
             ((logYphitemp2 != 0.) && (std::fabs((logYphitemp2 - logYphitemp1) / (logYphitemp2)) > adapt_error_phi)));

        if (this->full_comput)
        {
          Y_eq = Yeq(std::exp(x + 2. * dx));

          if (Y_eq > 0.)
            if (((!std::isnan(logYtemp2) && !std::isinf(logYtemp2)) &&
                 ((logYtemp2 < std::log(Y_eq)) || (std::fabs(1. - logYtemp2 / std::log(Y_eq)) < 0.01))) ||
                ((std::isnan(logYtemp2) || std::isinf(logYtemp2)) &&
                 ((logY + 2. * dx * k1 < std::log(Y_eq)) ||
                  (std::fabs(1. - (logY + 2. * dx * k1) / std::log(Y_eq)) < 0.001))))
            {
              logYtemp2 = logYtemp1 = std::log(Y_eq);
              testY = 0;
            }
        }

        testdx = std::fabs(dx / x) > 1.e-8;

        if (!std::isinf(logYtemp1) && std::exp(logYtemp1) == 0.)
          logYtemp1 = -1.e100;
        if (!std::isinf(logYphitemp1) && std::exp(logYphitemp1) == 0.)
          logYphitemp1 = -1.e100;
        if (!std::isinf(logYtemp2) && std::exp(logYtemp2) == 0.)
          logYtemp2 = -1.e100;
        if (!std::isinf(logYphitemp2) && std::exp(logYphitemp2) == 0.)
          logYphitemp2 = -1.e100;

        if (testY || testYphi)
        {
          if ((phi_model.get() != 0 && rho_phi != 0.) && (logYphitemp2 - logYphitemp1 != 0.) && logYphitemp2 > -1.e-100)
            dx *= 2. * 0.9 *
                  std::min(std::max(std::min(std::fabs(logYtemp2 * adapt_error / (logYtemp2 - logYtemp1)) *
                                                 !std::isnan(logYtemp2 / (logYtemp2 - logYtemp1)),
                                             std::fabs(logYphitemp2 * adapt_error_phi / (logYphitemp2 - logYphitemp1)) *
                                                 !std::isnan(logYphitemp2 / (logYphitemp2 - logYphitemp1))),
                                    0.3),
                           1.);
          else
            dx *= 2. * 0.9 *
                  std::min(std::max(std::fabs(logYtemp2 * adapt_error / (logYtemp2 - logYtemp1)) *
                                        !std::isnan(logYtemp2 / (logYtemp2 - logYtemp1)),
                                    0.3),
                           1.);
        }
      }

      if (this->full_comput)
      {
        if ((std::exp(logY) < (1. + delta) * Yeq(std::exp(x))) &&
            ((std::exp(logYb) > (1. + delta) * Yeq(std::exp(x + dx)))))
        {
          const real_t y1 = logY - std::log((1. + delta) * Yeq(std::exp(x)));
          const real_t x1 = x;
          const real_t y2 = logYb - std::log((1. + delta) * Yeq(std::exp(x + dx)));
          const real_t x2 = x + dx;
          const real_t y0 = 0.;
          const real_t x0 = (y0 - y1) * (x2 - x1) / (y2 - y1) + x1;

          this->Tfo = std::exp(x0);
        }
        else if ((std::exp(logYb) < (1. + delta) * Yeq(std::exp(x + dx))) &&
                 ((std::exp(logYtemp2) > (1. + delta) * Yeq(std::exp(x + 2. * dx)))))
        {
          const real_t y1 = logYb - std::log((1. + delta) * Yeq(std::exp(x + dx)));
          const real_t x1 = x + dx;
          const real_t y2 = logYtemp2 - std::log((1. + delta) * Yeq(std::exp(x + 2. * dx)));
          const real_t x2 = x + 2. * dx;
          const real_t y0 = 0.;
          const real_t x0 = (y0 - y1) * (x2 - x1) / (y2 - y1) + x1;

          this->Tfo = std::exp(x0);
        }
      }

      logY = logYtemp2;
      logYphi = logYphitemp2;

      Y = std::exp(logY);
#ifdef DEBUG
      //     fprintf(fptr, "%.4e\t%.4e\n", std::exp(x), Y);
      fYt << std::exp(x) << '\t' << Y << '\n';
#endif
      Yphi = std::exp(logYphi);

      if (Yphi < 1e-10 * Y)
      {
        Yphi = rho_phi = 0.;
        logYphi = -1.e100;
      }

      x = x + 2. * dx;

      if (phi_model.get() != 0 && rho_phi != 0. && (logYphitemp2 - logYphitemp1 != 0.))
        dx *= 2. * 0.9 *
              std::max(1., std::min(std::min(std::fabs(logYtemp2 * adapt_error / (logYtemp2 - logYtemp1)) *
                                                 !std::isnan(logYtemp2 / (logYtemp2 - logYtemp1)),
                                             std::fabs(Yphitemp2 * adapt_error_phi / (Yphitemp2 - Yphitemp1)) *
                                                 !std::isnan(Yphitemp2 / (Yphitemp2 - Yphitemp1))),
                                    2.));
      else
        dx *= 2. * 0.9 *
              std::max(1., std::min(std::fabs(logYtemp2 * adapt_error / (logYtemp2 - logYtemp1)) *
                                        !std::isnan(logYtemp2 / (logYtemp2 - logYtemp1)),
                                    2.));

      if (dx > x - std::log(Tfin))
        dx = x - std::log(Tfin);
    } while (T > Tfin);

#ifdef DEBUG
    //   fclose(fptr);
    //   fclose(sig);
    //   if(this->phi_model) fclose(phi0);
    //   printf("%d iterations were needed for the integration of the Boltzmann equation.\n\n",counter);

    fYt.close();
    fsig.close();
    if (fphi0)
      fphi0.close();
    std::cout << counter << " iterations were needed for the integration of the Boltzmann equation.\n\n";
#endif

    if ((this->mgravitino != 0.) && (this->mgravitino < this->relicmass))
      return 2.755e8 * this->mgravitino * Y;

    return 2.742e8 * this->relicmass * Y;
  }
  /*--------------------------------------------------------------*/

  real_t BoltzmannSolver::relic_density()
  /* Computes the relic density from the param and paramrelic structures, and Wefftab */
  {
    //   real_t omega=0.;
    //   if(this->solver==2) omega=relic_density_linsolver();
    //   else omega=relic_density_logsolver();
    //   return omega;
    return solver == 2 ? relic_density_linsolver() : relic_density_logsolver();
  }

  /*--------------------------------------------------------------*/
} // end of namespace __SPEC_LIB_NAME__
