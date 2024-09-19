#pragma once 
#ifndef __cplusplus
#define __cplusplus
#endif

#include "config.hpp"
#include "params_new.hpp"
#include <vector>

namespace scalar2to2
{

constexpr const static int NTABMAX=1000;

struct DataWithTemperature_t
{
  real_t T;
  real_t sqrtgstar;
  real_t heff;
  
  DataWithTemperature_t(){};
  
  DataWithTemperature_t(const real_t &T_in, 
                      const real_t &sqrtgstar_in,
                      const real_t &heff_in
                      )
    : T(T_in), sqrtgstar(sqrtgstar_in), heff(heff_in) {};
    
  void print(std::ostream &out=std::cout) const 
  {
    out << "T = " << T << ' ';
    out << "sqrtgstar = " << sqrtgstar << ' ';
    out << "heff = " << heff << '\n';
  };
  
};


struct Relicparam_t
{
    /* structure containing the cosmological model parameters */
  public:
    enum QCDeosModel : int {A=1, B, B2, B3, C, Bonn, IdealGas};

  protected:
    csl::InitSanitizer<int> model_eff { "model_eff" };
  
  private:
    std::array<DataWithTemperature_t,276> dataT;
    // Note that in the vector, temperatures will be stored in a decreasing order    
    
  public:
    // Original members from SuperIso Relic V. 4
    
    csl::InitSanitizer<int> entropy_model { "entropy_model" };     // Entropy model switch
    csl::InitSanitizer<int> energy_model { "energy_model" };       // Energy model switch

    csl::InitSanitizer<real_t> dd0 { "dd0" };             // Dark matter density
    csl::InitSanitizer<real_t> ndd { "ndd" };             // Dark matter density slope
    csl::InitSanitizer<real_t> Tdend { "Tdend" };         // Dark matter density cutoff temperature
    csl::InitSanitizer<real_t> Tddeq { "Tddeq" };         // Dark matter density equality temperature

    csl::InitSanitizer<real_t> sd0 { "sd0" };            // Dark entropy density
    csl::InitSanitizer<real_t> nsd { "nsd" };            // Dark entropy density slope
    csl::InitSanitizer<real_t> Tsend { "Tsend" };        // Dark entropy density cutoff temperature

    csl::InitSanitizer<real_t> Sigmad0 { "Sigmad0" };        // Dark entropy injection amplitude
    csl::InitSanitizer<real_t> nSigmad { "nSigmad" };        // Dark entropy injection slope
    csl::InitSanitizer<real_t> TSigmadend { "TSigmadend" };  // Dark entropy injection cutoff temperature

    csl::InitSanitizer<real_t> Sigmarad0 { "Sigmarad0" };      // Standard entropy injection amplitude
    csl::InitSanitizer<real_t> nSigmarad { "nSigmarad" };      // Standard entropy injection slope
    csl::InitSanitizer<real_t> TSigmaradend { "TSigmaradend" };// Standard entropy injection cutoff temperature

    csl::InitSanitizer<real_t> nt0 { "nt0" };            // Non-thermal production of relics amplitude
    csl::InitSanitizer<real_t> nnt { "nnt" };            // Non-thermal production of relics slope
    csl::InitSanitizer<real_t> Tnend { "Tnend" };        // Non-thermal production of relics cutoff temperature

    csl::InitSanitizer<int> coupd { "coupd" };           // Dark fluid coupling to plasma switch

    csl::InitSanitizer<real_t> quintn2 { "quintn2" };        // Effective quintessence model n2 parameter
    csl::InitSanitizer<real_t> quintn3 { "quintn3" };        // Effective quintessence model n3 parameter
    csl::InitSanitizer<real_t> quintn4 { "quintn4" };        // Effective quintessence model n4 parameter
    csl::InitSanitizer<real_t> quintT12 { "quintT12" };      // Effective quintessence model T12 parameter
    csl::InitSanitizer<real_t> quintT23 { "quintT23" };      // Effective quintessence model T23 parameter
    csl::InitSanitizer<real_t> quintT34 { "quintT34" };      // Effective quintessence model T34 parameter

    csl::InitSanitizer<int> phi_model { "phi_model" };             // Decaying scalar field model switch
    csl::InitSanitizer<real_t> eta_phi { "eta_phi" };              // Decaying scalar field model eta parameter
    csl::InitSanitizer<real_t> Gamma_phi { "Gamma_phi" };          // Decaying scalar field model decay rate
    csl::InitSanitizer<real_t> rhot_phi_Tmax { "rhot_phi_Tmax" };  // Decaying scalar field model energy density at Tmax
    csl::InitSanitizer<real_t> n_phi { "n_phi" };                  // Decaying scalar field model power law index
    csl::InitSanitizer<real_t> rhot_phi0 { "rhot_phi0" };          // Decaying scalar field model initial energy density
    csl::InitSanitizer<real_t> Tphi0 { "Tphi0" };                  // Decaying scalar field model initial temperature

    
    csl::InitSanitizer<real_t> T_RH { "T_RH" };                              // Temperature at the reheating epoch
    csl::InitSanitizer<real_t> Sigmatildestar { "Sigmatildestar" };          // Approximated cross section times relative velocity at freeze out
    csl::InitSanitizer<real_t> Sigmatildestar_max { "Sigmatildestar_max" };  // Maximum of the approximated cross section times relative velocity
    csl::InitSanitizer<real_t> Tstdstar_max { "Tstdstar_max" };              // Temperature where the cross section times relative velocity is maximum

    csl::InitSanitizer<real_t> mgravitino { "mgravitino" };      // Gravitino mass

    csl::InitSanitizer<real_t> relicmass { "relicmass" };        // Mass of the relic particles
    csl::InitSanitizer<int> scalar { "scalar" };             // Switch for the MSSM annihilation processes

    csl::InitSanitizer<int> solver { "solver" };                 // Switch for linear or logarithmic differential equation solver
    csl::InitSanitizer<int> beta_samples { "beta_samples" };     // Number of samples for the Coulomb log

    csl::InitSanitizer<real_t> Tfo { "Tfo" };                   // Freeze-out temperature
    csl::InitSanitizer<real_t> Tmax { "Tmax" };                 // Maximum temperature

    csl::InitSanitizer<int> full_comput { "full_comput" };      // Switch to deactivate the fast freeze-out temperature determination

    csl::InitSanitizer<int> use_table_rhoPD { "use_table_rhoPD" };                     // Switch for using the table of rho * PD
    csl::InitSanitizer<real_t> table_rhoPD[2][NTABMAX];  // Table of rho * PD
    csl::InitSanitizer<int> size_table_rhoPD { "size_table_rhoPD" };                   // Size of the table of rho * PD

    /*---------------------*/
    /* AlterBBN parameters */
    /*---------------------*/

    csl::InitSanitizer<int> err { "err" };             // Error code
    csl::InitSanitizer<int> failsafe { "failsafe" };   // Switch for the integration method
    csl::InitSanitizer<real_t> eta0 { "eta0" };        // Initial baryon to photon ratio
    csl::InitSanitizer<real_t> Nnu { "Nnu" };          // Number of neutrinos (e+- included)
    csl::InitSanitizer<real_t> dNnu { "dNnu" };        // Number of extra neutrinos (delta N_nu)
    csl::InitSanitizer<real_t> life_neutron { "life_neutron" };              // Neutron lifetime
    csl::InitSanitizer<real_t> life_neutron_error { "life_neutron_error" };  // Error on the neutron lifetime
    csl::InitSanitizer<real_t> xinu1 { "xinu1" };       // Degeneracy parameter for e- and neutrino
    csl::InitSanitizer<real_t> xinu2 { "xinu2" };       // Degeneracy parameter for muon and neutrino
    csl::InitSanitizer<real_t> xinu3 { "xinu3" };       // Degeneracy parameter for tau and neutrino
    csl::InitSanitizer<real_t> m_chi { "m_chi" };       // Mass of WIMP
    csl::InitSanitizer<real_t> g_chi { "g_chi" };       // Coupling of WIMP to SM particles
    csl::InitSanitizer<real_t> Tinit { "Tinit" };       // Initial temperature
    csl::InitSanitizer<real_t> Tnudec { "Tnudec" };     // Neutrino decoupling temperature
    csl::InitSanitizer<int> wimp { "wimp" };            // Switch to enable (1) / disable (0) WIMPs
    csl::InitSanitizer<int> SMC_wimp { "SMC_wimp" };    // WIMP coupling to SM particles. 1 for EM, 2 for neutrino, 3 for neutrino and equivalent neutrino
    csl::InitSanitizer<int> selfConjugate { "selfConjugate" };
    csl::InitSanitizer<int> fermion { "fermion" };
    csl::InitSanitizer<int> EM_coupled { "EM_coupled" };
    csl::InitSanitizer<int> neut_coupled { "neut_coupled" };
    csl::InitSanitizer<int> neuteq_coupled { "neuteq_coupled" };
    csl::InitSanitizer<real_t> chi2 { "chi2" };
    csl::InitSanitizer<int> nobs { "nobs" };
    csl::InitSanitizer<real_t> fierz { "fierz" };       // Fierz interference term from LQ sector
    csl::InitSanitizer<real_t> B_chi { "B_chi" };       // branching ratio of WIMP DM of mass m_p < m_chi < m_n to explain tau_n anomaly
    csl::InitSanitizer<real_t> rhob0 { "rhob0" };       // current baryon density
    csl::InitSanitizer<real_t> b_cdm_ratio { "b_cdm_ratio" }; // current ratio of baryon density to cold dark matter density

    // Methods for initialisation or costruction
    
    void setQCDeosModel(const int x);
    void init_heffgeff();
    void init_heffgeff(const int x);
    
    void init(const real_t &mrelic, const int x);
    
    Relicparam_t(const real_t &mrelic=-1., const int x=2);  
    virtual ~Relicparam_t(){};
    
    
    // Metod to print the elements
    void print_relicparam(std::ostream& out=std::cout) const ;
    virtual void print(std::ostream& out=std::cout)const{print_relicparam(out);};
    
    // The following methods need to be called after the constructor
    // Their purpose is to initialise this class with modified cosmological model
    // according to the same convention of SuperIso Relic v4
    void Init_cosmomodel_param(const real_t &eta, const real_t &Nnu_local, const real_t &dNnu_local, const real_t &life_neutron_local, const real_t &life_neutron_error_local, const real_t &xinu1_local, const real_t &xinu2_local, const real_t &xinu3_local);
    void Init_wimp(const real_t &mass_wimp, int EM_coupled_local, int neut_coupled_local, int neuteq_coupled_local, int fermion_local, int selfConjugate_local, const real_t &g_chi_local);
    void Init_dark_density(const real_t &dd0, const real_t &ndd, const real_t &T_end);
    void Init_dark_density2(const real_t &ndd, const real_t &Tddeq, const real_t &T_end);
    void Init_dark_coupling(int coupD);
    void Init_quintessence(const real_t &T12, const real_t &n2, const real_t &T23, const real_t &n3, const real_t &T34, const real_t &n4);
    void Init_dark_entropy(const real_t &sd0, const real_t &nsd, const real_t &T_end);
    void Init_dark_entropySigmaD(const real_t &Sigmad0, const real_t &nSigmad, const real_t &T_end);
    void Init_entropySigmarad(const real_t &Sigmarad0, const real_t &nSigmarad, const real_t &T_end);
    void Init_nonthermal(const real_t &nt0, const real_t &nnt, const real_t &T_end);
    void Init_gravitino(const real_t &mgravitino);
    void Init_scalarfield(const real_t &rhotilde_phi_Tmax, const real_t &Tmax, const real_t &T_RH, const real_t &eta_phi, const real_t &n_phi);
    void Init_dark_density_table(real_t table[2][NTABMAX], int nlines);
    void Init_neutron_decay(const real_t &tau, const real_t &tau_err, const real_t &fierz, const real_t &m_chi, const real_t &B_chi);


    // Methods to get thermodynamical quantities
    
    real_t getheff(const real_t &T);
    real_t getsqrtgstar(const real_t &T);
    real_t getgeff(const real_t &T);
    
    inline real_t heff(const real_t &T)  { return getheff(T);};
    inline real_t sqrtgstar(const real_t &T) { return getsqrtgstar(T);} ;
    inline real_t geff(const real_t &T) {return getgeff(T);};
    
    real_t dark_density(const real_t &T);
    real_t dark_density_pressure(const real_t &T);
    real_t sigma_entropy(const real_t &T);
    real_t dark_entropy(const real_t &T);
    real_t dark_entropy_derivative(const real_t &T);
    real_t dark_entropy_Sigmad(const real_t &T);
    real_t dark_entropy_Sigmarad(const real_t &T);
    real_t entropy_Sigmarad(const real_t &T);
    real_t nonthermal(const real_t &T);
    real_t neutdens(const real_t &Tnu);
    real_t neutdens_deriv(const real_t &Tnu);

};

real_t neutN(const real_t &T);

}// end of namespace scalar2to2

inline std::ostream& operator << (std::ostream& out, const scalar2to2::DataWithTemperature_t &var)
{
  var.print(out);
  return out;
}

inline std::ostream& operator << (std::ostream& out , const scalar2to2::Relicparam_t &var)
{
  var.print_relicparam(out);
  return out;
}
