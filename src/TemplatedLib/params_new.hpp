#pragma once

#ifndef __cplusplus
#define __cplusplus
#endif
#include "config.hpp"
#include "cparams.h"
#include "correspondance.hpp"
#include "marty/lha/lha.h"
#include <chrono>
#include <cstring>

// #define DEBUG

namespace scalar2to2
{

inline real_t get_value( mty::lha::LHAFileData data, const std::string block, const int N)
{
    return data.getValue(block, N).value();  
}

inline real_t get_value( mty::lha::LHAFileData data, const std::string block, const int N, const real_t pdgvalue)
{
    return data.getValue(block, N).value_or(pdgvalue);  
}

inline real_t get_value_mat( mty::lha::LHAFileData data, const std::string block, const int N)
{
    return data.getValue(block, N, N).value();
}

inline real_t get_value_mat( mty::lha::LHAFileData data, const std::string block, const int N, const real_t pdgvalue)
{
    return data.getValue(block, N, N).value_or(pdgvalue);
}

namespace mathematicalconstants
{
  constexpr const real_t pi = M_PI;
  constexpr const real_t zeta3= 3.1415926535897932385;
}

namespace units_conversion //Taken from the exact values of the PDG 2022
{
  constexpr const real_t c_in_m_s = 299792458.; // Speed of light in vacuo in m/s
  constexpr const real_t h_in_J_s = 6.62607015e-34;// Planck's constant in J*s

  constexpr const real_t hbar_in_J_s = h_in_J_s/2./mathematicalconstants::pi; // Reduced Planck's constant in J*s
  constexpr const real_t hbar_in_MeV_s = 6.582119569e-22; // Reduced Planck's constant in MeV*s

  constexpr const real_t hbarc_MeV_fm = 197.3269804; // hbar*c in MeV*fm
  
  constexpr const real_t hbarc2_Gev2_mbarn = 0.389379372;// (hbar*c)^2 in Gev^2 * mbarn
  
  constexpr const real_t e_C = 1.602176634e-19; // Elementary charge in Coulomb

  // Energy -> length^-1
  constexpr const real_t MeV_to_fmm1 = 1./hbarc_MeV_fm; // conversion from MeV to fm^-1
  constexpr const real_t GeV_to_pmm1 = MeV_to_fmm1;     // conversion from GeV to pm^-1
  constexpr const real_t GeV_to_cmm1 = GeV_to_pmm1*1.0e+10; // conversion from GeV to cm^-1
  constexpr const real_t GeV_to_mm1  = GeV_to_pmm1*1.0e+12; // conversion from GeV to m^-1

  // Energy^-2 -> length^2
  constexpr const real_t GeVm2_to_mbarn = hbarc2_Gev2_mbarn;// conversion from GeV^-2 to mbarn
  constexpr const real_t GeVm2_to_pbarn = GeVm2_to_mbarn*1.0e+9; // conversion from GeV^-2 to pbarn

  // Energy -> time^-1
  constexpr const real_t GeV_to_sm1  = GeV_to_mm1/c_in_m_s; // conversion from GeV to s^-1
  
  // Need to convert the GeV^-2 into cm^2 (surface) and mupliply by c in cm/s to have sigmav in cm^3/s
  constexpr const real_t GeVm2_to_cm3_over_s = c_in_m_s*1.0e+2/(GeV_to_cmm1*GeV_to_cmm1);  // conversion from GeV^-2 to cm^3/s
}

namespace pdg2016Value
{
    // Default values from PDG 2016 and their lha block code convention
    
    // 1/Alpha_em at the scale E = 0
    constexpr const real_t alpha_em_inv =  137.035999710;
    
    // Block SMINPUTS N=3
    constexpr const real_t alpha_str_Mz = 0.1181;
    
    // m_b(m_b) in the MS bar scheme
    // Block SMINPUTS N=5
    constexpr const real_t mb_mb = 4.18;
    
    // top quark pole mass 
    // Block SMINPUTS N=6    
    constexpr const real_t m_top_pole = 173.34;
    // From the kinematics of ttbar events
    
    // Block SMINPUTS N=7
    constexpr const real_t m_tau_pole = 1.77686;
    
    // electron pole mass
    // Block SMINPUTS N=11
    constexpr const real_t m_el_pole = 510.9989461e-6; /* electron mass in GeV */
    
    // muon pole mass
    // Block SMINPUTS N=13
    constexpr const real_t m_mu_pole = 0.105658;
    
    // m_d at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=21
    constexpr const real_t m_d2GeV = 4.7e-3;
    
    // m_u at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=22
    constexpr const real_t m_u2GeV = 2.2e-3;
    
    // m_s at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=23
    constexpr const real_t m_s2GeV = 9.6e-2;
    
    // charm m_c(m_c) mass in the MSbar scheme
    // Block SMINPUTS N=24
    constexpr const real_t m_c_m_c = 1.27;
    
    // W mass
    // Block MASS N=24
    constexpr const real_t m_W = 80.385;
    
    // Z mass
    // Block SMINPUTS N=4
    constexpr const real_t m_Z = 91.1876;
    
    
    constexpr const real_t theta_W = 0.490449;
    // CKM parameters 
    constexpr const real_t delta_wolf = 1.20146;
    constexpr const real_t V_ud = 0.974339;
    constexpr const real_t V_us = 0.225057;
    constexpr const real_t V_ub_mod = 0.00366661;
    constexpr const complex_t V_cd = {-0.224922, -0.000136891};
    constexpr const complex_t V_cs = { 0.97351, -3.16198e-05};
    constexpr const real_t V_cb = 0.0410782;
    constexpr const complex_t V_td = {0.00795853, -0.00332961};
    constexpr const complex_t V_ts = {-0.0403218, -0.00076909};
    constexpr const real_t V_tb = 0.999149;
  
    constexpr const real_t Gamma_Z = 2.4952;
    constexpr const real_t Gamma_W = 2.085;
    
    constexpr const real_t m_e = 510.9989461e-6; /* electron mass in GeV */
    constexpr const real_t k_B = 8.617330e-5;		/* Boltzmann's constant in GeV/GK */
    constexpr const real_t alphaem = 1./alpha_em_inv; /* 0.007297353 fine-structure constant */
    constexpr const real_t sigma_SB = 0.16449340668482282; /* Stefan-Boltzmann constant = pi^2/60 */
    
    constexpr const real_t Gn = 6.67428e-8;  /* Newton's constant in cm^3.g^-1.s^-2 */
    constexpr const real_t G = 6.708609142443796e-39; /* Gn*pow(m_to_GeV,3.)*pow(g_to_GeV,-1.)*pow(s_to_GeV,-2.) Newton constant in GeV^-2 */
    constexpr const real_t Mplanck = 1.2209102930946623e+19; /* in GeV, more precise definition than before */
    constexpr const real_t DMpn = 0.0012934; /* mass difference between neutron and proton in GeV */
    constexpr const real_t M_u = 1.660538921e-24; /* unit atomic mass M_u in g */
    constexpr const real_t zeta = 1.6103162253325862; /* 3*k_B/(2*c^2*Mu) in GeV^-1 */

    constexpr const real_t K_to_eV = k_B; /* conversion factor T(10**9 K) * K_to_GeV = T(GeV) or T(K) * K_to_eV = T(eV) */
    

    constexpr const real_t g_to_GeV = 5.60958884538932e+23; /* conversion factor M(g) * g_to_GeV = M(GeV) */
    constexpr const real_t kg_to_GeV = 5.60958884538932e+26; /* conversion factor M(kg) * kg_to_GeV = M(GeV) */
    
    constexpr const real_t m_to_GeV =  5.067730582705779e+15; /* conversion factor L(m) * m_to_GeV = L(GeV^-1) */
    constexpr const real_t cm_to_GeV = 5.067730582705779e+13; /* conversion factor L(cm) * cm_to_GeV = L(GeV^-1) */
    constexpr const real_t s_to_GeV  = 1.5192674078711377e+24; /* conversion factor t(s) * s_to_GeV = t(GeV^-1) */
    
}

namespace pdgValue
{
    // Default values from PDG 2022 and their lha block code convention
    
    /* In the block SM_INPUT 
     * a set of SM parameter is provided in order to compute 
     * the spectrum and other quantities 
     * 
     * In the block MASS
     * the spectrum (pole masses) has to be provided 
     * 
     * We list here the reference values taken from the PDG 2022
     * They are chosen if no alternative value is provided in the lha file
     */
    

  
    // 1/Alpha_em at the scale E = 0
    constexpr const real_t alpha_em_inv =  137.035999180;
    
    // 1/Alpha_em at the scale E = m_tau in the MS bar scheme
    constexpr const real_t alpha_em_Mtau_inv = 133.471;
    
    // 1/Alpha_em at the scale E = M_Z in the MS bar scheme
    // Block SMINPUTS N=1
    constexpr const real_t alpha_em_MZ_inv =  127.951 ;
    
    // Fermi constant in GeV-2
    // Block SMINPUTS N=2
    constexpr const real_t Gfermi =  1.1663788e-5;
    
    // Alpha_strong at the scale E = M_Z in the MS bar scheme
    // Block SMINPUTS N=3
    constexpr const real_t alpha_str_Mz = 0.1179;
    
    // Z pole mass
    // Block SMINPUTS N=4
    constexpr const real_t m_Z_pole = 91.1876-0.034;
    
    // M_Z(M_Z) in the MS bar scheme
    // Block MASS N=23
    constexpr const real_t m_Z_m_Z = 91.1876;
    // Experimental velue from Breit-Wigner distribution (M_Z(M_Z))
    
    // m_b(m_b) in the MS bar scheme
    // Block SMINPUTS N=5
    constexpr const real_t mb_mb = 4.18;
    
    // m_b pole computed from mb_mb
    // Block MASS N=5
    constexpr const real_t mb_pole = 4.78;
    
    // top quark pole mass 
    // Block SMINPUTS N=6    
    constexpr const real_t m_top_pole = 172.69;
    // From the kinematics of ttbar events
    
    // tau pole mass
    // Block SMINPUTS N=7
    constexpr const real_t m_tau_pole = 1.77686;
    
    // electron pole mass
    // Block SMINPUTS N=11
    constexpr const real_t m_el_pole = 5.1099895000e-4;
    
    // muon pole mass
    // Block SMINPUTS N=13
    constexpr const real_t m_mu_pole = 0.1056583755;
    
    // m_d at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=21
    constexpr const real_t m_d2GeV = 4.67e-3;
    
    // m_u at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=22
    constexpr const real_t m_u2GeV = 2.16e-3;
    
    // m_s at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=23
    constexpr const real_t m_s2GeV = 9.34e-2;
    
    // charm m_c(m_c) mass in the MSbar scheme
    // Block SMINPUTS N=24
    constexpr const real_t m_c_m_c = 1.27;
    
    // CKM parameters in the Wolfenstein parametrisation
    // Block VCKMIN N=1,2,3,4 respectively
    constexpr const real_t ckm_lambda = 0.22500;
    constexpr const real_t ckm_A      = 0.826;
    constexpr const real_t ckm_rhobar = 0.159;
    constexpr const real_t ckm_etabar = 0.348;
    

    // CKM matrix assignment from Wolfenstein values
    /***************************************************************/
    // const inline static complex_t I{0.0,1.0};

    const inline static real_t s12=ckm_lambda;
    const inline static real_t s23=ckm_A*ckm_lambda*ckm_lambda;
    
    const inline static complex_t rhopIeta={ckm_rhobar,ckm_etabar};
    
    const inline static real_t s13= std::abs(ckm_A*std::pow(ckm_lambda,3.)*(rhopIeta)*std::sqrt(1.-std::pow(ckm_A*ckm_lambda*ckm_lambda,2.))/std::sqrt(1.-ckm_lambda*ckm_lambda)/(1.-std::pow(ckm_A*ckm_lambda*ckm_lambda,2.)*(rhopIeta)));

    
    const inline static complex_t expid=
       ((complex_t)(ckm_A)*std::pow(ckm_lambda,3.)*
       (rhopIeta)*
       std::sqrt(1.-std::pow(ckm_A*ckm_lambda*ckm_lambda,2.))
       /std::sqrt(1.-ckm_lambda*ckm_lambda)
       /(1.-std::pow(ckm_A*ckm_lambda*ckm_lambda,2.)*(rhopIeta)))
       /s13;

    const inline static real_t c12=std::sqrt(1.-s12*s12);
    const inline static real_t c13=std::sqrt(1.-s13*s13);
    const inline static real_t c23=std::sqrt(1.-s23*s23);

    const inline static complex_t CKM[4][4] = 
     { {0    ,                            0,       0,     0},
       {0    ,c12*c13                     ,c12*c13,s13/expid},
       { 0   , -s12*c23-c12*s23*s13*expid ,c12*c23-s12*s23*s13*expid, s23*c13},
       {0,   s12*s23-c12*c23*s13*expid, -c12*s23-s12*c23*s13*expid,c23*c13}
     };
       

    const real_t V_ud=CKM[1][1].real();
    const real_t V_us=CKM[1][2].real();
    const real_t V_ub_mod=s13;
    const real_t delta_wolf = std::log(expid).imag();
    const complex_t V_cd=CKM[2][1];
    const complex_t V_cs=CKM[2][2];
    const real_t V_cb=CKM[2][3].real();
    const complex_t V_td=CKM[3][1];
    const complex_t V_ts=CKM[3][2];
    const real_t V_tb=CKM[3][3].real();
    /**************************************************************/
    
    constexpr const real_t m_h = 125.25;
    
    constexpr const real_t Gamma_Z = 2.4952;
    constexpr const real_t Gamma_W = 2.085;
    
    constexpr const real_t m_e = m_el_pole;
    // Boltzmann constant in eV/K
    constexpr const real_t k_B=8.617333262e-5; // Bolzmann constant in eV/K
    const real_t alphaem = 1./alpha_em_inv;
    
    constexpr const real_t Gn = 6.67430e-8;/* Newton's constant in cm^3.g^-1.s^-2 */
    constexpr const real_t G = 6.70883e-39; /* Gn*pow(m_to_GeV,3.)*pow(g_to_GeV,-1.)*pow(s_to_GeV,-2.) Newton constant in GeV^-2 */
    constexpr const real_t Mplanck = 1.220890e+19; /* in GeV, more precise definition than before */
    
    constexpr const real_t M_u = 1.6605306660e-24; /* unit atomic mass M_u in g */

    constexpr const real_t K_to_eV = k_B; /* conversion factor T(10**9 K) * K_to_GeV = T(GeV) or T(K) * K_to_eV = T(eV) */
    
    constexpr const real_t kg_to_eV = 5.609588603e+35;
    
    constexpr const real_t kg_to_GeV = kg_to_eV*1.0e-9;
    constexpr const real_t g_to_GeV= 1.0e-3*kg_to_eV*1.0e-9;

    constexpr const real_t hbarc=197.3269804; // hbar * c in MeV*fm
    constexpr const real_t fm_to_eVm1= 1./hbarc*1.0e-6;
    
    constexpr const real_t m_to_GeV = 1.0e+15*fm_to_eVm1*1.0e+9;
    constexpr const real_t cm_to_GeV = 1.0e-2*m_to_GeV;
    
    constexpr const real_t c_em = 299792458; // Speed of light in m/s
    const real_t s_to_GeV  = c_em*m_to_GeV; /* conversion factor t(s) * s_to_GeV = t(GeV^-1) */
  
    constexpr const real_t hbar = 6.582119569e-22; //hbar in Mev *s 
    const real_t sigma_SB = pdg2016Value::sigma_SB;
    const real_t zeta = pdg2016Value::zeta; /* 3*k_B/(2*c^2*Mu) in GeV^-1 */
    const real_t DMpn = pdg2016Value::DMpn; /* mass difference between neutron and proton in GeV */
}

struct darkpackparam_t
{
  public :
    /* Parameters necessary for the running */
    csl::InitSanitizer<real_t> Running_scale { "Running_scale"};
    //

    // Alpha_strong at the scale E = M_Z in the MS bar scheme
    // Block SMINPUTS N=3
    csl::InitSanitizer<real_t> alpha_str_Mz { "alpha_strong(M_Z)"};
    //

    // Z pole mass
    // Block SMINPUTS N=4
    csl::InitSanitizer<real_t> m_Z_pole { "m_Z,pole" };

    // m_b(m_b) in the MS bar scheme
    // Block SMINPUTS N=5
    csl::InitSanitizer<real_t> mbmb { "m_b(m_b)" };
    //

    // top quark pole mass
    // Block SMINPUTS N=6
    csl::InitSanitizer<real_t> mtop_pole { "m_top,pole" };
    //

    /* Parameters as read from lha file */

    // 1/Alpha_em at the scale E = M_Z in the MS bar scheme
    // Block SMINPUTS N=1
    csl::InitSanitizer<real_t> alpha_em_MZ_inv {"alpha_em_MZ_inv"};
    // Fermi constant in GeV-2
    // Block SMINPUTS N=2
    csl::InitSanitizer<real_t> Gfermi  {"G_F"};
    // tau pole mass
    // Block SMINPUTS N=7
    csl::InitSanitizer<real_t> m_tau_pole {"m_tau_pole"};
    // electron pole mass
    // Block SMINPUTS N=11
    csl::InitSanitizer<real_t> m_el_pole {"m_el_pole"};
    // muon pole mass
    // Block SMINPUTS N=13
    csl::InitSanitizer<real_t> m_mu_pole {"m_mu_pole"};
    // m_d at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=21
    csl::InitSanitizer<real_t> m_d2GeV {"m_d2GeV"};
    // m_u at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=22
    csl::InitSanitizer<real_t> m_u2GeV {"m_u2GeV"};
    // m_S at 2 GeV in the MSbar scheme
    // Block SMINPUTS N=22
    csl::InitSanitizer<real_t> m_s2GeV {"m_s2GeV"};
    // charm m_c(m_c) mass in the MSbar scheme
    // Block SMINPUTS N=24
    csl::InitSanitizer<real_t> m_c_m_c {"m_c_m_c"};

    // CKM parameters in the Wolfenstein parametrisation
    // Block VCKMIN N=1,2,3,4 respectively
    csl::InitSanitizer<real_t> ckm_lambda { "ckm_lambda"};
    csl::InitSanitizer<real_t> ckm_A      { "ckm_A"};
    csl::InitSanitizer<real_t> ckm_rhobar { "ckm_rhobar"};
    csl::InitSanitizer<real_t> ckm_etabar { "ckm_etabar"};
    
    // Boolean that determines if, in the instanciated Param_t 
    // variable, the masses_vector array has been initialised
    bool aremassesassigned=false;
    
    // Boolean that determines if, in the instanciated Param_t 
    // variable, the width_vector array has been initialised
    bool arewidthsassigned=false;
    
    // The masses_vector array contains the masses of all the 
    // particles. They appear in the same order as they are defined 
    // in the corr::Part_t enumeration in correspondance.h.
    std::array<csl::InitSanitizer<real_t>,corr::TOTAL_PARTICLES+1>  masses_vector;
    
    // The widths_vector array contains the widths of all the 
    // particles. They appear in the same order as they are defined 
    // in the corr::Part_t enumeration in correspondance.h.
    std::array<csl::InitSanitizer<real_t>,corr::TOTAL_PARTICLES+1>  widths_vector;

protected :
    // Integer value to label the lightest BSM particle, trough the 
    // enumeration corr::Part_t in correspondance.h.
    csl::InitSanitizer<int> part_lbsm { "part_lbsm"};
    
    // Floating-point value to store the value of the lightest 
    // BSM particle's mass
    csl::InitSanitizer<real_t> m_lbsm { "m_lbsm"};
    
    // Floating-point value to store the value of the lightest 
    // BSM particle's electric charge
    csl::InitSanitizer<real_t> em_charge_lbsm{ "q_lbsm"};
    
    // Integer value to store the value of the lightest 
    // BSM particle's spin-helicity degrees of freedom
    csl::InitSanitizer<int> g_lbsm { "g_lbsm"};
    
public :     
  
    inline void empty_init()
    {
      // Basic initilisation (no values assigned)
      aremassesassigned=false;
      arewidthsassigned=false;
      for(size_t i = 0 ; i < corr::TOTAL_PARTICLES+1 ; i++)
      {
        // Here we set the names of the csl::InitSanitizer variables
        // which are contained in arrays
        masses_vector[i].setName("m_"+corr::part_names[i]); 
        widths_vector[i].setName("Gamma_"+corr::part_names[i]);
      }
    };
    
    darkpackparam_t()
    {
     empty_init(); 
    }

    // This method returns the integer value in the enumeration 
    // "corr::Part_t" which corresponds to the lightest BSM 
    // particle
    inline int getLightestBSMpart() const
    { return part_lbsm;};

    // This method returns the floating-point value for the mass
    // of the lightest BSM particle
    inline real_t getLightestBSMmass() const
    { return m_lbsm;};

    // This method returns the floating-point value for the electric 
    // charge of the lightest BSM particle
    inline real_t getLightestBSMem_charge() const
    { return em_charge_lbsm;};

    // This method returns the integer value for the spin-helicity
    // degrees of freedom of the lightest BSM particle
    inline int getLightestBSMdof() const
    { return g_lbsm;};
    
    // Comparison operators
    bool compare_with_darkpackparam( const darkpackparam_t &) const  ;
    
    virtual bool operator==(const darkpackparam_t &) const  ;

    virtual bool operator!=(const darkpackparam_t &other) const 
    {
      return !(*this == other);
    };
    
    void print(std::ostream& out=std::cout) const ;
};
  
struct Param_t : public darkpackparam_t , param_t
{

public :
  
    void ReadSMparams(const std::string namefile);
  
    Param_t() : darkpackparam_t()
    {};    
    
    // This constructor instantiates a variable by reading the content
    // of a LHA file located in "filename". 
    // If the user does not create a leshouchesfrommarty.h/.cpp in 
    // the model-dependent libraries, the installation script 
    // does that, by creating an empty library
    Param_t(const std::string filename);
    
    // Procedure to automatically assign the elements in the 
    // "masses_vector" array, by calling corr::getMassFirst
    inline void AssignMassesVector()
    {
      using namespace corr;
      for(int i = 1 ; i <= TOTAL_PARTICLES ; i++)
      {
        masses_vector[i]= std::abs(getMassFirst(i, *this));
      }
      aremassesassigned=true;
    };
    
    // Procedure to determine the lightest BSM particle, and 
    // assign the relevant data member in the instance of this 
    // structure
    void setLightestBSMparticle() ;
    
    // Procedure to modify the content of a class instance with 
    // the one of another, by checking that the SM parameters 
    // and the DM candidate remain the same
    // Returns 0 for success
    unsigned int  change_values(const Param_t &other);
    
    
    // Comparison operators
    bool operator==(const Param_t &) const ;

    bool operator!=(const Param_t &other) const 
    {
      return !(*this == other);
    };

    void Print(std::ostream& out=std::cout) const ;
};

// The following class is defined to have an easy way of measuring time
class Timer
{
  std::chrono::_V2::system_clock::time_point start, end;
  
  std::chrono::duration<float> duration;
  
  public:
    Timer()
    {
      start = std::chrono::high_resolution_clock::now();
    }
    
    real_t getElapsed()
    {
      end = std::chrono::high_resolution_clock::now();
      duration = end-start;
      return static_cast<real_t>(duration.count());
    } 
      
    ~Timer()
    {
      end = std::chrono::high_resolution_clock::now();
      duration = end-start;
      std::cout <<"Timer took " << duration.count()*1000 << " ms" << std::endl;
    }
};

} // end of namespace scalar2to2



inline std::ostream& operator << (std::ostream& out, const scalar2to2::param_t &par)
{
  par.print(out); 
  return out;
}

inline std::ostream& operator << (std::ostream& out, const scalar2to2::Param_t &par)
{
  par.Print(out); 
  return out;
}
