#include "relicparam.hpp"
#include "macros.hpp"

// #define DEBUG
template <typename T>
static inline T SQUARE(const T x) {return x*x;}

namespace scalar2to2
{

void Relicparam_t::setQCDeosModel(const int x)
{
  model_eff=x;
}

void Relicparam_t::init_heffgeff()
{  
  int i;

#ifdef DEBUG
  std::cout << "Initialising with " << model_eff << '\n';
#endif
  
  switch(model_eff)
  {
    case 1:
    {
      const real_t tableA[276][3]=
      {
        #include "sgStar_heff/sgStar_heff_A.tab"
      };
      for(i=0;i<=275;i++) 
      {
        dataT[i]=DataWithTemperature_t(tableA[i][0], tableA[i][1], tableA[i][2]);
      }
      break;
    }
    case 2:
    {
      const real_t tableB[276][3]=
      {
        #include "sgStar_heff/sgStar_heff_B.tab"
      };
      for(i=0;i<=275;i++) 
      {
        dataT[i]=DataWithTemperature_t(tableB[i][0], tableB[i][1], tableB[i][2]);
      }
      break;
    }
    case 3:
    {
      const real_t tableB2[276][3]=
      {
       #include "sgStar_heff/sgStar_heff_B2.tab"
      };
      for(i=0;i<=275;i++) 
      {
        dataT[i]=DataWithTemperature_t(tableB2[i][0], tableB2[i][1], tableB2[i][2]);
      }
      break;
    }
    case 4:
    {
      const real_t tableB3[276][3]=
      {
        #include "sgStar_heff/sgStar_heff_B3.tab"
      };
      for(i=0;i<=275;i++) 
      {
        dataT[i]=DataWithTemperature_t(tableB3[i][0], tableB3[i][1], tableB3[i][2]);
      }
      break;
    }
    case 5:
    {
      const real_t tableC[276][3]=
      {
        #include "sgStar_heff/sgStar_heff_C.tab"
      };
      for(i=0;i<=275;i++) 
      {
        dataT[i]=DataWithTemperature_t(tableC[i][0], tableC[i][1], tableC[i][2]);
      }
      break;
    }
    case 6:
    {
      const real_t tableBonn[276][3]=
      {
        #include "sgStar_heff/sgStar_heff_Bonn.tab"
      };
      for(i=0;i<=275;i++) 
      {
        dataT[i]=DataWithTemperature_t(tableBonn[i][0], tableBonn[i][1], tableBonn[i][2]);
      }
      break;
    }
    default:
    {
      const real_t tableold[276][3]=
      {
        #include "sgStar_heff/sgStar_heff_old.tab"
      };
      for(i=0;i<=275;i++) 
      {
        dataT[i]=DataWithTemperature_t(tableold[i][0], tableold[i][1], tableold[i][2]);
      }
    }
  }

#ifdef DEBUG
  for( auto elem : dataT ) std::cout << elem;
#endif
  
  return;
}

void Relicparam_t::init_heffgeff(const int x)
{
  setQCDeosModel(x);
  init_heffgeff();
}

void Relicparam_t::init(const real_t &mrelic, const int x)
{
#ifdef DEBUG
  std::cout << "Calling init with relicmass value " << mrelic << std::endl;
#endif
  setQCDeosModel(x);
  init_heffgeff();
  
  // The rest of the default constructor is the correspondant of the function 
  //      void Init_cosmomodel(struct relicparam* paramrelic)
  // in SuperIso Relic V 4
  full_comput=0;
  
  solver=1; /* 1=logarithmic, 2=linear */
    
  failsafe=1;         // 0=fast, 1=precise (default), ... See stand_cosmo.c
  err=0;
  Tinit=27.;          // Starting at T = 27 x 10^9 K as default
  Tnudec=27.;           // Neutrino decoupling T = 27 x 10^9 K as default
  eta0=6.10e-10;      // Baryon-to-photon ratio (Planck 2015 results XIII)
  Nnu=3.046;          // Number of SM neutrinos, e+- reheating included
  dNnu=0.;            // Number of extra neutrino species (e.g. sterile neutrinos)
  life_neutron=880.2; // Neutron lifetime (PDG2018)
  life_neutron_error=1.0; // Neutron lifetime uncertainty (PDG2017)
  xinu1=0.; 
  xinu2=0.; 
  xinu3=0.; 
  beta_samples=50;  // how accurately to model n<-->p beta reactions

  b_cdm_ratio=0.02242/0.11933;   // current baryon to cold dark matter density ratio (Planck 2018 results VI)
  
  wimp=0;
  //m_chi=g_chi=SMC_wimp=selfConjugate=fermion=EM_coupled=neut_coupled=neuteq_coupled=0;
  m_chi=0;
  g_chi=0;
  SMC_wimp=0;
  selfConjugate=0;
  fermion=0;
  EM_coupled=0;
  neut_coupled=0;
  neuteq_coupled=0;

  fierz=0.;      // no fierz interface in the standard model 
  B_chi=0.;      // default is no branching to dark matter m_p < m_chi < m_n

  dd0=ndd=Tdend=Tddeq=0.;
  sd0=nsd=Tsend=0.;
  nt0=nnt=Tnend=0.;
  Sigmad0=nSigmad=TSigmadend=0.;
  Sigmarad0=nSigmarad=TSigmaradend=0.;
  coupd=0;
    
  mgravitino=0.;
  
  phi_model=0;
  eta_phi=Gamma_phi=n_phi=rhot_phi_Tmax=Tphi0=rhot_phi0=0.;

  entropy_model=1;
  energy_model=1;
  
//   mgravitino=relicmass=0.;
  mgravitino=0.;
  if(mrelic > 0.) relicmass=mrelic;
  
  full_comput=scalar=0;
  Tfo=0.;
  
  Tmax=10000.;
  
  chi2=0.;
  nobs=0;
    
  use_table_rhoPD=size_table_rhoPD=0;
}

Relicparam_t::Relicparam_t(const real_t &m, const int x)
{
#ifdef DEBUG
   std::cout << "Calling constructor with relicmass value " << m << std::endl;
#endif
   init(m,x);
}


void Relicparam_t::print_relicparam(std::ostream &out) const 
{
  out << "Printing the Relicparam_t structure:\n";
  out << (model_eff);
  out << (entropy_model);
  out << (energy_model);
  out << (dd0);
  out << (ndd);
  out << (Tdend);
  out << (Tddeq);
  out << (sd0);
  out << (nsd);
  out << (Tsend);
  out << (Sigmad0);
  out << (nSigmad);
  out << (TSigmadend);
  out << (Sigmarad0);
  out << (nSigmarad);
  out << (TSigmaradend);
  out << (nt0);
  out << (nnt);
  out << (Tnend);
  out << (coupd);
  out << (quintn2);
  out << (quintn3);
  out << (quintn4);
  out << (quintT12);
  out << (quintT23);
  out << (quintT34);
  out << (phi_model);
  out << (eta_phi);
  out << (Gamma_phi);
  out << (rhot_phi_Tmax);
  out << (n_phi);
  out << (rhot_phi0);
  out << (Tphi0);

  out << (T_RH);
  out << (Sigmatildestar);
  out << (Sigmatildestar_max);
  out << (Tstdstar_max);
  out << (mgravitino);
  out << (relicmass);
  out << (scalar);
  out << (solver);
  out << (beta_samples);
  out << (Tfo);
  out << (Tmax);
  out << (full_comput);
  out << (use_table_rhoPD);
  out << (size_table_rhoPD);
  for (int i = 0; i < size_table_rhoPD; ++i) 
  {
      out << (table_rhoPD[0][i]);
      out << (table_rhoPD[1][i]);
  }
  out << (err);
  out << (failsafe);
  out << (eta0);
  out << (Nnu);
  out << (dNnu);
  out << (life_neutron);
  out << (life_neutron_error);
  out << (xinu1);
  out << (xinu2);
  out << (xinu3);
  out << (m_chi);
  out << (g_chi);
  out << (Tinit);
  out << (Tnudec);
  out << (wimp);
  out << (SMC_wimp);
  out << (selfConjugate);
  out << (fermion);
  out << (EM_coupled);
  out << (neut_coupled);
  out << (neuteq_coupled);
  out << (chi2);
  out << (nobs);
  out << (fierz);
  out << (B_chi);
  out << (rhob0);
  out << (b_cdm_ratio);
  
}

real_t Relicparam_t::getheff(const real_t &T)
{
  if(T >= dataT.front().T) 
  {
#ifdef VERBOSE
      std::cerr << "Relicparam_t::getheff : you're trying to get heff for an too high temperature\n";
#endif
      return dataT.front().heff;
  }
  
  if( T <= dataT.back().T )
  {
#ifdef VERBOSE
      std::cerr << "Relicparam_t::getheff : you're trying to get heff for an too low temperature\n";
#endif
      return dataT.back().heff;
  }
  
  int ie = 1;
  
  while( T < dataT[ie].T ) ie++;
  
  const real_t logT=std::log(T);
  
  const real_t heff1=dataT[ie].heff;
  const real_t heff2=dataT[ie-1].heff;
  
  const real_t logT1=std::log(dataT[ie].T);
  const real_t logT2=std::log(dataT[ie-1].T);
      
  return (heff2-heff1)/(logT2-logT1)*(logT-logT1)+heff1;
}

real_t Relicparam_t::getsqrtgstar(const real_t &T)
{
  if(T >= dataT.front().T) 
  {
#ifdef VERBOSE
      std::cerr << "Relicparam_t::getsqrtgstar : you're trying to get sqrtgstar for an too high temperature\n";
#endif
      return dataT.front().sqrtgstar;
  }
  
  if( T <= dataT.back().T )
  {
#ifdef VERBOSE
      std::cerr << "Relicparam_t::getsqrtgstar : you're trying to get sqrtgstar for an too low temperature\n";
#endif
      return dataT.back().sqrtgstar;
  }
  
  int ie = 1;
  
  while( T < dataT[ie].T ) ie++;
  
  const real_t logT=std::log(T);
  
  const real_t sqrtgstar1=dataT[ie].sqrtgstar;
  const real_t sqrtgstar2=dataT[ie-1].sqrtgstar;
  
  const real_t logT1=std::log(dataT[ie].T);
  const real_t logT2=std::log(dataT[ie-1].T);
      
  return (sqrtgstar2-sqrtgstar1)/(logT2-logT1)*(logT-logT1)+sqrtgstar1;
}

real_t Relicparam_t::getgeff(const real_t &T)
{
  const real_t heff0 = getheff(T);
  
  return SQUARE(heff0/getsqrtgstar(T)*(1. + ( getheff(T*1.001) - getheff(T*0.999) )/0.006/heff0 ));
} 


#ifdef DEBUG
#undef DEBUG
#endif


void Relicparam_t::Init_cosmomodel_param(const real_t &eta, const real_t &Nnu_local, const real_t &dNnu_local, const real_t &life_neutron_local, const real_t &life_neutron_error_local, const real_t &xinu1_local, const real_t &xinu2_local, const real_t &xinu3_local)
/* modifies the values of the baryon-to-photon ratio eta, the number of SM neutrinos Nnu_local, extra neutrino species dNnu_local
 *  and the neutron lifetime life_neutron_local */
{
    this->eta0=eta;
    this->Nnu=Nnu_local;
    this->dNnu=dNnu_local;
    this->life_neutron=life_neutron_local;
    this->life_neutron_error=life_neutron_error_local;
    this->xinu1=xinu1_local;
    this->xinu2=xinu2_local;
    this->xinu3=xinu3_local;
    this->beta_samples=50;  // how accurately to model n<-->p beta reactions
    return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_wimp(const real_t &mass_wimp, int EM_coupled_local, int neut_coupled_local, int neuteq_coupled_local, int fermion_local, int selfConjugate_local, const real_t &g_chi_local)
/* modifies the parameters of an included light WIMP */
{
    this->m_chi=mass_wimp;
    this->g_chi=g_chi_local;
    this->fermion=fermion_local;
    this->EM_coupled=EM_coupled_local;
    this->neut_coupled=neut_coupled_local;
    this->neuteq_coupled=neuteq_coupled_local;
    this->wimp=1;
    this->selfConjugate=selfConjugate_local;
    return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_dark_density(const real_t &dd0_local, const real_t &ndd_local, const real_t &T_end)
{
  if(phi_model.get() != 0)
  {
    this->energy_model=0;
    return;
  }

  this->energy_model=1;
  this->dd0=dd0_local;
  this->ndd=ndd_local;
  this->Tdend=T_end;

  this->Tddeq=0.;

  this->use_table_rhoPD=this->size_table_rhoPD=0;

  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_dark_density2(const real_t &ndd_local, const real_t &Tddeq_local, const real_t &T_end)
{
  if(phi_model.get() != 0)
  {
    this->energy_model=0;
    return;
  }

  this->energy_model=2;
  this->ndd=ndd_local;
  this->Tddeq=Tddeq_local;
  this->Tdend=T_end;

  this->dd0=0.;

  this->use_table_rhoPD=this->size_table_rhoPD=0;
  
  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_dark_coupling(int coupD)
{
  this->coupd=coupD;

  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_quintessence(const real_t &T12, const real_t &n2, const real_t &T23, const real_t &n3, const real_t &T34, const real_t &n4)
{
  if(phi_model.get() != 0)
  {
    this->energy_model=0;
    return;
  }

  this->energy_model=3;
  this->dd0=0.;
  this->use_table_rhoPD=this->size_table_rhoPD=0;
  
  this->quintn2=n2;
  this->quintn3=n3;
  this->quintn4=n4;
  this->quintT12=T12;
  this->quintT23=T23;
  this->quintT34=T34;
  
  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_dark_entropy(const real_t &sd0_local, const real_t &nsd_local, const real_t &T_end)
{
  if(phi_model.get() != 0) return;

  this->sd0=sd0_local;
  this->nsd=nsd_local;
  this->Tsend=T_end;

  this->use_table_rhoPD=this->size_table_rhoPD=0;
  
  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_dark_entropySigmaD(const real_t &Sigmad0_local, const real_t &nSigmad_local, const real_t &T_end)
{
  if(phi_model.get() != 0) return;
  
  this->Sigmad0=Sigmad0_local;
  this->nSigmad=nSigmad_local;
  this->TSigmadend=T_end;

  this->use_table_rhoPD=this->size_table_rhoPD=0;
  
  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_entropySigmarad(const real_t &Sigmarad0_local, const real_t &nSigmarad_local, const real_t &T_end)
{
  if(phi_model.get() != 0) return;

  this->Sigmarad0=Sigmarad0_local;
  this->nSigmarad=nSigmarad_local;
  this->TSigmaradend=T_end;

  this->use_table_rhoPD=this->size_table_rhoPD=0;
  
  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_nonthermal(const real_t &nt0_local, const real_t &nnt_local, const real_t &T_end)
{
  if(phi_model.get() != 0) return;

  this->nt0=nt0_local;
  this->nnt=nnt_local;
  this->Tnend=T_end;

  this->use_table_rhoPD=this->size_table_rhoPD=0;
  
  return;
}
/*--------------------------------------------------------------*/

void Relicparam_t::Init_gravitino(const real_t &mgravitino_local)
{
  this->mgravitino=mgravitino_local;
  
  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_scalarfield(const real_t &rhotilde_phi_Tmax, const real_t &Tmax_local, const real_t &T_RH_local, const real_t &eta_phi_local, const real_t &n_phi_local)
{
  this->phi_model=1;
  this->full_comput=1;
  
  this->T_RH=T_RH_local;
  this->eta_phi=eta_phi_local;  /* = b/m_phi */
  this->Gamma_phi=std::sqrt(4.*std::pow(pi,3.)*getgeff(T_RH_local)/45.)*T_RH_local*T_RH_local/Mplanck;
  this->rhot_phi_Tmax=rhotilde_phi_Tmax;
  this->Tmax=Tmax_local;
  this->n_phi=n_phi_local;
  
  this->rhot_phi0=this->rhot_phi_Tmax;
  
  Init_dark_density(0.,0.,0.);
  Init_dark_density2(0.,0.,0.);
  Init_dark_entropy(0.,0.,0.);
  Init_dark_entropySigmaD(0.,0.,0.);
  Init_entropySigmarad(0.,0.,0.);
  Init_nonthermal(0.,0.,0.);
  Init_dark_coupling(0);
  
  this->use_table_rhoPD=this->size_table_rhoPD=0;

  return;
}

/*--------------------------------------------------------------*/

void Relicparam_t::Init_dark_density_table(real_t table[2][NTABMAX], int nlines)
{
  Init_dark_density(0.,0.,0.);
  Init_dark_density2(0.,0.,0.);
  Init_dark_entropy(0.,0.,0.);
  Init_dark_entropySigmaD(0.,0.,0.);
  Init_entropySigmarad(0.,0.,0.);
  Init_nonthermal(0.,0.,0.);
  
  this->use_table_rhoPD=1;
  
  int ie,je;
  
  this->size_table_rhoPD=nlines;
  for(ie=0;ie<=1;ie++) for(je=0;je<nlines;je++) this->table_rhoPD[ie][je]=table[ie][je];
    
  return;
}


/*--------------------------------------------------------------*/

void Relicparam_t::Init_neutron_decay(const real_t &tau, const real_t &tau_err, const real_t &fierz_local, const real_t &m_chi_local, const real_t &B_chi_local)
{
  this->life_neutron=tau;
  this->life_neutron_error=tau_err; // Neutron lifetime uncertainty (PDG2017)
  this->wimp=0;
  this->m_chi=m_chi_local;
  this->g_chi=0;
  this->SMC_wimp=0;
  this->selfConjugate=0;
  this->fermion=0;
  this->EM_coupled=0;
  this->neut_coupled=0;
  this->neuteq_coupled=0;
  this->fierz=fierz_local;      // no fierz_local interface in the standard model 
  this->B_chi=B_chi_local;      // default is no branching to dark matter m_p < m_chi_local < m_n
  this->beta_samples=1000;    // if we are gonna simulate beta decay, we're gonna simulate it.
}


/*--------------------------------------------------------------*/

real_t Relicparam_t::dark_density(const real_t &T)
{
  if(phi_model.get() != 0) return 0.;
  
  if(this->size_table_rhoPD>1&&this->use_table_rhoPD)
  {
    int ie=1;
  
    if(T<this->table_rhoPD[0][this->size_table_rhoPD-1])
    {
      ie=this->size_table_rhoPD-1;
    }
    else if(T<this->table_rhoPD[0][0]&&T>this->table_rhoPD[0][this->size_table_rhoPD-1])
    {
      while(T<this->table_rhoPD[0][ie]&&ie<this->size_table_rhoPD) ie++;
    }

    real_t logrhoD1,logrhoD2,rhoD,logT1,logT2,logT;
    logT=std::log(T);
    logrhoD1=std::log(this->table_rhoPD[1][ie]);
    logrhoD2=std::log(this->table_rhoPD[1][ie-1]);
    logT1=std::log(this->table_rhoPD[0][ie]);
    logT2=std::log(this->table_rhoPD[0][ie-1]);
    
    rhoD=exp((logrhoD2-logrhoD1)/(logT2-logT1)*(logT-logT1)+logrhoD1);
    
    return rhoD;
  }
  
  if(T<this->Tdend) return 0.;
  
  
  if(this->energy_model==3)
  {
    const real_t H0=67.8/3.0856e19; /* Hubble constant in second */  
    const real_t rho_Lambda=0.7*H0*H0/(8.*pi*Gn)/2.322e17; 
    
    if(T<=this->quintT12) return rho_Lambda;
    
    const real_t rho02=rho_Lambda;
    if(T<=this->quintT23) return rho02*std::pow(T/this->quintT12,this->quintn2);
    
    const real_t rho03=rho02*std::pow(this->quintT23/this->quintT12,this->quintn2);
    if(T<=this->quintT34) return rho03*std::pow(T/this->quintT23,this->quintn3);
    
    const real_t rho04=rho03*std::pow(this->quintT34/this->quintT23,this->quintn3);
    return rho04*std::pow(T/this->quintT34,this->quintn4);
  }
  
  if(this->energy_model==2)
  {
    if(this->Tddeq==0.) return 0.;
    
    const real_t geffT=getgeff(T);
    const real_t rhorad=pi*pi/30.*geffT*std::pow(T,4.);
    
    return rhorad*(getgeff(this->Tddeq)/geffT)*std::pow(getheff(T)/heff(this->Tddeq),this->ndd/3.)*std::pow(T/this->Tddeq,this->ndd);
  }
  
  if(this->energy_model==1)
  {
    if(this->dd0==0.) return 0.;
    
    const real_t rho_photon_1MeV=pi*pi/15.*1.e-12;
        
    return this->dd0*rho_photon_1MeV*std::pow(T/1.e-3,this->ndd);
  }
  
  return 0.;
}


real_t Relicparam_t::dark_density_pressure(const real_t &T)
{
  if(phi_model.get() != 0) return 0.;
    
  if(T<this->Tdend) return 0.;
  
  if(this->energy_model==1) if(this->dd0==0.) return 0.;
  
  if(this->energy_model==2) if(this->Tddeq==0.) return 0.;
  
  const real_t ddark_density_dT=(dark_density(T*1.001)-dark_density(T*0.999))/0.002/T;
  
  const real_t dentropy_dT=2.*pi*pi/45.*(getheff(T*1.001)*std::pow(T*1.001,3.)-heff(T*0.999)*std::pow(T*0.999,3.))/0.002/T+dark_entropy_derivative(T);
  
  const real_t entropy=2.*pi*pi/45.*getheff(T)*std::pow(T,3.)+dark_entropy(T);
  
  //return (this->ndd/3.-1.)*dark_density(T); /* outdated */
  return (ddark_density_dT-dentropy_dT/entropy*dark_density(T))*entropy/dentropy_dT;
}

/*--------------------------------------------------------------*/

real_t Relicparam_t::sigma_entropy(const real_t &T)
{
  if(phi_model.get() != 0) return 1.;

  if(this->Sigmarad0==0.) return 1.;
    
  real_t lnT,dlnT,Ttmp;
  int ie,nmax;
  real_t integ=0.;

  real_t heffT,geffT,darkdensitytilde,heffTdT,geffTdT,darkdensitytildeTdT,Htilde;
  real_t Sigmatildestar_local,dSigmatildestar_dT;
  
  nmax=10;
  
  lnT=std::log(1.e-15);
  
  dlnT=(std::log(T)-lnT)/nmax;
  
  for(ie=1;ie<nmax;ie++) 
  {
    lnT+=dlnT;
    Ttmp=std::exp(lnT);
    
    heffT=getheff(Ttmp);
    geffT=getgeff(Ttmp);
    darkdensitytilde=dark_density(Ttmp)/(pi*pi/30.*geffT*std::pow(Ttmp,4.));
    
    Htilde=std::sqrt(1.+darkdensitytilde); /*Htilde = H / std::sqrt(8 pi / 3 M_P^2) / rho_rad */
    
    Sigmatildestar_local=45.*std::sqrt(5.)/4./std::pow(pi,3.5)*Mplanck/heffT/std::sqrt(geffT)/std::pow(Ttmp,5.)/Htilde*entropy_Sigmarad(Ttmp);

    heffTdT=getheff(Ttmp*1.001);
    geffTdT=getgeff(Ttmp*1.001);
    darkdensitytildeTdT=dark_density(Ttmp*1.001)/(pi*pi/30.*geffTdT*std::pow(Ttmp*1.001,4.));
    
    dSigmatildestar_dT=((45.*std::sqrt(5.)/4./std::pow(pi,3.5)*Mplanck/heffTdT/std::sqrt(geffTdT)/std::pow(Ttmp*1.001,5.)/std::sqrt(1.+darkdensitytildeTdT)*entropy_Sigmarad(Ttmp*1.001))-Sigmatildestar_local)/0.001/Ttmp;
    
    integ+=-dSigmatildestar_dT/std::pow(1.-Sigmatildestar_local,2.)*std::log(heffT*std::pow(Ttmp,3.));
  }

  heffT=getheff(T);
  geffT=getgeff(T);
  darkdensitytilde=dark_density(T)/(pi*pi/30.*geffT*std::pow(T,4.));
  
  Htilde=std::sqrt(1.+darkdensitytilde); /*Htilde = H / std::sqrt(8 pi / 3 M_P^2) / rho_rad */
  
  Sigmatildestar_local=45.*std::sqrt(5.)/4./std::pow(pi,3.5)*Mplanck/heffT/std::sqrt(geffT)/std::pow(T,5.)/Htilde*entropy_Sigmarad(T);

  heffTdT=getheff(T*1.001);
  geffTdT=getgeff(T*1.001);
  darkdensitytildeTdT=dark_density(T*1.001)/(pi*pi/30.*geffTdT*std::pow(T*1.001,4.));
  
  dSigmatildestar_dT=((45.*std::sqrt(5.)/4./std::pow(pi,3.5)*Mplanck/heffTdT/std::sqrt(geffTdT)/std::pow(T*1.001,5.)/std::sqrt(1.+darkdensitytildeTdT)*entropy_Sigmarad(T*1.001))-Sigmatildestar_local)/0.001/T;
  
  integ+=-dSigmatildestar_dT/std::pow(1.-Sigmatildestar_local,2.)*std::log(heffT*std::pow(T,3.))/2.;
      
  integ*=dlnT;
  
  return std::exp(integ);
}


real_t Relicparam_t::dark_entropy(const real_t &T)
{
  if(phi_model.get() != 0) return 0.;

  if((this->sd0==0.)&&(this->Sigmad0==0.)) return 0.;
  
  if((this->Sigmad0==0.)&&(T<this->Tsend)) return 0.;
  
  if(this->Sigmad0==0.)
  {
    const real_t s_photon_1MeV=4.*pi*pi/45.*1.e-9;
  
    return this->sd0*s_photon_1MeV*std::pow(T/1.e-3,this->nsd);
  }
  else
  {
    real_t lnT,dlnT,Ttmp;
    int ie,nmax;
    real_t integ=0.;

    real_t heffT,geffT,darkdensitytilde,Htilde;
    real_t Sigmatildestar_local=0.;

    nmax=50;
    
    lnT=std::log(1.e-15);
    
    dlnT=(std::log(T)-lnT)/nmax;
    
    for(ie=1;ie<nmax;ie++) 
    {
      lnT+=dlnT;
      Ttmp=std::exp(lnT);
      
      heffT=getheff(Ttmp);
      geffT=getgeff(Ttmp);
      darkdensitytilde=dark_density(Ttmp)/(pi*pi/30.*geffT*std::pow(Ttmp,4.));
      
      Htilde=std::sqrt(1.+darkdensitytilde); /*Htilde = H / std::sqrt(8 pi / 3 M_P^2) / rho_rad */
      
      if(this->Sigmarad0!=0.) Sigmatildestar_local=45.*std::sqrt(5.)/4./std::pow(pi,3.5)*Mplanck/heffT/std::sqrt(geffT)/std::pow(Ttmp,5.)/Htilde*entropy_Sigmarad(Ttmp);
      
      integ+=getsqrtgstar(Ttmp)*dark_entropy_Sigmad(Ttmp)/Htilde/(1.-Sigmatildestar_local)/sigma_entropy(Ttmp)/std::pow(heffT*std::pow(Ttmp,3.),(2.-Sigmatildestar_local)/(1.-Sigmatildestar_local));
    }

    heffT=getheff(T);
    geffT=getgeff(T);
    darkdensitytilde=dark_density(T)/(pi*pi/30.*geffT*std::pow(T,4.));
    
    Htilde=std::sqrt(1.+darkdensitytilde); /*Htilde = H / std::sqrt(8 pi / 3 M_P^2) / rho_rad */
    
    if(this->Sigmarad0!=0.) Sigmatildestar_local=45.*std::sqrt(5.)/4./std::pow(pi,3.5)*Mplanck/heffT/std::sqrt(geffT)/std::pow(T,5.)/Htilde*entropy_Sigmarad(T);
    
    integ+=getsqrtgstar(T)*dark_entropy_Sigmad(T)/Htilde/(1.-Sigmatildestar_local)/sigma_entropy(T)/std::pow(heffT*std::pow(T,3.),(2.-Sigmatildestar_local)/(1.-Sigmatildestar_local))/2.;
        
    integ*=dlnT;

    return Mplanck*std::sqrt(45./4./std::pow(pi,3.))*std::pow(heffT*std::pow(T,3.),1./(1.-Sigmatildestar_local))*sigma_entropy(T)*integ;  
  }
}


real_t Relicparam_t::dark_entropy_derivative(const real_t &T)
{
  if(phi_model.get() != 0) return 0.;

  if((this->sd0==0.)&&(this->Sigmad0==0.)) return 0.;
  
  if((this->Sigmad0==0.)&&(T<this->Tsend)) return 0.;
  
  if(this->Sigmad0==0.)
  {
    return this->nsd*dark_entropy(T)/T;
  }
  else
  {
    const real_t heffT=getheff(T);
    const real_t geffT=getgeff(T);
//old unused variable
    const real_t rhoradT=pi*pi/30.*geffT*std::pow(T,4.);
    const real_t darkdensitytilde=dark_density(T)/rhoradT;
    const real_t Sigmarad=entropy_Sigmarad(T);
    
    const real_t Htilde=std::sqrt(1.+darkdensitytilde); /*Htilde = H / std::sqrt(8 pi / 3 M_P^2) / rho_rad */

    const real_t Sigmatildestar_local=45.*std::sqrt(5.)/4./std::pow(pi,3.5)*Mplanck/heffT/std::sqrt(geffT)/std::pow(T,5.)/Htilde*Sigmarad;

    return 3.*getsqrtgstar(T)/T/(1.-Sigmatildestar_local)/heffT
    *(std::sqrt(geffT)*dark_entropy(T)
    -std::sqrt(5.*Mplanck/4./std::pow(pi,3.))/T/T*dark_entropy_Sigmad(T)/Htilde);
  }
}


real_t Relicparam_t::dark_entropy_Sigmad(const real_t &T)
{
  if(phi_model.get() != 0) return 0.;

  if((this->sd0==0.)&&(this->Sigmad0==0.)) return 0.;
  
  if((this->Sigmad0==0.)&&(T<this->Tsend)) return 0.;
  
  if(this->Sigmad0==0.)
  {  
    const real_t heffT=getheff(T);
    const real_t geffT=getgeff(T);
//old unused variable
    const real_t rhoradT=pi*pi/30.*geffT*std::pow(T,4.);
    const real_t darkdensitytilde=dark_density(T)/rhoradT;
    const real_t Sigmarad=entropy_Sigmarad(T);
  
    const real_t Htilde=std::sqrt(1.+darkdensitytilde); /*Htilde = H / std::sqrt(8 pi / 3 M_P^2) / rho_rad */

    const real_t Sigmatildestar_local=45.*std::sqrt(5.)/4./std::pow(pi,3.5)*Mplanck/heffT/std::sqrt(geffT)/std::pow(T,5.)/Htilde*Sigmarad;

    return std::sqrt(4.*std::pow(pi,3.)/5.)/Mplanck*Htilde*T*T*(std::sqrt(geffT)*dark_entropy(T)-heffT/3./getsqrtgstar(T)*T*(1.-Sigmatildestar_local)*dark_entropy_derivative(T));
  }
  else
  {
    if(T<this->TSigmadend) return 0.;

    const real_t s_photon_1MeV=4.*pi*pi/45.*1.e-9;
  
    const real_t Sigma_photon_1MeV= 1./Mplanck*std::sqrt(8.*pi*pi*pi/5.)*(1.e-6)*s_photon_1MeV;
  
    return this->Sigmad0*Sigma_photon_1MeV*std::pow(T/1.e-3,this->nSigmad);
  }
}

/*--------------------------------------------------------------*/

real_t Relicparam_t::entropy_Sigmarad(const real_t &T)
{
  if( (this->phi_model).get() != 0) return 0.; //this->Gamma_phi*this->rho_phi/T;

  if((this->Sigmarad0==0.)||(T<this->TSigmaradend)) return 0.;
  
  const real_t s_photon_1MeV=4.*pi*pi/45.*1.e-9;
  
  const real_t Sigma_photon_1MeV=1./Mplanck*std::sqrt(8.*pi*pi*pi/5.)*(1.e-6)*s_photon_1MeV;
  
  return this->Sigmarad0*Sigma_photon_1MeV*std::pow(T/1.e-3,this->nSigmarad);
}

/*--------------------------------------------------------------*/

real_t Relicparam_t::nonthermal(const real_t &T)
{
  if(phi_model.get() != 0) return 0.; //this->eta_phi*this->Gamma_phi*this->rho_phi;

  if((this->nt0==0.)||(T<this->Tnend)) return 0.;
  
  return this->nt0*1.e-50*std::pow(T/1.e-3,this->nnt);
}

/*--------------------------------------------------------------*/

real_t Relicparam_t::neutdens(const real_t &Tnu)
/* Computes the neutrino density, including any effects from a neutrino degeneracy */
{
    if((this->xinu1==0.)&&(this->xinu2==0.)&&(this->xinu3==0.))
    {
        /* No degeneracy, relativistic approximation */
        return 2.*pi*pi/30.*7./8.*this->Nnu*std::pow(Tnu,4.);
    }

    int ie,je,n;
    real_t rho=0.;
    real_t xinu[4];
    real_t max1,max2,int1,int2;
    real_t x;

    xinu[1]=this->xinu1;
    xinu[2]=this->xinu2;
    xinu[3]=this->xinu3;

    /* SM neutrinos */
    for(ie=1;ie<=3;ie++)
    {
        /* The factor (this->Nnu/3.) includes extra DOF from non-rel. e+- and non-inst. nu decoupling */
        if(std::fabs(xinu[ie])<=0.03)
        {
            rho+=(this->Nnu/3.)*2.*pi*pi/30.*std::pow(Tnu,4.)*(7./8.+(15./(4*pi*pi))*xinu[ie]*xinu[ie]+
                                                               (15./(8.*std::pow(pi,4.)))*std::pow(xinu[ie],4.));
        }
        else if(std::fabs(xinu[ie])>=30.)
        {
            rho+=(this->Nnu/3.)*std::pow(Tnu,4.)/(8.*pi*pi)*std::pow(xinu[ie],4.)*(1.+12.*1.645/xinu[ie]/xinu[ie]);
        }
        else
        {
            /* Neutrinos */
            max1=(88.029+xinu[ie])*Tnu;
            int1=0.;
            n=50;
            for(je=1;je<=n-1;je++)
            {
                x=(real_t)je/(real_t)n*max1;
                int1+=1./(2.*pi*pi)*std::pow(x,3.)/(1.+std::exp(x/Tnu-xinu[ie]));
            }
            int1+=0.5*1./(2.*pi*pi)*std::pow(max1,3.)/(1.+std::exp(max1/Tnu-xinu[ie]));
            int1*=(this->Nnu/3.)*max1/(real_t)n;
            rho+=int1;

            /* Anti-neutrinos */
            max2=(88.029-xinu[ie])*Tnu;
            if(max2>0.)
            {
                int2=0.;
                n=50;
                for(je=1;je<=n-1;je++)
                {
                    x=(real_t)je/(real_t)n*max2;
                    int2+=1./(2.*pi*pi)*std::pow(x,3.)/(1.+std::exp(x/Tnu+xinu[ie]));
                }
                int2+=0.5/(2.*pi*pi)*std::pow(max2,3.)/(1.+std::exp(max2/Tnu+xinu[ie]));
                int2*=(this->Nnu/3.)*max2/(real_t)n;
                rho+=int2;
            }
        }
    }
    return rho;
}

/*--------------------------------------------------------------*/

real_t Relicparam_t::neutdens_deriv(const real_t &Tnu)
/* Computes the temperature (Tnu) derivative of the neutrino energy density */
{
    if((this->xinu1==0.)&&(this->xinu2==0.)&&(this->xinu3==0.))
    {
        return 7.*pi*pi/30.*this->Nnu*std::pow(Tnu,3.);
    }

    int ie,je,n;
    real_t drho=0.;
    real_t xinu[4];
    real_t max1,max2,int1,int2;
    real_t x;

    xinu[1]=this->xinu1;
    xinu[2]=this->xinu2;
    xinu[3]=this->xinu3;

    /* SM neutrinos */
    for(ie=1;ie<=3;ie++)
    {
        if(std::fabs(xinu[ie])<=0.03)
        {
            drho+=(this->Nnu/3.)*4.*pi*pi/15.*std::pow(Tnu,3.)*(7./8.+(15./(4*pi*pi))*xinu[ie]*xinu[ie]
                                                                +(15./(8.*std::pow(pi,4.)))*std::pow(xinu[ie],4.));
        }
        else if(std::fabs(xinu[ie])>=30.)
        {
            drho+=(this->Nnu/3.)*std::pow(Tnu,3.)/(2.*pi*pi)*std::pow(xinu[ie],4.)*(1.+12.*1.645/xinu[ie]/xinu[ie]);
        }
        else
        {
            max1=(88.029+xinu[ie])*Tnu;
            int1=0.;
            n=50;
            for(je=1;je<=n-1;je++)
            {
                x=(real_t)je/(real_t)n*max1;
                int1+=1./(2.*pi*pi)*std::pow(x,3.)/(1.+std::exp(x/Tnu-xinu[ie]));
            }
            int1+=0.5*1./(2.*pi*pi)*std::pow(max1,3.)/(1.+std::exp(max1/Tnu-xinu[ie]));
            int1*=(this->Nnu/3.)*4.*max1/Tnu/(real_t)n;
            drho+=int1;

            max2=(88.029-xinu[ie])*Tnu;
            if(max2>0.)
            {
                int2=0.;
                n=50;
                for(je=1;je<=n-1;je++)
                {
                    x=(real_t)je/(real_t)n*max2;
                    int2+=1./(2.*pi*pi)*std::pow(x,3.)/(1.+std::exp(x/Tnu+xinu[ie]));
                }
                int2+=0.5/(2.*pi*pi)*std::pow(max2,3.)/(1.+std::exp(max2/Tnu+xinu[ie]));
                int2*=(this->Nnu/3.)*4.*max2/Tnu/(real_t)n;
                drho+=int2;
            }
        }
    }
    return drho;
}

/*--------------------------------------------------------------*/

real_t neutN(const real_t &T)
/* Computes the round N(z) function of the neutrinos - Pisanti et al., 0705.0290, eq. (A24) */
{
  const real_t z=m_e/T;

  if(z>=4.) return 0.;

  real_t logNz=0.;
  
  const std::array<real_t,14> n={-10.21703221236002,61.24438067531452,-340.3323864212157,1057.2707914654834,-2045.577491331372,2605.9087171012848,-2266.1521815470196,1374.2623075963388,-586.0618273295763,174.87532902234145,-35.715878215468045,4.7538967685808755,-0.3713438862054167,0.012908416591272199};

  for(int ie=0;ie<14;ie++) logNz+=n[ie]*std::pow(z,ie);

  return std::exp(logNz);
}

} // end of namespace scalar2to2
