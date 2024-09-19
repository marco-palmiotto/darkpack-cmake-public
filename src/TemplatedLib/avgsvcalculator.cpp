#include "avgsvcalculator.hpp"
#define DEFAULT_T_LIMIT (1.0e-5*input.getLightestBSMmass())

#define OUTPATH "/home/palmiotto/Documents/MartyWorkspace/mssm-dev/RootWorkSpace/"

template <typename T>
static inline T SQUARE(const T x) {return x*x;}

template <typename T>
static inline T sq(const T x) {return x*x;}

static inline unsigned int kdelta(const int i, const int j) 
{
  return static_cast<unsigned int>(i==j);
}

namespace __SPEC_LIB_NAME__
{
  
inline real_t getContributionSingle(Process2to2 p, 
                       Param_t input,
                       const real_t Ecm, const real_t ctheta)
{
  return  p.get_process_dependent_contrib_dweff_dcos(input, Ecm, ctheta);
}

inline real_t getContribution(const std::shared_ptr<std::vector<Process2to2>> &vp, 
                       Param_t input,
                       const real_t Ecm, const real_t ctheta, const bool cut=false)
{
  // This function calculates the sum of all the process-dependent contributions to dWeff/dcostheta,
  // with no parallelisation, and taking into account possible cuts
  
  real_t result = 0.;
#ifdef DEBUG 
  size_t counter = 0; //Defining a counter for the number of processes
#endif
  for(auto p : *vp) 
  {
    if( cut 
      && p.getMass(1, input) + p.getMass(2,input) > 4*input.getLightestBSMmass() ) continue;// Ending cuts on Weff 
    result += p.get_process_dependent_contrib_dweff_dcos(input, Ecm, ctheta);
#ifdef DEBUG 
    counter++;
#endif
  }
  
#ifdef DEBUG 
  if(!cut && counter != (*vp).size()) std::cout << "Incorrect number of processes\n";
#endif
  
  return result;
}

void AvgSvCalculator::sort_bsm_spectrum()
{
  bsm_particles_sorted=corr::bsm_particles;
  
  std::sort(bsm_particles_sorted.begin(), bsm_particles_sorted.end(), 
            [&](const int a, const int b){ 
              return this->input.masses_vector[a]<this->input.masses_vector[b];
            });

#ifdef DEBUG
  std::cout << "bsm_particles_sorted is sorted:\n";
  for(size_t i = 0; i < bsm_particles_sorted.size(); i++)
    std::cout << i << " " << this->input.masses_vector[bsm_particles_sorted[i]] << " " <<
     corr::part_names[bsm_particles_sorted[i]] << "\n";
#endif
  
}

void AvgSvCalculator::init_vector(std::shared_ptr<std::vector<Process2to2>> &p_in)
{
#ifdef DEBUG
  std::cout << "AvgSvCalculator::init_vector: beginning\n";
#endif
  p_in->reserve(10000);
  
  for(auto i : corr::bsm_particles )
  {
    for(auto j : corr::bsm_particles )
    {
      for(auto k : corr::sm_particles )
      {
        for(auto l : corr::sm_particles  )
        {
          std::array<int,4> fields_list={i, j, k, l};
          for(bool bi : {true, false})
          {
            for(bool bj : {true, false})
            {
              for(bool bk : {true, false})
              {
                for(bool bl : {true, false})
                {                  
                  std::array<bool,4> anti_part_list={bi, bj, bk, bl};
                  p_in->emplace_back(fields_list, anti_part_list);     
                  if( !(p_in->back()).checkExistance() ) 
                  {
#ifdef DEBUG
                    std::cout << "The process does not exist, let's discard it...";
#endif
                    p_in->pop_back(); 
#ifdef DEBUG
                    std::cout << " done!\n";
#endif
                    continue;
                  }
#ifdef DEBUG
                  std::cout << (p_in->back()).getName() << " exists\nVerifying no process with the same key is listed\n";
#endif
                  bool isnew=true;
                  int pcount=p_in->size()-1;
                  
                  while(--pcount >=0 && isnew)
                  {
#ifdef DEBUG
                    std::cout << "count = " << pcount << std::endl;
                    std::cout << "process in vector is " << (*p_in)[pcount].getName() << std::endl;
                    const std::string key1=(p_in->back()).getKey();
                    const std::string key2=(*p_in)[pcount].getKey();
                    std::cout << "Evaluating isnew" << std::endl;
#endif
                    isnew=(p_in->back()).getKey()!= (*p_in)[pcount].getKey();
                  };
                  
                  if( !isnew ) 
                  { 
#ifdef DEBUG
                    std::cout << "The process is not new, let's discard it...";
#endif
                    p_in->pop_back(); 
#ifdef DEBUG
                    std::cout << " done!\n";
#endif
                    continue;
                  }
#ifdef DEBUG
                  std::cout << "The process is new, let's keep it... ";
#endif
                  (p_in->back()).setRunningData(&run);
                  (p_in->back()).setRunningExternal();
#ifdef DEBUG
                  std::cout << " done!\n";
#endif
                }
              }
            }
          }
        }
      }
    }
  }
#ifdef DEBUG
  std::cout << "Calling shrink_to_fit()\n";
#endif
  p_in->shrink_to_fit();
}

AvgSvCalculator::AvgSvCalculator(const Param_t &input_ext)
{
  input = input_ext;
  run.init(input);
  run.setHiggsLoopFactor(true);
  
  runningenabled=true;
  
  Weffcuts=false;
  p_ptr=std::make_shared<std::vector<Process2to2>>();
  T_lim_sigmav=DEFAULT_T_LIMIT;
  orderT=1;
  sort_bsm_spectrum();
  init_vector(p_ptr);
  wastherenosplitting=false;
  werealldeltaxinull=false;
}

AvgSvCalculator::AvgSvCalculator(const Param_t &input_ext, std::shared_ptr<std::vector<Process2to2>> p_ext)
{
  input = input_ext;
  run.init(input);
  run.setHiggsLoopFactor(true);
  
  runningenabled=true;
  
  Weffcuts=false;
  sort_bsm_spectrum();
  p_ptr=p_ext;
  orderT=1;
  T_lim_sigmav=DEFAULT_T_LIMIT;
  wastherenosplitting=false;
  werealldeltaxinull=false;
}

void AvgSvCalculator::changeInput(const Param_t &input_ext, const bool hasMassSpetrumChanged)
{
#ifdef DEBUG
  std::cout << "Entering AvgSvCalculator::changeInput\n";
#endif
  input = input_ext;
  run.init(input);
  run.setHiggsLoopFactor(true);
  
  if(hasMassSpetrumChanged )
  {
    sort_bsm_spectrum();
    
    g2_Wefftab.clear();
    sqrtStab.clear();
    pefftab.clear();
    
    setLimitingValues();
  }else
  {
    setLimitingValues();
    
    if( g2_Wefftab.size() != 0)
    {
      g2_Wefftab[0]=Weffderiv.at(0)*sq(input.getLightestBSMdof());
      for(size_t i = 1 ; i < sqrtStab.size() ; i++)
      {
        g2_Wefftab[i]= get_g2_Weff(sqrtStab[i]);
      }
    }
    
  }
  
  setTaylorCoeffSVx_nosplitting();
}

void AvgSvCalculator::setProcessList(std::shared_ptr<std::vector<Process2to2>> p_ext)
{
  p_ptr=p_ext;
}

bool AvgSvCalculator::operator== (const AvgSvCalculator &other) const
{
  bool temp=true;
  
  if(this->input!=other.input) 
  { 
    std::cout << "inputs are different\n";
    temp=false;
  }
  // run != other.run;
  if(this->p_ptr!=other.p_ptr) 
  {
    std::cout << "p_ptr are different\n";
    temp=false;
  }

  if(this->Weffcuts!=other.Weffcuts) 
  {
    std::cout << "Weffcuts are different\n";
    temp=false;
  }
  if(this->sqrtStab!=other.sqrtStab)
  {
    std::cout << "sqrtStab are different\n";
    temp=false;
  }
  if(this->pefftab!=other.pefftab) 
  {
    std::cout << "pefftab are different\n";
    temp=false;
  }
  if(this->g2_Wefftab!=other.g2_Wefftab)
  {
    std::cout << "g2_Wefftab are different\n";
    temp=false;
  }
  
  if(this->T_lim_sigmav != other.T_lim_sigmav)
  {
      std::cout << "T_lim_sigmav are different\n";
      temp = false;
  }
  
  if(this->orderT != other.orderT)
  {
      std::cout << "orderT are different\n";
      temp = false;
  }
  if(this->bsm_particles_sorted != other.bsm_particles_sorted)
  {
      std::cout << "bsm_particles_sorted are different\n";
      temp = false;
  }
  if(this->ytab != other.ytab)
  {
      std::cout << "ytab are different\n";
      temp = false;
  }
  if(this->Weffderiv != other.Weffderiv)
  {
      std::cout << "Weffderiv are different\n";
      temp = false;
  }
  if(this->Weffderiv_err != other.Weffderiv_err)
  {
      std::cout << "Weffderiv_err are different\n";
      temp = false;
  }
  if(this->TaylorCoeffSVT != other.TaylorCoeffSVT)
  {
      std::cout << "TaylorCoeffSVT are different\n";
      temp = false;
  }
  if(this->wastherenosplitting != other.wastherenosplitting)
  {
      std::cout << "wastherenosplitting are different\n";
      temp = false;
  }
  if(this->werealldeltaxinull != other.werealldeltaxinull)
  {
      std::cout << "werealldeltaxinull are different\n";
      temp = false;
  }
  if(this->runningenabled != other.runningenabled)
  {
      std::cout << "runningenabled are different\n";
      temp = false;
  }
  
// std::cin.get();
  return temp;
}

AvgSvCalculator::AvgSvCalculator(const AvgSvCalculator &other)
 : input(other.input),
   run(other.run),
   p_ptr(other.p_ptr),
   Weffcuts(other.Weffcuts),
   sqrtStab(other.sqrtStab), 
   pefftab(other.pefftab), 
   g2_Wefftab(other.g2_Wefftab),
   T_lim_sigmav(other.T_lim_sigmav),
   orderT(other.orderT),
   wastherenosplitting(false) ,
   werealldeltaxinull(false)
   /*
   TaylorCoeffSVT(other.TaylorCoeffSVT)
   ytab(other.ytab),
   Weffderiv(other.Weffderiv),
   Weffderiv_err(other.Weffderiv_err),
   bsm_particles_sorted(other.bsm_particles_sorted)
   */
   ,runningenabled(other.runningenabled)
{
   this->sort_bsm_spectrum();
}

AvgSvCalculator& AvgSvCalculator::operator=(const AvgSvCalculator &other)
{
  if( this == &other ) return *this;
  
  input=other.input;
  run=other.run;
  p_ptr=other.p_ptr;
  Weffcuts=other.Weffcuts;
  sqrtStab=other.sqrtStab; 
  pefftab=other.pefftab; 
  g2_Wefftab=other.g2_Wefftab;
  T_lim_sigmav=other.T_lim_sigmav;
  orderT=other.orderT;
  wastherenosplitting=false;
  werealldeltaxinull=false;
  this->sort_bsm_spectrum();

  /*
  ytab=other.ytab;
  Weffderiv=other.Weffderiv;
  Weffderiv_err=other.Weffderiv_err;
  TaylorCoeffSVT=other.TaylorCoeffSVT;
  */
  
  runningenabled=other.runningenabled;
  run.init(input);
  
  return *this;
}


void AvgSvCalculator::print_procs(std::ostream& out, const bool verbose) const
{
  out << "There are " << p_ptr->size() +1 << " processes in the list:\n";
  for(size_t i = 0 ; i< p_ptr->size() ; i++)
  {
    out << i+1 << " : ";
    if(verbose) out << (*p_ptr)[i];
    else out << (*p_ptr)[i].getName();
    out << '\n';
  }  
}

real_t AvgSvCalculator::get_g2_dWeff_dcos(const real_t & sqrtS, const real_t & cosine)
{
// This function calculates [g^2_LSP * dWeff/dcostheta]
#ifdef DEBUG
  std::cout << "Calling method for running\n";
#endif
  runAtScale(sqrtS);
#ifdef DEBUG
  std::cout << "Running performed\n" <<
               "Launching threads for computing each process' contribution\n";
  std::cout << "inputs are\n"; input.Print();
#endif
  std::vector<std::future<real_t>> vres;// I am obliged to use a dynamic vector because there's no way of knowing,
                                        // at compile time, how many processes we have to store
                                        // IDEA: Instantiating a static array of size 100 and reusing it?
  vres.reserve(p_ptr->size());
  
  for(auto proc : *p_ptr ) 
  {
    if(Weffcuts
      &&proc.getMass(1, input) + proc.getMass(2, input) > 4*getMassLBSM() ) continue;
    
    vres.emplace_back(std::async(std::launch::async, &getContributionSingle, proc, input, sqrtS, cosine));
  }
  real_t result = 0.;
  const real_t peff= 0.5*std::sqrt(SQUARE(sqrtS) - 4*SQUARE(getMassLBSM()) );
  const real_t prefac = 1./peff/sqrtS/(8.0*M_PI);
  for( size_t i = 0; i < vres.size() ; i++ )
  {
#ifndef DEBUG
    result += vres[i].get();
#else
    real_t temp= vres[i].get();
    std::cout << "Contribution from " << (*p_ptr)[i].getName() << " to dWeff/dcos("<< sqrtS<< ")="
              << prefac/SQUARE(input.getLightestBSMdof()) << " * " << temp << " = " << 
              prefac/SQUARE(input.getLightestBSMdof())*temp << '\n';
    result += temp;
#endif
  }
  
  return result*prefac;
}

// This function calculates [ g^2_LSP * Weff]
real_t AvgSvCalculator::get_g2_Weff(const real_t & sqrtS)
{
#ifdef DEBUG
  std::cout << "Calling method for running\n";
#endif
  runAtScale(sqrtS);
#ifdef DEBUG
  std::cout << "Running performed\n";
#endif
  
  auto lambda = [=](real_t costheta) mutable
  {
    return getContribution( p_ptr,input, sqrtS, costheta , Weffcuts);
  };
  
  real_t discr;
  
  real_t weff = advmath::integrate_gauss_comparative(-1., 1., lambda, 5.0e-3, &discr);
  
  if(discr > 5.0e-3) weff=advmath::integrate_trap(-1., 1., lambda, 5.0e-3, &discr);
  
  const real_t peff= 0.5*std::sqrt(SQUARE(sqrtS) - 4*SQUARE(getMassLBSM()) );
  const real_t prefac = 1./peff/sqrtS/(8.0*M_PI);
  
  return prefac*weff;
}

real_t AvgSvCalculator::getdWeff_dcos(const real_t & sqrtS, const real_t & cosine)
{
// This function calculates [ dWeff/dcostheta]
  return get_g2_dWeff_dcos(sqrtS, cosine)/SQUARE(input.getLightestBSMdof());
}

real_t AvgSvCalculator::getWeff(const real_t & sqrtS)
{
// This function calculates [ Weff]
  return get_g2_Weff(sqrtS)/SQUARE(input.getLightestBSMdof());
}

real_t AvgSvCalculator::getdWeff_dcos_fromW12(const real_t & sqrtS, const real_t & cosine)
{
// TESTING function
  const real_t peff =  0.5*std::sqrt(SQUARE(sqrtS) - 4*SQUARE(getMassLBSM()) );
#ifdef DEBUG
  std::cout << "Initialised peff\n";
#endif
  const real_t g1 = input.getLightestBSMdof();
  real_t sum = 0.;
  
  for(auto proc : *p_ptr ) 
  {
    if(Weffcuts
      &&proc.getMass(1, input) + proc.getMass(2, input) > 4*getMassLBSM() ) continue;
#ifdef DEBUG
    std::cout << "Computing the contribution of " << proc << std::endl;
#endif
    sum += proc.getDiffW12Contrib(input, sqrtS, cosine, true);
     
  }
   return 1./peff/SQUARE(g1)*sum;
}

real_t AvgSvCalculator::getWeff_fromW12(const real_t & sqrtS)
{
// TESTING function
// This method returs Weff computed from the integration getdWeff_dcos_fromW12, with 
// no optimisations
#ifdef DEBUG
  std::cout << "Calling method for running\n";
#endif
  runAtScale(sqrtS);
#ifdef DEBUG
  std::cout << "Running performed\n";
#endif
  
  auto lambda = [=](real_t costheta) mutable
  {
    return getdWeff_dcos_fromW12(sqrtS, costheta);
  };
  
  real_t discr;
  
  real_t weff = advmath::integrate_gauss_comparative(-1., 1., lambda, 5.0e-3, &discr);
  
  if(discr < 5.0e-3) weff=advmath::integrate_trap(-1., 1., lambda, 5.0e-3, &discr);
  
  return weff;
}

void AvgSvCalculator::setWeffcuts(bool x)
{
  Weffcuts = x;
}

void AvgSvCalculator::tabulateValues(const real_t &sqrtSmax, const size_t &Nmax)
{
  /* Computes a table of Weff*gDM^2 for different values of sqrtS */
  
#ifdef DEBUG
  std::cout << "Resetting vectors\n";
#endif
  sqrtStab.clear();
  pefftab.clear();
  g2_Wefftab.clear();
  
  sqrtStab.reserve(2*Nmax);
  pefftab.reserve(2*Nmax);
  g2_Wefftab.reserve(2*Nmax);
  

//Normal algorithm
//   Timer timer;
#ifdef DEBUG
#ifdef OUTPATH
  std::string nametab = static_cast<std::string>(OUTPATH)+static_cast<std::string>("weff_tab_new.dat");
#endif
#endif
  
  auto no_duplicates =  [] (const real_t x, const std::vector<real_t>::iterator itbegin, const std::vector<real_t>::iterator itend)
  {
    for( std::vector<real_t>::iterator it = itbegin ; it != itend ; it++ )
      if(advmath::getdiscrepancy(x, *it) < 1.0e-3) return false;
    return true;
  };
  
  // Note that sqrtSmin^2 has to be sqrtSmin^2 = 4*m_DM^2
  const real_t sqrtSmin = 2.*getMassLBSM();
  
#ifdef DEBUG  
  std::cout << "Adding mass values, sum of pair of masses values, mass +/- width/2\n";
#endif
  for( auto i = 1 ; i <= corr::TOTAL_PARTICLES ; i++ ) 
  {
#ifdef DEBUG  
  std::cout << " i = " << i << " part is " << corr::part_names[i] << std::endl;
  std::cout << "  temp = input.masses_vector[i]\n";
#endif
    real_t temp = input.masses_vector[i];
    if(  temp < sqrtSmax 
      && temp > sqrtSmin 
      && no_duplicates(temp, sqrtStab.begin(), sqrtStab.end() )) 
        sqrtStab.push_back(temp);
        
    // Adding sum of the each pair (m_i, m_j) 
    for( auto j = i ; j <= corr::TOTAL_PARTICLES ; j++)
    {
#ifdef DEBUG  
      std::cout << "  temp = input.masses_vector[i]+ input.masses_vector[j]\n";
#endif
      temp = input.masses_vector[i]+input.masses_vector[j];
      if(  temp < sqrtSmax 
        && temp > sqrtSmin 
        && no_duplicates(temp, sqrtStab.begin(), sqrtStab.end() )  ) 
          sqrtStab.push_back(temp);
    }
#ifdef DEBUG  
    std::cout << "  Neglecting very small widths\n";
#endif    
    if(input.widths_vector[i] < 1.0e-6) continue;  
#ifdef DEBUG  
    std::cout << "  temp = input.masses_vector[i] + 0.5*input.widths_vector[i];\n";
#endif   
    temp = input.masses_vector[i] + 0.5*input.widths_vector[i];  
    if(temp < sqrtSmax && temp > sqrtSmin && no_duplicates(temp, sqrtStab.begin(), sqrtStab.end() )  ) sqrtStab.push_back(temp);  
#ifdef DEBUG  
    std::cout << "  temp = input.masses_vector[i] - 0.5*input.widths_vector[i];\n";
#endif
    temp = input.masses_vector[i] - 0.5*input.widths_vector[i];   
    if(temp < sqrtSmax && temp > sqrtSmin && no_duplicates(temp, sqrtStab.begin(), sqrtStab.end() )  ) sqrtStab.push_back(temp); 
  } 
  
  sqrtStab.push_back(sqrtSmax);
//   sqrtStab.push_back(sqrtSmin);
  sqrtStab.push_back(sqrtSmin*1.0001);
#ifdef DEBUG  
  std::cout << "Adding: the 0.99, the 1.01 points, the halfpoints\n";
#endif
  size_t temp_size = sqrtStab.size();
  for( size_t i = 0 ; i < temp_size ; i++)
  {
    real_t temp = 0.99*sqrtStab[i];
    if(temp < sqrtSmax && temp > sqrtSmin && no_duplicates(temp, sqrtStab.begin(), sqrtStab.end() )  ) sqrtStab.push_back(temp);
    
    temp=1.01*sqrtStab[i];
    if(temp < sqrtSmax && temp > sqrtSmin && no_duplicates(temp, sqrtStab.begin(), sqrtStab.end() )  ) sqrtStab.push_back(temp);
    
    for( size_t j = i; j < temp_size ; j++)
    {
      temp = 0.5*( sqrtStab[i] + sqrtStab[j] );
      if(temp < sqrtSmax && temp > sqrtSmin && no_duplicates(temp, sqrtStab.begin(), sqrtStab.end() )  ) sqrtStab.push_back(temp);
    }
  }
  
#ifdef DEBUG  
  std::cout << "Sorting values\n";
#endif
  std::sort(sqrtStab.begin(), sqrtStab.end());
  
#ifdef DEBUG  
  std::cout << "Computing the relative peff and Weff for the default points\n";
#endif 
  for( real_t E : sqrtStab )
  {
    real_t temp = SQUARE(E) - SQUARE(sqrtSmin);
    if(temp < 0.)
    {
      std::cerr<< "Peff value for energy " << E << " is imaginary\n";
      exit(2);
    }
    temp = 0.5*std::sqrt(temp);
    if( pefftab.size() != 0 && temp < pefftab.back() ) 
    {
      std::cerr<< "Peff value "<< temp << " for energy " << E << " is lower than previus\n";
      exit(2);
    }
    // peff = 0.5 sqrt( s - 4*mlsp ) => 2 peff = sqrt () ; 4 peff2 = s - 4mlsp2 ; s = 4(peff2 + mlsp2)
    
    real_t temp2 = get_g2_Weff(E);
    if(temp2 < 0.)
    {
      std::cerr<< "Weff value for energy " << E << " is imaginary\n";
      exit(2);
    }
    
    pefftab.push_back(temp);
    g2_Wefftab.push_back(temp2);
  }
  
#ifdef DEBUG  
  for(size_t i = 0 ; i < pefftab.size() ; i++) std::cout << "i = " << i << " peff = " << pefftab[i] << std::endl;
  std::cout << "Adding mid-point values if needed\n";
#endif 
  temp_size = sqrtStab.size();
  
  if(temp_size < Nmax) 
  {
  #ifdef DEBUG 
    std::cout << "Adding mid-point values: list conversion\n";
  #endif 
    //   Let's use lists to be more efficient, then we'll re-convert to vectors
    std::list< std::pair<real_t,real_t> > list;
    
    
    #ifdef DEBUG  
      std::cout << "Adding to the list the limiting value\n";
    #endif 
    
    if(Weffderiv.size() ==0)
    {
      setLimitingValues();
    }
    
    list.push_back({0., Weffderiv.front()*sq(input.getLightestBSMdof())});
    
    
    for( size_t i = 0 ; i < pefftab.size() ; i++) list.push_back({pefftab[i], g2_Wefftab[i]});
    
    std::list< std::pair<real_t,real_t> >::iterator next, prev, current;

    prev= list.begin();
    current = std::next(prev);
    next=std::next(current);
    
    while( temp_size < Nmax && next != list.end())
    {
      // Using lambda values because they need to be evaluated each time 
      // With primitive types we would need an update after each change
      const real_t peffprev = prev->first ;
      auto peff = [&current]() { return current->first; };
      auto peffnext = [&next] () { return next->first; };

      const real_t weffprev = prev->second;
      auto weff = [&current]() { return current->second; };
      auto weffnext = [&next] () { return next->second; };
      
      auto cond1 = [&] () {
        return peff() - peffprev > pefftab.back()/Nmax;
//         return peff()- peffprev > 0.;
      };
      auto cond2 = [&] () 
      {
        return 
        std::abs(1. - (weffnext() - weff() ) / ( peffnext() - peff() ) * ( peff() - peffprev ) / ( weff() - weffprev ) ) > 0.05;
      };
      
      auto cond3 = [&] () {
        return peff()- peffprev > pefftab.back()/10./Nmax;
//         return peff()- peffprev > 0.;
      };
      auto cond4 = [&] () 
      {
        return
        std::abs((weffnext() - weff() ) / ( peffnext() - peff() )) > 
        10.*std::abs( ( g2_Wefftab.back() - g2_Wefftab.front() ) / (pefftab.back() - pefftab.front()) );
//         return true;
      };
      
//       const real_t Wefftolerance=1.0e-3;
      
      while(( temp_size < Nmax 
            && advmath::getdiscrepancy(peff(), peffprev) > 1.0e-3
            && ( ( cond1() && cond2() ) || ( cond3() && cond4() ) ) 
            ) || std::abs((weffprev - weff())/(weffprev + weff()))>0.8)
      {
  #ifdef DEBUG
        std::cout << "Current value is " << peff() << " encoded in " << current->first << std::endl;
        std::cout << "  prev = " << prev->first << " next = " << next->first << std::endl;
  #endif
        const real_t newpeff = 0.5*( peff() + peffprev ) ;
        if(advmath::getdiscrepancy(peffprev, newpeff)<1.0e-12||newpeff< peffprev) break;
        const real_t sqrtStemp = 2.*std::sqrt( SQUARE(newpeff) + SQUARE(getMassLBSM()) );
        const real_t newweff = get_g2_Weff(sqrtStemp);
        
        list.insert(current, {newpeff, newweff});
        temp_size++;
        current = std::next(prev);
        next = std::next(current);
  #ifdef DEBUG
        std::cout << "Added point " << temp_size << " with peff = " << current->first << std::endl;
        std::cout << "  prev = " << prev->first << " next = " << next->first << std::endl;
  #endif
  #ifdef DEBUG
        std::cout << " New next value is " << next->first << " new current is " << current->first << std::endl;
        if(peffnext() -peff() <= 0.)
        {
          std::cerr << "Problem in the order\n"; 
          exit(1);
        }
  #endif
      }
      prev++; current++;next++; 
    }


  #ifdef DEBUG  
    std::cout << "Putting everything into the vectors again\n";
  #endif 
    sqrtStab.clear();
    pefftab.clear();
    g2_Wefftab.clear();
    sqrtStab.reserve(2 * Nmax);
    pefftab.reserve(2 * Nmax);
    g2_Wefftab.reserve(2 * Nmax);
    
    #ifdef DEBUG
      std::cout << "Initial list contents:\n";
      for (const auto& entry : list) 
      {
        std::cout << "entry.first = " << entry.first << ", entry.second = " << entry.second << std::endl;
      }
    #endif

    for (auto entry : list)
    {
      #ifdef DEBUG
        std::cout << "Processing entry: entry.first = " << entry.first << ", entry.second = " << entry.second << std::endl;
      #endif

      if ( std::isnormal(entry.second) || entry.second==0. )
      {
        pefftab.push_back(entry.first);
        g2_Wefftab.push_back(entry.second);

        real_t sqrtStabValue = 2. * std::sqrt(SQUARE(pefftab.back()) + SQUARE(getMassLBSM()));
        sqrtStab.push_back(sqrtStabValue);

        #ifdef DEBUG  
          std::cout << "Added: pefftab.back() = " << pefftab.back() 
                    << ", g2_Wefftab.back() = " << g2_Wefftab.back() 
                    << ", sqrtStab.back() = " << sqrtStab.back() << std::endl;
        #endif
      }
      else
      {
        #ifdef DEBUG  
          std::cout << "Skipped: entry.first = " << entry.first 
                    << ", entry.second = " << entry.second << std::endl;
        #endif
      }
    }
  }

  sqrtStab.shrink_to_fit();
  g2_Wefftab.shrink_to_fit();
  pefftab.shrink_to_fit();
    
//     sqrtStab.clear();
//     pefftab.clear();
//     g2_Wefftab.clear();
//     sqrtStab.reserve(2*Nmax);
//     pefftab.reserve(2*Nmax);
//     g2_Wefftab.reserve(2*Nmax);
//     for( auto entry : list )
//     {
//       if( !std::isnan(entry.second) && !std::isinf(entry.second))
//       { // Filtering nan and inf
//         pefftab.push_back( entry.first );
//         g2_Wefftab.push_back( entry.second );
//         sqrtStab.push_back( 2.*std::sqrt( SQUARE(pefftab.back()) + SQUARE(getMassLBSM() ))) ;
//       }
//     }
//   } // End of the if statement to refill the vectors with additional elements
//   
//   sqrtStab.shrink_to_fit();
//   g2_Wefftab.shrink_to_fit();
//   pefftab.shrink_to_fit();
}

void AvgSvCalculator::print_g2_wefftable(std::ostream& out) const
{
/*   Prints out the g2_Wefftab on the output stream specified
 *   Fields are:
 *   peff sqrtS Weff
 *   separated by a \t
 */
  
  for(size_t i = 0 ; i != g2_Wefftab.size() ; i++)
    out <<  pefftab[i] << '\t'  << sqrtStab[i] << '\t' << g2_Wefftab[i] << '\n';
  
}

// Applying formula 65 from Edsjio-Gondolo "Neutralino relic desity including coannihilation"
// with no optimisations
real_t AvgSvCalculator::getAverageSigmav_coan_noopti(const real_t &T)
{
  real_t result = 0.;
  
  // Allocating g2_Wefftab if necessary
  if(g2_Wefftab.size() == 0) 
  {
    constexpr const real_t Tfreeze=25.;
    constexpr const real_t Beps=1.e-6;
    const real_t maxenergy=2.*getMassLBSM()-Tfreeze*std::log(Beps);
    tabulateValues(std::max(10.*getMassLBSM(), maxenergy), 
                   3000);
  }
  
//   real_t last_value = advmath::K1(sqrtStab.back()/T)*std::pow(pefftab.back(), 2)*Wefftab.back();
//   std::cout << "The last value of the integrand is " << last_value << std::endl;
//   std::cin.get();
  
  real_t side1, side2;
  
  side1 = advmath::K1(sqrtStab[0]/T)*std::pow(pefftab[0], 2)*g2_Wefftab[0];
  side2 = side1;
  
  for(size_t i = 1 ; i <= pefftab.size()-2 ; i++)
  {
   // Integrating with Riemann rule on the average of the size 
   side1 = side2; 
   side2 = advmath::K1(sqrtStab[i]/T)*std::pow(pefftab[i], 2)*g2_Wefftab[i];
   real_t deltap = pefftab[i] - pefftab[i-1];
   
   if( deltap <= 0) 
   {
     std::cout << "Error! deltap <= 0\n";
   // std::cin.get();
     continue;
   }
   result+=0.5*(side1+side2)*deltap;
  }
  // The dependency on g^2_DM is eliminated in the ratio
  result=result/T;
  
  real_t denominator = 0.;
  
  for(auto part : corr::bsm_particles)
  {
    denominator+= corr::part_hel_dof[part]*SQUARE(input.masses_vector[part])*advmath::K2(input.masses_vector[part]/T);
  }
  
  result = result/denominator/denominator;

  return result;
}

real_t AvgSvCalculator::getAverageSigmav_coan_siso(const real_t &T )
{
// This function computes the average annihilation cross-section by using the same algorithm 
// as SuperIso Relic 4.1

  constexpr const real_t Tfreeze=25.;
  constexpr const real_t Beps=1.e-6;
  const real_t maxenergy=2.*getMassLBSM()-Tfreeze*std::log(Beps);
  
  if(g2_Wefftab.size() == 0) tabulateValues(std::max(10.*getMassLBSM(),
                                                  maxenergy)
                                        , 600);
  
  const real_t T1 = 1.0e-4*getMassLBSM();

#ifdef DEBUG_OUTPUT
#ifdef OUTPATH
  std::string filename=static_cast<std::string>(OUTPATH)+static_cast<std::string>("sigmav_tab_full_new.dat");
  FILE *ftabptr=fopen(filename.c_str(), "a");
  if ( ftabptr == nullptr )
  {
      printf( "File sigmav_tab_full_new.dat is impossible to open for appending\n");
    // std::cin.get();
  }
#endif
#endif

  if( T < T1 ) 
  {
#ifdef DEBUG_OUTPUT
    std::cout << " T < T1 : T-T1 = " << T-T1 << std::endl;
    if(ftabptr!=nullptr)
    {
      // Let's print T, T1, peff, Weff, 
      //             T  T1 <sigmav>
      fprintf(ftabptr, "-1\t%.5e\t%.5e\t%.5e\t%.5e\n", T, T1, pefftab.front(), g2_Wefftab.front() );
      fprintf(ftabptr, "-1\t%.5e\t%.5e\t%.5e\n", T, T1, g2_Wefftab.front()/( 4. * SQUARE(getMassLBSM()) * input.getLightestBSMdof()) );
      fclose(ftabptr);
    }
#endif
    return g2_Wefftab.front()/( 4. * SQUARE(getMassLBSM()) * input.getLightestBSMdof() );
  }
  
  real_t numerator = 0., denominator = 0.;
  const size_t n_max= pefftab.size()-2;
  real_t C0;
  real_t A, B, KA, KB, C; 
  
  
  KA = advmath::K1exp(sqrtStab[0]/T, 2.*getMassLBSM()/T);
  A  = SQUARE(pefftab[0])*g2_Wefftab[0]*KA;
  KB = advmath::K1exp(sqrtStab[0+1]/T, 2.*getMassLBSM()/T);
  B  = SQUARE(pefftab[0+1])*g2_Wefftab[0+1]*KB;
  C0 = 0.5*(A+B)*(pefftab[0+1]-pefftab[0]);
#ifdef DEBUG_OUTPUT
  if(ftabptr != nullptr) fprintf(ftabptr, 
    "-2\t%.5e\t%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", 
            T, 0, pefftab[0], sqrtStab[0], A, KA, B, (A+B),(pefftab[0+1]-pefftab[0]),  C0);
  if( C0 != 0. && (!std::isnormal(C0) || !(C0>0.) ) )
  {
      std::string contributiontype = ( C0 < 0. ) ? "negative" : "non-positive" ;
      std::cout << "Warning in the numerator of <sigmav>:\n";
      std::cout << contributiontype << " contribution to <sigma v> for T = " << T << " GeV:" << std::endl;
      std::cout << " A = " << A << " B = " << B << " C = " << C0 << std::endl;
      numerator+=0.;
  }else
#endif
  numerator+=C0;
  
  for(size_t i = 1; i < n_max; i++)
  {
    KA = KB;
    A  = B;
    KB = advmath::K1exp(sqrtStab[i+1]/T, 2.*getMassLBSM()/T);
    B  = SQUARE(pefftab[i+1])*g2_Wefftab[i+1]*KB;
    C  = 0.5*(A+B)*(pefftab[i+1]-pefftab[i]);
//     if(C/C0 < 1.0e-6) break;
#ifdef DEBUG_OUTPUT
    if(ftabptr != nullptr) fprintf(ftabptr, 
      "-2\t%.5e\t%ld\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", 
      T, i, pefftab[i], sqrtStab[i], A, KA, B, (A+B),(pefftab[i+1]-pefftab[i]),  C);
    if( (C != 0.) && (!std::isnormal(C) || !(C>0.)) )
    {
        std::string contributiontype = ( C < 0. ) ? "negative" : "non-positive" ;
        std::cout << "Warning in the numerator of <sigmav>:\n";
        std::cout << contributiontype << " contribution to <sigma v> for T = " << T << " GeV:" << std::endl;
        std::cout << " A = " << A << " B = " << B << " C = " << C << std::endl;
        numerator+=0.;
    }else
#endif
    numerator+=C;
  }
  
  for(auto part : corr::bsm_particles)
  {
    C =  corr::part_hel_dof[part]*SQUARE(input.masses_vector[part])*advmath::K2exp( input.masses_vector[part]/T , getMassLBSM()/T);
#ifdef DEBUG
    if( !std::isnormal(C) || !(C>0) ) 
    {
        std::cout << "Warning in the denominator of <sigmav>: ";
        std::cout << " C = " << C << std::endl;
        std::cout << " part = " << corr::part_names[part] ;
        std::cout << " dof  = " << corr::part_hel_dof[part] ;
        std::cout << " mass = " << input.masses_vector[part];
        std::cout << " arg1 = " << input.masses_vector[part]/T;
        std::cout << " arg2 = " << getMassLBSM()/T;
        std::cout << " K2exp(arg1,arg2) = " << advmath::K2exp( input.masses_vector[part]/T , getMassLBSM()/T) << std::endl;;
    }
#endif
    denominator+= C;
  }

#ifdef DEBUG_OUTPUT
  fclose(ftabptr);
  filename=static_cast<std::string>(OUTPATH)+static_cast<std::string>("sigmav_tab_new.dat");
  FILE *fptr=fopen(filename.c_str(),"a");
  
  if(fptr==nullptr) 
  {
      printf("Impossible to open sigmav_tab_new.dat\n");
    // std::cin.get();
  }else
  {
    fprintf(fptr, "%.5e\t%.5e\t%.5e\t%.5e\n", T, numerator, denominator, std::sqrt(2.*T/getMassLBSM()));
    fclose(fptr);
  }

#endif
  // The square root compensates the use of K1,2exp at the place of K1,2
  denominator*=denominator*T*std::sqrt(2.*T/getMassLBSM()); 
  
  real_t result = numerator/denominator;
  
  return !std::isnan(result) ? result : 0.;
}

// This function 

real_t AvgSvCalculator::getAverageSigmav_coan_hightemp(const real_t &T)
{

  constexpr const real_t Tfreeze=25.;
  constexpr const real_t Beps=1.e-6;
  const real_t maxenergy=2.*getMassLBSM()-Tfreeze*std::log(Beps);
  
  if(g2_Wefftab.size() == 0) tabulateValues(std::max(10.*getMassLBSM(),
                                                  maxenergy)
                                        , 600);

#ifdef DEBUG_OUTPUT
  std::string filename=static_cast<std::string>(OUTPATH)+static_cast<std::string>("sigmav_tab_full_new.dat");
  FILE *ftabptr=fopen(filename.c_str(), "a");
  if ( ftabptr == nullptr )
  {
      printf( "File sigmav_tab_full_new.dat is impossible to open for appending\n");
    // std::cin.get();
  }
#endif
  
  real_t numerator = 0., denominator = 0.;
  const size_t n_max= pefftab.size()-2;
  real_t C0;
  real_t A, B, KA, KB, C; 
  real_t num_x1, num_x2;
  const real_t num_z = 2.*getMassLBSM();
  const real_t denom = T;
  
  num_x1=sqrtStab[0];
  num_x2=sqrtStab[0+1];
  
  KA = advmath::K1exp_numdenom(num_x1, num_z, denom);
  A  = SQUARE(pefftab[0])*g2_Wefftab[0]*KA;
  KB = advmath::K1exp_numdenom(num_x2, num_z, denom);
  B  = SQUARE(pefftab[0+1])*g2_Wefftab[0+1]*KB;
  C0 = 0.5*(A+B)*(pefftab[0+1]-pefftab[0]);
#ifdef DEBUG_OUTPUT
  if(ftabptr != nullptr) fprintf(ftabptr, 
    "-2\t%.5e\t%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", 
            T, 0, pefftab[0], sqrtStab[0], A, KA, B, (A+B),(pefftab[0+1]-pefftab[0]),  C0);
  if( C0 != 0. && (!std::isnormal(C0) || !(C0>0.) ) )
  {
      std::string contributiontype = ( C0 < 0. ) ? "negative" : "non-positive" ;
      std::cout << "Warning in the numerator of <sigmav>:\n";
      std::cout << contributiontype << " contribution to <sigma v> for T = " << T << " GeV:" << std::endl;
      std::cout << " A = " << A << " B = " << B << " C = " << C0 << std::endl;
      numerator+=0.;
  }else
#endif
  if(!std::isnan(C0) && !std::isinf(C0)) numerator+=C0;
  else return -1.;
  
  for(size_t i = 1; i < n_max; i++)
  {
    KA = KB;
    A  = B;
    
    num_x2=sqrtStab[i+1];
    
    KB = advmath::K1exp_numdenom(num_x2, num_z, denom);
    
    B  = SQUARE(pefftab[i+1])*g2_Wefftab[i+1]*KB;
    C  = 0.5*(A+B)*(pefftab[i+1]-pefftab[i]);
    // if(C/C0 < 1.0e-6) break;
    #ifdef DEBUG_OUTPUT
      if(ftabptr != nullptr) fprintf(ftabptr, 
        "-2\t%.5e\t%ld\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", 
        T, i, pefftab[i], sqrtStab[i], A, KA, B, (A+B),(pefftab[i+1]-pefftab[i]),  C);
      if( (C != 0.) && (!std::isnormal(C) || !(C>0.)) )
      {
          std::string contributiontype = ( C < 0. ) ? "negative" : "non-positive" ;
          std::cout << "Warning in the numerator of <sigmav>:\n";
          std::cout << contributiontype << " contribution to <sigma v> for T = " << T << " GeV:" << std::endl;
          std::cout << " A = " << A << " B = " << B << " C = " << C << std::endl;
          numerator+=0.;
      }else
    #endif
    if(!std::isnan(C) && !std::isinf(C)) numerator+=C;
    else return -1.;
  }
  // Now computing the denominator
  // Note that we consider only g_i (i.e. the degrees of freedom of the particle i)
  // and not g_i/g_1, since the numerator contains already g_1 ^2* Weff, 
  // and later in the square they compensate
  
  for(auto part : corr::bsm_particles)
  {
      num_x2=input.masses_vector[part];
      
      C =  1.0*corr::part_hel_dof[part]*SQUARE(input.masses_vector[part])*advmath::K2exp_numdenom( num_x2 , num_z/2., T);
      #ifdef DEBUG
          if( !std::isnormal(C) || !(C>0) ) 
          {
              std::cout << "Warning in the denominator of <sigmav>: ";
              std::cout << " C = " << C << std::endl;
              std::cout << " part = " << corr::part_names[part] ;
              std::cout << " dof  = " << corr::part_hel_dof[part] ;
              std::cout << " mass = " << input.masses_vector[part];
              std::cout << " arg1 = " << input.masses_vector[part]/T;
              std::cout << " arg2 = " << getMassLBSM()/T;
              std::cout << " K2exp(arg1,arg2) = " << advmath::K2exp( input.masses_vector[part]/T , getMassLBSM()/T) << std::endl;;
          }
      #endif
      if( std::isnan(C) || std::isinf(C) ) return -1.;
      //NOTE : this is not the denominator, but just the sum, that will be squared and normalised with the prefactor later on
      denominator+= C;
  }
  

#ifdef DEBUG_OUTPUT
  fclose(ftabptr);
  filename=static_cast<std::string>(OUTPATH)+static_cast<std::string>("sigmav_tab_new.dat");
  FILE *fptr=fopen(filename.c_str(),"a");
  
  if(fptr==nullptr) 
  {
      printf("Impossible to open sigmav_tab_new.dat\n");
    // std::cin.get();
  }else
  {
    fprintf(fptr, "%.5e\t%.5e\t%.5e\t%.5e\n", T, numerator, denominator, std::sqrt(2.*T/getMassLBSM()));
    fclose(fptr);
  }
#endif

  // Now getting the result  
  // The square root compensates the use of K1,2exp at the place of K1,2
  denominator*=denominator*T*std::sqrt(2.*T/getMassLBSM()); 
  real_t result = numerator/denominator;
  return result;
}

struct SecondArgument_t{
  gsl_function *fptr;
  real_t error;
};

void AvgSvCalculator::setLimitingValues()
{
    /* This function calculates Weff(m_DM) and its derivatives
     *
     */  
  
#ifdef DEBUG
    std::cout << "Calculating limits\n" ;
#endif
    
#ifdef DEBUG
    std::cout << "Resetting vectors\n" ;
#endif
    Weffderiv.clear();
    Weffderiv.reserve(orderT+1);
    
    constexpr const real_t epsilon=1.0e-12/55;
    real_t limit;
    
    do
    {
      real_t prefac=1.0001;
      limit=getWeff(2.*getMassLBSM()*prefac);  
      prefac*=prefac;
    }while(std::isnan(limit) || std::isinf(limit));
    
    
    real_t oldlimit=10*limit;
    constexpr const int maxcount=100;
    constexpr const real_t maxdiscrepancy=1.0e-6;
    
    // Defining a 2-arguments lambda function that returns 
    // Weff(y) instead of Weff(Ecm)
    // It needs a void* as 2nd argument, to be compatible with the gsl_function format.
    auto lambdweffofy = [](real_t y, void *myClassinstance)
    {
      AvgSvCalculator *c=(AvgSvCalculator*)(myClassinstance);
      return c->getWeff(2.*(c->getMassLBSM())*std::sqrt(y));//OK: Weff is already good for the g_DM^2 factor
    };
    
    // Calculating the limit of Weff for y->1
    for(int count=1; count <= maxcount && std::abs(limit-oldlimit)/std::abs(limit+oldlimit) >= maxdiscrepancy; count++)
    {
#ifdef DEBUG
      std::cout << " count = " << count << std::endl;
#endif
      oldlimit=limit;
      limit=lambdweffofy(1.+epsilon/count,this);
    }
#ifdef DEBUG
    std::cout << " limit = " << limit << " +/- " << limit*std::abs(limit-oldlimit)/std::abs(limit+oldlimit)/2.   
    << std::endl;
#endif
    
    // Storing the value of the limit and its error in the vector
    Weffderiv.push_back(limit);
    Weffderiv_err.push_back(limit*std::abs(limit-oldlimit)/std::abs(limit+oldlimit)/2.);
    
    // Defining a lambda functions to calculate the right derivative of its second
    // argument in y. In fact, we need to do that at the threshold of Weff.
    // It needs a void* as 2nd argument, to be compatible with the gsl_function format.
    auto lambda1weffofy = [](real_t y, void * arg2 )
    {
#ifdef DEBUG
      std::cout << " Evaluating right derivative in lambda\n";
#endif
      real_t result, deltax=1.0e-6, abserr;
      
      // We want to use the gsl library, therefore we need to define the proper gsl_function
      // associated to the second argument      
      gsl_function *f = ((SecondArgument_t *)(arg2))->fptr;
      
      // Computing the derivative
      [[maybe_unused]] int exit_code=gsl_deriv_forward(f, y, deltax, &result, &abserr);
      
      // Storing the error as a member of arg2
      ((SecondArgument_t *)(arg2))->error=abserr;
      
#ifdef DEBUG
      std::cout << " The value is " << result << " +/- " << abserr << std::endl;
#endif
      
//       if(result != 0. && abserr/result > 0.5)
//       {
//         std::cout << " Returning 0.\n";
//         return 0.;
//       }
      return result;
    };
    
    
    // To compute the derivatives, we use the previous lambda function recursively
    // We therefore create 3 arrays of size 30 (very unlikely one wants to go beyond the 30th derivative)
    // Each array contains:
    std::array<gsl_function,30> fgsl;       // The gsl_function for the function in input
    std::array<SecondArgument_t,30> second; // The second argument to be provided for the lambda function 
    std::array<real_t,30> results;          // The results of each calculation
    
#ifdef DEBUG
    std::cout << " Order of the derivative = 1\n";
#endif
    
    // Let us define the gsl_function fo the 1st derivative: 
    //  - the function is the lambda function we defined first 
    //  - the argument is the current instance of the class
    fgsl[0].function=lambdweffofy;
    fgsl[0].params=(void*)(this);
    
    // Let us define the SecondArgument_t needed to compute the 1st derivative:
    //  - the gsl_function is the one we just defined
    //  - we initiate the error at 0.
    second[0].fptr=&fgsl[0];
    second[0].error=0.;
    
    // We compute the 1st result, and put it in the vector
    results[0] = lambda1weffofy(1.,(void *)(&second[0]));
    
    Weffderiv.push_back(results[0]);
    Weffderiv_err.push_back(second[0].error);
    
    for(int count = 1; count < orderT; count ++)
    {
#ifdef DEBUG
    std::cout << " Order of the derivative = "<< count+1 << std::endl;
#endif
      // For the derivatives of order >=2, we need to define the gsl_function 
      // as the lambda function that computes the derivative, and its parameters
      // as the parameters used in the previous iteration (the one in which we computed 
      // such a derivative)
      fgsl[count].function=lambda1weffofy;
      fgsl[count].params=(void*)(&second[count-1]);
      
      
      // This is the same as the 1st derivative
      second[count].fptr= &fgsl[count];
      second[count].error=0.;
      results[count]=lambda1weffofy(1., (void*)(&second[count]));
      
      // Saving the results in the vectors
      Weffderiv.push_back(results[count]);
      Weffderiv_err.push_back(second[count].error);
      
      if(std::abs(second[count].error/results[count]) > 0.5)
      {
        orderT=count+1;
        std::cout << "Error too large for order " << count +1 <<": stopping at the order "<< orderT << " for the derivatives of Weff\n";
        break;
      }
    }

}

void AvgSvCalculator::setTaylorCoeffSVx_nosplitting()
{
#ifdef DEBUG
  std::cout << "Entering setTaylorCoeffSVx_nosplitting()\n";
#endif
  
  if(Weffderiv.size()==0) setLimitingValues();
  
  const real_t prefactor=std::pow(4. * SQUARE(getMassLBSM()),-1.);
  real_t tempresult=0.;
  
  TaylorCoeffSVT.clear();
    for(int i=0; i<= orderT; i++)
    {
      
      tempresult=0.;
      for(int j=0; j<= i ; j++)
      {
        tempresult+=CoefflowTsv[i][j]*Weffderiv[j];
        
#ifdef DEBUG      
        std::cout << "i = " << i << " j = " << j << " Cij = " << CoefflowTsv[i][j] 
        << " Wj = " << Weffderiv[j] << " res = " << tempresult << std::endl;
        assert(!std::isnan(Weffderiv[j]));
#endif
        
      }
      
      TaylorCoeffSVT.push_back(prefactor*tempresult);
  }
  
#ifdef DEBUG   
  std::cout << "prefactor = " << prefactor << std::endl;
  std::cout << "m1 = " << getMassLBSM() << " g1 = " << input.getLightestBSMdof() << std::endl;
  print_coefficients();
#endif
  
  
}

#ifdef EXPERIMENTAL__SPLITTING
void AvgSvCalculator::setTaylorCoeffSVx_splitting(const real_t deltax[], const int M_rho, const int M_eta)
{
// EXPERIMENTAL FEATURE: WE ARE TESTING THIS!
  
#ifdef DEBUG
  std::cout << "Entering setTaylorCoeffSVx_splitting()\n";
#endif
//   if(werealldeltaxinull) return;
  TaylorCoeffSVT.clear();
  
  /* If there is mass splitting,
     the Taylor coefficients depend on eta, which depends from 1/T. 
     So they need to be updated */
  
  real_t phitilde[NlowTsv_split+1]={0.};
  std::vector<real_t> ci; 
  real_t eta_o1=0.; // Parameter eta at the first order in Delta x_i
  
  /* We recall that eta_o1 = eta_tilde/rho_tilde
  
     I need to define the phitildes
     Firstly, I need the eta_o1
     it is defined as eta_tilde/rho tilde 
     I define the cis */
  
  ci.reserve(corr::SIZEPHYSICALBSM);
  ci.emplace_back(1.);
  for(size_t i = 1 ; i < ci.capacity() ; i++)ci.emplace_back(
    1.*(corr::part_hel_dof[bsm_particles_sorted[i]]/input.getLightestBSMdof()
      *SQUARE(input.masses_vector[bsm_particles_sorted[i]]/getMassLBSM()))
  );
  //  I define rho_tilde
  real_t rho_tilde=0.;
  for(int i = 0 ; i < M_rho ; i++) for(int j = i ; j < M_rho ; j++)rho_tilde+=(2-kdelta(i,j))*ci[i]*ci[j];
  
  
  // I define eta_tilde
  real_t eta_tilde=0.;
  for(int i = 0 ; i < M_eta ; i++) for(int j = i ; j < M_eta ; j++)eta_tilde-=(1.-0.5*kdelta(i,j))*ci[i]*ci[j]*(deltax[i]+deltax[j]);
  
  werealldeltaxinull=(eta_tilde==0.);
  
  eta_o1=eta_tilde/rho_tilde;
  
#ifdef DEBUG
  std::cout << "Preferred values found are rho_tilde = " << rho_tilde << " eta_o1 = " << eta_o1 << std::endl;
  
  std::cout << "Tabulating values...\n";
  
  for(size_t i = 0 ; i < ci.size() ; i++)
  {
    printf("c[%3ld] = %.3e\n", i, ci[i]);    
  }
  
  real_t rho_tilde_tab[corr::SIZEPHYSICALBSM][corr::SIZEPHYSICALBSM], 
         eta_o1_tab[corr::SIZEPHYSICALBSM][corr::SIZEPHYSICALBSM], 
         zerothcoeff[corr::SIZEPHYSICALBSM][corr::SIZEPHYSICALBSM];
         
  for(int M_rho_temp = 0 ; M_rho_temp < corr::SIZEPHYSICALBSM ; M_rho_temp++)
  {
    rho_tilde_tab[M_rho_temp][0] = 0.;
    
    for(int i = 0 ; i < M_rho_temp ; i++) 
      for(int j = i ; j < M_rho_temp ; j++)
        rho_tilde_tab[M_rho_temp][0]+=(2-kdelta(i,j))*ci[i]*ci[j];
    
    for(int M_eta_temp = 0; M_eta_temp < corr::SIZEPHYSICALBSM ; M_eta_temp++)
    {
      rho_tilde_tab[M_rho_temp][M_eta_temp] = rho_tilde_tab[M_rho_temp][0];
      
      eta_o1_tab[M_rho_temp][M_eta_temp]=0.;
      for(int i = 0 ; i < M_eta_temp ; i++) 
        for(int j = i ; j < M_eta_temp ; j++)
          eta_o1_tab[M_rho_temp][M_eta_temp]
           -=(1.-0.5*kdelta(i,j))*ci[i]*ci[j]*(deltax[i]+deltax[j]);
      eta_o1_tab[M_rho_temp][M_eta_temp]/=rho_tilde_tab[M_rho_temp][M_eta_temp];
      
      zerothcoeff[M_rho_temp][M_eta_temp]=std::pow(4.*rho_tilde_tab[M_rho_temp][M_eta_temp]*SQUARE(getMassLBSM()), -1)*Weffderiv[0]*(sv_beta[0]+sv_lambda[0]*eta_o1_tab[M_rho_temp][M_eta_temp]);
    } 
  }
  
  for(size_t i = 0 ; i < ci.size() ; i++)
  {
    printf("rho_tild[%3ld] = %.3e\n", i, rho_tilde_tab[i][0]);    
  }
  
  
  std::cout << "Tabulated values of the two for different combinations are \n";
  printf("%6s", "Mr\\Me");
  // First row
  for(int i = 0 ; i <corr::SIZEPHYSICALBSM ; i++)
    printf("%19d", i);
  printf("\n");
  
  // Next rows 
  for(int i = 0 ; i <corr::SIZEPHYSICALBSM ; i++)
  {
    printf("%6d", i);
    for(int j = 0 ; j < corr::SIZEPHYSICALBSM ; j++)
    {
      char tempstring[20];
      sprintf(tempstring, "(%.1e %.1e)", rho_tilde_tab[i][j], eta_o1_tab[i][j]);
      printf("%19s", tempstring);
    }
    printf("\n");
    
    printf("%6s", " ");
    for(int j = 0 ; j < corr::SIZEPHYSICALBSM ; j++)
    {
      char tempstring[20];
      sprintf(tempstring, "(%.1e)",zerothcoeff[i][j]);
      printf("%19s", tempstring);
    }
    printf("\n");
  }
#endif
  
  eta_o1=(eta_o1 > 0.5) ? 0. : eta_o1;
  
  // it is a linear relation in the variable beta (at the first order)
  for(int j = 0 ; j <= orderT; j++) phitilde[j]=sv_beta[j]+eta_o1*sv_lambda[j];
  
  // I define the prefactor
  const real_t prefactor=std::pow(4.*rho_tilde*SQUARE(getMassLBSM()), -1);
  
  for(int n = 0; n <=orderT ; n++)
  {
    real_t temp_coeff=0.;
#ifdef DEBUG 
    std::cout << "n = " << n << std::endl;
#endif
    for(int i = 0 ; i <=n ; i++)
      for(int j = 0 ; j <=n ; j++)
      {
        temp_coeff+=Weffderiv[i]*Qnij[n][i][j]*phitilde[j];
#ifdef DEBUG
        std::cout << "Weffderiv[" << i << "] = " << Weffderiv[i] << 
                   "\nphitilde["  << j << "] = " << phitilde[j] <<
                   "\nQ_[" << n << " " << i << " " << j << "] = " << Qnij[n][i][j] << 
                   "\ntemp_coeff = " << temp_coeff << 
                   "\nprefactor = " << prefactor <<
                   std::endl;
#endif
      }
    
    // Pushing back the whole coefficient in 1/x
    TaylorCoeffSVT.emplace_back(prefactor*temp_coeff);
  }
#ifdef DEBUG 

  std::cout << "rho_tilde = " << rho_tilde 
           << "\neta_tilde = " << eta_tilde
           << "\neta_o1 = " << eta_o1 
           << "\nm1 = " << getMassLBSM() << " g1 = " << input.getLightestBSMdof() << std::endl;
  print_coefficients();
// std::cin.get();
#endif
  wastherenosplitting=false;
}
#endif

void AvgSvCalculator::print_coefficients(std::ostream& out) const
{
  size_t count = 0; 
  
  for(count = 0 ; count < Weffderiv.size() ; count ++ )
  {
    out << "Weffderiv["<< count << "] = " << Weffderiv[count] << " +/- " << Weffderiv_err[count] <<    
    '\n';
  }
  count =0;
  
  for(auto elem : TaylorCoeffSVT)
  {
    out << "TaylorCoeffSVT["<< count << "] = " << elem << " GeV^-" << (2+count) <<  '\n';
    count++;
  }
}

real_t AvgSvCalculator::getAverageSigmav_coan_lowtemp(const real_t &T)
{
#ifdef DEBUG
  std::cout << "Entering getAverageSigmav_coan_lowtemp\n";
#endif
  real_t result=0. ;
  if(Weffderiv.size()==0)
  {
    setLimitingValues();
  }
  const real_t x = getMassLBSM()/T;
  setTaylorCoeffSVx_nosplitting();  
  
  for(int i=0; i<= orderT; i++)
  {
#ifdef DEBUG 
    std::cout  << " i = " << i << " TaylorCoeff = " << TaylorCoeffSVT[i] << std::endl;
    assert(std::isnormal(TaylorCoeffSVT[i]));
#endif
    result+=TaylorCoeffSVT[i]*std::pow(x, -i);
  }

  return result>0 ? result : 0.;
}

#ifdef EXPERIMENTAL__SPLITTING
real_t AvgSvCalculator::getAverageSigmav_coan_lowtemp_degenerate_mass(const real_t &T)
{
#ifdef DEBUG
  std::cout << "Entering getAverageSigmav_coan_lowtemp\n";
#endif
  real_t result=0. ;
  if(Weffderiv.size()==0)
  {
    setLimitingValues();
  }
  const real_t x = getMassLBSM()/T;
/*
  //Here I evaluate the mass splitting condition
  //  Firstly I have to assign the values of the deltax_i , and define some control variables 
  std::array<real_t, corr::SIZEPHYSICALBSM> deltax;
  bool aretherenans=false;
  int M_eta=0, M_rho=0;
  
  for(int i=0; i < corr::SIZEPHYSICALBSM; i++) 
  {
    deltax[i]=(getMass(bsm_particles_sorted[i])-getMassLBSM())/T;
    if(std::isinf(deltax[i])||std::isnan(deltax[i])) aretherenans=true;
  }
  
#ifdef DEBUG
  for(int i=0; i < corr::SIZEPHYSICALBSM; i++) 
  {
    std::cout << "deltax[" << i<< "]= " << deltax[i] << std::endl;
  }
#endif
  
  if(!aretherenans)
  {
#ifdef DEBUG
    std::cout << "All Delta x_i are well defined\n";
    std::cout << "They are " << corr::SIZEPHYSICALBSM << std::endl;
#endif
    // If all the elements of deltax_i are well-defined we can count how many particles 
    // contribute to betatilde. We do it in in the variable M
    const real_t thresholdx=1.0e-1, thresholdm=1.0e-2;
    
    // Counting M_eta
    for(M_eta=0 ;  M_eta<corr::SIZEPHYSICALBSM && deltax[M_eta]<thresholdx ; M_eta++) 
    {
#ifdef DEBUG
      std::cout << "They are " << corr::SIZEPHYSICALBSM << std::endl;
      if(deltax[M_eta]<thresholdx) std::cout << "Below the thereshoold\n";
      if(M_eta<corr::SIZEPHYSICALBSM) std::cout << "M_eta = " << M_eta << std::endl;
#endif
    }
    
    // Counting M_rho
    for(M_rho=0 ;  M_rho<corr::SIZEPHYSICALBSM && deltax[M_rho]*T<thresholdm ; M_rho++) 
    {
#ifdef DEBUG
      std::cout << "They are " << corr::SIZEPHYSICALBSM << std::endl;
      if(deltax[M_rho]<thresholdm) std::cout << "Below the thereshoold\n";
      if(M_rho<corr::SIZEPHYSICALBSM) std::cout << "M_rho = " << M_rho << std::endl;
#endif
    }
#ifdef DEBUG
    std::cout << "M_rho = "<< M_rho << std::endl;
#endif
    
    // Calling the splitting function for M_rho M_eta found
    if(M_rho>1||M_eta>1) setTaylorCoeffSVx_splitting(&deltax[0], M_rho, M_eta);
  }
  
  if((M_rho <=1 && M_eta<1)||aretherenans || TaylorCoeffSVT.size()==0)
*/ 
  setTaylorCoeffSVx_nosplitting();
  
  
  for(int i=0; i<= orderT; i++)
  {
#ifdef DEBUG 
    std::cout  << " i = " << i << " TaylorCoeff = " << TaylorCoeffSVT[i] << std::endl;
    assert(std::isnormal(TaylorCoeffSVT[i]));
#endif
    result+=TaylorCoeffSVT[i]*std::pow(x, -i);
  }

  return result>0 ? result : 0.;
}
#endif

real_t AvgSvCalculator::getAverageSigmav_coan(const real_t &T)
{
  if( T_lim_sigmav > 0.&& T < T_lim_sigmav ) return getAverageSigmav_coan_lowtemp(T);
  
  real_t result=getAverageSigmav_coan_hightemp(T) ;
  
  return (std::isnormal(result) && result >0.) ? result : getAverageSigmav_coan_lowtemp(T);
}


} // Enf of namespace __SPEC_LIB_NAME__
