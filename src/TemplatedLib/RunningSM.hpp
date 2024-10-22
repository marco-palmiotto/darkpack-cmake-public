#pragma once

#ifndef __cplusplus
#define __cplusplus
#endif

#include "config.hpp"
#include "params_new.hpp"
#include "correspondance.hpp"
// #include <ctime>
// #include <unistd.h>
// #include <iostream>
// #include <iomanip>
// #include <cstdlib>
// #include <fstream>
// #include <cassert>
// #include <thread>
// #include <future>


namespace __SPEC_LIB_NAME__
{

struct RunningSM{
  public:
    enum ParticlesList : short int { UP=0, DOWN, STRANGE, CHARM, BEAUTY, TOP, EL, MU, TAU, NUE, NUMU, NUTAU, GLUON, W, Z, PHOTON, HIGGS};
  
  private:
    
    std::vector<std::string> ParticlesNames = { "UP", "DOWN", "STRANGE", "CHARM", "BEAUTY", "TOP", "EL", "MU", "TAU", "NUE", "NUMU", "NUTAU", "GLUON", "W", "Z", "PHOTON", "HIGGS"};
    
    // Parameters needed for the running
    // In the comment after the ; there's the code in the SMINPUT block
    csl::InitSanitizer<real_t> alphas_MZ {"alphas_MZ"};// from 3
    csl::InitSanitizer<real_t> mass_Z {"mass_Z"};//4
    
    csl::InitSanitizer<real_t> m_d2GeV {"m_d2GeV"};//21
    csl::InitSanitizer<real_t> m_u2GeV {"m_u2GeV"};//22
    csl::InitSanitizer<real_t> m_s2GeV {"m_s2GeV"};//23
    
    csl::InitSanitizer<real_t> mc_mc {"mc_mc"};//24
    
    csl::InitSanitizer<real_t> mb_mb {"mb_mb"};//5
    
    csl::InitSanitizer<real_t> mass_top_pole {"mass_top_pole"};//6
    
    // Boolean parameter to determine if the light quark masses will or not be ran
    bool runlightquarks, runcharm, higgsloops;
    
    // It is useful to store the pole masses of the b and the c
    csl::InitSanitizer<real_t> mass_b_pole {"mass_b_pole"};
    csl::InitSanitizer<real_t> mass_c_pole {"mass_c_pole"};
    
    
    // With those two arrays, we store each quark mass as provided 
    // in the inputs at its scale
    // For the top, we'll have to call GetMtopMtop()
    csl::InitSanitizer<real_t> QuarkMass[6], QuarkMassQ[6];

    // Parameters for optimization: here we store each quark mass
    // as computed in the last call for the running, and each scale.
    
    csl::InitSanitizer<real_t> LastQuarkMass[6], LastQuarkMassQ[6];
      
    // Parametrs for optimization: here we store the last value
    // of alpha_strong and its energy, as computed in the running
    csl::InitSanitizer<real_t> alphas {"alphas"};//
    csl::InitSanitizer<real_t> Q_alphas {"Q_alphas"};//
    
    csl::InitSanitizer<real_t> Lambda3 {"Lambda3"};
    csl::InitSanitizer<real_t> alphasMZ_Lambda3 {"alphasMZ_Lambda3"};
    
    csl::InitSanitizer<real_t> Lambda4 {"Lambda4"};
    csl::InitSanitizer<real_t> alphasMZ_Lambda4 {"alphasMZ_Lambda4"};
    
    csl::InitSanitizer<real_t> Lambda5 {"Lambda5"};
    csl::InitSanitizer<real_t> alphasMZ_Lambda5 {"alphasMZ_Lambda5"};
    
    csl::InitSanitizer<real_t> Lambda6 {"Lambda6"};
    csl::InitSanitizer<real_t> alphasMZ_Lambda6 {"alphasMZ_Lambda6"};
    
  public:
    real_t AlphaStrong(real_t Q, real_t mtop_pole=-1.0); 
    // Normally you don't need to specify the top pole mass, it's read
    // from the class element    
    real_t GetMtopMtop(); 
    
    void setSanitNames()
    {
      std::array<std::string, 4> listofprefixes = {"Mass_", "Q_", "LastMass_", "LastQ_"};
      std::array<csl::InitSanitizer<real_t> *, 4>  listofarrays = {QuarkMass,  QuarkMassQ,  LastQuarkMass,  LastQuarkMassQ};
      for(size_t i = 0 ; i <4 ; i++)
      {
        for(size_t j = 0 ; j <6 ; j++) 
          listofarrays[i][j].setName(listofprefixes[i]+ParticlesNames[j]); 
      }
      
    };
    
    void init(const Param_t &input);
    RunningSM(const Param_t &input){setSanitNames();init(input);};
    RunningSM();
    ~RunningSM()
    {
#ifdef DEBUG
      std::cout << "Called ~RunningSM\n";    
#endif
    };
    
    bool operator== (const RunningSM & ) const;
    
    bool operator!= (const RunningSM &other) const
    {
      return !(*this==other);
    };
    
    void print(std::ostream & out=std::cout)const;
    
    inline void RunLightQuarks(bool light=true)
    {
      runlightquarks=light;
      runcharm= runlightquarks ? true : runcharm;
    };
    
    inline void RunCharmMass(bool runvalue=true)
    {
     runcharm=runvalue;
     runlightquarks= runcharm ? runlightquarks : false;
    };
    
    inline void setHiggsLoopFactor(bool x=true)
    {
      higgsloops=x;
    };
      
    real_t GetTopPoleMass()const{return mass_top_pole;};

    real_t GetMbMb()const{return mb_mb;};
      
    real_t RunQuarkMass(real_t quark_mass, const real_t & Qinit, const real_t & Qfin);
    real_t GetMcPole();
    real_t GetMbPole();
    real_t GetMcPole1Loop();
    real_t GetMbPole1Loop();
    real_t GetMcPole3Loops();
    real_t GetMb1S();
    real_t GetMcMcFromPole(const real_t & mcpole, const int loop);
    real_t GetMbMbFromPole(const real_t & mbpole, const int loop);
    
    real_t GetQuarkMass(enum ParticlesList part, const real_t &Qf);
    void HandleParamRunning(Param_t &input, const real_t &Q);
};

void TestRunningMasses( RunningSM * const run, Param_t &input);

} // End of namespace __SPEC_LIB_NAME__

inline std::ostream& operator << (std::ostream& out, const __SPEC_LIB_NAME__::RunningSM &var)
{
  var.print(out);
  return out;
}
