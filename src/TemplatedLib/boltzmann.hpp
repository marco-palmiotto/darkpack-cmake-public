#pragma once

#ifndef __cplusplus
#define __cplusplus
#endif

#include "avgsvcalculator.hpp"
#include "relicparam.hpp"

namespace __SPEC_LIB_NAME__
{

class BoltzmannSolver : public AvgSvCalculator, public Relicparam_t
{  
  private :
  
    real_t T_prev_boltzphi;
    real_t rhophi_prev_boltzphi;
    real_t T_prev_boltzphi_log;
    real_t rhophi_prev_boltzphi_log;
  
  public :    
    
    void reset_prevs();
    
    void setQCDeosModel(const int x)
    {
      Relicparam_t::init(getMassLBSM(), x);
      reset_prevs();
    };
    
    BoltzmannSolver(const Param_t &param, int x=2):AvgSvCalculator(param),Relicparam_t(getMassLBSM(),x)
    { 
      reset_prevs();
    };
    
    BoltzmannSolver(const AvgSvCalculator &sop, int x=2):AvgSvCalculator(sop),Relicparam_t(getMassLBSM(),x)
    { 
      reset_prevs();
    };
    
    BoltzmannSolver(const Param_t &param, std::shared_ptr<std::vector<Process2to2>> proclist, int x=2)
      : AvgSvCalculator(param, proclist), Relicparam_t(getMassLBSM(),x)
    {
      reset_prevs();
    }
  
    unsigned int changeInput(const Param_t & other)
    {
      AvgSvCalculator::changeInput(other);
      Relicparam_t::init(other.getLightestBSMmass(), model_eff);
      reset_prevs();
      return 0;
    };    
    
//     ~BoltzmannSolver() {};
    
    void print(std::ostream& out=std::cout) const override
    {
      print_procs(out); 
      print_g2_wefftable(out);
      print_relicparam(out);
    };
    
    real_t Yeq(const real_t &T);
    
    real_t Sigmatildestar_init(const real_t & T, const real_t & delta);
    
    real_t dYeq_dT(const real_t &T); /* derivative of Yeq */
    real_t boltzright(const real_t & T, const real_t & sv, const real_t & Y, const real_t & Yphi);
    void   setTfo(real_t delta);
    real_t boltzright_phi(const real_t & T, const real_t & Y, const real_t & Yphi);
    real_t boltzright_log(const real_t & T, const real_t & sv, const real_t & Y, const real_t & Yphi);
    real_t boltzright_phi_log(const real_t & T, const real_t & Y, const real_t & Yphi);
    real_t relic_density_linsolver();
    real_t relic_density_logsolver();
    real_t relic_density();
};

} // End of namespace __SPEC_LIB_NAME__

inline std::ostream& operator << (std::ostream& out , const __SPEC_LIB_NAME__::BoltzmannSolver &var)
{
  var.print(out);
  return out;
}
