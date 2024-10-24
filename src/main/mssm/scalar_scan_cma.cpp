#include <iostream>
#include <iomanip>
#include "dp_scalar2to2/config.hpp"
#include "dp_scalar2to2/correspondance.hpp" 
#include "dp_scalar2to2/leshouchesfrommarty.hpp"
#include "dp_scalar2to2/process.hpp"
#include "dp_scalar2to2/avgsvcalculator.hpp"
#include "dp_scalar2to2/boltzmann.hpp"
#include "libcmaes/cmaes.h"
#include <libcmaes/esoptimizer.h>
#include <libcmaes/cmastrategy.h>
#include <libcmaes/llogging.h>
#include <fstream>
using namespace scalar2to2;
using namespace advmath;
using namespace libcmaes;

const std::string PATHPLOTS=std::string(OUTPATH);
constexpr const double omega_h2_target=0.12;
constexpr const double omega_h2_err=0.00091;

std::string filename="";

static inline double sq(const double x)
{
  return x*x;
}

static inline double gaussian_shape(const double val, const double mean, const double stddev)
{
  return -std::exp( - sq(val-mean)/2./sq(stddev) );
}

double function_to_minimize(const double *parameters, const int npar=3)
{
  /* 
    BLOCK SCALARMODEL
      1  0   # g_d
      2  7.80000000E-1   # g_u    [0 -1]
      3  0   # g_l
      4  7.8000000E-1     # g_chi
      5  2.00000000E+03   # m_chi [100 GeV - 1 TeV]
      6  4.40000000E+03   # m_phi
  */
  // static unsigned int ncalls=0;
  static Param_t input(filename.c_str());
  input.g_d=0.;
  input.g_u=1.;
  input.refresh();
  static BoltzmannSolver boltz(input);

  // parameters in input will be
  input.g_chi = parameters[0];
  input.m_chi = parameters[1];
  input.m_phi = parameters[2];

  input.refresh();
  boltz.changeInput(input);
  const double omega=boltz.relic_density();
  return std::abs(omega - omega_h2_target);
}


int main(int argc, char *argv[])
{

  std::cout << "Running " << argv[0] << std::endl;

  if(argc < 2)
  {
    std::cout << "This function requires 1 argument\n";
    return 1;
  }
  filename=argv[1];
  
  std::cout << "Reading inputs from " << filename << std::endl;

  constexpr const double sigma=omega_h2_err;
  constexpr const int npar=3;
  // parameters in input will be
  // input.g_chi = parameters[0];
  // input.m_chi = parameters[1];
  // input.m_phi = parameters[2];
  std::vector<double> x0={0.4, 4.0e+3, 8.1e+3}; // Initial values
  double lbounds[npar]={1.e-4, 100. , 100. }, 
         ubounds[npar]={1.0  , 1.e+4, 1.e+4}; // arrays for lower and upper parameter bounds, respectively                

  constexpr const unsigned int nsols=4;
  double sols[nsols][npar]={
        {0.379524 , 999.867,  2200.91},
        {0.849991 , 320.633 , 153.001},
        {0.3772   , 1000.02 , 2199.99},
        {0.327959 , 4000.05 , 8100.03}
  }; 
  for(unsigned i=0 ; i < nsols ; i++)
  {
    std::cout << function_to_minimize(sols[i], npar) << std::endl;
  }

  GenoPheno<pwqBoundStrategy> gp(lbounds,ubounds,npar); // genotype / phenotype transform associated to bounds.  

  // -1 for automatically decided lambda, 0 is for random	seeding	of the internal generator. 
  CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0,sigma,-1,0,gp); 
                                               
  cmaparams.set_algo(aCMAES);
  FitFunc f = &function_to_minimize;
  CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(f, cmaparams);

  std::cout << "best solution: ";
  cmasols.print(std::cout,0,gp);
  std::cout << std::endl;
  std::cout << "Expected Distance from Minimum: " << cmasols.edm() << std::endl;
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";
  return cmasols.run_status();
}