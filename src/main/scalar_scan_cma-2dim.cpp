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

int main(int argc, char *argv[])
{

  std::cout << "Running " << argv[0] << '\n';

  constexpr const double sigma=omega_h2_err;
  constexpr const int npar=2;
  constexpr const int n_required_args=npar+4;
  constexpr const char *par_names[npar] = {"g_chi", "m_phi/m_chi"};

  constexpr const double init_seed_g_chi=0.3;

  double lbounds[npar]={1.e-12, 0.25}, 
         ubounds[npar]={1.0  , 5.0 }; // arrays for lower and upper parameter bounds, respectively  

    std::cout << "This function requires " << n_required_args << " arguments:\n"
                << " - The name of the input file\n"
                << " - The name of the output file\n"
                << " - The seed for the parameter g_chi\n"
                << " - The seed for the parameter m_phi/m_chi\n"
                << " - The initial value for m_chi\n";

  if(argc < n_required_args)
  {
    return 1;
  }
  filename=argv[1];
  
  std::cout << "The input file is " << filename << '\n';    
 
  char *end;
  const double init_m_chi = std::strtod(argv[argc-1],&end);
  constexpr const double end_m_chi=1.e+5;
  constexpr const int npoints=1000;
  const double stepsize=std::pow(end_m_chi/init_m_chi, 1./npoints);     

  Param_t input(filename.c_str());
  input.g_d=0.;
  input.g_u=1.;
  input.m_chi=init_m_chi;
  input.refresh();
  BoltzmannSolver boltz(input);

  auto function_to_minimize = [&] (const double *parameters, const int npar=2)
  {
    // parameters in input will be
    input.g_chi = parameters[0];
    input.m_phi = input.m_chi*parameters[1];
    input.refresh();
    boltz.changeInput(input);
    const double omega=boltz.relic_density();
    return std::abs(omega - omega_h2_target);
  };     

  std::vector<double> x0;
  x0.reserve(npar);
  for(int i = 0 ; i < npar; i++)
  {
    const int index_arg=i+3;
    const double value=strtod(argv[index_arg],&end);
    if(value < lbounds[i] || value > ubounds[i])
    {
      std::cerr << "Error in argument: " << index_arg << ", parameter " << i <<  ", " << par_names[i] << '\n'
                << "value " << value << " is out of the bound interval"
                << " [ " << lbounds[i] << ", " << ubounds[i] << " ].\n";
      return 2;
    }
    std::cout << "x0.size()=" << x0.size() << ", emplacing back " << value << '\n';
    x0.emplace_back(value);
  }

  if(x0.size() != npar)
  {
    std::cerr << "Initial parameters have wrong size: x0.size()=" << x0.size() << ", npar=" << npar << '\n';
    return 2;
  } 

  GenoPheno<pwqBoundStrategy> gp(lbounds,ubounds,npar); // genotype / phenotype transform associated to bounds.  
  std::cout << "FIXED: " << input.m_chi
            << "SEED: g_chi=" << x0[0] << " , m_phi/m_chi=" << x0[1] << '\n';
  CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0,sigma,-1,0,gp); // -1 for automatically \
decided lambda, 0 is for random	seeding	of the internal generator.                                              
  cmaparams.set_algo(aCMAES);
  FitFunc f = function_to_minimize;
  CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(f, cmaparams);

  std::cout << "seed is\n";
  for(int i=0; i < npar ; i++) std::cout << par_names[i] << " = " << x0[i] << '\n';
  std::cout << "\nbest solution: ";
  cmasols.print(std::cout,0,gp);
  std::cout << '\n';
  std::cout << "Expected Distance from Minimum: " << cmasols.edm() << '\n';
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";

  // Candidate bcand = cmasols.best_candidate();
  // double fmin = bcand.get_fvalue(); // min objective function value the optimizer converged to
  // std::vector<double> x_stdv = bcand.get_x(); // vector of objective function parameters at minimum.
  // const double* x_dptr = bcand.get_x_ptr(); // vector of objective function parameters at minimum, as C-style double array
  // Eigen::VectorXd x_ev = bcand.get_x_dvec(); // vector of objective function parameters at minimum, as Eigen vector
  // double edm = cmasols.edm(); // expected distance to the minimum.


  std::ofstream outfile{argv[2]};

  Candidate bcand = cmasols.best_candidate();
  double fmin = bcand.get_fvalue(); // min objective function value the optimizer converged to
  std::vector<double> x_stdv = bcand.get_x(); // vector of objective function parameters at minimum.
  const double* x_dptr = bcand.get_x_ptr(); // vector of objective function parameters at minimum, as C-style double array
  Eigen::VectorXd x_ev = bcand.get_x_dvec(); // vector of objective function parameters at minimum, as Eigen vector
  double edm = cmasols.edm(); // expected distance to the minimum.

  std::cout << "# m_chi   m_phi/m_chi    g_chi   pull-Oh2\nRESULT=";
  std::cout << input.m_chi.get() << '\t' << x_dptr[1] << '\t' << x_dptr[0] << '\t' << fmin << '\n';

  while(input.m_chi < end_m_chi)
  {
    input.m_chi = input.m_chi*stepsize;
    x0[0]=(fmin < 1.0e-3) ? x_dptr[0] : 0.3;
    x0[1]=(fmin < 1.0e-3) ? x_dptr[1] : 2.2;
    CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams1(x0,sigma,-1,0,gp); // -1 for automatically \
  decided lambda, 0 is for random	seeding	of the internal generator.                                              
    cmaparams1.set_algo(aCMAES);
    cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(f, cmaparams1);

    bcand = cmasols.best_candidate();
    fmin = bcand.get_fvalue(); // min objective function value the optimizer converged to
    x_stdv = bcand.get_x(); // vector of objective function parameters at minimum.
    x_dptr = bcand.get_x_ptr(); // vector of objective function parameters at minimum, as C-style double array
    x_ev = bcand.get_x_dvec(); // vector of objective function parameters at minimum, as Eigen vector
    edm = cmasols.edm(); // expected distance to the minimum.

    std::cout << "RESULT=" <<input.m_chi.get() << '\t' << x_dptr[1] << '\t' << x_dptr[0] << '\t' << fmin << '\n';
  }


  return cmasols.run_status();
}