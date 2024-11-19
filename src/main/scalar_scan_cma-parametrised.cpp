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

double function_to_minimize(const double *parameters, [[maybe_unused]] const int npar=3)
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

  std::cout << "Running " << argv[0] << '\n';

  constexpr const double sigma=omega_h2_err;
  constexpr const int npar=3;
  constexpr const int n_required_args=npar+3;
  constexpr const char *par_names[npar] = {"g_chi", "m_chi", "m_phi"};
  // parameters in input will be
  // input.g_chi = parameters[0];
  // input.m_chi = parameters[1];
  // input.m_phi = parameters[2];

  if(argc < n_required_args)
  {
    std::cout << "This function requires " << n_required_args << " arguments:\n"
              << " - The name of the input file\n"
              << " - The name of the output file\n"
              << " - The seed for the parameter g_chi\n"
              << " - The seed for the parameter m_chi\n"
              << " - The seed for the parameter m_phi";
    return 1;
  }
  filename=argv[1];
  
  std::cout << "The input file is " << filename << '\n';
  double lbounds[npar]={1.e-4, 1.   , 1. }, 
         ubounds[npar]={1.0  , 1.e+4, 1.e+5}; // arrays for lower and upper parameter bounds, respectively                

  std::vector<double> x0;
  x0.reserve(npar);
  for(int i = 0 ; i < npar; i++)
  {
    const int index_arg=i+3;
    char *end;
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
  CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0,sigma,-1,0,gp); // -1 for automatically decided lambda, 0 is for random	seeding	of the internal generator.                                              
  cmaparams.set_algo(aCMAES);
  FitFunc f = &function_to_minimize;
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
  if(!outfile)
  {
    std::cerr << "Impossible to open " << argv[2] << " exiting\n";
    return 0;
  }

  outfile << "seed is\n";
  for(int i=0; i < npar ; i++) outfile << par_names[i] << " = " << x0[i] << '\n';
  outfile << "\nbest solution: ";
  cmasols.print(outfile,0,gp);
  outfile << '\n';
  outfile << "Expected Distance from Minimum: " << cmasols.edm() << '\n';
  outfile << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";


  return cmasols.run_status();
}