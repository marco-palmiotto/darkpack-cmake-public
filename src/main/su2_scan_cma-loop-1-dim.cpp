#include <iostream>
#include <iomanip>
#include "libcmaes/cmaes.h"
#include <libcmaes/esoptimizer.h>
#include <libcmaes/cmastrategy.h>
#include <libcmaes/llogging.h>
#include <fstream>

#include "dp_su2f2to2/config.hpp"
#include "dp_su2f2to2/correspondance.hpp" 
#include "dp_su2f2to2/leshouchesfrommarty.hpp"
#include "dp_su2f2to2/process.hpp"
#include "dp_su2f2to2/avgsvcalculator.hpp"
#include "dp_su2f2to2/boltzmann.hpp"

using namespace su2f2to2;

using namespace advmath;
using namespace libcmaes;

const std::string PATHPLOTS=std::string(OUTPATH);
constexpr const double omega_h2_target=0.12;
constexpr const double omega_h2_err=0.00091;

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
  int read_args=0;
  std::cout << "Running " << argv[read_args++] << '\n';

  constexpr const double sigma=omega_h2_err;
  constexpr const int npar=1;
  constexpr const int n_required_args=8;
  constexpr const char *par_names[npar] = {"g_chi"}; // Name of the parameters to be optimised
  constexpr const double m_phi_over_mv=2.;

  [[maybe_unused]] constexpr const double init_seed_g_chi=0.3;

  double lbounds[npar]={1.e-6}, 
         ubounds[npar]={1.0  }; // arrays for lower and upper parameter bounds, respectively  

  std::cout << "This function requires " << n_required_args-1 << " arguments:\n"
              << " - The name of the input file\n"
              << " - The name of the output file\n"
              << " - The seed for the parameter g_f\n"
              << " - The value for the parameter m_chi/m_V\n"
              << " - The initial value for m_V\n"
              << " - The final value for m_V\n"
              << " - The desired numbewr of points\n";

  if(argc < n_required_args)
  {
    return 1;
  }
  
  Param_t input(argv[read_args]);
  std::cout << "The input file is " << argv[read_args++] << '\n';    
  
  std::ofstream outfile{argv[read_args]};
  if(!outfile)
  {
    std::cerr << "Impossible to open " << argv[read_args] << '\n';
    return 1;
  }
  std::cout << "Output file: " << argv[read_args++] << '\n';

  // Initialising x0
  std::vector<double> x0;
  x0.reserve(npar);
  char *end;
  for(int i = 0 ; i < npar; i++)
  {
    const double value=strtod(argv[read_args],&end);
    if(value < lbounds[i] || value > ubounds[i])
    {
      std::cerr << "Error in argument: " << read_args << ", parameter " << i <<  ", " << par_names[i] << '\n'
                << "value " << value << " is out of the bound interval"
                << " [ " << lbounds[i] << ", " << ubounds[i] << " ].\n";
      return 2;
    }
    read_args++;
    // std::cout << "x0.size()=" << x0.size() << ", emplacing back " << value << '\n';
    x0.emplace_back(value);
  }

  if(x0.size() != npar)
  {
    std::cerr << "Initial parameters have wrong size: x0.size()=" << x0.size() << ", npar=" << npar << '\n';
    return 2;
  } 

  std::cout << "SEED: ";
  for(size_t i =0; i < x0.size(); i++) std::cout << par_names[i] << '=' << x0[i] << ' ';
  std::cout << '\n';

  // Reading parameters from CLI arguments
  const double mchi_over_mv  = std::strtod(argv[read_args++],&end),
               init_m_V = std::strtod(argv[read_args++],&end),
               end_m_V  = std::strtod(argv[read_args++],&end);
  const int    npoints  = std::strtod(argv[read_args++],&end);

  // Computing the stepsize
  const double stepsize = std::pow(end_m_V/init_m_V, 1./npoints);   

  std::cout << "init_m_V=" << init_m_V << '\n'
            << "end_m_V=" << end_m_V  << '\n'
            << "npoints=" << npoints  << '\n';

  // The following vector will contain each value of m_V to be varied
  // Note: not stricly necessary now, but good to set it up for parallelising maybe
  std::vector<double> abscissae(npoints);
  for(int i=0 ; i <npoints ; i++) abscissae[i]=init_m_V*std::pow(stepsize,i);

  std::cout << "Abscissas are:";
  for(int i=0 ; i <npoints ; i++) std::cout << ' ' << abscissae[i] ;
  std::cout << '\n';
  
  // Function and paramters' first definition
  input.g_f=x0[0];
  input.m_V=abscissae[0];
  input.m_phi=2.*abscissae[0];
  input.m_chi_dm_1=mchi_over_mv*abscissae[0];
  input.m_chi_dm_2=mchi_over_mv*abscissae[0];
  input.refresh();
  BoltzmannSolver boltz(input);

  auto function_to_minimize = [&] (const double *parameters, [[maybe_unused]] const int size_array=2)
  {
    // parameters in input will be
    input.g_f= parameters[0];
    input.refresh();
    boltz.changeInput(input);
    const double omega=boltz.relic_density();
    return std::abs(omega - omega_h2_target);
  };     

  // Minimization algorithm

  GenoPheno<pwqBoundStrategy> gp(lbounds,ubounds,npar); // genotype / phenotype transform associated to bounds.  

  // -1 for automatically decided lambda, 0 is for random	seeding	of the internal generator.          
  CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0,sigma,-1,0,gp);                                     
  cmaparams.set_algo(aCMAES);
  FitFunc f = function_to_minimize;
  CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(f, cmaparams);

  std::cout << "seed is\n";
  for(int i=0; i < npar ; i++) std::cout << par_names[i] << " = " << x0[i] << '\n';
  std::cout << "\nbest solution: ";
  cmasols.print(std::cout,0,gp);
  std::cout << '\n';
  std::cout << "Expected Distance from Minimum: " << cmasols.edm() << '\n';
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n ";

  Candidate bcand = cmasols.get_best_seen_candidate();

  double fmin = bcand.get_fvalue(); // min objective function value the optimizer converged to
  Eigen::VectorXd bestparameters_dvec = gp.pheno(cmasols.get_best_seen_candidate().get_x_dvec()); // Applying pheno to the best candidate output
  double edm = cmasols.edm(); // expected distance to the minimum.
  int return_val=cmasols.run_status();
  
  std::cout << "# m_phi/m_V  m_chi/m_V m_V   g_f   pull-Oh2   EDM  exit_code  \nRESULT=";
  std::cout << m_phi_over_mv << '\t' <<  mchi_over_mv << '\t' << input.m_V.get() << '\t' << bestparameters_dvec[0] << '\t' << fmin << '\t' << edm << '\t' << return_val << '\n';

  outfile << "# m_phi/m_V  m_chi/m_V m_V  g_f   pull-Oh2   EDM  exit_code \n";
  outfile << m_phi_over_mv << '\t' <<  mchi_over_mv << '\t' << input.m_V.get() << '\t' << bestparameters_dvec[0] << '\t' << fmin  << '\t' << edm << '\t' << return_val << '\n';

  for(int i = 1 ; i < npoints; i++)
  {
    input.m_V = abscissae[i];
    input.m_phi= m_phi_over_mv*input.m_V;
    x0[0]=(fmin < 1.0e-3) ? bestparameters_dvec[0] : 0.1;
    CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams1(x0,sigma,-1,0,gp); // -1 for automatically decided lambda, 0 is for random	seeding	of the internal generator.                                              
    cmaparams1.set_algo(aCMAES);
    cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(f, cmaparams1);

    bcand = cmasols.get_best_seen_candidate();
    fmin = bcand.get_fvalue(); // min objective function value the optimizer converged to
    bestparameters_dvec = gp.pheno(bcand.get_x_dvec()); // Applying pheno to the best candidate output
    edm = cmasols.edm(); // expected distance to the minimum.

    std::cout << "RESULT=" 
              << m_phi_over_mv << '\t' 
              << mchi_over_mv << '\t' 
              << input.m_V.get()<< '\t' 
              << bestparameters_dvec[0] << '\t' 
              << fmin  << '\t' 
              << edm << '\t' 
              << return_val << '\n';

    outfile << m_phi_over_mv << '\t' 
            << mchi_over_mv << '\t' 
            << input.m_V.get()<< '\t' 
            << bestparameters_dvec[0] << '\t' 
            << fmin  << '\t' 
            << edm << '\t' 
            << return_val << '\n';
    return_val=cmasols.run_status();
  }


  return return_val;
  
}