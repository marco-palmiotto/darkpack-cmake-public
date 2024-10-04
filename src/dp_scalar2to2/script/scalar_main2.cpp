/* This is the first example file.
 * We start by including the standard files of the library
 * Then, we include all the new headers, in external_libs.h 
 */
#include <iostream>
#include <iomanip>
#include "config.hpp"
#include "correspondance.hpp" 
#include "leshouchesfrommarty.hpp"
#include "advmath.hpp"
#include "process.hpp"
#include "avgsvcalculator.hpp"
#include "boltzmann.hpp"

using namespace scalar2to2;
using namespace advmath;

scalar2to2::RunningSM run;

inline double sq(const double x)
{
    return x*x;
}

double gamma_phi_pred_onlyup(const Param_t &input)
{
  const double v_h = std::pow(2., -0.25)*std::pow(input.Gfermi, -0.5);
  
  std::vector<int> upquarks = {corr::u, corr::c, corr::t};
  double quarkcontrib = 0.;
  
  for( int quark : upquarks )
  {
    double mf= input.masses_vector[quark];
    double Y = std::sqrt(2.)*mf/v_h;
    quarkcontrib+= 3.*sq(Y*input.g_u)*input.m_phi/16./M_PI*std::pow( 1. - 4.*sq(mf/input.m_phi), 3./2.);
  }
  
  return quarkcontrib;
}

// Computation of Gamma_phi
double gamma_phi_pred(const Param_t &input)
{
//   const double v_h = std::pow(2., -0.25)*std::pow(input.Gfermi, -0.5);
  
  const double quarkcontrib = gamma_phi_pred_onlyup(input);
  
//   input.Print();
  
  const double chicontrib= sq(input.g_chi)*input.m_phi/8./M_PI*std::pow( 1. - 4.*sq(input.m_chi/input.m_phi), 3./2.);


//   printf("%s = %.3e ", "input.g_chi", double(input.g_chi));
//   
//   std::cout << "chicontrib=" << chicontrib << std::endl;
  
  const double strongcontrib=0.;
  
  return quarkcontrib + chicontrib + strongcontrib;
}

// Computation of the sum of the squared amplitudes from a given formula (the one I found)
double sum_Squaredampl_pred( Param_t &input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double C=sq(mf*gf*input.g_chi*input.e_em/(input.m_W*std::sin(input.theta_W)));
  // Defining s
  const double s = Ecm*Ecm;
  
  
  double result =Ncolors*C;
  
  // Writing the numerator of the result
  result*=(s-4*sq(input.m_chi))*(s-4*sq(mf));
  
  // Writing the denominator of the result
  result/=( sq(s-sq(input.m_phi)) + sq(input.m_phi*input.Gamma_phi));
  
  return result;
}

double dW_dcos_pred( Param_t &input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double M2 =sum_Squaredampl_pred(input, Ecm, Ncolors, gf, mf);
  
  const double s = sq(Ecm);
  const double p12= 0.5*std::sqrt(s-4.*sq(input.m_chi));
  const double p34= 0.5*std::sqrt(s-4.*sq(mf));
  const double peff = p12;
  
  double result = p12*p34*M2;
  result/= 8.*M_PI*sq(corr::part_hel_dof[corr::chi])*peff*Ecm;
 
  
  return result;
}

double dW_dcos_pred_nocoeff( Param_t &input, const double Ecm, const int Ncolors, const double gf, const double mf)
{
  run.HandleParamRunning(input, Ecm);
  // Defining prefactor
  const double M2 =sum_Squaredampl_pred(input, Ecm, Ncolors, gf, mf);
  
  const double s = sq(Ecm);
  const double p12= 0.5*std::sqrt(s-4.*sq(input.m_chi));
  const double p34= 0.5*std::sqrt(s-4.*sq(mf));
  
  return p12*p34*M2;
}


double Weff_pred_fromM2( Param_t &input, const double Ecm)
{
  run.HandleParamRunning(input, Ecm);
  
  std::cout << "\nIn Weff_pred_fromM2, " << input.masses_vector[corr::c] << " " << input.m_c << " " <<
               input.m_c_m_c << "\n";
  
  double sum=0.;
  
  const std::array<int, 9> 
        part({corr::u ,corr::c,  corr::t,
              corr::d, corr::s,  corr::b,
              corr::e, corr::mu, corr::tau}),
        part_color({3,3,3,
                    3,3,3,
                    1,1,1});
  const std::array<double, 9> 
        part_gf({input.g_u,input.g_u,input.g_u,
                input.g_d,input.g_d,input.g_d,
                input.g_l,input.g_l,input.g_l});
         
  size_t counter=0;
  
  for(auto i = 0 ; i <9; i++)
    {
      Process2to2 singleproc({corr::chi, corr::chi, part[i],part[i]},{false, true, false,true}); 
      if(singleproc.checkExistance())
      {
        double contrib=dW_dcos_pred(input, Ecm, part_color[i], part_gf[i], input.masses_vector[part[i]]);
        sum+=contrib;
        std::cout << "For " << singleproc.getName() << " dW/dcos(" << Ecm << ")= " << contrib << '\n';
        counter++;
      }
    }
    
  std::cout << "In the sum there were " << counter << " processes.\n";
  
  return 2*sum;
}

double Weff_formula( Param_t &input, const double Ecm)
{
  run.HandleParamRunning(input, Ecm);
  
//   std::cout << input;
  
//   std::cout << "Defining prefactor\n";
  const double C = sq(input.g_chi.get()*input.e_em/input.m_W/std::sin(input.theta_W));
  
//   std::cout << "Defining s\n";
  const double s = Ecm*Ecm;
  
  double result =C/(8.*M_PI*Ecm*sq(input.getLightestBSMdof()))*(s-4.*sq(input.m_chi));
  // Writing the denominator of the result
  result/=( sq(s-sq(input.m_phi)) + sq(input.m_phi*input.Gamma_phi));
  
  double sum=0.;
  
    const std::array<int, 9> part({corr::u, corr::c, corr::t,corr::d, corr::s, corr::b,corr::e, corr::mu, corr::tau}),
                        part_color({3,3,3,3,3,3,1,1,1});
    const std::array<double, 9> part_gf({input.g_u,input.g_u, input.g_u,input.g_d,input.g_d,input.g_d,input.g_l,input.g_l,input.g_l});
//   std::cout << "calculating the sum over all fermions\n";
  for(auto i = 0 ; i <9; i++)
    {
      Process2to2 singleproc({corr::chi, corr::chi, part[i],part[i]},{false, true, false,true}); 
      if(singleproc.checkExistance())
      {
        sum+=sq(1.0*input.masses_vector[part[i]]*part_gf[i])*part_color[i]*std::pow(s-4.*sq(input.masses_vector[part[i]]), 3./2.);
      }
    }
    
//   std::cout << "Multiplying the sum\n";
  result*=sum;
  
  return result;
}


int main(int argc, char ** argv)
{
 
    if(argc < 2)
    {
      std::cout << "This function requires 1 argument\n";
      return 1;
    }
    std::cout << std::setprecision(5) << std::scientific; 
    
    std::cout << "Reading inputs from " << argv[1] << std::endl;
    
    struct Param_t input(argv[1]);
    
    input.Print(); // Printing out the input structure
    
    run.init(input);
    run.RunCharmMass(false);
    
    std::cout << "v = " << 2*(input.m_W*std::sin(input.theta_W)/input.e_em);
    
    std::cout << "\nDefining setofproc\n";
    AvgSvCalculator allproc(input);
    allproc.setWeffcuts(false);
//     allproc.enablerunning(false);
    std::cout << "Class constructed\n";
    allproc.print();
    
    allproc.print_procs(std::cout, true);
    
    double T = 0.01*input.m_chi;
    
    // Computing single 2 to 2 sum of M^2
    double sqrts = 10.*input.m_chi;
    double costheta = 0.5;
    const std::array<int, 9> part({corr::u, corr::c, corr::t,corr::d, corr::s, corr::b,corr::e, corr::mu, corr::tau}),
                        part_color({3,3,3,3,3,3,1,1,1});
    const std::array<double, 9> part_gf({input.g_u,input.g_u, input.g_u,input.g_d,input.g_d,input.g_d,input.g_l,input.g_l,input.g_l});
    std::cout << std::setprecision(5) << std::scientific; 
    constexpr int width_field = 18;
    for(auto i = 0 ; i <9; i++)
    {
      Process2to2 singleproc({corr::chi, corr::chi, part[i],part[i]},{false, true, false,true}); 
      if(singleproc.checkExistance())
      {
        std::cout << '\n' << singleproc ;
        run.HandleParamRunning(input, sqrts);
        Param_t input_m=input;
        singleproc.setRunningData(&run);
        singleproc.setRunningExternal();
        double marty=singleproc.getSumSquaredAmpl(input_m, sqrts, costheta);
        double pred=sum_Squaredampl_pred(input, sqrts, part_color[i], part_gf[i], input.masses_vector[part[i]]);
        
        double martydweff=singleproc.get_g2_dweff_dcos(input_m, sqrts, 0.5)/sq(corr::part_hel_dof[corr::chi]);
        double predweff=dW_dcos_pred(input, sqrts, part_color[i], part_gf[i], input.masses_vector[part[i]]);
        
        double dweffdcoscontrib=singleproc.getDiffW12Contrib(input_m, sqrts, 0.5);
        double WeffContrib=singleproc.getTotalW12Contrib(input_m, sqrts);
        
            // Print the header row
        std::cout << std::setw(width_field) << ".getSumSquaredAmpl" << '\t' 
              << std::setw(width_field) << "sum_Squaredampl_pred" << '\t'
              << std::setw(width_field) << "1/2" << '\t'
              << std::setw(width_field) << ".get_g2_dweff_dcos/g2" << '\t'
              << std::setw(width_field) << "dW_dcos_pred" << '\t'
              << std::setw(width_field) << "4/5" << '\t'
              << std::setw(width_field) << ".getDiffW12Contrib" << '\t'
              << std::setw(width_field) << ".getTotalW12Contrib" << '\n';
        std::cout << std::setw(width_field) << marty << '\t' 
                  << std::setw(width_field) << pred << '\t'
                  << std::setw(width_field) << marty / pred << '\t'
                  << std::setw(width_field) << martydweff << '\t'
                  << std::setw(width_field) << predweff << '\t'
                  << std::setw(width_field) << martydweff / predweff << '\t'
                  << std::setw(width_field) << dweffdcoscontrib << '\t'
                  << std::setw(width_field) << WeffContrib << '\n';
                  
      }
    }
    
    run.HandleParamRunning(input, sqrts);
    std::cout << "dWeff/dcos = " << allproc.getdWeff_dcos(sqrts,0.) << std::endl; // With no multithreading
    std::cout << "2 x dWeff/dcos = " << 2.*allproc.getdWeff_dcos(sqrts,0.) << std::endl;
    std::cout << "Weff (num. integral) = " << allproc.getWeff(sqrts) << std::endl; // Uses multithreading
    

    std::cout << "Input for formula is\n"<< input;
    
    std::cout << "From formula = "   << Weff_formula(input, sqrts) << std::endl;
    double temp = Weff_pred_fromM2(input, sqrts);
    std::cout << "From formula M2= " << temp << std::endl;
    
    std::cout << "Ratio Weff/formula  = " << allproc.getWeff(sqrts)/Weff_formula(input, sqrts)<< std::endl;  
    
    
//     std::cout << "\nWeff_1234 = " << allproc.getWeff_fromW12(sqrts) << std::endl;
    
    double sigmav = allproc.getAverageSigmav(T);
    std::cout << " - <sigma v>( " << T << " ) = " << sigmav << " GeV^-2 = " << units_conversion::GeVm2_to_cm3_over_s*sigmav << " cm^3/s\n\n";
    
    
    auto sigmav_pred = [] (const double Temp, const int Nc, const int part_l,const double gf, Param_t &input_l) mutable
    {
        run.HandleParamRunning(input_l, input_l.getLightestBSMmass());
        //Higgs vev in the SM 
        const double v_h = std::pow(2., -0.25)*std::pow(input_l.Gfermi, -0.5);
        const double Y_f = std::sqrt(2.)*input_l.masses_vector[part_l]/v_h;
        // The temperature corresponds to the center of mass energy
        double numerator = Nc*3.
                           *sq(input_l.g_chi*gf*Y_f)
                           *input_l.m_chi
                           *Temp
                           *std::pow( 1. - sq(input_l.masses_vector[part_l]/input_l.m_chi), 3./2.);
        double denom = 8*M_PI
                       *(  sq(sq(input_l.m_phi) - 4.*sq(input_l.m_chi)) 
                         + sq(input_l.m_phi*input_l.Gamma_phi));
        
        return numerator/denom;     
    };
    
    
    double sigmavpred = 0.;
    
    for( auto count_part : {corr::e, corr::mu, corr::tau} )
      sigmavpred +=sigmav_pred(T, 1, count_part, input.g_l, input);
      
    std::cout << "Computed sigmav for leptons\n";
      
    for( auto count_part : {corr::d, corr::s, corr::b}) 
      sigmavpred += sigmav_pred(T, 3, count_part,input.g_d, input);
      
    std::cout << "Computed sigmav for d-type\n";
      
    for( auto count_part : {corr::u, corr::c, corr::t}) 
      sigmavpred += sigmav_pred(T, 3, count_part,input.g_u, input);
    
    std::cout << "Computed sigmav for u-type\n";
      

    std::cout << "The predicted sigmav is " << sigmavpred << " GeV^-2\n";
    std::cout << "The ratio is " << sigmavpred/sigmav << std::endl;
    
    double gamma_phi_up_trial = width_phi(input).real();
    std::cout << " Gamma_phi_up_computed = " << gamma_phi_up_trial ;
    std::cout << "\n Gamma_phi_up_pred = " << gamma_phi_pred(input);
    
    
    std::ofstream fout{"Weff.out"};
    if(!fout)
    {
      std::cerr << "Impossible to open Weff.out\n";
      exit(1);      
    }
    
    std::cout << "\n\nOpened Weff.out\n";
    
    fout.precision(17);
    
     
    std::vector<double> Ecmv;
    Ecmv.reserve(allproc.getWefftabsize());
    for(size_t i = 0 ; i < allproc.getWefftabsize() ; i++)
    {
      
      double peff, g2_weff;
      allproc.get_g2_WeffTabElement(i, peff, Ecmv[i], g2_weff);
      
      const double Ecm = Ecmv[i];
      
      run.HandleParamRunning(input, Ecm);
      
      fout << std::scientific 
           << Ecm
           << '\t' 
           << 0.5*Ecm*Ecm/allproc.getMassLBSM()/allproc.getMassLBSM()-1.
           << '\t' 
           << g2_weff/sq(corr::part_hel_dof[corr::chi]) 
           << '\t' 
           << Weff_formula(input, Ecm) 
           << '\t' 
           << allproc.getWeff(Ecm)       
           << '\n';
    }
    
    fout.close();
    
    fout.open("diffWeffTop.out");
    if(!fout)
    {
      std::cerr << "Impossible to open diffWeffTop.out\n";
      exit(1);      
    }
    
    // Testing a single process
    Process2to2 singleproc({corr::chi, corr::chi, corr::t,corr::t},{false, true, false,true}); 
    
    for(auto Ecm : Ecmv)
    {
      Param_t input_m=input;
      run.HandleParamRunning(input, Ecm);
      
      fout << std::scientific << 
              Ecm << '\t' 
           << 0.5*Ecm*Ecm/allproc.getMassLBSM()/allproc.getMassLBSM()-1. << '\t' 
           << singleproc.get_process_dependent_contrib_dweff_dcos(input_m, Ecm, 0.) << '\t' 
           << dW_dcos_pred(input, Ecm, part_color[corr::t], input.g_u, input.masses_vector[corr::t])<< '\n';
    }
    
    
    fout.close();
    double smallest_double = std::numeric_limits<double>::min();
    double smallest_positive_double = std::numeric_limits<double>::min() * std::numeric_limits<double>::epsilon();

    std::cout << "\nSmallest double: " << smallest_double << std::endl;
    std::cout << "Absolute value of the smallest positive double: " << smallest_positive_double << std::endl;
    
    
    
    fout.open("scalarsigmav.out");
    if(!fout)
    {
      std::cerr << "Impossible to open scalarsigmav.out\n";
      exit(1);      
    }
    
    std::cout << "\nOpened scalarsigmav.out\n";
    
    
    for(T=1.0e-3 ; T<= 100.*input.getLightestBSMmass() ; T*=1.1)
    {
      fout << std::scientific 
           << T/input.getLightestBSMmass()
           << '\t'
           << allproc.getAverageSigmav_coan_hightemp(T)
           << '\t'
           << allproc.getAverageSigmav_coan_lowtemp(T) 
           << '\t';
           
//       sigmavpred = 0.;
//       for( auto count_part : {corr::u, corr::c, corr::t, corr::d, corr::s, corr::b}) sigmavpred += sigmav_pred(T, 3, count_part, input);
//       
//       for( auto count_part : {corr::e, corr::mu, corr::tau} ) sigmavpred +=sigmav_pred(T, 1, count_part, input);
//       
        sigmavpred = 0.;
        for( auto count_part : {corr::e, corr::mu, corr::tau} ) sigmavpred +=sigmav_pred(T, 1, count_part, input.g_l, input);
        for( auto count_part : {corr::d, corr::s, corr::b}) sigmavpred += sigmav_pred(T, 3, count_part,input.g_d, input);
        for( auto count_part : {corr::u, corr::c, corr::t}) sigmavpred += sigmav_pred(T, 3, count_part,input.g_u, input);
        
    
      fout << sigmavpred
           << '\n';     
      
    }
      
    fout.close();
    
    
    // Computing derivatives of Weff
    
    allproc.print_coefficients();
    
    sigmavpred = 0.; T=1.;
    for( auto count_part : {corr::e, corr::mu, corr::tau} ) sigmavpred +=sigmav_pred(T, 1, count_part, input.g_l, input);
    for( auto count_part : {corr::d, corr::s, corr::b}) sigmavpred += sigmav_pred(T, 3, count_part,input.g_d, input);
    for( auto count_part : {corr::u, corr::c, corr::t}) sigmavpred += sigmav_pred(T, 3, count_part,input.g_u, input);;
    
    std::cout << "Taylor_pred_1 =  " << sigmavpred*input.getLightestBSMmass() << '\n'
              << "Ratio ourTaylor/FormulaTaylor  = " << allproc.get_TaylorCoeffSVT(1)/(sigmavpred*input.getLightestBSMmass()) 
              << '\n';

  
    
    fout.open("YYeq.dat");
    if(!fout)
    {
      std::cerr << "Impossible to open YYeq.dat\n";
      exit(1);      
    }
    
    
    BoltzmannSolver boltz(allproc);
    
    
    
    
    fout.close();
    
    return 0;
}
