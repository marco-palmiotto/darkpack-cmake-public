#include "mssm2to2.h"
#include "correspondance.hpp"
#include "advmath.hpp"
#include "params_new.hpp"
#include "input_manip.hpp"
#include "leshouches.h"
#include "leshouchesfromsuperiso.hpp"
#include "process.hpp"
#include "avgsvcalculator.hpp"
#include "RunningSM.hpp"
#include "avgsvcalculator.hpp"
using namespace mssm2to2;

// #define EMIN 419.009
// #define EMAX 20000
// #define DELTAE  0.1
// #define INC  150

#define MAX_ENERGY_RESOLUTION 4


#ifndef GEVTOPBARN
  #define GEVTOPBARN (0.389379*1.0e+9)
#endif

//Mathematical macros
#ifndef SQUARE
  #define SQUARE(X) (X*X)
#endif

//Output macros
#ifndef PRECISION
  #define PRECISION 10
#endif

#ifndef SPACE
  #define SPACE 20
#endif 

const std::string PATHPLOTS=std::string(OUTPATH)+"MARTY/" ;
const std::string PATHPLOTSRUN=std::string(OUTPATH)+"MARTY_run/" ;

#ifndef NUMBER_OF_CALLS
  #define NUMBER_OF_CALLS 100000
#endif

#ifndef NUMTHREADS
  #define NUMTHREADS 16
#endif

#ifdef gauss_order
#undef gauss_order
#endif
#define gauss_order 8

#define DEBUG

static std::string getSname(const Process2to2 &proc)
{
  //printdebug("Entering getSname");
  using namespace corr;
    const std::string listSisoSusy[SIZEPHYSICALBSM]={
        "o1", "o2", "o3", "o4",
        "c1", "c2", 
        "t1", "t2", 
        "b1", "b2",
        "dr", "dl",
        "ur", "ul",
        "cr", "cl",
        "sr", "sl",
        "er", "el",
        "mr", "ml",
        "l1", "l2",
        "ne", "nm", "nl",
        "go"};

    const std::string listSisoSM[SIZEPHYSICALSM]={
        "h", "hh", "hc", "h3",
        "w", "z", "g", "a", 
        "d", "u", "s", "c", "b", "t",
        "e", "m", "l",
        "ne", "nm", "nl"};
        
    const std::string listMartySusy[SIZEPHYSICALBSM]={
        "N_1", "N_2", "N_3", "N_4",
        "C_1", "C_2",
        "st_1", "st_2",
        "sb_1", "sb_2",
        "sd_R", "sd_L",
        "su_R", "su_L",
        "sc_R", "sc_L",
        "ss_R", "ss_L",
        "se_R", "se_L",
        "smu_R", "smu_L",
        "stau_1", "stau_2",
        "snu_e", "snu_mu", "snu_tau",
        "sG"};

    const std::string listMartySM[SIZEPHYSICALSM]={
        "h", "H0", "Hp", "A0",
        "W", "Z", "G", "A",
        "d", "u", "s", "c", "b", "t",
        "e", "mu", "tau",
        "nu_e", "nu_mu", "nu_tau"};
        
   std::ostringstream sout;    
   //printdebug("Arrays and stringstream initialised");
   for(int i : {0, 1})
   {
     int j=-1;
     do
     {
       j++;
       //printdebug("About to compare "+part_names[proc.getField(i)]+" and "+listMartySusy[j]);
     }while( part_names[proc.getField(i)] != listMartySusy[j] );
     sout << listSisoSusy[j] ;
     if( !proc.getMatter(i) ) sout << "bar";
   }
   for(int i : {2, 3})
   {
     int j=-1;
     do
     {
       j++;
       //printdebug("About to compare "+part_names[proc.getField(i)]+" and "+listMartySusy[j]);
     }while(part_names[proc.getField(i)] != listMartySM[j]);
     sout << listSisoSM[j] ;
     if( !proc.getMatter(i) ) sout << "bar";
   }

   std::string name = sout.str();
   return name;
}

AvgSvCalculator *allprocsptr;

std::vector<double> abscissas;


// For a single process, writes the NON-RUNNED 
// - Sum of the squared amplitudes 
// - (g_DM)^2 dW/dcos contribution 
void test_getSingle_ext(Process2to2 p, mssm2to2::Param_t input, bool wantrun=false)
{
  // Tests a single process trough the whole range of energies 
  std::vector<double> ctheta = {-0.3, 0.0, 0.3, 0.5};
  
  char tempstringc[2000];
  char pathplots[1000];
  if(wantrun) sprintf(pathplots, "%s", PATHPLOTSRUN.c_str());
  else sprintf(pathplots, "%s", PATHPLOTS.c_str());
  
  sprintf(tempstringc, "[ -d \"%s\" ] || mkdir -p \"%s\"", pathplots, pathplots);
  std::system(tempstringc);

  const std::string fileprefix = pathplots+static_cast<std::string>("dec_")+p.getMname();
  RunningSM *runptr;
  if(wantrun) runptr = new RunningSM(input);
  for( auto cosine : ctheta )
  {
    std::cout << "Computing " << p.getName() << " for costheta = " << cosine << std::endl;
      
    sprintf(tempstringc, "%.1f", cosine);
    const std::string filename=fileprefix+"_sumsquaredamp"+tempstringc+".dat";
    std::ofstream outf{filename};
    if (!(outf))
    {
        // Print an error and exit
        std::cerr << "Uh oh, the file " << filename << " could not be opened for writing!" << std::endl;
        do
        {
            sleep(0.5);
            outf.open(filename, std::ios::out);
        }while( !(outf) ) ;     
    }
    outf << std::setprecision(PRECISION) << std::scientific; 

    //Creating the file to write the results
    const std::string filenamet=fileprefix+"_contrib"+tempstringc+".dat";
    std::ofstream outft{filenamet};
    if (!(outft))
    {
      // Print an error and exit
      std::cerr << "Uh oh, the file " << filenamet << " could not be opened for writing!" << std::endl;
      do
      {
        sleep(0.5);
        outft.open(filenamet, std::ios::out);
      }while( !(outft) );
    }
    outft << std::setprecision(PRECISION) << std::scientific; 
    
    //Cycling on the energies
    for(double Ecm : abscissas)
    {
      if(wantrun) runptr->HandleParamRunning(input, Ecm);
      double SumSquaredAmpl = p.getSumSquaredAmpl(input, Ecm, cosine);
      // In SuperIso it is the fomula multiplied by 8 pi
      double contrib = p.get_g2_dweff_dcos(input, Ecm, cosine)*8*M_PI;
      if(SumSquaredAmpl <0.) continue;
      outf  << Ecm << "\t" << SumSquaredAmpl << std::endl;
      outft << Ecm << "\t" << contrib << std::endl;
    }    
  }
  if(wantrun) delete runptr;
}

// For each process in setofproc, writes the NON-RUNNED 
// - Sum of the squared amplitudes 
// - (g_DM)^2 dW/dcos contribution 
void test_ext_norun_local(AvgSvCalculator setofproc)
{
  std::cout << "For each process in setofproc, writes the NON-RUNNED\n";
  std::cout << "- Sum of the squared amplitudes \n";
  std::cout << "- (g_DM)^2 dW/dcos contribution \n";

  for( size_t i = 0 ; i < setofproc.getN() ; i++)
  {
    std::async(std::launch::async, test_getSingle_ext, setofproc.getProcess(i), setofproc.getInput(), false);
  }
}

void test_ext_yesrun_local(AvgSvCalculator setofproc)
{
  std::cout << "For each process in setofproc, writes the RUNNED\n";
  std::cout << "- Sum of the squared amplitudes \n";
  std::cout << "- (g_DM)^2 dW/dcos contribution \n";

  for( size_t i = 0 ; i < setofproc.getN() ; i++)
  {
    std::async(std::launch::async, test_getSingle_ext, setofproc.getProcess(i), setofproc.getInput(), true);
  }
}

int main(int argc, char ** argv)
{
 
    Timer t;
    if(argc < 1)
    {
      std::cout << "This function requires 1 argument\n";
      return 1;
    }
    std::cout << std::setprecision(5) << std::scientific; 
    
    std::cout << "Reading inputs from " << argv[1] << std::endl;
    
    struct Param_t input; // We need a Param_t struct to store data in the format we're going to use
    int err;
    superisosupport::ReadLHA(input, argv[1], &err);
//     input = readmodule::ReadLHA(argv[1]);
    if(err != 0) return err;
    
    std::cout << std::setprecision(5) << std::scientific; 
    input.Print(); 
    
    mssm2to2::Param_t input_ckm1 = input, input_mumass0 = input;
    
    input_manip::CKMIdentity(input_ckm1);
    
    input_mumass0.m_e = 0.;
    input_mumass0.m_mu = 0.;
    
    std::cout << "Reading abscissas from file\n";
    FILE * fptr;
    std::string tempstring = std::string(OUTPATH)+"Weff_siso.dat";
    fptr = fopen(tempstring.c_str(), "r");
    if(fptr == nullptr)
    {
      std::cout << "Impossible to open " << tempstring << std::endl;
      exit(1);
    }
    abscissas.clear();
    do
    {
      double sqrtS;
      
      char buffer[100];
      for(size_t i = 0 ; i < 12 ; i++) fgets(buffer, 100,fptr);
      fscanf(fptr, "%lf%*lf", &sqrtS);
      if(sqrtS < 0.99*input_ckm1.getLightestBSMmass() )
      {
        std::cerr << "Error: sqrtS too small\n";
        continue ;
      }
      if(abscissas.size()>0)
      {
        if(sqrtS < abscissas.back() )
        {
         std::cout << "You're trying to add a smaller value\n";
         fclose(fptr);
         exit(-1);          
        }
//         while(sqrtS-abscissas.back() > MAX_ENERGY_RESOLUTION)
//         {
//         double value =  abscissas.back()+MAX_ENERGY_RESOLUTION;
//         abscissas.push_back(value);
//         }
      }
      abscissas.push_back(sqrtS);
    }while(!feof(fptr));
    fclose(fptr);    
    
    std::cout << "Read " << abscissas.size() << " abscissas:\n";
    for(double Ecm : abscissas) std::cout << Ecm << std::endl;
    
    
    std::cout << "Defining setofproc\n";
    allprocsptr = new AvgSvCalculator(input_ckm1);
    std::cout << "Class constructed\n";
#ifdef DEBUG
      std::cout << "Defining namesptr\n";
#endif 
    FILE *namesptr;
    tempstring=std::string(OUTPATH)+"names_clean.dat";
#ifdef DEBUG
      std::cout << "Opening " << tempstring << "\n";
#endif 
    namesptr = fopen(tempstring.c_str(), "w");
    if(namesptr == nullptr)
    {
      std::cout << "Impossible to open names_clean.dat" << std::endl;
      return 1;
    }
#ifdef DEBUG
      std::cout << "success\n";
#endif 
    for(size_t i = 0 ; i < allprocsptr->getN() ; i++)
    {
#ifdef DEBUG
      std::cout << "Writing " ;
      std::cout << allprocsptr->getProcess(i).getMname() ;
      std::cout << " " ; 
      std::cout << getSname(allprocsptr->getProcess(i)) << std::endl;
#endif
      fprintf(namesptr, "%s %s\n", 
                ("dec_"+allprocsptr->getProcess(i).getMname()).c_str(), 
                getSname(allprocsptr->getProcess(i)).c_str()
            );
    }
    
    fclose(namesptr);
#ifdef DEBUG
      std::cout << "File closed\n";
#endif 
      
    RunningSM *run = new RunningSM(input);
    
    std::cout << "Total processes = " << allprocsptr->getN() << std::endl;
    
    if( allprocsptr->getN() == 0) return 1;
    
    AvgSvCalculator trial(*allprocsptr);
    
    delete allprocsptr;
    
    allprocsptr = &trial;
    
    Timer t1;
//     double elapsed_prev = 0.;
    for(size_t j = 0 ; j < allprocsptr->getN() ; j++)
    {
      std::cout << "Computing runned total xsec for " << allprocsptr->getProcess(j).getName() << " n = " << j << " / " << allprocsptr->getN() << std::endl;
      double elapsed = t1.getElapsed();
      double total = (elapsed)/(j+1)*(allprocsptr->getN());
      double remaining = total - elapsed;
      std::cout << "Elapsed time: " << (int)(elapsed)/60 <<" m " << (int)(elapsed)%60 << " s\n";
      std::cout << "Remaining time is: " << (int)(remaining)/60 <<" m " << (int)(remaining)%60 << " s\n";
      std::cout << "Total time is: " << (int)(total)/60 <<" m " << (int)(total)%60 << " s\n";
      Process2to2 proc = allprocsptr->getProcess(j);
      
      proc.setRunningData(run);
      proc.setRunningInternal();
      
//       test_getSingle_int(proc, input);
      
      std::cout << "Defing filename\n"; 
      std::string filename = PATHPLOTSRUN+static_cast<std::string>("dec_")+proc.getMname()+"_xsec_run.dat";
      std::cout << "Creating the output file\n";
      std::ofstream outfx;
      outfx.open(filename, std::ios::out);
      if (!(outfx))
      {
          // Print an error and exit
          std::cerr << "Uh oh, the file " << filename << " could not be opened for writing!" << std::endl;
          do
          {
              sleep(0.5);
              outfx.open(filename, std::ios::out);
          }while( !(outfx) );
      }
      outfx << std::setprecision(5) << std::scientific;
      
      std::cout << filename << " created\n";      
      
      for(double Ecm : abscissas)
      {
// #define DEBUG
// #ifdef DEBUG
//         std::cout << "Working at Ecm = " << Ecm <<std::endl;;
// #endif
        double err=0;
        double xsec = GEVTOPBARN*proc.getTotalCrossSection(input_mumass0, Ecm, &err);
        err*=xsec;
        outfx << Ecm << "\t" << xsec << '\t' << err << std::endl;
      }

      outfx.close();
    }//End of cross section calculation 
    
    std::cout << "Calling test_ext_norun\n";
    test_ext_norun_local(*allprocsptr);
    
    std::cout << "Calling test_ext_yesrun\n";
    test_ext_yesrun_local(*allprocsptr);
    
    tempstring=std::string(OUTPATH)+"Weff_new.dat";
    std::ofstream outweff{tempstring.c_str()};
    
    if(!outweff)
    {
      std::cerr << "Impossible to open Weff_new.dat\n";
      exit(1);      
    }
    std::cout << "File Weff_new.dat opened\n";
//     double mlsp = allprocsptr->getMassLSP();
    outweff << std::setprecision(5) << std::scientific; 
    for(size_t i = 0 ; i < abscissas.size() ; i++)
    {
      double sqrtS, value;
      int count = 0;
      do
      {
        sqrtS=abscissas[i]*(1.+1.0e-4*count++);
        value=allprocsptr->get_g2_Weff(sqrtS);
      }while(!std::isnormal(value));      
      outweff << sqrtS << "\t" << value << std::endl;
    }
      
    delete run;
    return 0;
}
