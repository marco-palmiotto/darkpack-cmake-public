#include "all_libs.hpp"
namespace mssm2to2{
std::string getSname(const Process2to2 &proc)
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
} // end of namespace mssm2to2
