// This function just generates the file smBsm.hpp
// The purpose of such a file is to be read in the 
// numerical library, in order for it to know which 
// particles are SM and which are BSM 
std::string generateCorrespondance(std::vector<Particle> psm, std::vector<Particle> pbsm, const std::string filename="smBsm.hpp") 
{        
    std::string folder_path="auxiliary_library/"+static_cast<std::string>(LIBNAME);
    std::string tempstring = "[ -d \""+folder_path+"\" ] || mkdir -p \""+folder_path+"\"";
//     std::cout << "tempstring is " << tempstring << std::endl;exit(0);
    
    std::system(tempstring.c_str());
  
    std::ofstream fileout{folder_path+static_cast<std::string>("/")+filename};
    
    fileout << "inline const int SIZEPHYSICALSM = " << psm.size()     << ";\n"; 
    fileout << "inline const int SIZEPHYSICALBSM = " << pbsm.size()  << ";\n";  
    fileout << "inline const int TOTAL_PARTICLES = " << pbsm.size()+psm.size() << ";\n"; 
    
    fileout << "inline const std::array<int, SIZEPHYSICALBSM> bsm_particles{";
    for( size_t i = 0 ; i < pbsm.size()-1; i++ )
    {
      fileout << "\n  " << pbsm[i]->getName() << ",";      
    }
      fileout << "\n  " << pbsm.back()->getName() << "\n};\n";
    
    
    
    
    fileout << "inline const std::array<int,SIZEPHYSICALSM> sm_particles{";
    for( size_t i = 0 ; i < psm.size()-1; i++ )
    {
      fileout << "\n  " << psm[i]->getName() << ",";      
    }
      fileout << "\n  " << psm.back()->getName() << "\n};\n";
    

    return filename;
}
