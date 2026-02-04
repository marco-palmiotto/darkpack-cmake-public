#include "marty/models/pmssm_lem.h"
#include "marty/sgl/sgl.h"
#include <filesystem>
#include <locale>
#include <marty.h>

#include "marty/ufo/electriccharge.h"

#define LIBNAME "mssm2to2" // Name of the output library
#define FULL_PARAM


// Output files and paths
#define DIMCHEPSTRING 30 // String dimension for calcHep names' output

#define CORRESPONDANCE // Defining this variable the creation of correspondance.hpp is enabled
#define SQUAREDAMP     // Defining this variable squared amplitudes are calculated and added to the
//                        library. Disabling this behaviour can be useful if the user just wants
//                        to generate the auxiliary files for the library
#define UPDATE_PART_DETAILS

// #define DISABLE_LIB_CREATION

// #define DEBUG

using namespace mty;
using namespace mty::experimental;
using namespace mty::sm_input;
using namespace csl;

using Process = std::vector<Insertion>;
using Processes = std::vector<Process>;


bool setShow = false;

#include "src/processName.cpp"

#include "src/generateCorrespondance.cpp"

#include "src/checkCPsymmetric.cpp"

#include "src/computeAndAddToLibFromList.cpp"


// ******* MSSM-specific functions and types

// Input files
#define INPUT_PROCESSES "data/processes_feynman.psiso" // File that contains the 2to2 processes relevant in SuperIso
#define UNITARY_PROCESSES                                                                                              \
  "data/processes_unitary.psiso"           // File that contains only the processes better suited in the unitary gauge
#define NAMES_CLEAN "data/names_clean.dat" // File containing the names in both MARTY and SuperIso convention
#define INPUT_NOSISO                                                                                                   \
  "" // File containing the processes you want to compute which are not included
     // in SuperIso

#define SIZEPYSICALSM   20 // Number of physical particles in the SM + Higgses
#define SIZEPYSICALSUSY 28 // Number of physical BSM particles (no Higgses)
// Output files and paths
#define DIMCHEPSTRING 30 // String dimension for calcHep names' output
#define NAMECHEPFILE                                                                                                   \
  "data/processes_chep.txt"  // File that will contain a correspondance among the names of the
                             // processes in the calcHep convention (1st column) and in the
                             // superiso convention (2nd column)
#define MAXNUMPROCESSES 5000 // Max number of processes to be stored in a Processes variable
#include "src/mssm/addFromFile.cpp"

typedef struct
{
  std::vector<Process> process;
  std::vector<std::string> namemarty, namesiso;
  std::vector<mty::Order> order;   // Options mty::Order::TreeLevel, mty::Order::OneLoop,
  std::vector<bool> leading_order; // true N1 N1 -> gamma gamma
  std::vector<mty::gauge::Type> Wgauge;
} Process2to2ToComputeVec;

std::vector<Process2to2ToCompute> convertToVecofP(Process2to2ToComputeVec input)
{
  std::vector<Process2to2ToCompute> output;
  Process2to2ToCompute temp;
  output.reserve(input.process.size());

  for (size_t i = 0; i < input.process.size(); i++)
  {
    temp.process = input.process[i];
    temp.order = input.order[i];
    temp.leading_order = input.leading_order[i];
    temp.Wgauge = input.Wgauge[i];

    output.push_back(temp);
  }

  return output;
}



int main()
{
  std::cout << "Defining the model\n";
  PMSSM_LEM mssm;


  std::cout << "\nSuppressing verbose output during the calculation of the amplitudes\n";
  mty::option::verboseAmplitude = false;
  std::cout << "\nChoosing not to call gafed to show Feynman diagrams during the computation\n";
  std::cout << "of an amplitude. They will be printed in a directory\n";
  setShow = false;
  std::cout << "Creating the library\n";
  mty::Library lib(LIBNAME);

  std::cout << "Adding the flag for using posix threads in the C++ files in the library\n";
  lib.addLibrary("-lpthread");

  std::cout << "\nListing all the Physical particles\n";
  std::vector<Particle> part_0 = mssm.getPhysicalParticles([&](Particle p) { return p->isPhysical(); });
  std::vector<Particle> part;
  for (size_t i = 0; i != part_0.size(); i++)
  {
    if (!(IsOfType<GhostBoson>(part_0[i]) || IsOfType<GoldstoneBoson>(part_0[i])))
      part.push_back(part_0[i]);
  }

  std::cout << "Re-fixing mass names:\n";
  for (size_t i = 0; i != part.size(); i++)
  {
    std::string tempstring;
    tempstring.clear();
    if (!part[i]->getMass()->getName().empty())
    {
      part[i]->getMass()->setName("m_" + part[i]->getName());
      tempstring = part[i]->getMass()->getName();
      std::cout << part[i]->getName() << " " << tempstring << std::endl;
    }
  }

  std::cout << "Generating Spectrum..." << std::endl;
  lib.generateSpectrum(mssm);
  std::cout << "Spectrum Generated" << std::endl;

  // The following replacement is necessary to have no symbol dependencies in the electric charge
  mssm.replace((1 + pow_s(tan_s(theta_W), 2)), 1 / pow_s(cos_s(theta_W), 2));

  std::cout << "Separating SM and BSM particles" << std::endl;
  std::vector<Particle> psm, pbsm;

  for (size_t i = 0; i != part.size(); i++)
  {
    if (part[i]->getName() == "h" || part[i]->getName() == "W" || part[i]->getName() == "Z" ||
        part[i]->getName() == "G" || part[i]->getName() == "A" || part[i]->getName() == "d" ||
        part[i]->getName() == "u" || part[i]->getName() == "s" || part[i]->getName() == "c" ||
        part[i]->getName() == "b" || part[i]->getName() == "t" || part[i]->getName() == "e" ||
        part[i]->getName() == "mu" || part[i]->getName() == "tau" || part[i]->getName() == "nu_e" ||
        part[i]->getName() == "nu_mu" || part[i]->getName() == "nu_tau" || part[i]->getName() == "H0" ||
        part[i]->getName() == "Hp" || part[i]->getName() == "A0")
    {
      psm.push_back(part[i]);
    }
    else
    {
      pbsm.push_back(part[i]);
    }
  }

  char tempstring_c[1000];

  sprintf(tempstring_c, "[ -d \"auxiliary_library/%s\" ] || mkdir -p \"auxiliary_library/%s\"", LIBNAME, LIBNAME);

  generateCorrespondance(psm, pbsm, "smBsm.hpp");

  std::cout << "\nSetting widths symbolically\n";
  // This first group of widths is formed by the widths that are usually nonzero
  Expr Gamma_h = constant_s("Gamma_h");
  Expr Gamma_H0 = constant_s("Gamma_H0");
  Expr Gamma_Hp = constant_s("Gamma_Hp");
  Expr Gamma_W = constant_s("Gamma_W");
  Expr Gamma_Z = constant_s("Gamma_Z");
  Expr Gamma_top = constant_s("Gamma_top");
  Expr Gamma_A0 = constant_s("Gamma_A0");

  std::cout << "Constants have been defined\n";

  mssm.getParticle("h")->setWidth(Gamma_h);
  mssm.getParticle("H0")->setWidth(Gamma_H0);
  mssm.getParticle("Hp")->setWidth(Gamma_Hp);
  mssm.getParticle("W")->setWidth(Gamma_W);
  mssm.getParticle("Z")->setWidth(Gamma_Z);
  mssm.getParticle("t")->setWidth(Gamma_top);
  mssm.getParticle("A0")->setWidth(Gamma_A0);

  std::cout << "Setting secondary widths\n";
  // This second group of widths is formed by the widths that are usually zero
  Expr Gamma_gluino = constant_s("Gamma_gluino"); // This is 0 since the gluino is massless
  Expr Gamma_c1 = constant_s("Gamma_c1");
  Expr Gamma_c2 = constant_s("Gamma_c2");
  Expr Gamma_o1 = constant_s("Gamma_o1");
  Expr Gamma_o2 = constant_s("Gamma_o2");
  Expr Gamma_o3 = constant_s("Gamma_o3");
  Expr Gamma_o4 = constant_s("Gamma_o4");

  Expr Gamma_t1 = constant_s("Gamma_t1");
  Expr Gamma_t2 = constant_s("Gamma_t2");
  Expr Gamma_b1 = constant_s("Gamma_b1");
  Expr Gamma_b2 = constant_s("Gamma_b2");
  Expr Gamma_tau1 = constant_s("Gamma_tau1");
  Expr Gamma_tau2 = constant_s("Gamma_tau2");
  Expr Gamma_ul = constant_s("Gamma_ul");
  Expr Gamma_ur = constant_s("Gamma_ur");
  Expr Gamma_dl = constant_s("Gamma_dl");
  Expr Gamma_dr = constant_s("Gamma_dr");
  Expr Gamma_el = constant_s("Gamma_el");
  Expr Gamma_er = constant_s("Gamma_er");
  Expr Gamma_cl = constant_s("Gamma_cl");
  Expr Gamma_cr = constant_s("Gamma_cr");
  Expr Gamma_sl = constant_s("Gamma_sl");
  Expr Gamma_sr = constant_s("Gamma_sr");
  Expr Gamma_ml = constant_s("Gamma_ml");
  Expr Gamma_mr = constant_s("Gamma_mr");
  Expr Gamma_nuel = constant_s("Gamma_nuel");
  Expr Gamma_numl = constant_s("Gamma_numl");
  Expr Gamma_nutaul = constant_s("Gamma_nutaul");

  mssm.getParticle("G")->setWidth(Gamma_gluino);
  mssm.getParticle("C_1")->setWidth(Gamma_c1);
  mssm.getParticle("C_2")->setWidth(Gamma_c2);
  mssm.getParticle("N_1")->setWidth(Gamma_o1);
  mssm.getParticle("N_2")->setWidth(Gamma_o2);
  mssm.getParticle("N_3")->setWidth(Gamma_o3);
  mssm.getParticle("N_4")->setWidth(Gamma_o4);

  mssm.getParticle("st_1")->setWidth(Gamma_t1);
  mssm.getParticle("st_2")->setWidth(Gamma_t2);
  mssm.getParticle("sb_1")->setWidth(Gamma_b1);
  mssm.getParticle("sb_2")->setWidth(Gamma_b2);
  mssm.getParticle("stau_1")->setWidth(Gamma_tau1);
  mssm.getParticle("stau_2")->setWidth(Gamma_tau2);
  mssm.getParticle("su_L")->setWidth(Gamma_ul);
  mssm.getParticle("su_R")->setWidth(Gamma_ur);
  mssm.getParticle("sd_L")->setWidth(Gamma_dl);
  mssm.getParticle("sd_R")->setWidth(Gamma_dr);
  mssm.getParticle("se_L")->setWidth(Gamma_el);
  mssm.getParticle("se_R")->setWidth(Gamma_er);
  mssm.getParticle("sc_L")->setWidth(Gamma_cl);
  mssm.getParticle("sc_R")->setWidth(Gamma_cr);
  mssm.getParticle("ss_L")->setWidth(Gamma_sl);
  mssm.getParticle("ss_R")->setWidth(Gamma_sr);
  mssm.getParticle("smu_L")->setWidth(Gamma_ml);
  mssm.getParticle("smu_R")->setWidth(Gamma_mr);
  mssm.getParticle("snu_e")->setWidth(Gamma_nuel);
  mssm.getParticle("snu_mu")->setWidth(Gamma_numl);
  mssm.getParticle("snu_tau")->setWidth(Gamma_nutaul);



  std::cout << "\nGetting electric charges: (Feynman rules are needed at this purpose)\n";
  std::vector<double> electric_charges;
  for (const auto& particle : part)
  {
    // Create the necessary data for the module to extract the charge
    ufo::ElectricChargeMessenger elec{
        particle.get(),              // Particle from which we need the charge
        mssm.getParticle("A").get(), // Vector boson carrying the charge
        mssm,                        // Model
        ufo::SMData{}                // Model conventions
    };

    // The charge is read in the vertices
    ufo::Charge electronCharge = ufo::getElectricCharge(elec);
    // We convert it to a symbolic expression and display it
    csl::Expr symbCharge = csl::intfraction_s(electronCharge.num, electronCharge.denom);
    double float_charge = electronCharge.num * 1. / electronCharge.denom;
    std::cout << particle->getName() << " has charge: " << symbCharge << " = " << float_charge << std::endl;
    electric_charges.push_back(float_charge);
  }

  std::cout << "\nGetting degrees of freedom:\n";
  std::vector<int> dofs;
  for (const auto& particle : part)
  {
    dofs.push_back(particle->getNDegreesOfFreedom() * (particle->isSelfConjugate() ? 1 : 2));
    std::cout << particle->getName() << " has dof: " << dofs.back() << std::endl;
  }


  std::cout << "Reading the processes from file (MSSM)\n";

  // To this purpose, we use the type

  Process2to2ToComputeVec vec1;

  std::system("[ -f data/processes_chep.txt ] && rm -f data/processes_chep.txt");

  std::cout << "  Adding processes of the kind SUSY+SUSY -> SM+SM not in the list from superIso\n";
  std::string listNosiso = INPUT_NOSISO;
  if (listNosiso.size() != 0)
  {
    std::cout << "  The list is not empty: adding processes.\n\n";
    if (addFromFile(mssm, vec1.process, vec1.namemarty, vec1.namesiso, listNosiso) < 0)
    {
      std::cout << "  Something went wrong\n";
      exit(1);
    }
  }
  else
  {
    std::cout << "  No process in the list, moving forward.\n\n";
  }

  // We need to know how many of them there are
  const size_t n_nosiso = vec1.process.size();

  std::cout << "  Adding now the processes from the list of superIso\n\n";
  std::string listsuperiso = INPUT_PROCESSES;
  int temp = addFromFile(mssm, vec1.process, vec1.namemarty, vec1.namesiso, listsuperiso);
  if (temp <= 0)
  {
    std::cout << "  Something went wrong\n";
    exit(1);
  }

  // For the sake of clarity we save the total number of processes in a new variable
  const size_t n_siso = vec1.process.size();
#ifdef DEBUG
  std::cout << "  n_siso - (temp + n_nosiso) = " << n_siso - (temp + n_nosiso) << std::endl;
#endif

  // Setting Feynman gauge to each process as default

  for (size_t i = 0; i < n_siso; i++)
  {
    vec1.order.push_back(mty::Order::TreeLevel);
    vec1.Wgauge.push_back(mty::gauge::Type::Feynman);
    vec1.leading_order.push_back(false);
  }

  std::cout << "  Adding now the processes to be computed in the unitary gauge\n\n";
  std::string listUnitary = UNITARY_PROCESSES;
  if (listUnitary.size() != 0)
  {
    if (addFromFile(mssm, vec1.process, vec1.namemarty, vec1.namesiso, listUnitary) < 0)
    {
      std::cout << "  Something went wrong\n";
      exit(1);
    }
  }

  const size_t n_unitary = vec1.process.size();

  // Setting Unitary gauge to each process in the rest of the list
  for (size_t i = 0; i < n_siso; i++)
  {
    vec1.order.push_back(mty::Order::TreeLevel);
    vec1.Wgauge.push_back(mty::gauge::Type::Unitary);
    vec1.leading_order.push_back(false);
  }

  std::vector<Process2to2ToCompute> list_of_processes = convertToVecofP(vec1);

  std::cout << "We just read\n";

  for (size_t i = 0; i < list_of_processes.size(); i++)
  {
    std::cout << i << " ";
    list_of_processes[i].printname();
  }

  std::vector<Process2to2ToCompute> list_of_processes_1to2;

  computeAndAddToLibFromList(mssm, lib, list_of_processes, list_of_processes_1to2, "smBsm.hpp");
#ifndef DISABLE_LIB_CREATION

  // Computation of the widths at the three TreeLevel
  Expr width_h = mssm.computeWidth(Order::TreeLevel, "h");
  Expr width_W = mssm.computeWidth(Order::TreeLevel, "W");
  Expr width_Z = mssm.computeWidth(Order::TreeLevel, "Z");
  Expr width_t = mssm.computeWidth(Order::TreeLevel, "t");
  Expr width_H0 = mssm.computeWidth(Order::TreeLevel, "H0");
  Expr width_Hp = mssm.computeWidth(Order::TreeLevel, "Hp");
  Expr width_A0 = mssm.computeWidth(Order::TreeLevel, "A0");

  const std::array<Expr, 7> expr_to_expand = {width_h, width_W, width_Z, width_t, width_H0, width_Hp, width_A0};

  for (Expr expression : expr_to_expand)
  {
    auto temp = csl::Evaluated(expression, csl::eval::abbreviation);
    csl::ForEachLeaf(temp,
                     [&](csl::Expr& sub)
                     {
                       if (csl::IsConstant(sub) || csl::IsVariable(sub))
                       {
                         sub->setValue(CSL_UNDEF);
                       }
                     });
  }

  lib.addFunction("width_h", width_h);
  lib.addFunction("width_W", width_W);
  lib.addFunction("width_Z", width_Z);
  lib.addFunction("width_t", width_t);
  lib.addFunction("width_H0", width_H0);
  lib.addFunction("width_Hp", width_Hp);
  lib.addFunction("width_A0", width_A0);

  #ifdef SQUAREDAMP
  std::cout << "Writing the library on disk\n";
  lib.print();
  #endif
#endif

  return 0;
}
