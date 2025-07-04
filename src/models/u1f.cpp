#include <marty.h>
#include <marty/ufo/electriccharge.h>

#include "src/def_su2f.cpp"

using namespace mty;
using namespace mty::experimental;
using namespace mty::sm_input;
using namespace csl;


class U1f_Model : public virtual SU2f_Model
{
  public:
  U1f_Model(const int model) : SU2f_Model(model, false), SM_Model(false)
  {
    // Deactivating charged bosons
    mty::ModelData::removeParticle("chi_p");
    mty::ModelData::removeParticle("V");
    mty::ModelData::removeParticle("chi_dm_2");

    SU2f_Model::refresh();
  }
};

#include "marty/ufo/electriccharge.h"

#define LIBNAME         "u1f2to2" // Name of the output library
#define MAXNUMPROCESSES 5000      // Max number of processes to be stored in a Processes variable

#define CORRESPONDANCE // Defining this variable the creation of correspondance.h is enabled
#define SQUAREDAMP     // Defining this variable squared amplitudes are calculated and added to the
//                        library. Disabling this behaviour can be useful if the user just wants
//                        to generate the auxiliary files for the library

#define FULL_PARAM
#define DISPLAYSQUAREDAMPLITUDESEXPRESSIONS
// #define DEBUG

template <typename T> inline void printdebug(const T& s1)
{
#ifdef DEBUG
  std::cout << s1 << std::endl;
#endif
  return;
};


using Process = std::vector<Insertion>;
using Processes = std::vector<Process>;


bool setShow = false;

#include "src/processName.cpp"

#include "src/generateCorrespondance.cpp"

#include "src/checkCPsymmetric.cpp"

#include "src/computeAndAddToLibFromList.cpp"


int main()
{
  U1f_Model my_model(23);

  std::cout << my_model;


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
  std::vector<Particle> part_0 = my_model.getPhysicalParticles([&](Particle p) { return p->isPhysical(); });
  std::vector<Particle> part;
  std::cout << "Printing Physical particles before filtering\n";
  for (size_t i = 0; i != part_0.size(); i++)
  {
    std::cout << part_0[i]->getName() << '\n';
    if (!(IsOfType<GhostBoson>(part_0[i]) || IsOfType<GoldstoneBoson>(part_0[i])))
      part.push_back(part_0[i]);
  }

  std::cout << "Printing Physical particles after filtering\n";
  for (auto elem : part)
    std::cout << elem->getName() << '\n';

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
  lib.generateSpectrum(my_model);
  std::cout << "Spectrum Generated" << std::endl;

  // The following replacement is necessary to have no symbol dependencies in the electric charge
  my_model.replace((1 + pow_s(tan_s(theta_W), 2)), 1 / pow_s(cos_s(theta_W), 2));

  std::cout << "Separating SM and BSM particles" << std::endl;
  std::vector<Particle> psm;
  std::vector<Particle> pbsm({my_model.getParticle("chi_dm_1")});
  for (size_t i = 0; i != part.size(); i++)
  {
    bool istoadd = true;
    for (size_t j = 0; j != pbsm.size() && istoadd; j++)
    {
      istoadd = part[i]->getName() != pbsm[j]->getName();
    }
    if (istoadd)
      psm.push_back(part[i]);
  }

  if (pbsm.size() > 0)
    generateCorrespondance(psm, pbsm);

  std::cout << "\nSetting widths symbolically\n";
  Expr Gamma_h = constant_s("Gamma_h");
  Expr Gamma_W = constant_s("Gamma_W");
  Expr Gamma_Z = constant_s("Gamma_Z");
  Expr Gamma_top = constant_s("Gamma_top");
  Expr Gamma_V3 = constant_s("Gamma_V_3");

  std::cout << "Constants have been defined\n";

  my_model.getParticle("h")->setWidth(Gamma_h);
  my_model.getParticle("W")->setWidth(Gamma_W);
  my_model.getParticle("Z")->setWidth(Gamma_Z);
  my_model.getParticle("t")->setWidth(Gamma_top);
  my_model.getParticle("V_3")->setWidth(Gamma_V3);


  std::cout << "\nGetting electric charges: (Feynman rules are needed at this purpose)\n";
  std::vector<double> electric_charges;
  for (const auto& particle : part)
  {
    // Create the necessary data for the module to extract the charge
    ufo::ElectricChargeMessenger elec{
        particle.get(),                  // Particle from which we need the charge
        my_model.getParticle("A").get(), // Vector boson carrying the charge
        my_model,                        // Model
        ufo::SMData{}                    // Model conventions
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
    std::cout << particle->getName() << " has dof: " << particle->getNDegreesOfFreedom() << std::endl;
  }

  std::vector<Process2to2ToCompute> list_of_processes;


  // Our list is formed by all the processes of the kind Chibar_{1,2,3} Chi_{1,2,3} -> Fbar F
  // where F is a SM lepton

  std::cout << "Building the list of SM leptons\n";

  for (size_t i1 = 0; i1 < pbsm.size(); i1++)
  {
    for (mty::Insertion state1 : std::array<mty::Insertion, 2>({Incoming(pbsm[i1]), Incoming(AntiPart(psm[i1]))}))
    {
      for (size_t i2 = i1; i2 < pbsm.size(); i2++)
      {
        for (mty::Insertion state2 : std::array<mty::Insertion, 2>({Incoming(pbsm[i2]), Incoming(AntiPart(pbsm[i2]))}))
        {
          for (size_t j1 = 0; j1 < psm.size(); j1++)
          {
            for (mty::Insertion state3 :
                 std::array<mty::Insertion, 2>({Outgoing(psm[j1]), Outgoing(AntiPart(psm[j1]))}))
            {
              for (size_t j2 = j1; j2 < psm.size(); j2++)
              {
                for (mty::Insertion state4 :
                     std::array<mty::Insertion, 2>({Outgoing(psm[j2]), Outgoing(AntiPart(psm[j2]))}))
                {
                  Process2to2ToCompute proc_c;

                  proc_c.process = {(state1), (state2), (state3), (state4)};
                  proc_c.order = mty::Order::TreeLevel;
                  proc_c.leading_order = false;
                  proc_c.Wgauge = mty::gauge::Type::Feynman;

                  list_of_processes.push_back(proc_c);
                }
              }
            }
          }
        }
      }
    }
  }

  {
    int i = 0;
    for (const auto& single_proc : list_of_processes)
    {
      const std::string procname = processName(single_proc.process);
      std::cout << "Process n " << i++ << " / " << list_of_processes.size() << " : " << procname << '\n';
    }
  }

  computeAndAddToLibFromList(my_model, lib, list_of_processes);


  // Computation of the widths at the three TreeLevel
  Expr width_h = my_model.computeWidth(Order::TreeLevel, "h");
  Expr width_W = my_model.computeWidth(Order::TreeLevel, "W");
  Expr width_Z = my_model.computeWidth(Order::TreeLevel, "Z");
  Expr width_t = my_model.computeWidth(Order::TreeLevel, "t");
  Expr width_V3 = my_model.computeWidth(Order::TreeLevel, "V_3");


  std::vector<Expr> expr_to_expand = {width_h, width_W, width_Z, width_t, width_V3};

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
  lib.addFunction("width_V_3", width_V3);
  lib.addDefaultParameter("g_s", false);

#ifdef SQUAREDAMP
  std::cout << "Writing the library on disk\n";
  lib.print();
#endif

  return 0;
}
