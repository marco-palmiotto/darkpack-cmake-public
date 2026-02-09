int computeAndAddToLibFromList(mty::Model& model, // model
                               mty::Library& lib, // output library
                               std::vector<Process2to2ToCompute> listofprocs,
                               std::vector<Process2to2ToCompute> listofprocs_1to2,
                               std::string nameSmBsmFile = "smBsm.hpp" // file to include in correspondance.hpp
)
{
  size_t count_converted = 0; // Counts the number of computed amplitudes: needed to track progress
  char tempstring_c[300];
  int temp;
  std::string prevname1, prevname2; // Strings that contain the previous names of the particles in the initial state.
                                    // It is used to regroup the Feynman diagrams in the output files

  int Sf34;              // Symmetry factor for the 3rd and 4th particle
  int C1234;             // CP symmetry factor
  csl::Expr comb_factor; // csl expression for the line above

  std::cout << "Creating directory to save Feynman diagrams\n\n";
  std::vector<std::shared_ptr<mty::wick::Graph>> graphs;

  sprintf(tempstring_c, "[ -e %s_graphs ] || mkdir %s_graphs", lib.getName().c_str(), lib.getName().c_str());
  std::system(tempstring_c);

  std::cout << "Listing all the Physical particles\n";
  std::vector<std::string> pnames, mnames; // particle name and its mass name
  std::vector<Particle> part_0 = model.getPhysicalParticles([&](Particle p) { return p->isPhysical(); });
  std::vector<Particle> part;
  for (size_t i = 0; i != part_0.size(); i++)
  {
    if (!(IsOfType<GhostBoson>(part_0[i]) || IsOfType<GoldstoneBoson>(part_0[i])))
      part.push_back(part_0[i]);
  }
  for (size_t i = 0; i != part.size(); i++)
  {
    if (!part[i]->getMass()->getName().empty())
    {
      mnames.push_back("static_cast<double>(input." + part[i]->getMass()->getName() + ")");
    }
    else
    {
      mnames.push_back("0.0");
    }
    pnames.push_back(part[i]->getName());
  }

  std::cout << "\nGetting electric charges:\n";
  std::vector<double> electric_charges;
  for (const auto& particle : part)
  {
    // Create the necessary data for the module to extract the charge
    ufo::ElectricChargeMessenger elec{
        particle.get(),               // Particle from which we need the charge
        model.getParticle("A").get(), // Vector boson carrying the charge
        model,                        // Model
        ufo::SMData{}                 // Model conventions
    };

    // The charge is read in the vertices
    ufo::Charge electronCharge = ufo::getElectricCharge(elec);
    // We convert it to a symbolic expression and display it
    csl::Expr symbCharge = csl::intfraction_s(electronCharge.num, electronCharge.denom);
    double float_charge = electronCharge.num * 1. / electronCharge.denom;
    std::cout << particle->getName() << " has charge: " << symbCharge << " = " << float_charge << std::endl;
    electric_charges.push_back(float_charge);
  }

  std::cout << "\nGetting total degrees of freedom:\n";
  std::vector<int> dofs;
  dofs.reserve(part.size());
  for (const auto& particle : part)
  {
    dofs.push_back(particle->getNDegreesOfFreedom() * (particle->isSelfConjugate() ? 1 : 2));
    std::cout << particle->getName() << " has dof: " << dofs.back() << std::endl;
  }

  std::cout << "\nGetting helicity degrees of freedom:\n";
  std::vector<int> h_dofs;
  h_dofs.reserve(part.size());
  for (const auto& particle : part)
  {
    h_dofs.push_back(particle->getNDegreesOfFreedom());
    std::cout << particle->getName() << " has helicity dofs: " << h_dofs.back() << std::endl;
  }

  std::cout << "\nGetting if a particle has an antipatticle:\n";
  std::vector<bool> isSelfConjugate;
  for (const auto& particle : part)
  {
    isSelfConjugate.push_back(particle->isSelfConjugate());
    std::cout << particle->getName() << (isSelfConjugate.back() ? "is self-conjugate" : "is not self-conjugate")
              << std::endl;
  }

#ifdef CORRESPONDANCE
  std::cout << "\nCreating correspondance file\n";
  std::ofstream correspondancefileh;

  //     std::cout << "Building the directory if needed\n";
  //     sprintf(tempstring_c, "[ -d \"%s/include\" ] || mkdir -p \"%s/include\" ; [ -d \"%s/src\" ] || mkdir -p
  //     \"%s/src\"",
  //             lib.getName().c_str(), lib.getName().c_str(),lib.getName().c_str(),lib.getName().c_str());
  //     std::system( tempstring_c );

  std::cout << "Building the directory if needed\n";
  sprintf(tempstring_c, "[ -d \"auxiliary_library/%s\" ] || mkdir -p \"auxiliary_library/%s\"", lib.getName().c_str(),
          lib.getName().c_str());
  std::system(tempstring_c);

  sprintf(tempstring_c, "auxiliary_library/%s/correspondance.hpp", lib.getName().c_str());
  correspondancefileh.open(tempstring_c);

  if (!correspondancefileh)
  {
    std::cerr << "Impossible to write the correspondance.hpp file\n";
    exit(1);
  }
  std::cout << "Creating correspondance.hpp file\n";
  correspondancefileh << "#pragma once\n\n";
  correspondancefileh << "#ifndef __cplusplus\n";
  correspondancefileh << "#define __cplusplus\n";
  correspondancefileh << "#endif\n\n";
  correspondancefileh << "#include \"" << lib.getName() << ".h\"\n";
  correspondancefileh << "#include \"clib_" << lib.getName() << ".h\"\n";
  correspondancefileh << "#include <unordered_map>\n";
  correspondancefileh << "namespace " << lib.getName() << "{\n";
  correspondancefileh
      << "// We are using char 32 (space) as possible separator, char 126 (tilde) as antiParticle tag\n";
  correspondancefileh << "#define EMPTYCHAR 32\n";
  correspondancefileh << "#define ANTICHAR 126\n";
  correspondancefileh << "using Cfptr_t = ccomplex_return_t (*) (cparam_t const *);\n";
  correspondancefileh << "using CXXfptr_t = complex_t (*) (param_t const *);\n";
  correspondancefileh << "namespace corr{\n";
  correspondancefileh << "  using Entry_t = std::tuple<Cfptr_t, short int, short int>;\n";
  correspondancefileh << "  using Entry_t_1to2 = std::tuple<CXXfptr_t, short int, short int>;\n";
  correspondancefileh << "  enum Part_t { \n";
  #ifdef DEBUG
  std::cout << "Writing particle names in enum Part_t\n";
  #endif
  correspondancefileh << "        " << pnames[0] << "=1,\n";
  #ifdef DEBUG
  std::cout << "i = 0\n";
  #endif
  for (size_t i = 1; i != part.size() - 1; i++)
  {
  #ifdef DEBUG
    std::cout << "i = " << i << " Writing " << pnames[i] << std::endl;
  #endif
    correspondancefileh << "        " << pnames[i] << ",\n";
  }
  correspondancefileh << "       " << pnames.back() << std::endl;
  correspondancefileh << "   };\n\n";

  correspondancefileh << "#include \"" << nameSmBsmFile << "\"\n";
  #ifdef DEBUG
  std::cout << "Writing particle names in part_names\n";
  #endif
  correspondancefileh << "inline const std::array<std::string,TOTAL_PARTICLES+1> part_names={\" \",\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    correspondancefileh << "        \"" << pnames[i] << "\",\n";
  }
  correspondancefileh << "       \"" << pnames.back() << '\"' << std::endl;
  correspondancefileh << "   };\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle statistic\n";
  #endif

  correspondancefileh << "inline const std::array<bool,TOTAL_PARTICLES+1> isboson={false, \n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    correspondancefileh << "        " << (part[i]->isBosonic() ? "true" : "false") << ",\n";
  }
  correspondancefileh << "       " << (part.back()->isBosonic() ? "true" : "false") << std::endl;
  correspondancefileh << "   };\n\n";


  #ifdef DEBUG
  std::cout << "Writing particle codes\n";
  #endif
  correspondancefileh << "inline const std::string mass_codes={EMPTYCHAR, \n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    correspondancefileh << "        EMPTYCHAR+" << pnames[i] << ",\n";
  }
  correspondancefileh << "       EMPTYCHAR+" << pnames.back() << std::endl;
  correspondancefileh << "   };\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle charge\n";
  #endif
  correspondancefileh << "inline const std::array<double,TOTAL_PARTICLES+1> part_charge={ 0,\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    correspondancefileh << "        " << electric_charges.at(i) << ",\n";
  }
  correspondancefileh << "        " << electric_charges.back() << "\n};\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle total degrees of freedom\n";
  #endif
  correspondancefileh << "inline const std::array<int,TOTAL_PARTICLES+1> part_tot_dof={ 0,\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    correspondancefileh << "        " << dofs.at(i) << ",\n";
  }
  correspondancefileh << "        " << dofs.back() << "\n};\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle helicity degrees of freedom\n";
  #endif
  correspondancefileh << "inline const std::array<int,TOTAL_PARTICLES+1> part_hel_dof={ 0,\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    correspondancefileh << "        " << h_dofs.at(i) << ",\n";
  }
  correspondancefileh << "        " << h_dofs.back() << "\n};\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle conjugation\n";
  #endif
  correspondancefileh << "inline const std::array<bool,TOTAL_PARTICLES+1> part_isSelfConj={ false,\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    correspondancefileh << "        " << (isSelfConjugate.at(i) ? "true" : "false") << ",\n";
  }
  correspondancefileh << "        " << (isSelfConjugate.back() ? "true" : "false") << "\n};\n\n";


  correspondancefileh << "\nextern const std::unordered_map<std::string, Entry_t> squaredampl;\n";
  correspondancefileh << "\nextern const std::unordered_map<std::string, Entry_t_1to2> squaredampl_1to2;\n";

  #ifdef DEBUG
  std::cout << "Writing getMassFirst\n";
  #endif
  correspondancefileh << "\ninline double getMassFirst(const int p, const param_t &input)\n";
  correspondancefileh << "{\n";

  correspondancefileh << "    switch (p)\n";
  correspondancefileh << "    {\n";

  for (size_t i = 0; i != part.size(); i++)
  {
    correspondancefileh << "        case(" << pnames[i] << ") : \n";
    correspondancefileh << "            return " << mnames[i] << ";\n";
  }
  correspondancefileh << "        default : \n";
  correspondancefileh << "            return -1.;\n";
  correspondancefileh << "    }\n";
  correspondancefileh << "}// End of getMassFirst\n";

  correspondancefileh << "}//End of namespace corr\n";
  correspondancefileh << "}//End of namespace " << lib.getName() << "\n";
  // correspondancefileh << "#endif\n";
  correspondancefileh.close();

  std::cout << "correspondance.h has been created succesfully.\n";

  #ifndef SQUAREDAMP
    #ifdef UPDATE_PART_DETAILS
  return 0;
    #endif
  #endif

  #ifndef DISABLE_LIB_CREATION

  std::cout << "Creating initialise_map.cpp file\n";
  sprintf(tempstring_c, "auxiliary_library/%s/initialise_map.cpp", lib.getName().c_str());
  correspondancefileh.open(tempstring_c);
  if (!correspondancefileh)
  {
    std::cerr << "Impossible to write the " << tempstring_c << " file\n";
    return 1;
  }

  correspondancefileh << "#include \"correspondance.hpp\"\n"
                      << "namespace " << lib.getName() << "::corr{\n"
                      << "const std::unordered_map<std::string, Entry_t> squaredampl{ \n";

  #endif // End of ifdef CORRESPONDANCE

  // Assigning the 1st values to the strings with the names
  // of the previous processes (starting with 0)
  prevname1 = listofprocs[0].process[0].getField()->getName();
  prevname2 = listofprocs[0].process[1].getField()->getName();

  const size_t lastproc = listofprocs.size();
  for (size_t i = 0; i != lastproc; i++)
  {
    count_converted++;
    std::string procname = processName(listofprocs[i].process);
    std::string nameSumSq = "sumSqAmpl_" + procname;
    std::string nameCombFac = "combFac_" + procname;
    std::cout << "Process n " << i + 1 << " / " << lastproc << "\n";

    // Adding the process in the library
    std::cout << "Computing " << nameSumSq << " ";
    // Setting W boson gauge
    model.getParticle("W")->setGaugeChoice(listofprocs[i].Wgauge);

    // Creating a lambda function that returns the amplitude at the desired order
    auto lambda_ampl = [&]()
    {
      if (!listofprocs[i].leading_order)
        return model.computeAmplitude(listofprocs[i].order, listofprocs[i].process);
      // If you want it at the leading order, let us compute the tree-level expression
      auto ampl_tree = model.computeAmplitude(mty::Order::TreeLevel, listofprocs[i].process);
      if (!ampl_tree.empty())
        return ampl_tree;
      auto ampl_one_loop = model.computeAmplitude(mty::Order::OneLoop, listofprocs[i].process);
      return ampl_one_loop;
    };

    auto ampl = lambda_ampl();
    if (ampl.empty())
    {
      std::cout << "... It seems this process is forbidden at One Loop order! Let's go on!\n";
      std::cerr << "Found forbidden process in the list : " << procname << std::endl;
      continue;
    }
    std::cout << std::endl;
    Expr s1;
    // Saving the single diagram
    SaveDiagrams(lib.getName() + "_graphs/" + procname + ".json", ampl);
    // Grouping diagrams and saving them
    if (listofprocs[i].process[0].getField()->getName() != prevname1 ||
        listofprocs[i].process[1].getField()->getName() != prevname2)
    {
      // If one of the two particles is different, I save the file and I start a new one
      SaveDiagrams(lib.getName() + "_graphs/" + prevname1 + "_" + prevname2 + ".json", graphs);
      prevname1 = listofprocs[i].process[0].getField()->getName();
      prevname2 = listofprocs[i].process[1].getField()->getName();
      graphs.clear();
    }
    auto process_graphs = ampl.obtainGraphs();
    graphs.insert(graphs.end(), process_graphs.begin(), process_graphs.end());

  #ifdef SQUAREDAMP
    s1 = model.computeSquaredAmplitude(ampl, false);

    if (s1 == CSL_0)
    {
      std::cout << "The squared amplitude is null\n";
      continue;
    }

    // You need to evaluate all the abbreviation to have the full set of parameters
    // as input data
    auto temp = csl::Evaluated(s1, csl::eval::abbreviation);
    #ifdef FULL_PARAM
    csl::ForEachLeaf(temp,
                     [&](csl::Expr& sub)
                     {
                       if (csl::IsConstant(sub) || csl::IsVariable(sub))
                       {
                         sub->setValue(CSL_UNDEF);
                       }
                     });
    #endif
    if (temp->dependsOn(dirac4.C_matrix.get()))
    {
      std::cerr << "Found conjugation matrix dependence for " << nameSumSq << std::endl;
      exit(2);
    }
    lib.addFunction(nameSumSq, s1);
    #ifdef DISPLAYSQUAREDAMPLITUDESEXPRESSIONS
    std::cout << "|M|^2 = " << DeepRefreshed(temp) << std::endl;
    #endif
  #endif
    // computing the combinatorial factor
    comb_factor = ampl.getKinematics().getDegreesOfFreedomFactor();
    if (!csl::IsNumerical(comb_factor))
    {
      std::cerr << "The number of degrees of freedom is not numerical!\n";
      exit(3);
    }
    comb_factor = Evaluated(comb_factor, csl::eval::all);

    std::string mass_string = " ";
    for (size_t j = 0; j < 4; j++)
    {
      Expr tempmass = listofprocs[i].process[j].getField()->getMass();
      std::string tempstring = (tempmass != CSL_0) ? "static_cast<double>(input." + tempmass->getName() + ")" : "0.0";
      mass_string += tempstring + ",";
    }

    // Determining the symmetry factor
    Sf34 =
        (((listofprocs[i].process[2].isOutgoingParticle() && listofprocs[i].process[3].isOutgoingParticle()) ||
          (listofprocs[i].process[2].isOutgoingAntiParticle() && listofprocs[i].process[3].isOutgoingAntiParticle())) &&
         listofprocs[i].process[2].getField()->getName() == listofprocs[i].process[3].getField()->getName())
            ? 2
            : 1;

    // Determining the CP symmetry factor
    C1234 = checkCPsymmetric(model, listofprocs[i].process) ? 2 : 1;
  #ifdef CORRESPONDANCE
    // Writing on the correspondance file
    correspondancefileh << "{ {";
    for (int j = 0; j <= 2; j++)
    {
      if (!listofprocs[i].process[j].isParticle())
        correspondancefileh << "ANTICHAR,";
      correspondancefileh << "EMPTYCHAR+corr::" << listofprocs[i].process[j].getField()->getName() << ",";
    }
    if (!listofprocs[i].process[3].isParticle())
      correspondancefileh << "ANTICHAR,";
    correspondancefileh << "EMPTYCHAR+corr::" << listofprocs[i].process[3].getField()->getName() << "},";
    // Writing the value field
    correspondancefileh << "{ ";
    // Writing the functions
    correspondancefileh << " &c_" << nameSumSq << ", ";
    correspondancefileh << comb_factor << ",";
    // Writing the number
    correspondancefileh << C1234;
    correspondancefileh << "} }";
    if (i != lastproc - 1)
      correspondancefileh << ",\n";
  #endif // End of ifdef CORRESPONDANCE
  } // Ends the cycle on the processes
  #ifdef CORRESPONDANCE
  correspondancefileh << "};\n}\n";
  correspondancefileh.close();
  #endif
  #ifndef SQUAREDAMP
    #ifdef UPDATE_PART_DETAILS
  return 0;
    #endif
  #endif

  #ifndef DISABLE_LIB_CREATION

  std::cout << "Creating initialise_map_1to2.cpp file\n";
  sprintf(tempstring_c, "auxiliary_library/%s/initialise_map_1to2.cpp", lib.getName().c_str());
  correspondancefileh.open(tempstring_c);
  if (!correspondancefileh)
  {
    std::cerr << "Impossible to write the " << tempstring_c << " file\n";
    return 1;
  }

  correspondancefileh << "#include \"correspondance.hpp\"\n"
                      << "namespace " << lib.getName() << "::corr{\n"
                      << "const std::unordered_map<std::string, Entry_t_1to2> squaredampl_1to2{ \n";

  #endif // End of ifdef CORRESPONDANCE

  // Assigning the 1st values to the strings with the names
  // of the previous processes (starting with 0)
  prevname1 = listofprocs_1to2[0].process[0].getField()->getName();
  prevname2 = listofprocs_1to2[0].process[1].getField()->getName();

  const size_t lastproc_1to2 = listofprocs_1to2.size();
  for (size_t i = 0; i != lastproc_1to2; i++)
  {
    count_converted++;
    std::string procname = processName(listofprocs_1to2[i].process);
    std::string nameSumSq = "sumSqAmpl_" + procname;
    std::string nameCombFac = "combFac_" + procname;
    std::cout << "Process n " << i + 1 << " / " << lastproc_1to2 << "\n";

    // Adding the process in the library
    std::cout << "Computing " << nameSumSq << " ";
    // Setting W boson gauge
    model.getParticle("W")->setGaugeChoice(listofprocs_1to2[i].Wgauge);

    // Creating a lambda function that returns the amplitude at the desired order
    auto lambda_ampl = [&]()
    {
      if (!listofprocs_1to2[i].leading_order)
        return model.computeAmplitude(listofprocs_1to2[i].order, listofprocs_1to2[i].process);
      // If you want it at the leading order, let us compute the tree-level expression
      auto ampl_tree = model.computeAmplitude(mty::Order::TreeLevel, listofprocs_1to2[i].process);
      if (!ampl_tree.empty())
        return ampl_tree;
      auto ampl_one_loop = model.computeAmplitude(mty::Order::OneLoop, listofprocs_1to2[i].process);
      return ampl_one_loop;
    };

    auto ampl = lambda_ampl();
    if (ampl.empty())
    {
      std::cout << "... It seems this process is forbidden at One Loop order! Let's go on!\n";
      std::cerr << "Found forbidden process in the list : " << procname << std::endl;
      continue;
    }
    std::cout << std::endl;
    Expr s1;
    // Saving the single diagram
    SaveDiagrams(lib.getName() + "_graphs/" + procname + ".json", ampl);
    // Grouping diagrams and saving them
    if (listofprocs_1to2[i].process[0].getField()->getName() != prevname1 ||
        listofprocs_1to2[i].process[1].getField()->getName() != prevname2)
    {
      // If one of the two particles is different, I save the file and I start a new one
      SaveDiagrams(lib.getName() + "_graphs/" + prevname1 + "_" + prevname2 + ".json", graphs);
      prevname1 = listofprocs_1to2[i].process[0].getField()->getName();
      prevname2 = listofprocs_1to2[i].process[1].getField()->getName();
      graphs.clear();
    }
    auto process_graphs = ampl.obtainGraphs();
    graphs.insert(graphs.end(), process_graphs.begin(), process_graphs.end());
  #ifdef SQUAREDAMP
    s1 = model.computeSquaredAmplitude(ampl, false);

    if (s1 == CSL_0)
    {
      std::cout << "The squared amplitude is null\n";
      continue;
    }

    // You need to evaluate all the abbreviation to have the full set of parameters
    // as input data
    auto temp = csl::Evaluated(s1, csl::eval::abbreviation);
    #ifdef FULL_PARAM
    csl::ForEachLeaf(temp,
                     [&](csl::Expr& sub)
                     {
                       if (csl::IsConstant(sub) || csl::IsVariable(sub))
                       {
                         sub->setValue(CSL_UNDEF);
                       }
                     });
    #endif
    if (temp->dependsOn(dirac4.C_matrix.get()))
    {
      std::cerr << "Found conjugation matrix dependence for " << nameSumSq << std::endl;
      exit(2);
    }
    lib.addFunction(nameSumSq, s1);
    #ifdef DISPLAYSQUAREDAMPLITUDESEXPRESSIONS
    std::cout << "|M|^2 = " << DeepRefreshed(temp) << std::endl;
    #endif
  #endif
    // computing the combinatorial factor
    comb_factor = ampl.getKinematics().getDegreesOfFreedomFactor();
    if (!csl::IsNumerical(comb_factor))
    {
      std::cerr << "The number of degrees of freedom is not numerical!\n";
      exit(3);
    }
    comb_factor = Evaluated(comb_factor, csl::eval::all);

    std::string mass_string = " ";
    for (size_t j = 0; j < 3; j++)
    {
      Expr tempmass = listofprocs_1to2[i].process[j].getField()->getMass();
      std::string tempstring = (tempmass != CSL_0) ? "static_cast<double>(input." + tempmass->getName() + ")" : "0.0";
      mass_string += tempstring + ",";
    }

    // Determining the symmetry factor
    Sf34 =
        (((listofprocs_1to2[i].process[1].isOutgoingParticle() &&
           listofprocs_1to2[i].process[2].isOutgoingParticle()) ||
          (listofprocs_1to2[i].process[1].isOutgoingAntiParticle() &&
           listofprocs_1to2[i].process[2].isOutgoingAntiParticle())) &&
         listofprocs_1to2[i].process[1].getField()->getName() == listofprocs_1to2[i].process[2].getField()->getName())
            ? 2
            : 1;

    // Determining the CP symmetry factor
    C1234 = 1;
  #ifdef CORRESPONDANCE
    // Writing on the correspondance file
    correspondancefileh << "{ {";
    for (int j = 0; j <= 1; j++)
    {
      if (!listofprocs_1to2[i].process[j].isParticle())
        correspondancefileh << "ANTICHAR,";
      correspondancefileh << "EMPTYCHAR+corr::" << listofprocs_1to2[i].process[j].getField()->getName() << ",";
    }
    if (!listofprocs_1to2[i].process[2].isParticle())
      correspondancefileh << "ANTICHAR,";
    correspondancefileh << "EMPTYCHAR+corr::" << listofprocs_1to2[i].process[2].getField()->getName() << "},";
    // Writing the value field
    correspondancefileh << "{ ";
    // Writing the functions
    correspondancefileh << " &" << nameSumSq << ", ";
    correspondancefileh << comb_factor << ",";
    // Writing the number
    correspondancefileh << C1234;
    correspondancefileh << "} }";
    if (i != lastproc_1to2 - 1)
      correspondancefileh << ",\n";
  #endif // End of ifdef CORRESPONDANCE
  } // Ends the cycle on the processes
  #ifdef CORRESPONDANCE
  correspondancefileh << "};\n}\n";
  correspondancefileh.close();
  #endif
#endif // ended #ifndef DISABLE_LIB_CREATION
  return 0;
}
