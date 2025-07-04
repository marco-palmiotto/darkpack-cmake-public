int computeAndAddToLibFromList(mty::Model& model, // model
                               mty::Library& lib, // output library
                               std::vector<Process2to2ToCompute> listofprocs,
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
  std::ofstream output_stream;

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
  output_stream.open(tempstring_c);

  if (!output_stream)
  {
    std::cerr << "Impossible to write the correspondance.hpp file\n";
    exit(1);
  }
  std::cout << "Creating correspondance.hpp file\n";
  output_stream << "#pragma once\n\n";
  output_stream << "#ifndef __cplusplus\n";
  output_stream << "#define __cplusplus\n";
  output_stream << "#endif\n\n";
  output_stream << "#include \"" << lib.getName() << ".h\"\n";
  output_stream << "#include \"clib_" << lib.getName() << ".h\"\n";
  output_stream << "#include <unordered_map>\n";
  output_stream << "namespace " << lib.getName() << "{\n";
  output_stream << "// We are using char 32 (space) as possible separator, char 126 (tilde) as antiParticle tag\n";
  output_stream << "#define EMPTYCHAR 32\n";
  output_stream << "#define ANTICHAR 126\n";
  output_stream << "using Cfptr_t = ccomplex_return_t (*) (cparam_t const *);\n";
  output_stream << "namespace corr{\n";
  output_stream << "  using Entry_t = std::tuple<Cfptr_t, short int, short int>;\n";
  output_stream << "  enum Part_t { \n";
  #ifdef DEBUG
  std::cout << "Writing particle names in enum Part_t\n";
  #endif
  output_stream << "        " << pnames[0] << "=1,\n";
  #ifdef DEBUG
  std::cout << "i = 0\n";
  #endif
  for (size_t i = 1; i != part.size() - 1; i++)
  {
  #ifdef DEBUG
    std::cout << "i = " << i << " Writing " << pnames[i] << std::endl;
  #endif
    output_stream << "        " << pnames[i] << ",\n";
  }
  output_stream << "       " << pnames.back() << std::endl;
  output_stream << "   };\n\n";

  output_stream << "#include \"" << nameSmBsmFile << "\"\n";
  #ifdef DEBUG
  std::cout << "Writing particle names in part_names\n";
  #endif
  output_stream << "inline const std::array<std::string,TOTAL_PARTICLES+1> part_names={\" \",\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    output_stream << "        \"" << pnames[i] << "\",\n";
  }
  output_stream << "       \"" << pnames.back() << '\"' << std::endl;
  output_stream << "   };\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle statistic\n";
  #endif

  output_stream << "inline const std::array<bool,TOTAL_PARTICLES+1> isboson={false, \n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    output_stream << "        " << (part[i]->isBosonic() ? "true" : "false") << ",\n";
  }
  output_stream << "       " << (part.back()->isBosonic() ? "true" : "false") << std::endl;
  output_stream << "   };\n\n";


  #ifdef DEBUG
  std::cout << "Writing particle codes\n";
  #endif
  output_stream << "inline const std::string mass_codes={EMPTYCHAR, \n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    output_stream << "        EMPTYCHAR+" << pnames[i] << ",\n";
  }
  output_stream << "       EMPTYCHAR+" << pnames.back() << std::endl;
  output_stream << "   };\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle charge\n";
  #endif
  output_stream << "inline const std::array<double,TOTAL_PARTICLES+1> part_charge={ 0,\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    output_stream << "        " << electric_charges.at(i) << ",\n";
  }
  output_stream << "        " << electric_charges.back() << "\n};\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle total degrees of freedom\n";
  #endif
  output_stream << "inline const std::array<int,TOTAL_PARTICLES+1> part_tot_dof={ 0,\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    output_stream << "        " << dofs.at(i) << ",\n";
  }
  output_stream << "        " << dofs.back() << "\n};\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle helicity degrees of freedom\n";
  #endif
  output_stream << "inline const std::array<int,TOTAL_PARTICLES+1> part_hel_dof={ 0,\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    output_stream << "        " << h_dofs.at(i) << ",\n";
  }
  output_stream << "        " << h_dofs.back() << "\n};\n\n";

  #ifdef DEBUG
  std::cout << "Writing particle conjugation\n";
  #endif
  output_stream << "inline const std::array<bool,TOTAL_PARTICLES+1> part_isSelfConj={ false,\n";
  for (size_t i = 0; i != part.size() - 1; i++)
  {
    output_stream << "        " << (isSelfConjugate.at(i) ? "true" : "false") << ",\n";
  }
  output_stream << "        " << (isSelfConjugate.back() ? "true" : "false") << "\n};\n\n";


  output_stream << "\nextern const std::unordered_map<std::string, Entry_t> squaredampl;\n";

  #ifdef DEBUG
  std::cout << "Writing getMassFirst\n";
  #endif
  output_stream << "\ninline double getMassFirst(const int p, const param_t &input)\n";
  output_stream << "{\n";

  output_stream << "    switch (p)\n";
  output_stream << "    {\n";

  for (size_t i = 0; i != part.size(); i++)
  {
    output_stream << "        case(" << pnames[i] << ") : \n";
    output_stream << "            return " << mnames[i] << ";\n";
  }
  output_stream << "        default : \n";
  output_stream << "            return -1.;\n";
  output_stream << "    }\n";
  output_stream << "}// End of getMassFirst\n";

  output_stream << "}//End of namespace corr\n";
  output_stream << "}//End of namespace " << lib.getName() << "\n";
  // correspondancefileh << "#endif\n";
  output_stream.close();

  std::cout << "correspondance.h has been created succesfully.\n";

  #ifndef SQUAREDAMP
    #ifdef UPDATE_PART_DETAILS
  return 0;
    #endif
  #endif

  #ifndef DISABLE_LIB_CREATION

  std::cout << "Creating initialise_map.cpp file\n";
  sprintf(tempstring_c, "auxiliary_library/%s/initialise_map.cpp", lib.getName().c_str());
  output_stream.open(tempstring_c);
  if (!output_stream)
  {
    std::cerr << "Impossible to write the " << tempstring_c << " file\n";
    return 1;
  }

  output_stream << "#include \"correspondance.hpp\"\n"
                << "namespace " << lib.getName() << "::corr{\n"
                << "const std::unordered_map<std::string, Entry_t> squaredampl{ \n";

  #endif // End of ifdef CORRESPONDANCE

  // Assigning the 1st values to the strings with the names
  // of the previous processes (starting with 0)
  prevname1 = listofprocs[0].process[0].getField()->getName();
  prevname2 = listofprocs[0].process[1].getField()->getName();

  const size_t lastproc_index = listofprocs.size();
  size_t i = 1;
  for (auto& single_proc : listofprocs)
  {
    count_converted++;
    const std::string procname = processName(single_proc.process);
    const std::string nameSumSq = "sumSqAmpl_" + procname;
    const std::string nameCombFac = "combFac_" + procname;
    std::cout << "Process n " << i << " / " << lastproc_index << "\n";

    // Adding the process in the library
    std::cout << "Computing " << nameSumSq << " ";
    // Setting W boson gauge
    model.getParticle("W")->setGaugeChoice(single_proc.Wgauge);

    // Creating a lambda function that returns the amplitude at the desired order
    auto lambda_ampl = [&]()
    {
      if (!single_proc.leading_order)
        return model.computeAmplitude(single_proc.order, single_proc.process);
      // If you want it at the leading order, let us compute the tree-level expression
      auto ampl_tree = model.computeAmplitude(mty::Order::TreeLevel, single_proc.process);
      if (!ampl_tree.empty())
        return ampl_tree;
      auto ampl_one_loop = model.computeAmplitude(mty::Order::OneLoop, single_proc.process);
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
    if (single_proc.process[0].getField()->getName() != prevname1 ||
        single_proc.process[1].getField()->getName() != prevname2)
    {
      // If one of the two particles is different, I save the file and I start a new one
      SaveDiagrams(lib.getName() + "_graphs/" + prevname1 + "_" + prevname2 + ".json", graphs);
      prevname1 = single_proc.process[0].getField()->getName();
      prevname2 = single_proc.process[1].getField()->getName();
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
      Expr tempmass = single_proc.process[j].getField()->getMass();
      std::string tempstring = (tempmass != CSL_0) ? "static_cast<double>(input." + tempmass->getName() + ")" : "0.0";
      mass_string += tempstring + ",";
    }

    // Determining the symmetry factor
    Sf34 = (((single_proc.process[2].isOutgoingParticle() && single_proc.process[3].isOutgoingParticle()) ||
             (single_proc.process[2].isOutgoingAntiParticle() && single_proc.process[3].isOutgoingAntiParticle())) &&
            single_proc.process[2].getField()->getName() == single_proc.process[3].getField()->getName())
               ? 2
               : 1;

    // Determining the CP symmetry factor
    C1234 = checkCPsymmetric(model, single_proc.process) ? 2 : 1;
  #ifdef CORRESPONDANCE
    // Writing on the correspondance file
    output_stream << "{ {";
    for (int j = 0; j <= 2; j++)
    {
      if (!single_proc.process[j].isParticle())
        output_stream << "ANTICHAR,";
      output_stream << "EMPTYCHAR+corr::" << single_proc.process[j].getField()->getName() << ",";
    }
    if (!single_proc.process[3].isParticle())
      output_stream << "ANTICHAR,";
    output_stream << "EMPTYCHAR+corr::" << single_proc.process[3].getField()->getName() << "},";
    // Writing the value field
    output_stream << "{ ";
    // Writing the functions
    output_stream << " &c_" << nameSumSq << ", ";
    output_stream << comb_factor << ",";
    // Writing the number
    output_stream << C1234;
    output_stream << "} }";
    if (i != lastproc_index - 1)
      output_stream << ",\n";
  #endif // End of ifdef CORRESPONDANCE
    i++;
  } // Ends the cycle on the processes
  #ifdef CORRESPONDANCE
  output_stream << "};\n}\n";
  output_stream.close();
  #endif

#endif // ended #ifndef DISABLE_LIB_CREATION
  return 0;
}
