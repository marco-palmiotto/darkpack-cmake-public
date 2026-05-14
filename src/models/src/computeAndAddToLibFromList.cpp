#define PARTICLEDECAYS



csl::Expr squaredAmplitudeToPartialWidth_new(csl::Expr const& squaredAmplitude,
                                             std::vector<Insertion> const& insertions,
                                             mty::Amplitude const& initialAmplitude, bool applyMassCondition)
{
  csl::Expr squared = csl::DeepCopy(squaredAmplitude);
  auto pi = initialAmplitude.getKinematics().getOrderedMomenta();
  csl::Index mu = MinkowskiIndex();
  csl::Expr s12_init = pi[0](mu) * pi[1](+mu);
  csl::Expr s13_init = pi[0](mu) * pi[2](+mu);
  csl::Expr s23_init = pi[1](mu) * pi[2](+mu);
  csl::Expr M = insertions[0].getField()->getMass();
  csl::Expr m1 = insertions[1].getField()->getMass();
  csl::Expr m2 = insertions[2].getField()->getMass();
  csl::Expr E1 = (M * M + m1 * m1 - m2 * m2) / (2 * M);
  csl::Expr E2 = (M * M + m2 * m2 - m1 * m1) / (2 * M);
  csl::Expr p2 = E1 * E1 - m1 * m1;
  csl::Expr s12 = csl::Abbrev::makeAbbreviation(M * E1);
  csl::Expr s13 = csl::Abbrev::makeAbbreviation(M * E2);
  csl::Expr s23 = csl::Abbrev::makeAbbreviation(E1 * E2 + p2);
  csl::Replace(squared, s12_init, s12);
  csl::Replace(squared, s13_init, s13);
  csl::Replace(squared, s23_init, s23);
  squared *= csl::sqrt_s(csl::abs_s(p2)) / (8 * CSL_PI * M * M);

  std::string name_1 = insertions[1].getField()->getName();
  std::string name_2 = insertions[2].getField()->getName();


  csl::Expr squared_approx = CSL_0;

  if ((name_1 == "Z" && name_2 == "Z") || (name_1 == "W" && name_2 == "W"))
  {
    squared_approx = csl::DeepCopy(squared);
    csl::Expr x = (m1 * m2) / (M * M);
    csl::Expr term_1 = (1 - 8 * x + 20 * csl::pow_s(x, 2)) / (csl::sqrt_s(4 * x - 1));
    csl::Expr term_2 = csl::acos_s((3 * x - 1) / (2 * csl::pow_s(x, 1.5)));
    csl::Expr term_3 = (1 - x) / (2 * x) * (2 - 13 * x + 47 * csl::pow_s(x, 2));
    csl::Expr term_4 = (1 - 6 * x + 4 * csl::pow_s(x, 2)) * csl::log_s(x);
    csl::Expr F = 3 * term_1 * term_2 - term_3 - 1.5 * term_4;

    squared_approx *= F;
  }
  if (applyMassCondition)
  {
    csl::Expr final =
        csl::booleanOperator_s(csl::BooleanOperator::GreaterThanOrEqualTo, M, m1 + m2, squared, squared_approx);
    return final;
  }

  return squared;
}



static std::vector<std::vector<mty::Insertion>> getIndependentDecays(std::vector<mty::Insertion> const& insertions)
{
  HEPAssert(insertions.size() == 3, mty::error::TypeError,
            "This function should only be used for 1->2 processes, " + std::to_string(insertions.size()) +
                " insertions given.") std::vector<std::vector<mty::Insertion>>
      res;
  res.reserve(4);
  res.push_back(insertions);
  if (!insertions[1].getField()->isSelfConjugate())
  {
    res.push_back({insertions[0], AntiPart(insertions[1]), insertions[2]});
    if (!insertions[2].getField()->isSelfConjugate())
    {
      if (insertions[1].getField() != insertions[2].getField())
        res.push_back({insertions[0], insertions[1], AntiPart(insertions[2])});
      res.push_back({insertions[0], AntiPart(insertions[1]), AntiPart(insertions[2])});
    }
  }
  else if (!insertions[2].getField()->isSelfConjugate())
  {
    res.push_back({insertions[0], insertions[1], AntiPart(insertions[2])});
  }
  return res;
}

csl::Expr computeWidth_new(mty::Model& model, mty::Insertion const& particle, std::vector<Particle> psm,
                           std::vector<Particle> pbsm)
{
  int order = mty::Order::TreeLevel;
  if (particle.getField()->getMass() == CSL_0)
  {
    // No width for massless particles
    return CSL_0;
  }
  csl::ScopedProperty silent(&mty::option::verboseAmplitude, false);
  std::vector<csl::Expr> contributions;
  contributions.reserve(100);
  const auto& physicalParticles = model.getPhysicalParticles(
      [](Particle const& p) { return !IsOfType<GoldstoneBoson>(p) && !IsOfType<GhostBoson>(p); });
  for (size_t i = 0; i != physicalParticles.size(); ++i)
  {
    const auto& p1 = physicalParticles[i];
    for (size_t j = i; j != physicalParticles.size(); ++j)
    {
      const auto& p2 = physicalParticles[j];

      auto insertions = getIndependentDecays({Incoming(particle), Outgoing(p1), Outgoing(p2)});
      for (const auto& ins : insertions)
      {
        order = mty::Order::TreeLevel;
        mty::Amplitude ampl = model.computeAmplitude(order, ins);

        if (ampl.empty())
        {
          /*
          auto it_1 = std::find(pbsm.begin(), pbsm.end(), p1);
          auto it_2 = std::find(pbsm.begin(), pbsm.end(), p2);
          if ((ins[1].getField()->getName() == "A" || ins[2].getField()->getName() == "A" ||
               ins[1].getField()->getName() == "G" || ins[2].getField()->getName() == "G") &&
              (it_1 != pbsm.end() || it_2 != pbsm.end()))
          {
            order = mty::Order::OneLoop;
            std::cout << "Tree-level amplitude is null for decay " << particle.getField()->getName() << " -> "
                      << p1->getName() << " + " << p2->getName() << ". Let's compute the one-loop amplitude.\n";
            ampl = model.computeAmplitude(order, ins);
          }
            */
          std::cout << "Tree-level amplitude is null for decay " << particle.getField()->getName() << " -> "
                    << p1->getName() << " + " << p2->getName() << ". Let's go to the next decay.\n";
          continue;
        }
        csl::Expr squared;
        squared = model.computeSquaredAmplitude(ampl);

        if (squared != 0)
        {
          std::string decaying_part = ins[0].getField()->getName();
          std::string prod_1 = ins[1].getField()->getName();
          std::string prod_2 = ins[2].getField()->getName();
          std::cout << "Found decay ";
          std::cout << decaying_part << " -> ";
          std::cout << prod_1 << " + ";
          std::cout << prod_2 << std::endl;
          contributions.push_back(squaredAmplitudeToPartialWidth_new(squared, ins, ampl, true));
        }
      }
    }
  }
  std::cout << contributions.size() << " independent " << particle.getField()->getName() << " decays found.\n";
  return csl::sum_s(contributions);
}



void getArguments(csl::Expr a, mty::Model& model, std::vector<Particle> psm, std::vector<Particle> pbsm,
                  mty::Library& lib, std::vector<std::string>& process_names, std::vector<Process>& processes,
                  std::vector<csl::Expr>& expressions)
{
  for (size_t i = 0; i < Size(a); i++)
  {
    std::cout << i << std::endl;
    std::string procname;
    Process proc;
    if (Size(a[i]) < 2)
      continue;
    std::string test_decay = toString(a[i][0]);
    if (test_decay.size() < 2)
    {
      getArguments(a[i], model, psm, pbsm, lib, process_names, processes, expressions);
      continue;
    }
    std::vector<std::string> terms;
    std::string test_prods = toString(a[i][1]);
    test_decay.erase(0, 2);
    mty::Insertion ins_decay(test_decay);
    if (test_prods[0] == '2')
    {
      test_prods.erase(0, 4);
      terms.push_back(test_prods);
      terms.push_back(test_prods);
      mty::Insertion ins_prod(test_prods);

      if (!ins_prod.getField()->isSelfConjugate())
      {
        proc = {Incoming(model.getParticle(test_decay)), Outgoing(AntiPart(test_prods)), Outgoing(test_prods)};
      }
      else
      {
        proc = {Incoming(model.getParticle(test_decay)), Outgoing(test_prods), Outgoing(test_prods)};
      }
      procname = "partWidth_" + processName(proc);

      if (process_names.size() > 0)
      {
        auto it_name = std::find(process_names.begin(), process_names.end(), procname);
        if (it_name != process_names.end())
        {
          std::cout << "bingas " << procname << std::endl;
          continue;
        }
      }
    }
    else
    {
      std::cout << "case 2" << std::endl;
      std::stringstream ss(test_prods);
      std::string term;
      while (std::getline(ss, term, '+'))
      {
        term.erase(remove(term.begin(), term.end(), ' '), term.end());
        term.erase(0, 2);
        terms.push_back(term);
      }
      if (terms.size() < 2)
        continue;
      mty::Insertion ins_term0(terms[0]), ins_term1(terms[1]);
      for (size_t i = 0; i < terms.size(); i++)
      {
        std::cout << terms[i] << std::endl;
      }
      proc = {Incoming(ins_decay), Outgoing(ins_term0), Outgoing(ins_term1)};
      procname = "partWidth_" + processName(proc);
      if (process_names.size() > 0)
      {
        auto it_name = std::find(process_names.begin(), process_names.end(), procname);
        if (it_name != process_names.end())
        {
          std::cout << "bingas " << procname << std::endl;
          continue;
        }
      }
    }

    std::cout << test_decay << " -> " << terms[0] << " + " << terms[1] << std::endl;
    std::cout << procname << std::endl;


    auto decay_part = model.getParticle(test_decay);
    auto part_1 = model.getParticle(terms[0]);
    auto part_2 = model.getParticle(terms[1]);
    auto it_1 = std::find(pbsm.begin(), pbsm.end(), part_1);
    auto it_2 = std::find(pbsm.begin(), pbsm.end(), part_2);
    if (it_1 != pbsm.end() || it_2 != pbsm.end() || decay_part == part_1 || decay_part == part_2)
      continue;
    else
    {

      auto temp_a = csl::Evaluated(a[i], csl::eval::abbreviation);
      csl::ForEachLeaf(temp_a,
                       [&](csl::Expr& sub)
                       {
                         if (csl::IsConstant(sub) || csl::IsVariable(sub))
                         {
                           sub->setValue(CSL_UNDEF);
                         }
                       });
      expressions.push_back(a[i]);
      process_names.push_back(procname);
      processes.push_back(proc);
    }
  }
}

void printArgument(csl::Expr a)
{
  for (size_t i = 0; i != Size(a); ++i)
  {
    std::cout << i << " : " << a[i] << std::endl;
  }
  // Also possible to call a->size()
}

int computeAndAddToLibFromList(mty::Model& model, // model
                               mty::Library& lib, // output library
                               std::vector<Process2to2ToCompute> listofprocs, Process decaying_parts,
                               Process decay_prods, std::vector<Particle> psm, std::vector<Particle> pbsm,
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
  std::ofstream correspondancefileh, widthfile;

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
  correspondancefileh << "using CXXfptr_t = complex_t (*) (const param_t&);\n";
  correspondancefileh << "namespace corr{\n";
  correspondancefileh << "  using Entry_t = std::tuple<Cfptr_t, short int, short int>;\n";
  correspondancefileh << "  using Entry_t_1to2 = std::tuple<CXXfptr_t>;\n";
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
  correspondancefileh << "\nextern const std::unordered_map<std::string, Entry_t> widths_map;\n";

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

#endif // end of #ifdef CORRESPONDANCE
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
    std::cout << nameSumSq << "\n";
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

      std::cout << "bing\n";
      auto ampl_tree = model.computeAmplitude(mty::Order::TreeLevel, listofprocs[i].process);
      if (!ampl_tree.empty())
        return ampl_tree;
      std::cout << "Computing the amplitude at One Loop order since the tree-level amplitude is null\n";
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



  // Assigning the 1st values to the strings with the names
  // of the previous processes (starting with 0)~
  model.getParticle("W")->setGaugeChoice(listofprocs[0].Wgauge);
  /*
  for (auto& particle : decaying_parts)
  {
    for (size_t i = 0; i != decay_prods.size(); ++i)
    {
      const auto& p1 = decay_prods[i];
      for (size_t j = 1; j != decay_prods.size(); ++j)
      {
        const auto& p2 = decay_prods[j];
        auto insertions = getIndependentDecays({Incoming(particle), Outgoing(p1), Outgoing(p2)});
        for (const auto& ins : insertions)
        {
          mty::Amplitude ampl = model.computeAmplitude(mty::Order::TreeLevel, ins);
          if (ampl.empty())
          {
            // ampl = computeAmplitude(mty::Order::OneLoop, ins);
            continue;
          }
          csl::Expr squared, partial;
          if (!ampl.empty())
            squared = model.computeSquaredAmplitude(ampl);
          else
            continue;
          if (squared != 0)
          {
            std::cout << "Found decay ";
            std::cout << ins[0].getField()->getName() << " -> ";
            std::cout << ins[1].getField()->getName() << " + ";
            std::cout << ins[2].getField()->getName() << std::endl;
            partial = model.squaredAmplitudeToPartialWidth(squared, ins, ampl);
          }
          std::string procname = processName(ins);
          std::string nameSumSq = "sumSqAmpl_" + procname;
          std::string namePartWidth = "partWidth_" + procname;
          lib.addFunction(nameSumSq, squared);
          lib.addFunction(namePartWidth, partial);
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
            Expr tempmass = ins[j].getField()->getMass();
            std::string tempstring =
                (tempmass != CSL_0) ? "static_cast<double>(input." + tempmass->getName() + ")" : "0.0";
            mass_string += tempstring + ",";
          }
          // Determining the symmetry factor
          Sf34 = (((ins[1].isOutgoingParticle() && ins[2].isOutgoingParticle()) ||
                   (ins[1].isOutgoingAntiParticle() && ins[2].isOutgoingAntiParticle())) &&
                  ins[1].getField()->getName() == ins[2].getField()->getName())
                     ? 2
                     : 1;

          // Determining the CP symmetry factor
          C1234 = 1;
  #ifdef CORRESPONDANCE
          // Writing on the correspondance file
          correspondancefileh << "{ {";
          for (int j = 0; j <= 1; j++)
          {
            if (!ins[j].isParticle())
              correspondancefileh << "ANTICHAR,";
            correspondancefileh << "EMPTYCHAR+corr::" << ins[j].getField()->getName() << ",";
          }
          if (!ins[2].isParticle())
            correspondancefileh << "ANTICHAR,";
          correspondancefileh << "EMPTYCHAR+corr::" << ins[2].getField()->getName() << "},";
          // Writing the value field
          correspondancefileh << "{ ";
          // Writing the functions
          correspondancefileh << " &" << nameSumSq << ", ";
          correspondancefileh << " &" << namePartWidth << ", ";
          correspondancefileh << comb_factor << ",";
          // Writing the number
          correspondancefileh << C1234;
          correspondancefileh << "} }";
          if (particle != decaying_parts.back() || i != decay_prods.size() - 1 || j != decay_prods.size() - 1)
            correspondancefileh << ",\n";
  #endif // End of ifdef CORRESPONDANCE
        } // Ends the cycle on the processes
  #ifdef CORRESPONDANCE
        correspondancefileh << "};\n}\n";
        correspondancefileh.close();
  #endif
      }
    }
  }
  */

  sprintf(tempstring_c, "auxiliary_library/%s/widths_map.cpp", lib.getName().c_str());
  widthfile.open(tempstring_c);
  widthfile << "#include \"correspondance.hpp\"\n"
            << "namespace " << lib.getName() << "::corr{\n"
            << "const std::unordered_map<std::string, Entry_t> widths_map{ \n";

  sprintf(tempstring_c, "auxiliary_library/%s/initialise_map_1to2.cpp", lib.getName().c_str());
  correspondancefileh.open(tempstring_c);
  correspondancefileh << "#include \"correspondance.hpp\"\n"
                      << "namespace " << lib.getName() << "::corr{\n"
                      << "const std::unordered_map<std::string, Entry_t_1to2> squaredampl_1to2{ \n";

  for (auto& particle : decaying_parts)
  {

    std::string part_name = particle.getField()->getName();
    std::string width_name = "width_" + part_name;
    csl::Expr width_part = computeWidth_new(model, part_name, psm, pbsm);
    csl::Expr width_true = model.computeWidth(Order::TreeLevel, part_name);

    printArgument(width_part);

    std::vector<std::string> process_names;
    std::vector<Process> processes;
    std::vector<csl::Expr> expressions;



    auto temp = csl::Evaluated(width_part, csl::eval::abbreviation);
    csl::ForEachLeaf(temp,
                     [&](csl::Expr& sub)
                     {
                       if (csl::IsConstant(sub) || csl::IsVariable(sub))
                       {
                         sub->setValue(CSL_UNDEF);
                       }
                     });
    auto temp_true = csl::Evaluated(width_true, csl::eval::abbreviation);
    csl::ForEachLeaf(temp_true,
                     [&](csl::Expr& sub)
                     {
                       if (csl::IsConstant(sub) || csl::IsVariable(sub))
                       {
                         sub->setValue(CSL_UNDEF);
                       }
                     });
    getArguments(width_part, model, psm, pbsm, lib, process_names, processes, expressions);
    lib.addFunction(width_name, width_true);

    widthfile << "{ {";
    widthfile << "EMPTYCHAR+corr::" << particle.getField()->getName() << ",";
    widthfile << "},";
    widthfile << "{ ";
    // Writing the functions
    widthfile << " &c_" << width_name << ", 1, 1";
    // Writing the number
    widthfile << "} }";

    if (particle != decaying_parts.back())
      widthfile << ",\n";

    for (int i = 0; i < process_names.size(); i++)
    {
      std::cout << process_names[i];
      lib.addFunction(process_names[i], expressions[i]);
  #ifdef CORRESPONDANCE
      // Writing on the correspondance file
      correspondancefileh << "{ {";
      for (int j = 0; j <= 1; j++)
      {
        if (!processes[i][j].isParticle())
          correspondancefileh << "ANTICHAR,";
        correspondancefileh << "EMPTYCHAR+corr::" << processes[i][j].getField()->getName() << ",";
      }
      if (!processes[i][2].isParticle())
        correspondancefileh << "ANTICHAR,";
      correspondancefileh << "EMPTYCHAR+corr::" << processes[i][2].getField()->getName() << "},";
      // Writing the value field
      correspondancefileh << "{ ";
      // Writing the functions
      correspondancefileh << " &" << process_names[i] << ", ";
      // Writing the number
      correspondancefileh << "} }";
      correspondancefileh << ",\n";
  #endif // End of ifdef CORRESPONDANCE
    } // Ends the cycle on the processes
  } // Ends the cycle on the particles
  widthfile << "};\n}\n";
  widthfile.close();
  #ifndef PARTICLEDECAYS
  prevname1 = listofprocs_1to2[0].process[0].getField()->getName();
  prevname2 = listofprocs_1to2[0].process[1].getField()->getName();

  const size_t lastproc_1to2 = listofprocs_1to2.size();
  for (size_t i = 0; i != lastproc_1to2; i++)
  {
    count_converted++;
    std::string procname = processName(listofprocs_1to2[i].process);
    std::string nameSumSq = "sumSqAmpl_" + procname;
    std::string namePartWidth = "partWidth_" + procname;
    std::cout << nameSumSq << "\n";
    std::string nameCombFac = "combFac_" + procname;
    std::cout << "Process n " << i + 1 << " / " << lastproc_1to2 << "\n";

    // Adding the process in the library
    std::cout << "Setting the process in the library\n";
    std::cout << "Computing " << nameSumSq << " ";
    // Setting W boson gauge
    model.getParticle("W")->setGaugeChoice(listofprocs_1to2[i].Wgauge);
    std::cout << "W gauge set\n";
    // Creating a lambda function that returns the amplitude at the desired order
    auto lambda_ampl = [&]()
    {
      std::cout << "bing\n";
      if (!listofprocs_1to2[i].leading_order)
      {
        return model.computeAmplitude(listofprocs_1to2[i].order, listofprocs_1to2[i].process);
        std::cout << "Amplitude computed at order " << listofprocs_1to2[i].order << "\n";
      }
      // If you want it at the leading order, let us compute the tree-level expression
      auto ampl_tree = model.computeAmplitude(mty::Order::TreeLevel, listofprocs_1to2[i].process);
      if (!ampl_tree.empty())
        return ampl_tree;
      std::cout << "Computing the amplitude at One Loop order since the tree-level amplitude is null\n";
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
    Expr s1, p1;

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
    p1 = model.squaredAmplitudeToPartialWidth(s1, listofprocs_1to2[i].process, ampl, true);
    if (s1 == CSL_0)
    {
      std::cout << "The squared amplitude is null\n";
      continue;
    }

    // You need to evaluate all the abbreviation to have the full set of parameters
    // as input data
    auto temp = csl::Evaluated(s1, csl::eval::abbreviation);
    auto temp_p1 = csl::Evaluated(p1, csl::eval::abbreviation);
      #ifdef FULL_PARAM
    csl::ForEachLeaf(temp,
                     [&](csl::Expr& sub)
                     {
                       if (csl::IsConstant(sub) || csl::IsVariable(sub))
                       {
                         sub->setValue(CSL_UNDEF);
                       }
                     });

    csl::ForEachLeaf(temp_p1,
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
    lib.addFunction(namePartWidth, p1);
      #ifdef DISPLAYSQUAREDAMPLITUDESEXPRESSIONS
    std::cout << "|M|^2 = " << DeepRefreshed(temp) << std::endl;
      #endif
      #ifdef DISPLAYPARTIALWIDTHEXPRESSIONS
    std::cout << "Partial width " << namePartWidth << " = " << DeepRefreshed(temp_p1) << std::endl;
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
    correspondancefileh << " &" << namePartWidth << ", ";
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
  #endif // End of ifndef PARTICLEDECAYS
#endif   // ended #ifndef DISABLE_LIB_CREATION
  return 0;
}
