int addFromFile(mty::Model& model,                       // Model
                Processes& processes,                    // Vector where I add the processes
                std::vector<std::string>& names,         // Name of the function in MARTY
                std::vector<std::string>& namesSuperiso, // Name of the process in SuperIso
                const std::string filename,              // Input file
                bool one_to_two = false                  // Flag for 1to2 processes (default false)
)
{

  std::vector<std::string> namesChep;
  char outstrings[MAXNUMPROCESSES][DIMCHEPSTRING];
  std::ofstream savechepnames;
  savechepnames.open(NAMECHEPFILE, std::ios::app);
  if (!savechepnames)
  {
    std::cout << "Impossible to open " << NAMECHEPFILE << "\n";
    exit(1);
  }

  std::ofstream savenamesclean;
  savenamesclean.open(NAMES_CLEAN, std::ios::out);
  if (!savenamesclean)
  {
    std::cout << "Impossible to open names_clean.dat\n";
    savechepnames.close();
    exit(1);
  }
  const std::string listSisoSusy[SIZEPYSICALSUSY] = {"o1", "o2", "o3", "o4", "c1", "c2", "t1", "t2", "b1", "b2",
                                                     "dr", "dl", "ur", "ul", "cr", "cl", "sr", "sl", "er", "el",
                                                     "mr", "ml", "l1", "l2", "ne", "nm", "nl", "go"};

  const std::string listSisoSM[SIZEPYSICALSM] = {"h", "hh", "hc", "h3", "w", "z", "g", "a",  "d",  "u",
                                                 "s", "c",  "b",  "t",  "e", "m", "l", "ne", "nm", "nl"};


  const std::string listChepSusy[SIZEPYSICALSUSY] = {
      "~o1", "~o2", "~o3", "~o4", "~1+", "~2+", "~t1", "~t2", "~b1", "~b2", "~dR", "~dL", "~uR", "~uL",
      "~cR", "~cL", "~sR", "~sL", "~eR", "~eL", "~mR", "~mL", "~l1", "~l2", "~ne", "~nm", "~nl", "~g"};

  const std::string listChepSM[SIZEPYSICALSM] = {"h", "H", "H+", "H3", "W+", "Z", "G", "A",  "d",  "u",
                                                 "s", "c", "b",  "t",  "e",  "m", "l", "ne", "nm", "nl"};

  const std::string listChepSusyAnti[SIZEPYSICALSUSY] = {
      "~o1", "~o2", "~o3", "~o4", "~1-", "~2-", "~T1", "~T2", "~B1", "~B2", "~DR", "~DL", "~UR", "~UL",
      "~CR", "~CL", "~SR", "~SL", "~ER", "~EL", "~MR", "~ML", "~L1", "~L2", "~Ne", "~Nm", "~Nl", "~g"};

  const std::string listChepSMAnti[SIZEPYSICALSM] = {"h", "H", "H-", "H3", "W-", "Z", "G", "A",  "D",  "U",
                                                     "S", "C", "B",  "T",  "E",  "M", "L", "Ne", "Nm", "Nl"};

  const std::string listMartySusy[SIZEPYSICALSUSY] = {
      "N_1",   "N_2",   "N_3",    "N_4",    "C_1",   "C_2",    "st_1",    "st_2", "sb_1", "sb_2",
      "sd_R",  "sd_L",  "su_R",   "su_L",   "sc_R",  "sc_L",   "ss_R",    "ss_L", "se_R", "se_L",
      "smu_R", "smu_L", "stau_1", "stau_2", "snu_e", "snu_mu", "snu_tau", "sG"};

  const std::string listMartySM[SIZEPYSICALSM] = {"h", "H0", "Hp", "A0", "W", "Z",  "G",   "A",    "d",     "u",
                                                  "s", "c",  "b",  "t",  "e", "mu", "tau", "nu_e", "nu_mu", "nu_tau"};

  FILE* source;
  source = fopen(filename.c_str(), "r");
  if (source == nullptr)
  {
    std::cout << "File " << filename << " cannot be opened for reading!" << std::endl;
    savechepnames.close();
    savenamesclean.close();
    return 1;
  }

  int nproc = 0;
  int particle = 0;
  int tempc = 0;
  std::string line, buffer, prevline = "";
  char tempstring[20];

  if (!one_to_two)
  {
    while (!feof(source))
    {
      fscanf(source, "%s", tempstring);
      line = tempstring;
      if (line.size() < 4 || line.size() > 20)
      {
        std::cout << "Line has invalid length" << std::endl;
        savechepnames.close();
        savenamesclean.close();
        fclose(source);
        return 0;
      }
      if (line == prevline)
        continue;
      prevline = line;
#ifdef DEBUG
      std::cout << "Reading line " << line << std::endl;
#endif
      particle = 0;
      size_t position = 0;
      int indexes[4];
      bool antiPart[4];
      size_t count = 0;
      bool perfectmatch = false;
      Process tempProcess;
      Particle tempParticle;

      buffer.clear();
      tempProcess.clear();
      // Reading the incoming state (SUSY)
      do
      {
        if (line.at(position) != '\0')
          buffer.push_back(line.at(position));
        else
          break;
        for (int i = 0; i < SIZEPYSICALSUSY; i++)
        {
          if (buffer == listSisoSusy[i])
          {
            indexes[particle] = i;
            if (line.size() >= position + 3 && line[position + 1] == 'b' && line[position + 2] == 'a' &&
                line[position + 3] == 'r')
            {
              antiPart[particle] = true;
              position += 3;
            }
            else
            {
              antiPart[particle] = false;
            }
            particle++;
            buffer.clear();
            i = SIZEPYSICALSUSY;
          }
        }
        position++;
      } while (particle < 2);

      // Reading the outgoing state (SM)
      do
      {
        if (line.at(position) != '\0')
          buffer.push_back(line.at(position));
        else
          break;
        for (int i = 0; i < SIZEPYSICALSM; i++)
        {
          if (buffer == listSisoSM[i])
          {
            // If particle is h*, checking a second character
            bool cond1 = (buffer.size() == 1 && buffer.at(0) == 'h' && line.size() > position + 1);
            if (cond1 && !(particle == 2 && line.size() == position + 1))
            {
              perfectmatch = false;
              if (particle == 2)
              {
                // Check if the rest of the string matches exactly with a particle
                std::string termination_buf;
                for (size_t newpos = position + 1; newpos < line.size(); newpos++)
                  termination_buf.push_back(line.at(newpos));
                for (int j = 0; j < SIZEPYSICALSM; j++)
                {
                  if (termination_buf == listSisoSM[j] || termination_buf == listSisoSM[j] + "bar")
                  {
                    // If there's a match
                    perfectmatch = true;
                    i = 0;
                    j = SIZEPYSICALSM;
                  }
                }
              }
              if (!perfectmatch)
                switch (line.at(position + 1))
                {
                case 'h':
                  i = 1;
                  position++;
                  break;
                case 'c':
                  i = 2;
                  position++;
                  break;
                case '3':
                  i = 3;
                  position++;
                  break;
                default:
                  break;
                }
            }
            indexes[particle] = i;
            if (line.size() >= position + 3 && line.at(position + 1) == 'b' && line.at(position + 2) == 'a' &&
                line.at(position + 3) == 'r')
            {
              antiPart[particle] = true;
              position += 3;
            }
            else
            {
              antiPart[particle] = false;
            }
            particle++;
            buffer.clear();
            i = SIZEPYSICALSM;
          }
        }
        position++;
      } while (particle < 4);

      count = 0;

      do
      {
        tempParticle = model.getParticle(listMartySusy[indexes[count]]);
        if (antiPart[count])
        {
          tempProcess.push_back(AntiPart(Incoming(tempParticle)));
        }
        else
        {
          tempProcess.push_back(Incoming(tempParticle));
        }
        count++;
      } while (count < 2);

      do
      {
        tempParticle = model.getParticle(listMartySM[indexes[count]]);
        if (antiPart[count])
        {
          tempProcess.push_back(AntiPart(Outgoing(tempParticle)));
        }
        else
        {
          tempProcess.push_back(Outgoing(tempParticle));
        }
        count++;
      } while (count < 4);
      sprintf(outstrings[nproc], "%s,%s->%s,%s",
              antiPart[0] ? listChepSusyAnti[indexes[0]].c_str() : listChepSusy[indexes[0]].c_str(),
              antiPart[1] ? listChepSusyAnti[indexes[1]].c_str() : listChepSusy[indexes[1]].c_str(),
              antiPart[2] ? listChepSMAnti[indexes[2]].c_str() : listChepSM[indexes[2]].c_str(),
              antiPart[3] ? listChepSMAnti[indexes[3]].c_str() : listChepSM[indexes[3]].c_str());

      namesChep.push_back(outstrings[nproc]);
      names.push_back(processName(tempProcess));
      processes.push_back(tempProcess);
      namesSuperiso.push_back(line);
      savechepnames << outstrings[nproc] << "\t" << line << std::endl;
      savenamesclean << processName(tempProcess) << " " << line << std::endl;
      nproc++;
    } // Ending cycle
    savechepnames.close();
    fclose(source);
  }
  else
  {
    while (!feof(source))
    {
      fscanf(source, "%s", tempstring);
      line = tempstring;
      if (line.size() < 3 || line.size() > 20)
      {
        std::cout << "Line has invalid length" << std::endl;
        savechepnames.close();
        savenamesclean.close();
        fclose(source);
        return 0;
      }
      if (line == prevline)
        continue;
      prevline = line;
#ifdef DEBUG
      std::cout << "Reading line " << line << std::endl;
#endif
      particle = 0;
      size_t position = 0;
      int indexes[3];
      bool antiPart[3];
      size_t count = 0;
      bool perfectmatch = false;
      Process tempProcess;
      Particle tempParticle;

      buffer.clear();
      tempProcess.clear();
      // Reading the incoming state (SM)
      do
      {
        if (line.at(position) != '\0')
          buffer.push_back(line.at(position));
        else
          break;
        for (int i = 0; i < SIZEPYSICALSM; i++)
        {
          if (buffer == listSisoSM[i])
          {
            // If particle is h*, checking a second character
            bool cond1 = (buffer.size() == 1 && buffer.at(0) == 'h' && line.size() > position + 1);
            if (cond1 && !(particle == 1 && line.size() == position + 1))
            {
              perfectmatch = false;
              if (particle == 1)
              {
                // Check if the rest of the string matches exactly with a particle
                std::string termination_buf;
                for (size_t newpos = position + 1; newpos < line.size(); newpos++)
                  termination_buf.push_back(line.at(newpos));
                for (int j = 0; j < SIZEPYSICALSM; j++)
                {
                  if (termination_buf == listSisoSM[j] || termination_buf == listSisoSM[j] + "bar")
                  {
                    // If there's a match
                    perfectmatch = true;
                    i = 0;
                    j = SIZEPYSICALSM;
                  }
                }
              }
              if (!perfectmatch)
                switch (line.at(position + 1))
                {
                case 'h':
                  i = 1;
                  position++;
                  break;
                case 'c':
                  if (line.at(position + 2) == 'b')
                    i = 0; // To distinguinsh between hc and h + cbar
                  else
                  {
                    i = 2;
                    position++;
                  }
                  break;
                case '3':
                  i = 3;
                  position++;
                  break;
                default:
                  break;
                }
            }
            indexes[particle] = i;
            if (line.size() >= position + 3 && line.at(position + 1) == 'b' && line.at(position + 2) == 'a' &&
                line.at(position + 3) == 'r')
            {
              antiPart[particle] = true;
              position += 3;
            }
            else
            {
              antiPart[particle] = false;
            }
            particle++;
            buffer.clear();
            i = SIZEPYSICALSM;
          }
        }
        position++;
      } while (particle < 1);

      // Reading the outgoing state (SM)
      do
      {
        if (line.at(position) != '\0')
          buffer.push_back(line.at(position));
        else
          break;
        for (int i = 0; i < SIZEPYSICALSM; i++)
        {
          if (buffer == listSisoSM[i])
          {
            // If particle is h*, checking a second character
            bool cond1 = (buffer.size() == 1 && buffer.at(0) == 'h' && line.size() > position + 1);
            if (cond1 && !(particle == 1 && line.size() == position + 1))
            {
              perfectmatch = false;
              if (particle == 1)
              {
                // Check if the rest of the string matches exactly with a particle
                std::string termination_buf;
                for (size_t newpos = position + 1; newpos < line.size(); newpos++)
                  termination_buf.push_back(line.at(newpos));
                for (int j = 0; j < SIZEPYSICALSM; j++)
                {
                  if (termination_buf == listSisoSM[j] || termination_buf == listSisoSM[j] + "bar")
                  {
                    // If there's a match
                    perfectmatch = true;
                    i = 0;
                    j = SIZEPYSICALSM;
                  }
                }
              }
              if (!perfectmatch)
                switch (line.at(position + 1))
                {
                case 'h':
                  i = 1;
                  position++;
                  break;
                case 'c':
                  i = 2;
                  position++;
                  break;
                case '3':
                  i = 3;
                  position++;
                  break;
                default:
                  break;
                }
            }
            indexes[particle] = i;
            if (line.size() >= position + 3 && line.at(position + 1) == 'b' && line.at(position + 2) == 'a' &&
                line.at(position + 3) == 'r')
            {
              antiPart[particle] = true;
              position += 3;
            }
            else
            {
              antiPart[particle] = false;
            }
            particle++;
            buffer.clear();
            i = SIZEPYSICALSM;
          }
        }
        position++;
      } while (particle < 3);


      count = 0;

      do
      {
        tempParticle = model.getParticle(listMartySM[indexes[count]]);
        if (antiPart[count])
        {
          tempProcess.push_back(AntiPart(Incoming(tempParticle)));
        }
        else
        {
          tempProcess.push_back(Incoming(tempParticle));
        }
        count++;
      } while (count < 1);

      do
      {
        tempParticle = model.getParticle(listMartySM[indexes[count]]);
        if (antiPart[count])
        {
          tempProcess.push_back(AntiPart(Outgoing(tempParticle)));
        }
        else
        {
          tempProcess.push_back(Outgoing(tempParticle));
        }
        count++;
      } while (count < 3);

      sprintf(outstrings[nproc], "%s->%s,%s",
              antiPart[0] ? listChepSusyAnti[indexes[0]].c_str() : listChepSusy[indexes[0]].c_str(),
              antiPart[1] ? listChepSusyAnti[indexes[1]].c_str() : listChepSusy[indexes[1]].c_str(),
              antiPart[2] ? listChepSMAnti[indexes[2]].c_str() : listChepSM[indexes[2]].c_str());

      namesChep.push_back(outstrings[nproc]);
      names.push_back(processName(tempProcess));
      processes.push_back(tempProcess);
      namesSuperiso.push_back(line);
      savechepnames << outstrings[nproc] << "\t" << line << std::endl;
      savenamesclean << processName(tempProcess) << " " << line << std::endl;
      nproc++;
    } // Ending cycle
    savechepnames.close();
    fclose(source);
  }
#ifdef DEBUG
  std::cout << "Read " << nproc << " processes\n";
#endif
  return nproc;
}