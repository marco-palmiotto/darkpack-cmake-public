// This function just generates the file susysm.h
// It is not called in the main because you just need
// to call it once. The file susysm.h is moreover provided already.
void generateCorrespondance()
{
  const std::string listMartySusy[SIZEPYSICALSUSY] = {
      "N_1",   "N_2",   "N_3",    "N_4",    "C_1",   "C_2",    "st_1",    "st_2", "sb_1", "sb_2",
      "sd_R",  "sd_L",  "su_R",   "su_L",   "sc_R",  "sc_L",   "ss_R",    "ss_L", "se_R", "se_L",
      "smu_R", "smu_L", "stau_1", "stau_2", "snu_e", "snu_mu", "snu_tau", "sG"};

  const std::string listMartySM[SIZEPYSICALSM] = {"h", "H0", "Hp", "A0", "W", "Z",  "G",   "A",    "d",     "u",
                                                  "s", "c",  "b",  "t",  "e", "mu", "tau", "nu_e", "nu_mu", "nu_tau"};

  std::ofstream susysm{"auxiliary_library/susysm.h"};
  susysm << "inline const std::vector<int> susy_particles{\n";
  for (int i = 0; i < SIZEPYSICALSUSY - 1; i++)
  {
    susysm << "  " << listMartySusy[i] << ",\n";
  }
  susysm << "  " << listMartySusy[SIZEPYSICALSUSY - 1];
  susysm << "};\n";

  susysm << "inline const std::vector<int> sm_particles{\n";
  for (int i = 0; i < SIZEPYSICALSM - 1; i++)
  {
    susysm << "  " << listMartySM[i] << ",\n";
  }
  susysm << "  " << listMartySM[SIZEPYSICALSUSY - 1];
  susysm << "};\n";

  susysm << "inline const int SIZEPYSICALSM = 20;\n";
  susysm << "inline const int SIZEPYSICALSUSY = 28;\n";
  susysm << "inline const int TOTAL_PARTICLES = 48;\n";
}
