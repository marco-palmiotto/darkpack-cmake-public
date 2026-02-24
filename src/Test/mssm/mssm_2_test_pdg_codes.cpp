#include "dp_mssm2to2/pdg_mass_codes.hpp"
#include <cassert>
#include <iomanip>
#include <iostream>

using namespace mssm2to2::corr;

int main()
{
  std::cout << "Testing PDG code functions...\n" << std::endl;

  // Test 1: Get PDG code by Part_t enumeration for all particles
  std::cout << "=== Test 1: get_PDG_code(Part_t) for all particles ===" << std::endl;
  int pass_count = 0;
  int fail_count = 0;

  for (int i = 1; i <= TOTAL_PARTICLES; ++i)
  {
    Part_t particle = static_cast<Part_t>(i);
    int pdg = get_PDG_code(particle);
    std::string particle_name = part_names[i];

    std::cout << std::left << std::setw(3) << i << std::setw(20) << particle_name << std::setw(10)
              << "PDG: " << std::setw(10) << pdg;

    // Verify that the name lookup gives the same PDG code
    try
    {
      int pdg_from_name = get_PDG_code(particle_name);
      if (pdg == pdg_from_name)
      {
        std::cout << "[PASS]" << std::endl;
        pass_count++;
      }
      else
      {
        std::cout << "[FAIL] Name lookup returned: " << pdg_from_name << std::endl;
        fail_count++;
      }
    }
    catch (const std::out_of_range&)
    {
      std::cout << "[FAIL] Name lookup threw exception" << std::endl;
      fail_count++;
    }
  }

  std::cout << "\nTest 1 Results: " << pass_count << " passed, " << fail_count << " failed" << std::endl;

  if (fail_count > 0)
  {
    return 1;
  }

  std::cout << std::endl;

  // Test 2: Test invalid particle name (should throw exception)
  std::cout << "=== Test 2: Invalid particle name (exception handling) ===" << std::endl;
  try
  {
    get_PDG_code("invalid_particle");
    std::cout << "[FAIL] Should have thrown std::out_of_range" << std::endl;
    return 1;
  }
  catch (const std::out_of_range&)
  {
    std::cout << "[PASS] Correctly threw std::out_of_range for invalid particle name" << std::endl;
  }

  std::cout << std::endl;

  // Test 3: Test get_PDG_code_backend(Part_t)
  std::cout << "=== Test 3: get_PDG_code_backend(Part_t) for all particles ===" << std::endl;
  pass_count = 0;
  fail_count = 0;

  for (int i = 1; i <= TOTAL_PARTICLES; ++i)
  {
    Part_t particle = static_cast<Part_t>(i);
    const int pdg_backend = get_PDG_code_backend(particle);
    const int pdg_regular = get_PDG_code(particle);

    if (pdg_backend == pdg_regular)
    {
      std::cout << "[PASS] " << part_names[i] << " (PDG: " << pdg_backend << ")" << std::endl;
      pass_count++;
    }
    else
    {
      std::cout << "[FAIL] " << part_names[i] << " - backend returned " << pdg_backend << ", regular returned "
                << pdg_regular << std::endl;
      fail_count++;
    }
  }

  std::cout << "\nTest 3 Results: " << pass_count << " passed, " << fail_count << " failed" << std::endl;

  if (fail_count > 0)
  {
    return 1;
  }

  std::cout << std::endl;

  // Test 4: Test get_particle_name_from_PDG_code(int)
  std::cout << "=== Test 4: get_particle_name_from_PDG_code(int) ===" << std::endl;
  pass_count = 0;
  fail_count = 0;

  for (int i = 1; i <= TOTAL_PARTICLES; ++i)
  {
    Part_t particle = static_cast<Part_t>(i);
    int pdg = get_PDG_code(particle);
    std::string expected_name = part_names[i];

    try
    {
      std::string retrieved_name = get_particle_name_from_PDG_code(pdg);
      if (retrieved_name == expected_name)
      {
        std::cout << "[PASS] PDG " << std::setw(8) << pdg << " -> " << retrieved_name << std::endl;
        pass_count++;
      }
      else
      {
        std::cout << "[FAIL] PDG " << pdg << " - expected '" << expected_name << "', got '" << retrieved_name << "'"
                  << std::endl;
        fail_count++;
      }
    }
    catch (const std::out_of_range&)
    {
      std::cout << "[FAIL] PDG " << pdg << " threw std::out_of_range" << std::endl;
      fail_count++;
    }
  }

  std::cout << "\nTest 4 Results: " << pass_count << " passed, " << fail_count << " failed" << std::endl;

  if (fail_count > 0)
  {
    return 1;
  }

  std::cout << std::endl;

  // Test 5: Test get_dp_code_from_PDG_code(int)
  std::cout << "=== Test 5: get_dp_code_from_PDG_code(int) ===" << std::endl;
  pass_count = 0;
  fail_count = 0;

  for (int i = 1; i <= TOTAL_PARTICLES; ++i)
  {
    Part_t particle = static_cast<Part_t>(i);
    int pdg = get_PDG_code(particle);

    try
    {
      Part_t retrieved_particle = get_dp_code_from_PDG_code(pdg);
      if (retrieved_particle == particle)
      {
        std::cout << "[PASS] PDG " << std::setw(8) << pdg << " -> " << part_names[i] << std::endl;
        pass_count++;
      }
      else
      {
        std::cout << "[FAIL] PDG " << pdg << " - expected Part_t " << static_cast<int>(particle) << ", got "
                  << static_cast<int>(retrieved_particle) << std::endl;
        fail_count++;
      }
    }
    catch (const std::out_of_range&)
    {
      std::cout << "[FAIL] PDG " << pdg << " threw std::out_of_range" << std::endl;
      fail_count++;
    }
  }

  std::cout << "\nTest 5 Results: " << pass_count << " passed, " << fail_count << " failed" << std::endl;

  if (fail_count > 0)
  {
    return 1;
  }

  std::cout << std::endl;

  // Test 6: Test antiparticle handling (negative PDG codes)
  std::cout << "=== Test 6: Antiparticle handling (negative PDG codes) ===" << std::endl;
  pass_count = 0;
  fail_count = 0;

  // Test a few particles with antiparticle codes
  std::vector<std::string> test_particles = {"e", "mu", "u", "d"};

  for (const auto& particle_name : test_particles)
  {
    try
    {
      int pdg = get_PDG_code(particle_name);
      int antiparticle_pdg = -pdg;

      std::string antiparticle_name = get_particle_name_from_PDG_code(antiparticle_pdg);
      Part_t antiparticle_code = get_dp_code_from_PDG_code(antiparticle_pdg);
      Part_t particle_code = get_dp_code_from_PDG_code(pdg);

      // Antiparticles should have same Part_t code as particles
      if (antiparticle_code == particle_code)
      {
        std::cout << "[PASS] " << particle_name << " (PDG: " << pdg << ") <-> " << antiparticle_name
                  << " (PDG: " << antiparticle_pdg << ")" << std::endl;
        pass_count++;
      }
      else
      {
        std::cout << "[FAIL] Antiparticle mismatch for " << particle_name << std::endl;
        fail_count++;
      }
    }
    catch (const std::out_of_range&)
    {
      std::cout << "[FAIL] Exception thrown for particle " << particle_name << std::endl;
      fail_count++;
    }
  }

  std::cout << "\nTest 6 Results: " << pass_count << " passed, " << fail_count << " failed" << std::endl;

  if (fail_count > 0)
  {
    return 1;
  }

  std::cout << std::endl;
  std::cout << "\n=== All tests passed! ===" << std::endl;
  return 0;
}
