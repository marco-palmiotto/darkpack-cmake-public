# Devcontainers

This folder contains the Dockerfiles for building
  1. The devcontainer images: these are the ones ending with `-dev`. They are the most complete. They contain everything needed for development:
      1.  The dependencies of `MARTY`
      2.  The dependencies of `CMA-ES`
      3.  Doxygen, clang-format
      4.  Some profiling tools
      5.  The installation of `MARTY` 
  2. The images for deployment and CI tests: these ends with `-deploy`. They are the most minimalistic ones, containing
      1.  The dependencies of `MARTY`
      2.  Doxygen, clang-format
      3.  The installation of `MARTY` 
  3. The images for launching a scan: these end with `-scan`. They are like the minimalistic ones, but with also the library `CMA-ES`. They contain:
      1.  The dependencies of `MARTY`
      2.  The dependencies of `CMA-ES`
      3.  Doxygen, clang-format
      4.  The installation of `MARTY` 