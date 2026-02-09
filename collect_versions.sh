#!/bin/bash
set -euo pipefail

pkgs=(gcc gcc-c++ gcc-gfortran gsl gsl-devel texlive ImageMagick dvipng coreutils qtchooser cmake git ninja-build gdb eigen3-devel tbb tbb-devel wget python3 jupyter-notebook doxygen git-clang-format pipx texlive-type1cm bc valgrind perf root pdftk)

for p in "${pkgs[@]}"; do
  echo -n "$p: "
  if command -v dnf >/dev/null 2>&1; then
    # prefer repoquery if available
    if command -v repoquery >/dev/null 2>&1; then
      repoquery --available --latest-limit=1 --qf '%{name}-%{version}-%{release}.%{arch}' "$p" || echo "(no available pkg)"
    else
      dnf info "$p" 2>/dev/null | awk -F: '/^Version/ {v=$2} /^Release/ {r=$2} END {gsub(/^ +| +$/,"",v); gsub(/^ +| +$/,"",r); if (v) print v"-"r; else print "(not found)"}'
    fi
  else
    echo "dnf not found"
  fi
done
