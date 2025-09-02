# Code policies and conventions

The main part of the source code has been written in more than 3 years, and the code policies were not defined yet. This document is an attempt of defining some conventions, and good practices to observe for contributing code. Ideally, the rest of the code will at some point respect these rules as well:

- Names that start with a capital letter, and/or end with `_t` are names of custom types. Examples: `relicparam_t`, `RunningSM`, `Param_t`.
  - As a consequence, name of variables, or functions (including methods), must begin with a lowercase letter.
- Avoid abbreviations in the names of variables, to make them self-explanatory.
  - The only exception to this might be the names of variables defined in the `MARTY` self-generated code, where, for instance, the mass names are defined as `m_<name of the particle>` instead of `mass_<name of the particle>`.
- All the names of variables composed by more than one word should separate the words with an `_`.
  - In some parts of the code , this convention is not respected, and I use several capital letters instead. At some point, I will change that. 
- While defining enumerations, do it in a namespace, or use an `enum class`.
- ... This document is a work in progress