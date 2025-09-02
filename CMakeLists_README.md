# `CMakeLists.txt` Documentation

## Options

### `BUILD_DOCUMENTATION`

Option to enable or disable documentation generation with Doxygen.

### `ENABLE_SUBMODULE_CMA`

Option to enable or disable the compilation of the CMA submodule and dependent libraries.

### `ENABLE_CLANG_FORMAT`

Option to enable or disable clang-format to modify files.

## Variables

### `MODELS_LIBS_NAMES`

List of model libraries to compile.

### `EXTERNAL_LIBS`

List of external libraries to be linked.

### `COMPILE_LINK_FLAGS`

List of compile options for linking external libraries.

### `LISTS_COMMANDS_TO_INSTALL_HEADERS_ONLY`

List of commands to install headers only.

### `EXCLUDE_PATTERNS`

Patterns to exclude from clang-format.

## Targets

### `install_headers_all`

Custom target for headers installation.

### `darkpack_doc`

Custom target for documentation generation with Doxygen.

### `format`

Custom target to run clang-format to change files.

### `format-check`

Custom target to run clang-format to check code style compliance.

## Subdirectories

### `submodules/CMA-ES`

Adds the `CMA-ES` submodule if `ENABLE_SUBMODULE_CMA` is enabled.

## Documentation Generation

If `BUILD_DOCUMENTATION` is enabled and Doxygen is found, a custom target for documentation generation is added.

## Clang-Format

If `ENABLE_CLANG_FORMAT` is enabled, custom targets for running clang-format are added.

## Functions

### `print`

Prints the values of the given variables.

### `list_all_models`

Lists all models in the `src` directory that start with `dp_`.

### `should_exclude`

Checks if a file should be excluded based on the defined patterns.