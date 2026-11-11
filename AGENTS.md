# Repository Guidelines

## Project Structure & Module Organization
`src/` contains the main codebase: C/C++ solver code in `src/c/`, MATLAB and Python APIs in `src/m/`, and language wrappers in `src/wrappers/`. Regression and integration tests live in `test/NightlyRun/` with language-specific files such as `test103.m`, `test252.py`, and `test229.js`. Use `examples/` for tutorial-scale workflows, `jenkins/` for CI configuration and reference build flags, `packagers/` for release packaging, and avoid editing `externalpackages/` unless you are intentionally updating vendored dependencies.

## Build, Test, and Development Commands
ISSM uses an Autotools/Automake build. Simply run `make -j 8 install` in the project root directory to compile the programs needed.
The tests should be written through the MATLAB interface, examples of which can be found in `examples` and `test`.

## Coding Style & Naming Conventions
Follow the style already present in the touched directory; there is no repo-wide formatter config checked in. C/C++ code generally uses tabs for indentation, brace-on-function style, and `CamelCase` for core types with lowercase `.cpp` filenames. Python uses 4-space indentation, snake_case methods, and class names like `sealevelmodel`. MATLAB and JavaScript tests keep compact procedural scripts. Preserve existing test numbering and naming: `test####.(m|py|js)`.

## Testing Guidelines
No need to carry out the unit tests provided in the repository. Instead only test the changes you make without worrying about the older functionalities

## Commit & Pull Request Guidelines
No need to handle git commit or pull requests. 

## Repository Specific Instructions 
This is a codebase for a finite element code used for numerical ice sheet modeling. All code written should reflect numerical rigor and also a degree of readability and interpretability for anyone who might want to examine the code. Provide ample documentation always to explain what the code is doing. Respect and capitalize on the usage of MPI where possible.