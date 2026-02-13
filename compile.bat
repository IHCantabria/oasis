
@echo off

rem ============================================================
rem OASIS - Windows build script
rem ============================================================
rem Prerequisites:
rem   - CMake on PATH
rem   - MSVC (run from "Developer Command Prompt" or "x64 Native Tools")
rem   - Libraries: Armadillo, HDF5, OpenBLAS (or MKL)
rem   - Optional: SuperLU, stl_reader
rem
rem If libraries are not auto-detected, create CMakeUserConfig.cmake
rem from CMakeUserConfig.cmake.template and set the paths.
rem ============================================================

rem DELETE PREVIOUS BUILD ---------------------------------
rem rmdir /Q /S build

rem CONFIGURE ---------------------------------------------
cmake -B./build/ -H./

rem CHOOSE COMPILATION OPTION -----------------------------
rem DEBUG COMPILATION ---------------------------------
rem cmake --build ./build
rem copy ".\build\Debug\OASIS.exe" ".\bin\oasis.exe"

rem RELEASE COMPILATION -------------------------------
cmake --build ./build --config Release

if not exist bin mkdir bin
copy ".\build\Release\OASIS.exe" ".\bin\oasis.exe"
echo Build complete: bin\oasis.exe