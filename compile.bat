
@echo off

rem ============================================================
rem OASIS - Windows build script
rem ============================================================
rem Prerequisites:
rem   - CMake on PATH
rem   - MSVC (run from "Developer Command Prompt" or "x64 Native Tools")
rem   - Libraries: Armadillo, HDF5, OpenBLAS (or MKL)
rem   - Optional: SuperLU, stl_reader, OpenFAST, FASTurbine_wrapper
rem
rem ALL paths are configured in CMakeUserConfig.cmake only.
rem Copy CMakeUserConfig.cmake.template and set paths for your system.
rem ============================================================

rem DELETE PREVIOUS BUILD ---------------------------------
rem rmdir /Q /S build

rem CONFIGURE ---------------------------------------------
rem Toolchain file is auto-detected from CMakeUserConfig.cmake.
cmake -B build -S .
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: CMake configure failed!
    exit /b 1
)

rem CHOOSE COMPILATION OPTION -----------------------------
rem DEBUG COMPILATION ---------------------------------
rem cmake --build build
rem copy ".\build\Debug\OASIS.exe" ".\bin\oasis.exe"

rem RELEASE COMPILATION -------------------------------
cmake --build build --config Release
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Build failed!
    exit /b 1
)

rem Copy executable + all runtime DLLs from build output to bin/
rem (vcpkg auto-copies its DLLs to build/Release/ via VCPKG_APPLOCAL_DEPS;
rem  CMakeLists.txt post-build commands copy FASTurbine + Intel DLLs there too)
if not exist bin mkdir bin
for %%f in (".\build\Release\*.exe" ".\build\Release\*.dll") do (
    copy "%%f" ".\bin\" >nul
)
echo Build complete: bin\oasis.exe