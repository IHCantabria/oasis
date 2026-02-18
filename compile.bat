
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
rem Paths are configured in CMakeUserConfig.cmake.
rem Copy CMakeUserConfig.cmake.template and set paths for your system.
rem ============================================================

rem --- User-configurable paths (match CMakeUserConfig.cmake) ---
set FASTURBINE_ROOT=E:\00_GIT_REPOS\FASTurbine_wrapper\install
set INTEL_ONEAPI_ROOT=C:\Program Files (x86)\Intel\oneAPI\compiler\latest

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

if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Build failed!
    exit /b 1
)

if not exist bin mkdir bin
copy ".\build\Release\OASIS.exe" ".\bin\oasis.exe"
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Failed to copy executable!
    exit /b 1
)

rem Copy required DLLs for OpenFAST (if present)
if exist "%FASTURBINE_ROOT%\lib\fasturbwrapper.dll" (
    copy "%FASTURBINE_ROOT%\lib\fasturbwrapper.dll" ".\bin\" >nul
    echo Copied fasturbwrapper.dll to bin\
)
rem Copy Intel oneAPI Fortran runtime DLLs (required by OpenFAST)
set INTEL_BIN=%INTEL_ONEAPI_ROOT%\bin
if exist "%INTEL_BIN%\libifcoremd.dll" (
    copy "%INTEL_BIN%\libifcoremd.dll" ".\bin\" >nul
    copy "%INTEL_BIN%\libifportmd.dll" ".\bin\" >nul
    copy "%INTEL_BIN%\libmmd.dll" ".\bin\" >nul
    copy "%INTEL_BIN%\svml_dispmd.dll" ".\bin\" >nul
    echo Copied Intel Fortran runtime DLLs to bin\
)
echo Build complete: bin\oasis.exe