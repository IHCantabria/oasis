@ECHO OFF
REM ============================================================================
REM OpenFAST Build Script for Windows
REM ============================================================================
REM This script uses CMake + Ninja + Intel oneAPI (ifx) to build OpenFAST.
REM It auto-detects CMake, Visual Studio, and Intel oneAPI installations.
REM
REM Prerequisites (auto-detected):
REM   - CMake (searched in PATH and "C:\Program Files\CMake\bin")
REM   - Visual Studio 2019+ with C/C++ workload
REM   - Intel oneAPI HPC Toolkit (provides ifx Fortran compiler + MKL)
REM
REM Usage:
REM   build_openfast.bat [build_type]
REM   build_type: Debug, Release, or clean (default: Release)
REM ============================================================================

REM --- Basic variables ---
SET "ROOT_DIR=%~dp0"
SET "BUILD_TYPE=%~1"
IF "%BUILD_TYPE%"=="" SET "BUILD_TYPE=Release"
SET "BUILD_DIR=%ROOT_DIR%build"
SET "INSTALL_DIR=%ROOT_DIR%install"

REM --- Handle clean ---
IF /I "%BUILD_TYPE%"=="clean" (
    IF EXIST "%BUILD_DIR%" rmdir /S /Q "%BUILD_DIR%"
    ECHO Build directory cleaned.
    PAUSE
    EXIT /B 0
)

REM --- Build options ---
SET "DOUBLE_PRECISION=ON"
SET "USE_DLL_INTERFACE=ON"
SET "BUILD_SHARED_LIBS=OFF"
SET "BUILD_FASTFARM=OFF"
SET "BUILD_OPENFAST_CPP_API=OFF"
SET "BUILD_TESTING=OFF"
SET "BUILD_DOCUMENTATION=OFF"
SET "FPE_TRAP_ENABLED=OFF"
SET "OPENMP=OFF"
SET "CMAKE_GENERATOR=Ninja"

REM ============================================================================
REM PHASE 1: Environment setup (linear calls, no complex control flow)
REM ============================================================================

REM --- 1. Add CMake to PATH if needed ---
cmake --version >nul 2>&1
IF %ERRORLEVEL% NEQ 0 SET "PATH=C:\Program Files\CMake\bin;%PATH%"

REM --- 2. Find and call vcvarsall.bat ---
SET "VCVARSALL="
SET "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
IF EXIST "%VSWHERE%" FOR /F "tokens=*" %%i IN ('"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath 2^>nul') DO SET "VS_INSTALL_PATH=%%i"
IF DEFINED VS_INSTALL_PATH SET "VCVARSALL=%VS_INSTALL_PATH%\VC\Auxiliary\Build\vcvarsall.bat"

REM Fallback paths if vswhere didn't find it
IF NOT DEFINED VCVARSALL IF EXIST "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvarsall.bat" SET "VCVARSALL=C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvarsall.bat"
IF NOT DEFINED VCVARSALL IF EXIST "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" SET "VCVARSALL=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat"
IF NOT DEFINED VCVARSALL IF EXIST "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" SET "VCVARSALL=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat"

IF NOT DEFINED VCVARSALL (
    ECHO ERROR: Visual Studio not found. Install VS with C++ workload.
    PAUSE
    EXIT /B 1
)

ECHO [1/4] Initializing Visual Studio: %VCVARSALL%
CALL "%VCVARSALL%" x64 >nul 2>&1

REM --- 3. Find and call Intel oneAPI setvars.bat ---
SET "ONEAPI_SETVARS="
IF EXIST "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" SET "ONEAPI_SETVARS=C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
IF NOT DEFINED ONEAPI_SETVARS IF EXIST "C:\Program Files\Intel\oneAPI\setvars.bat" SET "ONEAPI_SETVARS=C:\Program Files\Intel\oneAPI\setvars.bat"

IF NOT DEFINED ONEAPI_SETVARS (
    ECHO ERROR: Intel oneAPI not found. Install Intel oneAPI HPC Toolkit.
    PAUSE
    EXIT /B 1
)

ECHO [2/4] Initializing Intel oneAPI: %ONEAPI_SETVARS%
CALL "%ONEAPI_SETVARS%" --force
REM After setvars.bat returns, jump directly to build phase
GOTO :PHASE2

REM ============================================================================
REM PHASE 2: Build (all complex logic is here, after environment is set up)
REM ============================================================================
:PHASE2

REM --- Verify tools ---
ECHO(
ECHO ============================================================================
ECHO Verifying build tools...
ECHO ============================================================================

cmake --version >nul 2>&1
IF %ERRORLEVEL% NEQ 0 (
    ECHO ERROR: CMake not available.
    PAUSE
    EXIT /B 1
)
ECHO CMake:
cmake --version 2>&1 | FINDSTR "cmake version"

cl >nul 2>&1
ECHO MSVC: cl.exe found

ifx --version >nul 2>&1
IF %ERRORLEVEL% NEQ 0 (
    ECHO ERROR: Intel Fortran compiler ifx not available.
    ECHO Make sure Intel oneAPI HPC Toolkit is installed.
    PAUSE
    EXIT /B 1
)
ECHO Intel Fortran:
ifx --version 2>&1 | FINDSTR /I "ifx"

REM --- Find Ninja ---
ninja --version >nul 2>&1
IF %ERRORLEVEL% EQU 0 (
    SET "CMAKE_GENERATOR=Ninja"
    ECHO Ninja:
    ninja --version
    GOTO :ninja_done
)

REM Try to find Ninja bundled with Visual Studio
IF DEFINED VS_INSTALL_PATH IF EXIST "%VS_INSTALL_PATH%\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja\ninja.exe" (
    SET "PATH=%VS_INSTALL_PATH%\Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja;%PATH%"
    SET "CMAKE_GENERATOR=Ninja"
    ECHO Ninja: found in VS installation
    GOTO :ninja_done
)

ECHO Ninja not found, using NMake Makefiles
SET "CMAKE_GENERATOR=NMake Makefiles"

:ninja_done
ECHO Generator: %CMAKE_GENERATOR%
ECHO(

REM --- Setup build directory ---
ECHO ============================================================================
ECHO [3/4] Configuring OpenFAST with CMake...
ECHO ============================================================================

IF EXIST "%BUILD_DIR%\CMakeCache.txt" (
    ECHO Clearing old CMake cache...
    DEL /F /Q "%BUILD_DIR%\CMakeCache.txt" >nul 2>&1
    IF EXIST "%BUILD_DIR%\CMakeFiles" rmdir /S /Q "%BUILD_DIR%\CMakeFiles" >nul 2>&1
)

IF NOT EXIST "%BUILD_DIR%" MKDIR "%BUILD_DIR%"

ECHO   Build Type: %BUILD_TYPE%
ECHO   Generator:  %CMAKE_GENERATOR%
ECHO   Install:    %INSTALL_DIR%
ECHO(

CD /D "%BUILD_DIR%"

cmake -G "%CMAKE_GENERATOR%" ^
    -DCMAKE_BUILD_TYPE=%BUILD_TYPE% ^
    -DCMAKE_Fortran_COMPILER=ifx ^
    -DCMAKE_C_COMPILER=cl ^
    -DCMAKE_CXX_COMPILER=cl ^
    -DCMAKE_INSTALL_PREFIX="%INSTALL_DIR%" ^
    -DDOUBLE_PRECISION=%DOUBLE_PRECISION% ^
    -DUSE_DLL_INTERFACE=%USE_DLL_INTERFACE% ^
    -DBUILD_SHARED_LIBS=%BUILD_SHARED_LIBS% ^
    -DBUILD_FASTFARM=%BUILD_FASTFARM% ^
    -DBUILD_OPENFAST_CPP_API=%BUILD_OPENFAST_CPP_API% ^
    -DBUILD_TESTING=%BUILD_TESTING% ^
    -DBUILD_DOCUMENTATION=%BUILD_DOCUMENTATION% ^
    -DFPE_TRAP_ENABLED=%FPE_TRAP_ENABLED% ^
    -DOPENMP=%OPENMP% ^
    ..

IF %ERRORLEVEL% NEQ 0 (
    ECHO(
    ECHO ERROR: CMake configuration failed!
    ECHO   - Check compiler and library availability
    ECHO   - Try: build_openfast.bat clean
    PAUSE
    EXIT /B 1
)

ECHO(
ECHO Configuration OK.
ECHO(

REM --- Build ---
ECHO ============================================================================
ECHO [4/4] Building OpenFAST... (this may take several minutes)
ECHO ============================================================================

cmake --build . --config %BUILD_TYPE%

IF %ERRORLEVEL% NEQ 0 (
    ECHO(
    ECHO ERROR: Build failed! Check errors above.
    PAUSE
    EXIT /B 1
)

ECHO(
ECHO ============================================================================
ECHO Build completed successfully!
ECHO ============================================================================

REM --- Install ---
ECHO(
ECHO Installing...
cmake --build . --config %BUILD_TYPE% --target install

IF %ERRORLEVEL% NEQ 0 (
    ECHO WARNING: Installation failed. Binaries may still be in build directory.
) ELSE (
    ECHO(
    ECHO Installed to: %INSTALL_DIR%
    ECHO   Executables: %INSTALL_DIR%\bin
    ECHO   Libraries:   %INSTALL_DIR%\lib
)

ECHO(
ECHO ============================================================================
ECHO Done. Run: %INSTALL_DIR%\bin\openfast.exe [input_file]
ECHO ============================================================================

CD /D "%ROOT_DIR%"
PAUSE
EXIT /B 0
