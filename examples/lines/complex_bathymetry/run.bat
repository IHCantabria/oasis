@echo off
REM ==============================================================
REM OASIS - Run lines/complex_bathymetry example on Windows
REM ==============================================================
REM
REM Before running, generate the bathymetry mesh file:
REM   cd resources
REM   python generate_bathymetry.py
REM ==============================================================

set OASIS_EXE=..\..\..\bin\oasis.exe
set INPUT_DIR=.

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

if not exist "input\complex_bathymetry.dat" (
    echo ERROR: Bathymetry mesh not found at input\complex_bathymetry.dat
    echo Generate it first: cd resources ^& python generate_bathymetry.py
    exit /b 1
)

echo Running OASIS lines/complex_bathymetry...
"%OASIS_EXE%" "%INPUT_DIR%"
echo.
echo Exit code: %ERRORLEVEL%
exit /b %ERRORLEVEL%
