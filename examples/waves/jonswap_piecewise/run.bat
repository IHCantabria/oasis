@echo off
REM ==============================================================
REM OASIS - Run waves/jonswap_piecewise example on Windows
REM ==============================================================

set OASIS_EXE=..\..\..\bin\oasis.exe
set INPUT_DIR=.

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo Running OASIS waves/jonswap_piecewise...
"%OASIS_EXE%" "%INPUT_DIR%"
echo.
echo Exit code: %ERRORLEVEL%
