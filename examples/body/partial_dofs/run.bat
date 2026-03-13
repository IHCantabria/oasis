@echo off
REM ==============================================================
REM OASIS - Run body/partial_dofs example on Windows
REM ==============================================================

set OASIS_EXE=..\..\..\bin\oasis.exe
set INPUT_DIR=.

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo Running OASIS body/partial_dofs...
"%OASIS_EXE%" "%INPUT_DIR%"
echo.
echo Exit code: %ERRORLEVEL%
