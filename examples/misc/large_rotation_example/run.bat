@echo off
REM ==============================================================
REM OASIS - Run large_rotation_example on Windows
REM Runs two cases: simplified rotation (rotSimpFlag=1) and
REM full rotation dynamics (rotSimpFlag=0), then compares results.
REM ==============================================================

set OASIS_EXE=..\..\..\bin\oasis.exe

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ====== Case 1: Simplified rotation (rotSimpFlag=1) ======
"%OASIS_EXE%" case_simplified
if %ERRORLEVEL% NEQ 0 (
    echo FAILED: case_simplified
    exit /b 1
)
echo.

echo ====== Case 2: Full rotation dynamics (rotSimpFlag=0) ======
"%OASIS_EXE%" case_full
if %ERRORLEVEL% NEQ 0 (
    echo FAILED: case_full
    exit /b 1
)
echo.

echo Both cases completed successfully.
