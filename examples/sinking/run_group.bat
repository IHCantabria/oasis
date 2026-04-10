@echo off
REM ==============================================================
REM OASIS - Run all sinking examples
REM ==============================================================

set OASIS_EXE=..\..\bin\oasis.exe

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ============================================
echo   Running sinking examples
echo ============================================
echo.

set FAILED=0

echo Running progressive_flooding...
cd progressive_flooding
call run.bat
if %ERRORLEVEL% NEQ 0 set FAILED=1
cd ..

echo.
if %FAILED% NEQ 0 (
    echo Some sinking examples FAILED
    exit /b 1
) else (
    echo All sinking examples PASSED
    exit /b 0
)
