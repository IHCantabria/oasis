@echo off
REM ==============================================================
REM OASIS - Run all morison examples
REM ==============================================================

set OASIS_EXE=..\..\bin\oasis.exe

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ============================================
echo   Running morison examples
echo ============================================
echo.

set FAILED=0

echo Running wind_drag...
cd wind_drag
call run.bat
if %ERRORLEVEL% NEQ 0 set FAILED=1
cd ..

echo.
echo Running current_drag...
cd current_drag
call run.bat
if %ERRORLEVEL% NEQ 0 set FAILED=1
cd ..

echo.
echo Running wind_drag_yaml...
cd wind_drag_yaml
call run.bat
if %ERRORLEVEL% NEQ 0 set FAILED=1
cd ..

echo.
echo Running current_drag_yaml...
cd current_drag_yaml
call run.bat
if %ERRORLEVEL% NEQ 0 set FAILED=1
cd ..

echo.
if %FAILED% NEQ 0 (
    echo Some morison examples FAILED
    exit /b 1
) else (
    echo All morison examples PASSED
    exit /b 0
)
