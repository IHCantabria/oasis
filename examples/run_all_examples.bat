@echo off
REM ==============================================================
REM OASIS - Run all examples on Windows
REM ==============================================================

if not exist "..\bin\oasis.exe" (
    echo ERROR: OASIS executable not found at ..\bin\oasis.exe
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ============================================
echo   Running all OASIS examples
echo ============================================
echo.

set FAILED_TESTS=
set TOTAL_TESTS=0
set PASSED_TESTS=0

REM --- freq_wave_example ---
echo [1/4] Running freq_wave_example...
cd freq_wave_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: freq_wave_example
    set FAILED_TESTS=%FAILED_TESTS% freq_wave_example
) else (
    echo PASSED: freq_wave_example
    set /a PASSED_TESTS+=1
)
echo.

REM --- generic_example ---
echo [2/4] Running generic_example...
cd generic_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: generic_example
    set FAILED_TESTS=%FAILED_TESTS% generic_example
) else (
    echo PASSED: generic_example
    set /a PASSED_TESTS+=1
)
echo.

REM --- qtf_example ---
echo [3/4] Running qtf_example...
cd qtf_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: qtf_example
    set FAILED_TESTS=%FAILED_TESTS% qtf_example
) else (
    echo PASSED: qtf_example
    set /a PASSED_TESTS+=1
)
echo.

REM --- turbine_example ---
echo [4/4] Running turbine_example...
cd turbine_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: turbine_example
    set FAILED_TESTS=%FAILED_TESTS% turbine_example
) else (
    echo PASSED: turbine_example
    set /a PASSED_TESTS+=1
)
echo.

echo ============================================
if "%FAILED_TESTS%"=="" (
    echo   All examples completed successfully!
    echo   %PASSED_TESTS%/%TOTAL_TESTS% tests passed
    echo ============================================
    exit /b 0
) else (
    echo   Some examples FAILED
    echo   %PASSED_TESTS%/%TOTAL_TESTS% tests passed
    echo   Failed tests:%FAILED_TESTS%
    echo ============================================
    exit /b 1
)
