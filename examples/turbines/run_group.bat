@echo off
setlocal enabledelayedexpansion
REM ==============================================================
REM OASIS - Run all turbines examples on Windows
REM ==============================================================

set OASIS_EXE=..\..\bin\oasis.exe

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ============================================
echo   Running turbines examples
echo ============================================
echo.

set FAILED_TESTS=
set TOTAL_TESTS=0
set PASSED_TESTS=0
set TEST_NUM=0

for %%d in (fixed_turbine moored_turbine) do (
    set /a TEST_NUM+=1
    echo [!TEST_NUM!] Running turbines/%%d...
    cd %%d
    call run.bat
    set RESULT=!ERRORLEVEL!
    cd ..
    set /a TOTAL_TESTS+=1
    if !RESULT! NEQ 0 (
        echo FAILED: turbines/%%d
        set FAILED_TESTS=!FAILED_TESTS! turbines/%%d
    ) else (
        echo PASSED: turbines/%%d
        set /a PASSED_TESTS+=1
    )
    echo.
)

echo ============================================
if "!FAILED_TESTS!"=="" (
    echo   All turbines examples passed!
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo ============================================
    exit /b 0
) else (
    echo   Some turbines examples FAILED
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo   Failed tests:!FAILED_TESTS!
    echo ============================================
    exit /b 1
)
