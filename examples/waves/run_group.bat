@echo off
setlocal enabledelayedexpansion
REM ==============================================================
REM OASIS - Run all waves examples on Windows
REM ==============================================================

set OASIS_EXE=..\..\bin\oasis.exe

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ============================================
echo   Running waves examples
echo ============================================
echo.

set FAILED_TESTS=
set TOTAL_TESTS=0
set PASSED_TESTS=0
set TEST_NUM=0

for %%d in (regular jonswap jonswap_piecewise multidirectional timeseries frequency_domain) do (
    set /a TEST_NUM+=1
    echo [!TEST_NUM!] Running waves/%%d...
    cd %%d
    call run.bat
    set RESULT=!ERRORLEVEL!
    cd ..
    set /a TOTAL_TESTS+=1
    if !RESULT! NEQ 0 (
        echo FAILED: waves/%%d
        set FAILED_TESTS=!FAILED_TESTS! waves/%%d
    ) else (
        echo PASSED: waves/%%d
        set /a PASSED_TESTS+=1
    )
    echo.
)

echo ============================================
if "!FAILED_TESTS!"=="" (
    echo   All waves examples passed!
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo ============================================
    exit /b 0
) else (
    echo   Some waves examples FAILED
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo   Failed tests:!FAILED_TESTS!
    echo ============================================
    exit /b 1
)
