@echo off
setlocal enabledelayedexpansion
REM ==============================================================
REM OASIS - Run all owcs examples on Windows
REM ==============================================================

set OASIS_EXE=..\..\bin\oasis.exe

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ============================================
echo   Running owcs examples
echo ============================================
echo.

set FAILED_TESTS=
set TOTAL_TESTS=0
set PASSED_TESTS=0
set TEST_NUM=0

for %%d in (free_decay excitation hole_coupling turbine_coupling) do (
    set /a TEST_NUM+=1
    echo [!TEST_NUM!] Running owcs/%%d...
    cd %%d
    call run.bat
    set RESULT=!ERRORLEVEL!
    cd ..
    set /a TOTAL_TESTS+=1
    if !RESULT! NEQ 0 (
        echo FAILED: owcs/%%d
        set FAILED_TESTS=!FAILED_TESTS! owcs/%%d
    ) else (
        echo PASSED: owcs/%%d
        set /a PASSED_TESTS+=1
    )
    echo.
)

echo ============================================
if "!FAILED_TESTS!"=="" (
    echo   All owcs examples passed!
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo ============================================
    exit /b 0
) else (
    echo   Some owcs examples FAILED
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo   Failed tests:!FAILED_TESTS!
    echo ============================================
    exit /b 1
)
