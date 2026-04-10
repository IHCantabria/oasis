@echo off
setlocal enabledelayedexpansion
REM ==============================================================
REM OASIS - Run all misc examples on Windows
REM ==============================================================

set OASIS_EXE=..\..\bin\oasis.exe

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ============================================
echo   Running misc examples
echo ============================================
echo.

set FAILED_TESTS=
set TOTAL_TESTS=0
set PASSED_TESTS=0
set TEST_NUM=0

for %%d in (freq_wave_example generic_example qtf_example turbine_example elastic_anchor_example large_rotation_example) do (
    set /a TEST_NUM+=1
    echo [!TEST_NUM!] Running misc/%%d...
    cd %%d
    call run.bat
    set RESULT=!ERRORLEVEL!
    cd ..
    set /a TOTAL_TESTS+=1
    if !RESULT! NEQ 0 (
        echo FAILED: misc/%%d
        set FAILED_TESTS=!FAILED_TESTS! misc/%%d
    ) else (
        echo PASSED: misc/%%d
        set /a PASSED_TESTS+=1
    )
    echo.
)

echo ============================================
if "!FAILED_TESTS!"=="" (
    echo   All misc examples passed!
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo ============================================
    exit /b 0
) else (
    echo   Some misc examples FAILED
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo   Failed tests:!FAILED_TESTS!
    echo ============================================
    exit /b 1
)
