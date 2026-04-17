@echo off
setlocal enabledelayedexpansion
REM ==============================================================
REM OASIS - Run all body examples on Windows
REM ==============================================================

set OASIS_EXE=..\..\bin\oasis.exe

if not exist "%OASIS_EXE%" (
    echo ERROR: OASIS executable not found at %OASIS_EXE%
    echo Build OASIS first with compile.bat
    exit /b 1
)

echo ============================================
echo   Running body examples
echo ============================================
echo.

set FAILED_TESTS=
set TOTAL_TESTS=0
set PASSED_TESTS=0
set TEST_NUM=0

for %%d in (free_decay fixed partial_dofs imposed_motion radiation excitation_linear excitation_instantpos nonlinear_hs_flat nonlinear_hs_waves qtf fixed_yaml free_decay_yaml radiation_yaml excitation_linear_yaml imposed_motion_yaml nonlinear_hs_flat_yaml partial_dofs_yaml) do (
    set /a TEST_NUM+=1
    echo [!TEST_NUM!] Running body/%%d...
    cd %%d
    call run.bat
    set RESULT=!ERRORLEVEL!
    cd ..
    set /a TOTAL_TESTS+=1
    if !RESULT! NEQ 0 (
        echo FAILED: body/%%d
        set FAILED_TESTS=!FAILED_TESTS! body/%%d
    ) else (
        echo PASSED: body/%%d
        set /a PASSED_TESTS+=1
    )
    echo.
)

echo ============================================
if "!FAILED_TESTS!"=="" (
    echo   All body examples passed!
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo ============================================
    exit /b 0
) else (
    echo   Some body examples FAILED
    echo   !PASSED_TESTS!/!TOTAL_TESTS! tests passed
    echo   Failed tests:!FAILED_TESTS!
    echo ============================================
    exit /b 1
)
