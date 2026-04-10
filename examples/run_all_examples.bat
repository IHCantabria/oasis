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
set TEST_NUM=0

REM ==============================================================
REM Body behavior examples
REM ==============================================================

REM --- body/free_decay ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/free_decay...
cd body\free_decay
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/free_decay
    set FAILED_TESTS=%FAILED_TESTS% body/free_decay
) else (
    echo PASSED: body/free_decay
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/fixed ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/fixed...
cd body\fixed
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/fixed
    set FAILED_TESTS=%FAILED_TESTS% body/fixed
) else (
    echo PASSED: body/fixed
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/partial_dofs ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/partial_dofs...
cd body\partial_dofs
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/partial_dofs
    set FAILED_TESTS=%FAILED_TESTS% body/partial_dofs
) else (
    echo PASSED: body/partial_dofs
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/imposed_motion ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/imposed_motion...
cd body\imposed_motion
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/imposed_motion
    set FAILED_TESTS=%FAILED_TESTS% body/imposed_motion
) else (
    echo PASSED: body/imposed_motion
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/radiation ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/radiation...
cd body\radiation
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/radiation
    set FAILED_TESTS=%FAILED_TESTS% body/radiation
) else (
    echo PASSED: body/radiation
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/excitation_linear ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/excitation_linear...
cd body\excitation_linear
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/excitation_linear
    set FAILED_TESTS=%FAILED_TESTS% body/excitation_linear
) else (
    echo PASSED: body/excitation_linear
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/excitation_instantpos ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/excitation_instantpos...
cd body\excitation_instantpos
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/excitation_instantpos
    set FAILED_TESTS=%FAILED_TESTS% body/excitation_instantpos
) else (
    echo PASSED: body/excitation_instantpos
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/nonlinear_hs_flat ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/nonlinear_hs_flat...
cd body\nonlinear_hs_flat
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/nonlinear_hs_flat
    set FAILED_TESTS=%FAILED_TESTS% body/nonlinear_hs_flat
) else (
    echo PASSED: body/nonlinear_hs_flat
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/nonlinear_hs_waves ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/nonlinear_hs_waves...
cd body\nonlinear_hs_waves
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/nonlinear_hs_waves
    set FAILED_TESTS=%FAILED_TESTS% body/nonlinear_hs_waves
) else (
    echo PASSED: body/nonlinear_hs_waves
    set /a PASSED_TESTS+=1
)
echo.

REM --- body/qtf ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running body/qtf...
cd body\qtf
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: body/qtf
    set FAILED_TESTS=%FAILED_TESTS% body/qtf
) else (
    echo PASSED: body/qtf
    set /a PASSED_TESTS+=1
)
echo.

REM ==============================================================
REM Wave type examples
REM ==============================================================

REM --- waves/regular ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running waves/regular...
cd waves\regular
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: waves/regular
    set FAILED_TESTS=%FAILED_TESTS% waves/regular
) else (
    echo PASSED: waves/regular
    set /a PASSED_TESTS+=1
)
echo.

REM --- waves/jonswap ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running waves/jonswap...
cd waves\jonswap
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: waves/jonswap
    set FAILED_TESTS=%FAILED_TESTS% waves/jonswap
) else (
    echo PASSED: waves/jonswap
    set /a PASSED_TESTS+=1
)
echo.

REM --- waves/jonswap_piecewise ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running waves/jonswap_piecewise...
cd waves\jonswap_piecewise
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: waves/jonswap_piecewise
    set FAILED_TESTS=%FAILED_TESTS% waves/jonswap_piecewise
) else (
    echo PASSED: waves/jonswap_piecewise
    set /a PASSED_TESTS+=1
)
echo.

REM --- waves/multidirectional ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running waves/multidirectional...
cd waves\multidirectional
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: waves/multidirectional
    set FAILED_TESTS=%FAILED_TESTS% waves/multidirectional
) else (
    echo PASSED: waves/multidirectional
    set /a PASSED_TESTS+=1
)
echo.

REM --- waves/timeseries ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running waves/timeseries...
cd waves\timeseries
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: waves/timeseries
    set FAILED_TESTS=%FAILED_TESTS% waves/timeseries
) else (
    echo PASSED: waves/timeseries
    set /a PASSED_TESTS+=1
)
echo.

REM --- waves/frequency_domain ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running waves/frequency_domain...
cd waves\frequency_domain
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: waves/frequency_domain
    set FAILED_TESTS=%FAILED_TESTS% waves/frequency_domain
) else (
    echo PASSED: waves/frequency_domain
    set /a PASSED_TESTS+=1
)
echo.

REM ==============================================================
REM Mooring lines examples
REM ==============================================================

REM --- lines/single_line ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/single_line...
cd lines\single_line
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/single_line
    set FAILED_TESTS=%FAILED_TESTS% lines/single_line
) else (
    echo PASSED: lines/single_line
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/multi_line ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/multi_line...
cd lines\multi_line
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/multi_line
    set FAILED_TESTS=%FAILED_TESTS% lines/multi_line
) else (
    echo PASSED: lines/multi_line
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/joint_connection ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/joint_connection...
cd lines\joint_connection
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/joint_connection
    set FAILED_TESTS=%FAILED_TESTS% lines/joint_connection
) else (
    echo PASSED: lines/joint_connection
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/elastic_anchor ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/elastic_anchor...
cd lines\elastic_anchor
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/elastic_anchor
    set FAILED_TESTS=%FAILED_TESTS% lines/elastic_anchor
) else (
    echo PASSED: lines/elastic_anchor
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/prescribed_motion ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/prescribed_motion...
cd lines\prescribed_motion
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/prescribed_motion
    set FAILED_TESTS=%FAILED_TESTS% lines/prescribed_motion
) else (
    echo PASSED: lines/prescribed_motion
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/viscoelastic ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/viscoelastic...
cd lines\viscoelastic
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/viscoelastic
    set FAILED_TESTS=%FAILED_TESTS% lines/viscoelastic
) else (
    echo PASSED: lines/viscoelastic
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/tabulated_stiffness ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/tabulated_stiffness...
cd lines\tabulated_stiffness
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/tabulated_stiffness
    set FAILED_TESTS=%FAILED_TESTS% lines/tabulated_stiffness
) else (
    echo PASSED: lines/tabulated_stiffness
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/tension_symmetric ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/tension_symmetric...
cd lines\tension_symmetric
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/tension_symmetric
    set FAILED_TESTS=%FAILED_TESTS% lines/tension_symmetric
) else (
    echo PASSED: lines/tension_symmetric
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/seabed_contact ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/seabed_contact...
cd lines\seabed_contact
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/seabed_contact
    set FAILED_TESTS=%FAILED_TESTS% lines/seabed_contact
) else (
    echo PASSED: lines/seabed_contact
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/friction_isotropic ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/friction_isotropic...
cd lines\friction_isotropic
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/friction_isotropic
    set FAILED_TESTS=%FAILED_TESTS% lines/friction_isotropic
) else (
    echo PASSED: lines/friction_isotropic
    set /a PASSED_TESTS+=1
)
echo.

REM --- lines/friction_anisotropic ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running lines/friction_anisotropic...
cd lines\friction_anisotropic
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: lines/friction_anisotropic
    set FAILED_TESTS=%FAILED_TESTS% lines/friction_anisotropic
) else (
    echo PASSED: lines/friction_anisotropic
    set /a PASSED_TESTS+=1
)
echo.

REM ==============================================================
REM Wind turbine examples
REM ==============================================================

REM --- turbines/fixed_turbine ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running turbines/fixed_turbine...
cd turbines\fixed_turbine
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: turbines/fixed_turbine
    set FAILED_TESTS=%FAILED_TESTS% turbines/fixed_turbine
) else (
    echo PASSED: turbines/fixed_turbine
    set /a PASSED_TESTS+=1
)
echo.

REM --- turbines/moored_turbine ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running turbines/moored_turbine...
cd turbines\moored_turbine
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: turbines/moored_turbine
    set FAILED_TESTS=%FAILED_TESTS% turbines/moored_turbine
) else (
    echo PASSED: turbines/moored_turbine
    set /a PASSED_TESTS+=1
)
echo.

REM ==============================================================
REM Spring examples
REM ==============================================================

REM --- springs/pile_connector ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running springs/pile_connector...
cd springs\pile_connector
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: springs/pile_connector
    set FAILED_TESTS=%FAILED_TESTS% springs/pile_connector
) else (
    echo PASSED: springs/pile_connector
    set /a PASSED_TESTS+=1
)
echo.

REM --- springs/neoprene_connector ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running springs/neoprene_connector...
cd springs\neoprene_connector
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: springs/neoprene_connector
    set FAILED_TESTS=%FAILED_TESTS% springs/neoprene_connector
) else (
    echo PASSED: springs/neoprene_connector
    set /a PASSED_TESTS+=1
)
echo.

REM ==============================================================
REM Miscellaneous examples (pre-existing)
REM ==============================================================

REM --- misc/freq_wave_example ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running misc/freq_wave_example...
cd misc\freq_wave_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: misc/freq_wave_example
    set FAILED_TESTS=%FAILED_TESTS% misc/freq_wave_example
) else (
    echo PASSED: misc/freq_wave_example
    set /a PASSED_TESTS+=1
)
echo.

REM --- misc/generic_example ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running misc/generic_example...
cd misc\generic_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: misc/generic_example
    set FAILED_TESTS=%FAILED_TESTS% misc/generic_example
) else (
    echo PASSED: misc/generic_example
    set /a PASSED_TESTS+=1
)
echo.

REM --- misc/qtf_example ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running misc/qtf_example...
cd misc\qtf_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: misc/qtf_example
    set FAILED_TESTS=%FAILED_TESTS% misc/qtf_example
) else (
    echo PASSED: misc/qtf_example
    set /a PASSED_TESTS+=1
)
echo.

REM --- misc/turbine_example ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running misc/turbine_example...
cd misc\turbine_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: misc/turbine_example
    set FAILED_TESTS=%FAILED_TESTS% misc/turbine_example
) else (
    echo PASSED: misc/turbine_example
    set /a PASSED_TESTS+=1
)
echo.

REM --- misc/elastic_anchor_example ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running misc/elastic_anchor_example...
cd misc\elastic_anchor_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: misc/elastic_anchor_example
    set FAILED_TESTS=%FAILED_TESTS% misc/elastic_anchor_example
) else (
    echo PASSED: misc/elastic_anchor_example
    set /a PASSED_TESTS+=1
)
echo.

REM --- misc/large_rotation_example ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running misc/large_rotation_example...
cd misc\large_rotation_example
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: misc/large_rotation_example
    set FAILED_TESTS=%FAILED_TESTS% misc/large_rotation_example
) else (
    echo PASSED: misc/large_rotation_example
    set /a PASSED_TESTS+=1
)
echo.

REM ==============================================================
REM Multi-body examples
REM ==============================================================

REM --- multibody/shared_hydb ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running multibody/shared_hydb...
cd multibody\shared_hydb
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: multibody/shared_hydb
    set FAILED_TESTS=%FAILED_TESTS% multibody/shared_hydb
) else (
    echo PASSED: multibody/shared_hydb
    set /a PASSED_TESTS+=1
)
echo.

REM --- multibody/separate_hydb ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running multibody/separate_hydb...
cd multibody\separate_hydb
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: multibody/separate_hydb
    set FAILED_TESTS=%FAILED_TESTS% multibody/separate_hydb
) else (
    echo PASSED: multibody/separate_hydb
    set /a PASSED_TESTS+=1
)
echo.

REM ==============================================================
REM OWC examples
REM ==============================================================

REM --- owcs/free_decay ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running owcs/free_decay...
cd owcs\free_decay
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: owcs/free_decay
    set FAILED_TESTS=%FAILED_TESTS% owcs/free_decay
) else (
    echo PASSED: owcs/free_decay
    set /a PASSED_TESTS+=1
)
echo.

REM --- owcs/excitation ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running owcs/excitation...
cd owcs\excitation
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: owcs/excitation
    set FAILED_TESTS=%FAILED_TESTS% owcs/excitation
) else (
    echo PASSED: owcs/excitation
    set /a PASSED_TESTS+=1
)
echo.

REM --- owcs/hole_coupling ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running owcs/hole_coupling...
cd owcs\hole_coupling
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: owcs/hole_coupling
    set FAILED_TESTS=%FAILED_TESTS% owcs/hole_coupling
) else (
    echo PASSED: owcs/hole_coupling
    set /a PASSED_TESTS+=1
)
echo.

REM --- owcs/turbine_coupling ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running owcs/turbine_coupling...
cd owcs\turbine_coupling
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: owcs/turbine_coupling
    set FAILED_TESTS=%FAILED_TESTS% owcs/turbine_coupling
) else (
    echo PASSED: owcs/turbine_coupling
    set /a PASSED_TESTS+=1
)
echo.

REM ==============================================================
REM Sinking examples
REM ==============================================================

REM --- sinking/progressive_flooding ---
set /a TEST_NUM+=1
echo [%TEST_NUM%] Running sinking/progressive_flooding...
cd sinking\progressive_flooding
call run.bat
set RESULT=%ERRORLEVEL%
cd ..\..
set /a TOTAL_TESTS+=1
if %RESULT% NEQ 0 (
    echo FAILED: sinking/progressive_flooding
    set FAILED_TESTS=%FAILED_TESTS% sinking/progressive_flooding
) else (
    echo PASSED: sinking/progressive_flooding
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
