#!/bin/bash
# ==============================================================
# OASIS - Run all examples on Linux/cluster
# ==============================================================

echo "============================================"
echo "  Running all OASIS examples"
echo "============================================"
echo

FAILED_TESTS=""
TOTAL_TESTS=0
PASSED_TESTS=0
TEST_NUM=0

# ==============================================================
# Body behavior examples
# ==============================================================

# --- body/free_decay ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/free_decay..."
cd body/free_decay || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/free_decay"
    FAILED_TESTS="$FAILED_TESTS body/free_decay"
else
    echo "PASSED: body/free_decay"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/fixed ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/fixed..."
cd body/fixed || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/fixed"
    FAILED_TESTS="$FAILED_TESTS body/fixed"
else
    echo "PASSED: body/fixed"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/partial_dofs ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/partial_dofs..."
cd body/partial_dofs || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/partial_dofs"
    FAILED_TESTS="$FAILED_TESTS body/partial_dofs"
else
    echo "PASSED: body/partial_dofs"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/imposed_motion ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/imposed_motion..."
cd body/imposed_motion || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/imposed_motion"
    FAILED_TESTS="$FAILED_TESTS body/imposed_motion"
else
    echo "PASSED: body/imposed_motion"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/radiation ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/radiation..."
cd body/radiation || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/radiation"
    FAILED_TESTS="$FAILED_TESTS body/radiation"
else
    echo "PASSED: body/radiation"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/excitation_linear ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/excitation_linear..."
cd body/excitation_linear || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/excitation_linear"
    FAILED_TESTS="$FAILED_TESTS body/excitation_linear"
else
    echo "PASSED: body/excitation_linear"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/excitation_instantpos ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/excitation_instantpos..."
cd body/excitation_instantpos || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/excitation_instantpos"
    FAILED_TESTS="$FAILED_TESTS body/excitation_instantpos"
else
    echo "PASSED: body/excitation_instantpos"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/nonlinear_hs_flat ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/nonlinear_hs_flat..."
cd body/nonlinear_hs_flat || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/nonlinear_hs_flat"
    FAILED_TESTS="$FAILED_TESTS body/nonlinear_hs_flat"
else
    echo "PASSED: body/nonlinear_hs_flat"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/nonlinear_hs_waves ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/nonlinear_hs_waves..."
cd body/nonlinear_hs_waves || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/nonlinear_hs_waves"
    FAILED_TESTS="$FAILED_TESTS body/nonlinear_hs_waves"
else
    echo "PASSED: body/nonlinear_hs_waves"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- body/qtf ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running body/qtf..."
cd body/qtf || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: body/qtf"
    FAILED_TESTS="$FAILED_TESTS body/qtf"
else
    echo "PASSED: body/qtf"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# ==============================================================
# Mooring lines examples
# ==============================================================

# --- lines/single_line ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/single_line..."
cd lines/single_line || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/single_line"
    FAILED_TESTS="$FAILED_TESTS lines/single_line"
else
    echo "PASSED: lines/single_line"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/multi_line ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/multi_line..."
cd lines/multi_line || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/multi_line"
    FAILED_TESTS="$FAILED_TESTS lines/multi_line"
else
    echo "PASSED: lines/multi_line"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/joint_connection ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/joint_connection..."
cd lines/joint_connection || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/joint_connection"
    FAILED_TESTS="$FAILED_TESTS lines/joint_connection"
else
    echo "PASSED: lines/joint_connection"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/elastic_anchor ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/elastic_anchor..."
cd lines/elastic_anchor || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/elastic_anchor"
    FAILED_TESTS="$FAILED_TESTS lines/elastic_anchor"
else
    echo "PASSED: lines/elastic_anchor"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/prescribed_motion ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/prescribed_motion..."
cd lines/prescribed_motion || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/prescribed_motion"
    FAILED_TESTS="$FAILED_TESTS lines/prescribed_motion"
else
    echo "PASSED: lines/prescribed_motion"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/viscoelastic ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/viscoelastic..."
cd lines/viscoelastic || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/viscoelastic"
    FAILED_TESTS="$FAILED_TESTS lines/viscoelastic"
else
    echo "PASSED: lines/viscoelastic"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/tabulated_stiffness ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/tabulated_stiffness..."
cd lines/tabulated_stiffness || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/tabulated_stiffness"
    FAILED_TESTS="$FAILED_TESTS lines/tabulated_stiffness"
else
    echo "PASSED: lines/tabulated_stiffness"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/tension_symmetric ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/tension_symmetric..."
cd lines/tension_symmetric || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/tension_symmetric"
    FAILED_TESTS="$FAILED_TESTS lines/tension_symmetric"
else
    echo "PASSED: lines/tension_symmetric"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/seabed_contact ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/seabed_contact..."
cd lines/seabed_contact || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/seabed_contact"
    FAILED_TESTS="$FAILED_TESTS lines/seabed_contact"
else
    echo "PASSED: lines/seabed_contact"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/friction_isotropic ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/friction_isotropic..."
cd lines/friction_isotropic || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/friction_isotropic"
    FAILED_TESTS="$FAILED_TESTS lines/friction_isotropic"
else
    echo "PASSED: lines/friction_isotropic"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- lines/friction_anisotropic ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running lines/friction_anisotropic..."
cd lines/friction_anisotropic || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: lines/friction_anisotropic"
    FAILED_TESTS="$FAILED_TESTS lines/friction_anisotropic"
else
    echo "PASSED: lines/friction_anisotropic"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# ==============================================================
# Wind turbine examples
# ==============================================================

# --- turbines/fixed_turbine ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running turbines/fixed_turbine..."
cd turbines/fixed_turbine || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: turbines/fixed_turbine"
    FAILED_TESTS="$FAILED_TESTS turbines/fixed_turbine"
else
    echo "PASSED: turbines/fixed_turbine"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- turbines/moored_turbine ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running turbines/moored_turbine..."
cd turbines/moored_turbine || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: turbines/moored_turbine"
    FAILED_TESTS="$FAILED_TESTS turbines/moored_turbine"
else
    echo "PASSED: turbines/moored_turbine"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# ==============================================================
# Spring examples
# ==============================================================

# --- springs/pile_connector ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running springs/pile_connector..."
cd springs/pile_connector || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: springs/pile_connector"
    FAILED_TESTS="$FAILED_TESTS springs/pile_connector"
else
    echo "PASSED: springs/pile_connector"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- springs/neoprene_connector ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running springs/neoprene_connector..."
cd springs/neoprene_connector || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: springs/neoprene_connector"
    FAILED_TESTS="$FAILED_TESTS springs/neoprene_connector"
else
    echo "PASSED: springs/neoprene_connector"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# ==============================================================
# Miscellaneous examples (pre-existing)
# ==============================================================

# --- misc/freq_wave_example ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running misc/freq_wave_example..."
cd misc/freq_wave_example || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: misc/freq_wave_example"
    FAILED_TESTS="$FAILED_TESTS misc/freq_wave_example"
else
    echo "PASSED: misc/freq_wave_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- misc/generic_example ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running misc/generic_example..."
cd misc/generic_example || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: misc/generic_example"
    FAILED_TESTS="$FAILED_TESTS misc/generic_example"
else
    echo "PASSED: misc/generic_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- misc/qtf_example ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running misc/qtf_example..."
cd misc/qtf_example || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: misc/qtf_example"
    FAILED_TESTS="$FAILED_TESTS misc/qtf_example"
else
    echo "PASSED: misc/qtf_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- misc/turbine_example ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running misc/turbine_example..."
cd misc/turbine_example || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: misc/turbine_example"
    FAILED_TESTS="$FAILED_TESTS misc/turbine_example"
else
    echo "PASSED: misc/turbine_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- misc/elastic_anchor_example ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running misc/elastic_anchor_example..."
cd misc/elastic_anchor_example || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: misc/elastic_anchor_example"
    FAILED_TESTS="$FAILED_TESTS misc/elastic_anchor_example"
else
    echo "PASSED: misc/elastic_anchor_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- misc/large_rotation_example ---
TEST_NUM=$((TEST_NUM + 1))
echo "[$TEST_NUM] Running misc/large_rotation_example..."
cd misc/large_rotation_example || exit 1
sh run.sl .
RESULT=$?
cd ../..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: misc/large_rotation_example"
    FAILED_TESTS="$FAILED_TESTS misc/large_rotation_example"
else
    echo "PASSED: misc/large_rotation_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

echo "============================================"
if [ -z "$FAILED_TESTS" ]; then
    echo "  All examples completed successfully!"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "============================================"
    exit 0
else
    echo "  Some examples FAILED"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "  Failed tests:$FAILED_TESTS"
    echo "============================================"
    exit 1
fi
