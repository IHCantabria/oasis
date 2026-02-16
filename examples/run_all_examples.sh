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

# --- freq_wave_example ---
echo "[1/4] Running freq_wave_example..."
cd freq_wave_example || exit 1
sh run.sl .
RESULT=$?
cd ..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: freq_wave_example"
    FAILED_TESTS="$FAILED_TESTS freq_wave_example"
else
    echo "PASSED: freq_wave_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- generic_example ---
echo "[2/4] Running generic_example..."
cd generic_example || exit 1
sh run.sl .
RESULT=$?
cd ..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: generic_example"
    FAILED_TESTS="$FAILED_TESTS generic_example"
else
    echo "PASSED: generic_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- qtf_example ---
echo "[3/4] Running qtf_example..."
cd qtf_example || exit 1
sh run.sl .
RESULT=$?
cd ..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: qtf_example"
    FAILED_TESTS="$FAILED_TESTS qtf_example"
else
    echo "PASSED: qtf_example"
    PASSED_TESTS=$((PASSED_TESTS + 1))
fi
echo

# --- turbine_example ---
echo "[4/4] Running turbine_example..."
cd turbine_example || exit 1
sh run.sl .
RESULT=$?
cd ..
TOTAL_TESTS=$((TOTAL_TESTS + 1))
if [ $RESULT -ne 0 ]; then
    echo "FAILED: turbine_example"
    FAILED_TESTS="$FAILED_TESTS turbine_example"
else
    echo "PASSED: turbine_example"
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
