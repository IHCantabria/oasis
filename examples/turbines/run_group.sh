#!/bin/bash
# ==============================================================
# OASIS - Run all turbines examples on Linux/cluster
# ==============================================================

echo "============================================"
echo "  Running turbines examples"
echo "============================================"
echo

FAILED_TESTS=""
TOTAL_TESTS=0
PASSED_TESTS=0
TEST_NUM=0

for dir in fixed_turbine moored_turbine; do
    TEST_NUM=$((TEST_NUM + 1))
    echo "[$TEST_NUM] Running turbines/$dir..."
    cd "$dir" || exit 1
    sh run.sl .
    RESULT=$?
    cd ..
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    if [ $RESULT -ne 0 ]; then
        echo "FAILED: turbines/$dir"
        FAILED_TESTS="$FAILED_TESTS turbines/$dir"
    else
        echo "PASSED: turbines/$dir"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
    echo
done

echo "============================================"
if [ -z "$FAILED_TESTS" ]; then
    echo "  All turbines examples passed!"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "============================================"
    exit 0
else
    echo "  Some turbines examples FAILED"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "  Failed tests:$FAILED_TESTS"
    echo "============================================"
    exit 1
fi
