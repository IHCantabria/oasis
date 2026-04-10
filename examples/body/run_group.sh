#!/bin/bash
# ==============================================================
# OASIS - Run all body examples on Linux/cluster
# ==============================================================

echo "============================================"
echo "  Running body examples"
echo "============================================"
echo

FAILED_TESTS=""
TOTAL_TESTS=0
PASSED_TESTS=0
TEST_NUM=0

for dir in free_decay fixed partial_dofs imposed_motion radiation excitation_linear excitation_instantpos nonlinear_hs_flat nonlinear_hs_waves qtf; do
    TEST_NUM=$((TEST_NUM + 1))
    echo "[$TEST_NUM] Running body/$dir..."
    cd "$dir" || exit 1
    sh run.sl .
    RESULT=$?
    cd ..
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    if [ $RESULT -ne 0 ]; then
        echo "FAILED: body/$dir"
        FAILED_TESTS="$FAILED_TESTS body/$dir"
    else
        echo "PASSED: body/$dir"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
    echo
done

echo "============================================"
if [ -z "$FAILED_TESTS" ]; then
    echo "  All body examples passed!"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "============================================"
    exit 0
else
    echo "  Some body examples FAILED"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "  Failed tests:$FAILED_TESTS"
    echo "============================================"
    exit 1
fi
