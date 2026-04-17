#!/bin/bash
# ==============================================================
# OASIS - Run all solvers examples on Linux/cluster
# ==============================================================

echo "============================================"
echo "  Running solvers examples"
echo "============================================"
echo

FAILED_TESTS=""
TOTAL_TESTS=0
PASSED_TESTS=0
TEST_NUM=0

for dir in bdf1 bdfn_order2 bdfn_order4 esdirk46 bdf1_yaml; do
    TEST_NUM=$((TEST_NUM + 1))
    echo "[$TEST_NUM] Running solvers/$dir..."
    cd "$dir" || exit 1
    sh run.sl .
    RESULT=$?
    cd ..
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    if [ $RESULT -ne 0 ]; then
        echo "FAILED: solvers/$dir"
        FAILED_TESTS="$FAILED_TESTS solvers/$dir"
    else
        echo "PASSED: solvers/$dir"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
    echo
done

echo "============================================"
if [ -z "$FAILED_TESTS" ]; then
    echo "  All solvers examples passed!"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "============================================"
    exit 0
else
    echo "  Some solvers examples FAILED"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "  Failed tests:$FAILED_TESTS"
    echo "============================================"
    exit 1
fi
