#!/bin/bash
# ==============================================================
# OASIS - Run all owcs examples on Linux/cluster
# ==============================================================

echo "============================================"
echo "  Running owcs examples"
echo "============================================"
echo

FAILED_TESTS=""
TOTAL_TESTS=0
PASSED_TESTS=0
TEST_NUM=0

for dir in free_decay excitation hole_coupling turbine_coupling; do
    TEST_NUM=$((TEST_NUM + 1))
    echo "[$TEST_NUM] Running owcs/$dir..."
    cd "$dir" || exit 1
    sh run.sl .
    RESULT=$?
    cd ..
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    if [ $RESULT -ne 0 ]; then
        echo "FAILED: owcs/$dir"
        FAILED_TESTS="$FAILED_TESTS owcs/$dir"
    else
        echo "PASSED: owcs/$dir"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
    echo
done

echo "============================================"
if [ -z "$FAILED_TESTS" ]; then
    echo "  All owcs examples passed!"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "============================================"
    exit 0
else
    echo "  Some owcs examples FAILED"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "  Failed tests:$FAILED_TESTS"
    echo "============================================"
    exit 1
fi
