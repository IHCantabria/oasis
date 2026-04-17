#!/bin/bash
# ==============================================================
# OASIS - Run all winches examples on Linux/cluster
# ==============================================================

echo "============================================"
echo "  Running winches examples"
echo "============================================"
echo

FAILED_TESTS=""
TOTAL_TESTS=0
PASSED_TESTS=0
TEST_NUM=0

for dir in constant_tension horizontal_control constant_tension_yaml; do
    TEST_NUM=$((TEST_NUM + 1))
    echo "[$TEST_NUM] Running winches/$dir..."
    cd "$dir" || exit 1
    sh run.sl .
    RESULT=$?
    cd ..
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    if [ $RESULT -ne 0 ]; then
        echo "FAILED: winches/$dir"
        FAILED_TESTS="$FAILED_TESTS winches/$dir"
    else
        echo "PASSED: winches/$dir"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
    echo
done

echo "============================================"
if [ -z "$FAILED_TESTS" ]; then
    echo "  All winches examples passed!"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "============================================"
    exit 0
else
    echo "  Some winches examples FAILED"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "  Failed tests:$FAILED_TESTS"
    echo "============================================"
    exit 1
fi
