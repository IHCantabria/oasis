#!/bin/bash
# ==============================================================
# OASIS - Run all waves examples on Linux/cluster
# ==============================================================

echo "============================================"
echo "  Running waves examples"
echo "============================================"
echo

FAILED_TESTS=""
TOTAL_TESTS=0
PASSED_TESTS=0
TEST_NUM=0

for dir in regular jonswap jonswap_piecewise multidirectional timeseries frequency_domain; do
    TEST_NUM=$((TEST_NUM + 1))
    echo "[$TEST_NUM] Running waves/$dir..."
    cd "$dir" || exit 1
    sh run.sl .
    RESULT=$?
    cd ..
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    if [ $RESULT -ne 0 ]; then
        echo "FAILED: waves/$dir"
        FAILED_TESTS="$FAILED_TESTS waves/$dir"
    else
        echo "PASSED: waves/$dir"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
    echo
done

echo "============================================"
if [ -z "$FAILED_TESTS" ]; then
    echo "  All waves examples passed!"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "============================================"
    exit 0
else
    echo "  Some waves examples FAILED"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "  Failed tests:$FAILED_TESTS"
    echo "============================================"
    exit 1
fi
