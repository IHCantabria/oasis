#!/bin/bash
# ==============================================================
# OASIS - Run all lines examples on Linux/cluster
# ==============================================================

echo "============================================"
echo "  Running lines examples"
echo "============================================"
echo

FAILED_TESTS=""
TOTAL_TESTS=0
PASSED_TESTS=0
TEST_NUM=0

for dir in single_line multi_line joint_connection elastic_anchor prescribed_motion viscoelastic tabulated_stiffness tension_symmetric seabed_contact friction_isotropic friction_anisotropic; do
    TEST_NUM=$((TEST_NUM + 1))
    echo "[$TEST_NUM] Running lines/$dir..."
    cd "$dir" || exit 1
    sh run.sl .
    RESULT=$?
    cd ..
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    if [ $RESULT -ne 0 ]; then
        echo "FAILED: lines/$dir"
        FAILED_TESTS="$FAILED_TESTS lines/$dir"
    else
        echo "PASSED: lines/$dir"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
    echo
done

echo "============================================"
if [ -z "$FAILED_TESTS" ]; then
    echo "  All lines examples passed!"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "============================================"
    exit 0
else
    echo "  Some lines examples FAILED"
    echo "  $PASSED_TESTS/$TOTAL_TESTS tests passed"
    echo "  Failed tests:$FAILED_TESTS"
    echo "============================================"
    exit 1
fi
