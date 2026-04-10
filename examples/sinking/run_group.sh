#!/bin/bash
# ==============================================================
# OASIS - Run all sinking examples
# ==============================================================

echo "============================================"
echo "  Running sinking examples"
echo "============================================"
echo

FAILED=0

echo "Running progressive_flooding..."
cd progressive_flooding || exit 1
sh run.sl .
if [ $? -ne 0 ]; then FAILED=1; fi
cd ..

echo
if [ $FAILED -ne 0 ]; then
    echo "Some sinking examples FAILED"
    exit 1
else
    echo "All sinking examples PASSED"
    exit 0
fi
