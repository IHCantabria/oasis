#!/bin/bash
# ==============================================================
# OASIS - Run all morison examples
# ==============================================================

echo "============================================"
echo "  Running morison examples"
echo "============================================"
echo

FAILED=0

echo "Running wind_drag..."
cd wind_drag || exit 1
sh run.sl .
if [ $? -ne 0 ]; then FAILED=1; fi
cd ..

echo
echo "Running current_drag..."
cd current_drag || exit 1
sh run.sl .
if [ $? -ne 0 ]; then FAILED=1; fi
cd ..

echo
echo "Running wind_drag_yaml..."
cd wind_drag_yaml || exit 1
sh run.sl .
if [ $? -ne 0 ]; then FAILED=1; fi
cd ..

echo
echo "Running current_drag_yaml..."
cd current_drag_yaml || exit 1
sh run.sl .
if [ $? -ne 0 ]; then FAILED=1; fi
cd ..

echo
if [ $FAILED -ne 0 ]; then
    echo "Some morison examples FAILED"
    exit 1
else
    echo "All morison examples PASSED"
    exit 0
fi
