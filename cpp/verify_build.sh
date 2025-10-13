#!/bin/bash
# Quick verification script for C++ build

echo "========================================================================"
echo "Richerme Ion Analog C++ - Build Verification"
echo "========================================================================"
echo ""

# Navigate to cpp directory
cd "$(dirname "$0")"

# Check for Eigen3
echo "Step 1: Checking for Eigen3..."
if [ -d "/opt/homebrew/include/eigen3" ] || [ -d "/usr/local/include/eigen3" ]; then
    echo "  ✓ Eigen3 found"
else
    echo "  ✗ Eigen3 not found - run: brew install eigen"
    exit 1
fi
echo ""

# Create build directory
echo "Step 2: Setting up build directory..."
rm -rf build_test
mkdir -p build_test
echo "  ✓ Build directory created"
echo ""

# Run CMake
echo "Step 3: Running CMake configuration..."
cd build_test
cmake .. > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "  ✓ CMake configuration successful"
else
    echo "  ✗ CMake configuration failed"
    cmake ..
    exit 1
fi
echo ""

# Build
echo "Step 4: Compiling library and tests..."
make -j$(sysctl -n hw.ncpu) > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "  ✓ Compilation successful"
else
    echo "  ✗ Compilation failed"
    make
    exit 1
fi
echo ""

# Run tests
echo "Step 5: Running test suite..."
echo "------------------------------------------------------------------------"
./test_richerme
TEST_RESULT=$?
echo "------------------------------------------------------------------------"
echo ""

if [ $TEST_RESULT -eq 0 ]; then
    echo "Step 6: Running example..."
    echo "------------------------------------------------------------------------"
    ./example_basic
    echo "------------------------------------------------------------------------"
    echo ""
fi

# Summary
echo "========================================================================"
if [ $TEST_RESULT -eq 0 ]; then
    echo "SUCCESS: All tests passed!"
    echo "========================================================================"
    echo ""
    echo "You can now use the library:"
    echo "  - Tests: ./build_test/test_richerme"
    echo "  - Example: ./build_test/example_basic"
    echo "  - Link to your code: -lricherme_ion_analog"
    echo ""
else
    echo "FAILURE: Tests failed (see output above)"
    echo "========================================================================"
    exit 1
fi
