#!/bin/bash

# Clean build directory
rm -rf build
mkdir -p build
cd build

# Run CMake configuration
cmake ..

# Build
make -j$(sysctl -n hw.ncpu)

# Run tests
ctest --verbose

echo "Build complete!"
