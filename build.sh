#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Build script for moments++
# Author: Gustavo V. Barroso
# Created: 2025-09-24
#------------------------------------------------------------------------------

set -e  # Exit on error
set -u  # Treat unset variables as errors

#------------------------------------------------------------------------------
# Parse build type from command line
#------------------------------------------------------------------------------
BUILD_TYPE="${1:-Debug}"  # Default to Debug if not provided

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
BUILD_DIR="build"
INSTALL_PREFIX="$HOME/.local"

# Optional flags
NATIVE_BUILD="ON"
DEBUG="OFF" # turning ON/OFF compulation of blocks: ifdef DEBUG
NAKED_D="OFF"

echo "🔧 Build type: $BUILD_TYPE"

#------------------------------------------------------------------------------
# Step 1: Create build directory
#------------------------------------------------------------------------------
echo "📁 Creating build directory: $BUILD_DIR"
mkdir -p "$BUILD_DIR"

#------------------------------------------------------------------------------
# Step 2: Configure with CMake
#------------------------------------------------------------------------------
echo "⚙️ Configuring project with CMake..."
cmake -S . -B "$BUILD_DIR" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
  -DNativeBuild="$NATIVE_BUILD" \
  -DDEBUG="$DEBUG" \
  -DNAKED_D="$NAKED_D"

#------------------------------------------------------------------------------
# Step 3: Build with parallel jobs
#------------------------------------------------------------------------------
echo "🔨 Building project..."
NUM_CORES=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)
cmake --build "$BUILD_DIR" -- -j"$NUM_CORES"

#------------------------------------------------------------------------------
# Step 4: Install to user prefix
#------------------------------------------------------------------------------
echo "📦 Installing to $INSTALL_PREFIX..."
cmake --install "$BUILD_DIR"

echo "✅ Build and install complete."

#bash build.sh Release   # optimized build
#bash build.sh Debug     # debug symbols for Valgrind
#bash build.sh           # defaults to Debug

