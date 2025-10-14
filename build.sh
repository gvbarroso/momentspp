#!/usr/bin/env bash

#------------------------------------------------------------------------------
# Build script for moments++
# Author: Gustavo V. Barroso
# Created: 2025-09-24
#------------------------------------------------------------------------------

set -e  # Exit on error
set -u  # Treat unset variables as errors

#------------------------------------------------------------------------------
# Parse arguments
#------------------------------------------------------------------------------
BUILD_TYPE="Release"  # Default to Release
CLEAN_BUILD="OFF"

for arg in "$@"; do
  case "$arg" in
    Release|Debug)
      BUILD_TYPE="$arg"
      ;;
    --clean)
      CLEAN_BUILD="ON"
      ;;
    *)
      echo "❌ Unknown argument: $arg"
      echo "Usage: ./build.sh [Release|Debug] [--clean]"
      exit 1
      ;;
  esac
done

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
BUILD_DIR="build"
INSTALL_PREFIX="$HOME/.local"

# Optional flags
NATIVE_BUILD="ON"
DEBUG="OFF"
NAKED_D="OFF"

echo "🔧 Build type: $BUILD_TYPE"
echo "🧹 Clean build: $CLEAN_BUILD"

#------------------------------------------------------------------------------
# Step 1: Clean build directory if requested
#------------------------------------------------------------------------------
if [[ "$CLEAN_BUILD" == "ON" ]]; then
  echo "🧼 Removing existing build directory: $BUILD_DIR"
  rm -rf "$BUILD_DIR"
fi

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

#./build.sh                   # default to Release: optimized build
#./build.sh Debug             # debug symbols for Valgrind
#./build.sh Release --clean   # clean and rebuild Release
#./build.sh Debug --clean     # clean and rebuild debug
