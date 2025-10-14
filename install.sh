#!/bin/bash
set -e

echo "🚀 Starting moments++ installation..."

# Detect OS
OS="$(uname -s)"
echo "🧭 Detected OS: $OS"

# Step 1: Install dependencies
if [[ "$OS" == "Darwin" ]]; then
  echo "🍎 Installing dependencies via Homebrew..."
  brew install boost gsl mpfr gmp libomp yaml-cpp findutils cmake
elif [[ "$OS" == "Linux" ]]; then
  echo "🐧 Installing dependencies via apt..."
  sudo apt update
  sudo apt install -y \
    build-essential \
    cmake \
    libboost-iostreams-dev \
    libboost-random-dev \
    libboost-regex-dev \
    libgsl-dev \
    libmpfr-dev \
    libgmp-dev \
    libomp-dev \
    libyaml-cpp-dev \
    findutils
else
  echo "❌ Unsupported OS: $OS"
  exit 1
fi

# Step 2: Set environment variables
echo "🔧 Configuring environment..."
export PATH="$HOME/.local/bin:$PATH"
export DYLD_LIBRARY_PATH="$HOME/.local/lib:$DYLD_LIBRARY_PATH"  # macOS only
export LD_LIBRARY_PATH="$HOME/.local/lib:$LD_LIBRARY_PATH"      # Linux only
export LDFLAGS="-L$HOME/.local/lib"
export CPPFLAGS="-I$HOME/.local/include"

# Step 3: Build and install Bio++ core
echo "📦 Cloning and building bpp-core..."
git clone https://github.com/BioPP/bpp-core.git
cd bpp-core
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local -DBUILD_SHARED_LIBS=ON
make -j$(nproc || sysctl -n hw.logicalcpu)
make install
cd ../..

# Step 4: Build and install moments++ using build.sh
echo "🔨 Building moments++ via build.sh..."
./build.sh Release

echo "✅ Installation complete!"
echo "moments++ installed at: $HOME/.local/bin/momentspp"
echo "Run 'momentspp' from anywhere if your PATH includes ~/.local/bin"

# Optional: Refresh shell cache
hash -r

