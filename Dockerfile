# ───────────────────────────────────────────────────────────────
# Stage 1: Build momentspp and dependencies
# ───────────────────────────────────────────────────────────────
FROM debian:bookworm-slim AS build

ARG NCORES=1
ARG EIGEN_VERSION=3.4.0
ARG YAML_CPP_VERSION=0.8.0
ARG BPP_CORE_COMMIT=9f8d3e2300afb1d4a9e06c0a99ac3c7aeba03581

ENV DEBIAN_FRONTEND=noninteractive
ENV CMAKE_ARGS="-DCMAKE_INSTALL_PREFIX=/usr/local \
                -DCMAKE_PREFIX_PATH=/usr/local \
                -DCMAKE_INSTALL_RPATH=/usr/local \
                -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE"

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    curl \
    ca-certificates \
    libboost-iostreams-dev \
    libgsl-dev \
    libmpfr-dev \
    libgmp-dev \
    && rm -rf /var/lib/apt/lists/*

# ───────────────────────────────────────────────────────────────
# Install Eigen from official release tarball
RUN curl -fsSL -o /tmp/eigen.tar.gz \
      https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz && \
    mkdir -p /tmp/eigen-src && \
    tar -xzf /tmp/eigen.tar.gz -C /tmp/eigen-src --strip-components=1 && \
    cmake -B /tmp/eigen-src/build -S /tmp/eigen-src ${CMAKE_ARGS} && \
    cmake --build /tmp/eigen-src/build -j${NCORES} && \
    cmake --install /tmp/eigen-src/build && \
    rm -rf /tmp/eigen*

# ───────────────────────────────────────────────────────────────
# Install Bio++ core3
RUN git clone https://github.com/BioPP/bpp-core.git /tmp/bpp-core && \
    cd /tmp/bpp-core && git checkout ${BPP_CORE_COMMIT} && \
    cmake -B build -S . ${CMAKE_ARGS} && \
    cmake --build build -j${NCORES} && \
    cmake --install build && \
    rm -rf /tmp/bpp-core

# ───────────────────────────────────────────────────────────────
# Install yaml-cpp
RUN git clone https://github.com/jbeder/yaml-cpp.git /tmp/yaml-cpp && \
    cd /tmp/yaml-cpp && git checkout ${YAML_CPP_VERSION} && \
    cmake -B build -S . ${CMAKE_ARGS} && \
    cmake --build build -j${NCORES} && \
    cmake --install build && \
    rm -rf /tmp/yaml-cpp

# ───────────────────────────────────────────────────────────────
# Build momentspp
WORKDIR /opt/momentspp
COPY CMakeLists.txt .
COPY src/ src/

RUN cmake -B build -S . ${CMAKE_ARGS} && \
    cmake --build build -j${NCORES} && \
    cmake --install build

# Optional: remove static libs to reduce image size
RUN find /usr/local/lib -type f -name '*.a' -delete

# ───────────────────────────────────────────────────────────────
# Stage 2: Final runtime image
# ───────────────────────────────────────────────────────────────
FROM debian:bookworm-slim AS final

COPY --from=build /usr/local/lib /usr/local/lib
COPY --from=build /usr/local/bin /usr/local/bin

RUN apt-get update && apt-get install -y --no-install-recommends \
    libboost-iostreams1.74.0 \
    libgomp1 \
    libgsl27 \
    libmpfr6 \
    libgmp10 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["momentspp"]

