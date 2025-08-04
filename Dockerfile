FROM debian:bookworm-slim AS build

ARG NCORES=1

ARG EIGEN_VERSION=3.4
ARG BPP_CORE_VERSION=9f8d3e2300afb1d4a9e06c0a99ac3c7aeba03581
ARG BPP_SEQ_VERSION=fe03cf8c41508bf31199e052116e31d2e4beca02
ARG BPP_PHYL_VERSION=1b16ca9028f7a06fb15fc64ae82d5c57d970be64
ARG YAML_CPP_VERSION=0.8.0
ARG BPP_MPP_CMAKE_ARGS="\
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DCMAKE_PREFIX_PATH=/usr/local \
    -DCMAKE_INSTALL_RPATH=/usr/local \
    -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE \
    "

RUN \
    set -ex; \
    apt-get update; \
    apt-get install -y build-essential cmake git libboost-iostreams-dev libgsl-dev; \
    rm -rf /var/lib/apt/lists/*

RUN set -ex; \
    git clone https://gitlab.com/libeigen/eigen.git; \
    cd eigen; git checkout "${EIGEN_VERSION}"; \
    mkdir build; cd build; cmake ..; make -j "${NCORES}"; make install; \
    cd /; rm -rf eigen

RUN set -ex; \
    git clone https://github.com/BioPP/bpp-core.git; \
    cd bpp-core; git checkout "${BPP_CORE_VERSION}"; \
    mkdir build; cd build; cmake ${BPP_MPP_CMAKE_ARGS} ..; \
    make -j "${NCORES}"; make install; \
    cd /; rm -rf bpp-core

# Depends on bpp-core
RUN set -ex; \
    git clone https://github.com/BioPP/bpp-seq.git; \
    cd bpp-seq; git checkout "${BPP_SEQ_VERSION}"; \
    mkdir build; cd build; cmake ${BPP_MPP_CMAKE_ARGS} ..; \
    make -j "${NCORES}"; make install; \
    cd /; rm -rf bpp-seq

# Depends on bpp-core, bpp-seq, and eigen
RUN set -ex; \
    git clone https://github.com/BioPP/bpp-phyl.git; \
    cd bpp-phyl; git checkout "${BPP_PHYL_VERSION}"; \
    mkdir build; cd build; cmake ${BPP_MPP_CMAKE_ARGS} ..; \
    make -j "${NCORES}"; make install; \
    cd /; rm -rf bpp-phyl

# Depends on bpp-core, bpp-seq, and eigen
RUN set -ex; \
    git clone https://github.com/jbeder/yaml-cpp.git; \
    cd yaml-cpp; git checkout "${YAML_CPP_VERSION}"; \
    mkdir build; cd build; cmake ..; make -j "${NCORES}"; make install; \
    cd /; rm -rf yaml-cpp

# Change to the build dir
WORKDIR /opt/momentspp

# Add the necessary files for building momentspp
ADD src src
ADD CMakeLists.txt .

# Build momentspp
RUN set -ex; \
    mkdir build; cd build; cmake ${BPP_MPP_CMAKE_ARGS} ..; \
    make -j "${NCORES}"; make install

# Remove static libs before creating final build
RUN find /usr/local/lib -type f -name '*.a' -delete


FROM debian:bookworm-slim AS final
COPY --from=build /usr/local/lib /usr/local/lib
COPY --from=build /usr/local/bin /usr/local/bin

RUN \
    set -ex; \
    apt-get update; \
    apt-get install -y --no-install-recommends libboost-iostreams1.74.0 libgomp1 libgsl27; \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["momentspp"]
