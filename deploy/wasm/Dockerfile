# Dockerfile for Thermal-FIST WASM build and deployment
# Multi-stage build: Stage 1 builds, Stage 2 serves
#
# Usage: docker build -t thermal-fist-wasm .
#        docker run -p 8081:80 thermal-fist-wasm
#        Then open http://localhost:8081 in your browser

# =============================================================================
# Stage 1: Build environment
# =============================================================================
FROM ubuntu:22.04 AS builder

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    python3 \
    python3-pip \
    wget \
    curl \
    xz-utils \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Install aqtinstall for downloading Qt
RUN pip3 install aqtinstall

# Install Emscripten 4.0.7 (matching local dev environment for Qt 6.10.x)
WORKDIR /opt
RUN git clone https://github.com/emscripten-core/emsdk.git
WORKDIR /opt/emsdk
RUN ./emsdk install 4.0.7 && \
    ./emsdk activate 4.0.7

# Install Qt 6.10.1 for WebAssembly (Qt 6.7+ uses different aqtinstall syntax: all_os wasm)
# Also install host Qt (linux_gcc_64) which is required for cross-compilation tools
WORKDIR /opt
RUN aqt install-qt linux desktop 6.10.1 linux_gcc_64 && \
    aqt install-qt all_os wasm 6.10.1 wasm_multithread

# Set up environment variables
ENV EMSDK=/opt/emsdk
ENV EM_CONFIG=/opt/emsdk/.emscripten
ENV QT_DIR=/opt/6.10.1
ENV QT_HOST_PATH=/opt/6.10.1/gcc_64

# Clone Thermal-FIST and run initial cmake configure (cached)
WORKDIR /opt
RUN git clone --depth 1 -b devel https://github.com/vlvovch/Thermal-FIST.git
WORKDIR /opt/Thermal-FIST/build-wasm
RUN cd /opt/emsdk && . ./emsdk_env.sh && cd /opt/Thermal-FIST/build-wasm && \
    /opt/6.10.1/wasm_multithread/bin/qt-cmake .. \
        -DTHERMALFIST_WASM=ON \
        -DQT_HOST_PATH=/opt/6.10.1/gcc_64

# Pull latest changes (use --build-arg CACHEBUST=$(date +%s) to force refresh)
ARG CACHEBUST=1
RUN cd /opt/Thermal-FIST && git pull

# Incremental build - uses BuildKit cache mount to persist object files
RUN --mount=type=cache,target=/opt/Thermal-FIST/build-wasm \
    cd /opt/emsdk && . ./emsdk_env.sh && cd /opt/Thermal-FIST/build-wasm && \
    /opt/6.10.1/wasm_multithread/bin/qt-cmake .. \
        -DTHERMALFIST_WASM=ON \
        -DQT_HOST_PATH=/opt/6.10.1/gcc_64 && \
    cmake --build . --parallel $(nproc) && \
    cp -r bin /opt/Thermal-FIST/build-output

# =============================================================================
# Stage 2: Serve with Caddy
# =============================================================================
FROM caddy:alpine

# Copy built WASM files (from build-output since cache mount isn't visible to COPY)
COPY --from=builder /opt/Thermal-FIST/build-output/ /srv/

# Copy Caddyfile
COPY Caddyfile /etc/caddy/Caddyfile

# Expose HTTP port
EXPOSE 80

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD wget -q --spider http://localhost/ || exit 1
