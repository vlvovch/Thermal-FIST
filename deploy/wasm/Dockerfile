# Dockerfile for Thermal-FIST WASM build and deployment
# Multi-stage build: Stage 1 builds, Stage 2 serves
#
# Usage: docker build -t thermal-fist-wasm .
#        docker run -p 8080:80 thermal-fist-wasm
#        Then open http://localhost:8080 in your browser

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
    && rm -rf /var/lib/apt/lists/*

# Install aqtinstall for downloading Qt
RUN pip3 install aqtinstall

# Install Emscripten (Qt 6.6 requires Emscripten 3.1.37)
WORKDIR /opt
RUN git clone https://github.com/emscripten-core/emsdk.git
WORKDIR /opt/emsdk
RUN ./emsdk install 3.1.37 && \
    ./emsdk activate 3.1.37

# Install Qt 6.6 for WebAssembly (multi-threaded for better performance)
WORKDIR /opt
RUN aqt install-qt linux desktop 6.6.0 wasm_multithread

# Set up environment
ENV EMSDK=/opt/emsdk
ENV QT_DIR=/opt/6.6.0

# Clone Thermal-FIST
WORKDIR /opt
RUN git clone --depth 1 -b feature/wasm-support https://github.com/vlvovch/Thermal-FIST.git

# Build for WASM
WORKDIR /opt/Thermal-FIST/build-wasm
RUN . /opt/emsdk/emsdk_env.sh && \
    /opt/6.6.0/wasm_multithread/bin/qt-cmake .. -DTHERMALFIST_WASM=ON && \
    cmake --build . --parallel $(nproc)

# =============================================================================
# Stage 2: Serve with Caddy
# =============================================================================
FROM caddy:alpine

# Copy built WASM files
COPY --from=builder /opt/Thermal-FIST/build-wasm/bin/ /srv/

# Copy Caddyfile
COPY Caddyfile /etc/caddy/Caddyfile

# Expose HTTP port
EXPOSE 80

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD wget -q --spider http://localhost/ || exit 1
