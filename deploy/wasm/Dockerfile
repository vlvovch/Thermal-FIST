# Dockerfile for Thermal-FIST WASM build and deployment
# Uses pre-built base image with Emscripten + Qt (see Dockerfile.base)
#
# First time setup:
#   docker build -f Dockerfile.base -t thermal-fist-wasm-base .
#
# Build & deploy (fast — only rebuilds Thermal-FIST):
#   docker compose build
#   docker compose up -d
#
# Update (pulls latest code and rebuilds):
#   ./update_docker.sh

# =============================================================================
# Stage 1: Build Thermal-FIST
# =============================================================================
FROM thermal-fist-wasm-base AS builder

# Clone and build Thermal-FIST
# ARG CACHEBUST forces git clone to re-run (bypasses Docker layer cache)
ARG CACHEBUST=1
WORKDIR /opt
RUN git clone --depth 1 -b devel https://github.com/vlvovch/Thermal-FIST.git

WORKDIR /opt/Thermal-FIST/build-wasm
RUN cd /opt/emsdk && . ./emsdk_env.sh && cd /opt/Thermal-FIST/build-wasm && \
    /opt/6.10.1/wasm_multithread/bin/qt-cmake .. \
        -DTHERMALFIST_WASM=ON \
        -DQT_HOST_PATH=/opt/6.10.1/gcc_64 && \
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
