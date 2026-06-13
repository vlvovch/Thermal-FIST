#!/bin/bash
# Update Thermal-FIST WASM deployment
#
# First time: ./update_docker.sh --init  (builds base image + app, ~30 min)
# Updates:    ./update_docker.sh          (rebuilds app only, ~5 min)

set -e

cd "$(dirname "$0")"

if [ "$1" = "--init" ]; then
    echo "=== Building base image (Emscripten + Qt) ==="
    echo "This takes ~20-30 minutes but only needs to be done once."
    docker build -f Dockerfile.base -t thermal-fist-wasm-base .
    echo ""
fi

# Check if base image exists
if ! docker image inspect thermal-fist-wasm-base >/dev/null 2>&1; then
    echo "Error: Base image not found. Run with --init first:"
    echo "  ./update_docker.sh --init"
    exit 1
fi

echo "=== Building Thermal-FIST WASM ==="
# Build new image while old container keeps serving traffic
docker compose build --build-arg CACHEBUST=$(date +%s)
# Recreate container with new image (only seconds of downtime)
docker compose up -d
echo ""
echo "=== Done ==="
docker compose ps
