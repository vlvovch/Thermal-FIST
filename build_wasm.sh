#!/bin/bash
# Build Thermal-FIST for WebAssembly locally
#
# Usage:
#   ./build_wasm.sh              # multi-threaded build (default)
#   ./build_wasm.sh st           # single-threaded build
#   ./build_wasm.sh mt           # multi-threaded build
#   ./build_wasm.sh serve        # build multi-threaded and start local server
#   ./build_wasm.sh serve-st     # build single-threaded and start local server
#
# Prerequisites:
#   - Emscripten SDK
#   - Qt 6 with WebAssembly components

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ---- Configuration (edit these paths if needed) ----
EMSDK="${EMSDK:-/Users/vlvovch/Code/tools/emsdk}"
QT_VERSION="${QT_VERSION:-6.10.1}"
QT_DIR="${QT_DIR:-/Users/vlvovch/Qt/${QT_VERSION}}"
# ----------------------------------------------------

MODE="${1:-mt}"

case "$MODE" in
  st)          THREADING="singlethread" ; SERVE=false ;;
  mt)          THREADING="multithread"  ; SERVE=false ;;
  serve)       THREADING="multithread"  ; SERVE=true  ;;
  serve-st)    THREADING="singlethread" ; SERVE=true  ;;
  *)
    echo "Usage: $0 [st|mt|serve|serve-st]"
    exit 1
    ;;
esac

QT_WASM="${QT_DIR}/wasm_${THREADING}"
QT_HOST="${QT_DIR}/macos"
BUILD_DIR="${SCRIPT_DIR}/build-wasm-${THREADING}"

# Validate paths
if [ ! -f "${EMSDK}/emsdk_env.sh" ]; then
  echo "Error: emsdk not found at ${EMSDK}"
  echo "Set EMSDK environment variable to your emsdk directory"
  exit 1
fi

if [ ! -f "${QT_WASM}/bin/qt-cmake" ]; then
  echo "Error: Qt WASM kit not found at ${QT_WASM}"
  echo "Install Qt ${QT_VERSION} WebAssembly (${THREADING}) via Qt Online Installer"
  exit 1
fi

if [ ! -d "${QT_HOST}" ]; then
  echo "Error: Qt host tools not found at ${QT_HOST}"
  echo "Install Qt ${QT_VERSION} desktop kit (macOS) via Qt Online Installer"
  exit 1
fi

echo "=== Thermal-FIST WASM Build ==="
echo "  Threading: ${THREADING}"
echo "  Emscripten: ${EMSDK}"
echo "  Qt WASM:    ${QT_WASM}"
echo "  Qt Host:    ${QT_HOST}"
echo "  Build dir:  ${BUILD_DIR}"
echo ""

# Source Emscripten environment
source "${EMSDK}/emsdk_env.sh"

# Create build directory
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Configure
echo "=== Configuring ==="
"${QT_WASM}/bin/qt-cmake" "${SCRIPT_DIR}" \
    -DTHERMALFIST_WASM=ON \
    -DQT_HOST_PATH="${QT_HOST}"

# Build
echo ""
echo "=== Building ==="
cmake --build . --parallel $(sysctl -n hw.ncpu 2>/dev/null || nproc)

echo ""
echo "=== Build complete ==="
echo "Output: ${BUILD_DIR}/bin/"
ls -lh "${BUILD_DIR}/bin/"*.wasm "${BUILD_DIR}/bin/"*.html 2>/dev/null

# Serve if requested
if $SERVE; then
  echo ""
  echo "=== Starting local server ==="
  cd "${BUILD_DIR}/bin"
  if [ "$THREADING" = "multithread" ]; then
    echo "Multi-threaded: using emrun (provides required COOP/COEP headers)"
    echo "Open http://localhost:8080/QtThermalFIST.html"
    emrun --port 8080 QtThermalFIST.html
  else
    echo "Single-threaded: using Python HTTP server"
    echo "Open http://localhost:8080/QtThermalFIST.html"
    python3 -m http.server 8080
  fi
fi
