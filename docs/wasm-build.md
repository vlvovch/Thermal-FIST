# Building Thermal-FIST for WebAssembly

QtThermalFIST can be compiled to WebAssembly (WASM) to run directly in a web browser without installation.

## Prerequisites

### 1. Emscripten Toolchain

Install Emscripten by following the [official guide](https://emscripten.org/docs/getting_started/downloads.html):

```bash
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
```

### 2. Qt 6 for WebAssembly

Use the [Qt Online Installer](https://www.qt.io/download-qt-installer) and select the WebAssembly components:
- `Qt 6.x.x` → `WebAssembly (single-threaded)` and/or `WebAssembly (multi-threaded)`

The multi-threaded version provides better performance for calculations but requires additional server configuration (see below).

## Building

### Single-threaded Build

```bash
# Source the Emscripten environment
source /path/to/emsdk/emsdk_env.sh

# Configure with qt-cmake
mkdir build-wasm
cd build-wasm
/path/to/Qt/6.x.x/wasm_singlethread/bin/qt-cmake ../ -DTHERMALFIST_WASM=ON

# Build
cmake --build . --parallel
```

### Multi-threaded Build (Recommended)

For better performance during calculations:

```bash
source /path/to/emsdk/emsdk_env.sh

mkdir build-wasm-mt
cd build-wasm-mt
/path/to/Qt/6.x.x/wasm_multithread/bin/qt-cmake ../ -DTHERMALFIST_WASM=ON

cmake --build . --parallel
```

## Running Locally

The build produces files in `build-wasm/bin/`:
- `QtThermalFIST.html` - Main entry point
- `QtThermalFIST.js` - JavaScript glue code
- `QtThermalFIST.wasm` - WebAssembly binary
- `qtloader.js` - Qt's loader script

### Single-threaded Build

Any HTTP server works:

```bash
cd build-wasm/bin
python3 -m http.server 8080
```

Open `http://localhost:8080/QtThermalFIST.html` in your browser.

### Multi-threaded Build

Multi-threaded WASM requires `SharedArrayBuffer`, which needs specific HTTP headers. Use `emrun`:

```bash
cd build-wasm-mt/bin
emrun --port 8080 QtThermalFIST.html
```

Or configure your web server to send these headers:
```
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

## Browser Compatibility

| Browser | Single-threaded | Multi-threaded |
|---------|-----------------|----------------|
| Chrome/Edge | Yes | Yes |
| Firefox | Yes | Yes |
| Safari | Yes | Limited* |

*Safari has unreliable WASM threading support. The application automatically falls back to synchronous execution on Safari.

## WASM-specific Behavior

### File I/O
- **Loading files**: Uses browser file picker dialog
- **Saving files**: Triggers browser download
- No direct filesystem access (sandboxed environment)

### Threading
- Multi-threaded builds run calculations in background threads
- Single-threaded builds block the UI during calculations
- Safari always uses synchronous execution regardless of build type

### Embedded Resources
The WASM build includes embedded data files:
- Default particle list (PDG2020 with nuclei)
- Default experimental data (ALICE Pb-Pb 2.76 TeV)
- Lattice QCD reference data (Wuppertal-Budapest and HotQCD)

## Deployment

To deploy on a web server:

1. Copy the contents of `build-wasm/bin/` to your web server
2. Ensure the server sends correct MIME types:
   - `.wasm` → `application/wasm`
   - `.js` → `application/javascript`
3. For multi-threaded builds, configure COOP/COEP headers (see above)

## Troubleshooting

### "SharedArrayBuffer is not defined"
- You're running a multi-threaded build without proper headers
- Use `emrun` or configure COOP/COEP headers on your server

### Calculations freeze the UI
- This is expected for single-threaded builds or Safari
- Use a multi-threaded build with Chrome/Firefox for background calculations

### Build fails with "Cannot find Qt"
- Ensure you're using `qt-cmake` from the WASM kit, not the native Qt installation
- Verify Emscripten environment is sourced before running cmake

### Large download size
- WASM builds are typically 20-40 MB
- Consider enabling gzip/brotli compression on your web server
