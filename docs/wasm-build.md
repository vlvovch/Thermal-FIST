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

### Quick Deployment (Single-threaded)

1. Copy the contents of `build-wasm/bin/` to your web server
2. Ensure the server sends correct MIME types:
   - `.wasm` → `application/wasm`
   - `.js` → `application/javascript`

### Caddy Deployment (Recommended)

Caddy is simple and provides automatic HTTPS via Let's Encrypt.

**1. Install Caddy:**

```bash
# Debian/Ubuntu
sudo apt install -y debian-keyring debian-archive-keyring apt-transport-https curl
curl -1sLf 'https://dl.cloudsmith.io/public/caddy/stable/gpg.key' | sudo gpg --dearmor -o /usr/share/keyrings/caddy-stable-archive-keyring.gpg
curl -1sLf 'https://dl.cloudsmith.io/public/caddy/stable/debian.deb.txt' | sudo tee /etc/apt/sources.list.d/caddy-stable.list
sudo apt update
sudo apt install caddy
```

**2. Copy files:**

```bash
sudo mkdir -p /var/www/thermal-fist
sudo cp -r build-wasm/bin/* /var/www/thermal-fist/
```

**3. Create `/etc/caddy/Caddyfile`:**

```caddyfile
your-domain.com {
    root * /var/www/thermal-fist
    file_server

    # Required headers for multi-threaded WASM (SharedArrayBuffer)
    header {
        Cross-Origin-Opener-Policy "same-origin"
        Cross-Origin-Embedder-Policy "require-corp"
    }

    # Serve QtThermalFIST.html as index
    try_files {path} /QtThermalFIST.html
}
```

For local testing without a domain (HTTP only):

```caddyfile
:8080 {
    root * /var/www/thermal-fist
    file_server

    header {
        Cross-Origin-Opener-Policy "same-origin"
        Cross-Origin-Embedder-Policy "require-corp"
    }

    try_files {path} /QtThermalFIST.html
}
```

**4. Start Caddy:**

```bash
sudo systemctl restart caddy
```

Caddy automatically obtains and renews HTTPS certificates via Let's Encrypt.

### Automated Docker Build (Recommended)

A complete Docker setup that builds everything from source is available in `deploy/wasm/`:

```bash
cd deploy/wasm
docker compose build   # Takes 15-30 minutes
docker compose up -d
```

This automatically:
- Installs Emscripten and Qt for WASM
- Clones and builds Thermal-FIST
- Serves via Caddy with proper headers

See [deploy/wasm/README.md](../deploy/wasm/README.md) for details.

### Docker with Pre-built Files

If you've already built the WASM files locally:

**1. Create deployment directory:**

```bash
mkdir thermal-fist-web
cd thermal-fist-web
cp -r /path/to/build-wasm/bin/* .
```

**2. Create `Caddyfile`:**

```caddyfile
:80 {
    root * /srv
    file_server

    header {
        Cross-Origin-Opener-Policy "same-origin"
        Cross-Origin-Embedder-Policy "require-corp"
    }

    try_files {path} /QtThermalFIST.html
}
```

**3. Create `Dockerfile`:**

```dockerfile
FROM caddy:alpine
COPY Caddyfile /etc/caddy/Caddyfile
COPY . /srv/
RUN rm -f /srv/Caddyfile /srv/Dockerfile
```

**4. Build and run:**

```bash
docker build -t thermal-fist-web .
docker run -d -p 8080:80 --name thermal-fist thermal-fist-web
```

Access at `http://your-server:8080`

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
