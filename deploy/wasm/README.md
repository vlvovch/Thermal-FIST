# Thermal-FIST WASM Docker Deployment

This Docker setup builds QtThermalFIST for WebAssembly and serves it via Caddy.

## Quick Start

```bash
# Build the image (takes 15-30 minutes on first build)
docker compose build

# Run the container
docker compose up -d

# Access at http://localhost:8080
```

## Manual Build

```bash
# Build
docker build -t thermal-fist-wasm .

# Run
docker run -d -p 8080:80 --name thermal-fist-wasm thermal-fist-wasm
```

## Using with External Caddy

If you already have Caddy running on the server and want to proxy to this container:

1. Build the image:
   ```bash
   docker compose build
   ```

2. Run the container on an internal port:
   ```bash
   docker run -d -p 127.0.0.1:8081:80 --name thermal-fist-wasm thermal-fist-wasm
   ```

3. Add to your server's Caddyfile:
   ```caddyfile
   thermal-fist.example.com {
       reverse_proxy localhost:8081

       # Required for multi-threaded WASM
       header {
           Cross-Origin-Opener-Policy "same-origin"
           Cross-Origin-Embedder-Policy "require-corp"
       }
   }
   ```

## Build Arguments

The Dockerfile uses:
- **Emscripten**: 3.1.37 (required by Qt 6.6)
- **Qt**: 6.6.0 (wasm_multithread)
- **Thermal-FIST**: feature/wasm-support branch

To use a different branch, modify the Dockerfile:
```dockerfile
RUN git clone --depth 1 -b master https://github.com/vlvovch/Thermal-FIST.git
```

## Single-threaded Build

The default multi-threaded build requires COOP/COEP headers (already configured in Caddyfile).
For simpler deployment without special headers, use single-threaded:

1. Change the Qt installation in Dockerfile:
   ```dockerfile
   RUN aqt install-qt linux desktop 6.6.0 wasm_singlethread
   ```

2. Update the build command:
   ```dockerfile
   /opt/6.6.0/wasm_singlethread/bin/qt-cmake ..
   ```

## Troubleshooting

### Build fails with "Qt not found"
The aqtinstall tool may fail if Qt servers are slow. Retry the build:
```bash
docker compose build --no-cache
```

### Large image size
The builder stage is ~10GB due to Qt and Emscripten. The final image is much smaller (~50MB) since only the built WASM files are copied.

### Browser shows blank page
Check browser console for errors. Ensure your browser supports WebAssembly.
