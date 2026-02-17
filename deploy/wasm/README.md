# Thermal-FIST WASM Docker Deployment

This Docker setup builds QtThermalFIST for WebAssembly and serves it via Caddy.

## Quick Start

```bash
# First time: build base image (Emscripten + Qt, ~20-30 min, done once)
./update_docker.sh --init

# Access at http://localhost:8081
```

## Updating

When Thermal-FIST code changes, update the deployment without rebuilding the base image:

```bash
./update_docker.sh
```

This only rebuilds Thermal-FIST itself (~5 min), not the Emscripten/Qt toolchain.

## Manual Build

```bash
# Build base image (only needed once, or when changing Emscripten/Qt versions)
docker build -f Dockerfile.base -t thermal-fist-wasm-base .

# Build Thermal-FIST and deploy
docker compose build --build-arg CACHEBUST=$(date +%s)
docker compose up -d
```

## Using with External Caddy

If you already have Caddy running on the server and want to proxy to this container:

1. Build and run:
   ```bash
   ./update_docker.sh --init
   ```

2. Add to your server's Caddyfile:
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

## Build Details

- **Base image** (`Dockerfile.base`): Emscripten 4.0.7 + Qt 6.10.1 (gcc_64 host + wasm_multithread)
- **App image** (`Dockerfile`): clones and builds Thermal-FIST from the `devel` branch
- **Serving**: Caddy with COOP/COEP headers for multi-threaded WASM

To use a different branch, modify the `git clone` line in `Dockerfile`:
```dockerfile
RUN git clone --depth 1 -b master https://github.com/vlvovch/Thermal-FIST.git
```

## Troubleshooting

### Build fails with "Base image not found"
Run `./update_docker.sh --init` to build the base image first.

### Build fails with "Qt not found"
The aqtinstall tool may fail if Qt servers are slow. Rebuild the base image:
```bash
docker build -f Dockerfile.base --no-cache -t thermal-fist-wasm-base .
```

### Large image size
The base image is ~10GB due to Qt and Emscripten. The final serving image is much smaller (~50MB) since only the built WASM files are copied.

### Browser shows blank page
Check browser console for errors. Ensure your browser supports WebAssembly.
