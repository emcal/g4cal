#!/bin/bash
# Build calsim + calreco variants + adapter tests, then run the tests.
# Geant4 prefix defaults to /app/geant4; override with GEANT4_DIR.
set -e
cd "$(dirname "$0")/.."
: "${GEANT4_DIR:=/app/geant4}"
cmake -B build -DCMAKE_PREFIX_PATH="$GEANT4_DIR"
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
