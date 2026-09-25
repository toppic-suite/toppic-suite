#!/usr/bin/env bash
# Install the Ubuntu packages needed to build the TopPIC Suite from source.
#
# The package list mirrors the one in docs/compile/ubuntu_build.md; keep the
# two in step. Other third-party code (htslib, ProteoWizard, ONNX Runtime) is
# vendored under ext/ and needs no package.
#
# Usage: scripts/install_ubuntu_packages.sh
# Run as a user with sudo rights (or as root).

set -euo pipefail

if ! command -v apt-get >/dev/null 2>&1; then
    echo "error: apt-get not found; this script is for Ubuntu/Debian systems" >&2
    exit 1
fi

if [ "$(id -u)" -eq 0 ]; then
    SUDO=""
elif command -v sudo >/dev/null 2>&1; then
    SUDO="sudo"
else
    echo "error: this script must run as root or with sudo available" >&2
    exit 1
fi

# Build tools and libraries (Git LFS is needed for the model files under res/
# and the ONNX Runtime library under ext/onnx/).
PACKAGES=(
    build-essential
    cmake
    clang
    git
    git-lfs
    zlib1g-dev
    libsqlite3-dev
    libpugixml-dev
    libboost-filesystem-dev
    libboost-iostreams-dev
    libboost-thread-dev
    libboost-chrono-dev
    libboost-serialization-dev
    libboost-program-options-dev
    qt6-base-dev    # GUI tools (topfd_gui, toppic_gui, ...)
)

$SUDO apt-get update
$SUDO apt-get install -y "${PACKAGES[@]}"

echo "All build dependencies are installed."
echo "Next: 'git lfs install', then configure and build with CMake"
echo "(see docs/compile/ubuntu_build.md)."
