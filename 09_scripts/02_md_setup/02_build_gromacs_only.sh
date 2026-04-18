#!/bin/bash
# 02_build_gromacs_only.sh  —  CRYPTAD
# Resume script: PLUMED 2.9.2 already installed at $HOME/apps/plumed-2.9.2.
# This script only patches + builds GROMACS 2023.3.
#
# PRE-REQUISITE: run 01_build_gromacs_plumed.sh first, OR manually install
#   PLUMED 2.9.2 to $HOME/apps/plumed-2.9.2 and place the GROMACS tarball at:
#   $HOME/build_plumed_gmx/gromacs-2023.3.tar.gz
#
#   Verify checksum before submitting:
#     sha256sum gromacs-2023.3.tar.gz
#     (expected SHA-256 from https://manual.gromacs.org/documentation/2023.3/download.html)
#
# Submit from CRYPTAD project root:
#   sbatch 09_scripts/02_md_setup/02_build_gromacs_only.sh
#
#SBATCH --job-name=build_gmx_only
#SBATCH --partition=orfoz
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=03:00:00
#SBATCH --output=build_gmx_only_%j.out
#SBATCH --error=build_gmx_only_%j.err

set -euo pipefail

# ─── Pinned versions ──────────────────────────────────────────────────────────
# Both version strings are used to derive tarball names, install paths, and
# the plumed patch target.  Change here only; nowhere else in this script.
PLUMED_VERSION="2.9.2"
GMX_VERSION="2023.3"

# ─── Expected SHA-256 of GROMACS tarball ─────────────────────────────────────
# Obtain from https://manual.gromacs.org/documentation/2023.3/download.html
# and paste the full hex digest here.
GMX_SHA256="PASTE_SHA256_HERE_from_manual.gromacs.org/documentation/2023.3/download.html"

# ─── Modules ─────────────────────────────────────────────────────────────────
module purge
module load comp/gcc/12.3.0
module load lib/openmpi/4.1.6
module load lib/cuda/12.4
module load comp/cmake/3.31.1

# ─── Guard: CUDA_HOME must be exported by lib/cuda/12.4 ──────────────────────
: "${CUDA_HOME:?lib/cuda/12.4 module did not export CUDA_HOME.
  Run 'module show lib/cuda/12.4' to find the correct path,
  then set it explicitly: export CUDA_HOME=/path/to/cuda/12.4}"

# ─── Paths ───────────────────────────────────────────────────────────────────
BUILDDIR="$HOME/build_plumed_gmx"
PLUMED_INSTALL="$HOME/apps/plumed-${PLUMED_VERSION}"
GMX_INSTALL="$HOME/apps/gromacs-plumed-${GMX_VERSION}"
NPROC=56

# ─── Guard: PLUMED must already be installed ─────────────────────────────────
if [[ ! -f "$PLUMED_INSTALL/bin/plumed" ]]; then
    echo "[FAIL] PLUMED ${PLUMED_VERSION} not found at $PLUMED_INSTALL/bin/plumed" >&2
    echo "       Run 01_build_gromacs_plumed.sh first to build PLUMED." >&2
    exit 1
fi
if [[ ! -f "$PLUMED_INSTALL/lib/libplumedKernel.so" ]]; then
    echo "[FAIL] libplumedKernel.so not found — PLUMED install appears incomplete." >&2
    exit 1
fi

# ─── Activate PLUMED environment ─────────────────────────────────────────────
export PATH="$PLUMED_INSTALL/bin:$PATH"
export PLUMED_KERNEL="$PLUMED_INSTALL/lib/libplumedKernel.so"
export LD_LIBRARY_PATH="$PLUMED_INSTALL/lib:${LD_LIBRARY_PATH:-}"

echo "--- Versions ---"
gcc      --version | head -1
mpicc    --version | head -1
nvcc     --version | grep 'release'
cmake    --version | head -1
plumed   info 2>&1 | head -3
echo "----------------"

# ─── Tarball integrity check ─────────────────────────────────────────────────
verify_sha256() {
    local file="$1" expected="$2"
    if [[ "$expected" == PASTE_* ]]; then
        echo "[WARN] SHA-256 placeholder not filled for $file — skipping checksum." >&2
        return 0
    fi
    local actual
    actual=$(sha256sum "$file" | awk '{print $1}')
    if [[ "$actual" != "$expected" ]]; then
        echo "[FAIL] SHA-256 mismatch for $file" >&2
        echo "       expected: $expected" >&2
        echo "       actual:   $actual" >&2
        echo "       Retrieve the original from https://ftp.gromacs.org/gromacs/" >&2
        exit 1
    fi
    echo "[OK]   SHA-256 verified: $file"
}

cd "$BUILDDIR"
verify_sha256 "gromacs-${GMX_VERSION}.tar.gz" "$GMX_SHA256"

# ─── Patch + build GROMACS ───────────────────────────────────────────────────
rm -rf "gromacs-${GMX_VERSION}"             # ensure clean extraction

echo ""
echo "=== Extracting GROMACS ${GMX_VERSION} ==="
tar xzf "gromacs-${GMX_VERSION}.tar.gz"
cd "gromacs-${GMX_VERSION}"

echo ""
echo "=== Applying PLUMED patch ==="
plumed patch -p --shared -e "gromacs-${GMX_VERSION}"

echo ""
echo "=== Running CMake ==="
mkdir -p build && cd build

cmake .. \
  -DCMAKE_INSTALL_PREFIX="$GMX_INSTALL" \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DGMX_MPI=ON \
  -DGMX_THREAD_MPI=OFF \
  -DGMX_GPU=CUDA \
  -DCUDA_TOOLKIT_ROOT_DIR="$CUDA_HOME" \
  -DGMX_CUDA_TARGET_COMPUTE="80;86" \
  -DGMX_FFT_LIBRARY=fftw3 \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_SIMD=AVX2_256 \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_BUILD_TYPE=Release

echo ""
echo "=== Compiling (${NPROC} cores) ==="
make -j${NPROC}
make install

# ─── Build verification ───────────────────────────────────────────────────────
echo ""
echo "=== Verifying build ==="
if [[ ! -f "$GMX_INSTALL/bin/gmx_mpi" ]]; then
    echo "[FAIL] gmx_mpi not found at $GMX_INSTALL/bin/gmx_mpi" >&2
    exit 1
fi
"$GMX_INSTALL/bin/gmx_mpi" --version > /dev/null 2>&1 \
    || { echo "[FAIL] gmx_mpi --version exited non-zero" >&2; exit 1; }
echo "[OK]   gmx_mpi functional"

# ─── Build manifest ───────────────────────────────────────────────────────────
MANIFEST="$GMX_INSTALL/build_manifest.txt"
{
    echo "built_at:        $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "plumed_version:  $(plumed info 2>&1 | head -1)"
    echo "gromacs_version: $("$GMX_INSTALL/bin/gmx_mpi" --version 2>&1 | grep 'GROMACS version')"
    echo "gcc:             $(gcc --version | head -1)"
    echo "mpi:             $(mpicc --version | head -1)"
    echo "cuda:            $(nvcc --version | grep 'release')"
    echo "build_host:      $(hostname)"
    echo "gromacs_sha256:  $GMX_SHA256"
} > "$MANIFEST"
echo "Build manifest: $MANIFEST"

# ─── Write module-load helper ─────────────────────────────────────────────────
cat > "$HOME/load_gmx_plumed.sh" << 'EOF'
# Source this file to activate the PLUMED-patched GROMACS:
#   source $HOME/load_gmx_plumed.sh
module purge
module load comp/gcc/12.3.0
module load lib/openmpi/4.1.6
module load lib/cuda/12.4
export PATH=$HOME/apps/plumed-2.9.2/bin:$HOME/apps/gromacs-plumed-2023.3/bin:$PATH
export PLUMED_KERNEL=$HOME/apps/plumed-2.9.2/lib/libplumedKernel.so
export LD_LIBRARY_PATH=$HOME/apps/plumed-2.9.2/lib:$HOME/apps/gromacs-plumed-2023.3/lib:${LD_LIBRARY_PATH:-}
echo "Loaded: $(gmx_mpi --version 2>&1 | grep 'GROMACS version')"
echo "PLUMED: $(plumed info 2>&1 | head -3)"
EOF
chmod +x "$HOME/load_gmx_plumed.sh"

echo ""
echo "=== Build complete ==="
echo "Binary   : $GMX_INSTALL/bin/gmx_mpi"
echo "Manifest : $MANIFEST"
echo "Helper   : source \$HOME/load_gmx_plumed.sh"
"$GMX_INSTALL/bin/gmx_mpi" --version | head -10
