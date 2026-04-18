#!/bin/bash
# 01_build_gromacs_plumed.sh  —  CRYPTAD
# Compile PLUMED 2.9.2 + GROMACS 2023.3 (MPI + CUDA) on ARF orfoz partition.
#
# PRE-REQUISITE (run on arf-ui1 before submitting):
#   mkdir -p $HOME/build_plumed_gmx && cd $HOME/build_plumed_gmx
#   wget https://github.com/plumed/plumed2/releases/download/v2.9.2/plumed-2.9.2.tgz
#   wget https://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz
#
#   Verify checksums before submitting (obtain from release pages):
#     PLUMED: https://github.com/plumed/plumed2/releases/tag/v2.9.2
#     GROMACS: https://manual.gromacs.org/documentation/2023.3/download.html
#   sha256sum plumed-2.9.2.tgz gromacs-2023.3.tar.gz
#
# Submit from CRYPTAD project root:
#   sbatch 09_scripts/02_md_setup/01_build_gromacs_plumed.sh
#
#SBATCH --job-name=build_gmx_plumed
#SBATCH --partition=orfoz
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=04:00:00
#SBATCH --output=build_gmx_plumed_%j.out
#SBATCH --error=build_gmx_plumed_%j.err

set -euo pipefail

# ─── Pinned versions ──────────────────────────────────────────────────────────
# Changing these requires updating both the tarball names below AND the
# checksums in VERIFY_SHA256 so that the integrity check catches the mismatch.
PLUMED_VERSION="2.9.2"
GMX_VERSION="2023.3"

# ─── Expected SHA-256 of source tarballs ─────────────────────────────────────
# Obtain from the official release pages listed in the PRE-REQUISITE header.
# Paste the full hex digest here; the build will abort if they don't match.
PLUMED_SHA256="PASTE_SHA256_HERE_from_github.com/plumed/plumed2/releases/tag/v2.9.2"
GMX_SHA256="PASTE_SHA256_HERE_from_manual.gromacs.org/documentation/2023.3/download.html"

# ─── Modules ─────────────────────────────────────────────────────────────────
module purge
module load comp/gcc/12.3.0
module load lib/openmpi/4.1.6
module load lib/cuda/12.4
module load comp/cmake/3.31.1

echo "--- Compiler / MPI / CUDA versions ---"
gcc      --version | head -1
mpicc    --version | head -1
nvcc     --version | head -1
cmake    --version | head -1
echo "--------------------------------------"

# ─── Guard: CUDA_HOME must be set by lib/cuda/12.4 module ────────────────────
: "${CUDA_HOME:?lib/cuda/12.4 module did not export CUDA_HOME.
  Run 'module show lib/cuda/12.4' to find the correct path,
  then set it explicitly: export CUDA_HOME=/path/to/cuda/12.4}"

# ─── Paths ───────────────────────────────────────────────────────────────────
BUILDDIR=$HOME/build_plumed_gmx
PLUMED_INSTALL=$HOME/apps/plumed-${PLUMED_VERSION}
GMX_INSTALL=$HOME/apps/gromacs-plumed-${GMX_VERSION}
NPROC=56

cd "$BUILDDIR"

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
        echo "       Retrieve the original tarball from the official release page." >&2
        exit 1
    fi
    echo "[OK]   SHA-256 verified: $file"
}

verify_sha256 "plumed-${PLUMED_VERSION}.tgz"       "$PLUMED_SHA256"
verify_sha256 "gromacs-${GMX_VERSION}.tar.gz"       "$GMX_SHA256"

# ─── 1. Build PLUMED ─────────────────────────────────────────────────────────
echo ""
echo "=== Building PLUMED ${PLUMED_VERSION} ==="
rm -rf "plumed-${PLUMED_VERSION}"           # ensure clean extraction
tar xzf "plumed-${PLUMED_VERSION}.tgz"
cd "plumed-${PLUMED_VERSION}"

# Use -mavx2 (explicit, portable across orfoz nodes) instead of -march=native
# (which embeds the build node's microarchitecture and is not reproducible).
# AVX2 is consistent with -DGMX_SIMD=AVX2_256 used below for GROMACS.
./configure \
  --prefix="$PLUMED_INSTALL" \
  --enable-mpi \
  CC=mpicc CXX=mpicxx \
  CXXFLAGS="-O3 -mavx2"

make -j${NPROC}
make install

export PATH=$PLUMED_INSTALL/bin:$PATH
export PLUMED_KERNEL=$PLUMED_INSTALL/lib/libplumedKernel.so
export LD_LIBRARY_PATH=$PLUMED_INSTALL/lib:${LD_LIBRARY_PATH:-}
echo "PLUMED installed: $(plumed info 2>&1 | head -1)"

# ─── 2. Patch + build GROMACS ────────────────────────────────────────────────
echo ""
echo "=== Patching GROMACS ${GMX_VERSION} with PLUMED ==="
cd "$BUILDDIR"
rm -rf "gromacs-${GMX_VERSION}"             # ensure clean extraction
tar xzf "gromacs-${GMX_VERSION}.tar.gz"
cd "gromacs-${GMX_VERSION}"

# Apply PLUMED patch (--shared = runtime-loadable kernel)
plumed patch -p --shared -e "gromacs-${GMX_VERSION}"

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

make -j${NPROC}
make install

# ─── 3. Build verification ────────────────────────────────────────────────────
echo ""
echo "=== Verifying build ==="
if [[ ! -f "$GMX_INSTALL/bin/gmx_mpi" ]]; then
    echo "[FAIL] gmx_mpi not found at $GMX_INSTALL/bin/gmx_mpi" >&2
    exit 1
fi
if [[ ! -f "$PLUMED_INSTALL/lib/libplumedKernel.so" ]]; then
    echo "[FAIL] libplumedKernel.so not found" >&2
    exit 1
fi
# Quick functional check: gmx_mpi must print version without error
"$GMX_INSTALL/bin/gmx_mpi" --version > /dev/null 2>&1 \
    || { echo "[FAIL] gmx_mpi --version exited non-zero" >&2; exit 1; }
echo "[OK]   gmx_mpi functional"

# ─── 4. Write build manifest ─────────────────────────────────────────────────
MANIFEST="$GMX_INSTALL/build_manifest.txt"
{
    echo "built_at:        $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "plumed_version:  $(plumed info 2>&1 | head -1)"
    echo "gromacs_version: $("$GMX_INSTALL/bin/gmx_mpi" --version 2>&1 | grep 'GROMACS version')"
    echo "gcc:             $(gcc --version | head -1)"
    echo "mpi:             $(mpicc --version | head -1)"
    echo "cuda:            $(nvcc --version | grep 'release')"
    echo "build_host:      $(hostname)"
    echo "plumed_sha256:   $PLUMED_SHA256"
    echo "gromacs_sha256:  $GMX_SHA256"
} > "$MANIFEST"
echo "Build manifest: $MANIFEST"

# ─── 5. Write module-load helper ─────────────────────────────────────────────
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
echo "Binaries : $GMX_INSTALL/bin/gmx_mpi"
echo "Manifest : $MANIFEST"
echo "Helper   : source \$HOME/load_gmx_plumed.sh"
"$GMX_INSTALL/bin/gmx_mpi" --version | head -10
