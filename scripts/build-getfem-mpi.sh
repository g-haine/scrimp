#!/usr/bin/env bash
# Build the official GetFEM release against the *active* conda/MPI stack.
# The source tarball avoids tracking upstream development HEAD.
set -euo pipefail

: "${CONDA_PREFIX:?Activate the target conda/micromamba environment first}"
GETFEM_VERSION="${GETFEM_VERSION:-5.5}"
GETFEM_SOURCE="${GETFEM_SOURCE:-release}"  # release or git
GETFEM_REF="${GETFEM_REF:-master}"
GETFEM_PREFIX="${GETFEM_PREFIX:-$CONDA_PREFIX}"
MAKE_JOBS="${MAKE_JOBS:-4}"
BUILD_DIR="${GETFEM_BUILD_DIR:-$PWD/.build/getfem}"
mkdir -p "$BUILD_DIR"

export PATH="$CONDA_PREFIX/bin:$PATH"
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export CPPFLAGS="-I$CONDA_PREFIX/include ${CPPFLAGS:-}"
export LDFLAGS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib ${LDFLAGS:-}"
export CC="${CC:-mpicc}"
export CXX="${CXX:-mpicxx}"
export FC="${FC:-mpifort}"
PYTHON_BIN="$(command -v python)"

if [[ "$GETFEM_SOURCE" == release ]]; then
  file="$BUILD_DIR/getfem-$GETFEM_VERSION.tar.gz"
  url="https://download-mirror.savannah.gnu.org/releases/getfem/stable/getfem-$GETFEM_VERSION.tar.gz"
  if [[ ! -f "$file" ]]; then
    curl --fail --location --retry 3 --output "$file" "$url"
  fi
  if [[ -n "${GETFEM_SHA256:-}" ]]; then
    echo "$GETFEM_SHA256  $file" | sha256sum -c -
  else
    echo "WARNING: no GETFEM_SHA256 configured (experimental build)." >&2
  fi
  tar -xzf "$file" -C "$BUILD_DIR"
  src="$BUILD_DIR/getfem-$GETFEM_VERSION"
elif [[ "$GETFEM_SOURCE" == git ]]; then
  command -v git >/dev/null
  src="$BUILD_DIR/getfem-git"
  if [[ ! -d "$src/.git" ]]; then
    git clone https://git.savannah.nongnu.org/git/getfem.git "$src"
  fi
  git -C "$src" fetch --tags origin
  git -C "$src" checkout --detach "$GETFEM_REF"
  ( cd "$src" && bash autogen.sh )
else
  echo "Unsupported GETFEM_SOURCE=$GETFEM_SOURCE (release or git)" >&2
  exit 2
fi

echo "GetFEM=$GETFEM_VERSION source=$GETFEM_SOURCE"
echo "MPI=$("$CC" -show 2>/dev/null || command -v "$CC")"
echo "Python=$PYTHON_BIN NumPy=$(python -c 'import numpy; print(numpy.__version__)')"
cd "$src"
./configure \
  --prefix="$GETFEM_PREFIX" \
  --enable-python \
  --enable-paralevel=2 \
  --with-pic \
  PYTHON="$PYTHON_BIN" CC="$CC" CXX="$CXX" FC="$FC"

# Fail here rather than shipping an accidentally *serial* GetFEM.
if ! grep -E 'GETFEM_PARA_LEVEL[^0-9]*2|enable.paralevel[^0-9]*2' config.log >/dev/null; then
  echo "NOTE: verify configure's MPI level 2 from the runtime smoke." >&2
fi
make -j "$MAKE_JOBS"
make install
python - <<'PY'
import getfem
level = int(getfem.util_mpi_parallelism_level())
assert level == 2, f"GetFEM compiled without MPI level 2: {level}"
print("GetFEM MPI level:", level)
PY
