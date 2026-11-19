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
export CC="${GETFEM_CC:-mpicc}"
export CXX="${GETFEM_CXX:-mpicxx}"
export FC="${GETFEM_FC:-mpifort}"

# Autoconf/libtool can silently substitute "false" for an absent
# archiver and only fail after compiling all GetFEM C++ objects:
#   libtool: link: false cr .libs/libgetfem.a ...
# Check the complete native toolchain *before* downloading/compiling.
for tool in "$CC" "$CXX" "$FC" ar ranlib objdump nm file make; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: missing required GetFEM build tool '$tool'." >&2
    echo "Ubuntu: sudo apt-get install binutils file" >&2
    exit 2
  fi
done
export AR="$(command -v ar)"
export RANLIB="$(command -v ranlib)"
PYTHON_BIN="$(command -v python)"

if [[ "$GETFEM_SOURCE" == release ]]; then
  # Optional separate DOWNLOAD cache: Docker mounts this directory as a
  # BuildKit cache, while its actual compiled sources remain in BUILD_DIR.
  # The default preserves the original host-conda behavior.
  archive_dir="${GETFEM_ARCHIVE_DIR:-$BUILD_DIR}"
  mkdir -p "$archive_dir"
  file="$archive_dir/getfem-$GETFEM_VERSION.tar.gz"
  url="https://download-mirror.savannah.gnu.org/releases/getfem/stable/getfem-$GETFEM_VERSION.tar.gz"
  if [[ ! -f "$file" ]]; then
    # Do not leave a corrupted partial archive in the persistent cache
    # when the download is interrupted or a mirror fails.
    incomplete="${file}.part.${BASHPID}"
    trap 'rm -f "$incomplete"' EXIT
    curl --fail --location --retry 3 --output "$incomplete" "$url"
    mv "$incomplete" "$file"
    trap - EXIT
  fi
  if [[ -n "${GETFEM_SHA256:-}" ]]; then
    if ! echo "$GETFEM_SHA256  $file" | sha256sum -c -; then
      rm -f "$file"
      echo "ERROR: GetFEM archive SHA256 mismatch; corrupt cached archive removed." >&2
      exit 2
    fi
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

# Apply every approved GetFEM adjustment (including MPI timing verbosity)
# as a version-pinned standard unified diff, in numeric filename order.
# No special-purpose source-rewriting Python script is needed.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Apply small SCRIMP-maintained GetFEM 5.5 source patches before compiling.
# Keep patches separate from SCRIMP's solver, and fail rather than silently
# building an unsupported source. The original release archive is untouched.
if [[ "$GETFEM_VERSION" == "5.5" ]]; then
  if [[ "$GETFEM_SOURCE" == git && "$GETFEM_REF" != "v5.5" && "$GETFEM_REF" != "c242f7fe0989e9a90846c63eaccc6b495d774724" ]]; then
    echo "ERROR: SCRIMP GetFEM 5.5 patches require the v5.5 tag/commit, not GETFEM_REF=$GETFEM_REF" >&2
    exit 2
  fi
  command -v patch >/dev/null || {
    echo "ERROR: 'patch' is required to apply SCRIMP GetFEM source fixes." >&2
    exit 2
  }
  for source_patch in "$SCRIPT_DIR"/patches/getfem-5.5/*.patch; do
    [[ -f "$source_patch" ]] || continue
    echo "GetFEM source patch: $(basename "$source_patch")"
    if patch --batch --forward --silent --dry-run -p1 -d "$src" < "$source_patch"; then
      patch --batch --forward -p1 -d "$src" < "$source_patch"
    elif patch --batch --silent --dry-run -R -p1 -d "$src" < "$source_patch"; then
      echo "Already applied: $(basename "$source_patch")"
    else
      echo "ERROR: GetFEM patch cannot be applied cleanly: $source_patch" >&2
      exit 2
    fi
  done
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

# GetFEM 5.5 ships GNU Automake 1.16's obsolete "py-compile" helper.
# That helper still uses "import imp", removed in Python 3.12; the
# shared library and _getfem extension are already compiled by now.
# Override only Automake's optional byte-compilation stage while
# keeping all library, header and extension installation active.
# The command-line variable propagates to recursive make invocations.
# Compile real .py modules afterwards using the supported stdlib API.
echo "Installing GetFEM without legacy Automake py-compile (Python >=3.12)"
make install 'am__py_compile=:'
python - <<'PYCOMPILE'
import compileall
import pathlib
import sysconfig
package = pathlib.Path(sysconfig.get_path("platlib")) / "getfem"
if not package.is_dir():
    raise RuntimeError(f"GetFEM Python package missing after install: {package}")
if not compileall.compile_dir(str(package), quiet=1, force=True):
    raise RuntimeError(f"Python byte-compilation failed for {package}")
print("GetFEM Python byte-compiled:", package)
PYCOMPILE

python - <<'PY'
import getfem
level = int(getfem.util_mpi_parallelism_level())
assert level == 2, f"GetFEM compiled without MPI level 2: {level}"
print("GetFEM MPI level:", level)
PY
