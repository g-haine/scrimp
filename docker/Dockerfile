# syntax=docker/dockerfile:1.7
# SCRIMP v1.3.0 CANDIDATE. Do not overwrite scrimp-mpi:local until qualified.
# Native dependencies are intentionally built BEFORE copying SCRIMP sources:
# edits to Python, examples or docs only rebuild the final inexpensive layers.
FROM mambaorg/micromamba:2.10.0-ubuntu24.04

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Keep APT index and .deb downloads in BuildKit caches instead of image layers.
USER root
RUN --mount=type=cache,id=scrimp-v130-apt-archives,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,id=scrimp-v130-apt-lists,target=/var/lib/apt/lists,sharing=locked \
    rm -f /etc/apt/apt.conf.d/docker-clean \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates binutils file libgl1 libglu1-mesa libxrender1 libxext6 \
        libxcursor1 libxinerama1 libxi6 libxrandr2 libxfixes3 libsm6 libice6 \
        libxft2 libfontconfig1 libfreetype6 libxmu6 libxt6 libx11-6 libxcb1

# BuildKit cache mounts start as root-owned (even with mode=0777 on
# an existing mount). Three previous user steps hit permission failures:
# micromamba pkgs.lock, pip caching and the GetFEM tarball download.
# Fix BOTH newly created and existing named caches from a root step;
# avoid any cache prune and preserve already-downloaded dependencies.
RUN --mount=type=cache,id=scrimp-v130-mamba-pkgs,target=/opt/conda/pkgs,sharing=locked,mode=0775 \
    --mount=type=cache,id=scrimp-v130-pip-cache,target=/tmp/scrimp-pip-cache,sharing=locked,mode=0775 \
    --mount=type=cache,id=scrimp-v130-getfem-tarballs,target=/tmp/scrimp-getfem-tarballs,sharing=locked,mode=0775 \
    chown -R "$MAMBA_USER:$MAMBA_USER" \
        /opt/conda/pkgs /tmp/scrimp-pip-cache /tmp/scrimp-getfem-tarballs \
    && test -w /opt/conda/pkgs \
    && test -w /tmp/scrimp-pip-cache \
    && test -w /tmp/scrimp-getfem-tarballs

USER $MAMBA_USER

# Verify actual non-root write access before invoking micromamba, pip or
# curl. Root's test -w alone is not proof that mambauser can write.
RUN --mount=type=cache,id=scrimp-v130-mamba-pkgs,target=/opt/conda/pkgs,sharing=locked,mode=0775 \
    --mount=type=cache,id=scrimp-v130-pip-cache,target=/tmp/scrimp-pip-cache,sharing=locked,mode=0775 \
    --mount=type=cache,id=scrimp-v130-getfem-tarballs,target=/tmp/scrimp-getfem-tarballs,sharing=locked,mode=0775 \
    for directory in /opt/conda/pkgs /tmp/scrimp-pip-cache /tmp/scrimp-getfem-tarballs; do \
        test -w "$directory" || { echo "ERROR: non-writable BuildKit cache $directory" >&2; exit 2; }; \
    done

# Micromamba's default package cache is /opt/conda/pkgs. This cache
# persists across builds but is not copied into the resulting image.
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/scrimp-environment.yml
RUN --mount=type=cache,id=scrimp-v130-mamba-pkgs,target=/opt/conda/pkgs,sharing=locked,mode=0775 \
    micromamba create -y -f /tmp/scrimp-environment.yml

# Guard the actual resolved native builds: the previous wildcard '*complex*'
# silently pulled cuda13_complex plus >1 GiB of CUDA libraries on a CPU host.
# Keep complex PETSc/SLEPc and NumPy 2 without CUDA dependencies.
RUN micromamba run -n scrimp python -c "import json, pathlib; \
    base=pathlib.Path('/opt/conda/envs/scrimp/conda-meta'); \
    pkgs=[json.loads(p.read_text()) for p in base.glob('*.json')]; \
    builds={p['name']: p.get('build', '') for p in pkgs}; \
    assert all(builds.get(k, '').startswith('complex_') for k in ('petsc','slepc')), builds; \
    assert not ({'libcublas', 'libcusolver', 'libcusparse', 'libcufft'} & builds.keys()), 'CUDA runtime unexpectedly resolved'; \
    print('CPU complex PETSc/SLEPc:', builds['petsc'], builds['slepc'])"

# Catch missing native Gmsh wheel runtime dependencies before costly builds.
RUN micromamba run -n scrimp python -c "import gmsh; print('Gmsh native library loaded:', gmsh.__version__)"

ENV ENV_NAME=scrimp \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MPLBACKEND=Agg \
    PETSC_DIR=/opt/conda/envs/scrimp \
    SLEPC_DIR=/opt/conda/envs/scrimp \
    LD_LIBRARY_PATH=/opt/conda/envs/scrimp/lib \
    SCRIMP_GETFEM_MATRIX_MODE=distributed \
    SCRIMP_GETFEM_RHS_MODE=auto \
    SCRIMP_GETFEM_SOURCE_RHS_MODE=raw

# WORKDIR would otherwise be created as root on the Docker builder.
# Set ownership ONCE, before the expensive native build layers.
USER root
RUN install -d -m 0755 -o "$MAMBA_USER" -g "$MAMBA_USER" /workspace/scrimp
USER $MAMBA_USER
WORKDIR /workspace/scrimp

# PETSc/SLEPc Python wrappers are built against conda-forge's exact
# native COMPLEX libraries and NumPy 2. SCRIMP edits do not invalidate this.
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    scripts/install-petsc-python-bindings.sh scripts/install-petsc-python-bindings.sh
RUN --mount=type=cache,id=scrimp-v130-pip-cache,target=/tmp/scrimp-pip-cache,sharing=locked,mode=0775 \
    SCRIMP_PIP_USE_CACHE=1 PIP_CACHE_DIR=/tmp/scrimp-pip-cache \
    bash scripts/install-petsc-python-bindings.sh

# Declare GetFEM-only build args as late as possible: changing MAKE_JOBS
# or the pinned tarball must not invalidate the apt/Conda/PETSc layers.
ARG GETFEM_VERSION=5.5
ARG MAKE_JOBS=4
# SHA256 of official GetFEM 5.5 source tarball, checked by the build script.
# Reviewed 2026-11-19 Savannah release variant, compared to immutable v5.5
# c242f7fe... + all approved patches (see V130_DOCKER_BUILD.md).
ARG GETFEM_SHA256=0e042770b86152b09be4bfd0ee160b8d10d49c3f8a0b9b06d55f1ffaa1d30fc4

# Only GetFEM's build script and approved 0000..0004 patches enter the
# next heavy layer. A change to SCRIMP itself does not rebuild GetFEM.
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    scripts/build-getfem-mpi.sh scripts/build-getfem-mpi.sh
# Copy ONLY approved *.patch files here: README/experimental edits must
# not invalidate an expensive native GetFEM compilation layer.
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    scripts/patches/getfem-5.5/*.patch scripts/patches/getfem-5.5/
RUN --mount=type=cache,id=scrimp-v130-getfem-tarballs,target=/tmp/scrimp-getfem-tarballs,sharing=locked,mode=0775 \
    GETFEM_ARCHIVE_DIR=/tmp/scrimp-getfem-tarballs \
    GETFEM_VERSION="${GETFEM_VERSION}" GETFEM_SHA256="${GETFEM_SHA256}" \
    MAKE_JOBS="${MAKE_JOBS}" bash scripts/build-getfem-mpi.sh

# Explicit application COPY: no recursive repository copy.
# The strict .dockerignore allowlist also excludes ALL local simulation
# results/meshes/plots, even if users generate them under an unusual name.
# Docs and notebooks do not belong in the runtime candidate image.
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    setup.py pyproject.toml README.md LICENSE ./
COPY --chown=$MAMBA_USER:$MAMBA_USER scrimp/ scrimp/
COPY --chown=$MAMBA_USER:$MAMBA_USER examples/ examples/
COPY --chown=$MAMBA_USER:$MAMBA_USER sandbox/ sandbox/
COPY --chown=$MAMBA_USER:$MAMBA_USER tests/ tests/
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ scripts/

RUN python -m pip install --no-deps --no-build-isolation -e . \
    && python sandbox/smoke_modern_stack.py \
    && micromamba env export -n scrimp --explicit > /workspace/scrimp/docker-resolved-spec.txt

# Build metadata deliberately comes after GetFEM: specifying a new
# revision must not force an expensive native rebuild.
ARG SCRIMP_REVISION=unknown
LABEL org.opencontainers.image.source="https://github.com/g-haine/scrimp" \
      org.opencontainers.image.title="SCRIMP v1.3.0 candidate (unqualified)" \
      org.opencontainers.image.revision="${SCRIMP_REVISION}"

CMD ["/bin/bash"]
