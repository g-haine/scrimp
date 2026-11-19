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

USER $MAMBA_USER

# Micromamba's default package cache is /opt/conda/pkgs. This cache
# persists across builds but is not copied into the resulting image.
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/scrimp-environment.yml
RUN --mount=type=cache,id=scrimp-v130-mamba-pkgs,target=/opt/conda/pkgs,sharing=locked,mode=0777 \
    micromamba create -y -f /tmp/scrimp-environment.yml

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
RUN --mount=type=cache,id=scrimp-v130-pip-cache,target=/tmp/scrimp-pip-cache,sharing=locked,mode=0777 \
    SCRIMP_PIP_USE_CACHE=1 PIP_CACHE_DIR=/tmp/scrimp-pip-cache \
    bash scripts/install-petsc-python-bindings.sh

# Declare GetFEM-only build args as late as possible: changing MAKE_JOBS
# or the pinned tarball must not invalidate the apt/Conda/PETSc layers.
ARG GETFEM_VERSION=5.5
ARG MAKE_JOBS=4
# SHA256 of official GetFEM 5.5 source tarball, checked by the build script.
ARG GETFEM_SHA256=f6c09fb2edd112b6abe04e0485534f5e965d6f8f9ad853308b7120bfc6722ab5

# Only GetFEM's build script and approved 0000..0004 patches enter the
# next heavy layer. A change to SCRIMP itself does not rebuild GetFEM.
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    scripts/build-getfem-mpi.sh scripts/build-getfem-mpi.sh
# Copy ONLY approved *.patch files here: README/experimental edits must
# not invalidate an expensive native GetFEM compilation layer.
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    scripts/patches/getfem-5.5/*.patch scripts/patches/getfem-5.5/
RUN --mount=type=cache,id=scrimp-v130-getfem-tarballs,target=/tmp/scrimp-getfem-tarballs,sharing=locked,mode=0777 \
    GETFEM_ARCHIVE_DIR=/tmp/scrimp-getfem-tarballs \
    GETFEM_VERSION="${GETFEM_VERSION}" GETFEM_SHA256="${GETFEM_SHA256}" \
    MAKE_JOBS="${MAKE_JOBS}" bash scripts/build-getfem-mpi.sh

# Keep the complete SCRIMP source and its frequently edited examples
# OUTSIDE the native dependency layers, for a fast iterative rebuild.
COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace/scrimp/
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
