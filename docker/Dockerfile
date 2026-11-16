# Modern, single-node OpenMPI/GetFEM 5.5/SCRIMP MPI build.
# EXPERIMENTAL: not yet native-build validated. Build + tests must pass
# before promoting or publishing a GHCR release tag.
FROM mambaorg/micromamba:2.10.0-ubuntu24.04

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG GETFEM_VERSION=5.5
ARG MAKE_JOBS=4

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates binutils file libgl1 libglu1-mesa libxrender1 libxext6 \
    libxcursor1 libxinerama1 libxi6 libxrandr2 libxfixes3 libsm6 libice6 \
    libxft2 libfontconfig1 libfreetype6 libxmu6 libxt6 libx11-6 libxcb1 \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/scrimp-environment.yml
RUN micromamba create -y -f /tmp/scrimp-environment.yml \
    && micromamba clean -a -y

# Fail early on missing transitive shared libraries of the pip Gmsh wheel.
# This is intentionally BEFORE the expensive PETSc bindings and GetFEM C++
# compilation; missing libXcursor/libXft previously surfaced minutes later.
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

WORKDIR /workspace/scrimp
COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace/scrimp

# Compile the Python wrappers *against* NumPy 2 and conda's
# complex PETSc/SLEPc libraries (conda petsc4py currently pins NumPy 1).
RUN bash scripts/install-petsc-python-bindings.sh

# The same GetFEM source-build script also works in local conda.
RUN GETFEM_VERSION="${GETFEM_VERSION}" \
    MAKE_JOBS="${MAKE_JOBS}" \
    bash scripts/build-getfem-mpi.sh

RUN python -m pip install --no-deps --no-build-isolation -e . \
    && python sandbox/smoke_modern_stack.py \
    && micromamba env export -n scrimp --explicit > /workspace/scrimp/docker-resolved-spec.txt

CMD ["/bin/bash"]
