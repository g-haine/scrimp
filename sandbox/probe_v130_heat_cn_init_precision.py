"""Experimental tighter DAE initializer profile for the canonical CN heat pilot.

No edits to GetFEM operators, PDE, algebraic constraints, time scheme, or
accepted FE physics. The standalone canonical user example is invoked
with exactly the same CLI parameters as the baseline pilot. Only
documented SCRIMP/PETSc initializer options differ. This is NOT the
qualified default unless native results prove the original FE parity.
"""
from petsc4py import PETSc

OPTIONS = {
    "init_algebraic_max_passes": "8",
    "init_algebraic_accuracy": "1e-13",
    "init_algebraic_ksp_rtol": "1e-14",
    "init_algebraic_ksp_atol": "1e-14",
}
options = PETSc.Options()
for name, value in OPTIONS.items():
    options[name] = value

if __name__ == "__main__":
    import json
    print("EXPERIMENTAL HEAT INITIALIZATION PRECISION PROFILE "
          + json.dumps(OPTIONS), flush=True)
    from examples.heat import main
    main()
