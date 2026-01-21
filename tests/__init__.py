"""Test helpers for optional dependencies."""

import sys
import types

if "petsc4py" not in sys.modules:  # pragma: no cover - testing helper
    fake_petsc4py = types.ModuleType("petsc4py")

    def _init(*args, **kwargs):
        return None

    class _FakeComm:
        def getRank(self) -> int:
            return 0

    class _FakePETSc:
        COMM_WORLD = _FakeComm()

    fake_petsc4py.init = _init
    fake_petsc4py.PETSc = _FakePETSc()
    sys.modules["petsc4py"] = fake_petsc4py

if "getfem" not in sys.modules:  # pragma: no cover - testing helper
    fake_getfem = types.ModuleType("getfem")

    class _FakeModel:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    def _asm(*args, **kwargs):  # pragma: no cover - helper
        raise RuntimeError("getfem is not available in the test environment")

    fake_getfem.Model = _FakeModel
    fake_getfem.asm = _asm
    fake_getfem.util_trace_level = lambda level: None
    fake_getfem.util_warning_level = lambda level: None
    sys.modules["getfem"] = fake_getfem

if "numpy" not in sys.modules:  # pragma: no cover - testing helper
    fake_numpy = types.ModuleType("numpy")
    fake_numpy.array = lambda data, *args, **kwargs: data
    fake_numpy.zeros = lambda shape, *args, **kwargs: [0] * (shape if isinstance(shape, int) else 1)
    fake_numpy.ones = lambda shape, *args, **kwargs: [1] * (shape if isinstance(shape, int) else 1)
    fake_numpy.pi = 3.141592653589793
    sys.modules["numpy"] = fake_numpy

if "matplotlib" not in sys.modules:  # pragma: no cover - testing helper
    fake_matplotlib = types.ModuleType("matplotlib")
    fake_pyplot = types.SimpleNamespace(
        plot=lambda *args, **kwargs: None,
        figure=lambda *args, **kwargs: None,
        show=lambda *args, **kwargs: None,
    )
    fake_matplotlib.pyplot = fake_pyplot
    sys.modules["matplotlib"] = fake_matplotlib
    sys.modules["matplotlib.pyplot"] = fake_pyplot
