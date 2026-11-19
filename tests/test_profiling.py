# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Test opt-in profiling without importing GetFEM, PETSc or SCRIMP."""
import importlib.util
from pathlib import Path
import time
import unittest


FILE = Path(__file__).resolve().parents[1] / "scrimp/utils/profiling.py"
spec = importlib.util.spec_from_file_location("_scrimp_profiling", FILE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class Communicator:
    """Deterministic two-rank stand-in; reject unforeseen collectives."""
    def __init__(self):
        self.calls = 0

    def allgather(self, local):
        self.calls += 1
        return [local, local]


class TestTimingProfiler(unittest.TestCase):
    def test_off_is_collective_free(self):
        profiler = mod.TimingProfiler(enabled=False)
        comm = Communicator()
        with profiler.measure("main.petsc.ts_solve"):
            pass
        self.assertIsNone(profiler.report(comm))
        self.assertEqual(comm.calls, 0)

    def test_nested_timing_is_explicitly_inclusive(self):
        profiler = mod.TimingProfiler(enabled=True)
        comm = Communicator()
        with profiler.measure("main.petsc.ts_solve"):
            with profiler.measure("main.petsc.callback.ifunction"):
                time.sleep(0.001)
            with profiler.measure("main.getfem.rhs"):
                time.sleep(0.001)
        self.assertEqual(comm.calls, 0)  # never in callback/assembly
        result = profiler.report(comm)
        self.assertEqual(comm.calls, 1)
        self.assertEqual(result["mpi_size"], 2)
        events = result["events"]
        self.assertEqual(
            events["main.petsc.callback.ifunction"]["calls_rank_min"], 1
        )
        self.assertEqual(
            events["main.petsc.ts_solve"]["calls_rank_max"], 1
        )
        self.assertGreaterEqual(
            events["main.petsc.ts_solve"]["seconds_rank_min"],
            events["main.petsc.callback.ifunction"]["seconds_rank_max"]
        )
        self.assertGreaterEqual(
            events["main.petsc.non_callback_remainder"]["seconds_rank_min"],
            0.
        )
        profiler.reset()
        self.assertFalse(profiler.report(comm)["events"])

    def test_context_exception_still_records(self):
        profiler = mod.TimingProfiler(enabled=True)
        with self.assertRaisesRegex(ValueError, "sample"):
            with profiler.measure("error"):
                raise ValueError("sample")
        self.assertEqual(profiler._calls["error"], 1)


if __name__ == "__main__":
    unittest.main()
