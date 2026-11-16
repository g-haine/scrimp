# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Opt-in local wall-clock timings with collective MPI summaries.

No MPI collective occurs inside a measured region or a PETSc callback.
Only `report()` communicates, and callers MUST invoke it on all ranks.
Timings are inclusive; PETSc TS includes time spent in its Python
callbacks and the callbacks include GetFEM/bridge work. Do not sum these
nested measurements as if they were disjoint contributions.
"""
from collections import defaultdict
from contextlib import contextmanager, nullcontext
import os
import time


class TimingProfiler:
    """Per-DPHS, per-solve profiler; enabled with SCRIMP_PROFILE=1.

    Timers use monotonic `perf_counter`; they introduce no barriers.
    Disabled by default, and report() returns None without collectives
    when disabled. A report on an MPI communicator is collective.
    """

    def __init__(self, enabled=None):
        self.enabled = (
            os.environ.get("SCRIMP_PROFILE", "0") == "1"
            if enabled is None else bool(enabled)
        )
        self.reset()

    def reset(self):
        self._seconds = defaultdict(float)
        self._calls = defaultdict(int)

    def measure(self, label):
        """Return a context manager for a named, rank-local phase."""
        if not self.enabled:
            return nullcontext()
        return self._measure_enabled(label)

    @contextmanager
    def _measure_enabled(self, label):
        start = time.perf_counter()
        try:
            yield
        finally:
            self._seconds[label] += time.perf_counter() - start
            self._calls[label] += 1

    def report(self, mpi_comm):
        """Collect timing min/mean/max and call counts over MPI ranks.

        All ranks must enter. Values are *per rank*, not a sum of
        parallel wall durations. No extra MPI synchronization occurs
        in the assembly or solver callbacks.
        """
        if not self.enabled:
            return None

        local = {
            name: {"seconds": self._seconds[name], "calls": self._calls[name]}
            for name in self._seconds
        }
        # PETSc's TS.solve includes its Python callbacks. The remainder
        # contains PETSc internal solver, MUMPS factor/solve, waiting,
        # monitors and other native work -- NOT pure MUMPS time.
        for phase in ("initial", "main"):
            key = f"{phase}.petsc.ts_solve"
            if key in local:
                callbacks = sum(
                    self._seconds.get(f"{phase}.petsc.callback.{cb}", 0.)
                    for cb in ("ifunction", "ijacobian")
                )
                local[f"{phase}.petsc.non_callback_remainder"] = {
                    "seconds": max(0., self._seconds[key] - callbacks),
                    "calls": self._calls[key],
                }

        all_ranks = mpi_comm.allgather(local)
        events = {}
        for name in sorted(set().union(*(d.keys() for d in all_ranks))):
            values = [d.get(name, {"seconds": 0., "calls": 0})
                      for d in all_ranks]
            durations = [x["seconds"] for x in values]
            counts = [x["calls"] for x in values]
            events[name] = {
                "seconds_rank_min": min(durations),
                "seconds_rank_mean": sum(durations) / len(durations),
                "seconds_rank_max": max(durations),
                "calls_rank_min": min(counts),
                "calls_rank_max": max(counts),
            }
        return {
            "enabled": True,
            "mpi_size": len(all_ranks),
            "events": events,
            "notes": [
                "Timers are inclusive; nested phases must not be summed.",
                "MPI aggregation is per-event rank min/mean/max; no "
                "collectives were added to solver callbacks.",
                "petsc.non_callback_remainder is TS.solve minus Python "
                "callback time; it is NOT an isolated MUMPS measurement.",
                "For MUMPS factorization and PETSc event detail, use "
                "PETSc native logging (-log_view) in a separate run.",
            ],
        }
