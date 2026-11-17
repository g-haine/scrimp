# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Fast subprocess tests for the independent rank-partitioned binary writer."""
import importlib.util
import json
from pathlib import Path
import queue
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


SOURCE = Path(__file__).resolve().parents[1] / "sandbox" / "async_binary_trajectory.py"
SPEC = importlib.util.spec_from_file_location("_scrimp_binary_writer", SOURCE)
IO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IO)


class AsyncBinaryTrajectoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def _write_pair(self, dtype):
        process, socket_path = IO.start_writer(
            self.directory, 2, timeout=10.,
        )
        self.addCleanup(lambda: process.poll() is None and process.kill())
        first = IO.RankSender(
            socket_path, rank=0, ranks=2, global_size=11, lo=0, hi=6,
            dtype=dtype, max_queued=16,
        )
        second = IO.RankSender(
            socket_path, rank=1, ranks=2, global_size=11, lo=6, hi=11,
            dtype=dtype, max_queued=16,
        )
        self.addCleanup(lambda: first.close() if not first._closed else None)
        self.addCleanup(lambda: second.close() if not second._closed else None)
        initial = np.arange(11, dtype=dtype)
        for step in (0, 1, 3, 7):
            vector = initial + step
            first.enqueue(step, step * 0.025, vector[:6])
            second.enqueue(step, step * 0.025, vector[6:])
        first.close()
        second.close()
        self.assertEqual(process.wait(timeout=10), 0)
        manifest = self.directory / "trajectory_manifest.json"
        self.assertTrue(manifest.is_file())
        return IO.LazyTrajectory(manifest), initial

    def test_binary_writer_round_trip_real(self):
        trajectory, initial = self._write_pair(np.float64)
        self.assertEqual(len(trajectory), 4)
        self.assertEqual(trajectory.global_size, 11)
        np.testing.assert_array_equal(trajectory.steps, [0, 1, 3, 7])
        np.testing.assert_allclose(
            trajectory.times, [0., 0.025, 0.075, 0.175]
        )
        for i, step in enumerate((0, 1, 3, 7)):
            np.testing.assert_array_equal(trajectory[i], initial + step)
        np.testing.assert_array_equal(trajectory[-1], initial + 7)
        np.testing.assert_array_equal(
            trajectory[1:3], np.stack((initial + 1, initial + 3))
        )

    def test_binary_writer_round_trip_complex(self):
        trajectory, initial = self._write_pair(np.complex128)
        for i, step in enumerate((0, 1, 3, 7)):
            np.testing.assert_array_equal(trajectory[i], initial + step)
        self.assertEqual(trajectory[-1].dtype, np.complex128)

    def test_enqueue_rejects_nonlocal_shape_and_nans(self):
        process, path = IO.start_writer(self.directory, 1, timeout=10.)
        self.addCleanup(lambda: process.poll() is None and process.kill())
        sender = IO.RankSender(
            path, rank=0, ranks=1, global_size=3, lo=0, hi=3,
            dtype=np.float64,
        )
        try:
            with self.assertRaisesRegex(ValueError, "Rank data"):
                sender.enqueue(0, 0., np.arange(4, dtype=np.float64))
            with self.assertRaises(FloatingPointError):
                sender.enqueue(0, 0., np.array([1., float("nan"), 3.]))
            sender.enqueue(0, 0., np.arange(3, dtype=np.float64))
        finally:
            sender.close()
        self.assertEqual(process.wait(timeout=10), 0)

    def test_queue_overflow_reports_failure_instead_of_silent_loss(self):
        # Force put_nowait to raise, making this test independent of
        # actual OS socket speed and disk write timing.
        q = queue.Queue(maxsize=1)
        with patch.object(q, "put_nowait", side_effect=queue.Full):
            class FakeSender:
                _closed = False
                _error = None
                local_size = 1
                dtype = np.dtype("<f8")
                _pending = q
            with self.assertRaises(IO.AsyncOutputBackpressure):
                IO.RankSender.enqueue(
                    FakeSender(), 0, 0., np.array([1.], dtype=np.float64)
                )

    def test_mismatched_rank_sample_times_never_publish_manifest(self):
        proc, path = IO.start_writer(self.directory, 2, timeout=10.)
        self.addCleanup(lambda: proc.poll() is None and proc.kill())
        a = IO.RankSender(
            path, rank=0, ranks=2, global_size=4, lo=0, hi=2,
        )
        b = IO.RankSender(
            path, rank=1, ranks=2, global_size=4, lo=2, hi=4,
        )
        a.enqueue(0, 0., np.ones(2))
        b.enqueue(0, 1., np.ones(2))
        a.close()
        b.close()
        self.assertNotEqual(proc.wait(timeout=10), 0)
        self.assertFalse(
            (self.directory / "trajectory_manifest.json").exists()
        )


if __name__ == "__main__":
    unittest.main()
