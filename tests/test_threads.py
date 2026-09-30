"""The thread-count setting of pyscal3."""
import os
import subprocess
import sys

import pyscal3


def test_set_and_get_num_threads():
    before = pyscal3.get_num_threads()
    try:
        pyscal3.set_num_threads(3)
        assert pyscal3.get_num_threads() == 3
        pyscal3.set_num_threads(1)
        assert pyscal3.get_num_threads() == 1
        pyscal3.set_num_threads(None)
        assert pyscal3.get_num_threads() >= 1
    finally:
        pyscal3.set_num_threads(before)


def _threads_at_import(env):
    code = "import pyscal3; print(pyscal3.get_num_threads())"
    full = {k: v for k, v in os.environ.items() if k not in ("PYSCAL_NUM_THREADS", "OMP_NUM_THREADS")}
    full.update(env)
    out = subprocess.run([sys.executable, "-c", code], env=full, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return int(out.stdout.strip())


def test_default_from_environment():
    assert _threads_at_import({"PYSCAL_NUM_THREADS": "2"}) == 2
    assert _threads_at_import({"OMP_NUM_THREADS": "3"}) == 3
    assert _threads_at_import({"PYSCAL_NUM_THREADS": "2", "OMP_NUM_THREADS": "3"}) == 2
    assert _threads_at_import({}) >= 1
