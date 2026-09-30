"""
Number of threads used by pyscal3's C++ routines.

By default pyscal3 uses the CPUs this process may run on. The environment
variable ``PYSCAL_NUM_THREADS``, or else ``OMP_NUM_THREADS``, sets a different
default when pyscal3 is imported; :func:`set_num_threads` changes it at run
time. Results do not depend on the number of threads.
"""

import os

import pyscal3.csystem as pc


def _available_cpus():
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def _default_threads():
    for name in ("PYSCAL_NUM_THREADS", "OMP_NUM_THREADS"):
        value = os.environ.get(name, "").strip()
        if value:
            try:
                n = int(value.split(",")[0])
            except ValueError:
                continue
            if n > 0:
                return n
    return _available_cpus()


def set_num_threads(n=None):
    """Set the number of threads used by pyscal3.

    Parameters
    ----------
    n : int or None
        Number of threads. None (or a value below 1) restores the default:
        ``PYSCAL_NUM_THREADS``, else ``OMP_NUM_THREADS``, else the number of
        CPUs available to the process.
    """
    if n is None or int(n) < 1:
        n = _default_threads()
    pc.set_num_threads(int(n))


def get_num_threads():
    """Return the number of threads used by pyscal3."""
    return pc.get_num_threads()
