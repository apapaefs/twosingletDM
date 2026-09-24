"""Process and resource policy shared by single-host scan campaigns."""

import contextlib
import os
import signal
import time


THREAD_ENVIRONMENT = {
    name: "1" for name in (
        "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS", "TRSM_HIGHS_THREADS",
    )
}


def worker_environment(cache_dir=None):
    env = dict(os.environ, **THREAD_ENVIRONMENT)
    env["MPLBACKEND"] = "Agg"
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        env["MPLCONFIGDIR"] = str(cache_dir)
    return env


def available_cpus():
    """Respect Linux affinity, with portable fallbacks for older Python."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        counter = getattr(os, "process_cpu_count", os.cpu_count)
        return max(1, counter() or 1)


def job_limit(value):
    if value == "auto":
        return value
    try:
        number = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError("jobs must be a positive integer or 'auto'") from error
    if number <= 0:
        raise ValueError("jobs must be positive")
    return number


def effective_jobs(requested, remaining):
    cpus = available_cpus()
    return min(cpus, cpus if requested == "auto" else requested, remaining)


def group_alive(pgid):
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # Never treat an inaccessible process group as a stale owner.


def signal_group(pgid, signum):
    try:
        os.killpg(pgid, signum)
    except ProcessLookupError:
        pass


def stop_processes(processes, grace):
    """Only signal groups created by this supervisor, including descendants."""
    for process in processes:
        signal_group(process.pid, signal.SIGTERM)
    deadline = time.monotonic() + grace
    while time.monotonic() < deadline:
        for process in processes:
            process.poll()  # reap leaders even while descendants are exiting
        if not any(group_alive(process.pid) for process in processes):
            break
        time.sleep(0.05)
    for process in processes:
        signal_group(process.pid, signal.SIGKILL)
        process.wait()


@contextlib.contextmanager
def cancellation_signals():
    state = {"signal": None}

    def request_stop(signum, frame):
        state["signal"] = signum

    previous = {sig: signal.signal(sig, request_stop)
                for sig in (signal.SIGINT, signal.SIGTERM)}
    try:
        yield state
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)
