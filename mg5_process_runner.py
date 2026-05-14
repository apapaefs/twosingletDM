from generate_mg5_trsm_xsecs import ProcLocation, get_mg5_xsec


def selected_mg5_processes(processes_to_run, proc_location=None):
    """Return the configured MG5 processes that are available to this scan."""
    if proc_location is None:
        proc_location = ProcLocation

    selected = []
    for process in processes_to_run:
        if process not in proc_location:
            available = ", ".join(sorted(proc_location.keys()))
            print(
                f"Warning: MG5 process '{process}' is not available; "
                f"available processes are: {available}. Skipping."
            )
            continue
        if process in selected:
            print(f"Warning: MG5 process '{process}' was requested more than once. Skipping duplicate.")
            continue
        selected.append(process)
    return selected


def run_mg5_processes(
    processes_to_run,
    run_tag,
    lambdas,
    k1,
    k2,
    k3,
    m2,
    w2,
    m3,
    w3,
    energy,
    *,
    get_xsec=get_mg5_xsec,
    proc_location=None,
):
    """Run all selected MG5 processes and return `{process: cross_section_pb}`."""
    if len(processes_to_run) == 0:
        print("Warning: RunMG5 is True but MG5ProcessesToRun is empty; no MG5 processes will be run.")
        return {}

    mg5xsecs = {}
    for process in selected_mg5_processes(processes_to_run, proc_location=proc_location):
        print("Running MG5 process", process)
        mg5xsecs[process] = get_xsec(
            process,
            run_tag,
            lambdas,
            k1,
            k2,
            k3,
            m2,
            w2,
            m3,
            w3,
            ecm=energy,
        )
        print("MG5", process, "xsec [pb] =", mg5xsecs[process])
    return mg5xsecs
