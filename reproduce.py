"""Reproduce the finite checks reported in Appendix C, in a fresh directory."""
from __future__ import annotations
import argparse, hashlib, json, shutil, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parent
CORE = [
    ("ordinary_pieri_audit.py", []),
    ("ordinary_pieri_structural_audit.py", []),
    ("ordinary_resultant_colon_audit.py", []),
    ("lagrangian_independent.py", ["--max-n", "10"]),
    ("lagrangian_parameter_checks.py", []),
    ("lagrangian_full_q_flag.py", ["--max-n", "10"]),
    ("counterexample_certificate.py", []),
    ("negative_controls.py", []),
]
EXTENSIONS = [
    ("ordinary_new_full_conductor.py", ["ordinary_full_conductor.json"]),
    ("ordinary_new_lg_conductor.py", ["lagrangian_full_conductor.json"]),
    ("lagrangian_new_round2_conductors.py", []),
    ("ordinary_new_lg_ungraded_invariant.py", ["intrinsic_pairings.json"]),
    ("lagrangian_new_ordinary_not_p2xp2.py", []),
    ("lagrangian_new_round4_infinite_family.py", []),
    ("ordinary_new_k2_second_family.py", ["30", "two_row_pieri.json"]),
    ("lagrangian_new_even_offset_moment.py", []),
]

REVIEWED = [
    ("ordinary_added_criterion_audit.py", []),
    ("lagrangian_transition_exact.py", []),
    ("intersection_ring_exact.py", []),
    ("coefficient_intrinsic_exact.py", []),
    ("lg_fullflag_derivative_exact.py", []),
]

NEAR_TOP = [
    ("near_top_six_weights_exact.py", []),
    ("near_top_p2_symbolic.py", []),
]

def now():
    return datetime.now(timezone.utc).isoformat()

def digest(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", choices=["all", "core", "extensions", "reviewed", "near_top"], default="all")
    ap.add_argument("--program", help="Run one named program from the selected suite")
    ap.add_argument("--output", type=Path, help="New output directory; must not exist")
    ap.add_argument("--list", action="store_true", help="List programs and arguments without executing")
    args = ap.parse_args()
    suites = {"core": CORE, "extensions": EXTENSIONS, "reviewed": REVIEWED, "near_top": NEAR_TOP}
    jobs = CORE + EXTENSIONS + REVIEWED + NEAR_TOP if args.suite == "all" else suites[args.suite]
    if args.program:
        jobs = [job for job in jobs if job[0] == args.program]
        if not jobs:
            ap.error("The requested program is not in the selected suite")
    if args.list:
        for name, extra in jobs:
            print("python " + " ".join([name] + extra))
        return 0
    if not __debug__:
        raise RuntimeError("Assertions must be enabled; do not invoke Python with -O")
    out = args.output.resolve() if args.output else BASE / ("run_" + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f"))
    source = (BASE / "checks").resolve()
    if out == source or source in out.parents:
        raise ValueError("The output directory must be outside the checks source tree")
    out.mkdir(parents=True, exist_ok=False)
    workspace = out / "work"
    manifest = {"started_at": now(), "status": "RUNNING", "suite": args.suite,
                "program_filter": args.program, "source_sha256": {},
                "scope": "Finite exact arithmetic checks; not formal verification of the unrestricted theorems."}
    records = []
    current_program = None
    def save_manifest():
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf8")
    save_manifest()
    try:
        workspace.mkdir()
        for file in sorted(source.glob("*.py")):
            shutil.copy2(file, workspace / file.name)
        files = sorted(workspace.glob("*.py"))
        if not files:
            raise RuntimeError("The checks source directory is empty")
        source_hashes = {file.name: digest(file) for file in files}
        manifest["source_sha256"] = source_hashes
        save_manifest()
        probe_code = "import json,sys,sympy; print(json.dumps({'python':sys.version,'sympy':sympy.__version__,'ignore_environment':sys.flags.ignore_environment,'optimize':sys.flags.optimize}))"
        probe = subprocess.run([sys.executable, "-E", "-X", "utf8", "-u", "-c", probe_code],
                               cwd=workspace, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, check=False,
                               creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        (out / "environment_probe.log").write_text(probe.stdout, encoding="utf8")
        if probe.returncode:
            raise RuntimeError("The child environment probe failed; inspect environment_probe.log")
        environment = json.loads(probe.stdout)
        manifest.update(environment)
        save_manifest()
        if environment["sympy"] != "1.14.0":
            raise RuntimeError("The child environment requires sympy==1.14.0; install requirements.txt")
        if environment["optimize"] != 0:
            raise RuntimeError("Assertions are disabled in the child environment")
        for index, (name, extra) in enumerate(jobs, 1):
            current_program = name
            record = {"program": name, "arguments": extra, "started_at": now()}
            started = time.monotonic()
            logfile = out / (name + ".log")
            print(f"[{index}/{len(jobs)}] {name}", flush=True)
            with logfile.open("w", encoding="utf8") as stream:
                result = subprocess.run([sys.executable, "-E", "-X", "utf8", "-u", str(workspace / name)] + extra,
                                        cwd=workspace, stdout=stream, stderr=subprocess.STDOUT,
                                        check=False, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            record.update(exit_code=result.returncode, seconds=round(time.monotonic()-started, 3),
                          ended_at=now(), status="PASS" if result.returncode == 0 else "FAIL",
                          log=logfile.name)
            records.append(record)
            (out / "summary.json").write_text(json.dumps(records, indent=2), encoding="utf8")
            if result.returncode:
                manifest.update(ended_at=now(), status="FAIL", completed_programs=len(records), failed_program=name)
                save_manifest()
                print(f"FAILED: inspect {logfile}", file=sys.stderr)
                return 1
        assert source_hashes == {file.name: digest(file) for file in files}, "A verification source changed during the run"
        manifest.update(ended_at=now(), status="PASS", completed_programs=len(records))
        save_manifest()
        print(f"All {len(records)} selected programs passed. Results: {out}", flush=True)
        return 0
    except BaseException as error:
        interrupted = isinstance(error, KeyboardInterrupt)
        manifest.update(ended_at=now(), status="INTERRUPTED" if interrupted else "FAIL",
                        error=str(error), completed_programs=len(records), failed_program=current_program)
        save_manifest()
        print(f"Reproduction stopped: {error}", file=sys.stderr)
        return 130 if interrupted else 1

if __name__ == "__main__":
    raise SystemExit(main())
