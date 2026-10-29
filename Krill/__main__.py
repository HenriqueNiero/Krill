
import sys
import os
import pathlib
import argparse
import datetime
import subprocess
import json
import csv
import time
import threading
import platform
import socket
import re
import atexit
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor, as_completed

from tqdm import tqdm
from cprint import *

import prepareFastas
import ARTS_extractor
import AntiSMASH_extractor
import get_screening_results
import count_bases_and_ORFs
import run_AntiSMASH_and_ARTS
import get_dbs_metainfo
import build_charts
import run_bigscape
import run_bigscape_analysis
from bigscape_integration import integrate_bigscape_with_krill


# ---------------------------------------------------------------------
# Performance monitoring
# ---------------------------------------------------------------------

def read_linux_processes():
    """
    Read process information from /proc.

    Returns a dictionary keyed by PID. Each entry contains:
    PPID, CPU time in clock ticks, and resident memory in bytes.
    Linux-specific; unavailable process entries are ignored.
    """
    processes = {}

    try:
        entries = os.listdir("/proc")
    except OSError:
        return processes

    page_size = os.sysconf("SC_PAGE_SIZE")

    try:
        clock_ticks = os.sysconf("SC_CLK_TCK")
    except (ValueError, OSError):
        clock_ticks = 100

    for entry in entries:
        if not entry.isdigit():
            continue

        pid = int(entry)

        try:
            with open(f"/proc/{pid}/stat", "r") as handle:
                stat_line = handle.read()

            # The process name can contain spaces and parentheses.
            closing_parenthesis = stat_line.rfind(")")
            if closing_parenthesis == -1:
                continue

            fields = stat_line[closing_parenthesis + 2:].split()

            # Fields are indexed relative to field 3 of /proc/PID/stat.
            ppid = int(fields[1])
            utime = int(fields[11])
            stime = int(fields[12])
            rss_pages = int(fields[21])

            processes[pid] = {
                "ppid": ppid,
                "cpu_ticks": utime + stime,
                "rss_bytes": max(0, rss_pages) * page_size,
                "clock_ticks": clock_ticks,
            }

        except (OSError, ValueError, IndexError, ProcessLookupError):
            # A process may finish while /proc is being read.
            continue

    return processes


def get_process_tree(root_pid):
    """Return the PIDs belonging to root_pid and its descendants."""
    processes = read_linux_processes()

    if root_pid not in processes:
        return {root_pid}, processes

    children = {}

    for pid, info in processes.items():
        children.setdefault(info["ppid"], []).append(pid)

    tree = set()
    pending = [root_pid]

    while pending:
        pid = pending.pop()

        if pid in tree:
            continue

        tree.add(pid)
        pending.extend(children.get(pid, []))

    return tree, processes


def get_process_tree_usage(root_pid):
    """
    Return total CPU ticks, summed resident memory, and process count
    for the current process tree.

    Summed RSS can count shared memory more than once.
    """
    tree, processes = get_process_tree(root_pid)

    cpu_ticks = 0
    rss_bytes = 0
    count = 0
    clock_ticks = os.sysconf("SC_CLK_TCK")

    for pid in tree:
        info = processes.get(pid)

        if info is None:
            continue

        cpu_ticks += info["cpu_ticks"]
        rss_bytes += info["rss_bytes"]
        count += 1
        clock_ticks = info["clock_ticks"]

    return {
        "cpu_ticks": cpu_ticks,
        "rss_bytes": rss_bytes,
        "process_count": count,
        "clock_ticks": clock_ticks,
    }


def get_system_information(threads):
    """Collect hardware, operating-system, Python, and SLURM details."""
    info = {
        "timestamp": datetime.datetime.now().astimezone().isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "operating_system": platform.system(),
        "architecture": platform.machine(),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "krill_pid": os.getpid(),
        "requested_krill_threads": threads,
        "cpu_count_os": os.cpu_count(),
        "cpu_affinity_count": None,
        "cpu_affinity_list": None,
        "cpu_model": None,
        "total_memory_gb": None,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_job_name": os.environ.get("SLURM_JOB_NAME"),
        "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "slurm_ntasks": os.environ.get("SLURM_NTASKS"),
        "slurm_mem_per_cpu": os.environ.get("SLURM_MEM_PER_CPU"),
        "slurm_mem_per_node": os.environ.get("SLURM_MEM_PER_NODE"),
        "slurm_partition": os.environ.get("SLURM_JOB_PARTITION"),
        "container_indicators": [],
    }

    if hasattr(os, "sched_getaffinity"):
        try:
            affinity = sorted(os.sched_getaffinity(0))
            info["cpu_affinity_count"] = len(affinity)
            info["cpu_affinity_list"] = affinity
        except OSError:
            pass

    # Linux CPU model.
    try:
        with open("/proc/cpuinfo", "r") as handle:
            cpuinfo = handle.read()

        match = re.search(r"^model name\s*:\s*(.+)$", cpuinfo, re.MULTILINE)

        if match:
            info["cpu_model"] = match.group(1).strip()

    except OSError:
        pass

    # Linux total physical memory visible in /proc/meminfo.
    # In a container, this may differ from the memory limit assigned
    # to the container or SLURM job.
    try:
        with open("/proc/meminfo", "r") as handle:
            meminfo = handle.read()

        match = re.search(r"^MemTotal:\s+(\d+)\s+kB$", meminfo, re.MULTILINE)

        if match:
            info["total_memory_gb"] = round(
                int(match.group(1)) * 1024 / (1024 ** 3), 3
            )

    except OSError:
        pass

    if os.path.exists("/.dockerenv"):
        info["container_indicators"].append("Docker indicator: /.dockerenv")

    try:
        with open("/proc/1/cgroup", "r") as handle:
            cgroup = handle.read()

        if "singularity" in cgroup.lower() or "apptainer" in cgroup.lower():
            info["container_indicators"].append("Singularity/Apptainer cgroup indicator")

    except OSError:
        pass

    if os.path.exists("/.singularity.d"):
        info["container_indicators"].append("Singularity indicator: /.singularity.d")

    if os.path.exists("/.apptainer.d"):
        info["container_indicators"].append("Apptainer indicator: /.apptainer.d")

    return info


class PerformanceMonitor:
    """
    Optional KRILL performance monitor.

    Outputs:
      PerformanceReport/system_info.json
      PerformanceReport/stage_metrics.csv
      PerformanceReport/resource_samples.csv
      PerformanceReport/performance_summary.json
    """

    def __init__(self, output_path, enabled=False, interval=1.0, threads=1):
        self.enabled = enabled
        self.interval = interval
        self.threads = threads
        self.output_path = pathlib.Path(output_path) / "PerformanceReport"

        self.started_at = None
        self.finished_at = None
        self.start_monotonic = None
        self.sampler_thread = None
        self.stop_event = threading.Event()
        self.lock = threading.Lock()

        self.stage_rows = []
        self.sample_rows = []
        self.active_stage = None
        self.active_stage_started = None
        self.active_stage_row = None
        self.stage_number = 0
        self.last_cpu_ticks = None
        self.last_sample_time = None
        self.peak_rss_mb = 0.0
        self.peak_cpu_percent_cores = 0.0
        self.peak_cpu_percent_normalized = 0.0
        self.max_process_count = 0
        self.error_count = 0
        self.system_info = None
        self._finished = False

    def start(self):
        if not self.enabled:
            return

        self.output_path.mkdir(parents=True, exist_ok=True)

        self.started_at = datetime.datetime.now().astimezone()
        self.start_monotonic = time.monotonic()

        self.system_info = get_system_information(self.threads)

        self.system_info["monitoring_interval_seconds"] = self.interval
        self.system_info["performance_report_directory"] = str(self.output_path.resolve())

        self._write_json(
            self.output_path / "system_info.json",
            self.system_info
        )

        # Establish the CPU baseline before the first sampling interval.
        usage = get_process_tree_usage(os.getpid())
        self.last_cpu_ticks = usage["cpu_ticks"]
        self.last_sample_time = time.monotonic()

        self.sampler_thread = threading.Thread(
            target=self._sample_loop,
            name="krill-performance-monitor",
            daemon=True,
        )
        self.sampler_thread.start()

        print(f"[Performance] Monitoring enabled.")
        print(f"[Performance] Report directory: {self.output_path.resolve()}")
        print(f"[Performance] Sampling interval: {self.interval:g} second(s)")

    def _sample_loop(self):
        while not self.stop_event.wait(self.interval):
            try:
                self._take_sample()
            except Exception as exc:
                # Monitoring errors must not stop the biological pipeline.
                print(f"[Performance] Warning: resource sample failed: {exc}")

    def _take_sample(self):
        now = time.monotonic()
        usage = get_process_tree_usage(os.getpid())

        elapsed = now - self.start_monotonic
        interval_elapsed = now - self.last_sample_time

        if interval_elapsed <= 0:
            interval_elapsed = self.interval

        tick_delta = max(0, usage["cpu_ticks"] - self.last_cpu_ticks)
        cpu_seconds = tick_delta / usage["clock_ticks"]

        # CPU percentage in core-equivalents: 100% means one fully
        # occupied CPU core; multithreaded jobs can exceed 100%.
        cpu_percent_cores = 100.0 * cpu_seconds / interval_elapsed

        # Normalize to the requested KRILL thread count.
        cpu_percent_normalized = (
            cpu_percent_cores / self.threads
            if self.threads > 0
            else 0.0
        )

        rss_mb = usage["rss_bytes"] / (1024 ** 2)

        row = {
            "timestamp": datetime.datetime.now().astimezone().isoformat(),
            "elapsed_seconds": round(elapsed, 3),
            "stage": self.active_stage or "Pipeline",
            "cpu_percent_cores": round(cpu_percent_cores, 2),
            "cpu_percent_of_requested_threads": round(
                cpu_percent_normalized, 2
            ),
            "process_tree_rss_mb": round(rss_mb, 2),
            "process_count": usage["process_count"],
        }

        with self.lock:
            self.sample_rows.append(row)
            self.peak_rss_mb = max(self.peak_rss_mb, rss_mb)
            self.peak_cpu_percent_cores = max(
                self.peak_cpu_percent_cores, cpu_percent_cores
            )
            self.peak_cpu_percent_normalized = max(
                self.peak_cpu_percent_normalized, cpu_percent_normalized
            )
            self.max_process_count = max(
                self.max_process_count, usage["process_count"]
            )

        self.last_cpu_ticks = usage["cpu_ticks"]
        self.last_sample_time = now

    def stage_start(self, name):
        if not self.enabled:
            return

        with self.lock:
            # Close a previously open stage if control flow did not
            # close it normally.
            if self.active_stage is not None:
                self._close_active_stage("interrupted")

            self.stage_number += 1
            self.active_stage = name
            self.active_stage_started = time.monotonic()
            self.active_stage_row = {
                "stage_number": self.stage_number,
                "stage": name,
                "start_timestamp": datetime.datetime.now().astimezone().isoformat(),
                "duration_seconds": None,
                "status": "running",
                "error": "",
            }

    def stage_end(self, status="success", error=""):
        if not self.enabled:
            return

        with self.lock:
            self._close_active_stage(status, error)

    def _close_active_stage(self, status, error=""):
        if self.active_stage_row is None:
            return

        duration = time.monotonic() - self.active_stage_started

        self.active_stage_row["duration_seconds"] = round(duration, 3)
        self.active_stage_row["end_timestamp"] = (
            datetime.datetime.now().astimezone().isoformat()
        )
        self.active_stage_row["status"] = status
        self.active_stage_row["error"] = str(error)[:1000]

        self.stage_rows.append(self.active_stage_row)

        if status != "success":
            self.error_count += 1

        self.active_stage = None
        self.active_stage_started = None
        self.active_stage_row = None

    @contextmanager
    def stage(self, name):
        """Context manager for measuring a pipeline stage."""
        self.stage_start(name)

        try:
            yield
        except Exception as exc:
            self.stage_end(status="failed", error=repr(exc))
            raise
        else:
            self.stage_end(status="success")

    def run_stage(self, name, function, *args, **kwargs):
        """Run a function while recording its duration and status."""
        with self.stage(name):
            return function(*args, **kwargs)

    @staticmethod
    def _write_json(path, data):
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, default=str)

    @staticmethod
    def _write_csv(path, rows, fieldnames):
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def finish(self, final_status="success", error=""):
        if not self.enabled or self._finished:
            return

        self._finished = True
        self.finished_at = datetime.datetime.now().astimezone()

        # Stop sampling and wait for the sampler thread to finish.
        self.stop_event.set()

        if self.sampler_thread and self.sampler_thread.is_alive():
            self.sampler_thread.join(timeout=max(2.0, self.interval + 1.0))

        # Capture a final sample even if the run finished before one
        # complete sampling interval elapsed.
        try:
            self._take_sample()
        except Exception as exc:
            print(f"[Performance] Warning: final resource sample failed: {exc}")

        with self.lock:
            if self.active_stage is not None:
                self._close_active_stage(
                    "failed" if final_status != "success" else "interrupted",
                    error,
                )

            stages = list(self.stage_rows)
            samples = list(self.sample_rows)

        total_seconds = (
            time.monotonic() - self.start_monotonic
            if self.start_monotonic is not None
            else 0.0
        )

        stage_fields = [
            "stage_number",
            "stage",
            "start_timestamp",
            "end_timestamp",
            "duration_seconds",
            "status",
            "error",
        ]

        sample_fields = [
            "timestamp",
            "elapsed_seconds",
            "stage",
            "cpu_percent_cores",
            "cpu_percent_of_requested_threads",
            "process_tree_rss_mb",
            "process_count",
        ]

        try:
            self._write_csv(
                self.output_path / "stage_metrics.csv",
                stages,
                stage_fields,
            )

            self._write_csv(
                self.output_path / "resource_samples.csv",
                samples,
                sample_fields,
            )

            stage_durations = sorted(
                stages,
                key=lambda row: row["duration_seconds"] or 0,
                reverse=True,
            )

            summary = {
                "start_timestamp": (
                    self.started_at.isoformat() if self.started_at else None
                ),
                "end_timestamp": self.finished_at.isoformat(),
                "total_runtime_seconds": round(total_seconds, 3),
                "total_runtime_minutes": round(total_seconds / 60, 3),
                "final_status": final_status,
                "final_error": str(error)[:1000],
                "requested_krill_threads": self.threads,
                "sampling_interval_seconds": self.interval,
                "number_of_samples": len(samples),
                "number_of_stages": len(stages),
                "number_of_failed_or_interrupted_stages": self.error_count,
                "peak_process_tree_rss_mb": round(self.peak_rss_mb, 2),
                "peak_cpu_percent_cores": round(
                    self.peak_cpu_percent_cores, 2
                ),
                "peak_cpu_percent_of_requested_threads": round(
                    self.peak_cpu_percent_normalized, 2
                ),
                "maximum_process_count": self.max_process_count,
                "slowest_stages": [
                    {
                        "stage": row["stage"],
                        "duration_seconds": row["duration_seconds"],
                        "status": row["status"],
                    }
                    for row in stage_durations[:10]
                ],
                "measurement_notes": [
                    "CPU percentage in core-equivalents may exceed 100% for multithreaded workloads.",
                    "CPU percentage normalized to requested threads is an estimate, not a scheduler-level utilization measurement.",
                    "Memory is the sum of resident set sizes (RSS) for the KRILL process and descendants; shared memory can be counted more than once.",
                    "Short sampling intervals can miss brief CPU or memory peaks.",
                    "System memory from /proc/meminfo may describe the host rather than the SLURM job limit.",
                    "External tools are included only while they are descendants of the monitored KRILL process.",
                ],
            }

            self._write_json(
                self.output_path / "performance_summary.json",
                summary,
            )

            print("\n[Performance] Monitoring report completed.")
            print(f"[Performance] Total runtime: {total_seconds / 60:.2f} minutes")
            print(f"[Performance] Peak process-tree RSS: {self.peak_rss_mb:.2f} MB")
            print(f"[Performance] Peak CPU: {self.peak_cpu_percent_cores:.2f}% of one CPU core")
            print(f"[Performance] Reports: {self.output_path.resolve()}")

        except Exception as exc:
            # Do not hide the pipeline result because a report could not
            # be written, for example if the filesystem is full.
            print(f"[Performance] Warning: could not write report: {exc}")


# ---------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------

default_threads = len(os.sched_getaffinity(0)) if hasattr(
    os, "sched_getaffinity"
) else (os.cpu_count() or 1)

cprint.ok('#################\n### Krill 1.0 ###\n#################')
cprint.ok('\nPlease cite us!\n')

parser = argparse.ArgumentParser()

parser.add_argument(
    '-noprep',
    '--do_not_prepare_fasta_files',
    help='Rename fasta files, its headers and store changes in a CSV file for control [DEFAULT: TRUE]',
    action='store_true'
)

parser.add_argument(
    '-t',
    '--threads',
    help='Threads to use in analysis [DEFAULT: {}]'.format(default_threads),
    type=int,
    default=default_threads
)

parser.add_argument(
    "--bigscape",
    action="store_true",
    help="Run BiG-SCAPE2 after antiSMASH"
)

parser.add_argument(
    "--bigscape_cutoff",
    default="0.3",
    help="BiG-SCAPE cutoff (default 0.3)"
)

parser.add_argument(
    'PATH',
    help='Working path with fasta files',
    type=str
)

parser.add_argument(
    "--pfam-path",
    dest="pfam",
    default=None,
    help="Path to Pfam-A.hmm"
)

parser.add_argument(
    "--bigscape-analysis",
    action="store_true",
    default=False,
    help="Run the optional BiG-SCAPE category analysis on an existing BiG-SCAPE database"
)

parser.add_argument(
    "--bigscape-analysis-threshold",
    type=float,
    default=0.3,
    help="Distance threshold used for BiG-SCAPE singleton classification [DEFAULT: 0.3]"
)

parser.add_argument(
    "--bigscape-analysis-permutations",
    type=int,
    default=999,
    help="Number of permutations used for PERMANOVA/PERMDISP [DEFAULT: 999]"
)

parser.add_argument(
    "--krill_plots",
    action="store_true",
    default=False,
    help="Generate all Krill plots after completing the analysis"
)

# New performance-monitoring options.
parser.add_argument(
    "--performance-monitor",
    action="store_true",
    help="Monitor KRILL runtime, CPU, memory, and pipeline stage durations"
)

parser.add_argument(
    "--performance-interval",
    type=float,
    default=1.0,
    help="Resource sampling interval in seconds [DEFAULT: 1]"
)

args = parser.parse_args()

if args.threads < 1:
    parser.error("--threads must be at least 1")

if args.performance_interval <= 0:
    parser.error("--performance-interval must be greater than 0")

krill_input_path = pathlib.Path(args.PATH).resolve()

# Fixed BiG-SCAPE settings.
args.bigscape_env = "bigscape"
args.bigscape_mix = True
args.bigscape_classify = "category"
args.bigscape_include_singletons = True

root_path = os.path.abspath(args.PATH)
root_database = False
dbs = []

# Find database directories inside the input directory.
# Files such as taxonomy_krill are intentionally ignored.
for x in os.listdir(root_path):
    full_path = os.path.join(root_path, x)

    if os.path.isdir(full_path):
        dbs.append(os.path.abspath(full_path))

# If no database directories exist, the input directory is a single database.
if not dbs:
    root_database = True
    dbs.append(root_path)

args.PATH = root_path

monitor = PerformanceMonitor(
    output_path=args.PATH,
    enabled=args.performance_monitor,
    interval=args.performance_interval,
    threads=args.threads,
)

# Fallback in case the pipeline exits unexpectedly.
atexit.register(monitor.finish, "interrupted", "KRILL exited before normal completion")

pipeline_status = "success"
pipeline_error = ""

try:
    monitor.start()

    # Convert all file extensions (fa, fna, etc.) to "fasta".
    monitor.run_stage(
        "Convert FASTA extensions",
        prepareFastas.convert2fasta,
        args.PATH,
    )

    # Calculate total size of data to be analyzed.
    fastas = list(pathlib.Path(args.PATH).rglob("*.fasta"))
    total_size = int(sum(os.stat(x).st_size for x in fastas) / (10 ** 6))

    # Display a visual progress bar for the processing of the whole execution.
    with tqdm(
        desc='Total Databases Analysis progress in MB',
        total=total_size,
        unit='MB',
        colour='blue',
        position=0,
        leave=False
    ) as dbar:

        for db in dbs:
            db_name = os.path.basename(db)
            os.chdir(db)

            if not args.do_not_prepare_fasta_files:
                monitor.run_stage(
                    f"{db_name}: Prepare FASTA files",
                    prepareFastas.run,
                    db,
                    "fasta",
                )

            os.makedirs(os.path.join(db, 'AntiSMASH'), exist_ok=True)
            os.makedirs(os.path.join(db, 'ARTS'), exist_ok=True)
            os.makedirs(os.path.join(db, 'ReportOutput'), exist_ok=True)

            db_fastas = list(pathlib.Path(db).glob("*.fasta"))

            # Run AntiSMASH and ARTS concurrently for each FASTA file.
            with monitor.stage(f"{db_name}: AntiSMASH and ARTS jobs"):
                with tqdm(
                    desc=f'{db_name} BGCs Analysis',
                    total=len(db_fastas),
                    unit='file',
                    colour='blue',
                    position=1,
                    leave=False
                ) as pbar:

                    with ThreadPoolExecutor(max_workers=args.threads) as executor:
                        result_futures = [
                            executor.submit(
                                run_AntiSMASH_and_ARTS.run,
                                fasta_file,
                                db,
                            )
                            for fasta_file in db_fastas
                        ]

                        for future in as_completed(result_futures):
                            try:
                                result = future.result()
                                pbar.update(1)
                                dbar.update(int(result))
                            except Exception as exc:
                                cprint.fatal(exc, interrupt=False)
                                raise

            cprint.info('# Extracting AntiSMASH and ARTS results...')

            monitor.run_stage(
                f"{db_name}: Prepare AntiSMASH extraction",
                AntiSMASH_extractor.prepare_extraction,
                db,
            )

            monitor.run_stage(
                f"{db_name}: Prepare ARTS extraction",
                ARTS_extractor.prepare_extraction,
                db,
            )

            monitor.run_stage(
                f"{db_name}: Parse AntiSMASH protoclusters",
                AntiSMASH_extractor.protocluster_parse,
                db,
                args.threads,
            )

            monitor.run_stage(
                f"{db_name}: BLAST cluster analysis",
                AntiSMASH_extractor.get_clusters_blast,
                db,
                args.threads,
            )

            monitor.run_stage(
                f"{db_name}: ARTS overview",
                ARTS_extractor.ARTS_overview,
                db,
            )

            monitor.run_stage(
                f"{db_name}: Extract ARTS results",
                ARTS_extractor.ARTS_Results_Extraction,
                db,
            )

            cprint.info('# Screening AntiSMASH and ARTS results...')

            monitor.run_stage(
                f"{db_name}: Screening results",
                get_screening_results.run,
                db,
            )

            cprint.info('# Getting fasta info (Bases and ORFs count)...')

            monitor.run_stage(
                f"{db_name}: Count bases and ORFs",
                count_bases_and_ORFs.run,
                db,
                "fasta",
                args.threads,
            )

            os.chdir(args.PATH)

    os.makedirs(os.path.join(args.PATH, 'DBsReportOutput'), exist_ok=True)

    cprint.info(
        '# Getting metainfo (Databases MBases, ORFs, BGCs/Mbases and BGCs/ORFs)...'
    )

    monitor.run_stage(
        "Database metadata report",
        get_dbs_metainfo.get,
        args.PATH,
        "fasta",
        root_database,
    )

    if args.bigscape:
        cprint.info("# Running global BiG-SCAPE2 analysis...")

        monitor.run_stage(
            "BiG-SCAPE analysis",
            run_bigscape.run,
            path=args.PATH,
            threads=args.threads,
            cutoff=args.bigscape_cutoff,
            env=args.bigscape_env,
            pfam=args.pfam,
        )

        cprint.info("# Integrating BiG-SCAPE results with Krill...")

        monitor.run_stage(
            "Integrate BiG-SCAPE results with KRILL",
            integrate_bigscape_with_krill,
            path=args.PATH,
            cutoff=args.bigscape_cutoff,
        )

    if args.bigscape_analysis:
        cprint.info("# Running optional category analysis for BiG-SCAPE results...")

        monitor.run_stage(
            "BiG-SCAPE category analysis",
            run_bigscape_analysis.run,
            path=args.PATH,
            env=args.bigscape_env,
            threshold=args.bigscape_analysis_threshold,
            permutations=args.bigscape_analysis_permutations,
        )

    cprint.info('# Building charts...')

    monitor.run_stage(
        "Build KRILL charts",
        build_charts.build_charts,
        args.PATH,
    )

    if args.krill_plots:
        cprint.info('# Generating all Krill plots...')

        plots_input = (
            krill_input_path
            / "DBsReportOutput"
            / "DBs_BGCs_with_Hits_BiGSCAPE.tsv"
        )

        plots_script = (
            pathlib.Path(__file__).resolve().parent
            / "krill_all_plots.py"
        )

        if not plots_script.is_file():
            raise FileNotFoundError(
                f"Could not find the plotting script: {plots_script}"
            )

        if not plots_input.is_file():
            raise FileNotFoundError(
                "Cannot generate Krill plots because the integrated "
                f"BiG-SCAPE table was not found: {plots_input}\n"
                "Run Krill with --bigscape to generate the required table."
            )

        plots_output = krill_input_path / "KrillPlots"
        plots_output.mkdir(parents=True, exist_ok=True)

        command = [
            sys.executable,
            str(plots_script),
            str(plots_input),
            "-o",
            str(plots_output),
            "--input-path",
            str(krill_input_path),
        ]

        cprint.info("Running krill_all_plots.py...")
        cprint.info(f"Input table: {plots_input}")
        cprint.info(f"Output directory: {plots_output}")

        monitor.run_stage(
            "Generate KRILL plots",
            subprocess.run,
            command,
            check=True,
        )

        cprint.ok("All Krill plots generated successfully.")

except BaseException as exc:
    pipeline_status = "failed"
    pipeline_error = repr(exc)
    raise

finally:
    monitor.finish(pipeline_status, pipeline_error)

cprint.ok(
    '\nAll done. Any questions please contact: '
    'henrique.niero@lnbio.cnpem.br \nCheers!'
)
