#!/usr/bin/env python3
"""Run AntiSMASH and ARTS for one FASTA, safely resuming completed work."""

import datetime
import json
import os
import shutil
import subprocess
from pathlib import Path


def _input_signature(fasta_file):
    """Signature used to detect whether this FASTA changed since processing."""
    fasta = Path(fasta_file).resolve()
    stat = fasta.stat()
    return {
        "path": str(fasta),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def _archive_incomplete(output_dir):
    """Keep an incomplete previous run for debugging, then free its output path."""
    output_dir = Path(output_dir)
    if not output_dir.exists():
        return

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    archived = output_dir.with_name(output_dir.name + ".incomplete_" + timestamp)
    shutil.move(str(output_dir), str(archived))
    print(f"Archived incomplete output: {archived}")


def _write_completion_marker(output_dir, signature, required_output):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    marker = output_dir / ".krill_complete.json"
    temporary_marker = output_dir / ".krill_complete.json.tmp"

    payload = {
        "input": signature,
        "required_output": str(Path(required_output).resolve()),
        "completed_at": datetime.datetime.now().astimezone().isoformat(),
    }

    temporary_marker.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8"
    )

    temporary_marker.replace(marker)


def _is_complete(output_dir, signature, required_output):
    """
    A result is complete only if:
      1. The output directory exists.
      2. The required output exists and is non-empty.
      3. A valid KRILL completion marker exists.
      4. The marker matches the current input FASTA.
      5. The marker points to the expected output file.
    """
    output_dir = Path(output_dir)
    required_output = Path(required_output)
    marker = output_dir / ".krill_complete.json"

    if not output_dir.is_dir():
        return False

    if not required_output.is_file() or required_output.stat().st_size == 0:
        return False

    # No marker means that completion was never confirmed.
    # This includes interrupted runs and old results without markers.
    if not marker.is_file():
        return False

    try:
        saved = json.loads(marker.read_text(encoding="utf-8"))

        if saved.get("input") != signature:
            return False

        saved_output = saved.get("required_output")

        if saved_output != str(required_output.resolve()):
            return False

        return True

    except (OSError, json.JSONDecodeError, TypeError):
        return False


def _find_coretable(arts_output_dir):
    """ARTS versions can place coretable.tsv directly or in a nested folder."""
    matches = sorted(Path(arts_output_dir).rglob("coretable.tsv"))
    return next((p for p in matches if p.is_file() and p.stat().st_size > 0), None)


def run(fasta_file, db):
    project_dir = Path(__file__).resolve().parent
    fasta_file = Path(fasta_file).resolve()
    db = Path(db).resolve()
    name = fasta_file.stem
    signature = _input_signature(fasta_file)

    antismash_dir = db / "AntiSMASH" / name
    antismash_gbk = antismash_dir / f"{name}.gbk"
    antismash_script = project_dir / "krill_run_antismash.sh"

    # ---------------------------- AntiSMASH ----------------------------
    if _is_complete(antismash_dir, signature, antismash_gbk):
        print(f"[RESUME] AntiSMASH already complete: {name}")
    else:
        if antismash_dir.exists():
            _archive_incomplete(antismash_dir)
        antismash_dir.mkdir(parents=True, exist_ok=True)

        print(f"[RUN] AntiSMASH: {name}")
        subprocess.run(
            [
                "bash",
                str(antismash_script),
                str(fasta_file),
                "1",
                "AntiSMASH",
                str(Path("AntiSMASH") / name),
            ],
            cwd=str(db),
            check=True,
        )

        if not antismash_gbk.is_file() or antismash_gbk.stat().st_size == 0:
            raise RuntimeError(
                f"AntiSMASH returned successfully but expected GenBank output "
                f"was not created: {antismash_gbk}"
            )
        _write_completion_marker(antismash_dir, signature, antismash_gbk)

    # ------------------------------- ARTS ------------------------------

    arts_dir = db / "ARTS" / name
    arts_script = project_dir / "krill_run_arts.sh"
    coretable = _find_coretable(arts_dir)

    print(f"\n[ARTS CHECK] Database: {db.name} | Genome: {name}", flush=True)
    print(f"[ARTS CHECK] Output directory: {arts_dir}", flush=True)

    if coretable is not None:
        print(f"[ARTS CHECK] Existing coretable found: {coretable}", flush=True)
        print(
            f"[ARTS CHECK] Existing coretable size: {coretable.stat().st_size} bytes",
            flush=True
        )
    else:
        print("[ARTS CHECK] No existing coretable.tsv found.", flush=True)

    if coretable is not None and _is_complete(arts_dir, signature, coretable):

        print(
            f"[ARTS SKIPPED] Existing ARTS results reused: {name}",
            flush=True
        )

    else:

        if not antismash_gbk.is_file() or antismash_gbk.stat().st_size == 0:
            raise FileNotFoundError(
                f"Cannot start ARTS because AntiSMASH GenBank output is missing: {antismash_gbk}"
            )

        if arts_dir.exists():
            print(
                f"[ARTS WARNING] Existing incomplete ARTS directory found. "
                f"Archiving: {arts_dir}",
                flush=True
            )
            _archive_incomplete(arts_dir)

        arts_dir.mkdir(parents=True, exist_ok=True)

        hmm_ref = project_dir / "db" / "bacterial_metagenomic_resist_models.hmm"

        # THIS MESSAGE CONFIRMS THAT KRILL IS ABOUT TO EXECUTE ARTS.
        print("\n" + "=" * 70, flush=True)
        print(f"[ARTS EXECUTION] STARTING ARTS: {name}", flush=True)
        print(f"[ARTS EXECUTION] Database: {db.name}", flush=True)
        print(f"[ARTS EXECUTION] Input GenBank: {antismash_gbk}", flush=True)
        print(f"[ARTS EXECUTION] Output directory: {arts_dir}", flush=True)
        print(f"[ARTS EXECUTION] HMM database: {hmm_ref}", flush=True)
        print(f"[ARTS EXECUTION] Timestamp: {datetime.datetime.now().astimezone().isoformat()}", flush=True)
        print("=" * 70 + "\n", flush=True)

        subprocess.run(
            [
                "bash",
                str(arts_script),
                str(antismash_gbk),
                str(hmm_ref),
                str(arts_dir),
                "1",
            ],
            cwd=str(db),
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

        print(
            f"[ARTS EXECUTION] ARTS command finished successfully: {name}",
            flush=True
        )

        coretable = _find_coretable(arts_dir)

        if coretable is None:
            raise RuntimeError(
                f"ARTS returned successfully but no non-empty coretable.tsv "
                f"was found under: {arts_dir}"
            )

        _write_completion_marker(arts_dir, signature, coretable)

        print(
            f"[ARTS COMPLETED] Completion marker created: {arts_dir}",
            flush=True
        )

    return int(fasta_file.stat().st_size / (10 ** 6))

