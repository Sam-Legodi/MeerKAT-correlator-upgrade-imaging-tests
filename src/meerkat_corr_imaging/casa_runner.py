"""Finite CASA subprocess supervision; no silence-based timeouts."""
from __future__ import annotations
import json
import os
from pathlib import Path
import selectors
import signal
import subprocess
import time
import uuid


def supervise(command, env, timeout=21600, shutdown_timeout=30, heartbeat=30):
    process = subprocess.Popen(command, env=env, stdin=subprocess.DEVNULL,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               start_new_session=True)
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ)
    start = last = time.monotonic()
    completed_at = None
    stage_started = None
    stage = "startup"
    pending = ""
    stage_timeout = float(env.get("MCI_STAGE_TIMEOUT_SECONDS", timeout))
    finished = False
    try:
        while selector.get_map() or process.poll() is None:
            now = time.monotonic()
            result_path = env.get('MCI_RESULT_PATH')
            if result_path and Path(result_path).exists() and completed_at is None:
                completed_at = now
            if stage_started is not None and now-stage_started > stage_timeout:
                raise TimeoutError("CASA stage exceeded deadline: " + stage)
            if now-start > timeout or (completed_at is not None and now-completed_at > shutdown_timeout):
                raise TimeoutError('CASA exceeded runtime/shutdown deadline')
            if now-last >= heartbeat:
                print(f'[CASA] {stage}; elapsed {now-start:.0f}s', flush=True)
                last = now
            for key, _ in selector.select(timeout=0.2):
                data = os.read(key.fd, 65536)
                if data:
                    decoded = data.decode('utf-8', errors='replace')
                    print(decoded, end='', flush=True)
                    pending += decoded
                    while '\n' in pending:
                        line, pending = pending.split('\n', 1)
                        if line.startswith('[CAL TASK] Starting '):
                            stage = line.split('Starting ', 1)[1]
                            stage_started = time.monotonic()
                        elif line.startswith('[CAL TASK] Finished '):
                            stage_started = None
                            stage = 'between tasks'
                    pending = pending[-65536:]
                else:
                    selector.unregister(key.fileobj)
        if process.wait() != 0:
            raise RuntimeError(f'CASA exited with status {process.returncode}')
        finished = True
    finally:
        selector.close()
        if not finished or process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        process.stdout.close()


def run_calibration(command, env, reports_dir, timeout, shutdown_timeout):
    env = dict(env)
    run_id = uuid.uuid4().hex
    directory = Path(reports_dir) / 'calibration_results'
    directory.mkdir(parents=True, exist_ok=True)
    path = directory.resolve() / (run_id + '.json')
    env.update(MCI_RUN_ID=run_id, MCI_RESULT_PATH=str(path))
    try:
        supervise(command, env, timeout, shutdown_timeout)
    except Exception as exc:
        # A failed task may be followed by CASA shutdown/pipe trouble. Preserve
        # its authenticated result rather than replacing it with a timeout.
        if path.is_file():
            try:
                failure = json.loads(path.read_text())
            except (ValueError, OSError):
                failure = {}
            if (failure.get('run_id') == run_id and
                    failure.get('ms') == env['MCI_CAL_MSFILE'] and failure.get('error')):
                raise RuntimeError('CASA calibration failed: ' + str(failure['error']) +
                                   '; process supervision: ' + str(exc)) from exc
        raise
    if not path.is_file():
        raise RuntimeError('CASA exited without a calibration result')
    result = json.loads(path.read_text())
    if result.get('run_id') != run_id or result.get('ms') != env['MCI_CAL_MSFILE']:
        raise RuntimeError('CASA result identity mismatch')
    if not result.get('calibration_applied') or result.get('error'):
        raise RuntimeError('CASA calibration failed: ' + str(result.get('error', 'application not validated')))
    if env.get('MCI_QUALITY_ENABLED') == '1' and result.get('quality') != 'passed':
        raise RuntimeError('Requested calibration quality checks did not pass')
    print('[CAL/IMAGING] Validated result: ' + str(path), flush=True)
    return result
