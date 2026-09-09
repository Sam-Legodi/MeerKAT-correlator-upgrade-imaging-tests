from __future__ import annotations

import contextlib
import contextvars
import os
import re
import shlex
import subprocess
import sys
import time
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Mapping, Optional, Sequence, TextIO, TypeVar


_CURRENT_AUDIT: contextvars.ContextVar[Optional["StepAudit"]] = contextvars.ContextVar(
    "mci_current_step_audit", default=None
)
_COMPLETED_RE = re.compile(r"^Completed (?P<input>.+?) ->")
_FAILED_RE = re.compile(r"^FAILED (?P<input>.+?): (?P<detail>.*)$")
_ERROR_RE = re.compile(
    r"(?:\[ERROR\]|\b(?:ERROR|CRITICAL|FATAL|SEVERE)\s*:|\bFAILED\b|"
    r"\bSkip:\s+not found\b|Traceback \(most recent call last\))",
    re.IGNORECASE,
)
_COMMAND_RE = re.compile(r"^\[[^]]+\]\s+->\s+")
_SAFE_STEP_RE = re.compile(r"[^A-Za-z0-9._-]+")
_T = TypeVar("_T")


@dataclass(frozen=True)
class InputOutcome:
    status: str
    detail: str = ""


class StepInputsFailed(RuntimeError):
    """Raised after all independent inputs have been attempted."""

    def __init__(self, step_name: str, failures: Sequence[BaseException]) -> None:
        self.step_name = step_name
        self.failures = tuple(failures)
        super().__init__(f"{step_name}: {len(failures)} input(s) failed")


class _Tee:
    def __init__(self, terminal: TextIO, logfile: TextIO) -> None:
        self._terminal = terminal
        self._logfile = logfile

    def write(self, text: str) -> int:
        self._terminal.write(text)
        self._logfile.write(text)
        return len(text)

    def flush(self) -> None:
        self._terminal.flush()
        self._logfile.flush()

    def isatty(self) -> bool:
        return self._terminal.isatty()

    @property
    def encoding(self) -> Optional[str]:
        return self._terminal.encoding


class StepAudit:
    def __init__(self, step_name: str, log_path: Path) -> None:
        self.step_name = step_name
        self.log_path = log_path
        self.inputs: "OrderedDict[str, None]" = OrderedDict()
        self.outcomes: "OrderedDict[str, InputOutcome]" = OrderedDict()
        self.error_lines: list[str] = []
        self.skipped_reasons: list[str] = []
        self.duration_seconds: Optional[float] = None

    @staticmethod
    def _normal_path(value: str) -> Optional[str]:
        if not value or " + " in value or " :: " in value:
            return None
        try:
            return os.path.normcase(str(Path(value).expanduser().resolve(strict=False)))
        except (OSError, RuntimeError, ValueError):
            return None

    def _known_key(self, value: str) -> str:
        value = value.strip().strip("'\"")
        if value in self.inputs or value in self.outcomes:
            return value
        normal = self._normal_path(value)
        if normal is not None:
            for known in (*self.inputs.keys(), *self.outcomes.keys()):
                if self._normal_path(known) == normal:
                    return known
        return value

    def register_inputs(self, inputs: Iterable[str]) -> None:
        for value in inputs:
            label = str(value)
            if label:
                self.inputs.setdefault(label, None)

    def succeeded(self, input_label: str, detail: str = "") -> None:
        key = self._known_key(input_label)
        self.inputs.setdefault(key, None)
        if self.outcomes.get(key, InputOutcome("SUCCESS")).status != "FAILED":
            self.outcomes[key] = InputOutcome("SUCCESS", detail)

    def failed(self, input_label: str, detail: str) -> None:
        key = self._known_key(input_label)
        self.inputs.setdefault(key, None)
        self.outcomes[key] = InputOutcome("FAILED", detail.strip())

    def skipped(self, reason: str) -> None:
        self.skipped_reasons.append(reason)

    def observe_line(self, line: str) -> bool:
        """Parse known per-input markers and collect error-bearing log lines."""
        clean = line.strip()
        if not clean:
            return False
        matched_outcome = False
        completed = _COMPLETED_RE.match(clean)
        if completed:
            self.succeeded(completed.group("input"))
            matched_outcome = True
        failed = _FAILED_RE.match(clean)
        if failed:
            self.failed(failed.group("input"), failed.group("detail"))
            matched_outcome = True
        error_line = (
            not _COMMAND_RE.match(clean)
            and "falling back" not in clean.lower()
            and _ERROR_RE.search(clean)
        )
        if error_line:
            if clean not in self.error_lines:
                self.error_lines.append(clean)
            if failed is None:
                for input_label in tuple(self.inputs):
                    parts = [part.strip() for part in input_label.split(" + ")]
                    if any(part and part in clean for part in parts):
                        self.failed(input_label, clean)
        return matched_outcome

    def audit_log(self) -> None:
        """Read the completed step log and extract outcomes and error lines."""
        try:
            text = self.log_path.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            self.error_lines.append(f"Could not read audit log: {exc}")
            return
        for line in text.splitlines():
            self.observe_line(line)

    def _summary_lines(self, exception: Optional[BaseException]) -> list[str]:
        succeeded = [
            (label, outcome)
            for label, outcome in self.outcomes.items()
            if outcome.status == "SUCCESS"
        ]
        failed = [
            (label, outcome)
            for label, outcome in self.outcomes.items()
            if outcome.status == "FAILED"
        ]
        not_run = [label for label in self.inputs if label not in self.outcomes]
        state = "FAILED" if exception is not None or failed or self.error_lines else "SUCCESS"
        if self.skipped_reasons and not self.inputs and exception is None:
            state = "SKIPPED"

        lines = [
            "",
            f"[AUDIT] Step: {self.step_name}",
            f"[AUDIT] Status: {state}",
            (
                "[AUDIT] Inputs: "
                f"{len(succeeded)} succeeded, {len(failed)} failed, {len(not_run)} not run"
            ),
        ]
        if self.duration_seconds is not None:
            seconds = self.duration_seconds
            hours, remainder = divmod(int(seconds), 3600)
            minutes, whole_seconds = divmod(remainder, 60)
            lines.insert(3, f"[AUDIT] Duration: {hours:02d}:{minutes:02d}:{whole_seconds:02d} ({seconds:.3f} seconds)")
        if succeeded:
            lines.append("[AUDIT] Succeeded inputs:")
            lines.extend(f"  SUCCESS  {label}" for label, _ in succeeded)
        if failed:
            lines.append("[AUDIT] Failed inputs:")
            for label, outcome in failed:
                suffix = f" :: {outcome.detail}" if outcome.detail else ""
                lines.append(f"  FAILED   {label}{suffix}")
        if not_run:
            lines.append("[AUDIT] Inputs not run:")
            lines.extend(f"  NOT RUN  {label}" for label in not_run)
        if self.skipped_reasons:
            lines.append("[AUDIT] Skip reason(s):")
            lines.extend(f"  {reason}" for reason in self.skipped_reasons)
        if exception is not None:
            lines.append(f"[AUDIT] Step exception: {type(exception).__name__}: {exception}")
        if self.error_lines:
            lines.append(f"[AUDIT] Error-bearing log lines: {len(self.error_lines)}")
            lines.extend(f"  {line}" for line in self.error_lines[:20])
            if len(self.error_lines) > 20:
                lines.append(f"  ... {len(self.error_lines) - 20} additional line(s) in the log")
        lines.append(f"[AUDIT] Log: {self.log_path}")
        return lines

    def print_summary(self, exception: Optional[BaseException]) -> None:
        print("\n".join(self._summary_lines(exception)))


def current_audit() -> Optional[StepAudit]:
    return _CURRENT_AUDIT.get()


def register_inputs(inputs: Iterable[str]) -> None:
    audit = current_audit()
    if audit is not None:
        audit.register_inputs(inputs)


def record_success(input_label: str, detail: str = "") -> None:
    audit = current_audit()
    if audit is not None:
        audit.succeeded(input_label, detail)


def record_failure(input_label: str, error: BaseException | str) -> None:
    audit = current_audit()
    if audit is not None:
        audit.failed(input_label, str(error))


def record_skip(reason: str) -> None:
    audit = current_audit()
    if audit is not None:
        audit.skipped(reason)


def raise_for_failures(step_name: str, failures: Sequence[BaseException]) -> None:
    if failures:
        raise StepInputsFailed(step_name, failures)


def run_logged_command(
    command: Iterable[str],
    *,
    prefix: str,
    inputs: Iterable[str] = (),
    env: Optional[Mapping[str, str]] = None,
    mark_success: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a command, tee its output, and update the current per-step audit."""
    command = [str(token) for token in command]
    input_labels = [str(value) for value in inputs]
    audit = current_audit()
    if audit is not None:
        audit.register_inputs(input_labels)

    print(f"{prefix} ->", " ".join(shlex.quote(token) for token in command))
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=dict(env) if env is not None else None,
        )
    except OSError as exc:
        print(f"{prefix} ERROR: could not start command: {exc}")
        if audit is not None:
            for input_label in input_labels:
                audit.failed(input_label, str(exc))
            audit.error_lines.append(f"could not start command: {exc}")
        raise
    assert process.stdout is not None
    saw_outcome_marker = False
    for line in process.stdout:
        print(line, end="")
        if audit is not None:
            saw_outcome_marker = audit.observe_line(line) or saw_outcome_marker
    returncode = process.wait()

    if audit is not None:
        if returncode == 0 and mark_success:
            for input_label in input_labels:
                audit.succeeded(input_label)
        elif returncode != 0:
            detail = f"command exited with status {returncode}"
            # A child such as PyBDSF can report mixed per-input outcomes. Only
            # fill missing records; preserve the explicit Completed/FAILED lines.
            for input_label in input_labels:
                key = audit._known_key(input_label)
                if key not in audit.outcomes:
                    audit.failed(input_label, detail)
            if not input_labels and not saw_outcome_marker:
                audit.error_lines.append(detail)

    completed = subprocess.CompletedProcess(command, returncode)
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, command)
    return completed


def run_step_with_audit(
    step_name: str,
    reports_dir: str | Path,
    runner: Callable[[], _T],
) -> _T:
    """Run one pipeline step, then audit and summarize its combined log."""
    audit_dir = Path(reports_dir).expanduser() / "pipeline_audits"
    audit_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    safe_step = _SAFE_STEP_RE.sub("_", step_name).strip("_") or "step"
    log_path = audit_dir / f"{timestamp}_{safe_step}.log"
    audit = StepAudit(step_name, log_path)
    exception: Optional[BaseException] = None
    result: Optional[_T] = None

    with log_path.open("w", encoding="utf-8") as logfile:
        stdout_tee = _Tee(sys.stdout, logfile)
        stderr_tee = _Tee(sys.stderr, logfile)
        token = _CURRENT_AUDIT.set(audit)
        try:
            with contextlib.redirect_stdout(stdout_tee), contextlib.redirect_stderr(stderr_tee):
                started = time.monotonic()
                print(f"[AUDIT] Starting step: {step_name}", flush=True)
                try:
                    result = runner()
                except BaseException as exc:
                    exception = exc
                finally:
                    audit.duration_seconds = time.monotonic() - started
                logfile.flush()
                audit.audit_log()
                audit.print_summary(exception)
        finally:
            _CURRENT_AUDIT.reset(token)

    if exception is not None:
        raise exception
    return result  # type: ignore[return-value]
