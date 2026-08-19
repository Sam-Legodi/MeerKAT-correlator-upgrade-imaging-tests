from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from meerkat_corr_imaging import cli
from meerkat_corr_imaging.audit import (
    StepInputsFailed,
    record_success,
    run_logged_command,
    run_step_with_audit,
)
from meerkat_corr_imaging.config import Config, PathsCfg, Target
from meerkat_corr_imaging.steps import step2_vis_analysis


def _only_log(reports_dir: Path) -> Path:
    logs = list((reports_dir / "pipeline_audits").glob("*.log"))
    assert len(logs) == 1
    return logs[0]


def test_step_audit_reports_mixed_child_outcomes(tmp_path: Path) -> None:
    reports_dir = tmp_path / "reports"
    child = (
        "import sys; "
        "print('Completed input-a -> outputs'); "
        "print('FAILED input-b: synthetic failure', file=sys.stderr); "
        "raise SystemExit(3)"
    )

    def runner() -> None:
        run_logged_command(
            [sys.executable, "-c", child],
            prefix="[TEST]",
            inputs=["input-a", "input-b"],
        )

    with pytest.raises(subprocess.CalledProcessError):
        run_step_with_audit("mixed_step", reports_dir, runner)

    text = _only_log(reports_dir).read_text()
    assert "[AUDIT] Status: FAILED" in text
    assert "[AUDIT] Inputs: 1 succeeded, 1 failed, 0 not run" in text
    assert "SUCCESS  input-a" in text
    assert "FAILED   input-b :: synthetic failure" in text


def test_step_audit_marks_all_inputs_success_for_zero_exit(tmp_path: Path) -> None:
    reports_dir = tmp_path / "reports"

    run_step_with_audit(
        "successful_step",
        reports_dir,
        lambda: run_logged_command(
            [sys.executable, "-c", "print('all good')"],
            prefix="[TEST]",
            inputs=["input-a", "input-b"],
        ),
    )

    text = _only_log(reports_dir).read_text()
    assert "[AUDIT] Status: SUCCESS" in text
    assert "[AUDIT] Inputs: 2 succeeded, 0 failed, 0 not run" in text


def test_step_audit_detects_error_lines_even_with_zero_exit(tmp_path: Path) -> None:
    reports_dir = tmp_path / "reports"

    run_step_with_audit(
        "semantic_error_step",
        reports_dir,
        lambda: print("ERROR: result validation failed"),
    )

    text = _only_log(reports_dir).read_text()
    assert "[AUDIT] Status: FAILED" in text
    assert "Error-bearing log lines: 1" in text
    assert "ERROR: result validation failed" in text


def test_error_line_is_attributed_to_matching_input(tmp_path: Path) -> None:
    reports_dir = tmp_path / "reports"

    run_step_with_audit(
        "input_error_step",
        reports_dir,
        lambda: run_logged_command(
            [
                sys.executable,
                "-c",
                "print('Band tclean failed for input-a: synthetic failure')",
            ],
            prefix="[TEST]",
            inputs=["input-a", "input-b"],
        ),
    )

    text = _only_log(reports_dir).read_text()
    assert "[AUDIT] Inputs: 1 succeeded, 1 failed, 0 not run" in text
    assert "FAILED   input-a :: Band tclean failed for input-a" in text
    assert "SUCCESS  input-b" in text


def test_independent_inputs_continue_after_one_fails(monkeypatch, tmp_path: Path) -> None:
    reports_dir = tmp_path / "reports"
    cfg = Config(
        project_name="continue-test",
        paths=PathsCfg(
            interim_dir=str(tmp_path / "interim"),
            reports_dir=str(reports_dir),
        ),
        reference=Target(name="reference", ms_paths=["reference.ms"]),
        tests=[Target(name="test", ms_paths=["test.ms"])],
    )
    attempted: list[str] = []

    def fake_run(command, *, inputs):
        input_label = list(inputs)[0]
        attempted.append(input_label)
        if input_label == "reference.ms":
            raise subprocess.CalledProcessError(7, list(command))
        record_success(input_label)

    monkeypatch.setattr(step2_vis_analysis, "_run_cmd", fake_run)
    with pytest.raises(StepInputsFailed):
        run_step_with_audit(
            "visibility_qa",
            reports_dir,
            lambda: step2_vis_analysis.run(cfg),
        )

    assert attempted == ["reference.ms", "test.ms"]
    text = _only_log(reports_dir).read_text()
    assert "[AUDIT] Inputs: 1 succeeded, 1 failed, 0 not run" in text
    assert "FAILED   reference.ms" in text
    assert "SUCCESS  test.ms" in text


def test_cli_all_audits_each_step_immediately(monkeypatch, tmp_path: Path) -> None:
    cfg = Config(
        project_name="audit-cli-test",
        paths=PathsCfg(reports_dir=str(tmp_path / "reports")),
    )
    events: list[str] = []

    monkeypatch.setattr(cli, "load_config", lambda path: cfg)
    monkeypatch.setattr(cli, "ALL_STEPS", ("first", "second"))
    monkeypatch.setattr(
        cli,
        "STEP_RUNNERS",
        {
            "first": ("first_step", lambda config: events.append("ran first")),
            "second": ("second_step", lambda config: events.append("ran second")),
        },
    )

    def audited(step_name, reports_dir, runner):
        events.append(f"audit start {step_name}")
        runner()
        events.append(f"audit end {step_name}")

    monkeypatch.setattr(cli, "run_step_with_audit", audited)
    cli.main(["--config", "unused.yaml", "all"])

    assert events == [
        "audit start first_step",
        "ran first",
        "audit end first_step",
        "audit start second_step",
        "ran second",
        "audit end second_step",
    ]
