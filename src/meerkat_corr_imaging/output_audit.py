"""End-of-run output inventory, including subprocess products and reused files."""
from __future__ import annotations

import contextlib
import contextvars
from dataclasses import dataclass, field
from pathlib import Path

from .output_paths import draft_docx_path


CURRENT_RUN = contextvars.ContextVar("mci_output_run", default=None)


def absolute(path):
    return Path(path).expanduser().resolve()


@dataclass(frozen=True)
class Scope:
    path: Path
    reused: bool = True
    pattern: str | None = None


def output_scopes(cfg, command):
    """Watch output locations only; never traverse input MeasurementSets."""
    scopes = [Scope(absolute(cfg.paths.reports_dir), False)]
    scopes.append(Scope(Path.cwd(), False, "*.log"))
    paths = []
    targets = [cfg.reference, *cfg.tests]
    if command == "vis":
        paths.extend(Path(cfg.paths.interim_dir) / Path(ms).stem
                     for target in targets for ms in target.ms_paths)
    elif command == "cal":
        paths.extend(absolute(ms).parent / "images"
                     for target in targets for ms in target.ms_paths)
        pair = cfg.casa.paired_astrometry
        if pair.get("enabled") and pair.get("output_dir"):
            paths.append(pair["output_dir"])
        if cfg.extra.get("force_calibrate"):
            for target in targets:
                for ms in target.ms_paths:
                    path = absolute(ms)
                    scopes.append(Scope(path.parent, True, path.name.replace(".ms", "cal") + "_*"))
    elif command in {"src", "low_high_slice"}:
        from .steps.step4_low_high_slice import resolved_targets
        slices = resolved_targets(cfg)
        if command == "low_high_slice":
            paths.extend(p for target in slices for p in (target.low_image, target.high_image))
        else:
            from .pybdsf_srcfind import catalogue_paths, compute_base_name, expand_images
            images = [image for target in targets for image in target.images]
            if cfg.low_high_slice.add_to_source_finding:
                images.extend(str(p) for target in slices for p in (target.low_image, target.high_image))
            images.extend(cfg.extra.get("images_globs", []))
            for image in expand_images(images):
                base = compute_base_name(Path(image).stem, cfg.pybdsf.base_prefix)
                paths.append(catalogue_paths(image, base)[0].parent)
                # PyBDSF may also leave its work directory and log beside the image.
                for pattern in ("*.pybdsf", "*.log"):
                    scopes.append(Scope(absolute(image).parent, False, pattern))
    elif command == "xm":
        for entry in cfg.extra.get("xmatch_pairs", []):
            if len(entry) < 2:
                continue
            path = absolute(entry[2] if len(entry) >= 3 else
                            Path(cfg.paths.sky_xmatches_dir) / f"{Path(entry[0]).stem}_X_{Path(entry[1]).stem}.fits")
            if path.parent.name != "Sky-CrossMatches":
                path = path.parent / "Sky-CrossMatches" / path.name
            paths.extend((path, path.parent / "xmatch_pybdsf.log"))
        from .survey_xmatch import resolve_survey_jobs
        try:
            jobs = resolve_survey_jobs(cfg)
        except (ValueError, KeyError, TypeError):
            # The runner can still produce legacy matches before reporting an
            # invalid survey job; retain those output locations for its audit.
            jobs = []
        for job in jobs:
            path = absolute(job["output"])
            paths.extend((path, path.with_suffix(".diagnostics.json")))
    elif command == "pos":
        paths.extend(absolute(a["xmatch_table"]).parent / "crossmatched-positions"
                     for a in cfg.extra.get("positions", []) if a.get("xmatch_table"))
    elif command == "flux":
        analyses = cfg.extra.get("flux") or []
        if isinstance(analyses, dict):
            analyses = [analyses]
        for analysis in analyses:
            if not analysis.get("ref_low_xmatch"):
                continue
            directory = absolute(analysis["ref_low_xmatch"]).parent / "crossmatched-fluxes"
            paths.append(directory)
            if analysis.get("docx_name"):
                path = absolute(draft_docx_path(directory / analysis["docx_name"]))
                paths.extend((path, path.with_suffix(".pdf"), path.with_suffix(".json"),
                              path.with_suffix(".thermal_noise.csv")))
    elif command == "report":
        paths.append(Path(cfg.paths.reports_dir) / "imaging_verification")
        report = cfg.extra.get("verification_report") or {}
        if report.get("output_docx"):
            label = cfg.tests[0].name.removeprefix("CMC2_")
            path = absolute(draft_docx_path(report["output_docx"], label))
            paths.extend((path, path.with_suffix(".pdf"), path.with_name(path.stem + "_metrics.json")))
    elif command == "sensitivity":
        from .sensitivity_workflow import json_path
        path = json_path(cfg)
        if path:
            paths.append(path)
    return scopes + [Scope(absolute(p)) for p in paths]


def snapshot(scopes):
    """Stat signatures detect rewrites even when bytes and size are identical."""
    files, reused = {}, set()
    for scope in scopes:
        if scope.pattern:
            paths = (p for match in scope.path.glob(scope.pattern)
                     for p in (match.rglob("*") if match.is_dir() else [match]))
        else:
            paths = scope.path.rglob("*") if scope.path.is_dir() else [scope.path]
        for path in paths:
            if not path.is_file():
                continue
            stat = path.stat()
            files[path] = (stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size, stat.st_ino)
            if scope.reused:
                reused.add(path)
    return files, reused


@dataclass
class StepOutputs:
    name: str
    scopes: list[Scope]
    before: dict
    after: dict = field(default_factory=dict)
    reused: set = field(default_factory=set)
    audit: object = None
    exception: BaseException | None = None
    reports: set = field(default_factory=set)


class RunOutputAudit:
    def __init__(self, cfg):
        self.cfg = cfg
        self.steps = []
        self.active = None
        self.warnings = []

    def __enter__(self):
        self.token = CURRENT_RUN.set(self)
        return self

    def __exit__(self, *exc):
        CURRENT_RUN.reset(self.token)
        self.print_summary()

    def finish_reports(self):
        """Attribute late PDF exports to the step that registered each DOCX."""
        seen = set()
        for step in reversed(self.steps):
            paths = [report.with_suffix(".pdf") for report in step.reports - seen]
            seen.update(step.reports)
            after, reused = self._snapshot([Scope(p) for p in paths])
            for path in paths:
                step.after.pop(path, None)
            step.after.update(after)
            step.reused.update(reused)

    @contextlib.contextmanager
    def finalization(self):
        scopes = [Scope(absolute(self.cfg.paths.reports_dir) / "produced_reports.log")]
        before, _ = self._snapshot(scopes)
        step = StepOutputs("run_finalization", scopes, before)
        self.steps.append(step)
        try:
            yield
        except BaseException as exc:
            step.exception = exc
            raise
        finally:
            self.finish_reports()
            step.after, step.reused = self._snapshot(scopes)

    def _snapshot(self, scopes):
        try:
            return snapshot(scopes)
        except OSError as exc:
            self.warnings.append(f"Output inventory incomplete: {exc}")
            return {}, set()

    @contextlib.contextmanager
    def step(self, command, name, manifest):
        try:
            scopes = output_scopes(self.cfg, command)
        except Exception as exc:
            # Invalid configuration must still reach the runner's normal audit.
            self.warnings.append(f"Could not resolve all {name} outputs: {exc}")
            scopes = [Scope(absolute(self.cfg.paths.reports_dir), False)]
        before, _ = self._snapshot(scopes)
        step = StepOutputs(name, scopes, before)
        self.steps.append(step)
        self.active = step
        old_count = len(self._reports(manifest))
        try:
            yield
        except BaseException as exc:
            step.exception = exc
            raise
        finally:
            step.after, step.reused = self._snapshot(scopes)
            step.reports = set(self._reports(manifest)[old_count:])
            self.active = None

    @staticmethod
    def _reports(manifest):
        import json
        if not manifest.exists():
            return []
        return [absolute(json.loads(line)) for line in manifest.read_text().splitlines() if line.strip()]

    def record_step_audit(self, audit, exception):
        if self.active is not None:
            self.active.audit = audit
            self.active.exception = exception

    def _intermediate_parent(self, path):
        if path.suffix.lower() in {".log", ".docx", ".pdf"}:
            return None
        interim = absolute(self.cfg.paths.interim_dir)
        if path.is_relative_to(interim):
            return path.parent
        # CASA table internals and PyBDSF work directories are not deliverables.
        for parent in path.parents:
            if parent.suffix in {".ms", ".cal", ".image", ".model", ".residual", ".psf", ".pb", ".sumwt", ".mask", ".weight", ".pybdsf"}:
                return parent.parent
        return None

    def _output_lines(self, step):
        entries = set()
        for path, signature in step.after.items():
            status = ("NEW" if path not in step.before else
                      "OVERWRITTEN" if step.before[path] != signature else "UNCHANGED")
            if status == "UNCHANGED" and path not in step.reused:
                continue
            parent = self._intermediate_parent(path)
            entries.add((status, parent or path, bool(parent)))
        lines = ["[AUDIT]  Outputs:"]
        for status, path, intermediate in sorted(entries, key=lambda e: (str(e[1]), e[0])):
            suffix = " (intermediate folder)" if intermediate else ""
            lines.append(f"  {status:<11} {path}{suffix}")
        return lines if entries else lines + ["  (none)"]

    def print_summary(self):
        print("\n[AUDIT] Run output summary:")
        for step in self.steps:
            if step.audit is None:
                lines = [f"[AUDIT] Step: {step.name}", "[AUDIT] Inputs: (not applicable)"]
                lines.extend(self._output_lines(step))
                if step.exception:
                    lines.append(f"[AUDIT] Step exception: {type(step.exception).__name__}: {step.exception}")
            else:
                lines = step.audit._summary_lines(step.exception)
                index = next(i for i, line in enumerate(lines) if line.startswith("[AUDIT] Inputs:"))
                lines[index + 1:index + 1] = self._output_lines(step)
            print("\n".join(lines), flush=True)
        if not self.steps:
            print("[AUDIT]  Outputs: (none; no steps ran)", flush=True)
        for warning in self.warnings:
            print(f"[AUDIT] {warning}", flush=True)
