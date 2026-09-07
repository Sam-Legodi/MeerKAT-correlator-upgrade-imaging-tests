# Prompt for the next `meerkat_corr_imaging.cli` planning cycle

Work in planning mode. Inspect the repository, current configurations, representative FITS headers, generated catalogues, existing verification reports and tests before proposing changes. Do not implement the plan in this task. Ask every question whose answer would materially affect the design, especially the authoritative primary-beam model and its version.

Repository:

`/Users/samuel/SARAO/MK+/Correlator/MeerKAT-correlator-upgrade-imaging-tests`

The current pipeline now generates consolidated image-domain reports with a `report` step. Preserve the implemented acceptance definitions, like-for-like comparison safeguards, rigid astrometry model, CMC1 RMS reference, `draft_` DOCX naming, neutral third-person cautions, JSON metrics sidecars and render-based DOCX QA.

Plan the next implementation cycle around these goals:

1. Add frequency-aware primary-beam correction for every currently non-PB-corrected continuum image product, particularly the independent low- and high-band MFImage slices. Determine whether correction should be applied to the delivered MFImage cuboids before slicing, to extracted planes after slicing, or both. Make this an explicit design decision based on the product structure and FITS metadata.
2. Identify and require an authoritative MeerKAT L-band primary-beam model and version. Prefer a project-approved model or library already present in the environment. Do not invent coefficients. If the repository and metadata do not establish the model unambiguously, ask which model/version must be used.
3. Specify correct image and uncertainty handling: evaluate the beam at each product's effective frequency on its WCS grid; calculate `I_pbcor = I / PB`; define a configurable PB cutoff; mask pixels outside the cutoff rather than allowing unstable amplification; carry the correction state, model, frequency, cutoff and provenance in FITS headers and audit logs. Address how image RMS and uncertainty maps should be interpreted after correction.
4. Make PB correction state explicit pipeline metadata rather than inferring it only from filenames. Enforce reference/test compatibility at runtime and fail clearly on PB/non-PB mismatches. Preserve non-PB products and diagnostics where they remain useful, but make PB-corrected catalogues the standard for position and photometric comparisons.
5. Run PyBDSF independently on each corrected reference and test image. Review whether the present source/component catalogue selection is correct for the intended statistic, whether PyBDSF options such as shapelet fitting should be configurable for speed, and how minimum catalogue and quality-match counts should be represented in the report.
6. Extend tests with a synthetic beam and artificial point sources. Test central-flux invariance, off-axis flux recovery, cutoff masking, WCS/frequency handling, FITS provenance, idempotency, like-for-like enforcement, deterministic filenames, backward compatibility and failure messages. Include an end-to-end fixture that proves low/high PB-corrected products reach PyBDSF, cross-matching and reporting.
7. Review the 4kL high-band failure, which currently has only two catalogue rows, zero quality-selected matches and a measured/CMC1 RMS ratio of about 893. Plan diagnostics that distinguish severely corrupted input, source-finder failure, wrong plane/frequency selection and configuration errors without silently lowering quality criteria.
8. Retain the 5-arcsec cross-match search gate, or justify a different gate that still leaves the 1-arcsec positional acceptance statistic unbiased by selection. Preserve the scale-fixed translation-plus-one-rotation model because the expected release-to-release offset is global and rigid.
9. Plan visibility-domain work as a separate, data-dependent stage. When visibilities become available, calculate SDP flagging fraction by time, channel and correlation/baseline for auto- and cross-correlations; scan-averaged mean and RMS by polarisation, scan and antenna/baseline; and amplitude oscillations. Use a 51-channel, third-order Savitzky--Golay model over an RFI-free channel mask derived from flags. Apply the existing Pass thresholds: flagging below 20% and oscillation below 1%. Do not fabricate these metrics from calibration-report images.
10. Plan time/scan-dependent source-position and relative-flux comparisons only when per-scan images or suitable image products become available. State the data contract needed to enable this stage.
11. Assess whether calibration-report PDF ingestion can extract reliable structured values. Keep visual PDF review contextual unless numeric values can be traced to stable report data or verified extraction rules.
12. Keep reports succinct and evidence-led. Every unassessed checklist item must state the missing input. Every generated DOCX basename must begin with `draft_`. Render and visually inspect every final page. Update the implementation and wording patch logs.

The planning response must include:

- the current pipeline/data-flow map and the proposed insertion points;
- a file-by-file change plan;
- configuration/schema changes and migration behaviour;
- PB mathematics, cutoff and provenance decisions;
- test and validation matrix;
- operational risks, performance implications and rollback strategy;
- acceptance criteria for the implementation;
- a list of unanswered questions, led by the authoritative MeerKAT primary-beam model/version and the preferred cuboid-versus-slice correction point.

