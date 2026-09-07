# DOCX report wording patch log

Date: 2026-08-31

## Scope reviewed

- All 10 DOCX reports found under `/Users/samuel/SARAO/MK+/Correlator`.
- DOCX-generating code in `positions_analysis.py`, `flux_analysis.py`, and `vis_amp_analyze.py`.
- Future per-scan overlay wording in `positions_analysis.py`.

## Findings and patches

| Finding | Occurrences | Patch |
|---|---:|---|
| Duplicated `Pro tipPro tip` and informal second-person PB-correction advice | 6 astrometry DOCX reports and the positions generator | Replaced with a standalone, bold-labelled `Caution:` paragraph in neutral third-person language. The caution now calls for like-for-like PB-correction states and states the potential 2φ bias directly. |
| Informal or second-person astrometry guidance (`think ...`, `your diagnostics`, `NB`, `your reference beam`) | 6 astrometry DOCX reports and the positions generator | Recast as formal interpretation and reference-convention text without addressing the reader or writer. |
| Writer-directed θ guidance (`Pair this`, `What you should hope to see`, `your ellipses`, `you've got`) | 6 astrometry DOCX reports and the positions generator | Replaced with an objective expected-behaviour summary plus a separate third-person caution about apparently uniform θ distributions with large or field-angle-dependent ρ/Bmaj. |
| Other direct-reader or conversational wording (`Use this`, `you can read it`, `They’re handy`) | 6 astrometry DOCX reports and the positions generator | Rewritten as concise statements of purpose and interpretation. |
| First-person flux-method summary (`we analyse`) | 2 flux DOCX reports and the flux generator | Replaced with `this report analyses` wording. |
| Imperative or note-like wording in future per-scan overlay descriptions | Positions generator only | Recast as descriptive third-person report text. |
| Similar wording in visibility-amplitude reports | 2 visibility DOCX reports reviewed | No matches found; no changes required. |
| Page-edge clipping exposed by render QA | 2 total-flux figures in the 32kL flux DOCX report | Reduced those two embedded figures from 6.5 to 6.1 inches so their captions remain fully visible. No analytical content changed. |

## Files patched

- 6 `draft_*_astrometry.docx` reports in `crossmatched-positions`.
- 2 `draft_fluxcmp_*.docx` reports in `crossmatched-fluxes`.
- `src/meerkat_corr_imaging/positions_analysis.py`.
- `src/meerkat_corr_imaging/flux_analysis.py`.

## Verification

- Reopened all 8 edited DOCX files successfully after saving.
- Confirmed that body text, tables, headers, and footers in all 10 reviewed DOCX files contain no remaining first- or second-person report language or any of the identified informal phrases.
- Confirmed that report-generator string literals contain none of the identified report phrases or first-/second-person wording.
- Compiled the modified Python package successfully.
- Ran the DOCX output-path unit tests successfully (2 tests).
- Rendered all edited DOCX files and visually inspected every page for layout regressions.
