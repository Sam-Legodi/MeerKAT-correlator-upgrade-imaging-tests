from pathlib import Path
from unittest import TestCase

from meerkat_corr_imaging.output_paths import draft_docx_path


class DraftDocxPathTests(TestCase):
    def test_prefixes_only_the_basename(self):
        report = Path(
            "/pipeline/Sky-CrossMatches/crossmatched-positions/"
            "ref_x_other_astrometry.docx"
        )

        self.assertEqual(
            draft_docx_path(report),
            str(report.with_name("draft_ref_x_other_astrometry.docx")),
        )

    def test_does_not_duplicate_existing_prefix(self):
        report = Path("/pipeline/reports/draft_flux_summary.docx")

        self.assertEqual(draft_docx_path(report), str(report))
