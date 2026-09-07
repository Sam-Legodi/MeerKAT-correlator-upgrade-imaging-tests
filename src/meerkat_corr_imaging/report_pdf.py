"""Render pipeline DOCX contents using existing Python dependencies only.

This is a report renderer, not a Word layout engine: paragraphs, tables and
embedded raster figures retain document order with independent PDF pagination.
"""
from __future__ import annotations

from io import BytesIO
from pathlib import Path
import os
import tempfile

from docx import Document
from docx.oxml.ns import qn
from docx.table import Table
from docx.text.paragraph import Paragraph
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib import image as mpimage


class _Pages:
    width, height, margin = 595.0, 842.0, 42.0

    def __init__(self, pdf):
        self.pdf = pdf
        self.number = 0
        self.figure = None
        self.new_page()

    def new_page(self):
        if self.figure is not None:
            self.finish_page()
        self.figure = Figure(figsize=(self.width / 72, self.height / 72), dpi=72)
        self.canvas = FigureCanvasAgg(self.figure)
        self.renderer = self.canvas.get_renderer()
        self.y = self.height - self.margin
        self.number += 1

    def finish_page(self):
        self.figure.text(.5, .025, str(self.number), ha='center', fontsize=8)
        self.pdf.savefig(self.figure)
        self.figure.clear()

    def wrap(self, text, width, size=10, bold=False):
        font = FontProperties(family='DejaVu Sans', size=size,
                              weight='bold' if bold else 'normal')
        def fits(value):
            return self.renderer.get_text_width_height_descent(value, font, False)[0] <= width
        lines = []
        for raw in text.expandtabs(4).split('\n'):
            line = ''
            for word in raw.split():
                candidate = (line + ' ' + word).strip()
                if fits(candidate):
                    line = candidate
                    continue
                if line:
                    lines.append(line)
                    line = ''
                # Break long paths/identifiers too, so table cells cannot overflow.
                for char in word:
                    if line and not fits(line + char):
                        lines.append(line)
                        line = ''
                    line += char
            lines.append(line)
        return lines or ['']

    def text(self, value, x, y, size=10, bold=False):
        self.figure.text(x / self.width, y / self.height, value,
                         fontsize=size, fontfamily='DejaVu Sans',
                         fontweight='bold' if bold else 'normal', va='top',
                         parse_math=False)

    def paragraph(self, paragraph, keep_next_height=0):
        style = paragraph.style.name if paragraph.style is not None else ''
        bold = style.startswith('Heading') or style == 'Title'
        size = 18 if style == 'Title' else 13 if bold else 10
        line_height = size * 1.4
        lines = self.wrap(paragraph.text, self.width - 2 * self.margin, size, bold)
        if paragraph.text:
            if self.y - min(len(lines), 3) * line_height - keep_next_height < self.margin:
                self.new_page()
            for line in lines:
                if self.y - line_height < self.margin:
                    self.new_page()
                self.text(line, self.margin, self.y, size, bold)
                self.y -= line_height
            self.y -= 8
        # Pictures are embedded in paragraphs; preserve their order.
        for blip in paragraph._p.iter(qn('a:blip')):
            relation = blip.get(qn('r:embed'))
            if relation:
                self.picture(paragraph.part.related_parts[relation].blob)

    def picture(self, blob):
        image = mpimage.imread(BytesIO(blob))
        ih, iw = image.shape[:2]
        width = self.width - 2 * self.margin
        height = width * ih / iw
        limit = self.height - 2 * self.margin - 12
        if height > limit:
            width *= limit / height
            height = limit
        if self.y - height < self.margin:
            self.new_page()
        axes = self.figure.add_axes([
            (self.width - width) / 2 / self.width,
            (self.y - height) / self.height, width / self.width, height / self.height])
        # PDF/SVG backends embed the original raster with interpolation='none'.
        # Default interpolation resamples it to the 72-DPI layout canvas, losing
        # source plot detail. Keep page geometry in points and preserve pixels.
        axes.imshow(image, interpolation='none', resample=False)
        axes.axis('off')
        self.y -= height + 12

    def table(self, table):
        from matplotlib.patches import Rectangle
        if not table.rows:
            return
        count = len(table.columns)
        width = (self.width - 2 * self.margin) / count
        line_height = 12

        def row_lines(row, bold):
            return [self.wrap(cell.text, width - 10, 8, bold) for cell in row.cells]

        def draw_chunk(lines, start, count, bold):
            height = count * line_height + 10
            for index, cell_lines in enumerate(lines):
                x = self.margin + index * width
                self.figure.add_artist(Rectangle(
                    (x / self.width, (self.y - height) / self.height),
                    width / self.width, height / self.height,
                    transform=self.figure.transFigure, linewidth=.4,
                    edgecolor='#abb5be', facecolor='#e9eef3' if bold else 'white', zorder=0))
                for offset, line in enumerate(cell_lines[start:start + count]):
                    self.text(line, x + 5, self.y - 5 - offset * line_height, 8, bold)
            self.y -= height

        header = row_lines(table.rows[0], True)
        header_count = max(map(len, header))
        for index, row in enumerate(table.rows):
            lines = header if index == 0 else row_lines(row, False)
            total = max(map(len, lines))
            start = 0
            while start < total:
                capacity = int((self.y - self.margin - 10) // line_height)
                if capacity < 1:
                    self.new_page()
                    # Repeat a reasonably sized header above continued table rows.
                    if index and header_count * line_height + 10 < 150:
                        draw_chunk(header, 0, header_count, True)
                    capacity = int((self.y - self.margin - 10) // line_height)
                count = min(total - start, capacity)
                draw_chunk(lines, start, count, index == 0)
                start += count
        self.y -= 12


def render_docx_pdf(source: Path) -> Path:
    """Atomically write a PDF beside a DOCX without an external converter."""
    source = Path(source).resolve()
    document = Document(source)
    target = source.with_suffix('.pdf')
    fd, temporary = tempfile.mkstemp(prefix='.mci-pdf-', suffix='.pdf', dir=source.parent)
    os.close(fd)
    try:
        with PdfPages(temporary) as pdf:
            pages = _Pages(pdf)
            elements = list(document.element.body)
            for index, element in enumerate(elements):
                if element.tag == qn('w:p'):
                    paragraph = Paragraph(element, document)
                    reserve = 0
                    if paragraph.style.name.startswith('Heading') and index + 1 < len(elements):
                        reserve = 36
                        following = elements[index + 1]
                        for blip in following.iter(qn('a:blip')):
                            relation = blip.get(qn('r:embed'))
                            if relation:
                                image = mpimage.imread(BytesIO(document.part.related_parts[relation].blob))
                                reserve = min((pages.width - 2 * pages.margin) * image.shape[0] / image.shape[1],
                                              pages.height - 2 * pages.margin - 60) + 12
                                break
                    pages.paragraph(paragraph, reserve)
                elif element.tag == qn('w:tbl'):
                    pages.table(Table(element, document))
            pages.finish_page()
        os.replace(temporary, target)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return target
