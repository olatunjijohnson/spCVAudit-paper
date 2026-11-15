#!/usr/bin/env python3
"""Assemble the paper and supplementary material into one bookmarked PDF.

This optional presentation step needs PyMuPDF; scientific reproduction does not.
Original pages and their internal/external links are retained unchanged.
"""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parents[1]
PARTS = (
    ("Main paper", "manuscript/main.pdf"),
    ("Supplementary material", "manuscript/supplement.pdf"),
    ("S1: All reference comparisons", "figures/S1-reference-all.pdf"),
    ("S2: All prediction-model comparisons", "figures/S2-learners-all.pdf"),
    ("S3: Covariance fits", "figures/S3-covariance-fits.pdf"),
    ("S4: All settings with fitted covariance", "figures/S4-nested-all.pdf"),
    ("S5: Jura targets and mean models", "figures/S5-jura-targets.pdf"),
    ("Ten-fold sensitivity check", "exploratory/fold-count-report.pdf"),
    ("Ten-fold sensitivity: all settings", "exploratory/fold-count-comparison.pdf"),
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output, date_label="13 November 2026"):
    output = output.expanduser().resolve()
    receipt = output.with_suffix(".json")
    if output.exists() or receipt.exists():
        raise FileExistsError("Use a new output name; existing reading packs are preserved")
    result = fitz.open()
    cover = result.new_page(width=595.276, height=841.89)
    outline = [[1, "Reading guide", 1]]
    records = []
    for title, relative in PARTS:
        path = ROOT / "revision" / relative
        with fitz.open(path) as document:
            offset = len(result)
            result.insert_pdf(document, links=False, annots=True)
            # LaTeX uses named destinations for citations. Resolve these to
            # explicit page destinations so they survive document concatenation.
            for number, page in enumerate(document):
                for link in page.get_links():
                    target = {"kind": link["kind"], "from": link["from"]}
                    if link["kind"] in (fitz.LINK_GOTO, fitz.LINK_NAMED) and link.get("page", -1) >= 0:
                        point = link.get("to", fitz.Point(0, 0))
                        if link["kind"] == fitz.LINK_NAMED:
                            # Resolved named destinations use PDF coordinates;
                            # insert_link expects page coordinates with top origin.
                            point = point * document[link["page"]].transformation_matrix
                        target.update(kind=fitz.LINK_GOTO, page=link["page"] + offset,
                                      to=point, zoom=link.get("zoom", 0))
                    elif link["kind"] == fitz.LINK_URI:
                        target["uri"] = link["uri"]
                    else:
                        raise ValueError("Unsupported source link in " + relative)
                    result[offset + number].insert_link(target)
            outline.append([1, title, offset + 1])
            outline.extend([level + 1, text, page + offset]
                           for level, text, page in document.get_toc() if page > 0)
            records.append({"title": title, "source": str(path.relative_to(ROOT)),
                            "source_sha256": sha(path), "first_pdf_page": offset + 1,
                            "pages": len(document)})

    # Adding pages invalidates earlier Page handles; reload the cover now.
    cover = result[0]
    cover.insert_text((54, 73), "Spatial cross-validation", fontsize=24, fontname="hebo")
    cover.insert_text((54, 105), "Paper and supplementary material", fontsize=19, fontname="hebo")
    cover.insert_text((54, 135), f"{date_label}  |  Research draft", fontsize=11)
    text = ("Nearest-neighbour distances and prediction error in spatial cross-validation. "
            "Olatunji Johnson. The main paper is followed by the supplementary material, "
            "Figures S1-S5 and the ten-fold sensitivity analysis. Use the linked entries "
            "below or the PDF bookmarks to navigate.")
    assert cover.insert_textbox(fitz.Rect(54, 165, 540, 230), text, fontsize=11, lineheight=1.4) >= 0
    cover.draw_line((54, 251), (541, 251), color=(0.65, 0.65, 0.65), width=0.5)
    cover.insert_text((54, 275), "DOCUMENT", fontsize=10, fontname="hebo")
    cover.insert_text((481, 275), "PDF PAGES", fontsize=9, fontname="hebo")
    for i, record in enumerate(records):
        y = 307 + i * 31
        first = record["first_pdf_page"]
        last = first + record["pages"] - 1
        cover.insert_text((54, y), record["title"], fontsize=11)
        cover.insert_text((490, y), f"{first}-{last}" if first != last else str(first), fontsize=11)
        cover.insert_link({"kind": fitz.LINK_GOTO, "from": fitz.Rect(50, y - 14, 545, y + 5),
                           "page": first - 1, "to": fitz.Point(0, 0)})
    note = ("Page ranges refer to this combined PDF. Each document retains its original "
            "page numbers. The paper is a research draft.\n\n"
            "Code, data sources and reproduction instructions:\n"
            "github.com/olatunjijohnson/spCVAudit-paper\n\n"
            "The spCVAudit package and tutorials:\n"
            "olatunjijohnson.github.io/spCVAudit")
    assert cover.insert_textbox(fitz.Rect(54, 620, 540, 765), note, fontsize=10, lineheight=1.4) >= 0
    result.set_toc(outline, collapse=1)
    result.set_metadata({"title": "Nearest-neighbour distances and prediction error in spatial cross-validation",
                         "author": "Olatunji Johnson",
                         "subject": "Research draft and supplementary material"})
    output.parent.mkdir(parents=True, exist_ok=True)
    result.save(output, garbage=4, deflate=True)
    result.close()

    # Verify copied page text, dimensions and hyperlink destinations after writing.
    verified_pages = 0
    verified_links = 0
    with fitz.open(output) as merged:
        for record in records:
            offset = record["first_pdf_page"] - 1
            with fitz.open(ROOT / record["source"]) as original:
                for number, source_page in enumerate(original):
                    copied = merged[offset + number]
                    assert source_page.get_text() == copied.get_text(), record["source"]
                    assert source_page.rect == copied.rect, record["source"]
                    source_links, copied_links = source_page.get_links(), copied.get_links()
                    assert len(source_links) == len(copied_links), record["source"]
                    for before, after in zip(source_links, copied_links):
                        expected_kind = fitz.LINK_GOTO if before["kind"] == fitz.LINK_NAMED else before["kind"]
                        assert expected_kind == after["kind"]
                        if expected_kind == fitz.LINK_GOTO:
                            assert after["page"] == before["page"] + offset
                            expected_point = before.get("to", fitz.Point(0, 0))
                            if before["kind"] == fitz.LINK_NAMED:
                                expected_point = expected_point * original[before["page"]].transformation_matrix
                            assert abs(after["to"] - expected_point) < 1e-3
                        elif before["kind"] == fitz.LINK_URI:
                            assert before["uri"] == after["uri"]
                        verified_links += 1
                    verified_pages += 1
        assert len(merged) == verified_pages + 1
        assert len(merged[0].get_links()) == len(records)
        assert all(1 <= item[2] <= len(merged) for item in merged.get_toc())
        pages, bookmarks = len(merged), len(merged.get_toc())
    receipt.write_text(json.dumps({"created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "output": output.name, "sha256": sha(output), "pages": pages,
        "bookmarks": bookmarks, "copied_pages_verified": verified_pages,
        "original_links_verified": verified_links, "sources": records,
        "builder_sha256": sha(Path(__file__)), "pymupdf_version": fitz.VersionBind,
        "scope": "Document assembly and navigation checks; not independent scientific review"}, indent=2) + "\n")
    print(f"Created {pages}-page reading pack with {bookmarks} bookmarks; verified {verified_pages} copied pages and {verified_links} original links.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--date", default="13 November 2026", help="Date printed on the reading guide")
    args = parser.parse_args()
    build(args.output, args.date)
