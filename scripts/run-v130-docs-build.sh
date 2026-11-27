#!/usr/bin/env bash
# Build the *actual* public Sphinx HTML, LaTeX and PDF documentation.
# Non-destructive: write diagnostics and generated files to a newly
# allocated sandbox/outputs/v130-docs-build-<commit>-<id>/ directory.
# Sphinx/TeX may be installed locally or from your normal docs CI image.
set -Eeuo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

python tests/test_v130_docs_release_source.py
python -m sphinx --version >/dev/null || {
    echo "ERROR: install Sphinx and sphinx_rtd_theme in the build environment" >&2
    exit 2
}
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-docs-build-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
echo "SCRIMP v1.3.0 complete documentation build"
echo "Source commit: $(git rev-parse HEAD)"
echo "Build output: $out"
echo "PDF requested: ${SCRIMP_V130_DOC_PDF:-1}"
echo "Strict warnings: ${SCRIMP_V130_DOC_STRICT:-1}"

opts=(-n)
if [[ "${SCRIMP_V130_DOC_STRICT:-1}" == 1 ]]; then
    opts+=(-W --keep-going)
fi
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export PYTHONPATH="$root:${PYTHONPATH:-}"

python -m sphinx "${opts[@]}" -b html docs/source "$out/html" \
    2>&1 | tee "$out/html.log"
test -s "$out/html/index.html" || {
    echo "FAIL: no generated Sphinx HTML index" >&2
    exit 3
}
echo "PASS_HTML: $out/html/index.html"

python -m sphinx "${opts[@]}" -b latex docs/source "$out/latex" \
    2>&1 | tee "$out/latex.log"
test -s "$out/latex/scrimp.tex" || {
    echo "FAIL: no generated Sphinx LaTeX scrimp.tex" >&2
    exit 3
}
echo "PASS_LATEX: $out/latex/scrimp.tex"

pdf_built=0
if [[ "${SCRIMP_V130_DOC_PDF:-1}" == 1 ]]; then
    if ! command -v latexmk >/dev/null 2>&1; then
        echo "ERROR: latexmk and a usable LaTeX toolchain are required for PDF" >&2
        echo "HTML/LaTeX are retained in $out; no PDF PASS claimed." >&2
        exit 4
    fi
    (cd "$out/latex" && latexmk -pdf -interaction=nonstopmode \
        -halt-on-error -file-line-error scrimp.tex) \
        2>&1 | tee "$out/pdf.log"
    test -s "$out/latex/scrimp.pdf" || {
        echo "FAIL: no nonempty actual scrimp.pdf was produced" >&2
        exit 5
    }
    pdf_built=1
    echo "PASS_PDF: $out/latex/scrimp.pdf"
fi

python - "$out" "$pdf_built" <<'PY'
import json,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1])
pdf=bool(int(sys.argv[2]))
def count_warnings(p):
    if not p.exists(): return None
    return sum("WARNING:" in line or "WARNING " in line
               for line in p.read_text(errors="replace").splitlines())
report={
    "schema":"scrimp-v130-sphinx-full-html-latex-pdf-build-v1",
    "status":("PASS_HTML_LATEX_PDF_BUILT" if pdf else
              "PASS_HTML_LATEX_PDF_SKIPPED_BY_REQUEST"),
    "source":subprocess.check_output(
        ["git","rev-parse","HEAD"],text=True).strip(),
    "html":str(root/"html/index.html"),
    "latex":str(root/"latex/scrimp.tex"),
    "pdf":str(root/"latex/scrimp.pdf") if pdf else None,
    "html_warnings":count_warnings(root/"html.log"),
    "latex_warnings":count_warnings(root/"latex.log"),
    "pdf_log_available":(root/"pdf.log").is_file(),
    "strict_warnings_enforced":bool(int(
        __import__("os").environ.get("SCRIMP_V130_DOC_STRICT","1"))),
}
(root/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(report,indent=2),flush=True)
PY
echo "DOCUMENTATION EVIDENCE: $out/summary.json"
