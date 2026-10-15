#!/usr/bin/env python3
"""Generate the Starlight theory pages from the LaTeX manuscript.

The manuscript in ``docs/manuscript/`` is the canonical description of the BEACH
models. This script turns each chapter into a Starlight page with the same
section, equation, figure, table, and citation numbers as the PDF. Numbers are
counted here and checked against the ``.aux`` file of the PDF build, so the
manuscript must be built first (``make -C docs/manuscript``, or omit
``--skip-pdf``).

TikZ figures are compiled one by one with the manuscript preamble and converted
to SVG with ``pdftocairo``. Requirements: XeLaTeX with the manuscript fonts,
BibTeX, ``pdftocairo`` (poppler-utils), and pandoc 3 (``pypandoc-binary``).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_ROOT = REPO_ROOT / "docs" / "manuscript"
MAIN_TEX = MANUSCRIPT_ROOT / "beach_models.tex"
PDF_BUILD = MANUSCRIPT_ROOT / "build"
SITE_ROOT = REPO_ROOT / "docs-site"
CONTENT_DIR = SITE_ROOT / "src" / "content" / "docs" / "theory"
PUBLIC_DIR = SITE_ROOT / "public" / "manuscript"
WORK_DIR = REPO_ROOT / "build" / "manuscript-site"
SITE_BASE = "/BEACH"
PAGE_BASE = f"{SITE_BASE}/theory"
PUBLIC_BASE = f"{SITE_BASE}/manuscript"
GITHUB_EDIT_ROOT = "https://github.com/Nkzono99/BEACH/edit/main"
PDF_NAME = "beach_models.pdf"

# Preamble commands that the standalone figure documents need.
FIGURE_PREAMBLE = re.compile(
    r"\\(usepackage(\[[^\]]*\])?\{(fontspec|xeCJK|amsmath,amssymb,bm|graphicx,xcolor|tikz)\}"
    r"|set(main|sans|mono|CJKmain|CJKsans|CJKmono)font"
    r"|usetikzlibrary|input\{figures/styles\}|newcommand)"
)
# A page cropped to the picture with a 4 pt margin, without the standalone class.
FIGURE_DOCUMENT = r"""\documentclass{{article}}
{preamble}
\newsavebox\figurebox
\begin{{document}}
\sbox\figurebox{{{tikz}}}
\pdfpagewidth=\dimexpr\wd\figurebox+8pt\relax
\pdfpageheight=\dimexpr\ht\figurebox+\dp\figurebox+8pt\relax
\hoffset=-1in \voffset=-1in
\shipout\vbox{{\kern4pt\hbox{{\kern4pt\box\figurebox}}}}
\end{{document}}
"""
# Preamble macros that pandoc must expand.
MACRO_LINE = re.compile(r"^\\newcommand.*$", re.MULTILINE)

CJK = "\u3000-\u30ff\u4e00-\u9fff\uff00-\uffef"
CJK_BREAK = re.compile(rf"(?<=[{CJK}])\n[ \t]*(?=[{CJK}])")

# pandoc reads tabular but not tabularx; the X-type columns become plain columns.
TABULARX = re.compile(r"\\begin\{tabularx\}\{[^}]*\}\{((?:[^{}]|\{[^{}]*\})*)\}")

ENV_PATTERN = re.compile(
    r"\\(?P<cmd>chapter|section|subsection)(?P<star>\*?)\{"
    r"|\\begin\{(?P<env>equation|align|figure|table)\}"
)


@dataclass
class Label:
    number: str
    slug: str
    anchor: str


@dataclass
class Figure:
    key: str
    number: str
    tikz: str
    caption: str
    short_caption: str


@dataclass
class Chapter:
    number: int
    source: Path
    slug: str
    title: str
    body: str = ""
    figures: list[Figure] = field(default_factory=list)
    table_captions: dict[str, tuple[str, str]] = field(default_factory=dict)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--skip-pdf",
        action="store_true",
        help="Use the existing docs/manuscript/build outputs instead of running make.",
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Do not compile TikZ figures (pages still reference the SVG files).",
    )
    return parser.parse_args()


# ---------------------------------------------------------------- LaTeX helpers


def strip_comments(text: str) -> str:
    return re.sub(r"(?<!\\)%.*", "", text)


def read_group(text: str, start: int, open_char: str = "{", close_char: str = "}") -> tuple[str, int]:
    """Return the contents of the balanced group opening at ``text[start]`` and the index after it."""
    if text[start] != open_char:
        raise ValueError(f"expected {open_char!r} at offset {start}: {text[start:start + 40]!r}")
    depth = 0
    index = start
    while index < len(text):
        char = text[index]
        if char == "\\":
            index += 2
            continue
        if char == open_char:
            depth += 1
        elif char == close_char:
            depth -= 1
            if depth == 0:
                return text[start + 1 : index], index + 1
        index += 1
    raise ValueError(f"unbalanced group at offset {start}")


def find_env_end(text: str, env: str, start: int) -> int:
    """Return the offset just after the ``\\end{env}`` matching the environment body starting at ``start``."""
    pattern = re.compile(rf"\\(begin|end)\{{{re.escape(env)}\}}")
    depth = 1
    for match in pattern.finditer(text, start):
        depth += 1 if match.group(1) == "begin" else -1
        if depth == 0:
            return match.end()
    raise ValueError(f"unterminated environment {env}")


def extract_command(text: str, name: str) -> tuple[str | None, str | None, str]:
    """Remove the first ``\\name[opt]{arg}`` from ``text``; return (opt, arg, rest)."""
    match = re.search(rf"\\{name}(?![A-Za-z])", text)
    if not match:
        return None, None, text
    index = match.end()
    while index < len(text) and text[index].isspace():
        index += 1
    option = None
    if index < len(text) and text[index] == "[":
        option, index = read_group(text, index, "[", "]")
    argument, end = read_group(text, index)
    return option, argument, text[: match.start()] + text[end:]


def anchor_id(label: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", label.lower()).strip("-")


def marker(kind: str, value: str) -> str:
    """Return a plain-text placeholder that survives pandoc unchanged."""
    return f"ZZ{kind}zz{value.replace('-', 'X')}ZZ"


def split_rows(body: str) -> list[str]:
    """Split an align body at top-level ``\\\\`` separators."""
    rows: list[str] = []
    depth = 0
    current = 0
    index = 0
    while index < len(body):
        if body.startswith("\\begin{", index):
            depth += 1
        elif body.startswith("\\end{", index):
            depth -= 1
        elif depth == 0 and body.startswith("\\\\", index):
            rows.append(body[current:index])
            index += 2
            current = index
            continue
        index += 1
    rows.append(body[current:])
    return [row for row in rows if row.strip()]


# ---------------------------------------------------------------- inputs


def build_pdf() -> None:
    subprocess.run(["make", "-C", str(MANUSCRIPT_ROOT)], check=True)


def load_aux() -> dict[str, str]:
    aux = (PDF_BUILD / "beach_models.aux").read_text(encoding="utf-8")
    return {
        match.group(1): match.group(2)
        for match in re.finditer(r"^\\newlabel\{([^}]+)\}\{\{([^}]*)\}", aux, re.MULTILINE)
    }


def load_bibliography() -> list[tuple[str, str]]:
    bbl = (PDF_BUILD / "beach_models.bbl").read_text(encoding="utf-8")
    entries = []
    for chunk in re.split(r"\\bibitem", bbl)[1:]:
        match = re.match(r"\s*(?:\[(?:[^\[\]]|\{[^{}]*\})*\])?\{([^}]+)\}", chunk)
        if not match:
            raise ValueError(f"cannot parse bibitem: {chunk[:80]!r}")
        text = chunk[match.end() :].split("\\end{thebibliography}")[0]
        text = text.replace("\\newblock", " ").replace("\\penalty0", " ")
        text = re.sub(r"\\doi\{([^}]*)\}", r"doi: \1", text)
        entries.append((match.group(1), " ".join(text.split())))
    return entries


def chapter_sources() -> list[Path]:
    main = MAIN_TEX.read_text(encoding="utf-8")
    return [MANUSCRIPT_ROOT / f"{name}.tex" for name in re.findall(r"\\input\{(sections/[^}]+)\}", main)]


def preamble_macros() -> str:
    main = MAIN_TEX.read_text(encoding="utf-8")
    preamble = main.split("\\begin{document}")[0]
    return "\n".join(MACRO_LINE.findall(preamble)) + "\n"


def figure_preamble() -> str:
    main = MAIN_TEX.read_text(encoding="utf-8")
    preamble = strip_comments(main.split("\\begin{document}")[0])
    commands: list[str] = []
    for line in preamble.splitlines():
        if line.startswith("\\"):
            commands.append(line)
        elif commands and line.strip():
            commands[-1] += "\n" + line
    kept = [command for command in commands if FIGURE_PREAMBLE.match(command)]
    return "\n".join(kept).replace("\\input{figures/styles}", f"\\input{{{MANUSCRIPT_ROOT / 'figures' / 'styles'}}}")


def front_matter_sections() -> list[tuple[str, str]]:
    """Return the unnumbered front-matter chapters (title, LaTeX body) before the table of contents."""
    main = strip_comments(MAIN_TEX.read_text(encoding="utf-8"))
    front = main.split("\\mainmatter")[0].split("\\maketitle")[1].split("\\tableofcontents")[0]
    sections = []
    for match in re.finditer(r"\\(?:chapter|section)\*\{([^}]*)\}(.*?)(?=\\(?:chapter|section)\*\{|\Z)", front, re.DOTALL):
        body = re.sub(r"\\addcontentsline\{[^}]*\}\{[^}]*\}\{[^}]*\}", "", match.group(2))
        sections.append((match.group(1), body.strip()))
    return sections


# ---------------------------------------------------------------- numbering


def chapter_slug(number: int, source: Path) -> str:
    stem = re.sub(r"^\d+_", "", source.stem).replace("_", "-")
    return f"ch{number:02d}-{stem}"


def number_chapter(chapter: Chapter, labels: dict[str, Label]) -> None:
    """Number sections, equations, figures, and tables; replace them with pandoc-friendly LaTeX."""
    text = strip_comments(chapter.source.read_text(encoding="utf-8"))
    text = re.sub(
        r"\\input\{(figures/[^}]+)\}",
        lambda m: strip_comments((MANUSCRIPT_ROOT / f"{m.group(1)}.tex").read_text(encoding="utf-8")),
        text,
    )
    text = TABULARX.sub(lambda m: "\\begin{tabular}{" + m.group(1).replace("Y", "l").replace("X", "l") + "}", text)
    text = text.replace("\\end{tabularx}", "\\end{tabular}")
    text = re.sub(r"\\(begin|end)\{center\}", "", text)
    counters = {"section": 0, "subsection": 0, "equation": 0, "figure": 0, "table": 0}
    output: list[str] = []
    position = 0
    c = chapter.number

    def take_label(fragment: str) -> tuple[list[str], str]:
        found = re.findall(r"\\label\{([^}]+)\}", fragment)
        return found, re.sub(r"\\label\{[^}]+\}", "", fragment)

    def register(keys: list[str], number: str, anchor: str | None = None) -> str:
        anchors = ""
        for key in keys:
            identifier = anchor_id(key) if anchor is None else anchor
            labels[key] = Label(number, chapter.slug, identifier)
            anchors += marker("ANCHOR", identifier)
        return anchors

    while True:
        match = ENV_PATTERN.search(text, position)
        if not match:
            output.append(text[position:])
            break
        output.append(text[position : match.start()])

        if match.group("cmd"):
            command = match.group("cmd")
            title, end = read_group(text, match.end() - 1)
            following = re.match(r"\s*\\label\{([^}]+)\}", text[end:])
            keys = [following.group(1)] if following else []
            if following:
                end += following.end()
            if command == "chapter":
                chapter.title = title
                register(keys, str(c), "")
            elif match.group("star"):
                output.append(f"\\{command}{{{title}}}")
            else:
                if command == "section":
                    counters["section"] += 1
                    counters["subsection"] = 0
                    number = f"{c}.{counters['section']}"
                else:
                    counters["subsection"] += 1
                    number = f"{c}.{counters['section']}.{counters['subsection']}"
                anchors = register(keys, number)
                output.append(f"\\{command}{{{number}\u2003{title}{anchors}}}")
            position = end
            continue

        env = match.group("env")
        end = find_env_end(text, env, match.end())
        body = text[match.end() : end - len(f"\\end{{{env}}}")]

        if env == "equation":
            counters["equation"] += 1
            number = f"{c}.{counters['equation']}"
            keys, body = take_label(body)
            anchors = register(keys, number)
            output.append(f"{anchors}\\begin{{equation*}}{body.strip()}\\tag{{{number}}}\\end{{equation*}}")
        elif env == "align":
            rows = []
            anchors = ""
            for row in split_rows(body):
                keys, row = take_label(row)
                if re.search(r"\\(nonumber|notag)\b", row):
                    rows.append(re.sub(r"\\(nonumber|notag)\b", "", row).strip())
                    if keys:
                        raise ValueError(f"label on unnumbered row: {keys}")
                    continue
                counters["equation"] += 1
                number = f"{c}.{counters['equation']}"
                anchors += register(keys, number)
                rows.append(f"{row.strip()}\\tag{{{number}}}")
            output.append(f"{anchors}\\begin{{align*}}" + "\\\\\n".join(rows) + "\\end{align*}")
        elif env == "figure":
            counters["figure"] += 1
            number = f"{c}.{counters['figure']}"
            body = re.sub(r"^\s*\[[^\]]*\]", "", body)
            tikz_start = body.index("\\begin{tikzpicture}")
            tikz_end = find_env_end(body, "tikzpicture", tikz_start + len("\\begin{tikzpicture}"))
            short, caption, _ = extract_command(body, "caption")
            keys, _ = take_label(body)
            if len(keys) != 1:
                raise ValueError(f"figure {number} in {chapter.source.name} needs exactly one label")
            key = keys[0]
            register(keys, number, anchor_id(key))
            chapter.figures.append(
                Figure(key, number, body[tikz_start:tikz_end], caption or "", short or caption or "")
            )
            output.append(f"\n\n{marker('FIG', anchor_id(key))}\n\n")
        else:  # table
            counters["table"] += 1
            number = f"{c}.{counters['table']}"
            body = re.sub(r"^\s*\[[^\]]*\]", "", body)
            _, caption, body = extract_command(body, "caption")
            keys, body = take_label(body)
            body = body.replace("\\centering", "")
            identifier = anchor_id(keys[0]) if keys else f"tab-{c}-{counters['table']}"
            register(keys, number, identifier)
            chapter.table_captions[identifier] = (number, caption or "")
            output.append(f"\n\n{marker('TABCAP', identifier)}\n\n{body.strip()}\n\n")
        position = end

    chapter.body = "".join(output)


def check_numbers(labels: dict[str, Label], aux: dict[str, str]) -> None:
    mismatched = [
        f"{key}: page {label.number}, PDF {aux.get(key)}"
        for key, label in labels.items()
        if aux.get(key) != label.number
    ]
    missing = sorted(key for key in aux if key not in labels and not key.startswith(("ch:front", "sec:front")))
    if mismatched or missing:
        raise SystemExit(
            "manuscript numbering differs from the PDF build:\n  "
            + "\n  ".join(mismatched + [f"{key}: not found in chapters" for key in missing])
        )


# ---------------------------------------------------------------- cross references


def label_href(label: Label, current_slug: str) -> str:
    fragment = f"#{label.anchor}" if label.anchor else ""
    if label.slug == current_slug and fragment:
        return fragment
    return f"{PAGE_BASE}/{label.slug}.html{fragment}"


def resolve_references(text: str, labels: dict[str, Label], citations: dict[str, int], current_slug: str) -> str:
    def ref(match: re.Match[str]) -> str:
        label = labels.get(match.group(2))
        if label is None:
            raise SystemExit(f"undefined reference: {match.group(2)}")
        link = f"\\href{{{label_href(label, current_slug)}}}{{{label.number}}}"
        return f"({link})" if match.group(1) == "eqref" else link

    def cite_link(number: int, key: str) -> str:
        return f"\\href{{{PAGE_BASE}/references.html#ref-{anchor_id(key)}}}{{{number}}}"

    def cite(match: re.Match[str]) -> str:
        keys = [key.strip() for key in match.group(1).split(",")]
        missing = [key for key in keys if key not in citations]
        if missing:
            raise SystemExit(f"undefined citation: {missing}")
        # natbib sort&compress: three or more consecutive numbers become a range.
        numbers = sorted({citations[key]: key for key in keys}.items())
        runs: list[list[tuple[int, str]]] = []
        for item in numbers:
            if runs and item[0] == runs[-1][-1][0] + 1:
                runs[-1].append(item)
            else:
                runs.append([item])
        parts = []
        for run in runs:
            if len(run) >= 3:
                parts.append(f"{cite_link(*run[0])}--{cite_link(*run[-1])}")
            else:
                parts.extend(cite_link(*item) for item in run)
        return "[" + ", ".join(parts) + "]"

    text = re.sub(r"\\(ref|eqref)\{([^}]+)\}", ref, text)
    return re.sub(r"\\cite[pt]?\{([^}]+)\}", cite, text)


# ---------------------------------------------------------------- pandoc


def pandoc_path() -> str:
    if os.environ.get("PANDOC"):
        return os.environ["PANDOC"]
    try:
        import pypandoc

        return pypandoc.get_pandoc_path()
    except (ImportError, OSError):
        return shutil.which("pandoc") or "pandoc"


def to_markdown(latex: str, macros: str) -> str:
    # A line break between two Japanese characters is not a word space.
    latex = CJK_BREAK.sub("", latex)
    result = subprocess.run(
        [
            pandoc_path(),
            "--from=latex",
            "--to=gfm-tex_math_gfm+tex_math_dollars",
            "--wrap=none",
            "--shift-heading-level-by=1",
        ],
        input=macros + latex,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def finish_markdown(text: str) -> str:
    # remark-math renders display math only when the delimiters sit on their own lines.
    text = re.sub(r"\s*\$\$(.+?)\$\$\s*", lambda m: f"\n\n$$\n{m.group(1).strip()}\n$$\n\n", text, flags=re.DOTALL)
    text = re.sub(r"ZZANCHORzz([A-Za-z0-9]+)ZZ", lambda m: f'<span id="{m.group(1).replace("X", "-")}"></span>', text)
    return re.sub(r"\n{3,}", "\n\n", text).strip() + "\n"


def caption_markdown(caption: str, labels: dict[str, Label], citations: dict[str, int], slug: str, macros: str) -> str:
    latex = resolve_references(caption, labels, citations, slug)
    return finish_markdown(to_markdown(latex, macros)).strip()


# ---------------------------------------------------------------- figures


def compile_figures(chapters: list[Chapter]) -> None:
    """Compile each TikZ figure to PDF and convert it to SVG.

    A figure is recompiled only when its standalone source changed, so a run on a
    machine without ``pdftocairo`` can be finished by a second run elsewhere.
    """
    preamble = figure_preamble()
    work = WORK_DIR / "figures"
    work.mkdir(parents=True, exist_ok=True)
    target = PUBLIC_DIR / "figures"
    target.mkdir(parents=True, exist_ok=True)

    pdfs = []
    for chapter in chapters:
        for figure in chapter.figures:
            name = anchor_id(figure.key)
            source = work / f"{name}.tex"
            pdf = work / f"{name}.pdf"
            document = FIGURE_DOCUMENT.format(preamble=preamble, tikz=figure.tikz)
            if not pdf.exists() or not source.exists() or source.read_text(encoding="utf-8") != document:
                source.write_text(document, encoding="utf-8")
                pdf.unlink(missing_ok=True)
                subprocess.run(
                    ["xelatex", "-interaction=nonstopmode", "-halt-on-error", source.name],
                    cwd=work,
                    check=True,
                    stdout=subprocess.DEVNULL,
                )
            pdfs.append((pdf, target / f"{name}.svg"))

    if shutil.which("pdftocairo") is None:
        raise SystemExit(
            f"pdftocairo (poppler-utils) is required to convert the figures; the PDFs are in {work}"
        )
    for pdf, svg in pdfs:
        subprocess.run(["pdftocairo", "-svg", str(pdf), str(svg)], check=True)


# ---------------------------------------------------------------- pages


def frontmatter(title: str, description: str, order: int, label: str, source: str | None) -> str:
    lines = [
        "---",
        f"title: {json.dumps(title, ensure_ascii=False)}",
        f"description: {json.dumps(description, ensure_ascii=False)}",
    ]
    lines.append(f"editUrl: {GITHUB_EDIT_ROOT}/{source}" if source else "editUrl: false")
    lines += [
        "sidebar:",
        f"  order: {order}",
        f"  label: {json.dumps(label, ensure_ascii=False)}",
        "---",
        "",
    ]
    return "\n".join(lines)


def chapter_notice() -> str:
    return (
        ":::note\n"
        f"この章は[技術原稿の PDF]({PUBLIC_BASE}/{PDF_NAME})から自動生成しています。"
        "節・式・図・表・文献の番号は PDF と同じです。"
        f"設定と実行の手順は扱いません。\n:::\n\n"
    )


def render_chapter(chapter: Chapter, labels: dict[str, Label], citations: dict[str, int], macros: str) -> str:
    latex = resolve_references(chapter.body, labels, citations, chapter.slug)
    text = finish_markdown(to_markdown(latex, macros))

    for figure in chapter.figures:
        identifier = anchor_id(figure.key)
        caption = caption_markdown(figure.caption, labels, citations, chapter.slug, macros)
        alt = re.sub(r"[$\\{}]", "", figure.short_caption).strip()
        block = (
            f'<figure class="theory-figure" id="{identifier}">\n'
            f'<img src="{PUBLIC_BASE}/figures/{identifier}.svg" alt={json.dumps(alt, ensure_ascii=False)} loading="lazy" />\n'
            "<figcaption>\n\n"
            f"**図 {figure.number}**　{caption}\n\n"
            "</figcaption>\n</figure>"
        )
        text = text.replace(marker("FIG", identifier), block)

    for identifier, (number, caption_tex) in chapter.table_captions.items():
        caption = caption_markdown(caption_tex, labels, citations, chapter.slug, macros)
        text = text.replace(
            marker("TABCAP", identifier),
            f'<span id="{identifier}"></span>**表 {number}**　{caption}',
        )

    leftover = re.findall(r"ZZ[A-Z]+zz[A-Za-z0-9]+ZZ", text)
    if leftover:
        raise SystemExit(f"unresolved placeholders in {chapter.slug}: {leftover[:5]}")

    title = f"第 {chapter.number} 章　{chapter.title}"
    source = chapter.source.relative_to(REPO_ROOT).as_posix()
    head = frontmatter(title, f"BEACH 技術原稿 第 {chapter.number} 章", 10 + chapter.number, f"{chapter.number}. {chapter.title}", source)
    return head + chapter_notice() + text


def render_references(entries: list[tuple[str, str]], macros: str) -> str:
    lines = []
    for number, (key, latex) in enumerate(entries, start=1):
        body = finish_markdown(to_markdown(latex, macros)).strip().replace("\n", " ")
        lines.append(f'<span id="ref-{anchor_id(key)}"></span>[{number}] {body}\n')
    head = frontmatter("参考文献", "BEACH 技術原稿の参考文献", 99, "参考文献", "docs/manuscript/references.bib")
    return head + "\n".join(lines)


def render_index(
    chapters: list[Chapter], labels: dict[str, Label], citations: dict[str, int], macros: str
) -> str:
    parts = [
        f"[PDF をダウンロード]({PUBLIC_BASE}/{PDF_NAME})（{PDF_NAME}）\n",
        "BEACH の物理モデルと数値手法を、研究背景、支配方程式、離散化、計算手順、適用限界の順に説明する技術原稿です。"
        "設定ファイルの書き方や実行手順は、このサイトの使い方のページを参照してください。\n",
    ]
    for title, body in front_matter_sections():
        if re.search(r"[A-Za-z]{4,}", title) and not re.search(r"[\u3040-\u30ff\u4e00-\u9fff]", title):
            continue  # The English abstract duplicates the Japanese summary.
        parts.append(f"## {title}\n\n" + finish_markdown(to_markdown(resolve_references(body, labels, citations, "index"), macros)))
    parts.append("## 章\n")
    for chapter in chapters:
        parts.append(f"{chapter.number}. [{chapter.title}]({PAGE_BASE}/{chapter.slug}.html)")
    parts.append(f"\n[参考文献]({PAGE_BASE}/references.html)\n")
    head = frontmatter(
        "技術原稿: 物理モデルと数値アルゴリズム",
        "BEACH の物理モデルと数値手法の技術原稿。PDF と章ごとのページ。",
        0,
        "目次と PDF",
        "docs/manuscript/beach_models.tex",
    )
    return head + "\n".join(parts) + "\n"


def write_pages(pages: dict[str, str]) -> None:
    if CONTENT_DIR.exists():
        shutil.rmtree(CONTENT_DIR)
    CONTENT_DIR.mkdir(parents=True)
    for name, content in pages.items():
        (CONTENT_DIR / f"{name}.md").write_text(content, encoding="utf-8")


def main() -> int:
    args = parse_args()
    if not args.skip_pdf:
        build_pdf()

    macros = preamble_macros()
    labels: dict[str, Label] = {}
    chapters: list[Chapter] = []
    for number, source in enumerate(chapter_sources(), start=1):
        chapter = Chapter(number, source, chapter_slug(number, source), "")
        number_chapter(chapter, labels)
        chapters.append(chapter)
    check_numbers(labels, load_aux())

    entries = load_bibliography()
    citations = {key: number for number, (key, _) in enumerate(entries, start=1)}

    pages = {"index": render_index(chapters, labels, citations, macros), "references": render_references(entries, macros)}
    for chapter in chapters:
        pages[chapter.slug] = render_chapter(chapter, labels, citations, macros)

    PUBLIC_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(PDF_BUILD / PDF_NAME, PUBLIC_DIR / PDF_NAME)
    if not args.skip_figures:
        compile_figures(chapters)
    write_pages(pages)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
