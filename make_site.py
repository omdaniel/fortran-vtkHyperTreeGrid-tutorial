#!/usr/bin/env python3
"""Build the GitHub Pages version of TUTORIAL.md into docs/.

pandoc turns the Markdown into HTML (no highlighting of its own); this script
then highlights the Fortran and shell listings with a small free-form Fortran
lexer, so the colours are exactly the ones in site/style.css. Needs pandoc and
the Python standard library only.

    python3 make_site.py        # or: make site
"""
import html
import re
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "docs"
SITE = HERE / "site"

TITLE = "Naive adaptive quadtrees in modern Fortran"
SUBTITLE = "From depth-first to vtkHyperTreeGrid"
DESCRIPTION = ("Build an adaptive quadtree depth first, grow it into a forest of "
               "quadtrees, and lay it out breadth first, the way vtkHyperTreeGrid stores it.")

# --------------------------------------------------------------------------
# A small free-form Fortran lexer
# --------------------------------------------------------------------------

KEYWORDS = set("""
program module submodule contains end subroutine function recursive pure elemental
impure result call return if then else elseif endif do enddo while select case
default associate allocate deallocate interface abstract import use implicit none
only intrinsic error stop exit cycle print write read open close where elsewhere
forall block concurrent go to goto continue nullify inquire rewind backspace
format entry include generic final extends sequence enum enumerator
""".split())
TYPES = set("integer real logical character complex double precision type class procedure".split())
ATTRIBUTES = set("""
parameter allocatable intent pointer target optional save value dimension public
private protected volatile asynchronous contiguous bind external non_overridable
deferred nopass pass
""".split())
INTRINSICS = set("""
abs achar acos adjustl adjustr aimag aint all allocated anint any asin associated
atan atan2 bit_size btest ceiling char cmplx conjg cos cosh count cshift dble
digits dim dot_product dprod eoshift epsilon exp exponent floor fraction huge iachar
iand ibclr ibits ibset ichar ieor index int ior ishft ishftc kind lbound len len_trim
log log10 logical matmul max maxloc maxval merge min minloc minval mod modulo
move_alloc nint not null pack popcnt present product radix random_number
random_seed range repeat reshape scan selected_int_kind selected_real_kind shape
sign sin sinh size spacing spread sqrt sum tan tanh tiny transfer transpose trim
ubound unpack verify storage_size findloc norm2 iso_fortran_env real64 real32
int32 int64 output_unit error_unit
""".split())
NAME_AFTER = {"subroutine", "function", "module", "program", "call", "submodule"}

FORTRAN_TOKEN = re.compile(r"""
    (?P<co>![^\n]*)
  | (?P<st>'(?:[^'\n]|'')*'|"(?:[^"\n]|"")*")
  | (?P<dot>\.(?:true|false|not|and|or|eqv|neqv|eq|ne|lt|le|gt|ge)\.(?:_\w+)?)
  | (?P<nu>(?<![\w.])(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][+-]?\d+)?(?:_\w+)?)
  | (?P<word>[A-Za-z_]\w*)
  | (?P<op>::|=>|\*\*|==|/=|<=|>=|[-+*/=<>%&])
  | (?P<ws>\s+)
  | (?P<other>.)
""", re.VERBOSE | re.IGNORECASE)


def span(cls, text):
    text = html.escape(text, quote=False)
    return f'<span class="{cls}">{text}</span>' if cls else text


def highlight_fortran(code):
    out = []
    sig = []          # recent significant tokens, (kind, lowercased text)
    for m in FORTRAN_TOKEN.finditer(code):
        kind, text = m.lastgroup, m.group()
        low = text.lower()
        cls = None
        if kind == "co":
            cls = "co"
        elif kind == "st":
            cls = "st"
        elif kind == "dot":
            cls = "cn" if low.startswith((".true", ".false")) else "kw"
        elif kind == "nu":
            cls = "nu"
        elif kind == "word":
            prev = sig[-1] if sig else ("", "")
            prev2 = sig[-2] if len(sig) > 1 else ("", "")
            nxt = code[m.end():].lstrip(" ")
            if prev == ("op", "%"):
                cls = None            # a component name, even if it spells a keyword
            elif prev == ("word", "intent") or (prev == ("other", "(") and prev2 == ("word", "intent")):
                cls = "at" if low in ("in", "out", "inout") else None
            elif low in TYPES and (low not in ("type", "class", "procedure") or nxt.startswith("(")):
                cls = "ty"
            elif low in ("type",) and not nxt.startswith("("):
                cls = "kw"            # a type definition: "type :: qnode", "end type"
            elif low in KEYWORDS:
                cls = "kw"
            elif low in ATTRIBUTES:
                cls = "at"
            elif prev[0] == "word" and prev[1] in NAME_AFTER:
                cls = "nm"
            elif prev == ("op", "::") and prev2 == ("word", "type"):
                cls = "nm"
            elif prev == ("word", "type") and prev2 == ("word", "end"):
                cls = "nm"
            elif low in INTRINSICS and nxt.startswith("("):
                cls = "bu"
            elif low in INTRINSICS and prev == ("word", "use"):
                cls = "bu"
        elif kind == "op":
            cls = "op" if text in ("::", "=>", "%") else None
        out.append(span(cls, text))
        if kind not in ("ws", "co"):
            sig = (sig + [(kind if kind != "other" else "other", low)])[-3:]
    return "".join(out)


SHELL_LINE = re.compile(r"^(\s*)(\S+)(.*)$")


def highlight_shell(code):
    lines = []
    for line in code.split("\n"):
        if line.lstrip().startswith("#") or not line.strip():
            lines.append(span("co" if line.strip() else None, line))
            continue
        body, _, comment = line.partition("  #")
        m = SHELL_LINE.match(body)
        indent, cmd, rest = m.groups()
        rest = re.sub(r"(?<!\S)(-{1,2}[\w=.+-]+)",
                      lambda f: "\0at\0" + f.group(1) + "\0/\0", rest)
        parts = re.split(r"\0(at|/)\0", rest)
        rendered, cls = [], None
        for p in parts:
            if p == "at":
                cls = "at"
            elif p == "/":
                cls = None
            else:
                rendered.append(span(cls, p))
        text = indent + span("nm" if cmd.startswith("./") else "kw", cmd) + "".join(rendered)
        if comment:
            text += span("co", "  #" + comment)
        lines.append(text)
    return "\n".join(lines)


LANGS = {
    "fortran": ("Fortran", highlight_fortran),
    "bash": ("Shell", highlight_shell),
    "text": ("Output", lambda c: html.escape(c, quote=False)),
}

# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------

PRE = re.compile(r'<pre class="(\w+)"><code>(.*?)</code></pre>', re.DOTALL)
PRE_PLAIN = re.compile(r'<pre><code>(.*?)</code></pre>', re.DOTALL)


def code_block(match):
    lang, body = match.group(1), html.unescape(match.group(2))
    label, fn = LANGS.get(lang, (lang.title(), lambda c: html.escape(c, quote=False)))
    return (f'<div class="codeblock lang-{lang}"><span class="lang">{label}</span>'
            f'<pre class="code"><code>{fn(body)}</code></pre></div>')


def main():
    md = (HERE / "TUTORIAL.md").read_text(encoding="utf-8")
    md = re.sub(r"\A# [^\n]*\n+", "", md)            # the title goes in the masthead
    cmd = ["pandoc", "-f", "markdown", "-t", "html5", "--no-highlight",
           "--toc", "--toc-depth=2", "--template", str(SITE / "template.html"),
           "-M", f"title={TITLE}", "-M", f"subtitle={SUBTITLE}",
           "-M", f"pagetitle={TITLE}: {SUBTITLE.lower()}",
           "-M", f"description={DESCRIPTION}"]
    page = subprocess.run(cmd, input=md, capture_output=True, text=True, check=True).stdout

    page = PRE.sub(code_block, page)
    page = PRE_PLAIN.sub(lambda m: code_block(_Fake("text", m.group(1))), page)
    page = re.sub(r"<table>", '<div class="table-wrap"><table>', page)
    page = re.sub(r"</table>", "</table></div>", page)
    page = re.sub(r'<h2 id="([^"]+)">(.*?)</h2>',
                  r'<h2 id="\1"><a class="anchor" href="#\1" aria-hidden="true">#</a>\2</h2>', page)
    page = page.replace('src="figures/', 'loading="lazy" src="figures/')

    OUT.mkdir(exist_ok=True)
    (OUT / "index.html").write_text(page, encoding="utf-8")
    shutil.copy(SITE / "style.css", OUT / "style.css")
    shutil.copytree(HERE / "figures", OUT / "figures", dirs_exist_ok=True)
    (OUT / ".nojekyll").write_text("")
    print(f"wrote {OUT / 'index.html'}")


class _Fake:
    """Lets an un-tagged code block go through code_block() as plain output."""
    def __init__(self, lang, body):
        self._g = (None, lang, body)

    def group(self, i):
        return self._g[i]


if __name__ == "__main__":
    main()
