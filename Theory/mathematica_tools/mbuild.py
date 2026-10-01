"""Convert a marked-up source file into a list of Mathematica cells (cells.wl) for buildnb.wl.

Block markers (each on its own line):
  ::Title:: ::Subtitle:: ::Section:: ::Subsection:: ::Subsubsection::   one-line headings
  ::Text:: ::Note:: ::Verify::   paragraphs (blank line = new cell); inline $tex$, **bold**, *italic*, `code`
  ::Item::          bullets, one per line starting with "- "
  ::ItemNumbered::  numbered items, one per line
  ::Math:: ::BoxMath::  one display formula per line (mini-TeX); trailing \\tag{n} puts (n) in the margin
  ::Code::          input cell (Wolfram Language, ASCII only); evaluated by buildnb.wl
  ::Table::         rows on lines, columns separated by " | "; first row is the header
Usage: python3 mbuild.py source.src cells.wl
"""
import re, sys

NAMED = {
    'alpha': 'α', 'beta': 'β', 'gamma': 'γ', 'delta': 'δ', 'epsilon': 'ε', 'varepsilon': 'ε', 'zeta': 'ζ', 'eta': 'η',
    'theta': 'θ', 'kappa': 'κ', 'lambda': 'λ', 'mu': 'μ', 'nu': 'ν', 'xi': 'ξ', 'pi': 'π', 'rho': 'ρ', 'sigma': 'σ',
    'tau': 'τ', 'phi': 'φ', 'chi': 'χ', 'psi': 'ψ', 'omega': 'ω', 'Gamma': 'Γ', 'Delta': 'Δ', 'Theta': 'Θ',
    'Lambda': 'Λ', 'Xi': 'Ξ', 'Pi': 'Π', 'Sigma': 'Σ', 'Phi': 'Φ', 'Psi': 'Ψ', 'Omega': 'Ω',
    'nabla': '∇', 'partial': '∂', 'times': '×', 'cdot': '·', 'otimes': '⊗', 'le': '≤', 'ge': '≥', 'leq': '≤', 'geq': '≥',
    'approx': '≈', 'to': '→', 'infty': '∞', 'pm': '±', 'sum': '∑', 'int': '∫', 'propto': '∝', 'll': '≪', 'gg': '≫',
    'neq': '≠', 'ne': '≠', 'parallel': '∥', 'perp': '⊥', 'equiv': '≡', 'Rightarrow': '⇒', 'in': '∈', 'langle': '〈',
    'rangle': '〉', 'ldots': '…', 'cdots': '⋯', 'dots': '…', 'ell': 'l', 'circ': '°', 'sim': '∼', 'simeq': '≃',
    'lesssim': '≲', 'gtrsim': '≳', 'leftarrow': '←', 'leftrightarrow': '↔', 'mapsto': '↦',
    'quad': '\u2003', 'qquad': '\u2003\u2003', ',': '\u2009', ';': '\u2005', '!': '', '{': '{', '}': '}', '|': '‖',
    ' ': '\u2005', 'oint': '∮', 'Big': '', 'big': '', 'Bigg': '', 'bigg': '', 'lvert': '|', 'rvert': '|',
    'top': '\x01Transpose\x02', 'T': '\x01Transpose\x02', 'one': '\x01DoubleStruckOne\x02',
}
FUNCS = {'ln', 'exp', 'log', 'sin', 'cos', 'tan', 'arccos', 'arcsin', 'arctan', 'det', 'tr', 'max', 'min', 'Re', 'Im', 'diag'}
ARG1 = {'mathcal', 'mathbb', 'mathbf', 'boldsymbol', 'bf', 'mathrm', 'operatorname', 'rm', 'text', 'sf', 'dot', 'ddot', 'hat', 'widehat',
        'bar', 'tilde', 'vec', 'sqrt'}


def wl_str(s):
    out = []
    for ch in s:
        if ch == '\\': out.append('\\\\')
        elif ch == '"': out.append('\\"')
        elif ch == '\n': out.append('\\n')
        elif ch in '\x01\x02' or ord(ch) < 128: out.append(ch)
        else: out.append('\\:%04x' % ord(ch))
    return '"' + ''.join(out).replace('\x01', '\\[').replace('\x02', ']') + '"'


def row(items):
    return items[0] if len(items) == 1 else 'RowBox[{' + ', '.join(items) + '}]'


class TeX:
    def __init__(self, s): self.s, self.i, self.word = s, 0, False

    def peek(self): return self.s[self.i] if self.i < len(self.s) else None

    def skip(self):
        while self.peek() == ' ': self.i += 1

    def group(self):
        self.skip()
        if self.peek() == '{':
            self.i += 1; items = self.seq(end='}'); self.i += 1
            return items
        return [self.atom()]

    def matrix(self):
        self.skip(); assert self.peek() == '{'
        depth, j = 0, self.i
        while True:
            if self.s[j] == '{': depth += 1
            elif self.s[j] == '}':
                depth -= 1
                if depth == 0: break
            j += 1
        body = self.s[self.i + 1:j]; self.i = j + 1

        def split(text, sep):
            out, d, cur, k = [], 0, '', 0
            while k < len(text):
                if text[k] == '{': d += 1
                elif text[k] == '}': d -= 1
                if d == 0 and text.startswith(sep, k):
                    out.append(cur); cur = ''; k += len(sep); continue
                cur += text[k]; k += 1
            out.append(cur); return out
        rows = [split(r, '&') for r in split(body, '\\\\')]
        grid = '{' + ', '.join('{' + ', '.join(tex_boxes(e) if e.strip() else '""' for e in r) + '}' for r in rows) + '}'
        return f'RowBox[{{"(", GridBox[{grid}, ColumnSpacings -> 1, RowSpacings -> 0.8], ")"}}]'

    def atom(self):
        c = self.peek()
        if c == '\\':
            self.i += 1
            m = re.match(r'[A-Za-z]+|[,;!{}| ]', self.s[self.i:])
            if not m: raise ValueError('bad escape in ' + self.s)
            name = m.group(0); self.i += len(name)
            if name in ('left', 'right', 'big', 'Big', 'bigl', 'bigr', 'Bigl', 'Bigr'): return '""'
            if name == 'pmatrix': return self.matrix()
            if name in ('frac', 'tfrac', 'dfrac'):
                a = row(self.group()); b = row(self.group()); return f'FractionBox[{a}, {b}]'
            if name in ARG1:
                name = {'mathbf': 'bf', 'boldsymbol': 'bf', 'mathrm': 'rm', 'operatorname': 'rm', 'widehat': 'hat'}.get(name, name)
                self.skip()
                if name in ('bf', 'rm', 'text', 'sf'):
                    prev = self.word; self.word = name in ('rm', 'text')
                    g = self.group() if self.peek() == '{' else self.seq(end='}')
                    self.word = prev
                else:
                    g = self.group()
                r = row(g)
                return {'bf': f'StyleBox[{r}, FontWeight -> "Bold"]', 'rm': f'StyleBox[{r}, FontSlant -> "Plain"]',
                        'text': f'StyleBox[{r}, FontSlant -> "Plain"]', 'sf': f'StyleBox[{r}, FontFamily -> "Helvetica", FontSlant -> "Plain"]',
                        'dot': f'OverscriptBox[{r}, "."]', 'ddot': f'OverscriptBox[{r}, ".."]', 'hat': f'OverscriptBox[{r}, "^"]',
                        'bar': f'OverscriptBox[{r}, "_"]', 'tilde': f'OverscriptBox[{r}, "~"]', 'vec': f'OverscriptBox[{r}, "\\[RightVector]"]',
                        'sqrt': f'SqrtBox[{r}]',
                        'mathbb': wl_str('\x01DoubleStruckCapital' + g[0].strip('"') + '\x02'),
                        'mathcal': wl_str('\x01ScriptCapital' + g[0].strip('"') + '\x02')}[name]
            if name in NAMED: return wl_str(NAMED[name])
            if name in FUNCS: return f'StyleBox[{wl_str(name)}, FontSlant -> "Plain"]'
            raise ValueError('unknown command \\' + name + ' in ' + self.s)
        if c == '{':
            return row(self.group())
        if c.isalnum() or c == '.':
            m = re.match(r'[0-9]+(\.[0-9]+)?|' + ('[A-Za-z]+' if self.word else '[A-Za-z]'), self.s[self.i:])
            tok = m.group(0) if m else c
            self.i += len(tok); return wl_str(tok)
        self.i += 1
        return wl_str(c)

    def seq(self, end=None):
        items = []
        while self.peek() is not None and self.peek() != end:
            c = self.peek()
            if c in '_^':
                self.i += 1
                arg = row(self.group())
                base = items.pop() if items else '""'
                if base in ('")"', '"]"') and items:          # attach scripts to a whole bracketed group
                    opener = '"("' if base == '")"' else '"["'
                    depth, k = 0, len(items) - 1
                    while k >= 0:
                        if items[k] == base: depth += 1
                        elif items[k] == opener:
                            if depth == 0: break
                            depth -= 1
                        k -= 1
                    if k >= 0:
                        group = items[k:] + [base]; del items[k:]
                        base = 'RowBox[{' + ', '.join(group) + '}]'
                if self.peek() in ('^', '_') and self.peek() != c:
                    other = self.peek(); self.i += 1; arg2 = row(self.group())
                    sub, sup = (arg, arg2) if c == '_' else (arg2, arg)
                    items.append(f'SubsuperscriptBox[{base}, {sub}, {sup}]')
                else:
                    items.append(f'SubscriptBox[{base}, {arg}]' if c == '_' else f'SuperscriptBox[{base}, {arg}]')
            elif c == ' ':
                self.i += 1; items.append('" "')
            else:
                items.append(self.atom())
        return items


def tex_boxes(s):
    t = TeX(s.strip()); items = t.seq()
    assert t.i == len(t.s), 'unparsed TeX: ' + s
    return row(items) if items else '""'


def inline_math(s, bold=False):
    b = tex_boxes(s)
    if bold: b = f'StyleBox[{b}, FontWeight -> "Bold"]'
    return f'Cell[BoxData[FormBox[{b}, TraditionalForm]], FormatType -> TraditionalForm]'


TOKEN = re.compile(r'(\*\*[^*]+\*\*|\$[^$]+\$|\*[^*\s$][^*$]*\*|`[^`]+`)')


def text_data(s):
    parts = []
    for piece in TOKEN.split(s):
        if not piece: continue
        if piece.startswith('**'):
            for sub in re.split(r'(\$[^$]+\$)', piece[2:-2]):
                if not sub: continue
                parts.append(inline_math(sub[1:-1], bold=True) if sub.startswith('$') else f'StyleBox[{wl_str(sub)}, FontWeight -> "Bold"]')
        elif piece.startswith('$'):
            parts.append(inline_math(piece[1:-1]))
        elif piece.startswith('`'):
            parts.append(f'StyleBox[{wl_str(piece[1:-1])}, FontFamily -> "Source Code Pro", FontColor -> RGBColor[0.55, 0.1, 0.1]]')
        elif piece.startswith('*'):
            parts.append(f'StyleBox[{wl_str(piece[1:-1])}, FontSlant -> "Italic"]')
        else:
            parts.append(wl_str(piece))
    return 'TextData[{' + ', '.join(parts) + '}]'


def entry_boxes(s):
    parts = []
    for piece in re.split(r'(\$[^$]+\$)', s.strip()):
        if piece: parts.append(f'FormBox[{tex_boxes(piece[1:-1])}, TraditionalForm]' if piece.startswith('$') else wl_str(piece))
    return row(parts) if parts else '""'


def check_balance(s):
    stack, i, instr = [], 0, False
    pairs = {']': '[', '}': '{', ')': '('}
    while i < len(s):
        c = s[i]
        if instr:
            if c == '\\': i += 2; continue
            if c == '"': instr = False
        elif c == '"': instr = True
        elif s.startswith('(*', i): i = s.index('*)', i) + 2; continue
        elif c in '[{(': stack.append(c)
        elif c in ']})':
            assert stack and stack[-1] == pairs[c], f'unbalanced near: {s[max(0, i - 80):i + 20]}'
            stack.pop()
        i += 1
    assert not stack and not instr, 'unclosed ' + str(stack)


EXTRA = {'Note': ', Background -> RGBColor[0.96, 0.97, 1.0], CellFrame -> {{3, 0}, {0, 0}}, CellFrameColor -> RGBColor[0.35, 0.55, 0.85]',
         'Verify': ', Background -> RGBColor[0.985, 0.965, 0.92], CellFrame -> {{4, 0}, {0, 0}}, CellFrameColor -> RGBColor[0.85, 0.55, 0.15]'}


def build(src):
    blocks = re.split(r'^::(\w+)::[ \t]*\n?', src, flags=re.M)
    assert blocks[0].strip() == '', 'content before first marker'
    cells = []
    for style, body in zip(blocks[1::2], blocks[2::2]):
        if style in ('Title', 'Subtitle', 'Section', 'Subsection', 'Subsubsection'):
            cells.append(f'Cell[{text_data(body.strip())}, "{style}"]')
        elif style in ('Text', 'Note', 'Verify'):
            for para in re.split(r'\n\s*\n', body.strip()):
                cells.append(f'Cell[{text_data(" ".join(l.strip() for l in para.splitlines()))}, "Text"{EXTRA.get(style, "")}]')
        elif style == 'Item':
            items = []
            for l in body.strip().splitlines():
                if l.startswith('- '): items.append(l[2:].strip())
                else: items[-1] += ' ' + l.strip()
            cells += [f'Cell[{text_data(it)}, "Item"]' for it in items]
        elif style == 'ItemNumbered':
            cells += [f'Cell[{text_data(l.strip())}, "ItemNumbered"]' for l in body.strip().splitlines() if l.strip()]
        elif style in ('Math', 'BoxMath'):
            lines = [l for l in body.strip().splitlines() if l.strip()]
            for idx, l in enumerate(lines):
                opts = ''
                tag = re.search(r'\\tag\{([^}]*)\}\s*$', l)
                if tag:
                    l = l[:tag.start()]
                    opts += f', CellFrameLabels -> {{{{None, Cell[{wl_str("(" + tag.group(1) + ")")}, "Text", FontSize -> 12]}}, {{None, None}}}}'
                if style == 'BoxMath':
                    top, bot = int(idx == 0), int(idx == len(lines) - 1)
                    opts += f', CellFrame -> {{{{1, 1}}, {{{bot}, {top}}}}}, CellFrameColor -> GrayLevel[0.4], Background -> RGBColor[0.97, 0.98, 1.0]'
                cells.append(f'Cell[BoxData[FormBox[{tex_boxes(l)}, TraditionalForm]], "DisplayFormula"{opts}]')
        elif style == 'Code':
            code = body.strip('\n')
            assert all(ord(c) < 128 for c in code), 'non-ASCII in code: ' + code[:80]
            check_balance(code)
            cells.append(f'CodeCell[{wl_str(code)}]')
        elif style == 'Table':
            rows = [r.split(' | ') for r in body.strip().splitlines()]
            grid = '{' + ', '.join('{' + ', '.join(entry_boxes(c) for c in r) + '}' for r in rows) + '}'
            gb = (f'GridBox[{grid}, GridBoxAlignment -> {{"Columns" -> {{{{Left}}}}}}, GridBoxDividers -> {{"Columns" -> {{{{True}}}}, "Rows" -> {{{{True}}}}}}, '
                  f'GridBoxBackground -> {{"Rows" -> {{RGBColor[0.9, 0.93, 0.97], None}}}}, GridBoxSpacings -> {{"Columns" -> {{{{1}}}}, "Rows" -> {{{{0.6}}}}}}]')
            cells.append(f'Cell[TextData[{{Cell[BoxData[{gb}], FontFamily -> "Source Sans Pro"]}}], "Text"]')
        else:
            raise ValueError('unknown block ' + style)
    return 'cells = {\n' + ',\n'.join(cells) + '\n};\n'


if __name__ == '__main__':
    out = build(open(sys.argv[1], encoding='utf-8').read())
    check_balance(out)
    assert all(ord(c) < 128 for c in out)
    open(sys.argv[2], 'w', encoding='ascii').write(out)
    print('wrote', sys.argv[2], len(out), 'bytes')
