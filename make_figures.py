#!/usr/bin/env python3
"""Draw the README figures from the cells written by export_cells.

Standard library only. Reads cells_single.txt, cells_forest.txt and
viridis256.txt (the 256-entry viridis LUT, CC0, values as in matplotlib)
and writes SVG files into figures/.

    ./export_cells && python3 make_figures.py
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "figures"
FONT = "Helvetica, Arial, sans-serif"
MONO = "Menlo, Consolas, monospace"
MAX_LEVEL = 3


# ---------------------------------------------------------------- data

def read_cells(name):
    cells = []
    for line in (HERE / name).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        t, path, level, refined, x0, y0, h, value = line.split()
        cells.append(dict(tree=int(t), path=path[1:], level=int(level),
                          refined=refined == "1", x0=float(x0), y0=float(y0),
                          h=float(h), value=float(value)))
    return cells  # depth-first (pre-order) per tree, as written


def breadth_first(cells):
    """Breadth-first order within one tree: by depth, then by path."""
    return sorted(cells, key=lambda c: (c["level"], c["path"]))


LUT = [tuple(float(v) for v in line.split())
       for line in (HERE / "viridis256.txt").read_text().split("\n") if line.strip()]


def colour(v, lo=0.0, hi=1.0):
    k = round((min(max(v, lo), hi) - lo) / (hi - lo) * (len(LUT) - 1))
    r, g, b = LUT[k]
    return "#%02x%02x%02x" % (round(255 * r), round(255 * g), round(255 * b))


# ---------------------------------------------------------------- svg helpers

def svg(width, height, body, label):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img" aria-label="{label}" '
            f'font-family="{FONT}">\n'
            f'<rect width="{width}" height="{height}" fill="#ffffff"/>\n'
            + "\n".join(body) + "\n</svg>\n")


def text(x, y, s, size=12, anchor="middle", weight="normal", family=FONT, fill="#000000"):
    return (f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" text-anchor="{anchor}" '
            f'font-weight="{weight}" font-family="{family}" fill="{fill}">{s}</text>')


def rect(x, y, w, h, fill="none", stroke="#000000", sw=1.0):
    return (f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def line(x1, y1, x2, y2, sw=1.0, extra=""):
    return (f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
            f'stroke="#000000" stroke-width="{sw}" {extra}/>')


def circle(cx, cy, r, fill="#000000", sw=1.0):
    return (f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r}" fill="{fill}" '
            f'stroke="#000000" stroke-width="{sw}"/>')


class Frame:
    """Maps domain coordinates (y up) to pixels (y down)."""

    def __init__(self, ox, oy, scale, ymax):
        self.ox, self.oy, self.s, self.ymax = ox, oy, scale, ymax

    def x(self, x):
        return self.ox + self.s * x

    def y(self, y):
        return self.oy + self.s * (self.ymax - y)

    def cell(self, c, **kw):
        return rect(self.x(c["x0"]), self.y(c["y0"] + c["h"]), self.s * c["h"], self.s * c["h"], **kw)


def colourbar(x, y, w, h, label):
    """Vertical colour bar for [0, 1], low at the bottom."""
    body, n = [], 128
    for k in range(n):
        v = (k + 0.5) / n
        body.append(f'<rect x="{x}" y="{y + h * (1 - (k + 1) / n):.3f}" width="{w}" '
                    f'height="{h / n + 0.4:.3f}" fill="{colour(v)}" stroke="none"/>')
    body.append(rect(x, y, w, h, sw=1))
    for v in (0.0, 0.25, 0.5, 0.75, 1.0):
        ty = y + h * (1 - v)
        body.append(line(x + w, ty, x + w + 5, ty))
        body.append(text(x + w + 9, ty + 4, f"{v:.2f}", 11, "start"))
    body.append(f'<text x="{x + w + 48}" y="{y + h / 2}" font-size="12" text-anchor="middle" '
                f'transform="rotate(-90 {x + w + 48} {y + h / 2})">{label}</text>')
    return body


# ---------------------------------------------------------------- figures

def fig_child_order():
    s, ox, oy = 90, 70, 30
    body = [rect(ox, oy, 2 * s, 2 * s, sw=2),
            line(ox + s, oy, ox + s, oy + 2 * s, 1.25),
            line(ox, oy + s, ox + 2 * s, oy + s, 1.25)]
    for ic in range(4):
        ix, iy = ic % 2, ic // 2
        cx, cy = ox + s * (ix + 0.5), oy + s * (1.5 - iy)
        body.append(text(cx, cy + 2, str(ic), 30, weight="bold"))
        body.append(text(cx, cy + 24, f"(ix, iy) = ({ix}, {iy})", 11))
    # axes
    ax, ay = ox, oy + 2 * s + 24
    body += [f'<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
             f'markerHeight="7" orient="auto"><path d="M0 0L10 5L0 10z" fill="#000"/></marker></defs>',
             line(ax, ay, ax + 60, ay, 1.25, 'marker-end="url(#arr)"'),
             text(ax + 72, ay + 4, "x", 13, "start"),
             line(ox - 24, oy + 2 * s, ox - 24, oy + 2 * s - 60, 1.25, 'marker-end="url(#arr)"'),
             text(ox - 24, oy + 2 * s - 68, "y", 13)]
    body.append(text(ox + 2 * s + 30, oy + 80, "child index", 13, "start", "bold"))
    body.append(text(ox + 2 * s + 30, oy + 100, "ic = ix + 2·iy", 13, "start", family=MONO))
    body.append(text(ox + 2 * s + 30, oy + 122, "bit 0: −x or +x half", 12, "start"))
    body.append(text(ox + 2 * s + 30, oy + 140, "bit 1: −y or +y half", 12, "start"))
    return svg(460, 270, body, "Child numbering inside a refined cell")


def fig_refinement_steps(cells):
    s, gap, ox, oy = 200, 36, 20, 20
    body = []
    for k in range(MAX_LEVEL + 1):
        f = Frame(ox + k * (s + gap), oy, s, 1.0)
        shown = [c for c in cells if c["level"] == k or (c["level"] < k and not c["refined"])]
        for c in shown:
            split = c["level"] == k and c["refined"]
            body.append(f.cell(c, fill="#d0d0d0" if split else "#ffffff", sw=1.0))
        if k < MAX_LEVEL:
            for c in cells:
                if c["level"] != k:
                    continue
                x0, y0, h = c["x0"], c["y0"], c["h"]
                for px, py in ((x0, y0), (x0 + h, y0), (x0, y0 + h), (x0 + h, y0 + h),
                               (x0 + h / 2, y0 + h / 2)):
                    body.append(circle(f.x(px), f.y(py), 2.2, sw=0))
        body.append(rect(f.x(0), f.y(1), s, s, sw=2))
        caption = f"level {k}" if k < MAX_LEVEL else "level 3 (final)"
        body.append(text(f.x(0.5), oy + s + 22, caption, 13, weight="bold"))
        if k < MAX_LEVEL:
            n_test = sum(1 for c in cells if c["level"] == k)
            n_split = sum(1 for c in cells if c["level"] == k and c["refined"])
            body.append(text(f.x(0.5), oy + s + 40, f"{n_test} tested, {n_split} split", 12))
        else:
            n_leaf = sum(1 for c in cells if not c["refined"])
            body.append(text(f.x(0.5), oy + s + 40, f"{n_leaf} leaves, no further tests", 12))
    lx = ox
    ly = oy + s + 66
    body += [circle(lx + 6, ly - 4, 2.2, sw=0),
             text(lx + 16, ly, "the five samples tested in each cell (four corners and the centre)", 12, "start"),
             rect(lx, ly + 12, 12, 12, fill="#d0d0d0"),
             text(lx + 18, ly + 22, "spread of the samples > tol, so the cell splits into four", 12, "start")]
    width = ox * 2 + 4 * s + 3 * gap
    return svg(width, ly + 40, body, "Refinement of the single tree, level by level")


def fig_tree(cells):
    order = breadth_first(cells)
    bf_index = {c["path"]: i for i, c in enumerate(order)}
    dfs_index = {c["path"]: i for i, c in enumerate(cells)}
    by_path = {c["path"]: c for c in cells}
    slot, ox, oy, dy = 23, 92, 44, 84
    # leaf slots in depth-first order; parents centred over their children
    xpos, nxt = {}, [0]

    def place(path):
        c = by_path[path]
        if not c["refined"]:
            xpos[path] = ox + slot * (nxt[0] + 0.5)
            nxt[0] += 1
        else:
            for ic in range(4):
                place(path + str(ic))
            xpos[path] = 0.5 * (xpos[path + "0"] + xpos[path + "3"])
    place("")
    width_nodes = slot * nxt[0]
    ypos = lambda c: oy + dy * c["level"]

    r = 9.5
    body = []
    for c in cells:   # edges first
        if c["refined"]:
            for ic in range(4):
                ch = by_path[c["path"] + str(ic)]
                body.append(line(xpos[c["path"]], ypos(c) + r, xpos[ch["path"]], ypos(ch) - r, 0.8))
    for c in cells:
        x, y = xpos[c["path"]], ypos(c)
        refined = c["refined"]
        body.append(circle(x, y, r, "#000000" if refined else "#ffffff", 1.2))
        body.append(text(x, y + 3.3, str(bf_index[c["path"]]), 9.5, weight="bold",
                         fill="#ffffff" if refined else "#000000"))
        # depth-first index: under a leaf (nothing hangs below it), or above a
        # refined node on the side away from its incoming edge
        if not refined:
            body.append(text(x, y + r + 12, str(dfs_index[c["path"]]), 8.5, fill="#707070"))
        else:
            left = (x < xpos[c["path"][:-1]]) if c["path"] else True
            dx = -(r + 1) if left else (r + 1)
            body.append(text(x + dx, y - r - 2, str(dfs_index[c["path"]]), 8.5,
                             "end" if left else "start", fill="#707070"))

    # per-depth annotations: left = depth, right = count and descriptor bits
    rx = ox + width_nodes + 24
    for d in range(MAX_LEVEL + 1):
        row = [c for c in order if c["level"] == d]
        y = oy + dy * d
        body.append(text(18, y + 4, f"depth {d}", 12, "start", "bold"))
        body.append(text(rx, y - 2, f"{len(row)} " + ("cell" if len(row) == 1 else "cells"), 12, "start"))
        if d < MAX_LEVEL:
            bits = "".join("1" if c["refined"] else "0" for c in row)
            groups = " ".join(bits[i:i + 4] for i in range(0, len(bits), 4))
            body.append(text(rx, y + 14, groups, 11.5, "start", family=MONO))
        else:
            body.append(text(rx, y + 14, "no bits (deepest level)", 11.5, "start"))
    ly = oy + dy * MAX_LEVEL + 52
    body += [circle(24, ly - 4, 8, "#000000", 1.2), text(38, ly, "refined (bit 1)", 12, "start"),
             circle(156, ly - 4, 8, "#ffffff", 1.2), text(170, ly, "leaf (bit 0)", 12, "start"),
             text(270, ly, "number in a node: breadth-first (storage) index", 12, "start"),
             text(580, ly, "grey number beside it: depth-first (build) index", 12, "start", fill="#707070")]
    width = rx + 250
    return svg(int(width), int(ly + 20), body, "The single quadtree with depth-first and breadth-first indices")


def fig_mesh(cells, nx, ny, scale, label_trees, coloured, aria):
    pad_l, pad_t = 40, 16
    f = Frame(pad_l, pad_t, scale, ny)
    body = []
    leaves = [c for c in cells if not c["refined"]]
    for c in leaves:
        fill = colour(c["value"]) if coloured else "#ffffff"
        body.append(f.cell(c, fill=fill, sw=0.75))
    for j in range(ny):     # root cells
        for i in range(nx):
            body.append(rect(f.x(i), f.y(j + 1), scale, scale, sw=2.25))
    if label_trees:
        for j in range(ny):
            for i in range(nx):
                t = i + nx * j
                tx, ty = f.x(i) + 6, f.y(j + 1) + 6
                body.append(rect(tx, ty, 92, 20, fill="#ffffff", sw=1))
                body.append(text(tx + 46, ty + 14, f"tree {t} ({i},{j})", 11.5))
    for i in range(nx + 1):
        body.append(text(f.x(i), f.y(0) + 18, f"{i}", 12))
    for j in range(ny + 1):
        body.append(text(f.x(0) - 10, f.y(j) + 4, f"{j}", 12, "end"))
    width = pad_l + nx * scale + 20
    height = pad_t + ny * scale + 30
    if coloured:
        body += colourbar(width + 10, pad_t, 18, ny * scale, "f(x, y) at the cell centre")
        width += 100
    return svg(int(width), int(height), body, aria)


def fig_root_index(nx, ny):
    s, ox, oy = 96, 40, 20
    f = Frame(ox, oy, s, ny)
    body = []
    for j in range(ny):
        for i in range(nx):
            t = i + nx * j
            body.append(rect(f.x(i), f.y(j + 1), s, s, sw=2))
            body.append(text(f.x(i + 0.5), f.y(j + 0.5) + 2, f"{t}", 26, weight="bold"))
            body.append(text(f.x(i + 0.5), f.y(j + 0.5) + 22, f"(i, j) = ({i}, {j})", 11))
    for i in range(nx + 1):
        body.append(circle(f.x(i), f.y(0), 3, sw=0))
        body.append(text(f.x(i), f.y(0) + 20, f"{i}", 12))
    for j in range(1, ny + 1):
        body.append(circle(f.x(0), f.y(j), 3, sw=0))
    for j in range(ny + 1):
        body.append(text(f.x(0) - 10, f.y(j) + 4, f"{j}", 12, "end"))
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            body.append(circle(f.x(i), f.y(j), 3, sw=0))
    tx = ox + nx * s + 28
    body += [text(tx, oy + 34, "tree index", 13, "start", "bold"),
             text(tx, oy + 56, "t = i + nx·j", 13, "start", family=MONO),
             text(tx, oy + 80, f"nx = {nx}, ny = {ny} root cells", 12, "start"),
             text(tx, oy + 100, f"Dimensions = {nx + 1} × {ny + 1} × 1 points", 12, "start"),
             text(tx, oy + 128, "dots: coarse grid points", 12, "start"),
             text(tx, oy + 146, "(XCoordinates, YCoordinates)", 12, "start")]
    return svg(tx + 230, oy + ny * s + 34, body, "Root cells of the forest and their tree indices")


def main():
    OUT.mkdir(exist_ok=True)
    single = read_cells("cells_single.txt")
    forest = read_cells("cells_forest.txt")
    figures = {
        "child_order.svg": fig_child_order(),
        "refinement_steps.svg": fig_refinement_steps(single),
        "single_tree_diagram.svg": fig_tree(single),
        "single_tree_field.svg": fig_mesh(single, 1, 1, 360, False, True,
                                          "Leaves of the single tree coloured by f"),
        "root_indexing.svg": fig_root_index(4, 2),
        "forest_mesh.svg": fig_mesh(forest, 4, 2, 150, False, False,
                                    "Forest of eight quadtrees, leaf cells and tree numbering"),
        "forest_field.svg": fig_mesh(forest, 4, 2, 150, False, True,
                                     "Leaves of the forest coloured by f"),
    }
    for name, content in figures.items():
        (OUT / name).write_text(content)
        print("wrote", OUT / name)


if __name__ == "__main__":
    main()
