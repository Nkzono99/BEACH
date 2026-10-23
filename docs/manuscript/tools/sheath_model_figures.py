"""Write the TikZ sources of the numerical outer-sheath-model figures.

    python tools/sheath_model_figures.py figures/data/sheath_model figures

Reads the CSV files written by tools/sheath_model_data.py and emits plain TikZ (no pgfplots), so the
manuscript builds with the same XeLaTeX run as the hand-drawn figures. Only NumPy is required.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------------------------- helpers
def read_csv(path: Path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    out = {}
    for key in rows[0]:
        try:
            out[key] = np.array([float(r[key]) for r in rows])
        except ValueError:
            out[key] = np.array([r[key] for r in rows])
    return out


def num(value: float) -> str:
    text = f"{value:.3f}".rstrip("0").rstrip(".")
    return "0" if text in ("-0", "") else text


class Axes:
    """A rectangular plot frame mapping data to centimetres, with optional log10 axes."""

    def __init__(self, x0, y0, width, height, xlim, ylim, xlog=False, ylog=False):
        self.x0, self.y0, self.w, self.h = x0, y0, width, height
        self.xlim, self.ylim, self.xlog, self.ylog = xlim, ylim, xlog, ylog
        self.cmds: list[str] = []

    def _t(self, v, lim, log):
        v = np.asarray(v, dtype=float)
        a, b = lim
        if log:
            with np.errstate(divide="ignore", invalid="ignore"):
                v, a, b = np.log10(v), math.log10(a), math.log10(b)
        return (v - a) / (b - a)

    def X(self, x):
        return self.x0 + self.w * self._t(x, self.xlim, self.xlog)

    def Y(self, y):
        return self.y0 + self.h * self._t(y, self.ylim, self.ylog)

    def pt(self, x, y):
        return f"({num(float(self.X(x)))},{num(float(self.Y(y)))})"

    def frame(self, xticks, yticks, xlabel=None, ylabel=None, grid=True, ylabel_offset=1.05,
              xlabel_offset=0.75, xticklabels=True, yticklabels=True):
        c = self.cmds
        x1, y1 = self.x0 + self.w, self.y0 + self.h
        if grid:
            for v, _ in xticks:
                X = num(float(self.X(v)))
                c.append(f"\\draw[black!10,line width=0.4pt] ({X},{num(self.y0)})--({X},{num(y1)});")
            for v, _ in yticks:
                Y = num(float(self.Y(v)))
                c.append(f"\\draw[black!10,line width=0.4pt] ({num(self.x0)},{Y})--({num(x1)},{Y});")
        c.append(f"\\draw[line width=0.6pt] ({num(self.x0)},{num(self.y0)}) rectangle ({num(x1)},{num(y1)});")
        for v, label in xticks:
            X = num(float(self.X(v)))
            c.append(f"\\draw[line width=0.6pt] ({X},{num(self.y0)})--++(0,0.1);")
            if xticklabels and label is not None:
                c.append(f"\\node[below,font=\\scriptsize\\sffamily] at ({X},{num(self.y0)}) {{{label}}};")
        for v, label in yticks:
            Y = num(float(self.Y(v)))
            c.append(f"\\draw[line width=0.6pt] ({num(self.x0)},{Y})--++(0.1,0);")
            if yticklabels and label is not None:
                c.append(f"\\node[left,font=\\scriptsize\\sffamily] at ({num(self.x0)},{Y}) {{{label}}};")
        if xlabel:
            c.append(f"\\node[font=\\small\\sffamily] at ({num(self.x0 + self.w / 2)},"
                     f"{num(self.y0 - xlabel_offset)}) {{{xlabel}}};")
        if ylabel:
            c.append(f"\\node[rotate=90,font=\\small\\sffamily] at ({num(self.x0 - ylabel_offset)},"
                     f"{num(self.y0 + self.h / 2)}) {{{ylabel}}};")

    def clip_begin(self):
        self.cmds.append("\\begin{scope}")
        self.cmds.append(f"\\clip ({num(self.x0)},{num(self.y0)}) rectangle "
                         f"({num(self.x0 + self.w)},{num(self.y0 + self.h)});")

    def clip_end(self):
        self.cmds.append("\\end{scope}")

    def line(self, x, y, style, min_step=0.015):
        """Polyline through finite points; non-finite values split the line."""
        X, Y = self.X(x), self.Y(y)
        finite = np.isfinite(X) & np.isfinite(Y)
        segments, current = [], []
        for xi, yi, ok in zip(X, Y, finite):
            if ok:
                current.append((float(xi), float(yi)))
            elif current:
                segments.append(current)
                current = []
        if current:
            segments.append(current)
        for seg in segments:
            if len(seg) < 2:
                continue
            kept = [seg[0]]
            for p in seg[1:-1]:
                if math.hypot(p[0] - kept[-1][0], p[1] - kept[-1][1]) >= min_step:
                    kept.append(p)
            kept.append(seg[-1])
            coords = " ".join(f"({num(a)},{num(b)})" for a, b in kept)
            self.cmds.append(f"\\draw[{style}] plot coordinates {{{coords}}};")

    def hline(self, y, style):
        Y = num(float(self.Y(y)))
        self.cmds.append(f"\\draw[{style}] ({num(self.x0)},{Y})--({num(self.x0 + self.w)},{Y});")

    def vline(self, x, style):
        X = num(float(self.X(x)))
        self.cmds.append(f"\\draw[{style}] ({X},{num(self.y0)})--({X},{num(self.y0 + self.h)});")

    def text(self, x, y, label, opts="font=\\scriptsize\\sffamily"):
        self.cmds.append(f"\\node[{opts}] at {self.pt(x, y)} {{{label}}};")

    def raw(self, cmd):
        self.cmds.append(cmd)


def legend(x, y, entries, dy=0.42, length=0.6):
    """Legend rows: (style, label); x,y is the top-left corner in cm."""
    cmds = []
    for i, (style, label) in enumerate(entries):
        Y = y - i * dy
        cmds.append(f"\\draw[{style}] ({num(x)},{num(Y)})--({num(x + length)},{num(Y)});")
        cmds.append(f"\\node[right,font=\\scriptsize\\sffamily] at ({num(x + length)},{num(Y)}) {{{label}}};")
    return cmds


def figure(body: list[str], caption_short: str, caption: str, label: str, source: str) -> str:
    return "\n".join([
        f"% Generated by tools/sheath_model_figures.py from {source}; edit the generator, not this file.",
        "\\begin{figure}[htbp]",
        "\\centering",
        "\\begin{tikzpicture}[beach diagram]",
        *body,
        "\\end{tikzpicture}",
        f"\\caption[{caption_short}]{{{caption}}}",
        f"\\label{{{label}}}",
        "\\end{figure}",
        "",
    ])


COLORS = {"A": "BeachBlue", "B": "BeachOrange", "C": "BeachTeal"}
DASH = {"A": "solid", "B": "dashed", "C": "dash dot"}


def ticks(values, fmt=None):
    return [(v, fmt(v) if fmt else f"${num(v)}$".replace("$-", "$-")) for v in values]


def log_ticks(lo_exp, hi_exp):
    return [(10.0 ** k, f"$10^{{{k}}}$") for k in range(lo_exp, hi_exp + 1)]


# --------------------------------------------------------------------------------------------- figures
def fig_profiles(data: Path) -> str:
    meta = {m["branch"]: m for m in json.loads((data / "examples.json").read_text())}
    body = []
    top = Axes(0, 4.3, 10.5, 3.6, (0, 60), (-5, 5))
    bottom = Axes(0, 0, 10.5, 3.6, (0, 60), (-0.5, 1.75))
    top.frame(ticks(range(0, 61, 10)), ticks([-4, -2, 0, 2, 4]), None, "電位 $\\phi$ [V]", xticklabels=False)
    bottom.frame(ticks(range(0, 61, 10)), ticks([-0.5, 0, 0.5, 1, 1.5]), "壁からの距離 $z$ [m]",
                 "電場 $E_z$ [V/m]")
    for ax in (top, bottom):
        ax.clip_begin()
        ax.hline(0, "black!45,dashed,line width=0.5pt")
    for branch in "ABC":
        prof = read_csv(data / f"profile_{branch}.csv")
        style = f"{COLORS[branch]},{DASH[branch]},line width=1.1pt"
        top.line(prof["z_m"], prof["phi_v"], style)
        bottom.line(prof["z_m"], prof["e_v_m"], style)
    for ax in (top, bottom):
        ax.clip_end()
    a = meta["A"]
    top.raw(f"\\draw[BeachBlue,dotted,line width=0.6pt] {top.pt(a['z_min_m'], -5)}--{top.pt(a['z_min_m'], a['phim_v'])};")
    top.raw(f"\\node[beach dot,fill=BeachBlue] at {top.pt(a['z_min_m'], a['phim_v'])} {{}};")
    top.text(a["z_min_m"] + 0.8, a["phim_v"] + 0.9, "$z_m$, $\\phi_m$", "anchor=west,font=\\scriptsize\\sffamily,text=BeachBlue")
    body += top.cmds + bottom.cmds
    entries = [(f"{COLORS[b]},{DASH[b]},line width=1.1pt",
                f"Type {b}：$\\alpha={meta[b]['alpha_deg']:.0f}^\\circ$、$\\phi_0={meta[b]['phi0_v']:.2f}$ V")
               for b in "ABC"]
    body += legend(4.6, 7.55, entries, dy=0.45)
    caption = ("零電流条件で解いた Type A、B、C の電位と電場。条件は表\\ref{tab:sm_example_inputs}で、"
               "型ごとに太陽高度 $\\alpha$ を変えた。Type A は $z_m\\approx" + f"{a['z_min_m']:.1f}" +
               "$ m に負の極小 $\\phi_m=" + f"{a['phim_v']:.3f}" + "$ V を持ち、極小の手前で電場の向きが反転する。"
               "Type B は Debye 長 $\\lambda_{D,\\mathrm{pe}}\\approx1.4$ m 程度で減衰し、Type C は太陽風電子の Debye 長"
               " $\\lambda_{D,e}\\approx8.7$ m 程度でゆっくり 0 V に近づく。電子 drift は 0 とした。"
               "曲線は $|\\phi|$ が 2.2 mV（$10^{-3}T_{\\mathrm{pe}}$）に達した点で止め、その先に 0 V の線を補っていない。")
    return figure(body, "Type A、B、C の電位と電場の数値解", caption, "fig:sm_profiles", "profile_[ABC].csv")


def fig_populations(data: Path) -> str:
    prof = read_csv(data / "profile_A.csv")
    meta = {m["branch"]: m for m in json.loads((data / "examples.json").read_text())}["A"]
    z, phi = prof["z_m"], prof["phi_v"]
    U = -phi  # electron potential energy [eV], 0 upstream
    zm, um, u0 = meta["z_min_m"], -meta["phim_v"], -meta["phi0_v"]
    ax = Axes(0, 0, 11.5, 7.6, (0, 45), (-4.3, 1.6))
    ax.raw(f"\\fill[black!4] {ax.pt(0, -4.3)} rectangle {ax.pt(zm, 1.6)};")
    ax.frame(ticks(range(0, 46, 5)), ticks([-4, -3, -2, -1, 0, 1]), "壁からの距離 $z$ [m]",
             "電子のポテンシャルエネルギー $U=-e\\phi$ [eV]", grid=False)
    ax.clip_begin()
    ax.hline(0, "black!45,dashed,line width=0.5pt")
    ax.line(z, U, "BeachInk,line width=1.4pt")
    ax.clip_end()

    def turning(level, side):
        mask = z <= zm if side == "lower" else z >= zm
        zz, uu = z[mask], U[mask]
        if side == "lower":
            idx = np.where(uu >= level)[0][0]
        else:
            idx = np.where(uu <= level)[0][0]
        a, b = idx - 1, idx
        return float(zz[a] + (level - uu[a]) * (zz[b] - zz[a]) / (uu[b] - uu[a]))

    arrow = "-{Latex[length=2mm]},line width=1.0pt"
    # Solar-wind electrons enter from upstream (right).
    w_pass, w_refl = 1.15, 0.2
    ax.raw(f"\\draw[BeachBlue,{arrow}] {ax.pt(45, w_pass)}--{ax.pt(0.15, w_pass)};")
    ax.text(30, w_pass + 0.22, "太陽風電子：障壁を越えて壁で吸収", "font=\\scriptsize\\sffamily,text=BeachBlue")
    zt = turning(w_refl, "upper")
    ax.raw(f"\\draw[BeachBlue,dashed,{arrow}] {ax.pt(45, w_refl)}--{ax.pt(zt, w_refl)};")
    ax.raw(f"\\node[beach dot,fill=BeachBlue] at {ax.pt(zt, w_refl)} {{}};")
    ax.text(zt + 1.0, w_refl - 0.55, "太陽風電子：上側で反射", "anchor=west,font=\\scriptsize\\sffamily,text=BeachBlue")
    # Photoelectrons leave the wall.
    w_esc, w_cap = 0.62, -0.9
    ax.raw(f"\\draw[BeachOrange,{arrow}] {ax.pt(0.15, w_esc)}--{ax.pt(45, w_esc)};")
    ax.text(31, w_esc + 0.22, "光電子：障壁を越えて上流へ脱出", "font=\\scriptsize\\sffamily,text=BeachOrange")
    zc = turning(w_cap, "lower")
    ax.raw(f"\\draw[BeachRed,{arrow}] {ax.pt(0.15, w_cap)}--{ax.pt(zc, w_cap)};")
    ax.raw(f"\\node[beach dot,fill=BeachRed] at {ax.pt(zc, w_cap)} {{}};")
    ax.text(zc + 0.6, w_cap, "光電子：下側で折り返し、壁へ戻る", "anchor=west,font=\\scriptsize\\sffamily,text=BeachRed")
    # Reference levels.
    ax.raw(f"\\draw[BeachInk,dotted,line width=0.6pt] {ax.pt(0, um)}--{ax.pt(45, um)};")
    ax.raw(f"\\node[beach dot] at {ax.pt(zm, um)} {{}};")
    ax.text(zm + 1.0, 0.5 * (um + w_esc), f"障壁 $-e\\phi_m={um:.2f}$ eV", "anchor=west,font=\\scriptsize\\sffamily")
    ax.raw(f"\\node[beach dot] at {ax.pt(0, u0)} {{}};")
    ax.text(0.6, u0 + 0.05, f"壁 $-e\\phi_0={u0:.2f}$ eV", "anchor=west,font=\\scriptsize\\sffamily")
    ax.text(9.5, -3.3, "下側（壁〜極小）", "font=\\small\\sffamily")
    ax.text(30, -3.3, "上側（極小〜上流）", "font=\\small\\sffamily")
    caption = ("Type A の電子のポテンシャルエネルギー $U(z)=-e\\phi(z)$ と軌道の分類。黒の曲線は"
               "図\\ref{fig:sm_profiles}の Type A の数値解で、色の水平線は保存する全エネルギー $W$ の例である（$W$ の値は説明用）。"
               "電子は $U(z)=W$ の点で折り返す。太陽風電子は $W>-e\\phi_m$ なら壁に達し、"
               "そうでなければ上側で反射する。光電子は壁での法線エネルギーが $e(\\phi_0-\\phi_m)$ を超えれば脱出し、"
               "それ以下なら下側で折り返して壁に戻る。この二つの閾値が、各位置に存在する速度領域を決める。")
    return figure(ax.cmds, "Type A の電子のエネルギー図と軌道の分類", caption, "fig:sm_populations", "profile_A.csv")


def fig_densities(data: Path) -> str:
    prof = read_csv(data / "profile_A.csv")
    meta = {m["branch"]: m for m in json.loads((data / "examples.json").read_text())}["A"]
    z = prof["z_m"]
    ax = Axes(0, 0, 10.5, 6.0, (0, 60), (0.03, 100), ylog=True)
    ax.frame(ticks(range(0, 61, 10)), [(0.1, "$0.1$"), (1, "$1$"), (10, "$10$"), (100, "$100$")],
             "壁からの距離 $z$ [m]", "数密度 [cm$^{-3}$]")
    ax.clip_begin()
    ax.vline(meta["z_min_m"], "BeachInk,dotted,line width=0.6pt")
    series = [("n_i_m3", "BeachInk,line width=1.2pt", "イオン $n_i$"),
              ("n_e_pass_m3", "BeachBlue,line width=1.1pt", "太陽風電子（通過）"),
              ("n_e_refl_m3", "BeachCyan,dashed,line width=1.1pt", "太陽風電子（反射）"),
              ("n_pe_free_m3", "BeachOrange,line width=1.1pt", "光電子（脱出）"),
              ("n_pe_capt_m3", "BeachRed,dash dot,line width=1.1pt", "光電子（捕捉）")]
    for key, style, _ in series:
        values = prof[key] * 1e-6
        values = np.where(values > 0, values, np.nan)
        ax.line(z, values, style)
    ax.clip_end()
    ax.text(meta["z_min_m"], 140, "$z_m$", "font=\\scriptsize\\sffamily")
    body = ax.cmds + legend(5.3, 2.1, [(s, label) for _, s, label in series], dy=0.38)
    caption = ("Type A（$\\alpha=60^\\circ$）の粒子集団ごとの数密度。縦軸は対数で、点線は電位極小の位置 $z_m$ である。"
               "反射した太陽風電子は極小より上流側だけに、捕捉された光電子は極小より壁側だけに存在する。"
               "壁の近くでは捕捉された光電子が負電荷を担い、上流では通過・反射電子と脱出光電子の和がイオン密度と釣り合う。"
               "各密度は電位から速度積分で直接求めた値であり、粒子の標本ではない。")
    return figure(body, "Type A の粒子集団ごとの数密度", caption, "fig:sm_densities", "profile_A.csv")


def fig_first_integral(data: Path) -> str:
    d = read_csv(data / "first_integral_A.csv")
    meta = {m["branch"]: m for m in json.loads((data / "examples.json").read_text())}["A"]
    phim, phi0 = meta["phim_v"], meta["phi0_v"]
    lower, upper = d["side"] == "lower", d["side"] == "upper"
    w, h = 6.2, 3.3
    ax_rl = Axes(0, 4.6, w, h, (phim, phi0), (-50, 5))
    ax_el = Axes(8.0, 4.6, w, h, (phim, phi0), (0, 3))
    ax_ru = Axes(0, 0, w, h, (phim, 0), (-0.25, 0.05))
    ax_eu = Axes(8.0, 0, w, h, (phim, 0), (0, 2.5))
    xt_l = ticks([0, 1, 2, 3])
    xt_u = ticks([-0.3, -0.2, -0.1, 0])
    ax_rl.frame(xt_l, ticks([-40, -20, 0]), "$\\phi$ [V]", "$\\rho/e$ [cm$^{-3}$]", ylabel_offset=1.0,
                xlabel_offset=0.6)
    ax_el.frame(xt_l, ticks([0, 1, 2, 3]), "$\\phi$ [V]", "$E^2$ [V$^2$/m$^2$]", ylabel_offset=0.8,
                xlabel_offset=0.6)
    ax_ru.frame(xt_u, ticks([-0.2, -0.1, 0]), "$\\phi$ [V]", "$\\rho/e$ [cm$^{-3}$]", ylabel_offset=1.0,
                xlabel_offset=0.6)
    ax_eu.frame(xt_u, ticks([0, 1, 2]), "$\\phi$ [V]", "$E^2$ [$10^{-4}$ V$^2$/m$^2$]", ylabel_offset=0.8,
                xlabel_offset=0.6)
    for ax, mask, key, scale, color in ((ax_rl, lower, "rho_over_e_m3", 1e-6, "BeachRed"),
                                        (ax_el, lower, "e2_v2_m2", 1.0, "BeachRed"),
                                        (ax_ru, upper, "rho_over_e_m3", 1e-6, "BeachBlue"),
                                        (ax_eu, upper, "e2_v2_m2", 1e4, "BeachBlue")):
        ax.clip_begin()
        ax.hline(0, "black!45,dashed,line width=0.5pt")
        ax.line(d["phi_v"][mask], d[key][mask] * scale, f"{color},line width=1.2pt")
        ax.clip_end()
    # Shade the positive and negative areas of rho_upper: their integrals cancel (connection condition).
    phi_u, rho_u = d["phi_v"][upper], d["rho_over_e_m3"][upper] * 1e-6
    for sign, color in ((-1, "BeachBlue!20"), (1, "BeachOrange!35")):
        sel = sign * rho_u >= 0
        xs, ys = phi_u[sel], rho_u[sel]
        if len(xs) > 1:
            pts = " ".join(ax_ru.pt(float(x), float(y)) for x, y in zip(xs[::3], ys[::3]))
            ax_ru.cmds.insert(0, f"\\fill[{color}] {ax_ru.pt(float(xs[0]), 0)} -- plot coordinates {{{pts}}} -- "
                                 f"{ax_ru.pt(float(xs[-1]), 0)} -- cycle;")
    for ax, title in ((ax_rl, "下側（壁〜極小）の電荷密度"), (ax_el, "下側の $E^2(\\phi)$"),
                      (ax_ru, "上側（極小〜上流）の電荷密度"), (ax_eu, "上側の $E^2(\\phi)$")):
        ax.raw(f"\\node[above,font=\\small\\sffamily] at ({num(ax.x0 + ax.w / 2)},{num(ax.y0 + ax.h)}) {{{title}}};")
    ax_ru.text(-0.27, -0.045, "負", "font=\\scriptsize\\sffamily")
    ax_ru.text(-0.11, 0.025, "正", "font=\\scriptsize\\sffamily")
    ax_eu.text(-0.03, 0.35, "$E^2(0)=0$", "anchor=east,font=\\scriptsize\\sffamily")
    ax_el.text(phi0, 2.85, "$E_0^2$", "anchor=east,font=\\scriptsize\\sffamily")
    body = ax_rl.cmds + ax_el.cmds + ax_ru.cmds + ax_eu.cmds
    caption = ("Type A（$\\alpha=60^\\circ$）の電荷密度 $\\rho(\\phi)$ と第一積分 $E^2(\\phi)$。横軸は位置ではなく電位で、"
               "上段は壁から極小まで、下段は極小から上流までの区間である。二つの区間では粒子集団が異なるため、"
               "同じ電位でも $\\rho$ が異なる。$E^2$ は極小 $\\phi_m$ で 0 から始まり、上側では上流 $\\phi=0$ で再び 0 に戻る。"
               "これは上側の電荷密度の負の面積と正の面積が打ち消すこと（式\\eqref{eq:sm_connection}）に等しい。"
               "$E^2$ は数値求積で独立に計算した。")
    return figure(body, "Type A の電荷密度と第一積分", caption, "fig:sm_first_integral", "first_integral_A.csv")


def fig_vdf(data: Path) -> str:
    d = read_csv(data / "vdf_A.csv")
    meta = {m["branch"]: m for m in json.loads((data / "examples.json").read_text())}["A"]
    n_i = 8.7e6
    vte = math.sqrt(2 * 1.602176634e-19 * 12.0 / 9.1093837015e-31)
    scale = vte / n_i
    panels = [("wall", "壁 $z=0$"), ("lower", f"下側 $z={0.5 * meta['z_min_m']:.1f}$ m"), ("upper", "上側 $z=40$ m")]
    body = []
    w, h, gap = 4.1, 4.4, 0.5
    for k, (pos, title) in enumerate(panels):
        ax = Axes(k * (w + gap), 0, w, h, (-6, 6), (1e-3, 30), ylog=True)
        ax.frame(ticks([-4, -2, 0, 2, 4]), log_ticks(-3, 1), "$v_z$ [$10^6$ m/s]",
                 "$f(v_z)\\,v_{te}/n_i$" if k == 0 else None, yticklabels=(k == 0), ylabel_offset=1.05)
        ax.raw(f"\\node[above,font=\\small\\sffamily] at ({num(ax.x0 + w / 2)},{num(h)}) {{{title}}};")
        ax.clip_begin()
        ax.vline(0, "black!45,dashed,line width=0.5pt")
        for species, keys, styles in (
                ("swe", ("g1_s_m4", "g2_s_m4", "g3_s_m4"),
                 ("BeachBlue,line width=1.1pt", "BeachCyan,line width=1.1pt", "BeachCyan,line width=1.1pt")),
                ("phe", ("g1_s_m4", "g2_s_m4", "g3_s_m4"),
                 ("BeachOrange,line width=1.1pt", "BeachRed,line width=1.1pt", "BeachRed,line width=1.1pt"))):
            sel = (d["position"] == pos) & (d["species"] == species)
            v = d["vz_mps"][sel] * 1e-6
            for key, style in zip(keys, styles):
                g = d[key][sel] * scale
                ax.line(v, np.where(g > 0, g, np.nan), style)
        ax.clip_end()
        body += ax.cmds
    entries = [("BeachBlue,line width=1.1pt", "太陽風電子（通過、壁へ向かう）"),
               ("BeachCyan,line width=1.1pt", "太陽風電子（反射、往復）"),
               ("BeachOrange,line width=1.1pt", "光電子（脱出、上流へ向かう）"),
               ("BeachRed,line width=1.1pt", "光電子（捕捉、往復）")]
    body += legend(0.3, -1.25, entries[:2], dy=0.45) + legend(7.0, -1.25, entries[2:], dy=0.45)
    caption = ("Type A（$\\alpha=60^\\circ$）の電子の法線速度分布 $f(v_z)$。$v_z>0$ は上流向きで、縦軸は対数、"
               "$v_{te}=\\sqrt{2eT_e/m_e}$ と上流イオン密度 $n_i$ で規格化した。"
               "分布は上流の Maxwell 分布と壁の半 Maxwell 分布をエネルギー保存で写したもので、"
               "到達できない速度で 0 になる。壁では太陽風電子は $v_z<0$ だけ、光電子は $v_z>0$ だけを持つ。"
               "下側では捕捉光電子が $|v_z|$ の小さい領域を両方向に占め、上側では反射電子が同じ役割を担う。")
    return figure(body, "Type A の電子の速度分布", caption, "fig:sm_vdf", "vdf_A.csv")


def fig_alpha_sweep(data: Path) -> str:
    d = read_csv(data / "alpha_sweep.csv")
    top = Axes(0, 4.6, 9.5, 4.0, (0, 90), (-16, 6))
    bottom = Axes(0, 0, 9.5, 3.6, (0, 90), (0.1, 5), ylog=True)
    top.frame(ticks(range(0, 91, 15)), ticks([-15, -10, -5, 0, 5]), None, "電位 [V]", xticklabels=False)
    bottom.frame(ticks(range(0, 91, 15)), [(0.1, "$0.1$"), (0.3, "$0.3$"), (1, "$1$"), (3, "$3$")], "太陽高度 $\\alpha$ [deg]",
                 "数流束 $/(n_i v_{\\mathrm{sw}})$", ylabel_offset=1.0)
    top.clip_begin()
    top.hline(0, "black!45,dashed,line width=0.5pt")
    for branch in "ABC":
        sel = d["branch"] == branch
        alpha = d["alpha_deg"][sel]
        # Break the curve where consecutive solved elevations are not adjacent.
        phi0 = d["phi0_v"][sel].copy()
        breaks = np.where(np.diff(alpha) > 1.5)[0]
        x = np.insert(alpha, breaks + 1, np.nan)
        y0 = np.insert(phi0, breaks + 1, np.nan)
        top.line(x, y0, f"{COLORS[branch]},{DASH[branch]},line width=1.2pt")
        if branch == "A":
            ym = np.insert(d["phim_v"][sel], breaks + 1, np.nan)
            top.line(x, ym, "BeachBlue,dotted,line width=1.2pt")
    top.clip_end()
    flux = 8.7e6 * 468e3
    sel = d["branch"] == "A"
    bottom.clip_begin()
    for key, style in (("gamma_e", "BeachBlue,line width=1.2pt"), ("gamma_i", "BeachInk,line width=1.2pt"),
                       ("gamma_esc", "BeachOrange,line width=1.2pt"),
                       ("gamma_emit", "BeachOrange,dashed,line width=0.9pt")):
        bottom.line(d["alpha_deg"][sel], d[key][sel] / flux, style)
    bottom.clip_end()
    body = top.cmds + bottom.cmds
    body += legend(9.9, 8.3, [(f"{COLORS[b]},{DASH[b]},line width=1.2pt", f"Type {b} の $\\phi_0$") for b in "ABC"]
                   + [("BeachBlue,dotted,line width=1.2pt", "Type A の $\\phi_m$")], dy=0.46)
    body += legend(9.9, 3.3, [("BeachBlue,line width=1.2pt", "太陽風電子 $\\Gamma_e$"),
                               ("BeachInk,line width=1.2pt", "イオン $\\Gamma_i$"),
                               ("BeachOrange,line width=1.2pt", "脱出光電子 $\\Gamma_{\\mathrm{esc}}$"),
                               ("BeachOrange,dashed,line width=0.9pt", "放出光電子 $\\Gamma_{\\mathrm{emit}}$")],
                   dy=0.46)
    body.append("\\node[right,font=\\scriptsize\\sffamily] at (9.9,1.2) {下段は Type A の根};")
    caption = ("太陽高度 $\\alpha$ を $1^\\circ$ 刻みで変えたときの零電流根。上段は壁電位 $\\phi_0$ と Type A の極小 $\\phi_m$、"
               "下段は Type A の根での壁への数流束である。各型を個別に指定し、直前の高度の根を初期値に加えて解いた。"
               "曲線の切れ目は、有限の探索で物理解が見つからなかった高度である（不存在の証明ではない）。"
               "下段では $\\Gamma_e=\\Gamma_i+\\Gamma_{\\mathrm{esc}}$ が零電流条件を表し、放出の大部分は壁へ戻る。")
    return figure(body, "太陽高度に対する零電流根", caption, "fig:sm_alpha_sweep", "alpha_sweep.csv")


MAP_COLORS = {1: ("BeachBlue!70", "A"), 2: ("BeachOrange!75", "B"), 4: ("BeachTeal!70", "C"),
              3: ("BeachMagenta!65", "A+B"), 5: ("BeachCyan!70", "A+C"), 6: ("BeachRed!65", "B+C"),
              7: ("BeachInk!70", "A+B+C"), 0: ("black!8", "なし")}


def fig_type_map(data: Path) -> str:
    d = read_csv(data / "type_map.csv")
    alphas = np.unique(d["alpha_deg"])
    ratios = np.unique(d["source_ratio"])
    found = {(a, r): int(t) for a, r, t in zip(d["alpha_deg"], d["source_ratio"], d["types_found"])}
    ax = Axes(0, 0, 10.0, 6.0, (0, 90), (0.5, 16), ylog=True)
    da = alphas[1] - alphas[0]
    lr = np.log(ratios)
    dl = lr[1] - lr[0]
    present = set()
    cmds = []
    for j, r in enumerate(ratios):
        # Overlap neighbouring cells slightly so that anti-aliasing leaves no seams between rows.
        lo, hi = math.exp(lr[j] - dl / 2), math.exp(lr[j] + 0.62 * dl)
        lo, hi = max(lo, 0.5), min(hi, 16.0)
        i = 0
        while i < len(alphas):
            t = found[(alphas[i], r)]
            k = i
            while k + 1 < len(alphas) and found[(alphas[k + 1], r)] == t:
                k += 1
            a0, a1 = max(alphas[i] - da / 2, 0.0), min(alphas[k] + 0.62 * da, 90.0)
            present.add(t)
            cmds.append(f"\\fill[{MAP_COLORS[t][0]}] {ax.pt(a0, lo)} rectangle {ax.pt(a1, hi)};")
            i = k + 1
    ax.cmds += cmds
    ax.frame(ticks(range(0, 91, 15)), [(0.5, "$0.5$"), (1, "$1$"), (2, "$2$"), (4, "$4$"), (8, "$8$"), (16, "$16$")],
             "太陽高度 $\\alpha$ [deg]", "光電子源の強さ $n_{\\mathrm{pe,ref}}/n_i$", grid=False)
    r0 = 64e6 / 8.7e6
    for alpha, label in ((60, "A"), (20, "B"), (10, "C")):
        ax.raw(f"\\node[draw=white,fill=BeachInk,circle,inner sep=1.6pt,line width=0.6pt] at {ax.pt(alpha, r0)} {{}};")
        ax.text(alpha, r0 * 1.22, label, "font=\\scriptsize\\sffamily\\bfseries,text=white")
    ax.hline(r0, "white,dashed,line width=0.5pt")
    body = ax.cmds
    order = [1, 2, 4, 3, 5, 6, 7, 0]
    y = 5.8
    for t in order:
        if t not in present:
            continue
        color, label = MAP_COLORS[t]
        body.append(f"\\fill[{color}] (10.5,{num(y - 0.15)}) rectangle (10.9,{num(y + 0.15)});")
        body.append(f"\\node[right,font=\\scriptsize\\sffamily] at (10.95,{num(y)}) {{{label}}};")
        y -= 0.45
    caption = (f"太陽高度 $\\alpha$ と光電子源の強さ $n_{{\\mathrm{{pe,ref}}}}/n_i$ の平面で、零電流の物理解が見つかった型。"
               f"格子は {len(alphas)}$\\times${len(ratios)} 点で、各点で A、B、C を個別に指定し、既定の初期値だけから解いた。"
               "他の条件は表\\ref{tab:sm_example_inputs}と同じで、白い破線と点は図\\ref{fig:sm_profiles}の例の位置である。"
               "「A+B」などは同じ条件で複数の型が見つかったことを表し、色は解の安定性による選択を表さない。"
               "「なし」は有限の探索で物理解が見つからなかった点で、解の不存在を保証しない。")
    return figure(body, "零電流の物理解が見つかる型の分布", caption, "fig:sm_type_map", "type_map.csv")


def fig_upstream_band(data: Path) -> str:
    d = read_csv(data / "upstream_band.csv")
    meta = {m["drift"]: m for m in json.loads((data / "upstream_band.json").read_text())}
    ax = Axes(0, 0, 10.5, 5.6, (1e-4, 0.3), (-2.5, 6), xlog=True)
    ax.frame(log_ticks(-4, -1), ticks([-2, 0, 2, 4, 6]), "上流からの電位の深さ $-\\phi$ [V]",
             "$\\mathrm{sgn}(E^2)\\sqrt{|E^2|}$ [mV/m]")
    band = meta["normal"]["band_v"]
    ax.raw(f"\\fill[BeachOrange!12] {ax.pt(1e-4, -2.5)} rectangle {ax.pt(band, 6)};")
    ax.clip_begin()
    ax.hline(0, "black!45,dashed,line width=0.5pt")
    ax.hline(-2.0, "BeachRed,dashed,line width=0.7pt")
    for mode, style in (("zero", "BeachBlue,line width=1.2pt"), ("normal", "BeachOrange,line width=1.2pt")):
        sel = d["drift"] == mode
        e2 = d["e2_v2_m2"][sel]
        ax.line(d["depth_v"][sel], np.sign(e2) * np.sqrt(np.abs(e2)) * 1e3, style)
    ax.clip_end()
    deficit = meta["normal"]["deficit_v_m"] * 1e3
    ax.text(1.5e-4, -2.0 + 0.3, "$-E_{\\mathrm{floor}}=-2$ mV/m", "anchor=west,font=\\scriptsize\\sffamily,text=BeachRed")
    ax.text(band, 5.4, f"帯の縁 {band * 1e3:.0f} mV", "anchor=west,font=\\scriptsize\\sffamily")
    ax.text(3e-3, -deficit - 0.45, f"場の不足 {deficit:.2f} mV/m", "font=\\scriptsize\\sffamily,text=BeachOrange")
    z, n = meta["zero"], meta["normal"]
    body = ax.cmds + legend(0.3, 5.2, [
        ("BeachBlue,line width=1.2pt", f"drift 0（$\\phi_0={z['phi0_v']:.2f}$ V）"),
        ("BeachOrange,line width=1.2pt", f"drift あり（$\\phi_0={n['phi0_v']:.2f}$ V）")], dy=0.5)
    caption = ("Type A（$\\alpha=60^\\circ$）の上側の区間で、上流 $\\phi=0$ の近くの $E^2(\\phi)$。"
               "縦軸は符号を保った $\\sqrt{|E^2|}$、横軸は上流からの電位の深さ（対数）である。"
               "電子 drift 0 の根は $E^2\\propto\\phi^2$ で 0 に近づく。太陽風の法線速度で drift した電子では、"
               "$h\\log(1/h)$ 項（式\\eqref{eq:hlogh}）によって上流に接する帯（橙の網掛け）で $E^2<0$ になる。"
               "帯の中の最小値の大きさが場の不足であり、この例では下限 $E_{\\mathrm{floor}}=2$ mV/m を下回るので根を受理する。"
               "drift の有無で根そのもの（$\\phi_0$、$\\phi_m$）も変わる。")
    return figure(body, "上流の帯と場の不足の数値例", caption, "fig:sm_upstream_band", "upstream_band.csv")


FIGURES = {
    "sheath_model_profiles.tex": fig_profiles,
    "sheath_model_populations.tex": fig_populations,
    "sheath_model_densities.tex": fig_densities,
    "sheath_model_first_integral.tex": fig_first_integral,
    "sheath_model_vdf.tex": fig_vdf,
    "sheath_model_alpha_sweep.tex": fig_alpha_sweep,
    "sheath_model_type_map.tex": fig_type_map,
    "sheath_model_upstream_band.tex": fig_upstream_band,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--only", nargs="*", default=None)
    args = parser.parse_args()
    for name, make in FIGURES.items():
        if args.only and name not in args.only:
            continue
        (args.output / name).write_text(make(args.data))
        print("wrote", args.output / name)


if __name__ == "__main__":
    main()
