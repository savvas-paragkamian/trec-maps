"""From the raw plate render to a finished page: keyed background, vignette, overlays, PNG and vector PDF.

Everything is drawn with matplotlib in pixel coordinates of the render, so the PNG is
pixel-exact and the PDF embeds the raster once with all lines and text as vectors.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from scipy import ndimage

from .config import hex_to_rgb

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
SANS = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]
SERIF = ["Palatino", "Georgia", "Times New Roman", "DejaVu Serif"]


# ----------------------------------------------------------------------------- raster layer

def key_background(rgba: np.ndarray, depth: np.ndarray, background: str) -> np.ndarray:
    """Replace ray misses (depth NaN, grown by 1 px against the engine's grey backdrop) with the page
    colour and lay a soft shadow under the plate. Returns RGB uint8."""
    bg = hex_to_rgb(background) * 255.0
    out = rgba[..., :3].astype(np.float32)
    miss = ndimage.binary_dilation(~np.isfinite(depth), iterations=1)
    h, w = miss.shape
    shadow = ndimage.gaussian_filter((~miss).astype(np.float32), sigma=max(3.0, w / 320))
    shadow = np.roll(shadow, max(2, h // 90), axis=0) * 0.45
    page = bg[None, None, :] * (1.0 - shadow[..., None])
    out[miss] = page[miss]
    return np.clip(np.rint(out), 0, 255).astype(np.uint8)


def vignette(rgb: np.ndarray, strength: float = 0.12) -> np.ndarray:
    if strength <= 0:
        return rgb
    h, w = rgb.shape[:2]
    x = ((np.arange(w, dtype=np.float32) - w / 2) / (w / 2)) ** 2
    y = ((np.arange(h, dtype=np.float32) - h / 2) / (h / 2)) ** 2
    r = np.sqrt(x[None, :] + y[:, None])
    mask = np.clip(1.0 - strength * np.clip(r - 0.55, 0, None) ** 1.5 * 2.2, 0, 1).astype(np.float32)
    return np.clip(rgb.astype(np.float32) * mask[..., None], 0, 255).astype(np.uint8)


# ----------------------------------------------------------------------------- page

class Page:
    """A matplotlib figure the size of the render; sizes are given in reference pixels of a 3840 px wide frame."""

    def __init__(self, rgb: np.ndarray, dpi: int = 200):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        matplotlib.rcParams["pdf.fonttype"] = 42
        self.plt = plt
        self.h, self.w = rgb.shape[:2]
        self.dpi = dpi
        self.scale = self.w / 3840.0
        self.fig = plt.figure(figsize=(self.w / dpi, self.h / dpi), dpi=dpi)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.imshow(rgb, extent=(-0.5, self.w - 0.5, self.h - 0.5, -0.5), interpolation="none", zorder=0)
        self.ax.set_xlim(-0.5, self.w - 0.5)
        self.ax.set_ylim(self.h - 0.5, -0.5)
        self.ax.axis("off")

    def pt(self, ref_px: float) -> float:
        """Reference pixels (at 3840 wide) -> points at this figure's dpi."""
        return ref_px * self.scale * 72.0 / self.dpi

    def px(self, ref_px: float) -> float:
        return ref_px * self.scale

    def text(self, x, y, s, *, size, color, ha="left", va="top", serif=False, weight="normal", **kw):
        return self.ax.text(x, y, s, fontsize=self.pt(size), color=color, ha=ha, va=va, zorder=10,
                            fontfamily=SERIF if serif else SANS, fontweight=weight, **kw)

    def polyline(self, px, py, vis, *, color, width, occluded_alpha=0.35):
        """Visible runs solid, hidden runs dotted and faded."""
        px, py, vis = np.asarray(px), np.asarray(py), np.asarray(vis, bool)
        if len(px) < 2:
            return
        change = np.flatnonzero(vis[1:] != vis[:-1]) + 1
        starts = np.concatenate([[0], change])
        ends = np.concatenate([change, [len(px)]])
        for a, b in zip(starts, ends):
            b = min(b + 1, len(px))   # share the joint with the next run
            if vis[a]:
                self.ax.plot(px[a:b], py[a:b], color=color, lw=self.pt(width), solid_capstyle="round", zorder=5)
            else:
                self.ax.plot(px[a:b], py[a:b], color=color, lw=self.pt(width * 0.8), ls=(0, (1.5, 2.5)),
                             alpha=occluded_alpha, zorder=5)

    def marker(self, x, y, *, size, face, edge="white", edge_width=1.2, shape="o", alpha=1.0, zorder=7):
        self.ax.plot([x], [y], marker=shape, ms=self.pt(size), mfc=face, mec=edge, mew=self.pt(edge_width),
                     ls="none", alpha=alpha, zorder=zorder)

    def segment(self, x0, y0, x1, y1, *, color, width, ls="-", alpha=1.0, zorder=6):
        self.ax.plot([x0, x1], [y0, y1], color=color, lw=self.pt(width), ls=ls, alpha=alpha, zorder=zorder,
                     solid_capstyle="round")

    def gradient_bar(self, x, y, w, h, strip_rgb: np.ndarray, *, edge):
        """A horizontal colour bar: `strip_rgb` is (1, n, 3) sRGB floats, left = first value."""
        self.ax.imshow(strip_rgb, extent=(x, x + w, y + h, y), aspect="auto", interpolation="bilinear", zorder=8)
        from matplotlib.patches import Rectangle
        self.ax.add_patch(Rectangle((x, y), w, h, facecolor="none", edgecolor=edge, lw=self.pt(1.5), zorder=9))

    def save(self, out_png: Path, out_pdf: Path | None, *, title: str = "") -> None:
        out_png = Path(out_png)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        self.fig.savefig(out_png, dpi=self.dpi, facecolor="white")
        if out_pdf:
            self.fig.savefig(out_pdf, dpi=self.dpi, metadata={"Title": title, "Creator": "trec-maps / forge3d"})
        self.plt.close(self.fig)


# ----------------------------------------------------------------------------- furniture

def draw_furniture(page: Page, *, title: str, subtitle: str, credits: list[str], legend: dict, style: dict) -> None:
    """Title block top-left, elevation and depth bars bottom-left, credits bottom-right.

    legend = {"land": (strip, lo, hi), "sea": (strip, lo, hi), "note": str}
    """
    W, H = page.w, page.h
    text, muted, accent = style["text"], style["muted"], style["accent"]
    margin = page.px(110)
    page.text(margin, margin, title, size=104, color=text, serif=True, weight="bold", va="top")
    page.segment(margin, margin + page.px(135), margin + page.px(150), margin + page.px(135), color=accent, width=4)
    page.text(margin, margin + page.px(160), subtitle, size=40, color=muted, va="top")

    bar_w, bar_h = page.px(620), page.px(30)
    y0 = H - margin - page.px(150)
    x = margin
    if legend.get("classes"):
        from matplotlib.patches import Rectangle
        rows = 4
        cols = int(np.ceil(len(legend["classes"]) / rows))
        sw, row_h, col_w = page.px(40), page.px(44), page.px(430)
        top = y0 + bar_h - rows * row_h + page.px(8)
        page.text(x, top - page.px(44), legend.get("classes_title", "Land cover"), size=30, color=text, va="bottom")
        for i, (label, colour) in enumerate(legend["classes"]):
            cx, cy = x + (i // rows) * col_w, top + (i % rows) * row_h
            page.ax.add_patch(Rectangle((cx, cy), sw, sw * 0.7, facecolor=colour, edgecolor="none", zorder=9))
            page.text(cx + sw + page.px(14), cy + sw * 0.35, label, size=26, color=text, va="center")
        x += cols * col_w + page.px(80)
    else:
        strip, lo, hi = legend["land"]
        page.text(x, y0 - page.px(44), "Land elevation (m)", size=30, color=text, va="bottom")
        page.gradient_bar(x, y0, bar_w, bar_h, strip, edge=muted)
        for frac, value in ((0.0, lo), (0.5, (lo + hi) / 2), (1.0, hi)):
            page.text(x + frac * bar_w, y0 + bar_h + page.px(10), f"{value:.0f}", size=26, color=muted,
                      ha="center" if 0 < frac < 1 else ("left" if frac == 0 else "right"), va="top")
        x += bar_w + page.px(120)
    strip, lo, hi = legend["sea"]
    page.text(x, y0 - page.px(44), "Sea floor depth (m)", size=30, color=text, va="bottom")
    page.gradient_bar(x, y0, bar_w, bar_h, strip, edge=muted)
    for frac, value in ((0.0, lo), (0.5, (lo + hi) / 2), (1.0, hi)):
        page.text(x + frac * bar_w, y0 + bar_h + page.px(10), f"{abs(value):.0f}", size=26, color=muted,
                  ha="center" if 0 < frac < 1 else ("left" if frac == 0 else "right"), va="top")
    if legend.get("note"):
        page.text(margin, y0 + bar_h + page.px(70), legend["note"], size=26, color=muted, va="top")
    for i, (label, kind, colour) in enumerate(legend.get("symbols", [])):
        # under the credits, top right: the sky is free there
        sx, sy = W - margin - page.px(1050), margin + page.px(10) + (len(credits) + 1) * page.px(34) + i * page.px(46)
        if kind == "pin":
            page.segment(sx, sy - page.px(16), sx, sy + page.px(16), color=colour, width=3, ls=(0, (2, 2)))
            page.marker(sx, sy - page.px(16), size=12, face="white", edge=colour, edge_width=2)
            page.marker(sx, sy, size=14, face=colour, edge="white")
            page.marker(sx, sy + page.px(16), size=10, face="white", edge=colour, edge_width=2, shape="s")
        else:
            page.marker(sx, sy, size=16, face=colour, edge="white", shape=kind)
        page.text(sx + page.px(34), sy, label, size=26, color=text, va="center")

    y = margin + page.px(10)
    for line in credits:
        page.text(W - margin, y, line, size=24, color=muted, ha="right", va="top")
        y += page.px(34)
