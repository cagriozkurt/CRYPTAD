"""
Assemble 5-panel conformational composite figure (Figure 9).

Layout:
  Row 1:  (a) DCCM BIN1 BAR          |  (b) BAR curvature
  Row 2:  (c) DCCM PICALM ANTH       |  (d) PIP2-loop
  Row 3:  (e) ΔDCCM PICALM (3-panel, full width)

Panel labels (a)–(e): white text with dark halo, placed 28 px from left, 48 px
from top of each panel — offset far enough from the top y-axis tick ("385").

Usage:
  python3 09_scripts/10_figures/03_assemble_conformational.py
  python3 09_scripts/10_figures/03_assemble_conformational.py --project-root /path/to/CRYPTAD
"""

import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

CANVAS_W   = 5000
DPI        = 300
LABEL_SIZE = 72
LABEL_X    = 28
LABEL_Y    = 48
HALO_OFF   = 2


def load(path: Path) -> Image.Image:
    return Image.open(path).convert("RGBA")


def resize_w(img: Image.Image, target_w: int) -> Image.Image:
    w, h = img.size
    return img.resize((target_w, int(h * target_w / w)), Image.LANCZOS)


def draw_label(img: Image.Image, letter: str,
               x: int = LABEL_X, y: int = LABEL_Y,
               size: int = LABEL_SIZE) -> Image.Image:
    draw = ImageDraw.Draw(img)
    # Try system fonts in order; fall back to PIL default
    for font_path in [
        "/System/Library/Fonts/Helvetica.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    ]:
        try:
            font = ImageFont.truetype(font_path, size)
            break
        except Exception:
            font = None
    if font is None:
        font = ImageFont.load_default()
    text = f"({letter})"
    for dx in range(-HALO_OFF, HALO_OFF + 1):
        for dy in range(-HALO_OFF, HALO_OFF + 1):
            if dx != 0 or dy != 0:
                draw.text((x + dx, y + dy), text, fill=(20, 20, 20, 230), font=font)
    draw.text((x, y), text, fill=(255, 255, 255, 255), font=font)
    return img


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD conformational composite figure (Figure 9)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    conf_dir     = project_root / "03_pocket_analysis" / "conformational"
    fig_dir      = project_root / "06_figures"
    sub_dir      = project_root / "08_manuscript" / "submission" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    sub_dir.mkdir(parents=True, exist_ok=True)

    panels = {
        "a": conf_dir / "dccm_S1_BIN1_BAR.png",
        "b": conf_dir / "bar_curvature.png",
        "c": conf_dir / "dccm_S3_PICALM_ANTH.png",
        "d": conf_dir / "pip2_loop.png",
        "e": conf_dir / "dccm_delta_S3_PICALM.png",
    }

    for key, path in panels.items():
        if not path.exists():
            raise FileNotFoundError(f"Panel ({key}) not found: {path}")

    imgs  = {k: load(v) for k, v in panels.items()}
    col_w = CANVAS_W // 2

    a = resize_w(imgs["a"], col_w)
    b = resize_w(imgs["b"], col_w)
    row1_h = max(a.height, b.height)

    c = resize_w(imgs["c"], col_w)
    d = resize_w(imgs["d"], col_w)
    row2_h = max(c.height, d.height)

    e      = resize_w(imgs["e"], CANVAS_W)
    row3_h = e.height

    total_h = row1_h + row2_h + row3_h
    canvas  = Image.new("RGBA", (CANVAS_W, total_h), (255, 255, 255, 255))

    canvas.paste(a, (0,      0))
    canvas.paste(b, (col_w,  0))
    canvas.paste(c, (0,      row1_h))
    canvas.paste(d, (col_w,  row1_h))
    canvas.paste(e, (0,      row1_h + row2_h))

    draw_label(canvas, "a", x=LABEL_X,         y=LABEL_Y)
    draw_label(canvas, "b", x=col_w + LABEL_X, y=LABEL_Y)
    draw_label(canvas, "c", x=LABEL_X,         y=row1_h + LABEL_Y)
    draw_label(canvas, "d", x=col_w + LABEL_X, y=row1_h + LABEL_Y)
    draw_label(canvas, "e", x=LABEL_X,         y=row1_h + row2_h + LABEL_Y)

    rgb     = canvas.convert("RGB")
    out_png = fig_dir / "fig_conformational.png"
    rgb.save(out_png, dpi=(DPI, DPI))
    print(f"Saved: {out_png.relative_to(project_root)}  ({rgb.width}×{rgb.height} px)")

    sub     = resize_w(rgb, 3300)
    sub_png = sub_dir / "Figure_09_conformational.png"
    sub.save(sub_png, dpi=(DPI, DPI))
    print(f"Saved: {sub_png.relative_to(project_root)}  ({sub.width}×{sub.height} px)")


if __name__ == "__main__":
    main()
