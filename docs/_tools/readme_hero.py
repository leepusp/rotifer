"""Refresh the landing hero screenshot shown at the top of README.md.

The README shows a picture of the site's hero (docs/_static/readme/hero.png).
Run this after changing the hero, with the docs being served locally:

    cd docs && make livehtml          # or any server on the built site
    python docs/_tools/readme_hero.py http://127.0.0.1:8000/

Needs Google Chrome (or Chromium) on PATH and Pillow. The page is captured
at 2x, cropped to the hero (the rows whose left edge is the hero's ink
color), scaled to 1600px wide and saved as a 128-color PNG.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from PIL import Image

OUT = Path(__file__).resolve().parents[1] / "_static" / "readme" / "hero.png"
INK = (0x10, 0x14, 0x13)  # --rot-ink in _static/css/tokens.css
WIDTH = 1600


def chrome():
    for name in ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser"):
        if shutil.which(name):
            return name
    sys.exit("Chrome or Chromium not found on PATH")


def capture(url, dest):
    subprocess.run(
        [chrome(), "--headless=new", "--no-sandbox", "--disable-gpu", "--hide-scrollbars",
         "--virtual-time-budget=4000", "--force-device-scale-factor=2",
         "--window-size=1280,900", f"--screenshot={dest}", url],
        check=True, capture_output=True,
    )


def hero_rows(im):
    """First and last row whose left edge is the hero's ink color."""
    rows = [y for y in range(im.height)
            if max(abs(a - b) for a, b in zip(im.getpixel((4, y)), INK)) < 6]
    if not rows:
        sys.exit("hero not found in the screenshot")
    return rows[0], rows[-1]


def main(url):
    with tempfile.TemporaryDirectory() as tmp:
        shot = os.path.join(tmp, "page.png")
        capture(url, shot)
        im = Image.open(shot).convert("RGB")
    top, bottom = hero_rows(im)
    crop = im.crop((0, top, im.width, bottom + 1))
    crop = crop.resize((WIDTH, round(crop.height * WIDTH / crop.width)), Image.LANCZOS)
    crop = crop.quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    crop.save(OUT, optimize=True)
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8000/")
