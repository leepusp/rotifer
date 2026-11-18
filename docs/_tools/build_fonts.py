"""Convert the Chivo and Chivo Mono variable TTFs into the woff2 files the site serves.

The site self-hosts its fonts (no external requests), so the woff2 files are
committed under docs/_static/fonts. Run this script only to refresh them:

    pip install fonttools brotli
    python docs/_tools/build_fonts.py path/to/Chivo[wght].ttf path/to/ChivoMono[wght].ttf

Download the TTFs from https://github.com/google/fonts/tree/main/ofl/chivo
and https://github.com/google/fonts/tree/main/ofl/chivomono. Both are under
the SIL Open Font License (docs/_static/fonts/OFL.txt).

The weight axis is kept, so one file per family covers every weight. Glyphs
are subset to Latin and Latin Extended-A, which covers English and Portuguese.
"""

import sys
from pathlib import Path

from fontTools import subset
from fontTools.ttLib import TTFont

OUT = Path(__file__).resolve().parents[1] / "_static" / "fonts"

# Basic Latin, Latin-1 Supplement, Latin Extended-A, and general punctuation
# (dashes, quotes, bullet, ellipsis), plus the middle dot used in the footer.
UNICODES = "U+0000-00FF,U+0100-017F,U+2000-206F,U+20AC,U+2122,U+2190-2193,U+2212"

NAMES = {"Chivo": "chivo-variable.woff2", "Chivo Mono": "chivo-mono-variable.woff2"}


def family(font):
    # Name ID 16 is the typographic family ("Chivo"); ID 1 can carry the
    # default instance's style ("Chivo Medium") in variable fonts.
    names = font["name"]
    return names.getDebugName(16) or names.getDebugName(1)


def build(src):
    font = TTFont(src)
    name = family(font)
    if name not in NAMES:
        sys.exit(f"{src}: unexpected family {name!r}")
    opts = subset.Options()
    opts.flavor = "woff2"
    opts.layout_features = ["*"]
    opts.name_IDs = ["*"]
    opts.notdef_outline = True
    sub = subset.Subsetter(opts)
    sub.populate(unicodes=subset.parse_unicodes(UNICODES))
    sub.subset(font)
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / NAMES[name]
    font.flavor = "woff2"
    font.save(dest)
    print(f"{src} -> {dest} ({dest.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for path in sys.argv[1:]:
        build(path)
