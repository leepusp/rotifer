# Fonts

The site uses two typefaces, both served from this folder so no page
makes a request to a font service:

- `chivo-variable.woff2`: Chivo, for text and headings.
- `chivo-mono-variable.woff2`: Chivo Mono, for code, commands and
  identifiers.

Each file is a variable font with the full weight axis, subset to Latin
and Latin Extended-A. They are made from the TTFs in Google Fonts'
repository by `docs/_tools/build_fonts.py`; its docstring says where
to download them and how to run it.

Both fonts are by the Chivo Project Authors (Omnibus-Type) and are
licensed under the SIL Open Font License 1.1, in `OFL.txt`.
