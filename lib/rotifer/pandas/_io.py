"""Read text from, and write text to, a path, a file-like object or a string."""

import os


def _read_source(source):
    """Text from a file-like object, an existing single-line path, or the text itself."""
    if hasattr(source, 'read'):
        return source.read()
    if isinstance(source, (str, os.PathLike)) and '\n' not in str(source) and os.path.exists(source):
        with open(source) as fh:
            return fh.read()
    return source


def _write_sink(text, buf):
    """Return `text` when buf is None, otherwise write it to a path or file-like and return None."""
    if buf is None:
        return text
    if isinstance(buf, (str, os.PathLike)):
        with open(buf, 'w') as fh:
            fh.write(text)
        return None
    buf.write(text) # file-like object; caller closes it
    return None
