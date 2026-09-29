"""
Catalog collections of local data files.

This module provides :class:`FileCollection`, a small helper that
scans a directory tree for data files (sequence alignments, by
default) and exposes the matching paths as a dataframe.

Note: this module currently fails to import because it calls
``loadConfig`` without importing it; see ``docs/OPEN_QUESTIONS.md``
in the repository.
"""

import os
import hashlib
import pathlib
import pandas as pd
import rotifer
logger = rotifer.logging.getLogger(__name__)
_defaults = {
    "basedir": "/databases/fadb",
    "checksum": True,
    "pattern": "*.fa",
    "recursive": True,
    "ignore": [],
}
config = loadConfig(__name__, defaults = _defaults)

# Classes
# Loading our alignments
class FileCollection():
    """
    Index the data files found under a directory tree.

    Scans ``basedir`` for files matching one or more glob patterns
    and collects the results in a dataframe, one row per file. Extra
    keyword arguments describe companion files to look up for each
    entry, each becoming an additional column.

    Parameters
    ----------
    basedir : path-like, default '/databases/fadb'
        Root of the directory tree to scan.
    pattern : str or list of str, default '*.fa'
        Glob pattern(s) identifying the data files. The matched
        suffix is stripped from each file name to build the entry
        name.
    recursive : bool, default True
        Whether to descend into subdirectories.
    checksum : bool, default True
        Whether to add a ``checksum`` column.
    ignore : list of str, default []
        Entry names to skip.
    **kwargs : dict
        Companion file groups, one per resulting column. Each value
        is a dictionary with a ``pattern`` key and an optional
        ``basedir`` key, searched relative to the entry name. A
        group that matches nothing yields None, a single match
        yields the path itself and several matches yield a list.

    Attributes
    ----------
    df : pandas.DataFrame
        One row per data file, with the ``name`` and ``path``
        columns plus any companion columns named by ``kwargs``.

    Notes
    -----
    Passing ``checksum=True`` (the default) appends a ``checksum``
    column name without producing a matching value, so building the
    dataframe raises a length mismatch; ``hashlib`` is imported but
    never used. See ``docs/OPEN_QUESTIONS.md``.

    Examples
    --------
    >>> from rotifer.db.local.core import FileCollection
    >>> fc = FileCollection(basedir="/databases/fadb", checksum=False)  # doctest: +SKIP
    >>> fc.df.head()  # doctest: +SKIP
    """

    def __init__(
            self,
            basedir   = config["basedir"],
            pattern   = config["pattern"],
            recursive = config["recursive"],
            checksum  = config["checksum"],
            ignore    = config["ignore"],
            **kwargs
        ):
        if "kwargs" in config and config["kwargs"]:
            kwargs = { **config["kwargs"], **kwargs }
        if isinstance(pattern,str):
            pattern = [pattern]
        basedir = pathlib.Path(basedir)
        if not basedir.exists():
            logger.error(f'No directory named {basedir.name}')
            return None

        # Find alignments
        db = []
        alignments = dict()
        for extension in pattern:
            # Load data
            if recursive:
                it = basedir.rglob(extension)
            else:
                it = basedir.glob(extension)
            for x in it:
                if extension[0] == "*":
                    suffix = extension[1:]
                else:
                    suffix = extension
                basename = x.name.replace(f"{suffix}","")
                if basename in ignore:
                    continue

                # Load other paths
                otherfiles = []
                for colname, other in kwargs.items():
                    patt = os.path.join(basename, other["pattern"])
                    basepath = pathlib.Path(other["basedir"]) if "basedir" in other else basedir
                    if recursive:
                        files = list(basepath.rglob(patt))
                    else:
                        files = list(basepath.glob(patt))
                    files = [ y.as_posix() for y in files ]
                    if len(files) == 0:
                        files = None
                    elif len(files) == 1:
                        files = files[0]
                    otherfiles.append(files)

                # Add data to stack
                db.append((basename, x.as_posix(), *otherfiles))

        # Create internal dataframe for paths
        colnames = ['name','path']
        if checksum:
            colnames.append('checksum')
        db = pd.DataFrame(db, columns=colnames + list(kwargs.keys()))

        # Store data in object
        self.df = db

    def __setitem__(self,key,value):
        """
        Replace a row of the underlying dataframe.

        Parameters
        ----------
        key : label
            Row label, as accepted by :meth:`pandas.DataFrame.loc`.
        value : object
            Replacement row.
        """
        self.df.loc[key] = value

    def __getitem__(self,key):
        """
        Fetch a row of the underlying dataframe.

        Parameters
        ----------
        key : label
            Row label, as accepted by :meth:`pandas.DataFrame.loc`.

        Returns
        -------
        pandas.Series
            The row describing one indexed file.
        """
        return self.df.loc[key]

    def keys(self):
        """
        List the row labels of the underlying dataframe.

        Returns
        -------
        list
            The dataframe index, as a list.
        """
        return self.df.index.tolist()

# Is this library being used as a script?
if __name__ == '__main__':
    pass
