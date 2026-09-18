#!/usr/bin/env python3

# Copyright 2020 by Robson F. de Souza.  All rights reserved.
# This file is part of the Rotifer distribution and governed by 
# the "BSD 3-Clause License".
#
# Please see the LICENSE file that should have been included as part of this
# package.

# Import external modules
import os
import fileinput
from fileinput import *

_state = None

class FileInput(fileinput.FileInput):
    def __init__(self, files=None, inplace=False, backup='', *, mode='r', openhook=None, delete=False):
        self.files   = files if isinstance(files,list) else [ files ]
        self._delete = delete
        super(__class__, self).__init__(files=files, inplace=inplace, backup=backup, mode=mode, openhook=openhook)

    def close(self):
        super(__class__, self).close()
        if self._delete:
            for f in self.files:
                if os.path.exists(f):
                    os.remove(f)
        _state = None

    def __del__(self):
        self.close()

    def read(self, bufsize):
        """
        Ugly hack to make rotifer.ncbi.io.fileinput compatible with Bio.SeqIO
        Should be replaced with a _io.read() compatible implementation
        """
        return ""

def input(files=None, inplace=False, backup="", *, mode="r", openhook=None, delete=False):
    """Return an instance of the FileInput class, which can be iterated.

    The parameters are passed to the constructor of the FileInput class.
    The returned instance, in addition to being an iterator,
    keeps global state for the functions of this module,.
    """
    if fileinput._state and fileinput._state._file:
            raise RuntimeError("input() already active")
    _state = FileInput(files, inplace, backup, mode=mode, openhook=openhook, delete=delete)
    return _state

# Internal method to use when returning concatenated file
# streams (handlers) with fileinput

if __name__ == '__main__':
    pass
