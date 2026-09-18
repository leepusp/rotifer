#!/usr/bin/env python3

__version__ = "0.01"

### Import rotifer package
import sys
import os
sys.path.insert(0, os.path.join(os.getcwd(), '../..'))

### Import core cli
import rotifer.core.cli as corecli

### Other packages
import argparse
import sys

### create arguments for sequences


class action:
    '''
    Rparser custom actions!
    '''
    def openload(self, parser, largs):
        if sys.stdin.isatty():
            if len(largs) == 0:
                parser.error('No input!')
        elif not '-' in largs:
            largs.insert(0, sys.stdin)
        sys.stdin


### Sequence actions
