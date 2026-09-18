#!/usr/bin/env python3

import argparse
import sys
import copy as _copy
### Add arguments ###

def parser(description =None): 
    return argparse.ArgumentParser(add_help = False, description = description)


### General Actions ###

class action:
    '''
    Rparser custom actions!
    '''
    def openread(self,parser,largs):
        '''
        Open, read, and close a file. Remove duplicates
        '''
        # Process sys.stdin
        if sys.stdin.isatty():
            if len(largs) == 0:
                parser.error("No input!")
        elif not "-" in largs:
            largs.insert(0,sys.stdin)

        # Process sys.argv
        myset=set()
        for iohandle in largs:
            try:
                iohandle = open(iohandle, mode="r")
                myset=myset.union(iohandle.read().splitlines())
                iohandle.close
            except:
                try:
                    myset=myset.union(iohandle.read().splitlines())
                except:
                    if iohandle == "-":
                        if sys.stdin.isatty():
                            parser.error("No input!")
                        else:
                            myset=myset.union(sys.stdin.read().splitlines())
                    else:
                        myset.add(iohandle)
        return(myset)

    def openread2(self, parser, largs):
        '''
        Open, read, and close a file. Keep duplicates.
        '''
        # Process sys.stdin
        if sys.stdin.isatty():
            if len(largs) == 0:
                parser.error("No input!")
        elif not "-" in largs:
            largs.insert(0,sys.stdin)

        # Process sys.argv
        myset= []
        for iohandle in largs:
            try:
                iohandle = open(iohandle, mode="r")
                myset.extend(iohandle.read().splitlines())
                iohandle.close
            except:
                try:
                    myset.extend(iohandle.read().splitlines())
                except:
                    if iohandle == "-":
                        if sys.stdin.isatty():
                            parser.error("No input!")
                        else:
                            myset.extend(sys.stdin.read().splitlines())
                    else:
                        myset.extend(iohandle)
        return(myset)

    def openload(self, parser, largs):
        if sys.stdin.isatty():
            if len(largs) == 0:
                parser.error("No input!")
        elif not "-" in largs:
            largs.insert(0,sys.stdin)

        # Process sys.argv
        mylist = list()
        for iohandle in largs:
            try:
                iohandle = open(iohandle, mode="r")
                mylist.append(iohandle)
            except:
                try:
                    mylist.append(iohandle)
                except:
                    if iohandle == "-":
                        if sys.stdin.isatty():
                            parser.error("No input!")
                        else:
                            mylist.append(iohandle)
                    else:
                        mylist.append(iohandle)
        return mylist

    class autoload(argparse.Action):
        '''
        An action to decide file type, open, and load the file into memory, close after
        So far it decide if string, file, stdin
        Optional parameters:
            duplicates = Boolean [True/False]
        '''
        def __init__(self, option_strings, file_type = 'list', duplicates = True,
                      *args, **kwargs):
            self._add1 = file_type
            self._add2 = duplicates
            super(action.autoload, self).__init__(option_strings=option_strings,
                                           *args, **kwargs)
        def __call__(self, parser, namespace, values, option_string = None):
            #TODO insert new function for the var file type
            #print(self._add1)
            if self._add2 == True:
                f = action.openread(self,parser, values)
            if self._add2 == False:
                f = action.openread2(self,parser, values)
            setattr(namespace, self.dest, f)

    class autoopen(argparse.Action):
        '''
        An action to decide file type, and open.
        '''
        def __call__(self, parser, namespace, values, option_string = None):
            f = action.openload(self, parser, values)
            setattr(namespace, self.dest, f)

    class dict_options(argparse._AppendAction):
        '''
        An action to use with nargs = '1',
        will return a list of dict, separator is '='
        Usage in command line
        -arg my=1 -arg yours=2
        Result:
        args.arg = [{my : 1}, {yours : 2}]
        '''
        def __init__(self, option_strings, default = {}, *args, **kwargs):
            super(action.dict_options, self).__init__(option_strings = option_strings,
                                                    *args, **kwargs)
        def __call__(self, parser, namespace, values, option_string = None):
            items = argparse._copy.copy(argparse._ensure_value(namespace, self.dest, []))
            _source, _target = values[0].split('=')
            dd = {}
            dd[_source] = _target
            items.append(dd)
            setattr(namespace, self.dest, items)

    class extend(argparse.Action):

        def __call__(self, parser, namespace, values, option_string=None):
            items = getattr(namespace, self.dest) or []
            items.extend(values)
            setattr(namespace, self.dest, items)

###

def parseargs(add_help = True,parents = []):
    '''
    Combine multiple argparse.ArgumentParser
    '''
    if len(parents) == 1:
        return argparse.ArgumentParser(parents=parents, add_help= add_help).parse_args()
    else:
        s = argparse.ArgumentParser(parents=parents, add_help = add_help)
        return s.parse_args()
####

def _ensure_value(namespace, name, value):
    if getattr(namespace, name, None) is None:
        setattr(namespace, name, value)
    return getattr(namespace, name)
