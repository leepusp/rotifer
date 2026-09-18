#!/usr/bin/env python3

import os
import sys
import types
import typing
import shelve
from copy import deepcopy
from datetime import datetime as dt

from rotifer.core import config as CoreConfig
from rotifer.core import logger as rcl
logger = rcl.get_logger(__name__)

"""
Rotifer core functions
"""

__version__ = 0.15

def chunks(list_to_chunk, chunk_size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(list_to_chunk), chunk_size):
        yield list_to_chunk[i:i + chunk_size]

def who_is_calling(obj):
    return str(type(obj)).replace("<class '","").replace("'>","")

def get_some_info():
    import inspect
    # get the call frame of the calling method
    frame = inspect.currentframe().f_back
    try:
        # find the name of the first variable in the calling
        # function - which is hopefully the "self"
        codeobj = frame.f_code
        try:
            self_name = codeobj.co_varnames[0]
        except IndexError:
            return None

        # try to access the caller's "self"
        try:
            self_obj = frame.f_locals[self_name]
        except KeyError:
            return None

        # check if the calling function is really a method
        self_type = type(self_obj)
        func_name = codeobj.co_name

        # iterate through all classes in the MRO
        for cls in self_type.__mro__:
            # see if this class has a method with the name
            # we're looking for
            try:
                method = vars(cls)[func_name]
            except KeyError:
                continue

            # unwrap the method just in case there are any decorators
            try:
                method = inspect.unwrap(method)
            except ValueError:
                pass

            # see if this is the method that called us
            if getattr(method, '__code__', None) is codeobj:
                return self_type.__name__

        # if we didn't find a matching method, return None
        return None

    finally:
        # make sure to clean up the frame at the end to avoid ref cycles
        del frame

def open_compressed(filename, mode='r', encoding='utf8'):
    """
    Open a bzip2 or a gzip file as if it was uncompressed.

    Usage:
      import rotifer.core.functions as rcf
      fh = open_compressed(filename, "rb")

    Parameters:
      filename:
        Path to the compressed file
      mode: string, default 'rb'
        Access mode
      enconding: string, default 'utf-8'
        File enconding
    """
    ext = os.path.splitext(filename)[1]
    if (mode == 'r') or not mode:
        mode = 'rt'
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode, encoding=encoding)
    elif ext == '.bz2':
        import bz2
        return bz2.open(filename, mode, encoding=encoding)
    else:
        return open(filename, mode, encoding=encoding)

def vmsg(message = ''):
    """
    Output a message to stderr.

    Three options available:
    1-) Input a string
    2-) Input a list with one element
    3-) Input a list with n elements

    If a list with more than one element is passed,
    the first element will be printed with timelog the other lines will be printed in a new line.
    -----------
    INPUT:
    message: A message in the format, str or list.
    -----------
    OUTPUT:
    A message in the format for option (1) and (2)
    ## [DAY/MONTH Hour:Minute:Second] <message>
    or

    A message in the format for option (3)
    ## [DAY/MONTH Hour:Minute:Second] <message>
    ### <message>
    ### ...
    -----------
    Example:
    import rotifer.core.functions as cf

    if verbose:
        cd.msg('Small test')
    >> ## [12/03/18 09:50:54] Small test
    """

    now = dt.now().strftime('[%D %H:%M:%S]')
    if isinstance(message, list):
        if len(message) == 0:
            message = ''.join(message)
        else:
            m = message[0] +'\n'
            m2 = ''
            for e in message[1:]:
                m2 += '### ' + e + '\n'
            message = m+m2.rstrip('\n')
    sys.stderr.write('## {0} {1}\n'.format(now, message))
    sys.stderr.flush()

class rdoc:
    '''
    Input a file with rotifer corecli
    Ouput a markdown file
    The first parameter f is the script path
    the md object is created containing a raw markdown
    Later it can be printed to stdout using writer() function.
    ---------
    EXAMPLE:
    import rotifer.core.functions as cf
    md = cf.rdoc(os.path.basename(__file__))
    # Raw markdown
    print(md)

    # Pretty markdown
    md.writer()
    '''

    def __init__(self, f):
        self.f = f
        self.md = self._markdown()

    def _markdown(self):
        __version__ = 0
        __authors__ = ''
        self.f =open(self.f).read().splitlines()
        s2 = []
        for line in self.f:
            if line.startswith('from') or line.startswith('import'):
                s2.append(line)

        docs = []
        for x in range(len(self.f)):
            if self.f[x].startswith('__doc__'):
                for y in range(x+1, len(self.f)):
                    if "'''" in self.f[y]:
                        break
                    else:
                        docs.append(self.f[y])

        dc = {}
        for ele in docs:
            if ele.isupper():
                k = ele
                dc[k] = ''
            else:
                dc[k] += ele+'\n'

        s = []
        for x in range(len(self.f)):
            if 'def parse_cli()' in self.f[x]:
                while True:
                    if self.f[x].lstrip().startswith('return') or self.f[x].lstrip().startswith('args'):
                        break
                    else:
                        x+=1
                        s.append(self.f[x].strip())

        fi = exec('\n'.join(s2), locals(), globals())
        st = '\n'.join(s[:-1])
        exec(st, locals(), globals())
        ls_info = []

        for e in parser_merged._actions:
            long_alias = []
            short_alias = []
            default = []
            arg_type = []

            if e.default != '==SUPPRESS==' and \
                '--fun' not in e.option_strings and \
                '--configdump' not in e.option_strings and \
                '--configfile' not in e.option_strings:
                if not e.option_strings: # empty alias
                    long_alias.append(e.dest)
                    short_alias.append('')
                    default = e.default if not isinstance( e.default, type(None)) else 'None'
                    arg_type = e.type if not isinstance(e.type, type(None)) else 'None'
                    helper = e.help if not isinstance(e.type, type(None)) else ''
                    ls_info.append([long_alias,
                                    short_alias,
                                    default,
                                    arg_type,
                                    helper])

                else:
                    longest_opt = max([len(x) for x in e.option_strings])
                    for option in e.option_strings:
                        if len(option) == longest_opt or len(option)+1 == longest_opt:
                            long_alias.append(option)
                        else:
                            if option:
                                s_option = option
                            else:
                                s_option = ''
                            short_alias.append(s_option)
                    default = e.default if not isinstance( e.default, type(None)) else 'None'
                    arg_type = e.type if not isinstance(e.type, type(None)) else 'None'
                    helper = e.help if not isinstance(e.help, type(None)) else ''
                    ls_info.append([long_alias,
                                    short_alias,
                                    default,
                                    arg_type,
                                    helper])

        toprint = '## Rotifer\n\n'
        head_doc = ['NAME', 'SYNOPSIS', 'DESCRIPTION', 'AUTHOR']
        for descriptor in head_doc:
            if descriptor in dc.keys():
                toprint += '## {0}\n'.format(descriptor)
                toprint += dc[descriptor] +'\n'
            else:
                toprint += '## {0}\n'.format(descriptor)

        toprint += '''
## OPTIONS
'''
        for e in ls_info:
            aliases = ' '.join(e[0]) +', ' + ' '.join(e[1])
            aliases_msg = '**'+ aliases.rstrip(', ')+'**'+'\n'
            default_msg = 'default: '+ str(e[2]) +'\n'
            type_msg = e[3]
            helper = e[4]
            toprint += '\n'.join([aliases_msg, helper, '',default_msg, ''])


        toprint += '''
## Program options summary


|**Long name**|**Aliases**|**Type**|
|-----:|----:|:-----:|
'''

        for e in ls_info:
            long_name = ' '.join(e[0])
            short_name = ' '.join(e[1])
            # type_name = e[3]
            type_name = str(e[3].__class__.__name__) if str(e[3].__class__.__name__) != 'type' else 'None'
            toprint += '|{0}|{1}|{2}|\n'.format(long_name,short_name,type_name)

        toprint = toprint.rstrip('\n')
        return toprint

    def writer(self, output= 'markdown', theme = '641.4291'):
        '''
        '''
        output = self._checker(output)
        self.theme = theme
        switcher = {'md': self._mdwriter,
                    'raw': self._rawwriter,
                    'man': self._manwriter}
        switcher[output]()

    def _checker(self, output):
        if output in ['md', 'markdown']:
            return 'md'
        if output in ['raw']:
            return 'raw'
        if output in ['man']:
            return 'man'
    def _mdwriter(self):
        mdv.term_columns = 60
        print(mdv.main(self.md, theme = self.theme))

    def _rawwriter(self):
        print(self.md)
    def _manwriter(self):
        pass

class dotdict(dict):
    """A class implementing dot annotation access to dictionary attributes
    Instead of using dc['key'] dictionary values can be accessed using dc.key.
    -----------
    PARAMETERS:
    input a python dictionary
    -----------
    Example:
    import rotifer.core.functions as cf

    dc = dotdict({'key1': 'value1',
                  'key2': 'value2'})

    dc.key1
    dc.key2
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def _flatten(*args):
    '''
    Flat a list of arguments
    Up to two levels
    list of list [[]]
    ---------
    PARAMETERS:
    args: list, or singleton
        a list of input
    ---------
    RETURNS:
    A flatten list
    --------
    EXAMPLE:
    args = a
    >> [a]
    args = [a,b]
    >> [a,b]
    args = [[a],b]
    >> [a,b]
    args = [[a,b], [c]]
    >> [a,b,c]
    '''
    ls = []
    if isinstance(args, (list, tuple)):
        for ele1 in args:
            if isinstance(ele1, (list, tuple)):
                for ele2 in ele1:
                    if isinstance(ele2, (list, tuple)):
                        for ele3 in ele2:
                            if isinstance(ele3, (list, tuple)):
                                for ele4 in ele3:
                                    ls.append(ele4)
                            else:
                                ls.append(ele3)
                    else:
                        ls.append(ele2)
            else:
                ls.append(ele1)
    else:
        ls.append(args)
    return ls


def not_kwargs(dict_args, key, value):
    """ Small function to check if an element is not in a dictionary.
    It is useful to set a default value to functions that exceeds the
    number of 5 arguments
    """
    if key in dict_args.keys():
        return dict_args[key]
    return value

def findDataFiles(load, all=False):
    '''
    This routine locates data files under Rotifer's standard locations.

    User-supplied qualified paths, such as load='./config.yml', are
    accepted but this method is mostly useful when loading files from
    the following standard locations:

      1) The user folder ~/.rotifer/share
      2) The path set by the environment variable ROTIFER_DATA
      3) <rotifer installation path>/share/rotifer

    If the load parameter starts with the prefix ':' and subdirectory
    names are separated by dots ('.'), files will be loaded from the
    standard locations.
    
    Example:
    
      load=":taxonomy.taxonomy.txt"
      
    means
    
      ~/.rotifer/share/taxonomy/taxonomy.txt

    and the path above is used if it exists. Otherwise, this
    example will map to

      ${ROTIFER_DATA}/taxonomy/taxonomy.txt

    if the ROTIFER_DATA variable is set and a taxonomy.txt file
    exists in the taxonomy subdirectory under ${ROTIFER_DATA}.

    Finally, if neither of the paths above exist, the file will be
    searched for under rotifer's installation directory, subdirectory
    share/rotifer:

      <path to rotifer>/share/rotifer/taxonomy/taxonomy.txt
    
    Notice that files from the user standard location take precedence
    over files under rotifer's installation path.
    '''
    import rotifer

    # Make sure input is a list
    if not isinstance(load,list):
        load = [ load ]

    files = []
    for f in load:
        # Rename if needed
        if f[0:8] == "rotifer.":
            f = f.replace("rotifer.",":",1)

        # Standard files
        if f[0] != ":" and os.path.exists(f):
                files.append(f)

        # Process path
        g = []
        f = f.replace(':','')
        if os.path.sep in f:
            f = f.split(os.path.sep)
        else:
            f = f.split('.')
            if len(f) > 1:
                g = deepcopy(f)
                g[-2] = g[-2] + '.' + g[-1]
                del(g[-1])

        # User path
        user_path = CoreConfig['userDataDirectory']
        if os.path.exists(os.path.join(user_path, *f)):
            files.append(os.path.join(user_path, *f))
        elif len(g) > 0 and os.path.exists(os.path.join(user_path, *g)):
            files.append(os.path.join(user_path, *g))

        # Databases
        databases = rotifer.config['data']
        if os.path.exists(os.path.join(databases, *f)):
            files.append(os.path.join(databases, *f))
        elif len(g) > 0 and os.path.exists(os.path.join(databases, *g)):
            files.append(os.path.join(databases, *g))

        # System path
        system_path = CoreConfig['baseDataDirectory']
        if os.path.exists(os.path.join(system_path, *f)):
            files.append(os.path.join(system_path, *f))
        elif len(g) > 0 and os.path.exists(os.path.join(system_path, *g)):
            files.append(os.path.join(system_path, *g))

    # Return
    if files and not all:
        files = files[0]
    return files

def loadConfig(filepath, user_path=CoreConfig['userConfig'], system_path=CoreConfig['baseConfig'], defaults={}):
    '''
    This routine loads configuration parameters from YAML file(s).

    If the load parameter starts with ':' and subdirectory names
    are separated by dots ('.'), the YAML file is searched in the
    user_path directory or, secondly, in the system_path directory.

    Example:
    
      conf = loadConfig(":db.pfam")

    means

      <user_path>/db/pfam.yml

    if <user_path>/db/pfam.yml exists.

    Therefore, if user_path is "~/.rotifer/etc", this
    example translates to

      conf = "~/rotifer/etc/db/pfam.yaml"

    If the location above doesn't exist, the path

      <system_path>/db/pfam.yml
    
    is used, i.e. the user path always takes precedence over
    Rotifer's installation path so that user customizations
    always apply first.
    '''
    from yaml import load as yaml_load
    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader

    # Load defaults
    from copy import deepcopy
    config = deepcopy(defaults)

    # Parse input
    if filepath[0] != ":":
        if os.path.exists(filepath):
            try:
                loaded = yaml_load(open(filepath), Loader=Loader)
                config.update(loaded)
                return config
            except:
                logger.error(f'Error while parsing file {load}: {sys.exc_info[1]}')
                return config
        elif filepath[0:8] == "rotifer.":
            expand_load = filepath[8:]
        else:
            logger.error(f'File {filepath} not found.')
            return config
    else:
        expand_load = filepath.replace(':','')

    # System path
    if os.path.exists(system_path):
        try:
            config.update(yaml_search(expand_load, system_path))
        except:
            logger.error(f'Error while parsing file {os.path.join(system_path,expand_load)}: {sys.exc_info}', file=sys.stderr)

    # User path
    if os.path.exists(user_path):
        try:
            config.update(yaml_search(expand_load, user_path))
        except:
            logger.error(f'Error while parsing file {os.path.join(user_path,expand_load)}: {sys.exc_info}', file=sys.stderr)

    # Return
    return config

def yaml_search(string, path):
    from yaml import load as yaml_load
    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader

    # If a file matches path + exts, use it!
    exts = ['.yml','.yaml','.config']
    ls = string.split('.')
    for ext in exts:
        target = os.path.join(path,*ls) + ext
        if os.path.exists(target):
            return yaml_load(open(target), Loader=Loader)

    # This loop searches for entries within files
    # Do we need it?
    last = None
    data = {}
    query = ls[-1]
    while ls:

        # Check if any file or directory matches
        if os.path.isdir(path):
            contents = os.listdir(path)
            targets = [ ls[0] + x for x in exts ]
            targets.append(ls[0])
            for target in targets: # Search order: exts, ls[0] (directory)
                if target in contents:
                    path = os.path.join(path, target)
                    break

        # Load matching file that matches an extension
        if os.path.isfile(path) and [x for x in exts if path.endswith(x)]:
            if not data:
                data = yaml_load(open(path), Loader=Loader)
                #print(f'Loading file {path}', file=sys.stderr)
                last = ls[0]
            # Search for matching elements in data
            if isinstance(data, dict):
                if ls[0] in data.keys():
                    data = data[ls[0]]
                    last = ls[0]
            elif isinstance(data,list):
                try:
                    data = data[int(ls[0])]
                    last = ls[0]
                except:
                    pass

        # Move stack
        del ls[0]
    #while ls:

    # Return
    if query == last:
        return data
    else:
        return {}

def loadClasses(load, user_path=os.path.join(CoreConfig['user'],'lib'), system_path=os.path.join(CoreConfig['base'],'lib','rotifer')):
    '''
    This function list all classes and methods
    The main problem is the import
    '''
    import inspect

    classes = {}

    search_file = ''
    yml = ''

    old = os.getcwd()

    if load.startswith(':'):
        yml = loadConfig(load, user_path = user_path,
                         system_path = system_path)

    if yml:
        for path in yml:
            if not search_file:
                os.chdir(path)
                sys.path.insert(0, path)
                for e in os.listdir(path):
                    if e.endswith('.py'):
                        modname = e[:-3]
                        sub = __import__(modname)

                        methods_and_classes = [x for x in dir(sub) if not x.startswith('__')]
                        for method in methods_and_classes:
                            if inspect.isclass(getattr(sub, method)):
                                classes[method] = getattr(sub, method)
                            elif inspect.isfunction(getattr(sub, method)):
                                classes[method] = getattr(sub,method)
                            else:
                                pass

    else:
        if os.path.sep in load:
            load = load.replace(os.path.sep, '.')

        ls = load.split('.')
        if ls[0].startswith('~'):
            path = os.path.expanduser('~')

        else:
            path = os.path.join(os.path.sep, ls[0])

        ls = ls[1:]

        while ls:
            if os.path.isdir( path):
                for ele in os.listdir(path):
                    if ele.startswith(ls[0]):
                        if os.path.isdir(os.path.join(path,ele)):
                            path = os.path.join(path, ele)
                        else:
                            pass
            elif ls[0]+'.py' in os.listdir(path):
                search_file = (ls[0])
                break

            del ls[0]

        if not search_file:
            os.chdir(path)
            sys.path.insert(0, path)
            for e in os.listdir(path):
                if e.endswith('.py'):
                    modname = e[:-3]
                    sub = __import__(modname)

                    methods_and_classes = [x for x in dir(sub) if not x.startswith('__')]
                    for method in methods_and_classes:
                        if inspect.isclass(getattr(sub, method)):
                            classes[method] = getattr(sub, method)
                        elif inspect.isfunction(getattr(sub, method)):
                            classes[method] = getattr(sub,method)
                        else:
                            pass

        if search_file:
            while ls:
                if ls[0]+'.py' in os.listdir(path):
                    load_class = ls[0]
                    os.chdir(path)
                    sys.path.insert(0, path)
                    sub = __import__(load_class)
                else:
                    try:
                        classes[ls[0]] = getattr(sub, ls[0])
                    except:
                        pass

                del ls[0]

            if not classes:
                methods_and_classes = [x for x in dir(sub) if not x.startswith('__')]
                for method in methods_and_classes:
                    if inspect.isclass(getattr(sub, method)):
                        classes[method] = getattr(sub, method)
                    elif inspect.isfunction(getattr(sub, method)):
                        classes[method] = getattr(sub,method)
                    else:
                        pass

        os.chdir(old)

    return classes

def openread(largs):
    '''
    Open, read, and close a file. Remove duplicates
    '''

    # Process sys.argv
    myset=set()
    for iohandle in largs:
        try:
            stdi = iohandle.read().splitlines()
        except:
            stdi = ''
        if ':' in iohandle:
            try:
                _ = open(loadConfig(iohandle))
                myset = myset.union(_.read().splitlines())
            except:
                sys.stderr.write(f'File not found: {iohandle}\n')
                sys.exit()
        else:
            try:
                iohandle = open(iohandle, mode="r")
                myset=myset.union(iohandle.read().splitlines())
                iohandle.close
            except:
                if stdi:
                    myset=myset.union(stdi)
                elif iohandle == "-":
                    myset=myset.union(sys.stdin.read().splitlines())
                else:
                    myset.add(iohandle)
    _ = [x for x in myset if x != '']
    return(_)

def import_path(path):
    '''
    Import Python code from any input file, whatever its extension.
    '''
    import importlib
    module_name = os.path.basename(path).replace('-', '_')
    spec = importlib.machinery.SourceFileLoader(module_name, path)
    spec = importlib.util.spec_from_loader(module_name, spec)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules[module_name] = module
    return module

def md5(filename):
    '''
    Calculate MD5 checksum for a file.
    '''
    # Copied from https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    import hashlib
    hash_md5 = hashlib.md5()
    with open(filename,"rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def _fasta_ls(*args):
    '''
    Read fasta file or fasta sequence and output a list, each element is a line in the list.
    -----------
    PARAMETERS:
    args: a list containing a fasta sequence or file
    -----------
    RETURNS:
    List: each element is a line of the fasta file
    -----------
    EXAMPLE:

    '''
    fasta_res = []
    for arg in args:
        if isinstance(arg, (list, tuple)):
            fasta_res.extend(list(arg))
        else:
            try:
                fasta_res += open(arg).read().splitlines()
            except:
                fasta_res += arg.split('\n')
    return fasta_res

def save_session(filename, exclude=['In','Out','exit','quit'], exclude_type=tuple([types.ModuleType,types.FunctionType,types.MethodType,typing.Type]), namespace=None):
    myshelve = shelve.open(filename)
    if namespace is None:
        namespace = sys.modules['__main__'].__dict__
    with shelve.open(filename) as myshelve:
        for y in myshelve:
            if y not in namespace:
                del(myshelve[y])
        for x in namespace:
            if x.startswith('_') or (x in exclude) or isinstance(namespace[x],exclude_type):
                continue
            try:
                myshelve[x] = namespace[x]
            except Exception as e:
                print(f'Will not save {x} ({type(namespace[x])}): {e}', file=sys.stderr)
        myshelve.sync()

def restore_session(filename, exclude=[], namespace=None):
    myshelve = shelve.open(filename)
    if namespace is None:
        namespace = sys.modules['__main__'].__dict__
    with shelve.open(filename) as myshelve:
        for x in myshelve.keys():
            if x in exclude:
                continue
            try:
                namespace[x] = myshelve[x]
            except Exception as e:
                print(f'Failed to restore {x}: {e}', file=sys.stderr)
