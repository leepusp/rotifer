# This module's documentation
__doc__ = '''Rotifer's core configuration and setup

Paths are resolved so that Rotifer works from a git checkout, from a
pip/conda installation and inside containers or on shared servers.

Environment variables:

  ROTIFER_BASE
    Root of a Rotifer tree holding etc/rotifer and share/rotifer.
    Overrides the configuration and data files shipped with the package.
  ROTIFER_USER
    User directory (default: ~/.rotifer). If it is not writable, a
    private directory under the system's temporary directory is used.
  ROTIFER_CACHE
    Cache directory (see rotifer.config['cache']).
'''

# Dependencies
import os
import sys

def _writable(path):
    '''True if path is, or can be created as, a writable directory.'''
    path = os.path.abspath(path)
    while not os.path.exists(path):
        parent = os.path.dirname(path)
        if parent == path:
            return False
        path = parent
    return os.path.isdir(path) and os.access(path, os.W_OK | os.X_OK)

def _get_base_dirs():
    '''Return the installation root and the system config and data directories.'''
    from os.path import expanduser, dirname, realpath, join, isdir

    # Configuration and data files shipped with the package (see MANIFEST.in)
    packageDir = dirname(dirname(realpath(__file__)))
    packageData = join(packageDir, "_data")

    # 1) Explicit override
    if os.environ.get('ROTIFER_BASE'):
        baseDir = realpath(expanduser(os.environ['ROTIFER_BASE']))
        return baseDir, join(baseDir,"etc","rotifer"), join(baseDir,"share","rotifer")

    # 2) Git checkout: <root>/lib/rotifer/core/__init__.py
    root = dirname(dirname(packageDir))
    if isdir(join(root,"etc","rotifer")) and isdir(join(root,"share","rotifer")):
        return root, join(packageData,"etc"), join(packageData,"share")

    # 3) Installed package (pip, conda, containers)
    return sys.prefix, join(packageData,"etc"), join(packageData,"share")

def _get_user_dir():
    '''Return the user directory, falling back to a temporary one if needed.'''
    from os.path import expanduser, realpath, join
    import tempfile
    userDir = realpath(expanduser(os.environ.get('ROTIFER_USER') or join('~','.rotifer')))
    if not _writable(userDir):
        fallback = join(tempfile.gettempdir(), f"rotifer-{os.getuid() if hasattr(os,'getuid') else 'user'}")
        from rotifer.core import logger as rcl
        rcl.get_logger(__name__).warning(f'User directory {userDir} is not writable, using {fallback}')
        userDir = fallback
    return userDir

# Core configuration
def _get_core_config():
    from os.path import dirname, realpath, join

    # Base directories
    baseDir, baseConfig, baseDataDirectory = _get_base_dirs()
    userDir = _get_user_dir()

    # Main dictionary: this is the most important part of this module!!!!
    config = {
        'base': baseDir,
        'baseConfig': baseConfig,
        'baseDataDirectory': baseDataDirectory,
        'package': dirname(dirname(realpath(__file__))),
        'user': userDir,
        'userConfig': join(userDir,"etc"),
        'userDataDirectory': join(userDir,"share"),
    }

    # Libraries
    libdir = join(userDir,"lib","python")
    if os.path.isdir(libdir) and libdir not in sys.path:
        sys.path.insert(0,libdir)

    return config

config = _get_core_config()
