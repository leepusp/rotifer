__doc__ = """
rrsw: download and unpack mirrors of public databases.

A table of targets (``rrsw.yml`` in ROTIFER's ``etc/rotifer/rrsw``
directory, or ``rrsw.<host>.yml`` when present) lists, for each target,
how to download it and how to process it::

    - [ NAME    , METHOD , SERVER          , PATH                                  , PROCESS , OPTIONS ]
    - [ uniprot , ftp    , ftp.uniprot.org , pub/databases/uniprot/current_release , unpack  ,         ]

Downloads go to ``<source>/<server>/<path>`` and processed files to
``<processed>/<name>`` (both under ROTIFER's data directory by default).

Download methods
    ``ftp``: mirror an FTP tree. Only paths matching an ``include``
    regular expression, and none of the ``exclude`` ones, of
    ``<name>.yml`` are downloaded; files are fetched again when their
    size or date changed. ``rsync``: run rsync, with filter rules from
    ``<name>.rules`` when that file exists.

Processors
    ``unpack``: decompress (gzip, bzip2) or extract (tar) every
    downloaded file and link the others into a fresh directory, which
    then replaces the previous one (kept as ``<name>.<date>``). Shell
    scripts listed under ``callback`` in the target's options or
    ``<name>.yml`` run afterwards, inside the new directory, with
    arguments ``<test 0|1> <source> <destination>``.
    ``unpack_pfam``: unpack, then index HMM libraries with hmmpress.
    ``unpack_edirect``: unpack NCBI's EDirect and keep only the programs.

This is a port of the former Perl ``rrsw``. The ``unpack_kegg``,
``unpack_blast``, ``dump_blast`` and ``index_fasta`` processors were not
ported: no target used them, and KEGG no longer offers FTP downloads.
"""

import datetime
import ftplib
import os
import re
import shutil
import socket
import subprocess
import tempfile
import time
from urllib.parse import urlparse
import rotifer

logger = rotifer.logging.getLogger(__name__)

TABLE_COLUMNS = ['NAME', 'METHOD', 'SERVER', 'PATH', 'PROCESS', 'OPTIONS']

DEFAULT_RSYNC_OPTIONS = ['-arHv', '--delete', '--prune-empty-dirs', '--delete-excluded']

# Configuration

def config_dir():
    """
    Directory of rrsw's configuration files (``etc/rotifer/rrsw``).
    """
    return os.path.join(rotifer.config['baseConfig'], 'rrsw')

def default_table():
    """
    The host-specific table of targets, if present, else ``rrsw.yml``.
    """
    host = socket.gethostname().split('.')[0]
    for name in (f'rrsw.{host}.yml', 'rrsw.yml'):
        path = os.path.join(config_dir(), name)
        if os.path.isfile(path):
            return path
    return None

def _yaml_loader():
    import yaml
    class Loader(yaml.SafeLoader):
        pass
    # Former Perl configurations may tag regular expressions
    Loader.add_constructor('tag:yaml.org,2002:perl/regexp', lambda loader, node: loader.construct_scalar(node))
    return Loader

def load_yaml(path):
    """
    Load a YAML configuration file.
    """
    import yaml
    with open(path) as handle:
        return yaml.load(handle, Loader=_yaml_loader())

def load_table(path=None):
    """
    Load the table of targets: a list of rows, the first being the header.
    Missing trailing columns are filled with ``None``.
    """
    rows = load_yaml(path or default_table()) or []
    return [list(row) + [None] * (len(TABLE_COLUMNS) - len(row)) for row in rows]

def select_targets(table, names=(), exclude=()):
    """
    Rows of the table for the named targets (case-insensitive; ``all``
    selects every target), without excluded targets. The header is kept.
    """
    names = {name.upper() for name in names}
    exclude = {name.upper() for name in exclude}
    rows = [row for row in table if str(row[0]).upper() not in exclude]
    if names and 'ALL' not in names:
        rows = [row for row in rows if row[0] == 'NAME' or str(row[0]).upper() in names]
    return rows

def target_options(target, plugin, options=None, defaults=None):
    """
    Options of a plugin for a target: defaults, then the target's table
    options, then ``<plugin>.yml`` and ``<target>.yml`` (later wins).
    """
    merged = dict(defaults or {})
    merged.update(options or {})
    for name in (plugin, target):
        path = os.path.join(config_dir(), f'{name}.yml')
        if os.path.isfile(path):
            merged.update(load_yaml(path) or {})
    return merged

# Download methods

class Rsync:
    """
    Mirror a target with rsync.
    """
    def __init__(self, target, test=False, options=None):
        self.target, self.test = target, test
        options = options or {}
        self.arguments = options.get('rsync', DEFAULT_RSYNC_OPTIONS) if options else DEFAULT_RSYNC_OPTIONS

    def execute(self, origin, destination):
        command = [shutil.which('rsync') or 'rsync', *self.arguments, '-n' if self.test else '--no-motd']
        rules = os.path.join(config_dir(), f'{self.target}.rules')
        if os.path.isfile(rules):
            command += ['-f', f'. {rules}']
        command += [origin, destination]
        return subprocess.run(command).returncode == 0

_LISTING = re.compile(r'^([-dlcbps])[-rwxsStT]{9}\S*\s+\d+\s+\S+\s+\S+\s+(\d+)\s+'
                      r'(\w{3}\s+\d+\s+(?:\d{1,2}:\d{2}|\d{4}))\s+(.+)$')

def parse_listing(lines):
    """
    Parse a Unix-style FTP directory listing into ``(name, type)`` pairs,
    where type is ``d``, ``f`` or ``l <link target>``.
    """
    entries = []
    for line in lines:
        match = _LISTING.match(line)
        if not match:
            continue
        kind, name = match.group(1), match.group(4)
        if kind == 'l':
            name, _, linked = name.partition(' -> ')
            entries.append((name, f'l {linked}'))
        elif kind == 'd':
            entries.append((name, 'd'))
        elif kind == '-':
            entries.append((name, 'f'))
        else:
            entries.append((name, kind))
    return entries

class Ftp:
    """
    Mirror part of an FTP tree, selected by regular expressions.

    Options (from the target's table row and ``<target>.yml``):
    ``include`` and ``exclude`` lists of regular expressions matched
    against remote paths; ``config``: another YAML file to read them from.
    """
    def __init__(self, target, test=False, options=None):
        self.target, self.test = target, test
        options = {'ftp_Passive': 1, 'follow_symlinks': 1, **(options or {})}
        path = os.path.join(config_dir(), f'{target}.yml')
        if options.get('config') and os.path.isfile(options['config']):
            path = options['config']
        if not os.path.isfile(path):
            path = os.path.join(config_dir(), 'ftp.yml')
        if os.path.isfile(path):
            options.update(load_yaml(path) or {})
        self.include = [re.compile(pattern) for pattern in options.get('include') or []]
        self.exclude = [re.compile(pattern) for pattern in options.get('exclude') or []]
        self.passive = bool(options.get('ftp_Passive', 1))

    def _ignore(self, path):
        if not any(pattern.search(path) for pattern in self.include):
            return True
        return any(pattern.search(path) for pattern in self.exclude)

    def execute(self, origin, destination):
        """
        Mirror ``origin`` (``ftp://host/path/``) under ``destination/host``.
        Returns True if any file was (or, in test mode, would be) downloaded.
        """
        uri = urlparse(origin)
        destination = os.path.realpath(destination) if os.path.isdir(destination) else destination
        ftp = ftplib.FTP(uri.hostname, timeout=300)
        ftp.set_pasv(self.passive)
        ftp.login('anonymous', '-@anonymous')
        ftp.voidcmd('TYPE I')
        remoteroot = uri.path if uri.path == '/' else uri.path.strip('/')
        try:
            ftp.cwd(remoteroot)
        except ftplib.error_perm:
            logger.error(f"failed to cd to {remoteroot} on FTP server {uri.hostname}")
            return True
        remotepwd = ftp.pwd()
        if remotepwd != '/':
            remotepwd = remotepwd.lstrip('/')
        base = os.path.join(destination, uri.hostname)
        localroot = os.path.join(base, remoteroot)
        if remotepwd not in ('/', remoteroot):
            # The remote root is a symbolic link: mirror the real directory and link to it
            realroot = os.path.join(base, remotepwd)
            if not self.test:
                if os.path.islink(localroot):
                    linked = os.path.realpath(localroot)
                    if linked.startswith(f'{base}/'):
                        if linked[len(base) + 1:].lstrip('/') != remotepwd:
                            if os.path.exists(linked):
                                shutil.rmtree(linked)
                            os.unlink(localroot)
                    else:
                        os.unlink(localroot)
                elif os.path.exists(localroot):
                    shutil.rmtree(localroot)
                if not os.path.islink(localroot):
                    os.makedirs(os.path.dirname(localroot), exist_ok=True)
                    os.symlink(realroot, localroot)
            localroot = realroot
        if self.test:
            if not os.path.isdir(localroot):
                logger.info(f"creating local root directory {localroot}")
        else:
            os.makedirs(localroot, exist_ok=True)
        try:
            return self._process(ftp, base)
        finally:
            ftp.close()

    def _process(self, ftp, base, check=False):
        remotepwd = ftp.pwd()
        if check and self._ignore(remotepwd):
            return False
        localpwd = re.sub('//+', '/', base + remotepwd)
        status = False
        lines = []
        ftp.retrlines('LIST', lines.append)
        for name, kind in parse_listing(lines):
            if name in ('.', '..'):
                continue
            remotepath = re.sub('//+', '/', f'{remotepwd}/{name}')
            if self._ignore(remotepath):
                continue
            uri = f'ftp://{ftp.host}{remotepath}'
            local = os.path.join(localpwd, name)
            if kind == 'd':
                try:
                    ftp.cwd(name)
                except ftplib.error_perm:
                    continue
                changed = self._process(ftp, base, True)
                if changed and remotepath != ftp.pwd():
                    self._link(re.sub('//+', '/', f'{base}/{ftp.pwd()}'), local)
                ftp.cwd(remotepwd)
                status = status or changed
            elif kind.startswith('l '):
                linked = kind[2:]
                try:
                    contents = [entry for entry in ftp.nlst(name) if entry not in ('.', '..')]
                except ftplib.error_perm:
                    contents = []
                changed = True
                if contents and '/' in contents[0] and contents[0].split('/')[0] == name:
                    try:
                        ftp.cwd(linked)
                    except ftplib.error_perm:
                        continue
                    changed = self._process(ftp, base, True)
                    ftp.cwd(remotepwd)
                if changed:
                    self._link(re.sub('//+', '/', f'{base}/{linked}'), local)
                    status = status or (changed and not self.test)
            elif kind == 'f':
                state = 0
                if os.path.lexists(local):
                    if os.path.getsize(local) != ftp.size(name):
                        state += 1
                    remote_time = datetime.datetime.strptime(ftp.voidcmd(f'MDTM {name}').split()[1][:14],
                                                             '%Y%m%d%H%M%S')
                    if os.path.getmtime(local) < remote_time.replace(tzinfo=datetime.timezone.utc).timestamp():
                        state += 2
                    if not state:
                        continue
                if self.test:
                    logger.info(f"get ({state}) {uri}")
                    status = True
                else:
                    logger.info(f"retrieving ({state}) {uri}")
                    os.makedirs(localpwd, exist_ok=True)
                    with open(local, 'wb') as out:
                        ftp.retrbinary(f'RETR {name}', out.write)
                    status = True
            else:
                logger.warning(f"ignoring unknown file type: {remotepath} ({kind})")
        return status

    def _link(self, origin, local):
        if self.test:
            return
        os.makedirs(os.path.dirname(local), exist_ok=True)
        if os.path.islink(local) or os.path.isfile(local):
            os.unlink(local)
        elif os.path.isdir(local):
            shutil.rmtree(local)
        os.symlink(origin, local)

# Processors

class Unpack:
    """
    Unpack a downloaded tree into a fresh copy of the processed directory.
    """
    SUFFIXES = {'.tar.gz': 'tar', '.tgz': 'tar', '.tar': 'tar', '.gz': 'gzip', '.zip': 'unzip', '.bz2': 'bzip2'}
    DEFAULTS = {'tar': {'command_line_options': ['-xv']}, 'gzip': {'command_line_options': ['-dc']},
                'bzip2': {'command_line_options': ['-dc']}}
    NAME = 'unpack'

    def __init__(self, target, test=False, options=None):
        self.target, self.test = target, test
        self.options = target_options(target, self.NAME, options, self.DEFAULTS)
        self.suffixes = dict(self.SUFFIXES)

    def _command(self, path):
        for suffix in sorted(self.suffixes, key=len, reverse=True):
            if path.endswith(suffix):
                return self.suffixes[suffix], suffix
        return None, None

    def _callbacks(self):
        scripts = []
        for rule in self.options.get('callback') or []:
            for directory in (os.path.join(rotifer.config['user'], 'lib', 'rrsw'),
                              os.path.join(rotifer.config['package'], 'rrsw')):
                if os.path.isfile(os.path.join(directory, rule)):
                    rule = os.path.join(directory, rule)
                    break
            if not rule.endswith('.sh'):
                raise ValueError(f"{self.target}: callback {rule!r} is not a shell script "
                                 "(Perl callbacks are no longer supported)")
            scripts.append(rule)
        return scripts

    def execute(self, source, destination):
        """
        Unpack every file under ``source`` into a new ``destination``.
        """
        scripts = self._callbacks()
        source = source.rstrip('/') or '/'
        if self.test:
            workdir = destination  # nothing is written in test mode
        else:
            os.makedirs(destination, exist_ok=True)
            workdir = tempfile.mkdtemp(prefix=f'.{self.target}.', dir=os.path.dirname(os.path.abspath(destination)))
        files = [source]
        for root, directories, names in os.walk(source):
            for name in sorted(directories) + sorted(names):
                files.append(os.path.join(root, name))
        for path in files:
            self.unpack_each(source, path, workdir)
        for script in scripts:
            subprocess.run([script, '1' if self.test else '0', source, destination], cwd=workdir)
        if not self.test:
            date = time.strftime('%Y%m%d_%H%M')
            if os.path.isdir(f'{destination}.{date}'):
                time.sleep(60)
                date = time.strftime('%Y%m%d_%H%M')
            os.chmod(workdir, 0o775)
            os.rename(destination, f'{destination}.{date}')
            os.rename(workdir, destination)
        return True

    def unpack_each(self, source, path, workdir):
        if not (os.path.isfile(path) or os.path.islink(path)):
            return
        command, suffix = self._command(path)
        name = os.path.basename(path)
        if suffix and name.endswith(suffix) and name != suffix:
            name = name[:-len(suffix)]
        if command != 'tar':
            relative = os.path.relpath(os.path.dirname(path), source)
            if relative != '.':
                if not self.test:
                    os.makedirs(os.path.join(workdir, relative), exist_ok=True)
                name = os.path.join(relative, name)
        options = (self.options.get(command) or {}).get('command_line_options', []) if command else []
        executable = shutil.which(command) or command if command else None
        if command == 'tar':
            argv = [executable, '-f', path, *options]
            if self.test:
                print(' '.join(argv))
            else:
                subprocess.run(argv, cwd=workdir, stdout=subprocess.DEVNULL, check=True)
        elif command in ('gzip', 'bzip2'):
            argv = [executable, path, *options]
            if self.test:
                print(' '.join(argv) + f' >{name}')
            else:
                with open(os.path.join(workdir, name), 'wb') as out:
                    subprocess.run(argv, stdout=out, check=True)
        elif command:
            # As in the Perl version, other formats (zip) are not unpacked
            logger.warning(f"{self.target}: not unpacking {path} (unsupported format)")
        else:
            if self.test:
                print(f'ln -s {path} {name}')
            elif not os.path.lexists(os.path.join(workdir, name)):
                os.symlink(path, os.path.join(workdir, name))

class UnpackPfam(Unpack):
    """
    Unpack Pfam, then index its HMM libraries with ``hmmpress``.
    """
    NAME = 'unpack_pfam'

    def execute(self, source, destination):
        super().execute(source, destination)
        hmmpress = shutil.which('hmmpress') or 'hmmpress'
        for name in sorted(os.listdir(destination)):
            if name.endswith('.hmm'):
                if self.test:
                    print(f'{hmmpress} -f {name}')
                else:
                    subprocess.run([hmmpress, '-f', name], cwd=destination, check=True)
        return True

class UnpackEdirect(Unpack):
    """
    Unpack NCBI's EDirect and keep only its programs.
    """
    NAME = 'unpack_edirect'

    def execute(self, source, destination):
        super().execute(source, destination)
        if self.test:
            return True
        readme = os.path.join(destination, 'README')
        if os.path.isfile(readme):
            os.unlink(readme)
        edirect = os.path.join(destination, 'edirect')
        if os.path.isdir(edirect):
            for name in os.listdir(edirect):
                path = os.path.join(edirect, name)
                if os.path.isfile(path) and (name in ('Mozilla-CA.tar.gz', 'README', 'setup.sh', 'setup-deps.pl')
                                             or name.endswith('.go')):
                    os.unlink(path)
            for name in os.listdir(edirect):
                shutil.move(os.path.join(edirect, name), os.path.join(destination, name))
            os.rmdir(edirect)
        for name in os.listdir(destination):
            if name.endswith('.Linux'):
                os.chmod(os.path.join(destination, name), 0o755)
        return True

METHODS = {'ftp': Ftp, 'rsync': Rsync}
PROCESSORS = {'unpack': Unpack, 'unpack_pfam': UnpackPfam, 'unpack_edirect': UnpackEdirect}

def _enabled(name):
    return name is not None and str(name) != '' and str(name).lower() not in ('off', 'none')

def rrsw(targets, table=None, source=None, processed=None, download=True, process=True,
         exclude=(), test=False):
    """
    Download and process targets.

    Parameters
    ----------
    targets : list of str
        Target names (``all`` for every target in the table).
    table : str, optional
        Table of targets (default: see ``default_table``).
    source : str, optional
        Root of downloaded files (default: ``<data>/src``).
    processed : str, optional
        Root of processed files (default: ROTIFER's data directory).
    download, process : bool
        Whether to download and to process the targets.
    exclude : list of str
        Targets to skip.
    test : bool
        Only show what would be done.
    """
    datadir = rotifer.config['data']
    source = source or os.path.join(datadir, 'src')
    processed = processed or datadir
    for row in select_targets(load_table(table), targets, exclude):
        name, method, server, path, processor, options = row[:6]
        if name == 'NAME':
            continue
        options = options or {}
        origin = f'{method}://{server}/{path}/'
        destination = os.path.join(source, str(server), str(path), '')
        if not test:
            os.makedirs(destination, exist_ok=True)
        changed = not download
        if download and _enabled(method):
            if method not in METHODS:
                raise ValueError(f"{name}: unknown download method {method!r}")
            downloader = METHODS[method](name, test=test, options=options)
            changed = downloader.execute(origin, source if method == 'ftp' else destination)
        if changed and process and _enabled(processor):
            if processor not in PROCESSORS:
                raise ValueError(f"{name}: unknown processor {processor!r}")
            PROCESSORS[processor](name, test=test, options=options).execute(
                destination, os.path.join(processed, name))
