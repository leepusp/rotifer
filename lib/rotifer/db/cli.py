
#!/usr/bin/env python3

"""
Reusable command line argument sets for database tools.

This module defines argument parser builders shared by rotifer's
command line programs that query the laboratory's PostgreSQL genome
annotation database: :class:`sql` collects the connection options and
:class:`rneighbors` collects the options of the ``rneighbors``
program.
"""

__version__ = "0.01"

### Import rotifer package
import sys
import os
sys.path.insert(0, os.path.join(os.getcwd(), '../..'))

import rotifer.core.cli as corecli

### Other packages
import argparse

class sql:
    """
    Argument parser builder for SQL database connection options.

    The available options are ``--database``, ``--port``, ``--host``
    and ``--user``, with defaults matching the laboratory's internal
    PostgreSQL server.
    """

    def __init__(self):
         self.dc = {'database': {'short': '-d',
                        'long': '--database',
                        'default': 'rotifer',
                        'help': "Rotifer's database name"},

             'port': {'short': '-p',
                        'long': '--port',
                        'default': '5433',
                        'help': "Rotifer's database server port number"},
             'host': {'short': '-H',
                        'long': '--host',
                        'default': '10.1.1.1',
                        'help': "Rotifer's database server name or IP address"},
             'user': {'short': '-u',
                        'long': '--user',
                        'default': 'rotifer',
                        'help': "User name to connect to Rotifer's database"}
                    }
    def input(self, number = 1, exclude = [],
              input_order = ['database', 'port','host','user'],
              add_help = False):
        """
        Build an argument parser with the connection options.

        Parameters
        ----------
        number : int, default 1
            Unused, kept for interface compatibility.
        exclude : list of str, default []
            Option names to leave out of the parser.
        input_order : list of str, default ``['database', 'port', 'host', 'user']``
            Options to add, in order.
        add_help : bool, default False
            Whether the parser defines its own ``--help`` option.

        Returns
        -------
        argparse.ArgumentParser
            The configured parser.
        """
        self.parser = argparse.ArgumentParser(add_help = add_help)
        for choice in input_order:
            if choice not in exclude:
                self.parser.add_argument(self.dc[choice]['short'],
                                         self.dc[choice]['long'],
                                         default = self.dc[choice]['default'],
                                         help = self.dc[choice]['help'],
                                         type = str)
        return self.parser


class rneighbors:
    """
    Argument parser builder for the ``rneighbors`` program.

    The options cover the input file of protein accessions, the size
    of the neighborhood (``--above`` and ``--below``), filters by
    genomic accession or assembly, output format selection and
    logging destinations.
    """

    def __init__(self):
        self.dc = {'above': {'short': '-a',
                        'long': '--above',
                        'default': 3,
                        'help': 'Rows above: maximum number of neighbors upstream of target loci'},

             'below': {'short': '-b',
                        'long': '--below',
                        'default': 3,
                        'help': 'Rows below:  maximum number of neighbors downstream of target loci'},

             'debug': {'short': '-de',
                        'long': '--debug',
                        'default': 0,
                        'help': 'Set verbosity level for error/warning/info messages'},

             'outformat': {'short': '-of',
                        'long': '--outformat',
                        'default': 'gi2operons',
                        'help': '''Output format. Currently available formats are:
            1) gi2operons: visual display of neighborhoods as blocks of rows (default);
            2) table: tabular output;
            3) missing: list of queries not present in the SQL database
                                '''},

             'log': {'short': '',
                        'long': '--log',
                        'default': '',
                        'help': 'Choose destination of error, warning, info and debug messages'},

             'file': {'short': 'file',
                        'long': '',
                        'default': '',
                        'help': "Input file(s)"},

             'gacc': {'short': '-gacc',
                        'long': '--genomic_accession',
                        'default': [],
                        'help': "A list containing genomic_accessions"},

             'asm': {'short': '-asm',
                        'long': '--assembly',
                        'default': [],
                        'help': "A list containing assemblies"}}

    def input(self, number = 1, exclude = [],
              input_order = ['file','above', 'below', 'gacc', 'asm',
                              'debug','outformat', 'log'],
              add_help = False, no_open = False):
        """
        Build an argument parser with the ``rneighbors`` options.

        Parameters
        ----------
        number : int, default 1
            Unused, kept for interface compatibility.
        exclude : list of str, default []
            Option names to leave out of the parser.
        input_order : list of str, optional
            Options to add, in order. Defaults to every available
            option: ``file``, ``above``, ``below``, ``gacc``,
            ``asm``, ``debug``, ``outformat`` and ``log``.
        add_help : bool, default False
            Whether the parser defines its own ``--help`` option.
        no_open : bool, default False
            When True, the ``--genomic_accession`` and ``--assembly``
            options take plain values instead of files loaded through
            :mod:`rotifer.core.cli`.

        Returns
        -------
        argparse.ArgumentParser
            The configured parser.
        """
        self.parser = argparse.ArgumentParser(add_help = add_help)
        for choice in input_order:
            if choice not in exclude:
                if choice == 'file':
                    self.parser.add_argument(self.dc[choice]['short'],
                                             nargs = '*',
                                             action = corecli.action.autoload,
                                             help = self.dc[choice]['help'])
                elif choice == 'log':
                    self.parser.add_argument(
                                             self.dc[choice]['long'],
                                             default = '',
                                             help = self.dc[choice]['help'])
                elif choice == 'debug':
                    self.parser.add_argument(self.dc[choice]['short'],
                                             self.dc[choice]['long'],
                                             default = 0,
                                             type = int,
                                             help = self.dc[choice]['help'])

                elif choice == 'above' or choice =='below':
                    self.parser.add_argument(self.dc[choice]['short'],
                                             self.dc[choice]['long'],
                                             default = self.dc[choice]['default'],
                                             type = int,
                                             help = self.dc[choice]['help'])
                elif choice == 'asm' or choice == 'gacc':
                    if not no_open:
                        self.parser.add_argument(self.dc[choice]['short'],
                                                 self.dc[choice]['long'],
                                                 default = self.dc[choice]['default'],
                                                 help = self.dc[choice]['help'],
                                                 action = corecli.action.autoload)
                    if no_open:
                        self.parser.add_argument(self.dc[choice]['short'],
                                                 self.dc[choice]['long'],
                                                 default = self.dc[choice]['default'],
                                                 help = self.dc[choice]['help']
                                                )
                else:
                    self.parser.add_argument(self.dc[choice]['short'],
                                             self.dc[choice]['long'],
                                             default = self.dc[choice]['default'],
                                             help = self.dc[choice]['help'],
                                             type = str)
        return self.parser



