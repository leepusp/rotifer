"""
Query the local NCBI Taxonomy database managed by the ETE toolkit.

The ETE toolkit (ete3) keeps a SQLite copy of the NCBI Taxonomy
database under the user's home directory. This module exposes that
database through the rotifer cursor interface. The first use of ete3
downloads the taxonomy dump from NCBI; afterwards all queries are
local and no network connection is used. Call
:meth:`TaxonomyCursor.update_database` to refresh the local copy.
"""

# Import external modules
import os
import sys
import types
import socket
import typing
import numpy as np
import pandas as pd
from copy import deepcopy
from datetime import datetime, timedelta

# Ete3
from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa

# Rotifer
import rotifer
import rotifer.db.core
from rotifer import GlobalConfig
from rotifer.taxonomy.utils import lineage
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)
config = loadConfig(__name__, defaults = {
    "taxdump_file": None,
})

# Classes
class TaxonomyCursor(rotifer.db.core.BaseCursor):
    """
    Fetch NCBI Taxonomy data from the local ETE toolkit database.

    Phased out taxonomy identifiers are transparently replaced by
    their current equivalents, and the original identifiers are
    reported in the ``alternative_taxids`` column.

    Parameters
    ----------
    progress : bool, default False
        Whether to print a progress bar.

    Attributes
    ----------
    ete3 : ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
        The underlying ETE toolkit database handle.
    taxcols : list of str
        Columns of the returned dataframes: ``taxid``, ``organism``,
        ``superkingdom``, ``lineage``, ``classification`` and
        ``alternative_taxids``.

    See Also
    --------
    rotifer.db.ncbi.TaxonomyCursor : delegator that falls back to
        NCBI Entrez

    Examples
    --------
    >>> from rotifer.db.local import ete3
    >>> tc = ete3.TaxonomyCursor()  # doctest: +SKIP
    >>> t = tc.fetchall([2599])  # doctest: +SKIP
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ete3 = NCBITaxa()
        self.taxcols = ['taxid','organism','superkingdom','lineage','classification','alternative_taxids']

    def update_database(self):
        """
        Download the latest NCBI taxdump and rebuild the local
        database.

        The download is performed by the ETE toolkit in a temporary
        directory. This is the only method of this class that uses
        the network.
        """
        from tempfile import TemporaryDirectory
        cwd = os.getcwd()
        with TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            self.ete3.update_taxonomy_database()
            os.chdir(cwd)

    def getids(self,obj):
        """
        Extract taxonomy identifiers from taxonomy dataframes.

        Parameters
        ----------
        obj : pandas.DataFrame or list of pandas.DataFrame
            Taxonomy dataframes produced by the cursor.

        Returns
        -------
        set of str
            All identifiers in the ``taxid`` and
            ``alternative_taxids`` columns.
        """
        if not isinstance(obj,list):
            obj = [obj]
        ids = set()
        for o in obj:
            if "taxid" in o.columns:
                ids.update(o.taxid.astype(str))
            if "alternative_taxids" in o.columns:
                ids.update(o.alternative_taxids.astype(str))
        return ids

    def __getitem__(self, accessions):
        """
        Fetch taxonomy data for one or more taxids, dictionary
        style.

        Parameters
        ----------
        accessions : int, str or iterable
            NCBI Taxonomy identifiers.

        Returns
        -------
        pandas.DataFrame
            One row per taxon, with the columns listed in
            ``taxcols``. Identifiers not found in the database
            are registered as missing.
        """
        targets = self.parse_ids(accessions)
        targetsAsIntegers = pd.Series(list(targets)).astype(int)
        #logger.info(f'Loading {len(targets)} taxids from Ete3 local database')

        # Find replacements for phased-out taxids
        oldToNewDF = ",".join(targets)
        oldToNewDF = f'SELECT * from merged where taxid_old in ({oldToNewDF})'
        oldToNewDF = pd.read_sql(oldToNewDF, self.ete3.db)
        idmap = oldToNewDF.set_index('taxid_old').taxid_new.to_dict()
        targetsAsIntegers.replace(idmap, inplace=True)

        # Fetch data from database
        tl = self.ete3.get_lineage_translator(targetsAsIntegers.tolist())
        if len(tl) == 0:
            return pd.DataFrame(columns=self.taxcols)
        ti = pd.Series(tl.values()).explode().unique()
        ti = self.ete3.get_taxid_translator(ti)
        li = [[x,ti[x],np.nan,np.nan,"; ".join([ ti[y] for y in tl[x] if ti[y] not in ["root","cellular organisms"] ]),x] for x in tl ]
        li = pd.DataFrame(li, columns=self.taxcols)
        li['superkingdom'] = li.classification.str.split("; ", expand=True)[0]
        li['lineage'] = lineage(li.classification)

        # Shuffle alternative taxids and query accessions
        idmap = oldToNewDF.set_index('taxid_new').taxid_old.to_dict()
        li.alternative_taxids = li.alternative_taxids.replace(idmap)
        li = li.astype(str)

        # Register missing entries
        missing = targets - self.getids(li)
        if len(missing) > 0:
            self.update_missing(missing,"Accession not found in database")

        #logger.info(f'Loaded {len(targets.intersection(self.getids(li)))} taxids from Ete3 database')
        return li

    def fetchone(self,accessions):
        """
        Iterate over taxonomy records, one taxon at a time.

        Parameters
        ----------
        accessions : iterable
            NCBI Taxonomy identifiers.

        Yields
        ------
        pandas.DataFrame
            The rows of one taxon.
        """
        # Process accessions
        li = self.__getitem__(accessions)
        for taxid, rows in li.groupby('taxid'):
            yield rows.copy()

    def fetchall(self,accessions):
        """
        Fetch taxonomy data for all taxids as one dataframe.

        Parameters
        ----------
        accessions : iterable
            NCBI Taxonomy identifiers.

        Returns
        -------
        pandas.DataFrame
            One row per taxon, with the columns listed in
            ``taxcols``.
        """
        return self.__getitem__(accessions)
