"""
Cursors that delegate data access to other cursors.

A delegator cursor does not talk to any database itself. Instead, it
instantiates one backend cursor per data source and forwards queries
to them, so that identifiers missing from one source can be recovered
from the next. The mapping between backend names (such as ``entrez``
or ``mirror``) and the modules that implement them is read from the
``readers`` and ``writers`` entries of the calling module's ``config``
dictionary (see, for example, :mod:`rotifer.db.ncbi`).
"""

# External libraries
import types
from copy import deepcopy

# Rotifer
import rotifer.db.core
import rotifer
logger = rotifer.logging.getLogger(__name__)

class DelegatorCursor(rotifer.db.core.BaseCursor):
    """
    Base class for cursors that dispatch queries to backend cursors.

    Subclasses are expected to live in a module that defines a
    ``config`` dictionary with ``readers`` and ``writers`` entries
    mapping backend names to module paths. For each requested backend,
    the delegator imports the module and instantiates the class named
    like the delegator subclass itself.

    Parameters
    ----------
    readers : list of str, default []
        Names of the backend modules used to retrieve data, in order
        of preference.
    writers : list of str, default []
        Names of the backend modules used to store retrieved data.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, optional
        Number of attempts to download each batch, shared with every
        backend cursor that accepts it.
    batch_size : int, optional
        Number of accessions per batch, shared with every backend
        cursor that accepts it.
    threads : int, optional
        Number of simultaneous threads, shared with every backend
        cursor that accepts it.

    Attributes
    ----------
    cursors : dict
        Backend cursor instances, keyed by backend name.

    Notes
    -----
    Attributes listed in the subclass's ``_shared_attributes`` are
    propagated to every backend cursor whenever they are set on the
    delegator.
    """

    def __init__(self, readers=[], writers=[], progress=True, tries=None, batch_size=None, threads=None, *args, **kwargs):
        super().__init__(progress=progress, *args, **kwargs)
        self.readers = readers.copy()
        self.writers = writers.copy()
        self.tries = tries
        self.batch_size = batch_size
        self.threads = threads
        self.reset_cursors()

    @property
    def _cursor_modules(self):
        """
        Reader and writer modules named in the configuration.

        Returns
        -------
        dict
            Imported modules, keyed by the names listed under the
            ``readers`` and ``writers`` keys of the ``config``
            dictionary of the module that defines this cursor.

        Raises
        ------
        ValueError
            If that module has no ``config`` attribute, or its
            configuration lacks the ``readers`` or ``writers`` key.
        """
        import inspect
        import importlib
        mymodule = inspect.getmodule(self)
        cursor_modules = dict()

        # Check configuration
        try:
            myconfig = getattr(mymodule,'config')
        except:
            error = f'Module {mymodule.__name__} has no configuration! Blame the developer!'
            logger.error(error)
            raise ValueError(f'No attribute "config" in module {mymodule.__name__}')
        if 'readers' not in myconfig:
            error = f'Missing dictionary of reader modules in configuration of module {mymodule.__name__}'
            logger.error(error)
            raise ValueError(error)
        if 'writers' not in myconfig:
            error = f'Missing dictionary of writer modules in configuration of module {mymodule.__name__}'
            logger.error(error)
            raise ValueError(error)

        # Load modules
        for module in set(self.readers + self.writers):
            if "readers" in myconfig and module in myconfig['readers']:
                module_name = myconfig['readers'][module]
            elif "writers" in myconfig and module in myconfig['writers']:
                module_name = myconfig['writers'][module]
            else:
                error = f'Missing module name "{module}" for writers or readers in {mymodule.__name__}.config'
                logger.error(error)
                raise ValueError(error)
            try:
                cursor_modules[module] = importlib.import_module(module_name)
            except:
                logger.error(f'Unable to load module {module_name}: %s.', exc_info=1)
                raise ImportError(f'Unable to load module {module_name}')

        return cursor_modules

    def reset_cursors(self):
        """
        Instantiate or reinstantiate every backend cursor.

        The backend class is looked up in each backend module using
        the delegator's own class name. Attributes listed in
        ``_shared_attributes`` that are not ``None`` are passed to
        the backend constructors. Backends whose module does not
        define a matching class are skipped with an error message.
        """
        myname = str(type(self)).split("'")[1].split(".")[-1]
        if hasattr(self,"_shared_attributes"):
            kwargs = { x: getattr(self,x) for x in self._shared_attributes if not isinstance(getattr(self,x),types.NoneType) }
        else:
            kwargs = dict()
        self.cursors = dict()
        cursor_modules = self._cursor_modules
        for modulename in cursor_modules:
            module = cursor_modules[modulename]
            try:
                cursorClass = getattr(module,myname)
            except:
                logger.error(f'Module {module.__name__} does not define a {myname} class')
                continue
            self.cursors[modulename] = cursorClass(**kwargs)

    def __setattr__(self, name, value):
        """
        Set an attribute, propagating shared ones to the backends.

        Attributes named in ``_shared_attributes`` are also assigned
        on every backend cursor that already defines them, so that
        retuning the delegator keeps its backends in sync. ``None``
        values are never propagated.

        Parameters
        ----------
        name : str
            Attribute name.
        value : object
            Value to assign. Forwarded to the backends only when it
            is not None.
        """
        super().__setattr__(name, value)
        if hasattr(self,'cursors') and hasattr(self,'_shared_attributes') and name in self._shared_attributes:
            for cursor in self.cursors.values():
                if hasattr(cursor,name) and not isinstance(value,types.NoneType):
                    cursor.__setattr__(name,value)

class SequentialDelegatorCursor(DelegatorCursor):
    """
    Delegator that queries its backends one after the other.

    Backends listed in ``readers`` are tried in order and each one
    only receives the identifiers that previous backends could not
    resolve. Results retrieved by a reader are handed to every
    backend listed in ``writers`` for storage.

    Parameters
    ----------
    readers : list of str, default []
        Names of the backend modules used to retrieve data, in order
        of preference.
    writers : list of str, default []
        Names of the backend modules used to store retrieved data.
    progress : bool, default True
        Whether to print a progress bar.
    tries : int, optional
        Number of attempts to download each batch.
    batch_size : int, optional
        Number of accessions per batch.
    threads : int, optional
        Number of simultaneous threads.

    See Also
    --------
    rotifer.db.ncbi.SequenceCursor : a concrete sequential delegator
    """

    def __init__(self, readers=[], writers=[], progress=True, tries=None, batch_size=None, threads=None, *args, **kwargs):
        super().__init__(readers=readers, writers=writers, progress=progress, tries=tries, batch_size=batch_size, threads=threads, *args, **kwargs)

    def __getitem__(self, accessions, *args, **kwargs):
        """
        Fetch data for one or more entries, dictionary style.

        Parameters
        ----------
        accessions : str or list of str
            Database identifiers.

        Returns
        -------
        object
            The value returned by the delegated cursors: a single
            result object when a single accession yields a single
            result, otherwise a list of result objects.

        Examples
        --------
        >>> import rotifer.db.ncbi as ncbi
        >>> tc = ncbi.TaxonomyCursor(progress=False)
        >>> t = tc[2599]  # doctest: +SKIP
        """
        targets = self.parse_ids(accessions)

        # Call cursors
        data = []
        todo = deepcopy(targets)
        for i in range(0,len(self.readers)):
            if len(todo) == 0:
                break
            cursorName = self.readers[i]
            if cursorName in self.cursors:
                cursor = self.cursors[cursorName]
            else:
                continue
            result = cursor.__getitem__(todo, *args, **kwargs)
            found = self.getids(result, *args, **kwargs)
            done = targets.intersection(found)
            for j in range(0,i+1):
                c = self.readers[j]
                if c not in self.cursors:
                    continue
                self.cursors[c].remove_missing(done)
            self.remove_missing(done)
            self.update_missing(data=cursor._missing)
            if isinstance(result,types.NoneType):
                continue
            if isinstance(result, list):
                data.extend(result)
            else:
                data.append(result)
            todo = todo - done
        if len(targets) == 1 and len(data) == 1:
            data = data[0]
        return data

    def fetchone(self, accessions, *args, **kwargs):
        """
        Get a generator to fetch data for iteratively.

        Parameters
        ----------
        accessions: list of strings
          Database identifiers.

        Returns
        -------
        Generator of Pandas dataframes
        """
        targets = self.parse_ids(accessions)

        # Call cursors
        todo = deepcopy(targets)
        for i in range(0,len(self.readers)):
            if len(todo) == 0:
                break
            cursorName = self.readers[i]
            if cursorName in self.cursors:
                cursor = self.cursors[cursorName]
            else:
                continue
            for result in cursor.fetchone(todo, *args, **kwargs):
                found = self.getids(result, *args, **kwargs)
                done = todo.intersection(found)
                for j in range(0,i+1):
                    c = self.readers[j]
                    if c not in self.cursors:
                        continue
                    self.cursors[c].remove_missing(done)
                for j in self.writers:
                    if j not in self.cursors:
                        continue
                    self.cursors[j].insert(result)
                self.remove_missing(done)
                self.update_missing(data=cursor._missing)
                todo = todo - done
                yield result

    def fetchall(self, accessions, *args, **kwargs):
        """
        Fetch data for all accessions.

        Parameters
        ----------
        accessions: list of database identifiers
          Database identifiers.

        Returns
        -------
        Pandas dataframe
        """
        targets = self.parse_ids(accessions)
        stack = []
        for data in self.fetchone(targets, *args, **kwargs):
            if isinstance(data,list):
                stack.extend(data)
            else:
                stack.append(data)
        return stack

