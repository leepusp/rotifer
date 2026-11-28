"""Manifest v1: the JSON file that describes a saved session.

The manifest is `<session dir>/manifest.json`. It holds plain JSON types only,
so other tools can read it without importing ROTIFER. `manifest_schema()` is the
single definition of its structure.
"""

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

from rotifer.core.session._files import atomic_write_text

FORMAT = "rotifer-session"
VERSION = 1
MANIFEST_NAME = "manifest.json"

# Maps an old manifest version N to a function that turns a version-N dict
# into a version N+1 dict. Empty while v1 is the only version.
_UPGRADES = {}


class ManifestError(ValueError):
    """
    The manifest is not valid or has an unsupported version.

    See Also
    --------
    Manifest.from_dict : Raises it for an invalid manifest.
    """


@dataclass
class VariableEntry:
    """
    What the manifest records about one saved variable.

    Parameters
    ----------
    type : str
        Name of the variable's type.
    codec : str
        Name of the codec that wrote the files.
    codec_version : int
        Version of that codec.
    files : list, optional
        File names, relative to the session directory. Default is empty.
    bytes : int, default 0
        Total size of the files in bytes.
    sha256 : str, optional
        SHA-256 digest of the files. Default is an empty string.
    fingerprint : str or None, default None
        Fingerprint of the value, or None if it could not be computed.
    saved : str, optional
        Time the variable was saved. Default is an empty string.
    version_bound : bool, default False
        True if the files can only be read by the same Python version.
    summary : dict, optional
        Summary of the value, as returned by `describe`. Default is empty.
    meta : dict, optional
        Codec-specific metadata, from the codec's write function. Default
        is empty.

    See Also
    --------
    Manifest : The whole manifest of a session.
    """

    type: str
    codec: str
    codec_version: int
    files: list = field(default_factory=list)  # relative to the session dir
    bytes: int = 0
    sha256: str = ""
    fingerprint: str | None = None
    saved: str = ""
    version_bound: bool = False
    summary: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)  # codec-specific, from its write()

    def to_dict(self):
        """
        Convert the entry to a dict of plain JSON types.

        Returns
        -------
        dict
            One key per field.

        See Also
        --------
        from_dict : The inverse conversion.
        """
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        """
        Build an entry from a dict, such as one made by `to_dict`.

        Parameters
        ----------
        data : dict
            Field values by name. Keys that are not fields are ignored.

        Returns
        -------
        VariableEntry
            The entry.

        See Also
        --------
        to_dict : The inverse conversion.
        """
        names = cls.__dataclass_fields__
        return cls(**{k: v for k, v in data.items() if k in names})


@dataclass
class Manifest:
    """
    The whole manifest of a session.

    Parameters
    ----------
    revision : int, default 0
        Number that grows with each save.
    created : str, optional
        Time the session was created. Default is an empty string.
    updated : str, optional
        Time of the last save. Default is an empty string.
    environment : dict, optional
        Versions of Python, ROTIFER and the main packages. Default has
        "python" and "rotifer" set to None and empty dicts for "packages"
        and "tools".
    context : dict, optional
        Notebook context. Default is empty.
    variables : dict, optional
        Maps variable names to `VariableEntry` objects. Default is empty.
    failed : dict, optional
        Maps names of variables that could not be saved to an error string.
        Default is empty.
    format : str, default "rotifer-session"
        Format identifier.
    version : int, default 1
        Manifest version.

    See Also
    --------
    VariableEntry : What the manifest records about one variable.
    read_manifest : Read a manifest from a session directory.
    write_manifest : Write a manifest to a session directory.
    manifest_schema : JSON Schema of the manifest.
    """

    revision: int = 0
    created: str = ""
    updated: str = ""
    environment: dict = field(default_factory=lambda: {
        "python": None, "rotifer": None, "packages": {}, "tools": {}})
    context: dict = field(default_factory=dict)
    variables: dict = field(default_factory=dict)  # name -> VariableEntry
    failed: dict = field(default_factory=dict)  # name -> error string
    format: str = FORMAT
    version: int = VERSION

    def to_dict(self):
        """
        Convert the manifest to a dict of plain JSON types.

        Returns
        -------
        dict
            The manifest, with each `VariableEntry` converted to a dict.

        See Also
        --------
        from_dict : The inverse conversion.
        """
        return {
            "format": self.format,
            "version": self.version,
            "revision": self.revision,
            "created": self.created,
            "updated": self.updated,
            "environment": self.environment,
            "context": self.context,
            "variables": {k: v.to_dict() for k, v in self.variables.items()},
            "failed": dict(self.failed),
        }

    @classmethod
    def from_dict(cls, data):
        """
        Build a manifest from a dict, such as one made by `to_dict`.

        Older manifest versions are upgraded step by step to the current one.

        Parameters
        ----------
        data : dict
            The manifest as plain JSON types.

        Returns
        -------
        Manifest
            The manifest.

        Raises
        ------
        ManifestError
            If "format" is not "rotifer-session", "version" is missing or
            not an integer, "version" is newer than the supported version,
            or there is no upgrade from an older version.

        See Also
        --------
        to_dict : The inverse conversion.
        """
        data = dict(data)
        if data.get("format") != FORMAT:
            raise ManifestError(f"not a {FORMAT} manifest (format={data.get('format')!r})")
        version = data.get("version")
        if not isinstance(version, int):
            raise ManifestError(f"manifest version is missing or not an integer: {version!r}")
        if version > VERSION:
            raise ManifestError(
                f"manifest version {version} is newer than the supported "
                f"version {VERSION}; upgrade ROTIFER to read it")
        while version < VERSION:
            if version not in _UPGRADES:
                raise ManifestError(f"no upgrade from manifest version {version}")
            data = _UPGRADES[version](data)
            version += 1
            data["version"] = version
        return cls(
            revision=data.get("revision", 0),
            created=data.get("created", ""),
            updated=data.get("updated", ""),
            environment=data.get("environment") or {
                "python": None, "rotifer": None, "packages": {}, "tools": {}},
            context=data.get("context") or {},
            variables={k: VariableEntry.from_dict(v)
                       for k, v in (data.get("variables") or {}).items()},
            failed=dict(data.get("failed") or {}),
            format=data["format"],
            version=version,
        )


def read_manifest(directory):
    """
    Read `<directory>/manifest.json`. Never locks.

    Parameters
    ----------
    directory : str or pathlib.Path
        The session directory.

    Returns
    -------
    Manifest
        The manifest.

    Raises
    ------
    FileNotFoundError
        If the manifest file does not exist.
    json.JSONDecodeError
        If the file is not valid JSON.
    ManifestError
        If the content is not a valid manifest or has an unsupported
        version.

    See Also
    --------
    write_manifest : Write the manifest atomically.
    """
    path = Path(directory) / MANIFEST_NAME
    with open(path, encoding="utf-8") as handle:
        return Manifest.from_dict(json.load(handle))


def write_manifest(directory, manifest):
    """
    Write the manifest atomically.

    The file is `<directory>/manifest.json`. A reader sees either the old or
    the new file, never a partial one.

    Parameters
    ----------
    directory : str or pathlib.Path
        The session directory. It is created if it does not exist.
    manifest : Manifest
        The manifest to write.

    See Also
    --------
    read_manifest : Read a manifest from a session directory.
    """
    text = json.dumps(manifest.to_dict(), indent=2, sort_keys=False) + "\n"
    atomic_write_text(Path(directory) / MANIFEST_NAME, text)


def manifest_schema():
    """
    Return the JSON Schema (draft 2020-12) of manifest v1 as a dict.

    Returns
    -------
    dict
        The schema.

    See Also
    --------
    read_manifest : Read a manifest from a session directory.

    Examples
    --------
    >>> schema = manifest_schema()
    >>> schema["title"]
    'ROTIFER session manifest'
    >>> schema["required"][:3]
    ['format', 'version', 'revision']
    """
    nullable_str = {"type": ["string", "null"]}
    summary = {
        "type": "object",
        "properties": {
            "kind": {"type": "string"},
            "length": {"type": ["integer", "null"]},
            "shape": {"type": ["array", "null"], "items": {"type": "integer"}},
            "columns": {"type": ["array", "null"], "items": {"type": "string"}},
            "dtypes": {"type": ["object", "null"],
                       "additionalProperties": {"type": "string"}},
            "preview": {"type": "string"},
        },
    }
    variable = {
        "type": "object",
        "required": ["type", "codec", "codec_version", "files", "bytes",
                     "sha256", "fingerprint", "saved", "version_bound",
                     "summary", "meta"],
        "properties": {
            "type": {"type": "string"},
            "codec": {"type": "string"},
            "codec_version": {"type": "integer"},
            "files": {"type": "array", "items": {"type": "string"}},
            "bytes": {"type": "integer", "minimum": 0},
            "sha256": {"type": "string"},
            "fingerprint": nullable_str,
            "saved": {"type": "string"},
            "version_bound": {"type": "boolean"},
            "summary": summary,
            "meta": {"type": "object"},
        },
    }
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": "ROTIFER session manifest",
        "type": "object",
        "required": ["format", "version", "revision", "created", "updated",
                     "environment", "context", "variables", "failed"],
        "properties": {
            "format": {"const": FORMAT},
            "version": {"const": VERSION},
            "revision": {"type": "integer", "minimum": 0},
            "created": {"type": "string"},
            "updated": {"type": "string"},
            "environment": {
                "type": "object",
                "required": ["python", "rotifer", "packages", "tools"],
                "properties": {
                    "python": nullable_str,
                    "rotifer": nullable_str,
                    "packages": {"type": "object",
                                 "additionalProperties": nullable_str},
                    "tools": {"type": "object",
                              "additionalProperties": nullable_str},
                },
            },
            "context": {"type": "object"},
            "variables": {"type": "object", "additionalProperties": variable},
            "failed": {"type": "object",
                       "additionalProperties": {"type": "string"}},
        },
    }
