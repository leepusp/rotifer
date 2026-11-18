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
    """The manifest is not valid or has an unsupported version."""


@dataclass
class VariableEntry:
    """What the manifest records about one saved variable."""

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
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        names = cls.__dataclass_fields__
        return cls(**{k: v for k, v in data.items() if k in names})


@dataclass
class Manifest:
    """The whole manifest of a session."""

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
    """Read `<directory>/manifest.json`. Never locks."""
    path = Path(directory) / MANIFEST_NAME
    with open(path, encoding="utf-8") as handle:
        return Manifest.from_dict(json.load(handle))


def write_manifest(directory, manifest):
    """Write the manifest atomically."""
    text = json.dumps(manifest.to_dict(), indent=2, sort_keys=False) + "\n"
    atomic_write_text(Path(directory) / MANIFEST_NAME, text)


def manifest_schema():
    """Return the JSON Schema (draft 2020-12) of manifest v1 as a dict."""
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
