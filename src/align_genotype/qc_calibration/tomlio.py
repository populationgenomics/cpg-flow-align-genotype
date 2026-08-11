"""TOML reading and scalar writing for the calibration artifacts.

Reading uses ``tomllib`` (stdlib from 3.11, ``tomli`` on 3.10 - already in the resolved
dependency graph under a version marker, so nothing new is declared). Writing is
hand-rolled: manifests and specs are flat tables of scalars, which doesn't justify a
style-preserving TOML writer dependency.
"""

import re
import sys
from pathlib import Path
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - 3.10 only
    import tomli as tomllib

# Keys we can write unquoted. Metric keys come from Picard/samtools and cohort labels
# from Metamist dataset names, so both are already identifier-like; anything else is
# rejected at load time rather than silently mangled when we write the file back out.
BARE_KEY_RE = re.compile(r'^[A-Za-z0-9_-]+$')


def require_bare_key(key: str, kind: str) -> None:
    """Raise unless `key` can be written as an unquoted TOML key."""
    if not BARE_KEY_RE.match(key):
        raise ValueError(f'{kind} {key!r} is not a bare TOML key; expected only letters, digits, underscore or hyphen')


def loads(text: str) -> dict[str, Any]:
    """Parse a TOML document from a string."""
    return tomllib.loads(text)


def load_path(path: str | Path) -> dict[str, Any]:
    """Parse a TOML document from a local or cloud path."""
    from cpg_utils import to_path  # noqa: PLC0415

    with to_path(path).open('rb') as f:
        return tomllib.load(f)


def fmt_value(value: Any) -> str:
    """Render a Python scalar as a TOML value."""
    if isinstance(value, bool):  # bool subclasses int - must be checked first
        return 'true' if value else 'false'
    if isinstance(value, (int, float)):
        return repr(value)
    escaped = str(value).replace('\\', '\\\\').replace('"', '\\"')
    return f'"{escaped}"'


def fmt_kv(key: str, value: Any, quote_key: bool = False) -> str:
    """Render a ``key = value`` line. `quote_key` matches config_template.toml's style."""
    rendered_key = f'"{key}"' if quote_key else key
    return f'{rendered_key} = {fmt_value(value)}'
