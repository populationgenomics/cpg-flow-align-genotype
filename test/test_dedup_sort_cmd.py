"""Unit tests for align.dedup_sort_cmd, the dupblaster + coordinate-sort segment of the alignment pipe.

Config access is monkeypatched, so these run offline. The command is inspected with shlex (POSIX
word-splitting, no expansion) and syntax-checked with `bash -n`, which is how the job runs it.
"""

import shlex
import subprocess

import pytest

from align_genotype.jobs import align

# read-name layout for old FLOWCELL:LANE:TILE:X:Y names, e.g. HYYNKCCXX:7:1101:42581:0
READ_NAME_REGEX = r'regex:^(?<su>\w+:\d+):(?<tile>\d+):\d+:\d+'


@pytest.fixture
def set_extra_args(monkeypatch):
    def _set(extra_args: list[str]) -> None:
        def fake_config_retrieve(key: list[str], default: object = None) -> object:
            return extra_args if key == ['dupblaster', 'extra_args'] else default

        monkeypatch.setattr(align.config, 'config_retrieve', fake_config_retrieve)

    return _set


def test_dupblaster_spills_to_batch_tmpdir(set_extra_args):
    set_extra_args([])
    tokens = shlex.split(align.dedup_sort_cmd(8))
    assert tokens[tokens.index('--tmp-dir') + 1] == '$BATCH_TMPDIR'


def test_extra_args_reach_dupblaster_verbatim(set_extra_args):
    set_extra_args(['--read-name-format', READ_NAME_REGEX])
    cmd = align.dedup_sort_cmd(8)

    # -n parses without executing; the command is built from the test's own fixed args
    syntax_check = subprocess.run(  # noqa: S603
        ['bash', '-n', '-c', f'true {cmd}'],  # noqa: S607
        capture_output=True,
        text=True,
        check=False,
    )
    assert syntax_check.returncode == 0, syntax_check.stderr

    tokens = shlex.split(cmd)
    assert tokens[tokens.index('--read-name-format') + 1] == READ_NAME_REGEX
