#
# femtoanalysis/__init__.py
#

from typing import List

from pathlib import Path
from dataclasses import dataclass
from contextlib import contextmanager as _contextmanager


@_contextmanager
def lock_dir(path: Path, wait=False):
    from shutil import rmtree
    from time import sleep
    from os import getpid

    path = path.absolute()

    has_printed_warning = False
    while True:
        try:
            path.mkdir()
        except FileExistsError:
            if not wait:
                raise
            if not has_printed_warning:
                print("Waiting for directory lock.")
                has_printed_warning = True
        else:
            break

    try:
        (path / 'pid').write_text(f'{getpid()}')
        yield path
    except Exception:
        raise
    finally:
        rmtree(path)


@dataclass
class FemtoAnalysis:
    pass


@dataclass
class FemtoAnalysisJob:
    aliroot_version: str
    analyses: List[FemtoAnalysis]
