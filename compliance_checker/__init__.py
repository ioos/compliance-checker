from collections.abc import Generator
from contextlib import contextmanager
from tempfile import NamedTemporaryFile
from typing import BinaryIO

try:
    from ._version import __version__
except ImportError:
    __version__ = "unknown"


@contextmanager
def tempnc(data: BinaryIO) -> Generator[str, None, None]:
    """
    Create a temporary in-memory NetCDF file using a NamedTemporaryFile.
    Close the file automatically after scope is exited.

    Type aliasing and tempfile creation credit to @ocefpaf
    https://github.com/ioos/compliance-checker/pull/799#discussion_r411420587

    Parameters
    ----------
    data (bytes): raw bytes to store in the NamedTemporaryFile

    Returns
    -------
    context-managed generator
    """
    tmp = None
    try:
        tmp = NamedTemporaryFile(suffix=".nc", prefix="compliance-checker_")
        tmp.write(data)
        tmp.flush()
        yield tmp.name
    finally:
        if tmp is not None:
            tmp.close()
