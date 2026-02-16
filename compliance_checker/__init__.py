from collections.abc import Generator
from contextlib import contextmanager
from tempfile import NamedTemporaryFile
from typing import BinaryIO

from compliance_checker.runner import ComplianceChecker

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


def run_checker(
    ds_loc,
    checker_names,
    verbose=False,
    criteria="normal",
    output_format="text",
    output_filename="-",
    skip_checks=None,
    include_checks=None,
    options=None,
):
    """
    Simplified API to run compliance checks.

    Usage:
        from compliance_checker import run_checker
        passed, errors = run_checker('my_data.nc', 'cf')
    """
    if isinstance(checker_names, str):
        checker_names = [checker_names]

    return ComplianceChecker.run_checker(
        ds_loc=ds_loc,
        checker_names=checker_names,
        verbose=verbose,
        criteria=criteria,
        skip_checks=skip_checks,
        include_checks=include_checks,
        output_filename=output_filename,
        output_format=output_format,
        options=options,
    )
