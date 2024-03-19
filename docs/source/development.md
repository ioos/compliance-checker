# Developer Notes

## pytest

When running the python test suite, there may be test errors.  Certain
tests `record` responses to remote queries for information.  If tests
fail, they will appear to continue to fail as the queries are cached.

To perform tests using fresh queries from remote services, use
`pytest --disable-vcr`.  In certain cases, clearing the cache is
also advised, use `pytest --clear-cache`.
